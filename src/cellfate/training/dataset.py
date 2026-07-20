"""Turn dataset shards into training tensors (Document 3, S3).

Reads the Document-2 artefacts (shards + splits + scalers) and materialises the
rows of one split into a ``TensorDataset``. The frozen scalers are applied here
exactly as they will be at inference, so train/serve see identical inputs. ΔAge
is gated by ``age_mask`` (NaN ages become 0 and are ignored by the masked loss).
"""

from __future__ import annotations

import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset

from cellfate.common import io
from cellfate.common.errors import DataSourceError
from cellfate.common.io import ArtifactPaths
from cellfate.common.scalers import Scalers

# TensorDataset column order (used throughout the training package).
X_I, FP_I, DT_I, YC_I, YA_I, AM_I, DONOR_I = range(7)

# Donor (``cell_line``) string -> stable integer code. Module level so that every
# split of a run shares one encoding; ``setdefault`` keeps a donor's code fixed for
# the process lifetime. Consumed by the inner leave-one-donor-out in xdonor_calib.
#
# The CODE VALUES depend on the order donors are first seen, so they can differ between
# processes (a name always maps to the same code within one). That is safe only because every
# statistic pooled in xdonor_calib is order-invariant -- quantile, median, sort, and a mean over
# all rows. **Any new pooled statistic must also be order-invariant**, or `run()` stops being
# reproducible (tests/test_correctness.py::test_run_is_reproducible).
DONOR_VOCAB: dict[str, int] = {}


def _donor_code(name: object) -> int:
    return DONOR_VOCAB.setdefault(str(name), len(DONOR_VOCAB))


def load_split_tensors(paths: ArtifactPaths, scalers: Scalers, regime: str,
                       split: str) -> TensorDataset:
    """Materialise one split of one regime as a TensorDataset of 7 columns:
    (X, fingerprint, dose_time, y_cls, y_age, age_mask, donor)."""
    wanted = {cid for cid, sp in io.load_splits(paths, regime).items() if sp == split}
    xs, fps, dts, ycls, yage, amask, donors = [], [], [], [], [], [], []
    pert_width: int | None = None
    for shard in sorted(paths.shards_dir.glob("*.parquet")):
        arr = io.shard_to_numpy(io.read_shard(shard))
        # perturbation input: chem fingerprint OR TF-cocktail multi-hot (one per dataset)
        pert = arr["u_chem_fp"] if arr["u_chem_fp"] is not None else arr["u_tf_emb"]
        if pert is None:
            raise DataSourceError("perturbation features missing (neither u_chem_fp nor u_tf_emb)")
        pert_width = pert.shape[1]
        ids = arr["cell_id"]
        keep = np.fromiter((c in wanted for c in ids), bool, len(ids))
        if not keep.any():
            continue
        am = arr["age_mask"][keep].astype(bool)
        ya = np.where(am, np.asarray(arr["y_age"][keep], np.float32), 0.0).astype(np.float32)
        xs.append(scalers.transform_x(arr["X"][keep]))
        fps.append(pert[keep].astype(np.float32))
        dts.append(scalers.transform_dose_time(arr["dose_time"][keep]))
        ycls.append(arr["y_cls"][keep].astype(np.float32))
        yage.append(ya)
        amask.append(am.astype(np.float32))
        # donor identity, for the inner leave-one-donor-out that fits calibration
        donors.append(np.fromiter((_donor_code(c) for c in np.asarray(arr["cell_line"])[keep]),
                                  np.int64, int(keep.sum())))

    g = len(scalers.params.x_mean)
    if not xs:  # empty split -> empty (but well-shaped) tensors
        from cellfate.common.constants import N_DOSE_TIME, N_FINGERPRINT_BITS
        w = pert_width if pert_width is not None else N_FINGERPRINT_BITS
        z = torch.empty(0)
        return TensorDataset(torch.empty(0, g), torch.empty(0, w),
                             torch.empty(0, N_DOSE_TIME), torch.empty(0, 3), z, z,
                             torch.empty(0, dtype=torch.long))
    return TensorDataset(
        torch.from_numpy(np.vstack(xs)).float(),
        torch.from_numpy(np.vstack(fps)).float(),
        torch.from_numpy(np.vstack(dts)).float(),
        torch.from_numpy(np.vstack(ycls)).float(),
        torch.from_numpy(np.concatenate(yage)).float(),
        torch.from_numpy(np.concatenate(amask)).float(),
        torch.from_numpy(np.concatenate(donors)).long(),
    )


def loader(ds: TensorDataset, batch_size: int, shuffle: bool) -> DataLoader:
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle, drop_last=False)
