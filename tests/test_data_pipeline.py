"""End-to-end ETL test: run the orchestrator on injected synthetic sources and
verify every artefact validates against the cellfate.common contracts."""

from __future__ import annotations

import numpy as np
import pytest

from cellfate.common import io
from cellfate.common.constants import Split
from cellfate.common.io import ArtifactPaths
from cellfate.common.scalers import Scalers
from cellfate.common.schemas import Sample
from cellfate.data import DataConfig, QCConfig, SyntheticSource, run
from cellfate.data.splits import CONTROL_SCAFFOLD


def _sources():
    # One non-cancer source (aging path runs) + one cancer source (age masked).
    return [
        SyntheticSource(name="synth", n_lines=2, n_compounds=4,
                        n_cells_per_condition=8, n_filler_genes=150, seed=1),
        SyntheticSource(name="tahoe", n_lines=2, n_compounds=4,
                        n_cells_per_condition=8, n_filler_genes=150, seed=2),
    ]


def _cfg(tmp_path) -> DataConfig:
    return DataConfig(
        out=str(tmp_path),
        gene_panel=str(tmp_path / "gene_panel.json"),
        n_genes=160,
        qc=QCConfig(min_genes=5, max_mito_frac=0.5),
        label_tau=0.5,
        clock="random",
        deconfound=True,
        split_fracs=(0.6, 0.2, 0.1, 0.1),
        split_regimes=("scaffold", "cell_line", "both"),
        primary_regime="scaffold",
        seed=0,
    )


def _reconstruct_samples(arr: dict) -> list[Sample]:
    """Rebuild Sample objects from a shard's numpy view (proves schema validity)."""
    fp = arr["u_chem_fp"]
    out: list[Sample] = []
    for i in range(len(arr["cell_id"])):
        masked = bool(arr["age_mask"][i])
        out.append(Sample(
            cell_id=str(arr["cell_id"][i]),
            X=[float(v) for v in arr["X"][i]],
            u_modality=str(arr["u_modality"][i]),
            u_chem_fp=([int(v) for v in fp[i]] if fp is not None else None),
            dose_time=[float(v) for v in arr["dose_time"][i]],
            y_cls=[float(v) for v in arr["y_cls"][i]],
            y_age=(float(arr["y_age"][i]) if masked else None),
            age_mask=masked,
            sig_scores=[float(v) for v in arr["sig_scores"][i]],
            cell_line=str(arr["cell_line"][i]),
            pert_id=str(arr["pert_id"][i]),
            scaffold_id=str(arr["scaffold_id"][i]),
            source=str(arr["source"][i]),
        ))
    return out


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("dataset")
    cfg = _cfg(tmp)
    summary = run(cfg, sources=_sources())
    return tmp, cfg, summary


# --------------------------------------------------------------------------- #
def test_summary_and_shards_exist(built):
    tmp, cfg, summary = built
    paths = ArtifactPaths.of(tmp)
    shards = sorted(paths.shards_dir.glob("*.parquet"))
    assert len(shards) == 4                      # 2 sources x 2 lines
    assert summary["n_shards"] == 4
    assert summary["n_samples"] > 0
    assert summary["panel_size"] == 160
    # label distribution should span more than one outcome class
    assert len(summary["label_distribution"]) >= 2


def test_every_shard_row_validates_as_a_sample(built):
    tmp, _, summary = built
    paths = ArtifactPaths.of(tmp)
    total = 0
    for shard in sorted(paths.shards_dir.glob("*.parquet")):
        arr = io.shard_to_numpy(io.read_shard(shard))
        assert arr["X"].shape[1] == summary["panel_size"]
        samples = _reconstruct_samples(arr)   # raises if any row breaks the schema
        total += len(samples)
        # y_age is present exactly when age_mask is True, across the whole shard
        mask = arr["age_mask"].astype(bool)
        ya = np.asarray(arr["y_age"], dtype=float)
        assert np.isnan(ya[~mask]).all()
        assert np.isfinite(ya[mask]).all()
    assert total == summary["n_samples"]


def test_manifest_consolidates_and_counts_match(built):
    tmp, _, summary = built
    paths = ArtifactPaths.of(tmp)
    rows = io.manifest_rows(io.load_manifest(paths))
    assert len(rows) == summary["n_samples"]
    assert {r.cell_id for r in rows}.__len__() == summary["n_samples"]  # unique ids


def test_age_masking_follows_cancer_source(built):
    tmp, _, _ = built
    paths = ArtifactPaths.of(tmp)
    rows = io.manifest_rows(io.load_manifest(paths))
    tahoe = [r.age_mask for r in rows if r.source == "tahoe"]
    synth = [r.age_mask for r in rows if r.source == "synth"]
    assert tahoe and not any(tahoe)   # cancer source -> all masked
    assert synth and any(synth)       # aging source -> some valid ages


def test_all_three_regime_splits_written(built):
    tmp, cfg, summary = built
    paths = ArtifactPaths.of(tmp)
    cell_ids = {r.cell_id for r in io.manifest_rows(io.load_manifest(paths))}
    for regime in cfg.split_regimes:
        mapping = io.load_splits(paths, regime)
        assert set(mapping) == cell_ids     # every cell assigned in every regime


def test_scaffold_split_has_no_group_leakage(built):
    tmp, _, _ = built
    paths = ArtifactPaths.of(tmp)
    rows = io.manifest_rows(io.load_manifest(paths))
    mapping = io.load_splits(paths, "scaffold")
    cid_to_scaf = {r.cell_id: (r.scaffold_id or r.pert_id) for r in rows}
    by_scaffold: dict[str, set[str]] = {}
    for cid, sp in mapping.items():
        by_scaffold.setdefault(cid_to_scaf[cid], set()).add(sp)
    for scaf, sps in by_scaffold.items():
        assert len(sps) == 1, f"scaffold {scaf} leaked across {sps}"
    assert by_scaffold[CONTROL_SCAFFOLD] == {Split.TRAIN.value}


def test_scalers_fit_and_usable(built):
    tmp, _, summary = built
    paths = ArtifactPaths.of(tmp)
    scalers = Scalers.load(paths.scalers_file)
    # standardising a panel-width vector returns a finite same-shaped vector
    x = np.zeros(summary["panel_size"], dtype=np.float32)
    z = scalers.transform_x(x)
    assert z.shape == (summary["panel_size"],)
    assert np.isfinite(z).all()
    assert scalers.params.gene_panel_hash == summary["gene_panel_hash"]


def test_pipeline_is_resumable(built, tmp_path_factory):
    tmp, cfg, summary = built
    # Re-running with all chunks already done must be a no-op (no new shards),
    # and must still finalise without error.
    again = run(cfg, sources=_sources())
    paths = ArtifactPaths.of(tmp)
    assert again["n_samples"] == summary["n_samples"]
    assert again["n_shards"] == summary["n_shards"]
    assert len(sorted(paths.shards_dir.glob("*.parquet"))) == 4


# --------------------------------------------------------------------------- #
# Cell-cycle deconfounding is fit on TRAIN cells only (no eval leakage).       #
# --------------------------------------------------------------------------- #
def test_deconfounder_is_fit_on_train_cells_only(tmp_path, monkeypatch):
    """The train-only fit must recover the TRAIN cell-cycle slope and ignore the
    (deliberately opposite) slope carried by held-out cells."""
    from cellfate.common.schemas import ManifestRow
    from cellfate.data import build_dataset as bd

    # capture the re-applied y_age instead of touching real shards on disk
    captured: dict[str, np.ndarray] = {}
    monkeypatch.setattr(bd.io, "rewrite_shard_yage",
                        lambda p, y: captured.__setitem__(str(p), np.asarray(y)))

    n = 60
    cc = np.linspace(-1.0, 1.0, n)
    train_slope, held_slope = 2.0, -5.0
    aux_train = bd.ChunkAux(
        cell_ids=[f"tr_{i}" for i in range(n)], cell_line=np.array(["L"] * n),
        is_control=np.zeros(n, bool), d_age_raw=train_slope * cc, cc=cc,
        age_mask=np.ones(n, bool),
    )
    aux_held = bd.ChunkAux(
        cell_ids=[f"te_{i}" for i in range(n)], cell_line=np.array(["L"] * n),
        is_control=np.zeros(n, bool), d_age_raw=held_slope * cc, cc=cc,
        age_mask=np.ones(n, bool),
    )
    aux_by_sid = {"s_train": aux_train, "s_held": aux_held}

    def _rows(aux):
        return [ManifestRow(cell_id=c, cell_line="L", pert_id="p", scaffold_id="sc",
                            source="synth", age_mask=True, shard_id="s", row_idx=i)
                for i, c in enumerate(aux.cell_ids)]
    rows = _rows(aux_train) + _rows(aux_held)
    splits = {"scaffold": {**{c: Split.TRAIN.value for c in aux_train.cell_ids},
                           **{c: Split.TEST.value for c in aux_held.cell_ids}}}

    cfg = DataConfig(out=str(tmp_path), gene_panel=str(tmp_path / "panel.json"),
                     deconfound=True, primary_regime="scaffold")
    coef = bd._deconfound_train_only(cfg, ArtifactPaths.of(tmp_path), rows, splits, aux_by_sid)

    # slope reflects TRAIN only (~+2.0); a pooled fit over both arms would be ~ -1.5
    assert abs(coef[0] - train_slope) < 0.2
    # the identical transform was re-applied to every shard (train and held-out)
    assert set(captured) == {str(ArtifactPaths.of(tmp_path).shard_file(s))
                             for s in ("s_train", "s_held")}
    # and it removes the cell-cycle correlation on the held-out shard too
    held = captured[str(ArtifactPaths.of(tmp_path).shard_file("s_held"))]
    assert abs(np.polyfit(cc, held, 1)[0] - (held_slope - train_slope)) < 0.2


# --------------------------------------------------------------------------- #
# Crash-resume: ProgressTracker skips completed chunks and the resumed build   #
# still fits the train-only deconfounder correctly (via on-disk sidecars).     #
# --------------------------------------------------------------------------- #
class _CrashAfter(SyntheticSource):
    """A synthetic source that raises once on a chosen chunk, then behaves."""
    def __init__(self, *a, crash_cid=None, **k):
        super().__init__(*a, **k)
        self._crash_cid, self._crashed, self.fetched = crash_cid, False, []

    def fetch(self, chunk):
        if chunk["id"] == self._crash_cid and not self._crashed:
            self._crashed = True
            raise RuntimeError(f"simulated drop on {chunk['id']}")
        self.fetched.append(chunk["id"])
        return super().fetch(chunk)


def test_crash_then_resume_skips_done_and_stays_correct(tmp_path):
    root = tmp_path / "ds"
    def cfg():
        return DataConfig(out=str(root), gene_panel=str(root / "panel.json"), n_genes=120,
                          qc=QCConfig(min_genes=5, max_mito_frac=0.5),
                          split_regimes=("scaffold", "cell_line", "both"),
                          primary_regime="scaffold", seed=0)
    kw = dict(n_lines=3, n_compounds=6, n_cells_per_condition=8, seed=1)
    all_ids = [c["id"] for c in SyntheticSource(name="synth", n_scaffold_families=7, **kw).plan()]
    crash_on = all_ids[1]

    # run 1: crashes partway through
    s1 = _CrashAfter(name="synth", n_scaffold_families=7, crash_cid=crash_on, **kw)
    with pytest.raises(RuntimeError):
        run(cfg(), sources=[s1])
    paths = ArtifactPaths.of(root)
    from cellfate.common.progress import ProgressTracker
    done1 = {c for c in all_ids if ProgressTracker(paths.progress_file).is_done(c)}
    assert done1 and crash_on not in done1                     # some done, crashed one not

    # run 2: resume
    s2 = _CrashAfter(name="synth", n_scaffold_families=7, crash_cid=None, **kw)
    summary = run(cfg(), sources=[s2])
    assert set(s2.fetched).isdisjoint(done1)                   # done chunks were skipped
    assert set(s2.fetched) == set(all_ids) - done1             # only the rest re-fetched
    assert summary["n_shards"] == len(all_ids)                 # dataset complete
    assert not (paths.root / "_cc_cache").exists()             # sidecar cache cleaned up
    # deconfounder still control-relative after the resume
    ctrl = []
    for sh in sorted(paths.shards_dir.glob("*.parquet")):
        a = io.shard_to_numpy(io.read_shard(sh))
        m = a["age_mask"].astype(bool)
        isc = np.array([s == CONTROL_SCAFFOLD for s in a["scaffold_id"]])
        ctrl += [a["y_age"][i] for i in range(len(m)) if m[i] and isc[i]]
    assert abs(float(np.mean(ctrl))) < 1e-6
