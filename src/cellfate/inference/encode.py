"""Turn a ``Request`` into standardised, model-ready tensors.

Reuses the Document-1 ``Scalers`` and the Document-2 perturbation encoders so
neither of those packages has to change: SMILES -> Morgan fingerprint,
(dose, time) -> [log10 dose, log time], then the bundled scaler standardises X
and the dose/time vector.
"""

from __future__ import annotations

import numpy as np
import torch

from cellfate.common.constants import N_FINGERPRINT_BITS, Modality
from cellfate.common.errors import ContractViolation
from cellfate.common.scalers import Scalers
from cellfate.data.perturbation import encode_dose_time, encode_fingerprints


def _descriptor_to_fp(descriptor, modality: Modality) -> np.ndarray:
    """A single (2048,) float32 fingerprint from a SMILES string or a bit vector."""
    if modality != Modality.CHEM:
        raise ContractViolation(
            f"this bundle supports only chemical perturbations; got '{modality.value}'. "
            "Genetic/TF perturbations need a model trained with those input encoders."
        )
    if isinstance(descriptor, str):
        if not descriptor.strip():
            raise ContractViolation("empty SMILES descriptor")
        return encode_fingerprints([descriptor])[0].astype(np.float32)
    arr = np.asarray(descriptor, dtype=np.float32).ravel()
    if arr.shape != (N_FINGERPRINT_BITS,):
        raise ContractViolation(
            f"fingerprint descriptor must have {N_FINGERPRINT_BITS} bits; got shape {arr.shape}"
        )
    return arr


def encode_batch(reqs, scalers: Scalers):
    """(X, fp, dose_time) as standardised float32 tensors for a list of requests."""
    if not reqs:
        raise ContractViolation("empty request batch")
    n_genes = len(scalers.params.x_mean)
    X = np.asarray([r.X_raw for r in reqs], dtype=np.float32)
    if X.ndim != 2 or X.shape[1] != n_genes:
        raise ContractViolation(
            f"X_raw has {X.shape[1] if X.ndim == 2 else '?'} genes; bundle expects {n_genes}"
        )
    fp = np.stack([_descriptor_to_fp(r.u_descriptor, r.u_modality) for r in reqs])
    dt = encode_dose_time([r.dose_uM for r in reqs], [r.time_h for r in reqs])
    Xz = scalers.transform_x(X)
    dtz = scalers.transform_dose_time(dt)
    return (torch.from_numpy(Xz).float(),
            torch.from_numpy(fp).float(),
            torch.from_numpy(dtz).float())
