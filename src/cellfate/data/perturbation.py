"""Perturbation encoding (Document 2, S7).

Chemical perturbations become a Morgan fingerprint (rdkit, lazy) plus an encoded
dose/time vector. When rdkit is unavailable (tests / CI), a deterministic hashed
fingerprint stands in so the pipeline still runs; production runs should install
the ``data`` extra so real Morgan fingerprints and Bemis-Murcko scaffolds are
used.
"""

from __future__ import annotations

import hashlib
import warnings

import numpy as np

from cellfate.common import constants as C

_DOSE_FLOOR = 1e-4   # uM, used for vehicle/control (dose = 0)
_TIME_FLOOR = 1e-2   # h

_RDKIT_WARNED = False


def _have_rdkit() -> bool:
    try:
        import rdkit  # noqa: F401
        return True
    except ImportError:
        return False


def hashed_fingerprint(smiles: str, n_bits: int = C.N_FINGERPRINT_BITS,
                       n_features: int = 48) -> np.ndarray:
    """Deterministic folded fingerprint from a SMILES string (rdkit-free).

    Sets ~``n_features`` bits chosen by hashing; mimics a sparse folded Morgan
    fingerprint closely enough to exercise the pipeline. Not chemically
    meaningful -- a stand-in for environments without rdkit.
    """
    bits = np.zeros(n_bits, dtype=np.uint8)
    if not smiles:
        return bits  # vehicle / control: empty fingerprint
    for i in range(n_features):
        h = hashlib.sha256(f"{smiles}|{i}".encode()).digest()
        bits[int.from_bytes(h[:4], "big") % n_bits] = 1
    return bits


def morgan_fingerprint(smiles: str, n_bits: int = C.N_FINGERPRINT_BITS,
                       radius: int = C.MORGAN_RADIUS) -> np.ndarray:
    """Real Morgan fingerprint via rdkit; falls back to a hashed one if absent."""
    global _RDKIT_WARNED
    if not smiles:
        return np.zeros(n_bits, dtype=np.uint8)
    if not _have_rdkit():
        if not _RDKIT_WARNED:
            warnings.warn("rdkit not installed; using hashed fingerprints "
                          "(install the 'data' extra for real Morgan FPs).", stacklevel=2)
            _RDKIT_WARNED = True
        return hashed_fingerprint(smiles, n_bits)
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return hashed_fingerprint(smiles, n_bits)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros(n_bits, dtype=np.uint8)
    for bit in fp.GetOnBits():
        arr[bit] = 1
    return arr


def bemis_murcko_scaffold(smiles: str) -> str | None:
    """Bemis-Murcko scaffold SMILES via rdkit; ``None`` if unavailable/invalid."""
    if not smiles or not _have_rdkit():
        return None
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))


def encode_fingerprints(smiles: list[str]) -> np.ndarray:
    """Encode a list of SMILES into an (N, n_bits) uint8 matrix.

    Fingerprints are computed once per unique SMILES, then broadcast -- many
    cells share the same compound.
    """
    uniq = {s: morgan_fingerprint(s) for s in set(smiles)}
    return np.stack([uniq[s] for s in smiles]).astype(np.uint8)


def encode_dose_time(dose_uM, time_h) -> np.ndarray:
    """Encode dose/time as [log10(dose_uM), log(time_h)], (N, 2) float32.

    Vehicle/control rows (dose 0) are floored so the log is finite.
    """
    d = np.log10(np.maximum(np.asarray(dose_uM, dtype=np.float64), _DOSE_FLOOR))
    t = np.log(np.maximum(np.asarray(time_h, dtype=np.float64), _TIME_FLOOR))
    return np.stack([d, t], axis=1).astype(np.float32)


def resolve_scaffolds(smiles: list[str], pert_ids: list[str],
                      provided: list[str | None] | None = None) -> list[str]:
    """Pick a scaffold id per row for the leave-scaffold-out split.

    Preference: an explicitly provided scaffold id (e.g. precomputed or
    ``"CONTROL"``) > a freshly computed Bemis-Murcko scaffold > the compound id
    itself (so each compound is at least its own group). Never empty -- the
    Sample contract requires a scaffold id for chemical perturbations.
    """
    out: list[str] = []
    for i, (smi, pid) in enumerate(zip(smiles, pert_ids, strict=True)):
        prov = provided[i] if provided is not None else None
        if prov:
            out.append(prov)
            continue
        scaf = bemis_murcko_scaffold(smi)
        out.append(scaf or pid)
    return out


# --------------------------------------------------------------------------- #
# Transcription-factor (cocktail) encoding -- the TF/genetic modality.         #
# --------------------------------------------------------------------------- #
def factors_of(pert_id: str) -> tuple[str, ...]:
    """Resolve a cocktail id to its canonical TF symbols.

    Recognises named cocktails ("OSKM"), '+'-separated symbol lists
    ("POU5F1+SOX2+KLF4"), and single-letter shorthand ("OSK"). Controls /
    unknown ids resolve to no factors.
    """
    if not pert_id or pert_id.lower() in ("control", "vehicle", "dmso", "none"):
        return ()
    key = pert_id.strip()
    if key in C.TF_COCKTAILS:
        return C.TF_COCKTAILS[key]
    if "+" in key:
        syms = tuple(s.strip() for s in key.split("+"))
        return tuple(s for s in syms if s in C.TF_VOCAB)
    if key in C.TF_VOCAB:
        return (key,)
    # single-letter shorthand (e.g. OSKM) if every char is a known letter
    if all(ch in C.TF_LETTERS for ch in key):
        return tuple(dict.fromkeys(C.TF_LETTERS[ch] for ch in key))
    return ()


def encode_factor_sets(pert_ids: list[str], doses: list[float] | None = None) -> np.ndarray:
    """Multi-hot TF-vocabulary matrix (N, len(TF_VOCAB)), each factor scaled by dose.

    This is the ``u_tf_emb`` model input for TF perturbations: OSKM and OSK share
    columns, so the encoder learns factor identity and generalises across cocktails.
    """
    n = len(pert_ids)
    doses = doses if doses is not None else [1.0] * n
    vidx = {t: i for i, t in enumerate(C.TF_VOCAB)}
    out = np.zeros((n, len(C.TF_VOCAB)), dtype=np.float32)
    for r, (p, d) in enumerate(zip(pert_ids, doses, strict=True)):
        for tf in factors_of(p):
            out[r, vidx[tf]] = float(d) if d else 1.0
    return out
