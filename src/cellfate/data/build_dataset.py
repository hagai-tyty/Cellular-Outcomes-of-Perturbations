"""ETL orchestrator (Document 2, S14).

End-to-end flow, resumable at chunk granularity:

    plan -> [per chunk: fetch -> QC -> normalise -> signatures -> soft labels
             -> cell-cycle -> dAge (+mask) -> deconfound -> panel X -> encode
             -> assemble -> write shard + manifest part] -> consolidate manifest
    -> splits (all regimes) -> fit scalers on the primary regime's train -> summary

``run`` is plain-Python and dependency-injectable (pass ``sources`` / ``clock``)
so it is fully testable; ``cli`` is the Hydra entry point used in production.
"""

from __future__ import annotations

import shutil
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from cellfate.common import constants as C
from cellfate.common import io
from cellfate.common.constants import Modality, Split
from cellfate.common.errors import ConfigError
from cellfate.common.io import ArtifactPaths
from cellfate.common.logging import get_logger, log_event
from cellfate.common.scalers import Scalers
from cellfate.common.schemas import ManifestRow

from .aging import (
    AgingClock,
    LinearClock,
    delta_age,
    recenter_on_control_arrays,
)
from .assemble import assemble_samples
from .chunking import ProgressTracker, plan_all
from .labels import fate_labels
from .normalize import GenePanel, fit_gene_panel, normalize_counts, to_panel_matrix
from .perturbation import (
    encode_dose_time,
    encode_factor_sets,
    encode_fingerprints,
    resolve_scaffolds,
)
from .proliferation import cell_cycle_score, deconfound_age, fit_deconfounder
from .qc import QCConfig, apply_qc
from .signatures import signature_scores
from .sources import SOURCE_REGISTRY, DataSource
from .splits import make_splits

log = get_logger("cellfate.data")

# The fate-label markers are held OUT of the model panel (anti-circularity); the
# cell-cycle score is computed from the full profile, so no gene is forced in.
_PANEL_EXCLUDE: tuple[str, ...] = C.LABEL_HOLDOUT


@dataclass
class ChunkAux:
    """Per-cell material a chunk hands back so a **train-fit** cell-cycle
    deconfounder can be applied after splits are known (see ``process_chunk``)."""
    cell_ids: list[str]
    cell_line: np.ndarray
    is_control: np.ndarray
    d_age_raw: np.ndarray
    cc: np.ndarray
    age_mask: np.ndarray


@dataclass
class DataConfig:
    out: str
    gene_panel: str
    n_genes: int = C.DEFAULT_N_GENES
    qc: QCConfig = field(default_factory=QCConfig)
    label_tau: float = 1.0
    clock: str = "random"        # "random" (synthetic/smoke) or path to a fitted clock .json
    deconfound: bool = True
    modality: str = "chem"       # "chem" (fingerprint) or "tf" (OSKM-style factor cocktail)
    split_fracs: tuple[float, ...] = (0.7, 0.1, 0.1, 0.1)
    split_regimes: tuple[str, ...] = ("scaffold", "cell_line", "both")
    primary_regime: str = "scaffold"
    holdout_cell_lines: tuple[str, ...] = ()   # "holdout" regime: these cell lines -> test
    fate_test_line: str = ""                   # "line_holdout" regime: hold out a slice of this line
    fate_test_frac: float = 0.15               # fraction of fate_test_line -> test
    seed: int = 0
    panel_ref_chunks: int = 1
    scaler_max_cells: int = 200_000
    continue_on_error: bool = False
    source_specs: tuple[dict, ...] = ()


# --------------------------------------------------------------------------- #
# Per-chunk pipeline                                                          #
# --------------------------------------------------------------------------- #
def process_chunk(src, chunk, panel, clock: AgingClock, cfg: DataConfig):
    """Run the full transform for one chunk.

    Returns ``(samples, aux)`` where ``samples`` carry the *raw* control-relative
    ΔAge (no cell-cycle deconfounding applied here). ``aux`` is the per-cell
    material needed to re-apply a **train-fit** deconfounder later, or ``None``
    when deconfounding is off or the chunk has no age-valid cells. Fitting the
    deconfounder on train cells only requires the split assignment, which does
    not exist until every chunk has been read -- hence the deferred two-pass.
    """
    raw = apply_qc(src.fetch(chunk), cfg.qc)
    if len(raw.obs) == 0:
        return [], None

    norm = normalize_counts(raw.counts)
    sig = signature_scores(norm, raw.genes)
    y_cls = fate_labels(norm, raw.genes, raw.obs, cfg.label_tau)
    cc = cell_cycle_score(norm, raw.genes)
    x_panel = to_panel_matrix(norm, raw.genes, panel)

    # the clock consumes the FULL profile (its own gene panel), NOT the 2000-HVG
    # model input x_panel -- so aging genes filtered out of the HVG panel still
    # reach the clock. The model still trains on x_panel below.
    d_age, age_mask = delta_age(clock, norm, raw.genes, raw.obs, raw.source)
    cell_ids = raw.obs["cell_id"].tolist()
    aux: ChunkAux | None = None
    if cfg.deconfound and age_mask.any():
        aux = ChunkAux(
            cell_ids=cell_ids,
            cell_line=raw.obs["cell_line"].to_numpy().copy(),
            is_control=raw.obs["is_control"].to_numpy().astype(bool),
            d_age_raw=np.asarray(d_age, dtype=np.float64).copy(),
            cc=np.asarray(cc, dtype=np.float64).copy(),
            age_mask=age_mask.copy(),
        )

    smiles = raw.obs["smiles"].tolist()
    pert_ids = raw.obs["pert_id"].tolist()
    doses = raw.obs["dose_uM"].to_numpy()
    dose_time = encode_dose_time(doses, raw.obs["time_h"].to_numpy())
    is_tf = cfg.modality == "tf"
    if is_tf:
        # TF cocktail (OSKM...): multi-hot factor vector, not a fingerprint
        fingerprints = np.zeros((len(cell_ids), 0), dtype=np.uint8)
        tf_emb = encode_factor_sets(pert_ids, doses.tolist())
        scaffold_id = raw.obs["scaffold_id"].tolist()
    else:
        fingerprints = encode_fingerprints(smiles)
        tf_emb = None
        scaffold_id = resolve_scaffolds(smiles, pert_ids, raw.obs["scaffold_id"].tolist())
    samples = assemble_samples(
        cell_ids=cell_ids,
        x_panel=x_panel,
        fingerprints=fingerprints,
        tf_emb=tf_emb,
        modality=Modality.TF if is_tf else Modality.CHEM,
        dose_time=dose_time,
        y_cls=y_cls,
        y_age=d_age,
        age_mask=age_mask,
        sig_scores=sig,
        cell_line=raw.obs["cell_line"].tolist(),
        pert_id=pert_ids,
        scaffold_id=scaffold_id,
        source=raw.source,
    )
    return samples, aux


# --------------------------------------------------------------------------- #
# Panel + clock + source construction                                         #
# --------------------------------------------------------------------------- #
def load_or_fit_panel(cfg: DataConfig, work) -> GenePanel:
    panel_path = Path(cfg.gene_panel)
    if panel_path.exists():
        panel = GenePanel.load(panel_path)
        log_event(log, "panel.loaded", path=str(panel_path), n=len(panel), hash=panel.hash())
        return panel

    pooled, genes = [], None
    for src, chunk in work[: cfg.panel_ref_chunks]:
        raw = apply_qc(src.fetch(chunk), cfg.qc)
        pooled.append(normalize_counts(raw.counts))
        genes = raw.genes
    panel = fit_gene_panel(np.vstack(pooled), genes, n_top=cfg.n_genes, must_exclude=_PANEL_EXCLUDE)
    panel.save(panel_path)
    log_event(log, "panel.fit", path=str(panel_path), n=len(panel), hash=panel.hash())
    return panel


def build_clock(cfg: DataConfig, panel) -> AgingClock:
    """Resolve ``cfg.clock`` to a clock. Either ``'random'`` (explicit, for
    synthetic/smoke runs -- meaningless ages) or a path to a fitted weights JSON
    (see scripts/fit_clock.py). Anything else fails loud -- no silent fallback."""
    spec = str(cfg.clock)
    if spec == "random":
        return LinearClock.random(panel, seed=cfg.seed)
    if Path(spec).exists():
        return LinearClock.from_json(spec)
    raise ConfigError(
        f"clock={spec!r}: not 'random' and no weights file exists at that path. "
        "Fit a real clock on an age-labelled dataset (scripts/fit_clock.py, e.g. on "
        "GSE113957 human fibroblasts) and point clock: at the resulting .json, or set "
        "clock='random' for synthetic/smoke runs (its ages are not meaningful)."
    )


def build_sources(cfg: DataConfig) -> list[DataSource]:
    sources: list[DataSource] = []
    for spec in cfg.source_specs:
        spec = dict(spec)
        key = spec.pop("name")            # registry key (e.g. "synthetic")
        if "source_name" in spec:         # optional per-instance name override
            spec["name"] = spec.pop("source_name")
        if key not in SOURCE_REGISTRY:
            raise ValueError(f"unknown source {key!r}; have {list(SOURCE_REGISTRY)}")
        sources.append(SOURCE_REGISTRY[key](**spec))
    return sources


# --------------------------------------------------------------------------- #
# Orchestrator                                                                #
# --------------------------------------------------------------------------- #
def run(cfg: DataConfig, sources: list[DataSource] | None = None,
        clock: AgingClock | None = None) -> dict:
    """Execute the ETL and return a summary dict."""
    paths = ArtifactPaths.of(cfg.out)
    sources = sources if sources is not None else build_sources(cfg)
    if not sources:
        raise ValueError("no data sources configured")
    work = plan_all(sources)

    panel = load_or_fit_panel(cfg, work)
    clock = clock if clock is not None else build_clock(cfg, panel)
    tracker = ProgressTracker(paths.progress_file)

    label_counts: Counter[str] = Counter()
    n_age_labeled = 0

    for src, chunk in work:
        cid = chunk["id"]
        if tracker.is_done(cid):
            continue
        sid = io.sanitize_id(cid)
        try:
            samples, aux = process_chunk(src, chunk, panel, clock, cfg)
            if not samples:
                tracker.mark_done(cid, 0)
                continue
            io.write_shard(paths.shard_file(sid), samples)
            if aux is not None:
                # persist per-chunk cell-cycle data so a resumed run can still
                # fit the deconfounder on TRAIN cells only (survives crashes).
                _write_cc_sidecar(paths, sid, aux)
            io.write_manifest_part(
                paths, sid, [ManifestRow.from_sample(s, sid, j) for j, s in enumerate(samples)]
            )
            for s in samples:
                label_counts[C.IDX_TO_CLASS[int(np.argmax(s.y_cls))]] += 1
                n_age_labeled += int(s.age_mask)
            tracker.mark_done(cid, len(samples))
            log_event(log, "chunk.done", chunk=cid, n=len(samples))
        except Exception as exc:  # noqa: BLE001 - recorded for resume
            tracker.mark_failed(cid, repr(exc))
            log_event(log, "chunk.failed", chunk=cid, err=repr(exc))
            if not cfg.continue_on_error:
                raise

    io.consolidate_manifest(paths)
    rows = io.manifest_rows(io.load_manifest(paths))
    splits = make_splits(rows, tuple(cfg.split_fracs), tuple(cfg.split_regimes), cfg.seed,
                         holdout_cell_lines=tuple(cfg.holdout_cell_lines),
                         fate_test_line=cfg.fate_test_line, fate_test_frac=cfg.fate_test_frac)
    for regime, mapping in splits.items():
        io.write_splits(paths, regime, mapping)

    # Cell-cycle deconfounding: fit on the primary regime's TRAIN age-valid cells
    # only, then re-apply the *same* transform to every shard (no eval leakage).
    # Sidecars are read from disk so this is correct after a resumed build too.
    aux_by_sid = _load_cc_sidecars(paths) if cfg.deconfound else {}
    coef = _deconfound_train_only(cfg, paths, rows, splits, aux_by_sid)

    _fit_scalers(cfg, paths, panel, splits, coef)
    _clear_cc_cache(paths)

    summary = {
        "n_samples": len(rows),
        "n_shards": tracker.n_done,
        "n_age_labeled": n_age_labeled,
        "gene_panel_hash": panel.hash(),
        "panel_size": len(panel),
        "label_distribution": dict(label_counts),
        "regimes": list(splits),
        "primary_regime": cfg.primary_regime,
        "split_sizes": {sp: int(v) for sp, v in
                        Counter(splits[cfg.primary_regime].values()).items()},
        "created_at": time.time(),
    }
    io.write_json(paths.root / "dataset_summary.json", summary)
    log_event(log, "dataset.done", **{k: summary[k] for k in ("n_samples", "n_shards")})
    return summary


def _cc_cache_dir(paths) -> Path:
    d = paths.root / "_cc_cache"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write_cc_sidecar(paths, sid: str, aux: ChunkAux) -> None:
    """Persist a chunk's cell-cycle material next to its shard so a train-only
    deconfounder can be fit after splits -- and after a resumed build."""
    np.savez(
        _cc_cache_dir(paths) / f"{sid}.npz",
        cell_ids=np.asarray(aux.cell_ids).astype("U"),
        cell_line=np.asarray(aux.cell_line).astype("U"),
        is_control=np.asarray(aux.is_control, dtype=bool),
        d_age_raw=np.asarray(aux.d_age_raw, dtype=np.float64),
        cc=np.asarray(aux.cc, dtype=np.float64),
        age_mask=np.asarray(aux.age_mask, dtype=bool),
    )


def _load_cc_sidecars(paths) -> dict[str, ChunkAux]:
    out: dict[str, ChunkAux] = {}
    d = paths.root / "_cc_cache"
    if not d.exists():
        return out
    for p in sorted(d.glob("*.npz")):
        z = np.load(p, allow_pickle=False)
        out[p.stem] = ChunkAux(
            cell_ids=[str(c) for c in z["cell_ids"]], cell_line=z["cell_line"],
            is_control=z["is_control"], d_age_raw=z["d_age_raw"],
            cc=z["cc"], age_mask=z["age_mask"],
        )
    return out


def _clear_cc_cache(paths) -> None:
    d = paths.root / "_cc_cache"
    if d.exists():
        shutil.rmtree(d, ignore_errors=True)


def _deconfound_train_only(cfg, paths, rows, splits, aux_by_sid) -> tuple[float, float]:
    """Fit ``ΔAge ~ a*cc + b`` on the primary regime's TRAIN age-valid cells,
    then re-apply ``deconfound_age`` + control re-centring to *every* shard.

    Returns the fitted ``(a, b)`` (stored in the scalers as metadata). Fitting on
    train only -- and applying the identical transform to val/test -- is what
    keeps the cell-cycle correction leak-free, mirroring the scaler fit.
    """
    if not cfg.deconfound or not aux_by_sid:
        return (0.0, 0.0)

    # Integrity check: every age-valid cell in the manifest must have a cached
    # sidecar (they are written before a chunk is marked done, so a resumed build
    # still has them). A shortfall means the cache was deleted or corrupted.
    manifest_age_valid = sum(int(r.age_mask) for r in rows)
    aux_age_valid = sum(int(a.age_mask.sum()) for a in aux_by_sid.values())
    if manifest_age_valid > aux_age_valid:
        raise ConfigError(
            "cell-cycle deconfounding cache is incomplete (manifest has "
            f"{manifest_age_valid} age-valid cells but the sidecars cover only "
            f"{aux_age_valid}). Rebuild into a clean output directory."
        )

    train_ids = {cid for cid, sp in splits[cfg.primary_regime].items()
                 if sp == Split.TRAIN.value}
    d_tr, cc_tr = [], []
    for aux in aux_by_sid.values():
        for i, cell in enumerate(aux.cell_ids):
            if aux.age_mask[i] and cell in train_ids:
                d_tr.append(aux.d_age_raw[i])
                cc_tr.append(aux.cc[i])
    coef = (fit_deconfounder(np.asarray(d_tr), np.asarray(cc_tr))
            if len(d_tr) >= 2 else (0.0, 0.0))

    # Pass 2: re-apply the single train-fit transform to every shard.
    for sid, aux in aux_by_sid.items():
        d = deconfound_age(aux.d_age_raw, aux.cc, coef)
        m = aux.age_mask
        y = np.full(d.shape[0], np.nan, dtype=np.float64)
        if m.any():
            y[m] = recenter_on_control_arrays(d[m], aux.cell_line[m], aux.is_control[m])
        io.rewrite_shard_yage(paths.shard_file(sid), y)
    return coef


def _fit_scalers(cfg, paths, panel, splits, coef) -> None:
    """Fit normalisation on the primary regime's TRAIN rows only (no leakage)."""
    train_ids = {cid for cid, sp in splits[cfg.primary_regime].items()
                 if sp == Split.TRAIN.value}
    xs, dts = [], []
    for shard in sorted(paths.shards_dir.glob("*.parquet")):
        arr = io.shard_to_numpy(io.read_shard(shard))
        keep = np.array([c in train_ids for c in arr["cell_id"]], dtype=bool)
        if keep.any():
            xs.append(arr["X"][keep])
            dts.append(arr["dose_time"][keep])
    if not xs:
        raise ValueError("no TRAIN rows found to fit scalers on")
    x_train = np.vstack(xs)
    dt_train = np.vstack(dts)
    if len(x_train) > cfg.scaler_max_cells:
        rng = np.random.default_rng(cfg.seed)
        idx = rng.choice(len(x_train), cfg.scaler_max_cells, replace=False)
        x_train, dt_train = x_train[idx], dt_train[idx]
    coef = tuple(float(x) for x in coef)
    Scalers.fit(x_train, dt_train, panel, proliferation_coef=coef).save(paths.scalers_file)


# --------------------------------------------------------------------------- #
# CLI (Hydra)                                                                 #
# --------------------------------------------------------------------------- #
def _config_from_omegaconf(cfg) -> tuple[DataConfig, list[DataSource] | None]:
    """Map a composed Hydra config to a DataConfig (+ optional injected sources)."""
    d = cfg.data
    qc = QCConfig(
        min_genes=int(d.qc.min_genes),
        max_mito_frac=float(d.qc.max_mito_frac),
        max_counts=(None if d.qc.get("max_counts") is None else float(d.qc.max_counts)),
    )
    # source_specs drives source construction; fall back to `sources` (drop weight).
    specs = d.get("source_specs")
    if specs is None:
        specs = [{k: v for k, v in dict(s).items() if k != "weight"}
                 for s in d.get("sources", [])]
    return (
        DataConfig(
            out=d.out,
            gene_panel=d.gene_panel,
            n_genes=int(d.get("n_genes", cfg.get("model", {}).get("g", C.DEFAULT_N_GENES))),
            qc=qc,
            label_tau=float(d.labels.tau),
            clock=str(d.labels.clock),
            deconfound=bool(d.labels.deconfound_proliferation),
            split_fracs=tuple(d.splits.fracs),
            split_regimes=tuple(d.splits.regimes),
            primary_regime=str(d.splits.get("primary", "scaffold")),
            seed=int(cfg.seed),
            source_specs=tuple(dict(s) for s in specs),
        ),
        None,
    )


def cli() -> None:  # pragma: no cover - exercised in production via Hydra
    try:
        import hydra
        from omegaconf import DictConfig
    except ImportError as exc:
        raise ConfigError("hydra-core/omegaconf required for the CLI") from exc

    config_dir = str(Path(__file__).resolve().parents[3] / "configs")

    @hydra.main(version_base=None, config_path=config_dir, config_name="config")
    def _main(cfg: DictConfig) -> None:
        data_cfg, injected = _config_from_omegaconf(cfg)
        summary = run(data_cfg, sources=injected)
        log_event(log, "cli.done", **{k: summary[k] for k in ("n_samples", "n_shards")})

    _main()


if __name__ == "__main__":
    cli()
