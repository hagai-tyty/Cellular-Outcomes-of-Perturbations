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
    age_label_policy,
    census_warnings,
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
    # STAGE 1.5.3 C-I. Cells whose ΔAge is COMPUTABLE, which is not the same question as whose
    # label the age head is allowed to TRAIN on. It differs from ``age_mask`` by exactly the
    # ``dataset_policy`` exclusions (C-1's ``AGE_MASKED_DATASETS``): a cell withheld from training
    # still has a perfectly good ΔAge, so it must still inform the cell-cycle deconfounder and the
    # control re-centring. Using ``age_mask`` for those made ``y_age`` itself depend on the
    # training policy, which is what confounded step 6's first run -- see results/STEP6_REPORT.md.
    deconfound_mask: np.ndarray = None  # type: ignore[assignment]
    # ARM C (step-6 follow-up). Cells whose ΔAge label is to be SHUFFLED — the label-permutation
    # control that separates "HFF's labels carry information" from "75 labels is simply too few".
    # Empty unless `DataConfig.age_shuffle_datasets` is set.
    shuffle_mask: np.ndarray = None  # type: ignore[assignment]
    # ARM D. The stratum a cell belongs to, `f"{cell_line}|{time_h}"`. With
    # `age_shuffle_strata=True` the permutation runs WITHIN each stratum, so the between-stratum
    # trajectory (day -> ΔAge, rho = -0.905) survives while cell-level pairing is destroyed.
    stratum: np.ndarray = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.deconfound_mask is None:      # pre-C-I callers keep the old meaning
            self.deconfound_mask = self.age_mask
        if self.shuffle_mask is None:
            self.shuffle_mask = np.zeros(len(self.cell_ids), dtype=bool)
        if self.stratum is None:      # one global stratum == arm C's unstratified permutation
            self.stratum = np.full(len(self.cell_ids), "__all__", dtype=object)


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
    harmonize: bool = False                    # cross-modality control-anchoring + Gill Projection
    harmonize_ref_dataset: str = "gill_bulk"   # dataset whose scale the clock reads (projection ref)
    seed: int = 0
    panel_ref_chunks: int = 1
    scaler_max_cells: int = 200_000
    continue_on_error: bool = False
    source_specs: tuple[dict, ...] = ()
    # Stage 1.5.3 C-2: mask donors outside the clock's fitted `meta.age_range`.
    # OFF by default because turning it on MOVES LABELS -- N2 and N3 are donor_age 0 and the
    # shipped clock's range starts at 1.0, so 30 of the 75 non-HFF training labels would go.
    # Enabling it is a pre-registered change with its own bar, never a side effect.
    enforce_clock_age_range: bool = False
    # ARM C (step-6 follow-up): permute ΔAge labels among these datasets' cells instead of using or
    # masking them. Same cells, same label count, same value distribution -- only the cell<->label
    # pairing destroyed. Empty = OFF, so arms A and B are unaffected. Set to {"hff_sc"} for arm C.
    # The seed is recorded in `dataset_summary.json`; the permutation is global across shards.
    age_shuffle_datasets: frozenset[str] = frozenset()
    age_shuffle_seed: int = 0
    # ARM D: permute WITHIN (cell_line, time_h) strata instead of globally. Arm C's global shuffle
    # destroys the between-timepoint trajectory as well as cell-level pairing, so it cannot tell
    # real per-cell signal from a day-level artefact. Stratifying preserves the trajectory and
    # destroys only the within-stratum pairing, which is what separates them.
    age_shuffle_strata: bool = False
    # CHANGE C-7: reject bulk samples that are not transcriptomes (G1 library band, G2 dynamic
    # range) and mask ΔAge for any cell_line the rejection leaves with no controls.
    # OFF by default because turning it on MOVES LABELS -- `N2_Fib_Sendai_Exp2` is donor N2's
    # entire zero-point, so N2's 21 ΔAge labels go with it, and HFF's labels move ~3x in five
    # of six folds because that column also sits inside `sigma_gill`. Enabling it is a
    # pre-registered change with its own bar and its own snapshot, never a side effect.
    # ONE flag, deliberately: the gate alone would strip the control and leave the ΔAge
    # unmasked, which is the state B2' exists to forbid.
    bulk_integrity_gate: bool = False


def _clock_range(clock: AgingClock, cfg: DataConfig) -> tuple[float, float] | None:
    """The clock's fitted age range, but only when the config opts in. Stage 1.5.3 C-2.

    Two conditions, both required: the run must ASK for the rule (`enforce_clock_age_range`),
    and the clock must actually declare a range. A clock with no provenance metadata returns
    None and the rule stays off -- an unknown range must never be treated as "no limits".
    """
    if not getattr(cfg, "enforce_clock_age_range", False):
        return None
    return getattr(clock, "age_range", None)


# --------------------------------------------------------------------------- #
# Per-chunk pipeline                                                          #
# --------------------------------------------------------------------------- #
def process_chunk(src, chunk, panel, clock: AgingClock, cfg: DataConfig, harmonizer=None,
                  census: dict | None = None,
                  no_control_lines: frozenset[str] = frozenset()):
    """Run the full transform for one chunk.

    Returns ``(samples, aux)`` where ``samples`` carry the *raw* control-relative
    ΔAge (no cell-cycle deconfounding applied here). ``aux`` is the per-cell
    material needed to re-apply a **train-fit** deconfounder later, or ``None``
    when deconfounding is off or the chunk has no age-valid cells. Fitting the
    deconfounder on train cells only requires the split assignment, which does
    not exist until every chunk has been read -- hence the deferred two-pass.

    ``census`` (optional dict) is filled with the per-line baseline count and
    composition -- Stage 1.5.2's gate G-a. Recording only: ΔAge is bit-identical
    whether or not it is supplied.
    """
    raw = apply_qc(src.fetch(chunk), cfg.qc)
    if len(raw.obs) == 0:
        return [], None

    norm = normalize_counts(raw.counts)
    sig = signature_scores(norm, raw.genes)

    if harmonizer is not None:
        # Cross-modality harmonization: Z-score against this dataset's control stats
        # for the model input + fate labels; project into the reference (bulk) scale
        # for the frozen clock (Gill Projection). Everything downstream selects its
        # own genes from the harmonized matrix on the harmonizer's gene space.
        ds_id = str(raw.obs["dataset_id"].iloc[0])
        x_scaled = harmonizer.transform(norm, raw.genes, ds_id).astype(np.float32)
        x_clock = harmonizer.project_to_clock(x_scaled).astype(np.float32)
        hgenes = harmonizer.genes
        y_cls = fate_labels(x_scaled, hgenes, raw.obs, cfg.label_tau)
        x_panel = to_panel_matrix(x_scaled, hgenes, panel)
        cc = cell_cycle_score(norm, raw.genes)     # cell cycle stays on raw norm
        d_age, age_mask, age_reason = delta_age(clock, x_clock, hgenes, raw.obs, raw.source,
                                                census=census,
                                                clock_age_range=_clock_range(clock, cfg),
                                                lines_without_controls=no_control_lines,
                                                enforce_no_unmasked_fallback=cfg.bulk_integrity_gate)
    else:
        y_cls = fate_labels(norm, raw.genes, raw.obs, cfg.label_tau)
        cc = cell_cycle_score(norm, raw.genes)
        x_panel = to_panel_matrix(norm, raw.genes, panel)
        # the clock consumes the FULL profile (its own gene panel), NOT the 2000-HVG
        # model input x_panel -- so aging genes filtered out of the HVG panel still
        # reach the clock. The model still trains on x_panel below.
        d_age, age_mask, age_reason = delta_age(clock, norm, raw.genes, raw.obs, raw.source,
                                                lines_without_controls=no_control_lines,
                                                enforce_no_unmasked_fallback=cfg.bulk_integrity_gate,
                                                census=census,
                                                clock_age_range=_clock_range(clock, cfg))
    cell_ids = raw.obs["cell_id"].tolist()
    # STAGE 1.5.3 C-I: the SAME policy minus the dataset rule. `masked_datasets=frozenset()`
    # rather than reading `age_mask_reason`, because the policy assigns only the FIRST reason that
    # fires -- a cell both out-of-clock-range and dataset-masked reads "dataset_policy", and
    # inverting on the reason string would wrongly readmit it. Recomputing is exact and cheap.
    deconfound_mask, _ = age_label_policy(
        len(cell_ids), raw.source, raw.obs,
        masked_datasets=frozenset(), clock_age_range=_clock_range(clock, cfg))
    # ARM C: the shuffle targets are exactly the cells arm B would MASK -- computable, but excluded
    # by the dataset policy. Derived through the same `age_label_policy` call rather than by
    # matching on cell_line, so arm C permutes precisely the label set arm B withholds.
    shuffle_mask = np.zeros(len(cell_ids), dtype=bool)
    if getattr(cfg, "age_shuffle_datasets", frozenset()):
        kept, _ = age_label_policy(
            len(cell_ids), raw.source, raw.obs,
            masked_datasets=frozenset(cfg.age_shuffle_datasets),
            clock_age_range=_clock_range(clock, cfg))
        shuffle_mask = deconfound_mask & ~kept
    aux: ChunkAux | None = None
    # Guard on the COMPUTABLE mask: an all-HFF chunk has no trainable label under arm B, but its
    # ΔAge values are exactly what the deconfounder needs. Guarding on `age_mask` dropped those
    # chunks from the fit entirely.
    if cfg.deconfound and deconfound_mask.any():
        aux = ChunkAux(
            cell_ids=cell_ids,
            cell_line=raw.obs["cell_line"].to_numpy().copy(),
            is_control=raw.obs["is_control"].to_numpy().astype(bool),
            d_age_raw=np.asarray(d_age, dtype=np.float64).copy(),
            cc=np.asarray(cc, dtype=np.float64).copy(),
            age_mask=age_mask.copy(),
            deconfound_mask=deconfound_mask.copy(),
            shuffle_mask=shuffle_mask.copy(),
            # ARM D's stratum key. `time_h` is the reprogramming timepoint axis
            # (`sources.py:292`); pairing it with `cell_line` keeps the key correct if a future
            # source contributes more than one line to the shuffled set.
            stratum=np.array(
                [f"{cl}|{th}" for cl, th in zip(raw.obs["cell_line"].to_numpy(),
                                                raw.obs["time_h"].to_numpy(), strict=True)],
                dtype=object),
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
        age_mask_reason=age_reason,
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
    apply_source_flags(cfg, sources)
    return sources


def apply_source_flags(cfg: DataConfig, sources) -> None:
    """Push config flags that GOVERN a source onto the sources actually being used.

    Change C-7. Called from BOTH `build_sources` and `run`, because `run` accepts injected
    sources (`run_multi_local.py` passes them, and so does every test with a synthetic source)
    and then never calls `build_sources`. A flag applied only on the construction path would
    silently never reach a production build -- a gate that cannot fire, which is precisely the
    defect class C-7 exists to remove. Idempotent, so applying it twice is harmless.

    Set centrally rather than per-source in `source_specs` so the gate cannot be enabled for one
    source and not another, and cannot be switched independently of rule 4 that consumes it.
    """
    for src in sources:
        if hasattr(src, "bulk_integrity_gate"):
            src.bulk_integrity_gate = bool(cfg.bulk_integrity_gate)


def lines_without_controls(cfg: DataConfig, sources) -> frozenset[str]:
    """Union, over sources that can answer, of cell_lines with ZERO admissible controls.

    Change C-7, component B. Only a source knows which of its own samples are controls, and
    only the bulk gate can remove one, so this asks the sources instead of re-reading the
    corpus -- there is no second full pass. A source that does not implement
    ``lines_without_controls`` contributes nothing, so an un-censused line is never masked by
    accident: the set names lines KNOWN to have none, and silence is not evidence of absence.

    Empty when the gate is off, which is what keeps B4 (bit-identical when disabled) true.
    """
    if not cfg.bulk_integrity_gate:
        return frozenset()
    out: set[str] = set()
    for s in sources:
        fn = getattr(s, "lines_without_controls", None)
        if callable(fn):
            out |= set(fn())
    return frozenset(out)


# --------------------------------------------------------------------------- #
# Orchestrator                                                                #
# --------------------------------------------------------------------------- #
def fit_harmonizer(cfg: DataConfig, work):
    """Pre-pass: fit the cross-modality Harmonizer on TRAINING control cells only.

    'Training control' = a control cell whose cell_line is NOT held out (the
    holdout regime holds out whole donors, so this is decidable from cell_line
    alone, before the full split). Statistics are grouped by ``dataset_id`` and the
    held-out donor never contributes -- enforced by a hard assertion.
    """
    from .harmonize import Harmonizer

    heldout = set(cfg.holdout_cell_lines)
    controls: dict[str, list] = {}
    leaked = 0
    for src, chunk in work:
        raw = apply_qc(src.fetch(chunk), cfg.qc)
        if len(raw.obs) == 0:
            continue
        obs = raw.obs
        if "dataset_id" not in obs.columns:
            raise ValueError("harmonize=True requires sources that stamp dataset_id "
                             "(GSE242423SingleCellSource / GillReprogrammingSource)")
        is_ctrl = obs["is_control"].to_numpy().astype(bool)
        not_test = ~obs["cell_line"].isin(heldout).to_numpy()
        keep = is_ctrl & not_test
        leaked += int((is_ctrl & obs["cell_line"].isin(heldout).to_numpy()).sum())
        if not keep.any():
            continue
        norm = normalize_counts(raw.counts)
        for ds in np.unique(obs["dataset_id"].to_numpy()[keep]):
            m = keep & (obs["dataset_id"].to_numpy() == ds)
            controls.setdefault(str(ds), []).append((norm[m], raw.genes))
    if not controls:
        raise ValueError("harmonizer pre-pass found no training control cells")
    h = Harmonizer.fit(controls, ref_dataset=cfg.harmonize_ref_dataset)
    # Hard leakage assertion: not one held-out control cell entered the statistics.
    # (We simply never added them above; this documents and guarantees it.)
    log_event(log, "harmonizer.fit", datasets=sorted(controls), n_genes=len(h.genes),
              heldout=sorted(heldout), heldout_controls_excluded=leaked)
    return h


def run(cfg: DataConfig, sources: list[DataSource] | None = None,
        clock: AgingClock | None = None) -> dict:
    """Execute the ETL and return a summary dict."""
    paths = ArtifactPaths.of(cfg.out)
    sources = sources if sources is not None else build_sources(cfg)
    # CHANGE C-7. Applied HERE, not only in `build_sources`, because callers may INJECT sources
    # -- `run_multi_local.py` does, and so does every test that hands in a synthetic source. A
    # flag set only on the construction path would silently never reach a production build,
    # which is exactly the class of defect this gate exists to remove. Setting it on whatever
    # sources `run` actually uses makes the config's flag mean the same thing either way.
    apply_source_flags(cfg, sources)
    if not sources:
        raise ValueError("no data sources configured")
    work = plan_all(sources)

    panel = load_or_fit_panel(cfg, work)
    clock = clock if clock is not None else build_clock(cfg, panel)
    # CHANGE C-7, component B. Computed BEFORE the harmonizer and independently of it:
    # `fit_harmonizer` is `if cfg.harmonize`, so a predicate hosted there would silently not
    # exist in harmonize=False builds (the arm B/C/D probes, any single-dataset build) -- a
    # data-integrity invariant that evaporates when a flag is off is a guard that cannot fire.
    no_control_lines = lines_without_controls(cfg, sources)
    if no_control_lines:
        log_event(log, "c7.no_control_lines", lines=sorted(no_control_lines))
    harmonizer = fit_harmonizer(cfg, work) if cfg.harmonize else None
    if harmonizer is not None:
        harmonizer.to_json(paths.bundle_dir / "harmonization.json"
                           if paths.bundle_dir.exists() else f"{cfg.out}/harmonization.json")
    tracker = ProgressTracker(paths.progress_file)

    label_counts: Counter[str] = Counter()
    n_age_labeled = 0
    baseline_census: dict = {}      # gate G-a: what each ΔAge zero-point actually rests on

    for src, chunk in work:
        cid = chunk["id"]
        if tracker.is_done(cid):
            continue
        sid = io.sanitize_id(cid)
        try:
            chunk_census: dict = {}
            samples, aux = process_chunk(src, chunk, panel, clock, cfg, harmonizer,
                                         census=chunk_census,
                                         no_control_lines=no_control_lines)
            # Key by chunk AND line, never by line alone. `cell_line` is NOT unique across
            # chunks -- HFF spans 45 of them (verify_stage1_5_results.json) -- so a plain
            # `.update()` keyed on the line silently kept 1 record and discarded 44, for the
            # dataset that carries ~99.8% of the age labels. The whole point of G-a is that a
            # baseline problem in ANY chunk stays visible.
            for line, rec in chunk_census.items():
                baseline_census[f"{cid}::{line}"] = {**rec, "chunk_id": cid, "cell_line": line}
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
            log_event(log, "chunk.done", chunk=cid, n=len(samples), baseline=chunk_census)
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
        # Gate G-a: what every ΔAge zero-point rests on, and what is wrong with it. Persisted so
        # a run can be audited without re-reading the raw data -- Stage 1.5 made `n=0` visible;
        # this makes `n=1` and cross-batch baselines visible too.
        "baseline_census": baseline_census,
        "baseline_warnings": census_warnings(baseline_census),
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
        deconfound_mask=np.asarray(aux.deconfound_mask, dtype=bool),
        shuffle_mask=np.asarray(aux.shuffle_mask, dtype=bool),
        stratum=np.asarray(aux.stratum).astype("U"),
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
            # Tolerant read: a sidecar written before C-I has no `deconfound_mask`. Falling back
            # to `age_mask` reproduces the pre-C-I behaviour exactly rather than crashing a
            # RESUMED build -- the same lesson as `rewrite_shard_yage`'s backfill.
            deconfound_mask=(z["deconfound_mask"] if "deconfound_mask" in z.files
                             else z["age_mask"]),
            shuffle_mask=(z["shuffle_mask"] if "shuffle_mask" in z.files
                          else np.zeros(len(z["age_mask"]), dtype=bool)),
            stratum=(z["stratum"].astype(object) if "stratum" in z.files
                     else np.full(len(z["age_mask"]), "__all__", dtype=object)),
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
    # STAGE 1.5.3 C-I: fit on cells whose ΔAge is COMPUTABLE, not on cells whose label the age
    # head may train on. Those are different questions, and conflating them made `y_age` depend on
    # `AGE_MASKED_DATASETS` -- so step 6's two arms had different TARGET VARIABLES, not just
    # different label counts. See results/STEP6_REPORT.md section 3, confound C-I.
    for aux in aux_by_sid.values():
        for i, cell in enumerate(aux.cell_ids):
            if aux.deconfound_mask[i] and cell in train_ids:
                d_tr.append(aux.d_age_raw[i])
                cc_tr.append(aux.cc[i])
    coef = (fit_deconfounder(np.asarray(d_tr), np.asarray(cc_tr))
            if len(d_tr) >= 2 else (0.0, 0.0))

    # Pass 2: re-apply the single train-fit transform to every shard. `deconfound_mask` again --
    # the control RE-CENTRING is part of computing ΔAge, so it must not narrow to the training
    # subset either. Which cells the loss uses is governed by the `age_mask` COLUMN, written
    # separately by `assemble_samples`; this function only decides the VALUE of `y_age`.
    ys: dict[str, np.ndarray] = {}
    for sid, aux in aux_by_sid.items():
        d = deconfound_age(aux.d_age_raw, aux.cc, coef)
        m = aux.deconfound_mask
        y = np.full(d.shape[0], np.nan, dtype=np.float64)
        if m.any():
            y[m] = recenter_on_control_arrays(d[m], aux.cell_line[m], aux.is_control[m])
        ys[sid] = y

    ys = _shuffle_age_labels(cfg, aux_by_sid, ys)

    for sid, y in ys.items():
        io.rewrite_shard_yage(paths.shard_file(sid), y)
    return coef


def _shuffle_age_labels(cfg, aux_by_sid: dict, ys: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """ARM C — permute ΔAge labels among the targeted cells, GLOBALLY across shards. Pure-ish.

    The label-permutation control for step 6's follow-up. Arm A trains on 33,688 labels, arm B on
    75; the ranking gap between them is confounded between *"HFF's labels carry information"* and
    *"75 labels is simply too few to learn from"*. Arm C holds label VOLUME at arm A's level and
    destroys only the cell<->label PAIRING, so the two explanations separate:

        C ranks like A -> the gain was volume / trunk regularisation; the labels are uninformative
        C ranks like B -> the labels carry real signal despite the artefact

    Three properties this must have, and why each is deliberate:

      * It runs **after** the deconfounder fit and the control re-centring, so the fitted
        coefficient and every ΔAge VALUE are bit-identical to arm A. Only the assignment moves.
      * It permutes **globally across shards**, not within each chunk. A within-chunk shuffle would
        leave the between-chunk structure intact -- chunks are timepoint-homogeneous, so per-chunk
        mean ΔAge would survive and the control would be far weaker than it looks.
      * It is **deterministic** given `age_shuffle_seed`, and the shard iteration order is sorted so
        the permutation does not depend on dict ordering.

    A no-op unless `cfg.age_shuffle_datasets` is set, so arms A and B are untouched by its presence.
    """
    if not getattr(cfg, "age_shuffle_datasets", frozenset()):
        return ys
    stratified = bool(getattr(cfg, "age_shuffle_strata", False))
    # Group the target slots by stratum. Unstratified (arm C) is the same code with every slot in
    # one group, so the two arms cannot drift apart through duplicated logic.
    groups: dict[str, list[tuple[str, int]]] = {}
    for sid in sorted(aux_by_sid):                    # sorted: order must not depend on dict order
        aux = aux_by_sid[sid]
        y = ys[sid]
        for i in np.flatnonzero(aux.shuffle_mask):
            if not np.isnan(y[i]):                    # only real labels take part
                key = str(aux.stratum[i]) if stratified else "__all__"
                groups.setdefault(key, []).append((sid, int(i)))

    seed = int(getattr(cfg, "age_shuffle_seed", 0))
    total = moved = singleton = 0
    # One generator, consumed in sorted-key order, so the whole permutation is reproducible from
    # the seed alone and does not depend on dict insertion order.
    rng = np.random.default_rng(seed)
    for key in sorted(groups):
        slots = groups[key]
        total += len(slots)
        if len(slots) < 2:
            singleton += len(slots)     # nothing to permute against; label stays put
            continue
        vals = np.array([ys[sid][i] for sid, i in slots], dtype=np.float64)
        perm = rng.permutation(len(vals))
        for (sid, i), v in zip(slots, vals[perm], strict=True):
            ys[sid][i] = v
        moved += int((perm != np.arange(len(perm))).sum())
    if total < 2:
        log.warning("age label shuffle requested but only %d target label(s) found; no-op", total)
        return ys
    log_event(log, "age_labels.shuffled", n=total, moved=moved, seed=seed,
              stratified=stratified, n_strata=len(groups), singleton_strata=singleton)
    return ys


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
