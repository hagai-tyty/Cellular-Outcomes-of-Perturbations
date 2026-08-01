"""STAGE 1.5 §9 — is the Fleischer clock BROKEN on this data, or just MIS-APPLIED / OUT-OF-DOMAIN?

    python experiments/diag_clock_validity.py                       # defaults: D:\\Gill
    python experiments/diag_clock_validity.py "D:\\Gill" "D:\\GSE113957"

READ-ONLY. Writes `diag_clock_validity_results.json`. Nothing is rebuilt or refitted; `src/` is not
touched. This runs BEFORE any of the four fix options (replace clock / replace target / narrow /
more data), because all four assume the instrument is fundamentally broken — and that has not been
established. Getting THIS wrong is the expensive mistake: abandoning a working target, or shipping
a broken one. So every axis is measured independently, with a pre-registered interpretation.

WHAT M1/E1/E1b ESTABLISHED, AND THE HOLE IN IT (STAGE_1_5 §7-§8)
---------------------------------------------------------------
M1: the clock did not separate the age extremes (contrast 11.8 yr vs a 53 yr gap). E1b: predicted
age rose during reprogramming (p≈0.045, marginal). Conclusion drawn: "ΔAge target unvalidated."

But three confounds were never ruled out, and each can produce those exact failures WITHOUT the
clock being wrong about aging:

  (H1 APPLICATION)  predict_age sums `w_g * x_g` only over genes present in the data
                    (`weights.get(g, 0.0)`). The clock has 33,155 genes. If the Gill matrix is
                    missing a large fraction, most of the model is silently dropped and every
                    prediction collapses toward the intercept (72.4). The M1 ages (36-99) DO
                    cluster around 72.4 -- the fingerprint of a partial clock.

  (H2 OUT-OF-RANGE) the clock was fit on ages [1, 96]. M1's "young" anchor was the two NEONATAL
                    donors (age 0), below that range; N2 read 98.7. Among IN-RANGE adults
                    (Y1,Y2 ~32 vs O1,O2 =53) the day-0 contrast is ~18 yr for a ~21 yr true gap --
                    i.e. the clock DOES track in-domain adult fibroblast age. M1 may have failed by
                    anchoring on donors the clock cannot read.

  (H3 OUT-OF-DOMAIN) the clock was fit on fibroblasts. Days 0-15 of OSKM reprogramming are cells
                    LEAVING fibroblast identity. If the "age rises" signal is driven by
                    OSKM/cell-cycle/pluripotency genes the clock happens to weight, that is the
                    clock reading cell-STATE, not aging -- fixable by domain restriction or a
                    reprogramming-aware target, not proof the aging axis is unreadable.

Each hypothesis has a decisive check below. The point is not to make the clock pass -- it is to
locate the failure precisely enough that the RIGHT fix is obvious and the WRONG retreat is avoided.

DISCIPLINE. Pure verdict functions, separated from I/O and unit-tested on every branch. Real-data
wiring reuses the production path (GillReprogrammingSource -> normalize_counts -> LinearClock) so
what is measured is what the pipeline actually does. Nothing is tuned; there are no free knobs.
"""

from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)


# --------------------------------------------------------------------------- #
# Pure logic — data-free, fully unit-tested; nothing below imports repo data.  #
# --------------------------------------------------------------------------- #
DONOR_AGE: dict[str, float] = {"N2": 0.0, "N3": 0.0, "Y1": 29.0, "Y2": 35.0, "O1": 53.0, "O2": 53.0}
CLOCK_CV_MAE = 12.26879346460328
CLOCK_AGE_RANGE = (1.0, 96.0)

# Coverage bars (fraction of the clock's total |weight| present in the data). These are about the
# MODEL, not gene counts: a clock can miss half its genes yet keep most of its predictive mass, or
# lose it. Bars pre-registered here.
COVERAGE_OK = 0.90            # >= this: application is not the problem
COVERAGE_CRIPPLED = 0.70      # <  this: the clock is running on a fraction of itself -> H1

# Marker sets for the directional attribution (H3). HGNC symbols. Curated, not tuned on results.
OSKM_PLURIPOTENCY = frozenset({
    "POU5F1", "SOX2", "KLF4", "MYC", "MYCL", "NANOG", "LIN28A", "LIN28B", "DPPA4", "DPPA2",
    "ZFP42", "SALL4", "TDGF1", "DNMT3B", "TERT", "UTF1", "ESRG", "PRDM14"})
CELL_CYCLE = frozenset({
    "MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "CCNB2", "CCNA2", "CCNE2", "CENPF", "BIRC5",
    "AURKA", "AURKB", "CDC20", "UBE2C", "BUB1", "PLK1", "TYMS", "RRM2", "MCM2", "MCM6"})
SENESCENCE_AGING = frozenset({           # genes a REAL aging signal would move; a positive control
    "CDKN2A", "CDKN1A", "GLB1", "SERPINE1", "LMNB1", "IGFBP3", "IGFBP5", "TP53", "B2M", "IL6"})


def weighted_coverage(clock_weights: dict[str, float], data_genes: list[str]) -> dict:
    """How much of the clock actually gets EVALUATED on this data.

    `predict_age` sums only over genes present in `data_genes`, so the operative quantity is not the
    gene COUNT but the fraction of total |weight| that survives. Also reports how concentrated the
    weight is (top-k share), because a clock whose signal lives in a few genes fails differently
    from one spread thin.
    """
    present = set(data_genes)
    items = [(g, abs(w)) for g, w in clock_weights.items()]
    tot = sum(w for _g, w in items) or 1.0
    kept = sum(w for g, w in items if g in present)
    n_overlap = sum(1 for g, _w in items if g in present)
    absw = np.array(sorted((w for _g, w in items), reverse=True))
    cum = np.cumsum(absw) / tot
    # smallest k carrying 50% / 90% of the weight -> how concentrated the clock is
    k50 = int(np.searchsorted(cum, 0.50) + 1)
    k90 = int(np.searchsorted(cum, 0.90) + 1)
    return {
        "n_clock_genes": len(clock_weights),
        "n_data_genes": len(present),
        "n_overlap": int(n_overlap),
        "frac_genes_present": float(n_overlap / max(len(clock_weights), 1)),
        "frac_abs_weight_present": float(kept / tot),
        "weight_k50": k50, "weight_k90": k90,
    }


def coverage_verdict(frac_abs_weight: float) -> dict:
    """H1. Is the clock running on enough of itself to be trusted?"""
    if frac_abs_weight >= COVERAGE_OK:
        s, r = "OK", f"{frac_abs_weight:.1%} of the clock's weight is present; application is not the problem"
    elif frac_abs_weight < COVERAGE_CRIPPLED:
        s, r = "CRIPPLED", (f"only {frac_abs_weight:.1%} of the clock's weight survives gene matching; "
                            "the clock is running on a fraction of itself -> APPLICATION defect (H1), "
                            "recoverable by fixing gene mapping, NOT by replacing the clock")
    else:
        s, r = "DEGRADED", (f"{frac_abs_weight:.1%} of the clock's weight is present -- materially "
                            "reduced; predictions are partially intercept-driven and must be read with care")
    return {"status": s, "frac_abs_weight_present": float(frac_abs_weight), "reason": r}


def intercept_dominance(pred_ages: list[float], intercept: float) -> dict:
    """H1 corollary. If predictions barely move off the intercept, the clock is effectively dead on
    this data regardless of why."""
    a = np.asarray([x for x in pred_ages if np.isfinite(x)], float)
    if len(a) < 2:
        return {"status": "CANNOT_VERIFY", "n": int(len(a)), "reason": "need >=2 predictions"}
    spread = float(a.std(ddof=1))
    off = float(np.mean(np.abs(a - intercept)))
    # a healthy clock over a real age span moves many years; near-intercept clustering is the tell
    status = "DEAD_NEAR_INTERCEPT" if spread < 3.0 else "MOVES"
    return {"status": status, "pred_sd_years": spread, "mean_abs_offset_from_intercept": off,
            "intercept": float(intercept),
            "reason": f"predictions have SD {spread:.1f} yr around a {intercept:.1f} yr intercept"}


def in_range_age_tracking(pred: dict[str, float], chrono: dict[str, float],
                          age_range: tuple[float, float] = CLOCK_AGE_RANGE) -> dict:
    """H2. Does the clock track age among donors WITHIN its fitted range?

    Splits donors by whether chrono age is in `age_range`. Reports the young/old contrast among
    in-range donors and the Spearman over them. This is the test M1 should arguably have run: a
    clock cannot be blamed for misreading ages it was never fit on.
    """
    in_range = {d: p for d, p in pred.items()
                if d in chrono and age_range[0] <= chrono[d] <= age_range[1]}
    out_range = {d: p for d, p in pred.items()
                 if d in chrono and not (age_range[0] <= chrono[d] <= age_range[1])}
    if len(in_range) < 3:
        return {"status": "CANNOT_VERIFY", "n_in_range": len(in_range),
                "out_of_range_donors": {d: chrono.get(d) for d in out_range},
                "reason": f"only {len(in_range)} donors inside {age_range}; cannot test in-range tracking"}
    ages = np.array([chrono[d] for d in in_range])
    vals = np.array([in_range[d] for d in in_range])
    # young/old GROUPS split at the median age, not single min/max donors -- at n=4 a group split
    # uses every in-range donor and is far less swingy than picking two extremes.
    med = float(np.median(ages))
    young = vals[ages < med]
    old = vals[ages > med]
    if len(young) == 0 or len(old) == 0:          # median tie left one side empty -> fall back
        young, old = vals[ages == ages.min()], vals[ages == ages.max()]
    contrast = float(old.mean() - young.mean())
    true_gap = float(ages[ages > med].mean() - ages[ages < med].mean()) if (
        (ages > med).any() and (ages < med).any()) else float(ages.max() - ages.min())
    # rank tracking within range
    ro = np.argsort(np.argsort(vals))
    rt = np.argsort(np.argsort(ages))
    rho = float(np.corrcoef(ro, rt)[0, 1]) if len(in_range) > 2 else float("nan")
    tracks = contrast > 0 and rho > 0
    return {
        "status": "TRACKS_IN_RANGE" if tracks else "NO_IN_RANGE_TRACKING",
        "n_in_range": len(in_range), "in_range_contrast_years": contrast,
        "true_gap_years": true_gap, "spearman_in_range": rho,
        "out_of_range_donors": {d: chrono.get(d) for d in out_range},
        "out_of_range_predictions": {d: float(out_range[d]) for d in out_range},
        "reason": (f"in-range contrast {contrast:+.1f} yr for a {true_gap:.0f} yr true gap, "
                   f"Spearman {rho:+.2f} over {len(in_range)} in-range donors; "
                   f"{len(out_range)} donors are OUT of range {age_range}"),
    }


def denominator_sensitivity(age_full: list[float], age_restricted: list[float]) -> dict:
    """H1 corollary. CP10k divides by the library size over whatever genes are in the matrix. The
    clock was fit with CP10k over ITS gene set. Comparing predictions normalised over the full data
    gene set vs over the clock-overlap set bounds how much the denominator mismatch distorts."""
    a = np.asarray(age_full, float)
    b = np.asarray(age_restricted, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 2:
        return {"status": "CANNOT_VERIFY", "reason": "need >=2 paired predictions"}
    shift = float(np.mean(np.abs(a - b)))
    status = "SENSITIVE" if shift > CLOCK_CV_MAE / 2 else "STABLE"
    return {"status": status, "mean_abs_shift_years": shift,
            "reason": f"switching the CP10k gene space moves predicted age by {shift:.1f} yr on average"}


def reproduction_verdict(mae: float | None, n: int) -> dict:
    """The gold check (H1 vs real failure). Apply the clock through THIS pipeline to fibroblasts of
    KNOWN age (its own training set if available). If MAE is near the clock's own CV error, the
    pipeline application is clean and any Gill failure is genuinely domain/data. If MAE blows up,
    the pipeline corrupts the clock and the whole escalation is partly an artefact."""
    if mae is None:
        return {"status": "SKIPPED", "reason": "no known-age fibroblast set available (pass its dir)"}
    if mae <= CLOCK_CV_MAE * 1.5:
        s, r = "REPRODUCES", (f"MAE {mae:.1f} yr on n={n} known-age fibroblasts, within 1.5x the "
                              f"clock's own CV error ({CLOCK_CV_MAE:.1f}) -> pipeline application is clean")
    elif mae <= CLOCK_CV_MAE * 3:
        s, r = "DEGRADED", (f"MAE {mae:.1f} yr -- worse than the clock's CV error but not random; "
                            "application is partly degraded")
    else:
        s, r = "BROKEN", (f"MAE {mae:.1f} yr on its OWN domain -> the pipeline application corrupts the "
                          "clock; M1/E1 inherit that and the escalation is (partly) an artefact")
    return {"status": s, "mae_years": float(mae), "n": int(n), "reason": r}


def attribute_direction(contrib_by_gene: dict[str, float]) -> dict:
    """H3. During reprogramming, decompose the age CHANGE into gene contributions and ask which
    biology drives the (positive, 'older') part. contrib_by_gene[g] = w_g * (x_g[late]-x_g[early]),
    so it sums to the total predicted age change. If OSKM/cell-cycle genes dominate the positive
    contribution, the clock is reading cell-STATE, not aging (out-of-domain)."""
    if not contrib_by_gene:
        return {"status": "CANNOT_VERIFY", "reason": "no gene contributions supplied"}
    total = float(sum(contrib_by_gene.values()))
    pos = {g: c for g, c in contrib_by_gene.items() if c > 0}
    pos_tot = float(sum(pos.values())) or 1.0

    def share(genes):
        return float(sum(c for g, c in pos.items() if g in genes) / pos_tot)

    s_oskm, s_cc, s_sen = share(OSKM_PLURIPOTENCY), share(CELL_CYCLE), share(SENESCENCE_AGING)
    confound = s_oskm + s_cc
    top = sorted(contrib_by_gene.items(), key=lambda kv: kv[1], reverse=True)[:15]
    if confound >= 0.30:
        s, r = "OUT_OF_DOMAIN_CONFOUND", (
            f"{confound:.0%} of the 'age rises' signal comes from OSKM/pluripotency ({s_oskm:.0%}) + "
            f"cell-cycle ({s_cc:.0%}) genes -- the clock is reading reprogramming cell-STATE, not "
            "aging. A domain-restricted or reprogramming-aware target addresses this (H3).")
    elif s_sen >= 0.20:
        s, r = "AGING_GENES_DRIVE_IT", (
            f"senescence/aging genes carry {s_sen:.0%} of the rise -- the wrong-direction signal is "
            "in genuine aging genes, which is harder to explain away as domain shift")
    else:
        s, r = "DIFFUSE", ("the rise is spread across genes with no dominant category; inconclusive "
                           "on H3, weight the other checks")
    return {"status": s, "share_oskm_pluripotency": s_oskm, "share_cell_cycle": s_cc,
            "share_senescence_aging": s_sen, "confound_share": confound,
            "total_change_years": total, "top_contributors": [(g, float(c)) for g, c in top],
            "reason": r}


def decide(coverage: dict, in_range: dict, reproduction: dict, attribution: dict) -> dict:
    """Fold the checks into the one conclusion that picks the fix. Ordered by decisiveness: an
    application defect (H1) makes everything else moot and is the recoverable, best-case answer."""
    if coverage["status"] == "CRIPPLED" or reproduction["status"] == "BROKEN":
        return {"action": "FIX_APPLICATION",
                "reason": "the clock is mis-applied on this data (H1): gene coverage and/or its own-"
                          "domain reproduction fail. Fix gene mapping / normalisation and RE-RUN "
                          "M1/E1 before any talk of replacing the clock. ΔAge is likely recoverable "
                          "as-is. This is the best case and it is NOT what the escalation assumed."}
    if in_range["status"] == "TRACKS_IN_RANGE" and attribution["status"] == "OUT_OF_DOMAIN_CONFOUND":
        return {"action": "TARGET_RECOVERABLE_DOMAIN_FIX",
                "reason": "the clock DOES track age on in-range fibroblasts, and the reprogramming "
                          "'wrong direction' is out-of-domain cell-state (H2+H3). ΔAge stays; the fix "
                          "is to restrict the clock to its domain or move to a reprogramming-aware "
                          "rejuvenation target (option B). NOT a retreat."}
    if in_range["status"] == "TRACKS_IN_RANGE":
        return {"action": "IN_DOMAIN_OK_INVESTIGATE_REPROGRAMMING",
                "reason": "the clock tracks in-range adult fibroblast age (H2 confirmed), so 'the "
                          "clock can't read age' is too strong. The reprogramming-phase behaviour "
                          "needs the attribution check to finish before choosing between a domain fix "
                          "and a new target."}
    if coverage["status"] == "OK" and reproduction["status"] == "REPRODUCES":
        return {"action": "GENUINE_CLOCK_LIMITATION",
                "reason": "application is clean (coverage OK, reproduces on its own domain) yet the "
                          "clock does not track age even in range -> this is a real clock/data "
                          "limitation. Options A (reprogramming-validated clock) or B (validated "
                          "signature) are now justified; more data (D) may be required to choose."}
    return {"action": "INCONCLUSIVE",
            "reason": "the checks do not yet converge; do not pick a fix. Most likely the gold "
                      "reproduction check was skipped (no known-age set) -- supply one and re-run."}


def bars() -> list[dict]:
    """Pre-registered thresholds and what each decides (ground rule §5b)."""
    return [
        {"id": "H1-coverage", "bar": f"frac of clock |weight| present >= {COVERAGE_OK:.0%} = OK, "
         f"< {COVERAGE_CRIPPLED:.0%} = CRIPPLED", "decides": "is the clock even fully applied"},
        {"id": "H1-reproduction", "bar": f"MAE on known-age fibroblasts <= {CLOCK_CV_MAE*1.5:.0f} yr "
         "= REPRODUCES", "decides": "does the pipeline apply the clock correctly on its own domain"},
        {"id": "H2-in-range", "bar": "positive young->old contrast AND Spearman > 0 among in-range "
         "donors", "decides": "does the clock track age it was actually fit to read"},
        {"id": "H3-attribution", "bar": "OSKM+cell-cycle share of the 'age rises' signal >= 30% "
         "= out-of-domain confound", "decides": "is the reprogramming reversal cell-state, not aging"},
    ]


# --------------------------------------------------------------------------- #
# Real-data wiring (imports repo machinery only when actually run)            #
# --------------------------------------------------------------------------- #
def _load_gill(gill_dir: str):
    root = Path(__file__).resolve().parents[1]
    for p in (root, root / "local_runners", root / "src"):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))
    from run_multi_local import discover_gill  # type: ignore

    from cellfate.data.aging import LinearClock
    from cellfate.data.sources import GillReprogrammingSource

    expr, series = discover_gill(gill_dir)
    src = GillReprogrammingSource(expr_tsv=expr, series_matrix=series)
    src._load()
    clock = LinearClock.from_json(root / "configs" / "clocks" / "fleischer_clock.json")
    return src, clock, root


def _predict(clock, mat_linear_rpm: np.ndarray, genes: list[str]) -> np.ndarray:
    from cellfate.data.normalize import normalize_counts
    return clock.predict_age(normalize_counts(mat_linear_rpm), genes)


def run_gill_checks(gill_dir: str) -> dict:
    """H1 coverage, H2 in-range tracking, denominator sensitivity, intercept dominance, and H3
    directional attribution — all from the Gill data + the frozen clock."""
    src, clock, _root = _load_gill(gill_dir)
    genes = list(src._genes)

    # ---- H1: gene coverage ----
    cov = weighted_coverage(clock.weights, genes)
    cov["verdict"] = coverage_verdict(cov["frac_abs_weight_present"])

    # ---- day-0 baseline age per donor (production path) ----
    from cellfate.data.normalize import normalize_counts  # noqa: F401 (kept explicit for clarity)
    baseline_pred: dict[str, float] = {}
    baseline_pred_restricted: dict[str, float] = {}
    clock_genes = set(clock.weights)
    keep_idx = [i for i, g in enumerate(genes) if g in clock_genes]
    genes_restricted = [genes[i] for i in keep_idx]
    for chunk in src.plan():
        donor = chunk["cell_line"]
        cols = [c for c in src._rpm.columns
                if c in src._meta and src._meta[c]["donor"] == donor
                and str(src._meta[c].get("ctype", "")).strip().lower() != "ipsc"]
        base_cols = [c for c in cols if src._meta[c]["day"] == 0.0]
        if not base_cols:
            continue
        mat = src._rpm[base_cols].to_numpy(dtype=np.float64).T
        baseline_pred[donor] = float(np.mean(_predict(clock, mat, genes)))
        # denominator sensitivity: CP10k over the clock-overlap gene space only
        mat_r = mat[:, keep_idx]
        baseline_pred_restricted[donor] = float(np.mean(_predict(clock, mat_r, genes_restricted)))

    # ---- H2: in-range tracking + out-of-range flagging ----
    in_range = in_range_age_tracking(baseline_pred, DONOR_AGE)
    intercept = float(clock.intercept)
    dom = intercept_dominance(list(baseline_pred.values()), intercept)
    den = denominator_sensitivity([baseline_pred[d] for d in baseline_pred],
                                  [baseline_pred_restricted[d] for d in baseline_pred])

    # ---- H3: directional attribution over the reprogramming phase (day 0 -> ~15) ----
    attribution = _attribute_reprogramming(src, clock, genes, day_max=15.0)

    return {"coverage": cov, "baseline_pred": baseline_pred,
            "baseline_pred_restricted": baseline_pred_restricted,
            "in_range": in_range, "intercept_dominance": dom,
            "denominator_sensitivity": den, "attribution": attribution}


def _attribute_reprogramming(src, clock, genes: list[str], day_max: float) -> dict:
    """Per-gene contribution to the age CHANGE from day 0 to the last reprogramming-phase day,
    pooled across donors. contrib[g] = w_g * mean_donor(x_g[late] - x_g[day0])."""
    from cellfate.data.normalize import normalize_counts

    w = np.array([clock.weights.get(g, 0.0) for g in genes], dtype=np.float64)
    deltas = []
    for chunk in src.plan():
        donor = chunk["cell_line"]
        cols = [c for c in src._rpm.columns
                if c in src._meta and src._meta[c]["donor"] == donor
                and str(src._meta[c].get("ctype", "")).strip().lower() != "ipsc"]
        base = [c for c in cols if src._meta[c]["day"] == 0.0]
        late = [c for c in cols if 0.0 < src._meta[c]["day"] <= day_max]
        if not base or not late:
            continue
        xb = normalize_counts(src._rpm[base].to_numpy(dtype=np.float64).T).mean(axis=0)
        xl = normalize_counts(src._rpm[late].to_numpy(dtype=np.float64).T).mean(axis=0)
        deltas.append(xl - xb)
    if not deltas:
        return {"status": "CANNOT_VERIFY", "reason": "no donor had both day-0 and reprogramming-phase samples"}
    mean_delta = np.mean(deltas, axis=0)
    contrib = {g: float(w[i] * mean_delta[i]) for i, g in enumerate(genes) if w[i] != 0.0}
    return attribute_direction(contrib)


def run_reproduction(known_dir: str | None) -> dict:
    """Gold check: apply the clock through this pipeline to fibroblasts of KNOWN age.

    Best-effort loader for a GEO series-matrix + expression pair (the clock's own GSE113957, or any
    fibroblast-aging set). Skips cleanly if absent or unparseable — the Gill checks stand alone."""
    if not known_dir or not Path(known_dir).exists():
        return reproduction_verdict(None, 0)
    try:

        from cellfate.data.aging import LinearClock
        from cellfate.data.normalize import normalize_counts
        root = Path(__file__).resolve().parents[1]
        clock = LinearClock.from_json(root / "configs" / "clocks" / "fleischer_clock.json")
        ages, mat, genes = _load_known_age_fibroblasts(known_dir)
        if ages is None or len(ages) < 5:
            return {"status": "SKIPPED",
                    "reason": f"could not parse >=5 known-age samples from {known_dir}; format not recognised"}
        pred = np.asarray(clock.predict_age(normalize_counts(mat), genes), float)
        ages = np.asarray(ages, float)
        mae = float(np.mean(np.abs(pred - ages)))
        v = reproduction_verdict(mae, len(ages))
        v["age_range_tested"] = [float(np.min(ages)), float(np.max(ages))]
        # Secondary, robust to absolute-scale drift (annotation-release/gene-set differences between
        # this NCBI matrix and what fit_clock used): does predicted age RANK the samples by true age?
        from scipy.stats import spearmanr
        v["spearman_pred_vs_age"] = float(spearmanr(pred, ages).correlation)
        v["pearson_pred_vs_age"] = float(np.corrcoef(pred, ages)[0, 1])
        # coverage of the clock on THIS matrix -- if low, a poor MAE is application, not the clock
        cov = weighted_coverage(clock.weights, genes)
        v["weighted_coverage_here"] = cov["frac_abs_weight_present"]
        v["reason"] += (f" | Spearman(pred,age)={v['spearman_pred_vs_age']:+.2f}, "
                        f"weighted coverage {cov['frac_abs_weight_present']:.0%}")
        return v
    except Exception as exc:  # noqa: BLE001
        return {"status": "SKIPPED", "reason": f"reproduction check errored: {exc!r}"[:200]}


# --- pure parsing helpers for the NCBI-generated GSE113957 layout (unit-tested) --------------- #
def parse_age_value(cell: str) -> float:
    """`'age: 30'` / `'30'` / `'30 years'` -> 30.0; unparseable -> NaN. Pure."""
    import re
    if cell is None:
        return float("nan")
    tail = cell.split(":", 1)[1] if ":" in cell else cell
    m = re.search(r"(\d+(?:\.\d+)?)", tail)
    return float(m.group(1)) if m else float("nan")


def series_gsm_to_age(series_paths: list[str]) -> dict[str, float]:
    """Merge one or more GEO series-matrix files into {GSM: age}.

    Reads `!Sample_geo_accession` (the GSM ids, column-aligned) and the
    `!Sample_characteristics_ch1` row whose values start with 'age'. Two platforms => two files,
    unioned. Pure w.r.t. logic; only reads the files it is given.
    """
    import gzip
    import re
    out: dict[str, float] = {}
    for path in series_paths:
        gsm: list[str] = []
        age_row: list[str] | None = None
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("!Sample_geo_accession"):
                    gsm = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                elif line.startswith("!Sample_characteristics_ch1"):
                    vals = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                    if vals and re.match(r"\s*age\b", vals[0], re.I):
                        age_row = vals
        if gsm and age_row and len(gsm) == len(age_row):
            for g, a in zip(gsm, age_row, strict=True):
                v = parse_age_value(a)
                if np.isfinite(v):
                    out[g] = v
    return out


def dedup_symbols_highest_total(symbols: list[str], counts: np.ndarray) -> tuple[list[str], np.ndarray]:
    """Collapse duplicate gene symbols keeping the row with the highest TOTAL count.

    `counts` is (n_genes, n_samples). Matches the clock's `dedup: highest_expressed` and, crucially,
    guarantees the symbol list handed to `predict_age` is UNIQUE -- otherwise a symbol appearing
    twice would have its weight counted twice. Pure.
    """
    totals = counts.sum(axis=1)
    best: dict[str, int] = {}
    for i, s in enumerate(symbols):
        if s not in best or totals[i] > totals[best[s]]:
            best[s] = i
    keep = [best[s] for s in dict.fromkeys(symbols) if s in best]  # stable, one row per symbol
    keep_syms = [symbols[i] for i in keep]
    return keep_syms, counts[keep]


def _load_known_age_fibroblasts(known_dir: str):
    """Load the NCBI-generated GSE113957 counts, map GeneID->Symbol, attach GSM ages.

    Expects in `known_dir`:
      * `*raw_counts*NCBI*.tsv(.gz)`  -- GeneID rows x GSM columns, integer counts
      * `*annot*.tsv(.gz)`            -- GeneID -> Symbol table
      * one or more `*series_matrix*` -- GSM -> age (both platforms)
    Returns (ages[n], counts[n_samples, n_genes] in the SAME order, gene_symbols) so the caller can
    run the production `normalize_counts` -> clock path. (None, None, None) if the files are absent.
    """
    import glob

    import pandas as pd

    def _first(pattern):
        hits = glob.glob(str(Path(known_dir) / pattern))
        return hits[0] if hits else None

    counts_f = _first("*raw_counts*NCBI*.tsv*") or _first("*raw_counts*.tsv*")
    annot_f = _first("*annot*.tsv*")
    series = glob.glob(str(Path(known_dir) / "*series_matrix*"))
    if not counts_f or not annot_f or not series:
        return None, None, None

    gsm_age = series_gsm_to_age(series)                              # {GSM: age}
    annot = pd.read_csv(annot_f, sep="\t", usecols=["GeneID", "Symbol"], dtype={"GeneID": str})
    id2sym = dict(zip(annot["GeneID"].astype(str), annot["Symbol"].astype(str), strict=True))

    cdf = pd.read_csv(counts_f, sep="\t", dtype={"GeneID": str}).set_index("GeneID")
    # keep only GSM columns that have a known age, in a fixed order
    sample_cols = [c for c in cdf.columns if c in gsm_age]
    if len(sample_cols) < 5:
        return None, None, None
    ages = np.array([gsm_age[c] for c in sample_cols], dtype=float)

    symbols = [id2sym.get(gid, "") for gid in cdf.index.astype(str)]
    counts = cdf[sample_cols].to_numpy(dtype=np.float64)            # genes x samples
    keep = [i for i, s in enumerate(symbols) if s]                  # drop unmapped GeneIDs
    symbols = [symbols[i] for i in keep]
    counts = counts[keep]
    symbols, counts = dedup_symbols_highest_total(symbols, counts)  # unique symbols
    return ages, counts.T, symbols                                  # (n_samples, n_genes)


def main() -> int:
    gill_dir = sys.argv[1] if len(sys.argv) > 1 else r"D:\Gill"
    known_dir = sys.argv[2] if len(sys.argv) > 2 else r"D:\GSE113957"
    print("STAGE 1.5 §9 — CLOCK VALIDITY: broken, mis-applied, or out-of-domain?\n")
    print("  PRE-REGISTERED BARS (ground rule §5b):")
    for b in bars():
        print(f"    {b['id']}: {b['bar']}\n        decides: {b['decides']}")
    print()

    try:
        g = run_gill_checks(gill_dir)
    except Exception as exc:  # noqa: BLE001
        print(f"  [!] Gill checks could not run ({exc!r}). Pass the dir as arg 1.")
        return 1
    repro = run_reproduction(known_dir)

    cov, inr = g["coverage"], g["in_range"]
    print("  H1 — GENE COVERAGE (is the clock fully applied?)")
    print(f"     {cov['n_overlap']}/{cov['n_clock_genes']} clock genes present "
          f"({cov['frac_genes_present']:.0%}); {cov['frac_abs_weight_present']:.0%} of |weight| kept "
          f"-> {cov['verdict']['status']}")
    print(f"     {cov['verdict']['reason']}")
    print("\n  H1 — OWN-DOMAIN REPRODUCTION (gold check)")
    print(f"     {repro['status']}: {repro['reason']}")
    print("\n  H2 — IN-RANGE AGE TRACKING")
    print("     baselines: " + ", ".join(f"{d}({DONOR_AGE.get(d)})={p:.0f}"
                                          for d, p in g["baseline_pred"].items()))
    print(f"     {inr['status']}: {inr['reason']}")
    print(f"\n  H1 — intercept dominance: {g['intercept_dominance']['status']} "
          f"({g['intercept_dominance']['reason']})")
    print(f"  H1 — CP10k denominator: {g['denominator_sensitivity']['status']} "
          f"({g['denominator_sensitivity']['reason']})")
    print("\n  H3 — DIRECTIONAL ATTRIBUTION (why does age rise in reprogramming?)")
    att = g["attribution"]
    print(f"     {att['status']}: {att.get('reason', '')}")
    if "top_contributors" in att:
        print("     top gene contributors to the rise: "
              + ", ".join(f"{gname}{c:+.2f}" for gname, c in att["top_contributors"][:8]))

    decision = decide(cov["verdict"], inr, repro, att)
    print(f"\n  ==> ACTION: {decision['action']}\n      {decision['reason']}")

    out = {"script": "diag_clock_validity", "utc": datetime.now(UTC).isoformat(timespec="seconds"),
           "gill_dir": gill_dir, "known_age_dir": known_dir, "bars": bars(),
           "coverage": cov, "reproduction": repro, "in_range": inr,
           "intercept_dominance": g["intercept_dominance"],
           "denominator_sensitivity": g["denominator_sensitivity"],
           "attribution": att, "baseline_pred": g["baseline_pred"],
           "baseline_pred_restricted": g["baseline_pred_restricted"], "decision": decision}
    _RESULTS / "diag_clock_validity_results.json".write_text(json.dumps(out, indent=2, default=str),
                                                        encoding="utf-8")
    print("\n  wrote diag_clock_validity_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
