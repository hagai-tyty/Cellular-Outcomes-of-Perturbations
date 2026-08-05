"""What IS HFF's ΔAge label? Within-timepoint identity decomposition. (GSE242423)

    python experiments/diag_hff_label_identity.py "D:\\GSE242423\\GSE242423"

READ-ONLY. Writes `results/diag_hff_label_identity_results.json`. `src/` untouched, no labels move.

THE QUESTION
------------
HFF supplies **33,613 of the project's 33,688 age labels (99.8 %)**, and what those labels MEAN is
unmeasured. Two facts sit in tension:

  * arm C: permuting them collapsed `rank_model_dage` 0.9476 -> 0.5765, **5.4x the entire A-B gap**.
    So they carry structure the model exploits -- "of unknown provenance" (STEP6 report).
  * G-c step 1: the trajectory has the SHAPE of real rejuvenation (rho -0.905 vs day, stable under
    leave-one-timepoint-out, reaching -24.0 yr against methylation's verified -24.1) -- but the
    identity artefact produces that same monotone shape, so shape cannot separate them.

This measures the provenance directly:

    WITHIN a single timepoint, how much of a cell's ΔAge is explained by its IDENTITY markers?

**Within-timepoint is the whole point.** Between timepoints ΔAge and pluripotency move together
trivially -- every cell is further along -- so a pooled correlation is uninformative by construction.
Within a timepoint, day is held constant and the question becomes: among cells at the SAME point in
the protocol, does the one that looks more pluripotent read younger?

WHY THIS ONE CAN RESOLVE WHEN NOTHING ELSE HAS
----------------------------------------------
Every dead end in this project is donor-limited: 3 donors, 6 folds, MDE = 1.049 x SD. HFF is the one
place that is not. **~4,800 cells per timepoint, 8 timepoints.** SE(rho) ~ 1/sqrt(n) ~ 0.014, which
makes any bar in [0.1, 0.9] trivially resolvable. This is the only well-powered question left in the
data, and nobody has asked it.

PRE-REGISTERED BEFORE THE FIRST RUN
-----------------------------------
  primary   Spearman(y_age, pluripotency) computed WITHIN each timepoint
  bar       |rho| >= 0.50 in at least 6 of 8 timepoints  =>  IDENTITY_DOMINATED
  also      R^2 of y_age on [pluripotency, somatic identity] within timepoint -- the fraction of the
            label that is identity rather than anything else
  reported  every timepoint, always. No timepoint is selected on its value.

  IDENTITY_DOMINATED   the label is an identity readout at cell level. Arm C's "exploitable
                       structure" is then pluripotency structure, and calling it age is unjustified.
  NOT_DOMINATED        ΔAge carries within-timepoint structure that identity does not explain. That
                       residual is the candidate age signal and is what arm C detected.

**Neither outcome moves a label.** This measures what the existing labels are.

HONEST LIMIT, STATED UP FRONT
-----------------------------
A real rejuvenation signal would ALSO correlate with pluripotency at cell level, because a cell
further along is both more pluripotent and more rejuvenated. So a strong correlation does not by
itself prove artefact. What the R^2 adds is the quantity that does matter: if identity explains
essentially ALL of the label's within-timepoint variance, there is nothing left for age to be.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

RHO_BAR = 0.50            # |rho| at or above this counts a timepoint as identity-dominated
MIN_TIMEPOINTS = 6        # ... in at least this many of the 8
IPSC_DAY = 21.0           # a cell-type change, not a timepoint on the trajectory -- excluded
MIN_GENES = 500           # the source's own empty-droplet gate
MAX_CELLS = 6000          # per timepoint; the recorded build kept ~4,800
# `cells_per_run=None` means ONE chunk holding every cell, densified: ~48k cells x 36,601 genes is
# ~7 GB and the process simply thrashes. The source's own docstring says this parameter "bounds peak
# RAM"; batching costs nothing here because every statistic is computed per timepoint AFTER the
# chunks are concatenated, so where the batch boundaries fall cannot change a result.
CELLS_PER_RUN = 4000


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- #
# Pure logic — unit-tested with no repo data present                           #
# --------------------------------------------------------------------------- #
def _ranks(v: np.ndarray) -> np.ndarray:
    order = np.argsort(v, kind="mergesort")
    r = np.empty(len(v), float)
    r[order] = np.arange(1, len(v) + 1)
    for val in np.unique(v):
        m = v == val
        if m.sum() > 1:
            r[m] = r[m].mean()
    return r


def spearman(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3:
        return float("nan")
    a, b = _ranks(x), _ranks(y)
    a, b = a - a.mean(), b - b.mean()
    den = np.sqrt((a**2).sum() * (b**2).sum())
    return float(a @ b / den) if den else float("nan")


def r2_on(y: np.ndarray, X: np.ndarray) -> float:
    """Fraction of y's variance explained by X (intercept included). Ranks in, so it is monotone-
    robust and comparable with the Spearman beside it."""
    y = _ranks(np.asarray(y, float))
    Z = np.column_stack([np.ones(len(y))] + [_ranks(np.asarray(c, float)) for c in np.asarray(X).T])
    beta, *_ = np.linalg.lstsq(Z, y, rcond=None)
    resid = y - Z @ beta
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return float(1.0 - (resid**2).sum() / ss_tot) if ss_tot else float("nan")


def decide(per_day: dict, bar: float = RHO_BAR, need: int = MIN_TIMEPOINTS) -> dict:
    """Pre-registered rule: identity-dominated iff |rho| >= bar in at least `need` timepoints."""
    hits = [d for d, v in per_day.items() if np.isfinite(v["rho_pluri"]) and abs(v["rho_pluri"]) >= bar]
    ok = len(hits) >= need
    return {
        "verdict": "IDENTITY_DOMINATED" if ok else "NOT_DOMINATED",
        "n_timepoints_over_bar": len(hits),
        "n_timepoints": len(per_day),
        "timepoints_over_bar": sorted(hits),
        "reason": (f"|rho(y_age, pluripotency)| >= {bar} within timepoint in {len(hits)} of "
                   f"{len(per_day)} timepoints (rule: >= {need})"),
    }


# --------------------------------------------------------------------------- #
# Real-data wiring                                                             #
# --------------------------------------------------------------------------- #
def build(hff_dir: Path):
    """Per-cell (day, y_age, pluripotency, somatic) for HFF, through the pipeline's own loader."""
    from cellfate.common import constants as C
    from cellfate.data.aging import LinearClock
    from cellfate.data.normalize import normalize_counts
    from cellfate.data.qc import QCConfig, apply_qc
    from cellfate.data.sources import GSE242423SingleCellSource

    dcv = _load("dcv", ROOT / "experiments" / "diag_clock_validity.py")
    genes_file = next(hff_dir.glob("*genes.tsv.gz"))
    samples = []
    for mtx in sorted(hff_dir.glob("*.matrix.mtx.gz")):
        label = mtx.name.split(".")[0].split("_")[-1]
        bc = mtx.with_name(mtx.name.replace(".matrix.mtx.gz", ".barcodes.tsv.gz"))
        if bc.exists():
            samples.append({"matrix": str(mtx), "barcodes": str(bc), "label": label})
    src = GSE242423SingleCellSource(samples, str(genes_file), min_genes=MIN_GENES,
                                    max_cells_per_sample=MAX_CELLS,
                                    cells_per_run=CELLS_PER_RUN)
    clock = LinearClock.from_json(ROOT / "configs" / "clocks" / "fleischer_clock.json")
    qc = QCConfig(max_mito_frac=0.20, min_genes=MIN_GENES)

    days, ages, pluri, somatic = [], [], [], []
    lib, ngene, mito = [], [], []      # technical covariates: is the residual just depth?
    for chunk in src.plan():
        raw = apply_qc(src.fetch(chunk), qc)
        if len(raw.obs) == 0:
            continue
        norm = normalize_counts(raw.counts)
        age = clock.predict_age(norm, raw.genes)
        gi = {g: i for i, g in enumerate(raw.genes)}
        pi = [gi[g] for g in dcv.OSKM_PLURIPOTENCY if g in gi]
        si = [gi[g] for g in C.DEFAULT_SIGNATURES["safe"] if g in gi]
        cnt = np.asarray(raw.counts.todense() if hasattr(raw.counts, "todense")
                         else raw.counts, dtype=np.float64)
        tot = cnt.sum(axis=1)
        mt = [i for i, g in enumerate(raw.genes) if str(g).upper().startswith("MT-")]
        lib.append(tot)
        ngene.append((cnt > 0).sum(axis=1).astype(float))
        mito.append(cnt[:, mt].sum(axis=1) / np.maximum(tot, 1.0) if mt else np.zeros(len(tot)))
        days.append(raw.obs["time_h"].to_numpy(dtype=float) / 24.0)
        ages.append(np.asarray(age, float))
        pluri.append(norm[:, pi].mean(axis=1) if pi else np.zeros(len(age)))
        somatic.append(norm[:, si].mean(axis=1) if si else np.zeros(len(age)))
    return {"day": np.concatenate(days), "age": np.concatenate(ages),
            "pluri": np.concatenate(pluri), "somatic": np.concatenate(somatic),
            "lib": np.concatenate(lib), "ngene": np.concatenate(ngene),
            "mito": np.concatenate(mito), "n_pi": len(pi), "n_si": len(si)}


def main() -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    if len(sys.argv) < 2:
        print(__doc__)
        return 2

    D = build(Path(sys.argv[1]))
    keep = D["day"] != IPSC_DAY
    day, age, pluri, somatic = D["day"][keep], D["age"][keep], D["pluri"][keep], D["somatic"][keep]
    lib, ngene, mito = D["lib"][keep], D["ngene"][keep], D["mito"][keep]
    n_pi, n_si = D["n_pi"], D["n_si"]
    print(f"\n[shape before statistic] {len(age)} cells, "
          f"{len(np.unique(day))} timepoints; {n_pi} pluripotency genes, {n_si} somatic genes")
    print(f"  PRE-REGISTERED: |rho| >= {RHO_BAR} within timepoint in >= {MIN_TIMEPOINTS} of them "
          f"=> IDENTITY_DOMINATED\n")

    per_day = {}
    print(f"  {'day':>5} {'n':>6} {'rho(age,pluri)':>15} {'rho(age,somatic)':>17} {'R2 ident':>12}"
          f" {'rho(lib)':>9} {'rho(ng)':>8} {'R2 tech':>9} {'R2 both':>9}")
    for d in sorted(np.unique(day)):
        m = day == d
        if m.sum() < 30:
            continue
        rp = spearman(age[m], pluri[m])
        rs = spearman(age[m], somatic[m])
        r2 = r2_on(age[m], np.column_stack([pluri[m], somatic[m]]))
        # TECHNICAL: is the 84-98% that identity does not explain simply sequencing depth?
        rl = spearman(age[m], lib[m])
        rg = spearman(age[m], ngene[m])
        rm = spearman(age[m], mito[m])
        r2t = r2_on(age[m], np.column_stack([lib[m], ngene[m], mito[m]]))
        r2b = r2_on(age[m], np.column_stack([pluri[m], somatic[m], lib[m], ngene[m], mito[m]]))
        per_day[str(int(d))] = {"n": int(m.sum()), "rho_pluri": rp, "rho_somatic": rs,
                                "r2_identity": r2, "mean_age": float(age[m].mean()),
                                "rho_lib": rl, "rho_ngene": rg, "rho_mito": rm,
                                "r2_technical": r2t, "r2_identity_plus_technical": r2b}
        print(f"  {int(d):5d} {int(m.sum()):6d} {rp:15.3f} {rs:17.3f} {r2:12.3f}"
              f" {rl:9.3f} {rg:8.3f} {r2t:9.3f} {r2b:9.3f}")

    v = decide(per_day)
    # the between-timepoint number, for contrast only -- it is uninformative by construction
    pooled = spearman(age, pluri)
    out = {"script": "diag_hff_label_identity", "utc": datetime.now(UTC).isoformat(),
           "n_cells": int(len(age)), "rho_bar": RHO_BAR, "min_timepoints": MIN_TIMEPOINTS,
           "n_pluripotency_genes": n_pi, "n_somatic_genes": n_si,
           "per_timepoint": per_day, "pooled_rho_pluri_DESCRIPTIVE_ONLY": pooled,
           "verdict": v}
    (_RESULTS / "diag_hff_label_identity_results.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print(f"\n  pooled rho(age, pluripotency) = {pooled:+.3f}  [DESCRIPTIVE ONLY -- between-"
          f"timepoint movement is shared by construction]")
    print(f"\n  VERDICT: {v['verdict']}  -- {v['reason']}")
    print("  wrote results/diag_hff_label_identity_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
