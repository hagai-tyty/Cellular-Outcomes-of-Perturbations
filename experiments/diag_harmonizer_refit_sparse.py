"""STAGE 1.5.6 OPTION 2 — re-fit the harmonizer on the SPARSE gene space, and measure.

    python experiments/diag_harmonizer_refit_sparse.py "D:\\Gill" "D:\\GSE242423\\GSE242423"

READ-ONLY. Writes `results/diag_harmonizer_refit_sparse_results.json`. `src/` untouched, no
labels move.

WHAT "RE-FIT ON THE SPARSE SPACE" ACTUALLY CHANGES — read `harmonize.py` before assuming
--------------------------------------------------------------------------------------
`Harmonizer.fit` computes, per dataset:

    mu    = MG.mean(axis=0)                  # harmonize.py:110   PER GENE
    sigma = MG.std(axis=0)                   # harmonize.py:111   PER GENE
    floor = float(np.median(sigma))          # harmonize.py:112   OVER THE GENE SET
    sigma = np.maximum(sigma, floor)         # harmonize.py:113

`mu_g` and `sigma_g` are **per-gene** statistics: restricting the gene set does not move a
single one of them. The gene set enters in exactly two places, and only two:

  1. **the variance floor** — `median(sigma)` is a SET-level statistic, so it moves
  2. **the admissible mask** — `genes_G` is the intersection of per-dataset sets with mean
     control expression >= 0.1 (`harmonize.py:88, 91`); genes outside it are not in the
     harmonizer's space at all, so the clock never sees them

**So option 2 is, mechanically, a change to the variance floor and nothing else.** If the
top-100 clock genes all sit above both floors already, re-fitting is a no-op and option 2
cannot work. That is a real possible outcome and this script must be able to return it.

A DISCREPANCY IN STEPS 1c/1d THAT THIS SCRIPT FIXES
---------------------------------------------------
`diag_harmonization_gain.py` floored sigma at the median over the **shared clock genes**. The
pipeline floors at the median over the **full admissible space** (`genes_G`), which is a
different and much larger set, and it applies an expression floor that the clock-gene set
never had. Steps 1c/1d therefore measured a *near-pipeline* gain, not the pipeline's own.
This script computes all three so the numbers can be reconciled rather than quietly replaced:

    A  PIPELINE     floor = median sigma over genes_G, admissibility applied   <- the truth
    B  CLOCK-SPACE  floor = median sigma over clock genes, no expr floor       <- what 1c/1d did
    C  SPARSE REFIT floor = median sigma over the top-k clock genes            <- option 2

THE PRE-REGISTERED BAR, WRITTEN BEFORE THIS RAN
------------------------------------------------
Step 1d found the top-100 clock has a LARGER gain than the dense clock (2.769 vs 2.152), so
sparsification made HFF worse. **Option 2 succeeds only if, under regime C, the top-100 gain
is no larger than the dense clock's gain in the same regime.** Bringing the gain down a
little while it still exceeds the dense clock's does not rescue §1 — it would still mean
sparsifying degrades HFF, which is 99.8 % of the age labels.

A FAIL-OPEN DEFAULT IN 1c/1d, FOUND AND REMOVED HERE
-----------------------------------------------------
`diag_harmonization_gain.py` picked Gill's controls with

    g_day = [float(m.group(1)) if (m := re.search(r"_d(\\d+)_", s)) else 0.0 for s in samples]
    g_ctrl = g_norm[g_day == 0.0]                     # <- anything UNPARSEABLE becomes a control

Checked: 118 of 124 Gill sample names carry `_dNN_`, and **not one carries `_d0_`**. The 6 that
do not parse are `N2/N3/O1/O2/Y1/Y2_Fib_Sendai_Exp2` — the day-0 dermal fibroblasts, which are
exactly what `sources.py:417` defines as `is_control` for this source. **So 1c/1d selected the
right 6 rows, but by accident: the fail-open default happened to coincide with the truth.** Here
the controls are selected by the pipeline's own definition and the count is asserted, so a rename
upstream fails loudly instead of silently pooling the whole time course into sigma_gill.

ONE EXACTNESS NOTE
-------------------
ΔAge is linear in expression, so the mean over day-14 cells of `(x_i - base) . w . ratio`
equals `(mean_day14(x) - base) . w . ratio`. This script therefore streams per-gene sums
instead of holding a 50k x 20k matrix. That is an identity, not an approximation.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

KS = (50, 100, 150, 300, 1000, None)
EPS = 1e-6            # harmonize.py:25
EXPR_FLOOR = 0.1      # DEFAULT_EXPR_FLOOR, harmonize.py:26
DENSE_GAIN_1D = 2.152  # regime B, for reconciliation
TOP100_GAIN_1D = 2.769


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def apply_floor(sigma_raw: np.ndarray, floor_space: np.ndarray) -> tuple[np.ndarray, float]:
    """`harmonize.py:112-113` — floor at the median sigma over the FITTED gene set.

    `floor_space` is the index set the harmonizer was fitted on. Changing it is the whole of
    what option 2 does.
    """
    floor = float(np.median(sigma_raw[floor_space]))
    return np.maximum(sigma_raw, floor), floor


def main() -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    if len(sys.argv) < 3:
        print(__doc__)
        return 2

    from cellfate.data.normalize import normalize_counts
    from cellfate.data.qc import QCConfig, apply_qc
    from cellfate.data.sources import GSE242423SingleCellSource
    dv = _load("dv", ROOT / "experiments" / "diag_dage_variants.py")
    hli = _load("hli", ROOT / "experiments" / "diag_hff_label_identity.py")

    clock = json.loads((ROOT / "configs" / "clocks" / "fleischer_clock.json").read_text("utf-8"))
    W = {k: float(v) for k, v in clock["weights"].items()}

    # ---- Gill: the reference. Bulk, small enough to hold ---------------------------------- #
    g_samples, g_genes, g_lin = dv.load_gill(Path(sys.argv[1]))
    g_norm = normalize_counts(g_lin)
    # `sources.py:417` -- is_control is the day-0 "Dermal fibroblast" of each donor. Select on that
    # definition, not on a regex whose failure mode is to call an unparsed sample a control.
    ctrl_mask = np.array([re.search(r"_d\d+_", s) is None and "_Fib_" in s for s in g_samples])
    n_ctrl = int(ctrl_mask.sum())
    if n_ctrl < 3:      # MIN_REPLICATES, harmonize.py:27
        print(f"FATAL: {n_ctrl} Gill control samples matched '_Fib_' with no _dNN_ token; "
              f"sigma_gill is undefined. Sample names may have changed: {g_samples[:4]}")
        return 1
    unparsed = int(sum(re.search(r"_d\d+_", s) is None for s in g_samples))
    if unparsed != n_ctrl:
        print(f"FATAL: {unparsed} samples carry no _dNN_ token but only {n_ctrl} look like "
              f"fibroblast baselines -- the control set is ambiguous, refusing to guess.")
        return 1
    G_ctrl = np.asarray(g_norm[ctrl_mask], dtype=np.float64)
    gill_mean, gill_sd = G_ctrl.mean(axis=0), G_ctrl.std(axis=0)
    print(f"\n[gill]  {len(g_samples)} samples, {n_ctrl} controls "
          f"(day-0 dermal fibroblasts: "
          f"{[s for s, m in zip(g_samples, ctrl_mask, strict=True) if m][:3]} ...)")

    # ---- HFF: stream per-gene sums. Linearity makes this exact, not an approximation ------ #
    hff = Path(sys.argv[2])
    genes_file = next(hff.glob("*genes.tsv.gz"))
    samples = []
    for mtx in sorted(hff.glob("*.matrix.mtx.gz")):
        bc = mtx.with_name(mtx.name.replace(".matrix.mtx.gz", ".barcodes.tsv.gz"))
        if bc.exists():
            samples.append({"matrix": str(mtx), "barcodes": str(bc),
                            "label": mtx.name.split(".")[0].split("_")[-1]})
    src = GSE242423SingleCellSource(samples, str(genes_file), min_genes=hli.MIN_GENES,
                                    max_cells_per_sample=hli.MAX_CELLS,
                                    cells_per_run=hli.CELLS_PER_RUN)
    qc = QCConfig(max_mito_frac=0.20, min_genes=hli.MIN_GENES)

    h_genes = None
    n0 = n14 = 0
    s0 = q0 = s14 = None
    for chunk in src.plan():
        raw = apply_qc(src.fetch(chunk), qc)
        if len(raw.obs) == 0:
            continue
        h_genes = raw.genes
        norm = np.asarray(normalize_counts(raw.counts), dtype=np.float64)
        d = raw.obs["time_h"].to_numpy(dtype=float) / 24.0
        m0 = d == 0.0
        if m0.any():
            X = norm[m0]
            s0 = X.sum(0) if s0 is None else s0 + X.sum(0)
            q0 = (X * X).sum(0) if q0 is None else q0 + (X * X).sum(0)
            n0 += int(m0.sum())
        m14 = d == 14.0
        if m14.any():
            X = norm[m14]
            s14 = X.sum(0) if s14 is None else s14 + X.sum(0)
            n14 += int(m14.sum())
    if not n0 or not n14:
        print("FATAL: no day-0 controls or no day-14 cells")
        return 1
    hff_mean = s0 / n0
    hff_sd = np.sqrt(np.maximum(q0 / n0 - hff_mean ** 2, 0.0))   # ddof=0, matches np.std
    hff_d14 = s14 / n14
    print(f"[hff]   {n0} day-0 controls, {n14} day-14 cells")

    # ---- gene spaces -------------------------------------------------------------------- #
    gi_g = {g: i for i, g in enumerate(g_genes)}
    gi_h = {g: i for i, g in enumerate(h_genes)}
    both = sorted(set(gi_g) & set(gi_h))
    ig_b = np.array([gi_g[g] for g in both])
    ih_b = np.array([gi_h[g] for g in both])

    adm = (gill_mean[ig_b] >= EXPR_FLOOR) & (hff_mean[ih_b] >= EXPR_FLOOR)
    genes_G = [g for g, a in zip(both, adm, strict=True) if a]
    igG = np.array([gi_g[g] for g in genes_G])
    ihG = np.array([gi_h[g] for g in genes_G])
    print(f"[space] {len(both)} genes in both datasets -> {len(genes_G)} admissible "
          f"(mean control expr >= {EXPR_FLOOR} in BOTH)")

    sgG, shG = gill_sd[igG], hff_sd[ihG]
    wG = np.array([W.get(g, 0.0) for g in genes_G])
    deltaG = hff_d14[ihG] - hff_mean[ihG]
    clock_mask = wG != 0.0
    print(f"[clock] {int(clock_mask.sum())} clock genes survive admissibility "
          f"(of {len(set(both) & set(W))} present in both)")

    # how many of the top-k clock genes are DROPPED by the expression floor?
    shared_nofloor = [g for g in both if g in W]
    w_nf = np.array([W[g] for g in shared_nofloor])
    order_nf = np.argsort(-np.abs(w_nf))
    inG = set(genes_G)

    order = np.argsort(np.where(clock_mask, -np.abs(wG), np.inf))   # clock genes first, by |w|
    n_clock = int(clock_mask.sum())

    # ---- regime A (pipeline) and regime C (option 2: refit on the sparse space) ---------- #
    all_G = np.arange(len(genes_G))
    per_k: dict[str, dict] = {}
    print(f"\n  {'k':>6} {'dropped':>8} | {'A gain':>8} {'A day14':>9} | {'C gain':>8} "
          f"{'C day14':>9} | {'C floor bites':>14}")
    for k in KS:
        kk = n_clock if k is None else min(k, n_clock)
        kept = order[:kk]
        n_dropped = kk - sum(1 for i in order_nf[:kk] if shared_nofloor[i] in inG)

        row = {"n_kept": int(kk), "n_top_k_dropped_by_expr_floor": int(n_dropped)}
        for tag, space in (("A_pipeline", all_G), ("C_sparse_refit", kept)):
            sg, fg = apply_floor(sgG, space)
            sh, fh = apply_floor(shG, space)
            ratio = sg / (sh + EPS)
            direct = float(deltaG[kept] @ wG[kept])
            harm = float((deltaG[kept] * ratio[kept]) @ wG[kept])
            row[tag] = {
                "floor_gill": fg, "floor_hff": fh,
                "day14_direct": direct, "day14_harmonized": harm,
                "gain": harm / direct if direct else float("nan"),
                "median_ratio_kept": float(np.median(ratio[kept])),
                "n_kept_clipped_gill": int((sgG[kept] < fg).sum()),
                "n_kept_clipped_hff": int((shG[kept] < fh).sum()),
            }
        a, c = row["A_pipeline"], row["C_sparse_refit"]
        bites = f"{c['n_kept_clipped_gill']}g/{c['n_kept_clipped_hff']}h"
        print(f"  {str(k):>6} {n_dropped:8d} | {a['gain']:8.3f} {a['day14_harmonized']:9.2f} | "
              f"{c['gain']:8.3f} {c['day14_harmonized']:9.2f} | {bites:>14}")
        per_k[str(k)] = row

    # ---- regime B: reproduce 1c/1d so the three can be reconciled ------------------------ #
    ig_s = np.array([gi_g[g] for g in shared_nofloor])
    ih_s = np.array([gi_h[g] for g in shared_nofloor])
    sg_s, _ = apply_floor(gill_sd[ig_s], np.arange(len(shared_nofloor)))
    sh_s, _ = apply_floor(hff_sd[ih_s], np.arange(len(shared_nofloor)))
    ratio_s = sg_s / (sh_s + EPS)
    delta_s = hff_d14[ih_s] - hff_mean[ih_s]
    regime_b = {}
    for k in (100, None):
        kept = order_nf[:len(shared_nofloor) if k is None else k]
        dd = float(delta_s[kept] @ w_nf[kept])
        dh = float((delta_s[kept] * ratio_s[kept]) @ w_nf[kept])
        regime_b[str(k)] = {"day14_direct": dd, "day14_harmonized": dh,
                            "gain": dh / dd if dd else float("nan")}
    print(f"\n  [regime B reconciliation] dense gain {regime_b['None']['gain']:.3f} "
          f"(1d said {DENSE_GAIN_1D}), top100 gain {regime_b['100']['gain']:.3f} "
          f"(1d said {TOP100_GAIN_1D})")

    # ---- the pre-registered verdict ------------------------------------------------------ #
    c100 = per_k["100"]["C_sparse_refit"]["gain"]
    cdense = per_k["None"]["C_sparse_refit"]["gain"]
    a100 = per_k["100"]["A_pipeline"]["gain"]
    adense = per_k["None"]["A_pipeline"]["gain"]
    passed = abs(c100) <= abs(adense)
    print(f"\n  PIPELINE (A): dense gain {adense:.3f}, top100 gain {a100:.3f}")
    print(f"  OPTION 2 (C): top100 gain {c100:.3f}  (dense, refit on its own space: {cdense:.3f})")
    print(f"  BAR: |top100 gain after refit| <= |dense gain| = {abs(adense):.3f}")
    verdict = "OPTION_2_WORKS" if passed else "OPTION_2_FAILS"
    if per_k["100"]["C_sparse_refit"]["n_kept_clipped_gill"] == 0 and \
            per_k["100"]["C_sparse_refit"]["n_kept_clipped_hff"] == 0:
        verdict += "__REFIT_IS_A_NO_OP"

    out = {"script": "diag_harmonizer_refit_sparse", "utc": datetime.now(UTC).isoformat(),
           "n_genes_both": len(both), "n_genes_admissible": len(genes_G),
           "n_clock_genes_admissible": n_clock,
           "n_clock_genes_in_both": len(shared_nofloor),
           "n_day0_controls": n0, "n_day14_cells": n14,
           "per_k": per_k, "regime_B_reconciliation": regime_b,
           "bar": "abs(top100 gain, refit) <= abs(dense gain, pipeline)",
           "verdict": verdict}
    (_RESULTS / "diag_harmonizer_refit_sparse_results.json").write_text(
        json.dumps(out, indent=2), "utf-8")
    print(f"\n  => {verdict}")
    print("  wrote results/diag_harmonizer_refit_sparse_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
