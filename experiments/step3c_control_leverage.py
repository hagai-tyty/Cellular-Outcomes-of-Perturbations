"""
STAGE 1.5.6 STEP 3c — is the degenerate control the carrier of the fold spread?

    python experiments/step3c_control_leverage.py "D:/Gill" "D:/GSE242423"

READ-ONLY. Writes `results/step3c_control_leverage_results.json`. `src/` untouched, no build
touched, no label moved. Pre-registration: `plans/STAGE_1_5_6_SPARSE_CLOCK.md` §5.8.

THE TEST (§5.8, unchanged)
--------------------------
On the **O1 fold** -- July's reference, `d_O1 = -24.023`, and one of the five folds that INCLUDE
the degenerate control -- refit the harmonizer FIVE times, each dropping ONE of its five control
samples (N2, N3, O2, Y1, Y2), and recompute HFF's day-14 dAge each time.

`MIN_REPLICATES = 3` (`harmonize.py:27`), so four controls remain legal in every arm.

**The four healthy drops are the built-in negative control.** Dropping ANY control moves
`sigma_gill`. The claim under test is not "N2's removal moves the number" but "**N2's removal is an
OUTLIER among the five**" -- measured in the same run, by the same code, on the same fold.

THE RECONSTRUCTION
------------------
`Harmonizer.transform` z-scores against the cell's own dataset, `project_to_clock` reverses into
the reference scale (`harmonize.py:117-131`), so for two cells of the SAME dataset the projected
difference is

    x'_1 - x'_2  =  (x_1 - x_2) * sigma_ref / (sigma_d + EPS)

and therefore, over the harmonizer's gene space with the frozen clock weights `w`:

    d_hat  =  SUM_g  delta_g * ( sigma_gill,g / (sigma_hff,g + EPS) ) * w_g
    delta_g = mean_day14(x_hff,g) - mean_day0(x_hff,g)      in log1p-CP10k space

**`delta` is COMMON to all six arms.** Only `sigma_gill` changes between them, so every
`Delta_k = d_hat^(-k) - d_hat^(5)` is exact regardless of how HFF's cells were sampled. That
robustness is why B1 is the primary bar and B2 the secondary one.

WHAT IS HELD FIXED, AND WHY IT IS STATED
-----------------------------------------
The gene set and `sigma_hff` are held at the O1 fold's SHIPPED values
(`runs/cellfate_multi/harmonization.json`). HFF's controls are fold-invariant
(`fit_harmonizer`: `keep = is_ctrl & ~cell_line.isin(heldout)`, decidable from `cell_line` alone),
so `sigma_hff` genuinely does not move; the gene set is held by convention, matching §5.6's
T3-only rung, so that this measures the control's effect THROUGH `sigma_gill` and not through the
admissible mask. §5.5 already showed the mask is not a carrier (F = 0.957 at maximum leverage).

BAR (§5.8, pre-registered)
---------------------------
  B1 OUTLIER (primary)   |Delta_N2| is the LARGEST of the five AND >= 2x the second largest
  B2 MAGNITUDE (primary) gap closed A = Delta_N2 / (d_N2 - d_O1) >= 0.70,  denominator +16.671 yr
  B3 DIRECTION (gate)    Delta_N2 > 0 -- removing the contaminant must move dAge TOWARD zero.
                         The wrong sign falsifies the mechanism outright.

Stated in advance: "O1 minus N2's control" is NOT the N2 fold. The N2 fold holds out N2 and
therefore INCLUDES O1's control, on a different admissible set. Exact reproduction of -7.352 is
NOT predicted and is not the bar.

BRANCHES: B1^B2^B3 ATTRIBUTED | B1^B3 without B2 PARTIAL | B1 fails GENERIC | B3 fails FALSIFIED
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "step3c_control_leverage_results.json"

SHIPPED = ROOT / "runs" / "cellfate_multi" / "harmonization.json"
EPS = 1e-8
DAY14, DAY0 = 14.0, 0.0
D_O1_OBSERVED = -24.023          # §4.7, the O1 fold's recorded HFF day-14
D_N2_OBSERVED = -7.352           # §4.7, the N2 fold's
GAP = D_N2_OBSERVED - D_O1_OBSERVED      # +16.671 yr
O1_FOLD_CONTROLS = ["N2", "N3", "O2", "Y1", "Y2"]     # O1 is held out, so its own control is absent
B2_MIN_GAP_CLOSED = 0.70
B1_OUTLIER_RATIO = 2.0


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def gill_controls(gill_dir: Path, genes_G: list[str]) -> dict[str, np.ndarray]:
    """Each donor's day-0 `_Fib_` control, normalised, aligned to the harmonizer's gene space."""
    from cellfate.data.normalize import normalize_counts
    dv = _load("dv", ROOT / "experiments" / "diag_dage_variants.py")
    samples, genes, lin = dv.load_gill(gill_dir)
    norm = np.asarray(normalize_counts(lin), dtype=np.float64)
    idx = {g: i for i, g in enumerate(genes)}
    cols = np.array([idx.get(g, -1) for g in genes_G])
    if (cols < 0).any():
        raise SystemExit(f"{int((cols < 0).sum())} shipped genes missing from the Gill matrix")
    out = {}
    for i, s in enumerate(samples):
        if re.search(r"_d\d+_", s) is None and "_Fib_" in s:
            out[s.split("_")[0]] = norm[i, cols]
    return out


def hff_delta(hff_dir: Path, genes_G: list[str]) -> tuple[np.ndarray, int, int]:
    """Per-gene mean(day14) - mean(day0) for HFF, log1p-CP10k, on the harmonizer's gene space."""
    from cellfate.data.normalize import normalize_counts
    from cellfate.data.qc import QCConfig, apply_qc
    from cellfate.data.sources import GSE242423SingleCellSource
    hli = _load("hli", ROOT / "experiments" / "diag_hff_label_identity.py")

    genes_file = next(hff_dir.glob("*genes.tsv.gz"))
    samples = []
    for mtx in sorted(hff_dir.glob("*.matrix.mtx.gz")):
        bc = mtx.with_name(mtx.name.replace(".matrix.mtx.gz", ".barcodes.tsv.gz"))
        if bc.exists():
            samples.append({"matrix": str(mtx), "barcodes": str(bc),
                            "label": mtx.name.split(".")[0].split("_")[-1]})
    src = GSE242423SingleCellSource(samples, str(genes_file), min_genes=hli.MIN_GENES,
                                    max_cells_per_sample=hli.MAX_CELLS,
                                    cells_per_run=hli.CELLS_PER_RUN)
    qc = QCConfig(max_mito_frac=0.20, min_genes=hli.MIN_GENES)

    n = len(genes_G)
    s14 = np.zeros(n)
    s0 = np.zeros(n)
    c14 = c0 = 0
    for chunk in src.plan():
        raw = apply_qc(src.fetch(chunk), qc)
        if len(raw.obs) == 0:
            continue
        day = raw.obs["time_h"].to_numpy(dtype=float) / 24.0
        m14 = np.isclose(day, DAY14)
        m0 = np.isclose(day, DAY0)
        if not (m14.any() or m0.any()):
            continue
        norm = np.asarray(normalize_counts(raw.counts), dtype=np.float64)
        idx = {g: i for i, g in enumerate(raw.genes)}
        cols = np.array([idx.get(g, -1) for g in genes_G])
        take = cols >= 0
        if m14.any():
            blk = norm[m14][:, cols[take]]
            s14[take] += blk.sum(axis=0)
            c14 += int(m14.sum())
        if m0.any():
            blk = norm[m0][:, cols[take]]
            s0[take] += blk.sum(axis=0)
            c0 += int(m0.sum())
    if c14 == 0 or c0 == 0:
        raise SystemExit(f"HFF stream produced day14={c14} day0={c0} cells")
    return s14 / c14 - s0 / c0, c14, c0


def sigma_floored(mat: np.ndarray) -> tuple[np.ndarray, float]:
    """`harmonize.py:111-113` -- per-gene std over control observations, floored at its median."""
    sig = mat.std(axis=0)
    floor = float(np.median(sig))
    return np.maximum(sig, floor), floor


def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    gill_dir, hff_dir = Path(sys.argv[1]), Path(sys.argv[2])
    if not SHIPPED.exists():
        print(f"FATAL: {SHIPPED} missing -- 3c gates on the shipped O1 harmonizer.")
        return 1

    ship = json.loads(SHIPPED.read_text("utf-8"))
    genes_G = list(ship["genes"])
    sig_hff = np.asarray(ship["stats"]["hff_sc"]["sigma"], float)
    sig_gill_ship = np.asarray(ship["stats"]["gill_bulk"]["sigma"], float)

    clock = json.loads((ROOT / "configs" / "clocks" / "fleischer_clock.json").read_text("utf-8"))
    W = {k: float(v) for k, v in clock["weights"].items()}
    w = np.array([W.get(g, 0.0) for g in genes_G])

    print("\n" + "=" * 78)
    print("STEP 3c — leave-one-CONTROL-out on the O1 fold's harmonizer")
    print("=" * 78)
    print(f"gene space {len(genes_G)}   clock genes present {int((w != 0).sum())}")

    ctrl = gill_controls(gill_dir, genes_G)
    missing = [d for d in O1_FOLD_CONTROLS if d not in ctrl]
    if missing:
        print(f"FATAL: controls not found: {missing}")
        return 1
    print(f"Gill controls found: {sorted(ctrl)}   O1 fold uses {O1_FOLD_CONTROLS}")

    # G0 — the five-control fit must reproduce the shipped sigma_gill
    base_mat = np.stack([ctrl[d] for d in O1_FOLD_CONTROLS])
    sig5, floor5 = sigma_floored(base_mat)
    unclamped = sig5 > floor5 * (1 + 1e-12)
    rel = np.abs(sig5[unclamped] - sig_gill_ship[unclamped]) / np.maximum(
        np.abs(sig_gill_ship[unclamped]), 1e-12)
    g0 = float(np.median(rel)) if rel.size else float("nan")
    g0_pass = bool(g0 < 1e-6)
    print(f"\n[G0] five-control fit vs shipped sigma_gill on {int(unclamped.sum())} unclamped "
          f"genes: median rel.err {g0:.2e}  -> {'PASS' if g0_pass else 'FAIL'}")
    if not g0_pass:
        print("   G0 failed; nothing downstream is read.")
        return 1

    delta, n14, n0 = hff_delta(hff_dir, genes_G)
    print(f"[HFF] day-14 cells {n14}, day-0 cells {n0}")

    def dhat(sig_gill: np.ndarray) -> float:
        return float((delta * (sig_gill / (sig_hff + EPS)) * w).sum())

    d_base = dhat(sig5)
    arms = {}
    for k in O1_FOLD_CONTROLS:
        keep = [d for d in O1_FOLD_CONTROLS if d != k]
        sig, fl = sigma_floored(np.stack([ctrl[d] for d in keep]))
        arms[k] = {"d_hat": dhat(sig), "floor_gill": fl, "n_controls": len(keep)}
        arms[k]["delta_yr"] = arms[k]["d_hat"] - d_base

    print(f"\n  baseline (all 5 controls)  d_hat = {d_base:+.3f} yr"
          f"   [recorded d_O1 = {D_O1_OBSERVED:+.3f}]")
    print("\n" + render_table(
        ["dropped", "n ctrl", "floor_gill", "d_hat", "Delta (yr)"],
        [[k, str(arms[k]["n_controls"]), f"{arms[k]['floor_gill']:.5f}",
          f"{arms[k]['d_hat']:+.3f}", f"{arms[k]['delta_yr']:+.3f}"]
         for k in O1_FOLD_CONTROLS], aligns=["l", "r", "r", "r", "r"]))

    mags = {k: abs(arms[k]["delta_yr"]) for k in O1_FOLD_CONTROLS}
    ordered = sorted(mags, key=lambda k: -mags[k])
    d_n2 = arms["N2"]["delta_yr"]
    second = mags[ordered[1]]
    b1 = bool(ordered[0] == "N2" and mags["N2"] >= B1_OUTLIER_RATIO * second)
    a = d_n2 / GAP
    b2 = bool(a >= B2_MIN_GAP_CLOSED)
    b3 = bool(d_n2 > 0)

    print("\n  ranked |Delta|: " + "  ".join(f"{k} {mags[k]:.3f}" for k in ordered))
    print(f"\n  B1 OUTLIER   N2 largest and >= {B1_OUTLIER_RATIO}x the second "
          f"({mags['N2']:.3f} vs {second:.3f}, ratio "
          f"{mags['N2']/second if second else float('inf'):.2f}x)  -> {'PASS' if b1 else 'FAIL'}")
    print(f"  B2 MAGNITUDE gap closed A = {a:+.3f} (bar >= {B2_MIN_GAP_CLOSED})"
          f"                          -> {'PASS' if b2 else 'FAIL'}")
    print(f"  B3 DIRECTION Delta_N2 = {d_n2:+.3f} > 0"
          f"                                     -> {'PASS' if b3 else 'FAIL'}")

    if not b3:
        verdict = "FALSIFIED"
    elif not b1:
        verdict = "GENERIC"
    elif b1 and b2:
        verdict = "ATTRIBUTED"
    else:
        verdict = "PARTIAL"

    print(f"\n  VERDICT: {verdict}")
    msg = {
        "ATTRIBUTED": "the degenerate control IS the carrier. Step 3b's ladder is unnecessary; "
                      "the fix is a data fix, leakage-free and NOT donor-blocked.",
        "PARTIAL": "N2's removal is special but carries <70% of the gap. 3b runs on the residue, "
                   "with the contaminant now a known term rather than a hypothesis.",
        "GENERIC": "every control drop moves it comparably -- the contaminant is NOT the carrier. "
                   "§5.7 stands as a data defect with its own owner; 3b runs as written.",
        "FALSIFIED": "removing the contaminant moved dAge the WRONG WAY. §5.7's mechanism is "
                     "wrong. Record it, do not rescue it; 3b runs as written.",
    }[verdict]
    print(f"  -> {msg}")
    print("\n  SCOPE: 'O1 minus N2's control' is NOT the N2 fold -- the N2 fold holds out N2 and")
    print("  therefore INCLUDES O1's control, on a different admissible set. Exact reproduction")
    print(f"  of {D_N2_OBSERVED:+.3f} was not predicted and is not the bar.")

    OUT.write_text(json.dumps({
        "script": "step3c_control_leverage", "gene_space": len(genes_G),
        "g0": {"median_rel_err": g0, "n_unclamped": int(unclamped.sum()), "pass": g0_pass},
        "hff": {"n_day14": n14, "n_day0": n0},
        "d_hat_baseline": d_base, "observed": {"d_O1": D_O1_OBSERVED, "d_N2": D_N2_OBSERVED,
                                               "gap": GAP},
        "arms": arms, "ranked_abs_delta": ordered,
        "bars": {"B1_outlier": b1, "B2_magnitude": b2, "B3_direction": b3,
                 "gap_closed_A": a, "outlier_ratio": mags["N2"] / second if second else None},
        "verdict": verdict, "reading": msg,
    }, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
