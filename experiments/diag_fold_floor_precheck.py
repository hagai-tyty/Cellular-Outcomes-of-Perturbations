"""
STAGE 1.5.6 step 3b — PRECHECK: can T1 (mask) or T2 (variance floor) carry the fold spread?

    python experiments/diag_fold_floor_precheck.py

READ-ONLY. Writes `results/diag_fold_floor_precheck_results.json`. `src/` untouched, no build
touched, no label moved. Reads only committed artefacts: each fold's `harmonization.json`,
the frozen clock, and `results/repro_hff_signature_armA_results.json`.

WHY THIS EXISTS
---------------
§5.3 proposed a three-term ablation ladder (T1 mask / T2 variance floor / T3 sigma_gill) behind a
full per-fold reconstruction. §5.4 then measured, on `runs/cellfate_multi/harmonization.json` (the
O1 fold's shipped harmonizer), that the variance floor is not a tail-trim but the transform's
CENTRE:

    floor_gill 0.15821   floor_hff 0.42388   ratio 0.3732
    1848 of 5328 genes (34.7%) clamped in BOTH datasets -> ratio identical to 12 decimals
    median ratio over all 5328 genes = 0.3732 = the floor constant itself

and predicted **T2 > T3**: if `floor_gill/floor_hff` shifts between folds, those 1848 ratios move
in LOCKSTEP -- a coherent, non-averaging perturbation -- whereas T3's per-gene sigma noise is
independent across genes and should largely cancel in a weighted sum over thousands.

§5.4's own observation is that this is checkable from **two scalars per fold**, before any
reconstruction exists. That is what this script does. It is a FALSIFIER, not the ladder.

WHAT IS COMPUTED
----------------
Per fold, from that fold's own `harmonization.json`:

    floor_ds = min(sigma_ds)        exact: sigma = max(sigma_raw, median(sigma_raw)), so the
                                    post-floor minimum IS the floor whenever anything clamped
    R_f      = floor_gill / floor_hff              <- T2's lockstep constant
    C_f      = sum_{g in G_f} |w_g| / sum_g |w_g|  <- T1's clock-weight coverage

MAXIMUM-LEVERAGE TEST (the reason a scalar can falsify a mechanism)
-------------------------------------------------------------------
The closed form is `d = sum_g delta_g * r_g * w_g`. Split the genes into the set B clamped in both
datasets -- whose ratio is exactly `R_f` -- and the rest:

    d_f  =  R_f * sum_{g in B} delta_g w_g   +   sum_{g not in B} delta_g r_g w_g

`d` is affine in `R_f` with slope `sum_B delta w`. Setting the second term to zero gives T2 its
**largest possible** leverage, and then `d_f = d_O1 * R_f / R_O1`. If even that prediction misses
the recorded `d_f` by more than the spread being explained, **T2 cannot be the carrier at any
leverage fraction.** The same construction bounds T1 using `C_f`.

READ:
  - a term's max-leverage prediction tracks d_f  -> it survives; the ladder must measure it
  - it misses by more than the spread            -> ELIMINATED, and the reconstruction can skip it

WHAT A NEGATIVE RESULT DOES *NOT* ESTABLISH
--------------------------------------------
`R_f` and `C_f` are SCALARS. Eliminating them eliminates the lockstep-constant channel and the
total-coverage channel. Floor effects that act through *which* genes get clamped, and mask effects
that act through gene IDENTITY rather than summed weight, are not scalar and are NOT tested here --
they live inside the reconstruction. This script narrows the ladder; it cannot close it.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

_RESULTS = Path(__file__).resolve().parents[1] / "results"
REPO = Path(__file__).resolve().parents[1]
OUT = _RESULTS / "diag_fold_floor_precheck_results.json"

DONORS = ["N2", "N3", "O1", "O2", "Y1", "Y2"]
REF_FOLD = "O1"            # the fold July's reference and §5.4's measurement both come from
ARM = "armA"
# A term is a plausible CARRIER only if, at maximum leverage, it removes >20% of the spread.
# Between that and F >= 1 it is not eliminated but cannot be the explanation either -- a
# distinction a two-state verdict would hide.
CARRIER_F = 0.80


def clock_abs_weights() -> dict[str, float]:
    c = json.loads((REPO / "configs" / "clocks" / "fleischer_clock.json").read_text("utf-8"))
    w = c.get("coefficients") or c.get("weights")
    if not isinstance(w, dict):
        raise SystemExit(f"unexpected clock format: keys={list(c)[:8]}")
    return {g: abs(float(v)) for g, v in w.items()}


def fold_row(donor: str, gw: dict[str, float], total_w: float) -> dict:
    h = json.loads((REPO / f"cellfate_loocv_{donor}_{ARM}" / "harmonization.json")
                   .read_text("utf-8"))
    genes = h["genes"]
    sg = np.asarray(h["stats"]["gill_bulk"]["sigma"], float)
    sh = np.asarray(h["stats"]["hff_sc"]["sigma"], float)
    fg, fh = float(sg.min()), float(sh.min())
    cg, ch = np.isclose(sg, fg), np.isclose(sh, fh)
    both = cg & ch
    ratio = sg / sh
    cover = sum(gw[g] for g in genes if g in gw) / total_w
    return {
        "donor": donor,
        "n_genes": len(genes),
        "floor_gill": fg,
        "floor_hff": fh,
        "R_floor_ratio": fg / fh,
        "frac_clamped_gill": float(cg.mean()),
        "frac_clamped_hff": float(ch.mean()),
        "frac_clamped_both": float(both.mean()),
        "n_clamped_both": int(both.sum()),
        "ratio_spread_within_both": float(ratio[both].max() - ratio[both].min()),
        "median_ratio_all": float(np.median(ratio)),
        "median_ratio_clamped_neither": float(np.median(ratio[~cg & ~ch])),
        "C_clock_weight_coverage": float(cover),
    }


def max_leverage(term: np.ndarray, d: np.ndarray, ref: int) -> dict:
    """Largest effect this scalar term could have: d_f = d_ref * term_f / term_ref."""
    pred = d[ref] / term[ref] * term
    resid = d - pred
    spread_d = float(d.max() - d.min())
    spread_r = float(resid.max() - resid.min())
    return {
        "predicted": pred.tolist(),
        "residual": resid.tolist(),
        "spread_d": spread_d,
        "spread_residual": spread_r,
        "F": spread_r / spread_d if spread_d else float("nan"),
        "worst_miss_yr": float(np.abs(resid).max()),
    }


def main() -> int:
    from scipy.stats import spearmanr

    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    repro = json.loads((_RESULTS / "repro_hff_signature_armA_results.json").read_text("utf-8"))
    d14 = repro["day14_by_fold"]

    gw = clock_abs_weights()
    total_w = float(sum(gw.values()))
    rows = [fold_row(dn, gw, total_w) for dn in DONORS]

    print("\n" + "=" * 78)
    print("PRECHECK — can the variance floor (T2) or the mask (T1) carry the 16.67 yr spread?")
    print("=" * 78)
    print(f"clock: {len(gw)} genes, sum|w| = {total_w:.3f}")

    print("\n" + render_table(
        ["fold", "genes", "floor_gill", "floor_hff", "R = ratio", "both%", "|w|-cover", "day-14"],
        [[r["donor"], str(r["n_genes"]), f"{r['floor_gill']:.5f}", f"{r['floor_hff']:.5f}",
          f"{r['R_floor_ratio']:.4f}", f"{r['frac_clamped_both']:.1%}",
          f"{r['C_clock_weight_coverage']:.4f}", f"{d14[r['donor']]:+.3f}"] for r in rows],
        aligns=["l"] + ["r"] * 7))

    ref = DONORS.index(REF_FOLD)
    d = np.array([d14[dn] for dn in DONORS])
    R = np.array([r["R_floor_ratio"] for r in rows])
    C = np.array([r["C_clock_weight_coverage"] for r in rows])

    out: dict = {"script": "diag_fold_floor_precheck", "ref_fold": REF_FOLD,
                 "per_fold": rows, "day14": {dn: d14[dn] for dn in DONORS}}

    for name, term, label in [("T2_floor", R, "T2  variance floor (R = floor_gill/floor_hff)"),
                              ("T1_mask", C, "T1  mask (clock-weight coverage)")]:
        ml = max_leverage(term, d, ref)
        rho = float(spearmanr(term, d).correlation)
        # Three states, because "not eliminated" and "is a carrier" are different claims.
        # F is the fraction of the spread SURVIVING the term's maximum possible leverage.
        if ml["F"] >= 1.0:
            state = "ELIMINATED"          # explains nothing, or makes the spread worse
        elif ml["F"] > CARRIER_F:
            state = "NOT A CARRIER"       # explains <20% even at max leverage
        else:
            state = "SURVIVES"            # the ladder must measure it
        ml["spearman_with_day14"] = rho
        ml["state"] = state
        ml["explains_frac_at_max_leverage"] = 1.0 - ml["F"]
        out[name] = ml

        print(f"\n--- {label} ---")
        print(f"   spread {term.min():.4f} .. {term.max():.4f}   ({term.max()/term.min():.2f}x)")
        print(f"   Spearman(term, day-14) = {rho:+.4f}")
        print(render_table(
            ["fold", "term", "day-14", "MAX-leverage pred", "residual"],
            [[dn, f"{term[i]:.4f}", f"{d[i]:+.3f}", f"{ml['predicted'][i]:+.3f}",
              f"{ml['residual'][i]:+.3f}"] for i, dn in enumerate(DONORS)],
            aligns=["l", "r", "r", "r", "r"]))
        print(f"   spread(day-14) {ml['spread_d']:.3f} yr   "
              f"spread(residual) {ml['spread_residual']:.3f} yr   F = {ml['F']:.3f}")
        print(f"   worst miss {ml['worst_miss_yr']:.3f} yr   "
              f"explains {100 * ml['explains_frac_at_max_leverage']:+.1f}% of the spread "
              "at MAX leverage")
        print(f"   -> {state}")

    surv = [k for k in ("T1_mask", "T2_floor") if out[k]["state"] == "SURVIVES"]
    out["surviving_scalar_terms"] = surv
    out["verdict"] = ("NO_SCALAR_TERM_IS_A_CARRIER" if not surv
                      else "SURVIVES:" + ",".join(surv))

    print("\n" + "=" * 78)
    print(f"   VERDICT: {out['verdict']}")
    if not surv:
        print("   Neither scalar channel can produce the spread at any leverage fraction.")
        print("   §5.4's T2 > T3 prediction is NOT supported: the fold with the anomalous floor")
        print("   ratio (Y1) has normal labels, and the fold that collapses (N2) has a normal one.")
        print("   Within harmonization this leaves T3 (per-gene sigma_gill). Outside it, the")
        print("   deconfounder -- which A2's downgrade left unbounded -- is now in scope.")
    print("   SCOPE: scalars only. Floor effects acting through WHICH genes clamp, and mask")
    print("   effects acting through gene IDENTITY, are not scalar and are NOT tested here.")
    print("   This narrows the ladder; it does not close it.")

    OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
