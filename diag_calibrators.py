"""Compare calibrator families HONESTLY on the cross-donor pool, without touching the graded folds.

Run 3 shipped a union-fitted logit-Platt and missed `fate_ece` (0.249 vs the <=0.169 bar). Four
questions were left open, and none of them needs a retrain:

  Q1 FAMILY    is `sigmoid(a*logit p + b)` the wrong shape? `scorecard._platt` uses
               `LogisticRegression` on RAW p, and a synthetic probe of the two disagreed with the
               measured numbers -- so this has to be settled on real data, not by simulation.
  Q2 FIT SET   the shipped fit was 97.7% in-distribution (103 of 4509 rows cross-donor). Does
               fitting on the pool alone actually beat it OUT of sample?
  Q3 CAPACITY  can ANY monotone map reach the bar here? Isotonic fitted in-sample is the ceiling:
               nothing that preserves ranking can do better, so if the ceiling is above 0.169 the
               problem is not the calibrator and no amount of family-shopping will fix it.
  Q4 EFFECTIVE-n  cells within a donor share that donor's offset. If the between-donor variance
               dominates, the pool's effective n is nearer 5 than 103 and a 2-parameter fit is
               being estimated from ~5 points -- which bounds how well ANY choice can generalise.

METHOD. Every out-of-sample number is leave-one-donor-out WITHIN the pool: fit on 4 donors,
score the 5th, pool the held-out predictions. That is the same shift the deployed calibrator
faces, and it never reads `test_*`. The in-sample column is printed beside it only to show the
size of the optimism -- never as a result.

    python diag_calibrators.py            # reads diag_dump/ written by dump_diag_bundle.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from cellfate.common.calibration import platt_safe                  # noqa: E402
from cellfate.common.constants import SAFE_IDX                      # noqa: E402
from cellfate.training.calibrate import fit_platt_binary            # noqa: E402

DUMP = Path("diag_dump")
BAR = 0.169          # STAGE_1_CALIBRATION.md section 3: >=40% drop from 0.281
FLOOR = 0.078        # measured ECE estimator floor at this geometry; see the run-3 notebook entry


def ece(p, y, bins: int = 10) -> float:
    """Verbatim the estimator scorecard.py grades with (scorecard.py:112). Do not 'improve' it --
    a different binning here would silently stop being the quantity under the bar."""
    p, y = np.asarray(p, float), np.asarray(y, float)
    edges = np.linspace(0, 1, bins + 1)
    e = 0.0
    for i in range(bins):
        hi = edges[i + 1] if i < bins - 1 else 1.0 + 1e-9
        m = (p >= edges[i]) & (p < hi)
        if m.sum():
            e += (m.sum() / len(p)) * abs(p[m].mean() - y[m].mean())
    return float(e)


# --------------------------------------------------------------- families ------- #
def _fit_identity(p, y):
    return lambda q: np.asarray(q, float)


def _fit_logit_platt(p, y):
    """What Stage 1 ships."""
    a, b = fit_platt_binary(p, y)
    return lambda q: platt_safe(q, a, b)


def _fit_logistic_on_p(p, y):
    """What scorecard._platt / every T8.2 number used: LogisticRegression on the RAW probability."""
    from sklearn.linear_model import LogisticRegression
    if not (0 < np.sum(y >= 0.5) < len(y)):
        return lambda q: np.asarray(q, float)
    lr = LogisticRegression(max_iter=1000).fit(np.asarray(p).reshape(-1, 1), (y >= 0.5).astype(int))
    return lambda q: lr.predict_proba(np.asarray(q, float).reshape(-1, 1))[:, 1]


def _fit_isotonic(p, y):
    """The capacity ceiling: the best any RANK-PRESERVING map can do."""
    from sklearn.isotonic import IsotonicRegression
    ir = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0).fit(p, y)
    return lambda q: ir.predict(np.asarray(q, float))


FAMILIES = {
    "identity (no calibration)": _fit_identity,
    "logit-Platt  [SHIPPED]": _fit_logit_platt,
    "logistic-on-p  [T8.2]": _fit_logistic_on_p,
    "isotonic": _fit_isotonic,
}


def lodo_scores(p, y, donor, fit) -> tuple[float, float]:
    """(out-of-sample ECE by leave-one-donor-out, in-sample ECE). The pair is the point:
    their gap IS the optimism, and a family that only wins in-sample has not won."""
    oof = np.empty_like(np.asarray(p, float))
    for d in np.unique(donor):
        te = donor == d
        if (~te).sum() < 3 or len(np.unique(y[~te] >= 0.5)) < 2:
            oof[te] = p[te]                       # cannot fit without both classes; pass through
            continue
        oof[te] = fit(p[~te], y[~te])(p[te])
    return ece(oof, y), ece(fit(p, y)(p), y)


def effective_n(p, y, donor) -> dict:
    """Q4. Split the calibration error into between-donor and within-donor variance.

    If a donor's cells share that donor's offset, the pool carries far less independent
    information than its row count suggests, and the usual n=103 intuition is wrong.
    """
    err = np.asarray(y, float) - np.asarray(p, float)
    groups = [err[donor == d] for d in np.unique(donor)]
    k = len(groups)
    means = np.array([g.mean() for g in groups])
    within = float(np.mean([g.var(ddof=0) for g in groups]))
    between = float(means.var(ddof=0))
    icc = between / (between + within) if (between + within) > 0 else float("nan")
    n_bar = float(np.mean([len(g) for g in groups]))
    # Kish: independent rows are discounted by the design effect 1 + (n_bar - 1) * ICC
    deff = 1.0 + (n_bar - 1.0) * icc if np.isfinite(icc) else float("nan")
    return {"donors": k, "n": int(len(err)), "between": between, "within": within,
            "icc": icc, "n_eff": len(err) / deff if deff and np.isfinite(deff) else float("nan")}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dump", default=str(DUMP))
    args = ap.parse_args()
    root = Path(args.dump)
    files = sorted(root.glob("*.npz"))
    if not files:
        print(f"No .npz in {root}/. Run `python dump_diag_bundle.py` first (or pass --dump).")
        return 1

    agg: dict[str, list[float]] = {k: [] for k in FAMILIES}
    agg_in: dict[str, list[float]] = {k: [] for k in FAMILIES}
    eff_rows, sat_rows = [], []

    for f in files:
        d = np.load(f, allow_pickle=False)
        if "pool_probs_mean" not in d or "pool_donor_id" not in d:
            print(f"  [!] {f.stem}: no pool or no donor labels -- skipped")
            continue
        p = d["pool_probs_mean"][:, SAFE_IDX]
        y = d["pool_targets"][:, SAFE_IDX]
        donor = d["pool_donor_id"]

        print(f"\n  === fold {f.stem} ===  pool n={len(p)}  donors={len(np.unique(donor))}")
        print(f"  {'family':<28}{'ECE (LODO)':>12}{'ECE (in-sample)':>18}{'optimism':>11}")
        print("  " + "-" * 67)
        for name, fit in FAMILIES.items():
            try:
                oos, ins = lodo_scores(p, y, donor, fit)
            except Exception as exc:  # noqa: BLE001
                print(f"  {name:<28}{'failed':>12}   {exc!r}"[:100])
                continue
            agg[name].append(oos)
            agg_in[name].append(ins)
            flag = "  <-- ceiling" if name == "isotonic" else ""
            print(f"  {name:<28}{oos:>12.3f}{ins:>18.3f}{oos - ins:>11.3f}{flag}")

        eff = effective_n(p, y, donor)
        eff_rows.append(eff)
        print(f"  effective n: {eff['n']} rows over {eff['donors']} donors -> "
              f"ICC={eff['icc']:.3f}, n_eff={eff['n_eff']:.1f}")
        sat_rows.append(float(np.mean(p > 0.99)))

    if not any(agg.values()):
        print("\nNothing scored. The dump has no pool rows with donor labels.")
        return 1

    print(f"\n\n  MEAN OVER {len(agg['isotonic'])} FOLDS   (bar {BAR:.3f}; "
          f"estimator floor {FLOOR:.3f})\n")
    print(f"  {'family':<28}{'ECE (LODO)':>12}{'vs bar':>10}")
    print("  " + "-" * 50)
    for name in FAMILIES:
        if not agg[name]:
            continue
        m = float(np.mean(agg[name]))
        print(f"  {name:<28}{m:>12.3f}{('PASS' if m <= BAR else 'miss'):>10}")

    # Q3 must be read OUT of sample. Isotonic in-sample is ~0 by construction -- it interpolates
    # the fitting points -- so it measures flexibility, never attainability. The informative
    # quantity is how the most flexible monotone family does on donors it did not see.
    iso = float(np.mean(agg["isotonic"])) if agg["isotonic"] else float("nan")
    iso_in = float(np.mean(agg_in["isotonic"])) if agg_in["isotonic"] else float("nan")
    best = min((float(np.mean(v)) for v in agg.values() if v), default=float("nan"))
    print(f"\n  CAPACITY (Q3): best out-of-sample over all families = {best:.3f}; "
          f"isotonic LODO = {iso:.3f}.")
    print(f"  (isotonic in-sample is {iso_in:.3f} -- it interpolates, so it bounds nothing.)")
    if best <= BAR:
        print(f"    -> a monotone calibrator DOES reach the bar ({BAR:.3f}) out of sample on this")
        print("       pool. The miss is a choice of family/fitting-set, not a capacity limit.")
    else:
        print(f"    -> NO family tried reaches {BAR:.3f} out of sample, and isotonic is the most")
        print("       flexible monotone map there is. The miss is then NOT a calibrator choice --")
        print("       re-calibrating cannot fix it and the next change has to be elsewhere.")

    if eff_rows:
        icc = float(np.nanmean([e["icc"] for e in eff_rows]))
        neff = float(np.nanmean([e["n_eff"] for e in eff_rows]))
        print(f"\n  EFFECTIVE n (Q4): mean ICC {icc:.3f}, n_eff {neff:.1f} of "
              f"{eff_rows[0]['n']} rows.")
        print(f"    -> {'donor offsets dominate; the fit sees ~%.0f independent points, not %d'% (neff, eff_rows[0]['n']) if icc > 0.3 else 'rows are close to independent; n_eff is near the row count'}")
    if sat_rows:
        print(f"\n  SATURATION: mean {np.mean(sat_rows):.1%} of pool rows have P(safe) > 0.99.")

    print("\n  Every ECE (LODO) above is out-of-sample by leave-one-donor-out WITHIN the pool.")
    print("  No number here was computed from test_* -- the graded folds are untouched.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
