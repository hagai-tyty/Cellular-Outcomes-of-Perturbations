"""STAGE 1.5.6 — the ΔAge LEDGER: truth vs expectation vs actual, per condition, per variant.

    python experiments/diag_dage_ledger.py "D:\\Gill" "D:\\GSE165178" "D:\\GSE165177" "D:\\GSE165179"

READ-ONLY. Writes `results/DAGE_LEDGER.md`, `results/diag_dage_ledger_results.json` and
`results/dage_ledger.csv`. `src/` untouched, no labels move.

WHAT THIS RECORDS, AND WHY IN THIS FORM
---------------------------------------
Every earlier test in this arc reported a single summary statistic per variant, which is why a
negative result never said WHERE it broke. This records, for every condition and every variant:

    TRUTH     the methylation ΔAge in years  -- the instrument 1.5.1 validated (negative control
              inert at +0.5/-2.4 yr, dose-response p = 0.0001, SNR 3.4)
    EXPECTED  what the pipeline SHOULD produce if the RNA clock worked -- i.e. the truth
    ACTUAL    what the RNA clock variant actually produced
    ERROR     actual - truth, in years

so a failure can be read per cell rather than inferred from a correlation.

THE REPLICATE FIX
-----------------
The previous run scored 30 Sendai rows drawn from 22 independent methylation samples, and 90
transient rows from 68 conditions -- exp1/exp2 replicates were left unaveraged, so the same
methylation value backed two rows. That is pseudo-replication: it inflates apparent significance
without adding information. **Here every modality is averaged to one row per (donor, arm, day)
condition BEFORE anything is scored**, which is the unit-of-analysis rule 1.5.1 established.

ΔAge, AND WHERE IT CANNOT BE FORMED
-----------------------------------
ΔAge = age(condition) - age(matched control), matched WITHIN (donor, day).

  transient  GSE165179 carries `negative_control` and `negative_control_intermediate` arms, so a
             matched control exists at every day and ΔAge is computable on both modalities.
  Sendai     GSE165178's 22 samples are days 9/11/15, markers CD13/SSEA4 -- **no untreated control
             at any day**. ΔAge cannot be formed on the methylation side, so that arm is recorded
             on ABSOLUTE age and labelled as such rather than quietly given a substitute baseline.

Both Horvath clocks are always reported, whichever way they fall.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
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

ADULT_AGE = 20.0
CTRL_ARMS = ("negative_control", "negative_control_intermediate")


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- #
# Pure logic — unit-tested with no repo data present                           #
# --------------------------------------------------------------------------- #
def anti_trafo(x: float, adult_age: float = ADULT_AGE) -> float:
    """Horvath's inverse transform: linear above `adult_age`, exponential below."""
    return (1 + adult_age) * math.exp(x) - 1 if x < 0 else (1 + adult_age) * x + adult_age


def average_by(keys: list, values: np.ndarray) -> tuple[list, np.ndarray, list[int]]:
    """Collapse rows sharing a key to their mean. Returns (unique_keys, means, n_replicates).

    This is the replicate fix: exp1/exp2 are repeats of ONE condition, so scoring them as two
    independent rows is pseudo-replication.
    """
    order: dict = {}
    for i, k in enumerate(keys):
        order.setdefault(k, []).append(i)
    uk = sorted(order)
    return uk, np.array([values[order[k]].mean() for k in uk]), [len(order[k]) for k in uk]


def delta_vs_control(keys: list[tuple], vals: np.ndarray, ctrl_arms=CTRL_ARMS):
    """ΔAge against the matched control at the SAME (donor, day). None where no control exists.

    Keys are (donor, arm, day). Matching within donor AND day is what makes this a controlled
    contrast rather than a time trend.
    """
    ctrl: dict = {}
    for (d, a, day), v in zip(keys, vals, strict=True):
        if a in ctrl_arms:
            ctrl.setdefault((d, day), []).append(v)
    out = []
    for (d, _a, day), v in zip(keys, vals, strict=True):
        c = ctrl.get((d, day))
        out.append(v - float(np.mean(c)) if c else None)
    return out


def score(truth: list, actual: list) -> dict:
    """MAE, bias, Spearman and sign agreement over the conditions where truth exists."""
    pairs = [(t, a) for t, a in zip(truth, actual, strict=True)
             if t is not None and a is not None and np.isfinite(t) and np.isfinite(a)]
    if len(pairs) < 3:
        return {"n": len(pairs), "mae": None, "bias": None, "rho": None, "sign_agree": None}
    t = np.array([p[0] for p in pairs], float)
    a = np.array([p[1] for p in pairs], float)
    nz = (t != 0) | (a != 0)
    return {
        "n": len(pairs),
        "mae": float(np.mean(np.abs(a - t))),
        "bias": float(np.mean(a - t)),
        "rho": _spear(a, t),
        "sign_agree": float(np.mean(np.sign(a[nz]) == np.sign(t[nz]))) if nz.any() else None,
    }


def _spear(x: np.ndarray, y: np.ndarray) -> float:
    def rk(v):
        o = np.argsort(v, kind="mergesort")
        r = np.empty(len(v), float)
        r[o] = np.arange(1, len(v) + 1)
        for u in np.unique(v):
            m = v == u
            if m.sum() > 1:
                r[m] = r[m].mean()
        return r
    a, b = rk(x) - rk(x).mean(), rk(y) - rk(y).mean()
    den = math.sqrt(float((a**2).sum() * (b**2).sum()))
    return float(a @ b / den) if den else float("nan")


# --------------------------------------------------------------------------- #
# Real-data wiring                                                             #
# --------------------------------------------------------------------------- #
def main() -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    if len(sys.argv) < 5:
        print(__doc__)
        return 2
    gill, meth_sen, rna_tr, meth_tr = (Path(a) for a in sys.argv[1:5])

    from cellfate.data.normalize import normalize_counts
    dv = _load("dv", ROOT / "experiments" / "diag_dage_variants.py")
    dvm = _load("dvm", ROOT / "experiments" / "diag_dage_variants_meth.py")
    dma = _load("dma", ROOT / "experiments" / "diag_methylation_anchor.py")
    m2a = _load("m2a", ROOT / "experiments" / "diag_m2a_calibratability.py")
    dcv = _load("dcv", ROOT / "experiments" / "diag_clock_validity.py")

    clock = json.loads((ROOT / "configs" / "clocks" / "fleischer_clock.json").read_text("utf-8"))
    W = {k: float(v) for k, v in clock["weights"].items()}
    b0 = float(clock.get("intercept", 0.0))

    def meth_years(meth_dir: Path) -> dict:
        bpath = next(meth_dir.glob("*Matrix_processed*.txt.gz"))
        out = {}
        for cfile, cname in dma.CLOCKS:
            cl = json.loads((ROOT / "configs" / "clocks" / f"{cfile}.json").read_text("utf-8"))
            Wm = {k: float(v) for k, v in cl["weights"].items()}
            samples, betas = dma.load_betas(bpath, set(Wm))
            out[cname] = {s: anti_trafo(float(dma.linear_predictor(betas[s], Wm)[0]))
                          for s in samples}
        return out

    ledger: list[dict] = []
    scores: dict = {}

    # ------------------------------------------------------------------ TRANSIENT ------ #
    samples, genes, lin = m2a.load_rna(rna_tr)
    norm = normalize_counts(lin)
    gi = {g: i for i, g in enumerate(genes)}
    plu = norm[:, [gi[g] for g in dcv.OSKM_PLURIPOTENCY if g in gi]].mean(axis=1)
    cc = norm[:, [gi[g] for g in dcv.CELL_CYCLE if g in gi]].mean(axis=1)
    variants = dvm.build_variants(norm, genes, W, plu, cc, dv)
    mt = meth_years(meth_tr)
    anyc = next(iter(mt))

    keep = [i for i, s in enumerate(samples) if s in mt[anyc] and m2a.parse_title(s)]
    keys = [(p["donor"], p["arm"], p["day"]) for i in keep
            if (p := m2a.parse_title(samples[i]))]
    uk, _, nrep = average_by(keys, np.zeros(len(keys)))

    meth_by_clock = {}
    for cname, table in mt.items():
        _, mv, _ = average_by(keys, np.array([table[samples[i]] for i in keep], float))
        meth_by_clock[cname] = (mv, delta_vs_control(uk, mv))

    rna_by_variant = {}
    for name, v in variants.items():
        _, av, _ = average_by(keys, (v[keep] + b0) if name not in ("resid_pluri", "resid_cc", "resid_both")
                              else v[keep])
        rna_by_variant[name] = (av, delta_vs_control(uk, av))

    for j, (d, a, day) in enumerate(uk):
        row = {"arm_set": "transient", "donor": d, "condition": a, "day": day,
               "n_replicates": nrep[j], "is_control": a in CTRL_ARMS}
        for cname, (mv, md) in meth_by_clock.items():
            tag = "sb" if "skin" in cname else "mt"
            row[f"TRUTH_meth_age_{tag}"] = round(float(mv[j]), 3)
            row[f"TRUTH_meth_dage_{tag}"] = (None if md[j] is None else round(float(md[j]), 3))
        for name, (av, ad) in rna_by_variant.items():
            row[f"ACTUAL_rna_age_{name}"] = round(float(av[j]), 3)
            row[f"ACTUAL_rna_dage_{name}"] = (None if ad[j] is None else round(float(ad[j]), 3))
            for cname, (_mv, md) in meth_by_clock.items():
                tag = "sb" if "skin" in cname else "mt"
                row[f"ERROR_{name}_{tag}"] = (
                    None if (ad[j] is None or md[j] is None) else round(float(ad[j] - md[j]), 3))
        ledger.append(row)

    for name, (_av, ad) in rna_by_variant.items():
        for cname, (_mv, md) in meth_by_clock.items():
            tag = "sb" if "skin" in cname else "mt"
            scores.setdefault("transient", {}).setdefault(name, {})[tag] = score(md, ad)

    # ------------------------------------------------------------------- SENDAI -------- #
    s2, g2, l2 = dv.load_gill(gill)
    n2 = normalize_counts(l2)
    gi2 = {g: i for i, g in enumerate(g2)}
    plu2 = n2[:, [gi2[g] for g in dcv.OSKM_PLURIPOTENCY if g in gi2]].mean(axis=1)
    cc2 = n2[:, [gi2[g] for g in dcv.CELL_CYCLE if g in gi2]].mean(axis=1)
    var2 = dvm.build_variants(n2, g2, W, plu2, cc2, dv)
    ms = meth_years(meth_sen)
    anys = next(iter(ms))
    key2 = [re.sub(r"_Sendai_Exp\d+$", "", s) for s in s2]
    keep2 = [i for i, k in enumerate(key2) if k in ms[anys]]
    kk = [key2[i] for i in keep2]
    uk2, _, nrep2 = average_by(kk, np.zeros(len(kk)))

    meth2 = {}
    for cname, table in ms.items():
        _, mv, _ = average_by(kk, np.array([table[key2[i]] for i in keep2], float))
        meth2[cname] = mv
    rna2 = {}
    for name, v in var2.items():
        _, av, _ = average_by(kk, (v[keep2] + b0) if name not in
                              ("resid_pluri", "resid_cc", "resid_both") else v[keep2])
        rna2[name] = av

    for j, k in enumerate(uk2):
        parts = k.split("_")
        row = {"arm_set": "sendai", "donor": parts[0], "condition": parts[-1],
               "day": parts[1], "n_replicates": nrep2[j], "is_control": False}
        for cname, mv in meth2.items():
            tag = "sb" if "skin" in cname else "mt"
            row[f"TRUTH_meth_age_{tag}"] = round(float(mv[j]), 3)
            row[f"TRUTH_meth_dage_{tag}"] = None      # no untreated control in GSE165178
        for name, av in rna2.items():
            row[f"ACTUAL_rna_age_{name}"] = round(float(av[j]), 3)
            row[f"ACTUAL_rna_dage_{name}"] = None
            for cname, mv in meth2.items():
                tag = "sb" if "skin" in cname else "mt"
                row[f"ERROR_{name}_{tag}"] = round(float(av[j] - mv[j]), 3)
        ledger.append(row)
    for name, av in rna2.items():
        for cname, mv in meth2.items():
            tag = "sb" if "skin" in cname else "mt"
            scores.setdefault("sendai_ABSOLUTE", {}).setdefault(name, {})[tag] = \
                score(list(mv), list(av))

    out = {"script": "diag_dage_ledger", "utc": datetime.now(UTC).isoformat(),
           "n_conditions": {"transient": len(uk), "sendai": len(uk2)},
           "note_sendai": "ABSOLUTE age only - GSE165178 has no untreated control at any day",
           "scores": scores}
    (_RESULTS / "diag_dage_ledger_results.json").write_text(json.dumps(out, indent=2), "utf-8")
    cols = sorted({k for r in ledger for k in r})
    with (_RESULTS / "dage_ledger.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(ledger)

    print(f"\nconditions: transient {len(uk)} (from {len(keep)} samples), "
          f"sendai {len(uk2)} (from {len(keep2)} samples)")
    for armset, per in scores.items():
        print(f"\n=== {armset} ===")
        print(f"  {'variant':>12} {'clock':>4} {'n':>4} {'MAE':>9} {'bias':>9} {'rho':>7} {'sign':>6}")
        for name, byc in per.items():
            for tag, s in byc.items():
                if s["mae"] is None:
                    continue
                print(f"  {name:>12} {tag:>4} {s['n']:4d} {s['mae']:9.2f} {s['bias']:+9.2f} "
                      f"{s['rho']:+7.3f} {s['sign_agree']:6.2f}"
                      if s["sign_agree"] is not None else
                      f"  {name:>12} {tag:>4} {s['n']:4d} {s['mae']:9.2f} {s['bias']:+9.2f} "
                      f"{s['rho']:+7.3f}      -")
    print("\nwrote results/dage_ledger.csv, results/diag_dage_ledger_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
