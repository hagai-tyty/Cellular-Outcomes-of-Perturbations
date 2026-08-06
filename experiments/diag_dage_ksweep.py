"""STAGE 1.5.6 step 3 — how sparse should the clock be? A k-sweep, scored in BOTH arms.

    python experiments/diag_dage_ksweep.py "D:\\Gill" "D:\\GSE165178" "D:\\GSE165177" "D:\\GSE165179"

READ-ONLY. Writes `results/diag_dage_ksweep_results.json`. `src/` untouched, no labels move.

WHY
---
The ledger found something the correlation-only tests had hidden for the whole arc: restricting the
Fleischer clock to its largest-|weight| genes does not merely help, it removes a **systematic bias**.

    raw     MAE 16.61 yr   bias -14.10   rho +0.703   sign 0.62     (multi-tissue, transient)
    top100  MAE  5.36 yr   bias  -1.61   rho +0.835   sign 0.94

The dense clock over-reports rejuvenation by ~15 years. That is not noise -- a bias that large and
that consistent is the signature of thousands of near-zero weights each contributing a little drift,
which is exactly the p >> n regime the clock was fitted in (33,155 genes from 133 samples).

**rho_partial could never have seen this.** It partials out donor, and it is scale- and shift-free,
so a uniform -15 yr offset is invisible to it. MAE and bias against a real instrument are not.

WHAT THIS SWEEPS
----------------
k in {10, 20, 50, 100, 150, 200, 300, 500, 1000, 2000, 5000, all}: keep only the k largest-|weight|
genes, zero the rest, and score ΔAge against methylation ΔAge.

  transient   68 conditions, ΔAge vs ΔAge (GSE165179 has matched negative-control arms)
  Sendai      22 conditions, ABSOLUTE age only (GSE165178 has no untreated control at any day)

Replicates are averaged to one row per condition BEFORE scoring, in both modalities.

SCORED ON FOUR THINGS, NOT ONE
------------------------------
  MAE    years of disagreement with the reference instrument
  bias   mean signed error -- the quantity that exposed this, and the one a correlation cannot see
  rho    Spearman, for ordering
  sign   fraction of conditions where both modalities agree on the DIRECTION of the change

A k that wins on MAE while losing on sign is not a fix, so all four are reported for every k and
both clocks, and the choice is made on agreement across the two arms rather than on one number.
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

KS = (10, 20, 50, 100, 150, 200, 300, 500, 1000, 2000, 5000, None)   # None = all genes


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


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
    dl = _load("dl", ROOT / "experiments" / "diag_dage_ledger.py")
    dma = _load("dma", ROOT / "experiments" / "diag_methylation_anchor.py")
    m2a = _load("m2a", ROOT / "experiments" / "diag_m2a_calibratability.py")

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
            out[cname] = {s: dl.anti_trafo(float(dma.linear_predictor(betas[s], Wm)[0]))
                          for s in samples}
        return out

    res: dict = {"script": "diag_dage_ksweep", "utc": datetime.now(UTC).isoformat(), "arms": {}}

    # ---------------------------------------------------------------- transient ---- #
    s1, g1, l1 = m2a.load_rna(rna_tr)
    n1 = normalize_counts(l1)
    mt = meth_years(meth_tr)
    anyc = next(iter(mt))
    keep = [i for i, s in enumerate(s1) if s in mt[anyc] and m2a.parse_title(s)]
    keys = [(p["donor"], p["arm"], p["day"]) for i in keep if (p := m2a.parse_title(s1[i]))]
    uk, _, _ = dl.average_by(keys, np.zeros(len(keys)))
    truth = {}
    for cname, tbl in mt.items():
        _, mv, _ = dl.average_by(keys, np.array([tbl[s1[i]] for i in keep], float))
        truth[cname] = dl.delta_vs_control(uk, mv)

    wv1 = np.array([W.get(g, 0.0) for g in g1])
    ord1 = np.argsort(-np.abs(wv1))
    res["arms"]["transient"] = {"n_conditions": len(uk), "mode": "dAge vs dAge", "k": {}}
    for k in KS:
        sub = wv1.copy()
        if k is not None:
            sub = np.zeros_like(wv1)
            sub[ord1[:k]] = wv1[ord1[:k]]
        _, av, _ = dl.average_by(keys, (n1[keep] @ sub) + b0)
        ad = dl.delta_vs_control(uk, av)
        res["arms"]["transient"]["k"][str(k)] = {
            c: dl.score(truth[c], ad) for c in truth}

    # ------------------------------------------------------------------- sendai ---- #
    s2, g2, l2 = dv.load_gill(gill)
    n2 = normalize_counts(l2)
    ms = meth_years(meth_sen)
    anys = next(iter(ms))
    key2 = [re.sub(r"_Sendai_Exp\d+$", "", s) for s in s2]
    keep2 = [i for i, kk in enumerate(key2) if kk in ms[anys]]
    kk2 = [key2[i] for i in keep2]
    uk2, _, _ = dl.average_by(kk2, np.zeros(len(kk2)))
    truth2 = {}
    for cname, tbl in ms.items():
        _, mv, _ = dl.average_by(kk2, np.array([tbl[key2[i]] for i in keep2], float))
        truth2[cname] = list(mv)

    wv2 = np.array([W.get(g, 0.0) for g in g2])
    ord2 = np.argsort(-np.abs(wv2))
    res["arms"]["sendai"] = {"n_conditions": len(uk2), "mode": "ABSOLUTE age", "k": {}}
    for k in KS:
        sub = wv2.copy()
        if k is not None:
            sub = np.zeros_like(wv2)
            sub[ord2[:k]] = wv2[ord2[:k]]
        _, av, _ = dl.average_by(kk2, (n2[keep2] @ sub) + b0)
        res["arms"]["sendai"]["k"][str(k)] = {c: dl.score(truth2[c], list(av)) for c in truth2}

    for arm, blk in res["arms"].items():
        print(f"\n=== {arm.upper()}  ({blk['mode']}, n = {blk['n_conditions']}) ===")
        for cname in truth:
            tag = "skin&blood" if "skin" in cname else "multi-tissue"
            print(f"\n  {tag}")
            print(f"  {'k':>6} {'MAE':>9} {'bias':>9} {'rho':>8} {'sign':>7}")
            for k in KS:
                s = blk["k"][str(k)][cname]
                if s["mae"] is None:
                    continue
                sg = f"{s['sign_agree']:7.2f}" if s["sign_agree"] is not None else "      -"
                print(f"  {str(k):>6} {s['mae']:9.2f} {s['bias']:+9.2f} {s['rho']:+8.3f}{sg}")

    (_RESULTS / "diag_dage_ksweep_results.json").write_text(json.dumps(res, indent=2), "utf-8")
    print("\nwrote results/diag_dage_ksweep_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
