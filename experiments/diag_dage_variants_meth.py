"""STAGE 1.5.6 step 2 — score every ΔAge variant against METHYLATION, in BOTH Gill arms.

    python experiments/diag_dage_variants_meth.py "D:\\Gill" "D:\\GSE165178" "D:\\GSE165177" "D:\\GSE165179"

READ-ONLY. Writes `results/diag_dage_variants_meth_results.json`. `src/` untouched, no labels move.

WHY, AND WHY NOT THE FACS OUTCOME
---------------------------------
Step 1 scored the same variants on predicting the FACS sort marker and **no variant helped** -- but
the baseline (day + pluripotency) was already at AUC 0.9221, leaving 0.078 of headroom, and SSEA4 is
itself a pluripotency surface marker. That test had a ceiling, so it could not separate "the label
adds nothing" from "there was nothing left to add".

Methylation has no such ceiling. It is a **different molecular layer**, it is the instrument 1.5.1
validated (negative control inert at +0.5/-2.4, dose-response p = 0.0001), and pluripotency does not
trivially predict it -- M-2a measured rho_partial 0.309 / 0.517 *after* partialling exactly that.

TWO ARMS, WHICH IS THE POINT
----------------------------
    Sendai      GSE165176 RNA  x  GSE165178 methylation   (join on donor_day_marker)
    transient   GSE165177 RNA  x  GSE165179 methylation   (join on exact title)

Independent experiments, different protocols, partly different donors. **A variant that wins in one
arm and not the other is fitting noise.** Consistency across the two is the finding; a single arm's
ranking is not.

THE ESTIMAND, AND WHY ABSOLUTE AGE IS THE RIGHT THING TO SCORE
--------------------------------------------------------------
rho_partial(variant, methylation | pluripotency, donor) -- 1.5.4's estimand, unchanged.

Scoring ABSOLUTE age rather than ΔAge is not a shortcut: ΔAge is `age - (a per-donor constant)`, and
partialling donor removes exactly that constant. So rho_partial(ΔAge, . | donor) is *identical* to
rho_partial(age, . | donor), and using absolute age avoids needing a control arm in the methylation
data -- which the Sendai methylation does not have (no day 0 among its 22 samples).

Both Horvath clocks are scored and both are always reported, whichever way they fall.
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


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def build_variants(norm: np.ndarray, genes: list[str], W: dict[str, float],
                   plu: np.ndarray, cc: np.ndarray, dv) -> dict[str, np.ndarray]:
    """The same nine definitions step 1 swept, computed identically."""
    wv = np.array([W.get(g, 0.0) for g in genes])
    ranked = np.argsort(-np.abs(wv))
    cov = dv.covered_weight_fraction(genes, W)
    out = {"raw": norm @ wv, "covnorm": (norm @ wv) * (1.0 / cov if cov else 1.0)}
    for k in (100, 500, 2000):
        sub = np.zeros_like(wv)
        sub[ranked[:k]] = wv[ranked[:k]]
        out[f"top{k}"] = norm @ sub
    out["ranknorm"] = dv.rank_normalise(norm) @ wv
    out["resid_pluri"] = dv.residualise(out["raw"], plu)
    out["resid_cc"] = dv.residualise(out["raw"], cc)
    out["resid_both"] = dv.residualise(out["raw"], np.column_stack([plu, cc]))
    return out


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
    dlc = _load("dlc", ROOT / "experiments" / "diag_learned_clock.py")
    dma = _load("dma", ROOT / "experiments" / "diag_methylation_anchor.py")
    m2a = _load("m2a", ROOT / "experiments" / "diag_m2a_calibratability.py")
    dcv = _load("dcv", ROOT / "experiments" / "diag_clock_validity.py")

    clock = json.loads((ROOT / "configs" / "clocks" / "fleischer_clock.json").read_text("utf-8"))
    W = {k: float(v) for k, v in clock["weights"].items()}

    def meth_lp(meth_dir: Path, matrix: str) -> dict:
        """Methylation linear predictor per sample, for both clocks.

        Deliberately does NOT call `dma.load_series`: that parser is hard-coded to GSE165179's
        characteristic names ("length of reprogramming phase (days)") and raises on GSE165178, which
        spells its fields differently. Nothing here needs that metadata — the sample names come from
        the beta matrix header and the join key is the title, so the series matrix is not consulted
        at all.
        """
        bpath = meth_dir / matrix
        out = {}
        for cfile, cname in dma.CLOCKS:
            cl = json.loads((ROOT / "configs" / "clocks" / f"{cfile}.json").read_text("utf-8"))
            Wm = {k: float(v) for k, v in cl["weights"].items()}
            samples, betas = dma.load_betas(bpath, set(Wm))
            out[cname] = {s: float(dma.linear_predictor(betas[s], Wm)[0]) for s in samples}
        return out

    arms = {}

    # ---- SENDAI: RNA title is `<donor>_<day>_<marker>_Sendai_ExpN`; meth title drops the suffix -- #
    samples, genes, lin = dv.load_gill(gill)
    norm = normalize_counts(lin)
    gi = {g: i for i, g in enumerate(genes)}
    plu = norm[:, [gi[g] for g in dcv.OSKM_PLURIPOTENCY if g in gi]].mean(axis=1)
    cc = norm[:, [gi[g] for g in dcv.CELL_CYCLE if g in gi]].mean(axis=1)
    var_sen = build_variants(norm, genes, W, plu, cc, dv)
    ms = meth_lp(meth_sen, next(meth_sen.glob("*Matrix_processed*.txt.gz")).name)
    key = [re.sub(r"_Sendai_Exp\d+$", "", s) for s in samples]
    anyc = next(iter(ms))
    idx = [i for i, k in enumerate(key) if k in ms[anyc]]
    arms["sendai"] = {"idx": np.array(idx), "keys": [key[i] for i in idx],
                      "var": var_sen, "plu": plu, "meth": ms,
                      "donor": np.array([samples[i].split("_")[0] for i in idx])}

    # ---- TRANSIENT: exact title match, then average exp replicates per (donor, arm, day) -------- #
    samples2, genes2, lin2 = m2a.load_rna(rna_tr)
    norm2 = normalize_counts(lin2)
    gi2 = {g: i for i, g in enumerate(genes2)}
    plu2 = norm2[:, [gi2[g] for g in dcv.OSKM_PLURIPOTENCY if g in gi2]].mean(axis=1)
    cc2 = norm2[:, [gi2[g] for g in dcv.CELL_CYCLE if g in gi2]].mean(axis=1)
    var_tr = build_variants(norm2, genes2, W, plu2, cc2, dv)
    mt = meth_lp(meth_tr, next(meth_tr.glob("*Matrix_processed*.txt.gz")).name)
    anyt = next(iter(mt))
    idx2 = [i for i, s in enumerate(samples2) if s in mt[anyt] and m2a.parse_title(s)]
    arms["transient"] = {"idx": np.array(idx2), "keys": [samples2[i] for i in idx2],
                         "var": var_tr, "plu": plu2, "meth": mt,
                         "donor": np.array([samples2[i].split("_")[0] for i in idx2])}

    out = {"script": "diag_dage_variants_meth", "utc": datetime.now(UTC).isoformat(), "arms": {}}
    for arm, D in arms.items():
        i = D["idx"]
        print(f"\n=== {arm.upper()}  n = {len(i)} joined samples, "
              f"{len(set(D['donor']))} donors ===")
        blk = {}
        for cname, table in D["meth"].items():
            y = np.array([table[k] for k in D["keys"]], float)
            print(f"\n  {cname}")
            print(f"  {'variant':>12} {'rho_partial':>12}")
            for name, v in D["var"].items():
                r = dlc.partial_spearman(v[i], y, D["plu"][i], D["donor"])
                blk.setdefault(name, {})[cname] = r
                print(f"  {name:>12} {r:12.3f}")
        rank = sorted(blk, key=lambda k: -min(blk[k].values()))
        out["arms"][arm] = {"n": len(i), "n_donors": len(set(D["donor"])),
                            "per_variant": blk, "ranked_by_worst_clock": rank}
        print(f"\n  ranked (by the WORSE clock, so a split cannot win): {', '.join(rank)}")

    a, b = out["arms"]["sendai"]["ranked_by_worst_clock"], out["arms"]["transient"]["ranked_by_worst_clock"]
    agree = a[0] == b[0]
    out["consistency"] = {"sendai_best": a[0], "transient_best": b[0], "same_winner": agree,
                          "spearman_of_rankings": float(dlc.spearman(
                              [a.index(v) for v in a], [b.index(v) for v in a]))}
    print(f"\n=== CONSISTENCY ===\n  sendai best: {a[0]}   transient best: {b[0]}   "
          f"same winner: {agree}")
    print(f"  rank agreement across arms (Spearman): "
          f"{out['consistency']['spearman_of_rankings']:+.3f}")
    (_RESULTS / "diag_dage_variants_meth_results.json").write_text(json.dumps(out, indent=2), "utf-8")
    print("  wrote results/diag_dage_variants_meth_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
