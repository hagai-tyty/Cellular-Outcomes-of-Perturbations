"""STAGE 1.5.2 §10 step 1 + M-2b — do RNA and methylation agree on the SSEA4 − CD13 CONTRAST?

    python experiments/diag_m2b_contrast_agreement.py                 # §10 step 1 only (shape)
    python experiments/diag_m2b_contrast_agreement.py --run           # + M-2b

READ-ONLY. Writes `diag_m2b_contrast_agreement_results.json`. `src/` untouched.

WHY THIS SERIES AND NOT THE OTHER ONE
-------------------------------------
M-2a ran on GSE165177 × GSE165179 because that pairing has 68 conditions and a real untreated
`negative_control` arm — the better instrument for the *verdict*. But **M-2b is defined on a
geometry only GSE165178 × GSE165176 has**: §5's `Δ = mean(SSEA4) − mean(CD13)` per (donor, day).
The transient series has no SSEA4/CD13 sort at all. GSE165178 is also, per §3, "the only series that
joins the samples the model actually trains on", so it is the one that can speak to the *labels*.

§10 STEP 1 — SHAPE BEFORE STATISTIC, WITH AN ABORT
--------------------------------------------------
"Load GSE165178. Print and record: sample count, donor roster, day grid, arm counts, probe count,
CpG coverage per clock, and the join result against GSE165176. **Abort if the join is not 22/22.**
No statistic is computed in this step." §9-R5 is why: the 22/22 join was verified on GEO *titles*,
and titles can mislead. This checks it against the actual downloaded matrices, and refuses to
proceed on a partial join rather than quietly analysing whatever matched.

THE BAR
-------
Sign agreement, **≥ 7/11** — moved from the registered ≥ 8/11 by the §6 freeze on 2026-07-31,
because 8/11 was UNRESOLVABLE (93.0%) and §5b requires an unresolvable bar to move *before* the run.
§6's own honest note applies and is repeated at the point of use: **this is a weaker test than
registered**, and any pass must be reported with that attached.

ρ(Δ_rna, Δ_meth) is reported alongside but is NOT a bar — at 11 pairs it was never registered as
one, and inventing one now would be choosing a criterion after seeing the geometry.

PRE-COMMITTED EXPECTATION (§5, restated so it cannot be revised afterwards)
--------------------------------------------------------------------------
"Methylation says reprogramming cells are YOUNGER (1.5.1: −24 to −27 yr). RNA has been reading
treated cells OLDER (+36.5 yr on non-responders). **Disagreement is the live hypothesis, and it is a
decisive result, not a failure of the stage.**"
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import re
import sys
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)


ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

ANNOT_COLS = 12
EXPECTED_JOIN = 22          # §3, verified on titles by REV FINAL §8.3
SIGN_BAR = 7                # §6 freeze 2026-07-31: moved from 8 to its usable_bar
COVERAGE_BAR = 0.90         # §9-R4
ADULT_SLOPE = 21.0

SSEA4, CD13 = "SSEA4", "CD13"


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


# --------------------------------------------------------------------------- #
# Pure logic — data-free, unit-tested                                          #
# --------------------------------------------------------------------------- #
def join_key(title: str) -> str | None:
    """`{donor}_{day}_{marker}` — the shared key. Pure.

    GSE165178 titles ARE the key (`Y2_d11_SSEA4`); GSE165176 titles are the key plus a batch suffix
    (`Y2_d11_SSEA4_Sendai_Exp1`). Anything not matching the shape returns None rather than being
    silently coerced, so a format surprise shows up as an unmatched sample instead of a wrong join.
    """
    m = re.match(r"^([A-Za-z]\d)_d(\d+)_(SSEA4|CD13)(?:_.*)?$", title.strip())
    return f"{m.group(1)}_d{m.group(2)}_{m.group(3)}" if m else None


def parse_key(key: str) -> dict:
    d, day, marker = key.split("_")
    return {"donor": d, "day": int(day[1:]), "marker": marker}


def verify_join(meth_titles: list[str], rna_titles: list[str], expected: int = EXPECTED_JOIN) -> dict:
    """§10 step 1's abort condition. Pure — takes titles, not files."""
    mk = {}
    for t in meth_titles:
        if (k := join_key(t)):
            mk.setdefault(k, []).append(t)
    rk = {}
    for t in rna_titles:
        if (k := join_key(t)):
            rk.setdefault(k, []).append(t)
    matched = sorted(set(mk) & set(rk))
    unmatched_meth = sorted(set(mk) - set(rk))
    info = [parse_key(k) for k in matched]
    return {
        "n_meth_samples": len(meth_titles), "n_meth_keyed": len(mk),
        "n_rna_samples": len(rna_titles), "n_rna_keyed": len(rk),
        "n_matched": len(matched), "expected": expected,
        "matched_keys": matched, "unmatched_meth": unmatched_meth,
        "donors": sorted({i["donor"] for i in info}),
        "days": sorted({i["day"] for i in info}),
        "arm_counts": dict(Counter(i["marker"] for i in info)),
        "rna_replicates_per_key": {k: len(v) for k, v in sorted(rk.items()) if k in set(matched)},
        "verdict": "OK" if len(matched) == expected and not unmatched_meth else "ABORT",
    }


def pair_contrasts(values: dict[str, float], keys: list[str]) -> list[dict]:
    """Δ = mean(SSEA4) − mean(CD13) per (donor, day). Pure.

    Replicates are AVERAGED, not treated as independent — 1.5.1's unit-of-analysis rule, and the fix
    for the pairing bug that once silently dropped 6 of 9 pairs.
    """
    grid: dict[tuple[str, int], dict[str, list[float]]] = {}
    for k in keys:
        p = parse_key(k)
        if k in values:
            grid.setdefault((p["donor"], p["day"]), {}).setdefault(p["marker"], []).append(values[k])
    out = []
    for (donor, day), arms in sorted(grid.items()):
        if SSEA4 in arms and CD13 in arms:
            t, c = float(np.mean(arms[SSEA4])), float(np.mean(arms[CD13]))
            out.append({"donor": donor, "day": day, "ssea4": t, "cd13": c, "delta": t - c})
    return out


def sign_agreement(d_rna: list[float], d_meth: list[float]) -> dict:
    """How many pairs move the same direction in both modalities? Pure."""
    pairs = [(a, b) for a, b in zip(d_rna, d_meth, strict=True)
             if np.isfinite(a) and np.isfinite(b)]
    agree = sum(1 for a, b in pairs if (a > 0) == (b > 0))
    return {"n_pairs": len(pairs), "n_agree": agree,
            "frac": agree / len(pairs) if pairs else float("nan")}


def agreement_by_day(rows: list[dict]) -> dict:
    """Sign agreement split by reprogramming day. Pure. DESCRIPTIVE — never a criterion.

    Added after seeing the first M-2b output, and labelled as such. It is here because the
    headline 7/11 is not distributed evenly across the grid, and a single pooled number hides
    which timepoints carry it. §4's confound applies to M-2b exactly as it does to M-2a: both
    modalities move with reprogramming progress, so agreement concentrated at the timepoints
    where BOTH effects are large is agreement about the day axis, not about age.
    """
    by_day: dict[int, dict] = {}
    for r in rows:
        d = by_day.setdefault(r["day"], {"n": 0, "agree": 0, "rna": [], "meth": []})
        d["n"] += 1
        d["agree"] += int(r["same_sign"])
        d["rna"].append(r["delta_rna_years"])
        d["meth"].append(r["delta_meth_years"])
    return {str(k): {"n": v["n"], "n_agree": v["agree"],
                     "mean_rna": float(np.mean(v["rna"])),
                     "mean_meth": float(np.mean(v["meth"]))}
            for k, v in sorted(by_day.items())}


def m2b_verdict(sign: dict, rho: float, bar: int = SIGN_BAR) -> dict:
    """What M-2b licenses. Pure. §7 reads this together with M-2a, never alone (§9-R1)."""
    if sign["n_pairs"] < 2:
        return {"status": "CANNOT_VERIFY", "reason": f"{sign['n_pairs']} usable pairs"}
    ok = sign["n_agree"] >= bar
    # Sitting exactly ON the bar is a FRAGILE pass: one pair flipping fails it. This project has
    # been burned three times by hairline margins (E1b by 0.009, D2 by 0.014, M-2a by 0.016), so
    # it is named in the status rather than left for the reader to notice.
    fragile = ok and sign["n_agree"] == bar
    return {
        "status": ("AGREE_FRAGILE" if fragile else "AGREE") if ok else "DISAGREE",
        "fragile": fragile,
        "n_agree": sign["n_agree"], "n_pairs": sign["n_pairs"], "bar": bar, "rho": rho,
        "reason": (f"{sign['n_agree']}/{sign['n_pairs']} pairs agree in sign against a bar of "
                   f"{bar}/{sign['n_pairs']} -> "
                   + ("the two modalities agree on the direction of the contrast"
                      if ok else
                      "the two modalities DISAGREE on the direction of the contrast, which §5 "
                      "pre-committed as the live hypothesis and a decisive result")),
        "caveat": " ".join(x for x in [
            ("bar moved 8/11 -> 7/11 by the §6 freeze because 8/11 was UNRESOLVABLE (93.0%). "
             "This is a WEAKER test than registered." if ok else ""),
            ("It also lands EXACTLY on the moved bar -- one pair flipping would fail it."
             if fragile else ""),
        ] if x),
    }


# --------------------------------------------------------------------------- #
# Data wiring                                                                  #
# --------------------------------------------------------------------------- #
def titles_of(series_matrix: Path) -> list[str]:
    with gzip.open(series_matrix, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                return [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
    return []


def rna_ages(gill_dir: Path) -> dict[str, float]:
    """Clock ages per GSE165176 sample, keyed by join key (replicates kept separate)."""
    import pandas as pd

    from cellfate.data.aging import LinearClock
    from cellfate.data.normalize import normalize_counts

    expr = next(gill_dir.glob("*Log2_RPM_Sendai*.txt.gz"))
    df = pd.read_csv(expr, sep="\t", low_memory=False)
    cols = list(df.columns)
    mat = df.set_index(cols[0])[cols[ANNOT_COLS:]]
    lin = np.power(2.0, mat.to_numpy(dtype=np.float64)) - 1.0
    lin[lin < 0] = 0.0
    sym = [str(s) for s in mat.index]
    order = np.argsort(-lin.sum(axis=1))
    seen, keep = set(), []
    for i in order:
        if sym[i] not in seen:
            seen.add(sym[i])
            keep.append(i)
    keep.sort()
    genes = [sym[i] for i in keep]
    norm = normalize_counts(lin[keep, :].T)
    clock = LinearClock.from_json(ROOT / "configs" / "clocks" / "fleischer_clock.json")
    ages = clock.predict_age(norm, genes)
    out: dict[str, list[float]] = {}
    for t, a in zip(mat.columns, ages, strict=True):
        if (k := join_key(str(t))):
            out.setdefault(k, []).append(float(a))
    return {k: float(np.mean(v)) for k, v in out.items()}


def meth_linpreds(meth_dir: Path) -> dict[str, dict]:
    """Linear predictor per GSE165178 sample, per clock, keyed by join key. Plus coverage.

    Returns the LINEAR PREDICTOR: `anti_trafo` is linear above adult_age, so a DIFFERENCE of two
    samples is `21 x (lp_t - lp_c)` with the intercept cancelling exactly (REV FINAL §4.3). M-2b is
    entirely differences, so no intercept is needed anywhere here.
    """
    dma = _load("diag_methylation_anchor", "experiments/diag_methylation_anchor.py")
    bpath = meth_dir / "GSE165178_Matrix_processed_sendai.txt.gz"
    out = {}
    for cfile, cname in dma.CLOCKS:
        W = {k: float(v) for k, v in json.loads(
            (ROOT / "configs" / "clocks" / f"{cfile}.json").read_text(encoding="utf-8")
        )["weights"].items()}
        samples, betas = dma.load_betas(bpath, set(W))
        vals, cov = {}, []
        for s in samples:
            if (k := join_key(str(s))):
                v, n = dma.linear_predictor(betas[s], W)
                vals.setdefault(k, []).append(v)
                cov.append(n)
        out[cname] = {"n_cpg": len(W),
                      "coverage_frac": (float(np.mean(cov)) / len(W)) if cov else 0.0,
                      "lp": {k: float(np.mean(v)) for k, v in vals.items()}}
    return out


def main() -> int:
    pos = [a for a in sys.argv[1:] if not a.startswith("--")]
    meth_dir = Path(pos[0] if pos else r"D:\GSE165178")
    gill_dir = Path(pos[1] if len(pos) > 1 else r"D:\Gill")

    print("STAGE 1.5.2 §10 step 1 — SHAPE BEFORE STATISTIC (read-only)\n")
    mt = titles_of(meth_dir / "GSE165178_series_matrix.txt.gz")
    rt = titles_of(gill_dir / "GSE165176_series_matrix.txt.gz")
    j = verify_join(mt, rt)
    print(f"  GSE165178 (methylation): {j['n_meth_samples']} samples, {j['n_meth_keyed']} keyed")
    print(f"  GSE165176 (RNA, the training series): {j['n_rna_samples']} samples, "
          f"{j['n_rna_keyed']} keyed")
    print(f"  JOIN on donor_day_marker: {j['n_matched']}/{j['expected']}  -> {j['verdict']}")
    print(f"  donors {j['donors']}   days {j['days']}   arms {j['arm_counts']}")
    reps = Counter(j["rna_replicates_per_key"].values())
    print(f"  RNA replicates per matched key: {dict(reps)}")
    if j["unmatched_meth"]:
        print(f"  [!] unmatched methylation samples: {j['unmatched_meth']}")

    out = {"script": "diag_m2b_contrast_agreement",
           "utc": datetime.now(UTC).isoformat(timespec="seconds"), "join": j}

    if j["verdict"] != "OK":
        print("\n  ==> ABORT (§10 step 1). The join is not 22/22, so no statistic is computed.")
        (_RESULTS / "diag_m2b_contrast_agreement_results.json").write_text(
            json.dumps(out, indent=2, default=str), encoding="utf-8")
        return 1

    print("\n  reading methylation matrix (clock CpGs only)...")
    meth = meth_linpreds(meth_dir)
    for cname, blk in meth.items():
        ok = blk["coverage_frac"] >= COVERAGE_BAR
        blk["coverage_verdict"] = "OK" if ok else "DEGRADED"
        print(f"     {cname:<28} {blk['n_cpg']} CpGs, coverage {blk['coverage_frac']:.1%} "
              f"[{blk['coverage_verdict']}]")
    out["methylation_coverage"] = {k: {"n_cpg": v["n_cpg"], "coverage_frac": v["coverage_frac"],
                                       "coverage_verdict": v["coverage_verdict"]}
                                   for k, v in meth.items()}

    if "--run" not in sys.argv:
        print("\n  §10 step 1 complete and PASSED. No statistic computed.")
        print("  M-2b itself is gated on G-a (§0). Re-run with --run once G-a is in.")
        (_RESULTS / "diag_m2b_contrast_agreement_results.json").write_text(
            json.dumps(out, indent=2, default=str), encoding="utf-8")
        return 0

    # ---- M-2b -------------------------------------------------------------------- #
    print("\nM-2b — do the two modalities agree on the SSEA4 - CD13 contrast?\n")
    print(f"  bar: sign agreement >= {SIGN_BAR}/11 (moved from 8/11 by the §6 freeze -- WEAKER "
          "than registered)")
    print("  rho(delta_rna, delta_meth) is reported but is NOT a bar.\n")
    age_rna = rna_ages(gill_dir)
    rna_pairs = pair_contrasts(age_rna, j["matched_keys"])
    from scipy.stats import spearmanr
    out["clocks"] = {}
    for cname, blk in meth.items():
        meth_pairs = pair_contrasts(blk["lp"], j["matched_keys"])
        idx = {(p["donor"], p["day"]): p for p in meth_pairs}
        rows, dr, dm = [], [], []
        for p in rna_pairs:
            q = idx.get((p["donor"], p["day"]))
            if q is None:
                continue
            d_m = q["delta"] * ADULT_SLOPE          # lp difference -> years, intercept-free
            rows.append({"donor": p["donor"], "day": p["day"],
                         "delta_rna_years": p["delta"], "delta_meth_years": d_m,
                         "same_sign": (p["delta"] > 0) == (d_m > 0)})
            dr.append(p["delta"])
            dm.append(d_m)
        sign = sign_agreement(dr, dm)
        rho = float(spearmanr(dr, dm).correlation) if len(dr) >= 4 else float("nan")
        v = m2b_verdict(sign, rho)
        v["pairs"] = rows
        v["mean_delta_rna_years"] = float(np.mean(dr)) if dr else float("nan")
        v["mean_delta_meth_years"] = float(np.mean(dm)) if dm else float("nan")
        v["by_day"] = agreement_by_day(rows)
        out["clocks"][cname] = v
        print(f"  === {cname} ===   {len(rows)} matched (donor, day) pairs")
        print(f"     mean delta  RNA {v['mean_delta_rna_years']:+7.2f} yr   "
              f"METH {v['mean_delta_meth_years']:+7.2f} yr")
        for r in rows:
            print(f"        {r['donor']} d{r['day']:<3} RNA {r['delta_rna_years']:+8.2f}   "
                  f"METH {r['delta_meth_years']:+7.2f}   "
                  f"{'agree' if r['same_sign'] else 'DISAGREE'}")
        print(f"     sign agreement {sign['n_agree']}/{sign['n_pairs']} (bar >= {SIGN_BAR})   "
              f"rho {rho:+.3f}   -> {v['status']}")
        print("     BY DAY (descriptive — §4's confound applies to M-2b too):")
        for d, b in v["by_day"].items():
            print(f"        day {d:>2}:  {b['n_agree']}/{b['n']} agree   "
                  f"mean RNA {b['mean_rna']:+7.2f}   mean METH {b['mean_meth']:+7.2f}")
        if v.get("caveat"):
            print(f"     [!] {v['caveat']}")
        print()

    st = {v["status"] for v in out["clocks"].values()}
    out["verdict_m2b"] = st.pop() if len(st) == 1 else "SPLIT"
    print(f"  ==> M-2b VERDICT: {out['verdict_m2b']}"
          f"{'  -- §6: a criterion met on one clock and not the other is SPLIT' if out['verdict_m2b'] == 'SPLIT' else ''}")

    (_RESULTS / "diag_m2b_contrast_agreement_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("\n  wrote diag_m2b_contrast_agreement_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
