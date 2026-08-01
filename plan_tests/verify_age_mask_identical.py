"""STAGE 1.5.3 steps 1-3 — the bit-identity gate: did any label move?

    python plan_tests/verify_age_mask_identical.py --capture   # BEFORE any change
    python plan_tests/verify_age_mask_identical.py --verify    # after each step

READ-ONLY with respect to `src/`. Writes `results/verify_age_mask_identical_results.json`.

WHAT THIS GATES
---------------
`STAGE_1_5_3_EXECUTE.md` §5 says steps 1-3 "cannot move a number", and §6 requires that be shown by
`np.array_equal` on real data rather than argued from the diff. C-1's whole design rests on it: the
new policy is a no-op at its defaults, so ΔAge and `age_mask` must come out **bit-identical**, and
the new `age_mask_reason` must be all-`None`.

This is the same discipline gate G-a shipped under, and the same one `STAGE_1_5_2_LABEL_ANCHOR.md`
§8.2 demands of any label change: *"If any ΔAge moves, the change is wrong -- revert, do not
rationalise."*

THE BAR, AND WHY ITS RESOLVABILITY IS ARGUED RATHER THAN SIMULATED
------------------------------------------------------------------
    max |ΔAge_after - ΔAge_before| == 0.0 exactly, AND age_mask identical elementwise.

`REF_GROUND_RULES.md` §5b asks whether a system meeting the intent EXACTLY would pass. Here the
intent is "the computation is unchanged", the computation is deterministic, and the comparison is
exact equality -- so a correct change passes with probability **1 by construction**. Simulating a
null would measure nothing.

**The meaningful check runs in the other direction: can this bar DETECT a violation?** A gate that
cannot fail is not a gate (the `verify_1a` lesson, which graded PASS on a warning it had itself
printed). So `--verify` also runs a self-test that injects a perturbation of 1 ULP into a copy of
the captured baseline and confirms the comparison rejects it. That self-test runs on every
invocation, and the script **aborts** if the bar cannot catch its own violation.

GEOMETRY
--------
All six Gill donors (bulk, fast, and the only cells with `donor_age` today) plus one HFF chunk from
GSE242423 (single-cell, ~980 cells, the 99.7% population). `--capture` stores a SHA-256 of the
concatenated arrays plus the raw values, so a later comparison cannot be fooled by a re-ordering.
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)
_BASELINE = _RESULTS / "verify_age_mask_identical_baseline.json"
_REPORT = _RESULTS / "verify_age_mask_identical_results.json"

GILL_DIR = r"D:\Gill"
GSE_DIR = r"D:\GSE242423"


# --------------------------------------------------------------------------- #
# Pure logic — data-free, unit-tested                                          #
# --------------------------------------------------------------------------- #
def fingerprint(d_age: np.ndarray, age_mask: np.ndarray) -> str:
    """Order-sensitive SHA-256 of (ΔAge, age_mask). Pure.

    Hashing the raw bytes of float64 means a change of one ULP changes the digest, which is the
    point: `allclose` would hide exactly the drift this gate exists to catch.
    """
    h = hashlib.sha256()
    h.update(np.ascontiguousarray(d_age, dtype=np.float64).tobytes())
    h.update(np.ascontiguousarray(age_mask, dtype=np.bool_).tobytes())
    return h.hexdigest()


def compare(before: dict, after: dict) -> dict:
    """Did anything move? Pure. Exact equality, never a tolerance."""
    keys = sorted(set(before) | set(after))
    rows, worst = [], 0.0
    ok = True
    for k in keys:
        b, a = before.get(k), after.get(k)
        if b is None or a is None:
            rows.append({"chunk": k, "status": "MISSING",
                         "detail": f"{'after' if a is None else 'before'} has no entry"})
            ok = False
            continue
        db, da = np.asarray(b["d_age"], float), np.asarray(a["d_age"], float)
        mb, ma = np.asarray(b["age_mask"], bool), np.asarray(a["age_mask"], bool)
        if db.shape != da.shape or mb.shape != ma.shape:
            rows.append({"chunk": k, "status": "SHAPE_CHANGED",
                         "detail": f"{db.shape}->{da.shape}"})
            ok = False
            continue
        d_ident = bool(np.array_equal(db, da))
        m_ident = bool(np.array_equal(mb, ma))
        hash_ident = b["sha256"] == a["sha256"]
        delta = float(np.max(np.abs(db - da))) if db.size else 0.0
        worst = max(worst, delta)
        # `age_mask_reason` must be entirely absent-or-None while the policies are off
        reasons = a.get("reasons_not_none", 0)
        row_ok = d_ident and m_ident and hash_ident and reasons == 0
        ok &= row_ok
        rows.append({"chunk": k, "n": int(db.size), "max_abs_delta": delta,
                     "d_age_identical": d_ident, "age_mask_identical": m_ident,
                     "sha256_identical": hash_ident, "reasons_not_none": reasons,
                     "status": "OK" if row_ok else "MOVED"})
    return {"verdict": "IDENTICAL" if ok else "MOVED", "max_abs_delta": worst,
            "n_chunks": len(keys), "rows": rows}


def self_test() -> dict:
    """Can the bar detect a violation? If not, it is not a gate. Pure.

    Three injected faults, each of the kind this gate exists to catch.
    """
    base = {"c": {"d_age": [1.0, 2.0, 3.0], "age_mask": [True, True, False],
                  "sha256": fingerprint(np.array([1.0, 2.0, 3.0]),
                                        np.array([True, True, False])),
                  "reasons_not_none": 0}}
    checks = {}

    # 1 ULP on one value -- the smallest possible float change
    d = [1.0, np.nextafter(2.0, 3.0), 3.0]
    pert = {"c": {**base["c"], "d_age": d,
                  "sha256": fingerprint(np.array(d), np.array([True, True, False]))}}
    checks["one_ulp"] = compare(base, pert)["verdict"] == "MOVED"

    # a flipped mask bit with ΔAge untouched
    m = [True, False, False]
    pert = {"c": {**base["c"], "age_mask": m,
                  "sha256": fingerprint(np.array([1.0, 2.0, 3.0]), np.array(m))}}
    checks["mask_flip"] = compare(base, pert)["verdict"] == "MOVED"

    # a reason string appearing while the policies are supposed to be off
    checks["reason_appeared"] = compare(
        base, {"c": {**base["c"], "reasons_not_none": 1}})["verdict"] == "MOVED"

    # and the control: an unchanged copy must PASS, or the gate is useless the other way
    checks["unchanged_passes"] = compare(base, base)["verdict"] == "IDENTICAL"

    return {"checks": checks, "all_pass": all(checks.values())}


# --------------------------------------------------------------------------- #
# Data wiring                                                                  #
# --------------------------------------------------------------------------- #
def _delta_age_call(clock, norm, genes, obs, source):
    """Call `delta_age` through whichever signature is current.

    Before C-1 it returns 2 values; after C-1 it returns 3. This gate has to run on BOTH sides
    of that change -- that is its entire job -- so it adapts instead of assuming.
    """
    from cellfate.data.aging import delta_age
    out = delta_age(clock, norm, genes, obs, source)
    if len(out) == 3:
        d, mask, reasons = out
        return d, mask, sum(1 for r in reasons if r is not None)
    d, mask = out
    return d, mask, 0


def collect() -> dict:
    """ΔAge + age_mask for all six Gill donors and one HFF chunk."""
    from cellfate.data.aging import LinearClock
    from cellfate.data.normalize import normalize_counts
    from cellfate.data.sources import GillReprogrammingSource

    clock = LinearClock.from_json(ROOT / "configs" / "clocks" / "fleischer_clock.json")
    out: dict = {}

    gill = GillReprogrammingSource(
        str(next(Path(GILL_DIR).glob("*Log2_RPM_Sendai*.txt.gz"))),
        str(Path(GILL_DIR) / "GSE165176_series_matrix.txt.gz"))
    for ch in gill.plan():
        raw = gill.fetch(ch)
        norm = normalize_counts(raw.counts)
        d, m, n_reasons = _delta_age_call(clock, norm, raw.genes, raw.obs, raw.source)
        out[ch["id"]] = {"d_age": [float(x) for x in d], "age_mask": [bool(x) for x in m],
                         "sha256": fingerprint(d, m), "reasons_not_none": n_reasons}
        print(f"   {ch['id']:<28} n={len(d):>4}  sha {out[ch['id']]['sha256'][:12]}")

    # one HFF chunk -- the 99.7% population, and the only source that exercises C-3
    try:
        sys.path.insert(0, str(ROOT / "local_runners"))
        from run_multi_local import discover_gse  # type: ignore

        from cellfate.data.sources import GSE242423SingleCellSource
        samples, genes_file = discover_gse(GSE_DIR)
        src = GSE242423SingleCellSource(samples, genes_file, cell_line="HFF", min_genes=500,
                                        max_cells_per_sample=200, cells_per_run=None, seed=0)
        ch = src.plan()[0]
        raw = src.fetch(ch)
        norm = normalize_counts(raw.counts)
        d, m, n_reasons = _delta_age_call(clock, norm, raw.genes, raw.obs, raw.source)
        out[ch["id"]] = {"d_age": [float(x) for x in d], "age_mask": [bool(x) for x in m],
                         "sha256": fingerprint(d, m), "reasons_not_none": n_reasons}
        print(f"   {ch['id']:<28} n={len(d):>4}  sha {out[ch['id']]['sha256'][:12]}")
    except Exception as exc:  # noqa: BLE001 -- recorded, never silently skipped
        print(f"   [!] HFF chunk NOT captured: {exc!r}")
        print("       The Gill donors still gate ΔAge; C-3's stamping is covered by unit tests.")
    return out


def main() -> int:
    mode = "--capture" if "--capture" in sys.argv else "--verify"
    print(f"STAGE 1.5.3 bit-identity gate  [{mode}]\n")

    st = self_test()
    print("  SELF-TEST -- can this bar detect a violation?")
    for k, v in st["checks"].items():
        print(f"     {k:<20} {'yes' if v else 'NO'}")
    if not st["all_pass"]:
        print("\n  ==> ABORT: the gate cannot catch its own violation, so it is not a gate.")
        return 1
    print("     -> the bar can fail, so a PASS means something.\n")

    if mode == "--capture":
        data = collect()
        _BASELINE.write_text(json.dumps(
            {"utc": datetime.now(UTC).isoformat(timespec="seconds"), "chunks": data},
            indent=2), encoding="utf-8")
        print(f"\n  captured {len(data)} chunks -> {_BASELINE.name}")
        print("  Run --verify after EACH of steps 1, 2 and 3.")
        return 0

    if not _BASELINE.exists():
        print(f"  [!] no baseline at {_BASELINE}. Run --capture BEFORE changing src/.")
        return 1
    before = json.loads(_BASELINE.read_text(encoding="utf-8"))["chunks"]
    after = collect()
    res = compare(before, after)
    print(f"\n  {'chunk':<30}{'n':>6}{'max|Δ|':>12}  status")
    print("  " + "-" * 62)
    for r in res["rows"]:
        print(f"  {r['chunk']:<30}{r.get('n', 0):>6}{r.get('max_abs_delta', float('nan')):>12.2e}"
              f"  {r['status']}")
    print(f"\n  ==> {res['verdict']}   max|Δ| = {res['max_abs_delta']:.2e}")
    if res["verdict"] != "IDENTICAL":
        print("      STAGE_1_5_2 §8.2: if any ΔAge moves, the change is WRONG -- revert, "
              "do not rationalise.")
    (_REPORT).write_text(json.dumps(
        {"script": "verify_age_mask_identical",
         "utc": datetime.now(UTC).isoformat(timespec="seconds"),
         "self_test": st, **res}, indent=2), encoding="utf-8")
    print(f"\n  wrote {_REPORT.name}")
    return 0 if res["verdict"] == "IDENTICAL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
