"""STAGE 1.5.1 REV FINAL §6.3 / §10.6 — are O1 and O2 the SAME PHYSICAL DONORS across the series?

    python experiments/diag_donor_identity.py                       # pre-register only
    python experiments/diag_donor_identity.py --run [178dir] [179dir]

READ-ONLY. Writes `diag_donor_identity_results.json`. `src/` untouched.

THE CLAIM THIS TESTS
--------------------
`STAGE_1_5_1_REV_FINAL.md` §10.6 lists this under "claims deliberately NOT validated":

    "O1/O2 are the same physical donors in GSE165176 and GSE165179 — ❌ not verifiable from the
     metadata — matched by label and age only."

**"Not verifiable from the METADATA" is true. It is verifiable from the DATA.** Both GSE165178 and
GSE165179 are Illumina methylation arrays, and methylation carries a genotype fingerprint: a subset
of CpG probes overlaps common SNPs and is driven by the donor's genome rather than by cell state.
That is the standard basis for array-based sample identity checking.

WHY THIS IS THE RIGHT PAIR OF SERIES
------------------------------------
§10.6 names GSE165176 (RNA) vs GSE165179 (methylation), which cannot be compared directly — one has
no genotype signal to speak of. But the identity question is really about the DONOR LABELS, and
those propagate:

    GSE165176 (Sendai RNA)  <-> GSE165178 (Sendai array)     same experiment, 22/22 join (verified)
    GSE165177 (transient RNA) <-> GSE165179 (transient array) same experiment, 90 pairs (verified)

So "is Sendai-O1 the same person as transient-O1?" is answerable **methylation-to-methylation**, and
the answer transfers to §10.6's question by the two within-experiment joins already established.

THE DESIGN, AND ITS BUILT-IN CONTROLS
-------------------------------------
The rosters overlap only partially, and that asymmetry is the control:

    GSE165178 (query)   O1  O2  Y1  Y2
    GSE165179 (target)  O1  O2  O3

* **O1 and O2 are the true positives** — each should match its own label best.
* **Y1 and Y2 have NO counterpart in GSE165179**, and **O3 has none in GSE165178**. They cannot match
  correctly, so they measure what a spurious match looks like. Without them a high correlation
  everywhere (batch, array chemistry, cell type) would be indistinguishable from identity.

PROBE SELECTION CANNOT SEE THE ANSWER
-------------------------------------
Probes are chosen using **GSE165179 alone**, and within it using **only the untreated
`Negative control fibroblast` arm** — never GSE165178, never a treated sample. Selection is by
between-donor F-statistic: variance of donor means over pooled within-donor variance. Genotype-driven
probes score high on that by construction; probes that track cell state do not, because within this
arm there is no treatment to track.

The selected set is then **checked for cell-state stability inside GSE165179** (untreated vs
transiently reprogrammed, same donors) before it is used across series. A probe set that moved with
reprogramming would answer a different question.
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

N_PROBES = 5000          # genotype-like panel size; the published rs panels are 59-65, so this is
                         # deliberately generous — a small hand-picked set could be cherry-picked
N_SIM = 20000
RNG = np.random.default_rng(0)

NC_ARM = "Negative control fibroblast"
TR_ARM = "Transiently reprogrammed fibroblast"
SHARED_LABELS = ("O1", "O2")      # the labels present in BOTH series -- the only testable ones
STABILITY_BAR = 0.95              # cross-arm correlation the panel must hold inside GSE165179


# --------------------------------------------------------------------------- #
# Pure logic — data-free, unit-tested                                          #
# --------------------------------------------------------------------------- #
def between_donor_f(beta: np.ndarray, donors: list[str]) -> np.ndarray:
    """Per-probe F = var(donor means) / mean(within-donor var). Pure.

    High F means the probe separates donors and is stable inside a donor — the behaviour of a
    genotype-driven probe. Probes constant across everything score 0 and are filtered out by the
    epsilon in the denominator rather than by a special case.
    """
    b = np.asarray(beta, float)
    groups = sorted(set(donors))
    d = np.asarray(donors)
    means = np.stack([b[:, d == g].mean(axis=1) for g in groups], axis=1)
    within = np.stack([b[:, d == g].var(axis=1, ddof=1) if (d == g).sum() > 1
                       else np.zeros(b.shape[0]) for g in groups], axis=1).mean(axis=1)
    return means.var(axis=1, ddof=1) / (within + 1e-6)


def assign(query: dict[str, np.ndarray], target: dict[str, np.ndarray]) -> dict:
    """Match each query donor to its best-correlating target donor. Pure.

    Returns per-query the ranked correlations and the MARGIN (best minus runner-up). The margin is
    what separates "this is the same person" from "everything correlates because it is the same
    array chemistry" -- a query with no true counterpart should have a small margin however high
    its top correlation is.
    """
    tnames = sorted(target)
    out = {}
    for q, qv in sorted(query.items()):
        cors = {t: float(np.corrcoef(qv, target[t])[0, 1]) for t in tnames}
        ranked = sorted(cors.items(), key=lambda kv: -kv[1])
        margin = ranked[0][1] - ranked[1][1] if len(ranked) > 1 else float("nan")
        out[q] = {"correlations": cors, "best": ranked[0][0], "best_r": ranked[0][1],
                  "runner_up": ranked[1][0] if len(ranked) > 1 else None,
                  "margin": margin, "has_true_counterpart": q in tnames,
                  "correct": (q in tnames) and ranked[0][0] == q}
    return out


def identity_verdict(assignment: dict, shared: tuple[str, ...] = SHARED_LABELS) -> dict:
    """Does the genotype evidence support the same-donor assumption? Pure.

    Two conditions, both pre-registered:
      1. every shared label matches ITSELF (the positive claim);
      2. the shared labels' margins EXCEED every no-counterpart donor's margin (the control) --
         otherwise condition 1 could be satisfied by chance at a rate of 1 in 3 per donor.
    """
    testable = [q for q in shared if q in assignment]
    if not testable:
        return {"status": "CANNOT_VERIFY", "reason": "no shared donor labels between the series"}
    correct = [q for q in testable if assignment[q]["correct"]]
    ctrl = [v["margin"] for q, v in assignment.items() if not v["has_true_counterpart"]]
    shared_margins = [assignment[q]["margin"] for q in testable]
    sep = (min(shared_margins) > max(ctrl)) if ctrl else None
    if len(correct) == len(testable) and sep:
        status = "SAME_DONORS"
        reason = (f"all {len(testable)} shared labels match themselves, and their margins "
                  f"(min {min(shared_margins):.4f}) exceed every no-counterpart donor's "
                  f"(max {max(ctrl):.4f}). The label match is genotype-backed, not assumed.")
    elif len(correct) == len(testable):
        status = "SAME_DONORS_WEAK"
        reason = ("all shared labels match themselves, but their margins do not separate from the "
                  "no-counterpart controls, so the match is not distinguishable from chance at "
                  "1-in-3 per donor.")
    elif correct:
        status = "INCONSISTENT"
        reason = (f"only {len(correct)} of {len(testable)} shared labels match themselves — the "
                  "label mapping is not reliable and §6.3's assumption fails for at least one.")
    else:
        status = "DIFFERENT_DONORS"
        reason = ("no shared label matches itself. The donor labels do NOT denote the same people "
                  "across the two series, and every cross-series statement keyed on them is void.")
    return {"status": status, "reason": reason, "n_testable": len(testable),
            "n_correct": len(correct), "shared_margins": shared_margins,
            "control_margins": ctrl, "margins_separate": sep}


def sim_random_assignment(n_shared: int, n_targets: int, n_sim: int = N_SIM) -> np.ndarray:
    """How often does a panel with NO identity signal get every shared label right by chance?"""
    hits = RNG.integers(0, n_targets, size=(n_sim, n_shared)) == 0
    return hits.all(axis=1).astype(float)


# --------------------------------------------------------------------------- #
# Data wiring                                                                  #
# --------------------------------------------------------------------------- #
def _series_meta(path: Path) -> dict:
    """{title: donor} from a GEO series matrix. Donor is the leading token of the title."""
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                t = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                return {x: x.split("_")[0] for x in t}
    return {}


def _cell_types(path: Path) -> dict:
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        titles, ct = [], None
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                v = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                if v and v[0].split(":")[0].strip() == "cell type":
                    ct = [x.split(":", 1)[1].strip() for x in v]
        return dict(zip(titles, ct, strict=True)) if ct else {}


def _sample_cols(header: list[str]) -> list[tuple[int, str]]:
    """Column indices of real samples; `Detection Pval` columns are interleaved and skipped."""
    return [(i, h) for i, h in enumerate(header[1:], start=1) if h.strip() != "Detection Pval"]


def stream_betas(path: Path, want_cols: list[str] | None = None,
                 want_probes: set[str] | None = None):
    """Yield (probe, {sample: beta}) from a processed matrix, reading only the columns asked for."""
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        rdr = csv.reader(fh)
        keep = _sample_cols(next(rdr))
        if want_cols is not None:
            wanted = set(want_cols)
            keep = [(i, h) for i, h in keep if h in wanted]
        for row in rdr:
            if not row:
                continue
            p = row[0].strip()
            if want_probes is not None and p not in want_probes:
                continue
            vals = {}
            for i, s in keep:
                if i < len(row):
                    try:
                        vals[s] = float(row[i])
                    except ValueError:
                        pass
            if len(vals) == len(keep):
                yield p, vals


def select_panel(m179: Path, meta179: dict, ct179: dict) -> tuple[list[str], list[str], list[str]]:
    """Pick the genotype-like panel from GSE165179's UNTREATED arm only. Never sees GSE165178."""
    nc = sorted(t for t, c in ct179.items() if c == NC_ARM)
    donors = [meta179[t] for t in nc]
    print(f"  selecting the panel from {len(nc)} untreated samples "
          f"({len(set(donors))} donors), GSE165179 only")
    probes, rows = [], []
    for p, v in stream_betas(m179, want_cols=nc):
        probes.append(p)
        rows.append([v[s] for s in nc])
    beta = np.asarray(rows, float)
    print(f"  scanned {beta.shape[0]} probes")
    f = between_donor_f(beta, donors)
    # a probe that never moves cannot fingerprint anything, however high its F ratio looks
    rng_ok = (beta.max(axis=1) - beta.min(axis=1)) > 0.2
    f = np.where(rng_ok, f, -np.inf)
    idx = np.argsort(-f)[:N_PROBES]
    return [probes[i] for i in sorted(idx)], nc, donors


def profiles(path: Path, panel: set[str], meta: dict, only: list[str] | None = None
             ) -> tuple[dict[str, np.ndarray], list[str]]:
    """Mean beta per donor over the panel, in a fixed probe order."""
    cols = only if only is not None else list(meta)
    acc: dict[str, list] = defaultdict(list)
    order = []
    for p, v in stream_betas(path, want_cols=cols, want_probes=panel):
        order.append(p)
        for s, b in v.items():
            acc[meta[s]].append(b)
    n = len(order)
    out = {}
    for d in acc:
        cols_d = [s for s in cols if meta[s] == d]
        arr = np.asarray(acc[d], float).reshape(n, len(cols_d))
        out[d] = arr.mean(axis=1)
    return out, order


def main() -> int:
    print("STAGE 1.5.1 REV FINAL §6.3 — are the donor labels the same PEOPLE across series?\n")
    print("  PHASE 1: pre-registration. No beta value is read in this phase.\n")
    sim = sim_random_assignment(len(SHARED_LABELS), 3)
    p_chance = float(sim.mean())
    print(f"  a panel with NO identity signal gets both shared labels right "
          f"{p_chance:.1%} of the time (1 in {1/max(p_chance,1e-9):.0f}).")
    print("  That is why a correct assignment ALONE is not the bar. The pre-registered bar is:")
    print("     (1) every shared label matches ITSELF, AND")
    print("     (2) the shared labels' margins exceed every no-counterpart donor's margin.")
    print(f"  Panel stability inside GSE165179 (untreated vs reprogrammed) must be "
          f">= {STABILITY_BAR:.2f}.\n")

    out: dict = {"script": "diag_donor_identity",
                 "utc": datetime.now(UTC).isoformat(timespec="seconds"),
                 "preregistration": {"p_both_correct_by_chance": p_chance,
                                     "n_probes": N_PROBES, "stability_bar": STABILITY_BAR,
                                     "shared_labels": list(SHARED_LABELS)}}

    if "--run" not in sys.argv:
        print("  Pre-registration only. Re-run with --run to measure.")
        Path("diag_donor_identity_results.json").write_text(
            json.dumps(out, indent=2, default=str), encoding="utf-8")
        return 0

    pos = [a for a in sys.argv[1:] if not a.startswith("--")]
    d178 = Path(pos[0] if pos else r"D:\GSE165178")
    d179 = Path(pos[1] if len(pos) > 1 else r"D:\GSE165179")
    m178 = d178 / "GSE165178_Matrix_processed_sendai.txt.gz"
    m179 = d179 / "GSE165179_Matrix_processed_transient.txt.gz"

    print("  PHASE 2: measurement.\n")
    meta178, meta179 = _series_meta(d178 / "GSE165178_series_matrix.txt.gz"), \
        _series_meta(d179 / "GSE165179_series_matrix.txt.gz")
    ct179 = _cell_types(d179 / "GSE165179_series_matrix.txt.gz")
    meta179 = {t: d for t, d in meta179.items() if d != "iPSC"}

    panel, _nc, _don = select_panel(m179, meta179, ct179)
    pset = set(panel)
    print(f"  panel: {len(panel)} probes\n")

    # --- stability: does the panel move with reprogramming INSIDE GSE165179? --------- #
    nc = sorted(t for t, c in ct179.items() if c == NC_ARM and t in meta179)
    tr = sorted(t for t, c in ct179.items() if c == TR_ARM and t in meta179)
    p_nc, o1 = profiles(m179, pset, meta179, only=nc)
    p_tr, o2 = profiles(m179, pset, meta179, only=tr)
    assert o1 == o2, "probe order differs between the two passes"
    stab = {d: float(np.corrcoef(p_nc[d], p_tr[d])[0, 1]) for d in sorted(p_nc) if d in p_tr}
    print("  panel stability inside GSE165179 (untreated vs transiently reprogrammed, same donor):")
    for d, r in stab.items():
        print(f"     {d}: r = {r:.4f}")
    stable = min(stab.values()) >= STABILITY_BAR
    out["panel_stability"] = {"per_donor": stab, "min": min(stab.values()),
                              "bar": STABILITY_BAR, "verdict": "STABLE" if stable else "MOVES"}
    if not stable:
        print(f"\n  [!] panel min r {min(stab.values()):.4f} < {STABILITY_BAR} — the panel tracks "
              "cell state, so it cannot arbitrate identity. Reporting and stopping.")
        Path("diag_donor_identity_results.json").write_text(
            json.dumps(out, indent=2, default=str), encoding="utf-8")
        return 0

    # --- the cross-series assignment ------------------------------------------------- #
    tgt, ordt = profiles(m179, pset, meta179, only=sorted(meta179))
    qry, ordq = profiles(m178, pset, meta178, only=sorted(meta178))
    common = [p for p in ordt if p in set(ordq)]
    it = {p: i for i, p in enumerate(ordt)}
    iq = {p: i for i, p in enumerate(ordq)}
    tgt = {d: v[[it[p] for p in common]] for d, v in tgt.items()}
    qry = {d: v[[iq[p] for p in common]] for d, v in qry.items()}
    print(f"\n  {len(common)} panel probes present in BOTH series")
    print(f"  query  (GSE165178, Sendai):    {sorted(qry)}")
    print(f"  target (GSE165179, transient): {sorted(tgt)}\n")

    a = assign(qry, tgt)
    out["assignment"] = a
    out["n_common_probes"] = len(common)
    hdr = sorted(tgt)
    print(f"  {'query':<8}" + "".join(f"{t:>10}" for t in hdr) + f"{'best':>8}{'margin':>10}")
    print("  " + "-" * (8 + 10 * len(hdr) + 18))
    for q in sorted(a):
        row = a[q]
        cells = "".join(f"{row['correlations'][t]:>10.4f}" for t in hdr)
        flag = "  <- no counterpart" if not row["has_true_counterpart"] else (
            "  OK" if row["correct"] else "  <- MISMATCH")
        print(f"  {q:<8}{cells}{row['best']:>8}{row['margin']:>10.4f}{flag}")

    v = identity_verdict(a)
    out["verdict"] = v
    print(f"\n  ==> §6.3 VERDICT: {v['status']}\n      {v['reason']}")

    Path("diag_donor_identity_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("\n  wrote diag_donor_identity_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
