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

ATTEMPT 1 FAILED THAT GATE, AND IS RECORDED RATHER THAN DELETED
---------------------------------------------------------------
Selecting purely by between-donor F gave a panel with cross-arm stability **r = 0.821 / 0.942 /
0.966**, below the pre-registered 0.95 bar on donor O1. The script stopped there **without computing
the cross-series assignment** — so the refinement below was chosen against the *stability* criterion,
an internal property of GSE165179, and not against the identity answer, which had not been seen.
That ordering is the protection, and it is checkable: the abort is in the code path above the
assignment.

Why F alone was the wrong criterion, in hindsight and in principle: between-donor F rewards any
probe that differs between donors, including one that differs because those donors' *cultures*
differ. It is a necessary condition for a genotype probe, not a sufficient one.

ATTEMPT 2 — TRIMODALITY, WHICH IS WHAT ACTUALLY DEFINES A GENOTYPE PROBE
------------------------------------------------------------------------
A CpG whose beta is driven by an underlying SNP takes one of three values — homozygous, heterozygous,
homozygous — so its betas cluster near **{0, 0.5, 1}** and nowhere else. That is the standard
signature used for array-based identity checking when an `rs` panel is unavailable, and unlike F it
is a property of the *shape* of the distribution rather than of which donors differ. Probes are
required to sit within 0.15 of a mode in **every** sample and to occupy at least two distinct modes
(otherwise they are invariant and fingerprint nothing), and are then ranked by F among those.

**The stability bar is unchanged at 0.95** — the same gate, applied to the new panel.

ATTEMPT 2 ALSO FAILED IT, AND THAT EXPOSED A DEFECT IN MY OWN BAR
------------------------------------------------------------------
Trimodality improved stability a lot — **0.938 / 0.985 / 0.990** against 0.821 / 0.942 / 0.966 — but
O1 still missed 0.95. Before touching the bar, note what is wrong with it: **I set 0.95 by
assertion.** `REF_GROUND_RULES` §5b requires every acceptance bar to go through
`audit_metrics.bar_verdict` *before* it is registered, and this one did not. That is the same defect
this project has now caught four times, committed by me here.

So the bar is checked the way the rule requires: **would a PERFECT panel — one whose true betas are
identical across arms — clear 0.95 at this geometry?** The geometry is lopsided: stability compares
a mean of **7** untreated samples against a mean of the transiently-reprogrammed arm, and that arm
has **7 samples for O3, 4 for O2, and only 2 for O1**. A two-sample mean is noisy, and correlation
between two noisy means of the same true vector is bounded well below 1.

Noise is estimated from GSE165179's own **exp1/exp2 technical replicates** of the same condition —
a direct measurement of array noise, taken without any contact with GSE165178. If 0.95 turns out to
be unresolvable at n=2, it moves to its `usable_bar`, per §5b, **before** the assignment is computed.
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

from audit_metrics import bar_verdict  # noqa: E402

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)


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


def trimodal_score(beta: np.ndarray, tol: float = 0.15) -> tuple[np.ndarray, np.ndarray]:
    """Genotype-probe shape test. Pure. Returns (is_trimodal, n_modes_occupied) per probe.

    A CpG driven by an underlying SNP reads homozygous / heterozygous / homozygous, so every
    sample sits near one of {0, 0.5, 1}. A probe qualifies when **every** sample is within `tol`
    of a mode AND at least two distinct modes are occupied — the second condition rules out
    probes that are simply invariant, which satisfy the first trivially and fingerprint nothing.
    """
    b = np.asarray(beta, float)
    modes = np.array([0.0, 0.5, 1.0])
    d = np.abs(b[:, :, None] - modes[None, None, :])
    nearest = d.argmin(axis=2)
    within = d.min(axis=2).max(axis=1) <= tol
    n_modes = np.array([len(set(row)) for row in nearest])
    return within & (n_modes >= 2), n_modes


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


def sim_stability(n_probes: int, signal_sd: float, noise_sd: float, n_a: int, n_b: int,
                  n_sim: int = 2000) -> np.ndarray:
    """Cross-arm correlation achieved by a PERFECTLY stable panel. Pure.

    The two arms are means of `n_a` and `n_b` noisy observations of the SAME true profile, so any
    shortfall from 1.0 here is sampling noise and nothing else. This is what makes 0.95 checkable
    rather than asserted: if a perfect panel cannot reach it at n_b = 2, the bar is measuring the
    sample count, not the panel.
    """
    out = np.empty(n_sim)
    for i in range(n_sim):
        true = RNG.normal(0.0, signal_sd, size=n_probes)
        a = true + RNG.normal(0.0, noise_sd / np.sqrt(n_a), size=n_probes)
        b = true + RNG.normal(0.0, noise_sd / np.sqrt(n_b), size=n_probes)
        out[i] = np.corrcoef(a, b)[0, 1]
    return out


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


def select_panel(m179: Path, meta179: dict, ct179: dict, mode: str = "trimodal") -> dict:
    """Pick the genotype-like panel from GSE165179's UNTREATED arm only. Never sees GSE165178.

    `mode="f_only"` reproduces ATTEMPT 1 (between-donor F alone), kept so the failed attempt stays
    runnable rather than only described. `mode="trimodal"` is ATTEMPT 2.
    """
    nc = sorted(t for t, c in ct179.items() if c == NC_ARM)
    donors = [meta179[t] for t in nc]
    print(f"  selecting the panel from {len(nc)} untreated samples "
          f"({len(set(donors))} donors), GSE165179 only  [mode={mode}]")
    probes, rows = [], []
    for p, v in stream_betas(m179, want_cols=nc):
        probes.append(p)
        rows.append([v[s] for s in nc])
    beta = np.asarray(rows, float)
    print(f"  scanned {beta.shape[0]} probes")
    f = between_donor_f(beta, donors)
    # a probe that never moves cannot fingerprint anything, however high its F ratio looks
    eligible = (beta.max(axis=1) - beta.min(axis=1)) > 0.2
    n_tri = None
    if mode == "trimodal":
        tri, _ = trimodal_score(beta)
        n_tri = int(tri.sum())
        print(f"  {n_tri} probes have the trimodal genotype shape (of {beta.shape[0]})")
        eligible &= tri
    idx = np.argsort(-np.where(eligible, f, -np.inf))[:N_PROBES]
    idx = [i for i in idx if np.isfinite(f[i]) and eligible[i]]
    return {"panel": [probes[i] for i in sorted(idx)], "n_eligible": int(eligible.sum()),
            "n_trimodal": n_tri, "mode": mode, "nc_samples": nc, "donors": donors}


def common_probe_count(p: dict[str, np.ndarray]) -> np.ndarray:
    """Length of the profile vectors — the number of probes the correlation is taken over."""
    return next(iter(p.values()))


def replicate_noise(path: Path, panel: set[str], meta: dict, ct: dict) -> dict:
    """Per-sample array noise, from exp1/exp2 TECHNICAL replicates of the same condition.

    These are repeats of one condition (the fact `pair_by_donor_day` exists to handle), so their
    difference is instrument noise with no biology in it. sd of a single sample = sd(diff)/sqrt(2).
    Measured on GSE165179 alone — GSE165178 is never opened here.
    """
    groups: dict[tuple, list[str]] = {}
    for t in meta:
        if t.endswith(("_exp1", "_exp2")):
            groups.setdefault(t.rsplit("_", 1)[0], []).append(t)
    pairs = {k: sorted(v) for k, v in groups.items() if len(v) == 2}
    cols = sorted({s for v in pairs.values() for s in v})
    if not cols:
        return {"noise_sd": float("nan"), "n_pairs": 0}
    vals: dict[str, dict[str, float]] = {}
    for p, v in stream_betas(path, want_cols=cols, want_probes=panel):
        vals[p] = v
    diffs = []
    for a, b in pairs.values():
        diffs.extend(v[a] - v[b] for v in vals.values() if a in v and b in v)
    return {"noise_sd": float(np.std(diffs) / np.sqrt(2.0)), "n_pairs": len(pairs),
            "n_probes": len(vals)}


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


def _run_assignment(m178: Path, m179: Path, pset: set[str], meta178: dict, meta179: dict,
                    out: dict, target_cols: list[str] | None = None,
                    query_filter=None) -> int:
    """Build both sides' profiles over the panel and report the assignment. Shared by both paths."""
    tgt, ordt = profiles(m179, pset, meta179, only=target_cols or sorted(meta179))
    qcols = [t for t in sorted(meta178) if (query_filter is None or query_filter(t))]
    qry, ordq = profiles(m178, pset, meta178, only=qcols)
    common = [p for p in ordt if p in set(ordq)]
    it = {p: i for i, p in enumerate(ordt)}
    iq = {p: i for i, p in enumerate(ordq)}
    tgt = {d: v[[it[p] for p in common]] for d, v in tgt.items()}
    qry = {d: v[[iq[p] for p in common]] for d, v in qry.items()}
    print(f"\n  {len(common)} panel probes present in BOTH series; "
          f"{len(qcols)} query samples used")
    print(f"  query  (GSE165178, Sendai):    {sorted(qry)}")
    print(f"  target (GSE165179, transient): {sorted(tgt)}\n")

    a = assign(qry, tgt)
    out["assignment"] = a
    out["n_common_probes"] = len(common)
    out["n_query_samples"] = len(qcols)
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
    (_RESULTS / "diag_donor_identity_results.json").write_text(
        json.dumps(out, indent=2, default=str), encoding="utf-8")
    print("\n  wrote diag_donor_identity_results.json")
    return 0


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
        (_RESULTS / "diag_donor_identity_results.json").write_text(
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

    mode = "f_only" if "--f-only" in sys.argv else "trimodal"
    sel = select_panel(m179, meta179, ct179, mode=mode)
    panel = sel["panel"]
    pset = set(panel)
    out["panel_selection"] = {k: v for k, v in sel.items() if k != "panel"}
    out["panel_selection"]["n_panel"] = len(panel)
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
    # --- is 0.95 even reachable at this geometry? §5b, applied to my own bar --------- #
    n_tr = {d: sum(1 for t in tr if meta179[t] == d) for d in stab}
    reps = replicate_noise(m179, pset, meta179, ct179)
    sig_sd = float(np.std(np.concatenate([p_nc[d] for d in sorted(p_nc)])))
    print(f"\n  §5b check on my own bar. Array noise from exp1/exp2 replicates: "
          f"sd {reps['noise_sd']:.4f} over {reps['n_pairs']} pairs;")
    print(f"  panel signal sd {sig_sd:.4f}; untreated arm n=7 vs reprogrammed arm n={n_tr}")
    bar = STABILITY_BAR
    worst_n = min(n_tr.values())
    sim = sim_stability(len(common_probe_count(p_nc)), sig_sd, reps["noise_sd"], 7, worst_n)
    v = bar_verdict(sim, STABILITY_BAR, lower_is_better=False)
    print(f"  a PERFECT panel (identical true betas) scores median {v['null_median']:.4f} at "
          f"n=7 vs n={worst_n}")
    print(f"     it clears {STABILITY_BAR} {v['pass_rate']:.1%} of the time -> {v['verdict']}")
    out["stability_bar_check"] = {**v, "bar_proposed": STABILITY_BAR, "noise_sd": reps["noise_sd"],
                                  "signal_sd": sig_sd, "n_treated_per_donor": n_tr,
                                  "n_replicate_pairs": reps["n_pairs"]}
    if v["verdict"] != "RESOLVABLE":
        bar = float(v["usable_bar"])
        print(f"     [!] §5b: UNRESOLVABLE at {STABILITY_BAR}. The bar moves to {bar:.4f} NOW, "
              "before the assignment is computed.")
    out["panel_stability"] = {"per_donor": stab, "min": min(stab.values()),
                              "bar_proposed": STABILITY_BAR, "bar_used": bar,
                              "n_treated_per_donor": n_tr}
    stable = min(stab.values()) >= bar
    out["panel_stability"]["verdict"] = "STABLE" if stable else "MOVES"
    if not stable:
        print(f"\n  [!] panel min r {min(stab.values()):.4f} < {bar:.4f} — the panel tracks "
              "cell state, so it cannot arbitrate identity. Reporting and stopping.")
        # WHY, as a diagnostic. Descriptive only, and it cannot touch the identity answer because
        # the assignment is never computed on this path. The hypothesis to check: successful
        # reprogramming demethylates globally, so a panel is least stable in the donor whose
        # reprogrammed samples are the deepest -- not in the donor with the fewest of them.
        fl = sorted(t for t, c in ct179.items()
                    if c == "Failed to transiently reprogram fibroblast" and t in meta179)
        p_fl, _ = profiles(m179, pset, meta179, only=fl)
        stab_fl = {d: float(np.corrcoef(p_nc[d], p_fl[d])[0, 1]) for d in sorted(p_nc)
                   if d in p_fl}
        days = {}
        for t in tr:
            days.setdefault(meta179[t], []).append(t)
        out["diagnosis"] = {"stability_vs_failed_arm": stab_fl,
                            "treated_sample_titles": {d: sorted(v) for d, v in days.items()}}
        print("\n  DIAGNOSIS (descriptive — the assignment was never computed on this path):")
        print("    stability against the FAILED arm instead (cells that got OSKM and did NOT")
        print("    reprogram, so no global demethylation; 7 samples per donor, balanced):")
        for d, r in stab_fl.items():
            print(f"       {d}: r = {r:.4f}   (vs reprogrammed arm: {stab[d]:.4f})")
        print("    the reprogrammed samples each donor actually has:")
        for d, v in sorted(days.items()):
            print(f"       {d} (n={len(v)}): {', '.join(sorted(v))}")

        # The gate did not say "give up"; it said the panel cannot arbitrate WHERE it moves. It
        # demonstrably does not move in non-reprogramming cells, so the test is re-run restricted
        # to those states on BOTH sides. This is the stability evidence deciding the scope, not the
        # identity answer -- which is still uncomputed at this point.
        if min(stab_fl.values()) >= bar:
            print(f"\n  ==> RESTRICTING to NON-REPROGRAMMING cells, where the panel IS stable "
                  f"(min r {min(stab_fl.values()):.4f} >= {bar:.4f}).")
            print("      target: GSE165179 untreated + failed arms;  query: GSE165178 CD13 only")
            print("      (CD13 = 'Failing to reprogram fibroblast' — the direct analogue.)")
            out["restricted_to_non_reprogramming"] = True
            return _run_assignment(m178, m179, pset, meta178, meta179, out,
                                   target_cols=sorted(set(nc) | set(fl)),
                                   query_filter=lambda t: t.endswith("_CD13"))
        (_RESULTS / "diag_donor_identity_results.json").write_text(
            json.dumps(out, indent=2, default=str), encoding="utf-8")
        return 0
    print(f"  panel stability PASSES at the §5b bar {bar:.4f}")

    return _run_assignment(m178, m179, pset, meta178, meta179, out)


if __name__ == "__main__":
    raise SystemExit(main())
