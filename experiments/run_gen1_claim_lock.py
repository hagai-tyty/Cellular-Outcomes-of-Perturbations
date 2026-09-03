"""Generation-1 claim lock.

Executes `GEN1_CLAIM_LOCK_V1.md`: freeze the abstract-level allowed and forbidden claims.

Every previous stage produced a number. This one produces sentences, and sentences are where a
project of this kind actually fails -- not in the statistics, but in the abstract, where "six
observed experimental conditions in one melanoma line" becomes "treatments in cancer" and a
detection proxy becomes "response". The evidence lock fixed what the claims are made of; this fixes
what may be said about it.

An abstract assembled entirely from permitted fragments can still read as a broader claim than any
of them. So listing allowed sentences is not enough: this also fires a corpus of tempting,
forbidden sentences at the same scanner Stage 26 used on the shipped tool, and fails if any goes
undetected. One instrument, one standard, for the code and for the prose.

Nothing here fits anything or produces a number. It refuses if the evidence lock does not verify.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import run_gen1_evidence_lock as EL  # noqa: E402
import run_stage26_scope_lock as S26  # noqa: E402

RESULTS = ROOT / "results"
OUT = RESULTS / "claim_lock"
OUT.mkdir(parents=True, exist_ok=True)

PLANS = ROOT / "plans" / "(newer)practical plans"
PLAN = PLANS / "GEN1_CLAIM_LOCK_V1.md"
SHIP_PLAN = PLANS / "STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md"
SHIP_PLAN_DIGEST = "8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48"

EVIDENCE_HANDOFF = RESULTS / "gen1_handoff_to_claim_lock.json"
EVIDENCE_MANIFEST = RESULTS / "evidence_lock" / "GEN1_EVIDENCE_MANIFEST.json"
CLAIM_INPUT = RESULTS / "evidence_lock" / "GEN1_CLAIM_LOCK_INPUT.json"

CLAIMS_JSON = OUT / "GEN1_CLAIM_LOCK.json"
CLAIMS_MD = OUT / "GEN1_CLAIMS.md"
ADVERSARIAL_JSON = OUT / "claim_lock_adversarial.json"
DIGEST_JSON = OUT / "GEN1_CLAIM_DIGEST.json"
HANDOFF_JSON = RESULTS / "gen1_handoff_to_manuscript.json"


# Wall-clock timings say nothing about content, and one of them sat inside a file this digest
# covers. The consequence was not cosmetic: re-running `--stage all` with nothing changed produced
# a DIFFERENT claim digest every time, so the value quoted in the manuscript and the README was
# valid for exactly one execution, and anyone reproducing the stage would conclude something had
# moved. A lock that cannot survive its own stage being re-run is not a lock.
VOLATILE_KEYS = {"runtime_seconds", "runtime_minutes", "total_runtime_seconds"}


def _strip_volatile(obj):
    if isinstance(obj, dict):
        return {k: _strip_volatile(v) for k, v in obj.items() if k not in VOLATILE_KEYS}
    if isinstance(obj, list):
        return [_strip_volatile(v) for v in obj]
    return obj


def stable_sha256(p: Path) -> str:
    """Hash CONTENT, not timing. JSON is normalised with volatile fields removed and keys sorted;
    everything else keeps the canonical-LF rule the evidence lock uses."""
    if p.suffix == ".json":
        payload = _strip_volatile(json.loads(p.read_text(encoding="utf-8")))
        blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(blob).hexdigest()
    return EL.canonical_lf_sha256(p)

def write_json(p: Path, obj: dict) -> dict:
    stamped = _strip_volatile({**obj, "module_sha256":
                               hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    p.write_text(json.dumps(stamped, indent=2) + "\n", encoding="utf-8")
    return stamped


def _j(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


# The evidence lock hashes every Gen-1 artifact -- 54 when this plan was written, 62 once the
# figures, source-data export and release bundle were added -- and none of them is this
# stage's output. That is correct
# layering, since re-locking to include them would make CL-A circular. But it leaves the document
# the manuscript is written from with no identity at all. So the claim lock hashes itself, by the
# same canonical-LF rule, and the manuscript binds to two numbers: the evidence digest for what the
# claims are made of, the claim digest for what may be said about it.
CLAIM_LOCK_FILES = [
    "plans/(newer)practical plans/GEN1_CLAIM_LOCK_V1.md",
    "experiments/run_gen1_claim_lock.py",
    "tests/test_gen1_claim_lock.py",
    "results/claim_lock/GEN1_CLAIMS.md",
]




def claim_digest() -> tuple[str, dict]:
    """The digest lives OUTSIDE the verdict JSON it covers -- otherwise it would hash itself."""
    per = {}
    for rel in CLAIM_LOCK_FILES + ["results/claim_lock/GEN1_CLAIM_LOCK.json"]:
        p = ROOT / rel
        per[rel] = stable_sha256(p) if p.exists() else "MISSING"
    canonical = "\n".join(f"{k}  {per[k]}" for k in sorted(per))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest(), per


# =============================================================================================== #
# The mandatory scope qualifiers (plan §2.1)
#
# A claim quoted without its qualifier is a claim the lock did not grant. These are not decoration;
# they are the difference between what was measured and what a reader will assume was measured.
# =============================================================================================== #
QUALIFIERS = {
    "system": "one BRAF-V600E melanoma cell line, WM989 (GSE279162), 1,401 lineage-traced clones",
    "vocabulary": "six observed experimental conditions: Acid, Cisplatin, CoCl2, Dabrafenib, "
                  "Doxorubicin, Trametinib",
    "outcome": "an observed post-treatment clone DETECTION proxy -- not death, sensitivity, "
               "resistance, clinical response or patient benefit",
    "evaluation": "clone-held-out, five outer folds, frozen before any result existed",
    "replication": "NONE. Independent biological replication is Generation 2 and was not "
                   "performed.",
}


# =============================================================================================== #
# CL-D — the adversarial abstract corpus (plan §4)
#
# Declared before it is run, so it cannot be trimmed after seeing which entries the scanner misses.
# A missed entry is a FAILED STAGE: it means the instrument cannot see a claim the manuscript could
# make. Each carries the nearest permitted phrasing, because a writer told only "no" learns nothing.
# =============================================================================================== #
ADVERSARIAL_CLAIMS: list[tuple[str, str, str]] = [
    ("generalisation", "The model generalises to new treatments.",
     "The model is supported for the six observed experimental conditions; anything else returns "
     "UNSUPPORTED_TREATMENT."),
    ("generalisation", "Our approach works for any drug.",
     "Our approach was evaluated on six observed experimental conditions in one cell line."),
    ("cross_system", "We predict outcomes in cancer cells.",
     "We predict a clone-detection outcome in one BRAF-V600E melanoma cell line, WM989."),
    ("cross_system", "Validated on other cell lines.",
     "Evaluated clone-held-out within WM989; no other cell line was tested."),
    ("cross_system", "The model predicts patient response.",
     "The model predicts an observed post-treatment clone-detection proxy, not clinical response."),
    ("clinical", "The tool identifies the best treatment for each clone.",
     "The tool orders six observed experimental conditions for each clone; this is "
     "experimental-condition selection, not treatment recommendation."),
    ("clinical", "This supports clinical decision-making.",
     "This is a research tool for the WM989 experiment and is not a clinical tool."),
    ("clinical", "A therapeutic recommendation tool for melanoma.",
     "A research predictor for six observed conditions in one melanoma cell line."),
    ("causal", "We estimate the causal effect of each treatment.",
     "We estimate a predictive association between pretreatment state and a later observed "
     "outcome, under clone-held-out evaluation."),
    ("calibration", "The model outputs a calibrated probability of death.",
     "The model outputs an uncalibrated future_detection_score; no calibration was frozen or "
     "passed."),
    ("replication", "Independently replicated in an external cohort.",
     "No independent biological replication was performed; clone-held-out folds and two endpoint "
     "families are not replication."),
    ("uniformity", "The interaction helps uniformly across all six conditions.",
     "Four of six conditions carry meaningful interaction; Cisplatin is negligible on C1 and "
     "Doxorubicin is negative on both endpoints."),
    ("role_a", "The finding was confirmed in a second system.",
     "A separately reconstructed reprogramming system gave positive but underpowered supporting "
     "evidence; its confirmation gate FAILED."),
    ("role_a", "Rewind confirms the result.",
     "Rewind supports the result without confirming it; gate 18.3 FAILED at 0.64, audited ~0.45."),
    ("single_cell", "The tool scores a single cell.",
     "The tool scores a clone-level pseudobulk profile; a single cell is not an equivalent input."),
]


# =============================================================================================== #
# CL-D.2 — the prose extension (plan §4.2)
#
# The Stage-26 scanner was calibrated on the tool's own surfaces: code, a model card, a JSON
# schema. The adversarial corpus above shows that is not enough for PROSE. Fired at the unmodified
# Stage-26 patterns, three of fifteen forbidden sentences walked straight through:
#
#   "We predict outcomes in cancer cells."         `cancer` is REVIEW_ONLY in Stage 26, not gating
#   "The finding was confirmed in a second system." no pattern for confirmation without the word
#                                                   Rewind or Role-A next to it
#   "The tool scores a single cell."                the pattern reads `score a single cell`; the
#                                                   sentence says `scores`
#
# The tool's documentation would never say "cancer cells". An abstract easily would.
#
# These patterns are ADDED here, never substituted, and `run_stage26_scope_lock.py` is not touched:
# it is a locked evidence artifact, and §3 of this plan permits adding to the forbidden set but
# never subtracting. The combined set is then turned back on the already-locked shipped surfaces,
# because a stricter instrument that is only ever pointed at new text is not an instrument.
# =============================================================================================== #
PROSE_PATTERNS: dict[str, list[str]] = {
    "2_cross_cell_line_or_patient": [
        r"\bcancer\b", r"\btumou?r\b", r"a second cell line", r"across cell lines",
        r"melanoma patients?"],
    "8_confirmed_role_a": [
        r"confirmed in a second", r"confirmed in another", r"confirmed in an independent",
        r"a second system", r"independent confirmation", r"confirms the (?:result|finding)"],
    "9_single_cell_equivalence": [
        r"scores? a single cell", r"scoring a single cell", r"single[- ]cell (?:input|scoring)",
        r"per[- ]cell (?:score|prediction)"],
}


def combined_patterns() -> dict[str, list[str]]:
    return {k: v + PROSE_PATTERNS.get(k, []) for k, v in S26.FORBIDDEN_CLAIMS.items()}


# Stage 26 excuses a forbidden phrase when a negation token appears anywhere within +-160
# characters. On code and a model card that is fine. On PROSE it is nearly toothless, because prose
# is full of legitimate negations -- and a proximity rule cannot tell which clause they govern:
#
#   "The model is not calibrated for abundance, and outputs a calibrated probability of death."
#   "We make no claim about dosing; the tool identifies the best treatment for each clone."
#   "This was not replicated internally, but was independently replicated in an external cohort."
#
# All three make a plainly forbidden claim. All three went undetected under the window rule.
#
# So prose is scanned CLAUSE-SCOPED: a negation excuses a hit only when it sits in the same clause.
# Same twelve tokens, same nine claims -- a tighter scope, not a different instrument.
# A newline is NOT a clause boundary. Treating it as one split "...and not a\nclinical
# recommendation." in the shipped predictor's docstring and reported a negated sentence as a
# forbidden claim. Whitespace is normalised first so the rule is immune to line wrapping.
CLAUSE_SPLIT = re.compile(r"[.;:]|,\s+(?:and|but|while|whereas|yet|although|though)\b")


# A word that negates itself. `uncalibrated` contains `calibrated`, and the Stage-26 pattern has no
# word boundary -- so "outputs an uncalibrated score" was flagged as claiming a calibrated
# probability. Excusing by prefix is narrower than adding a negation token, which would loosen the
# gate everywhere; the prefix must be contiguous with the match, so "run calibrated" is untouched.
NEGATING_PREFIXES = ("un", "non", "non-", "de")


def prose_scan(text: str, patterns: dict[str, list[str]]) -> list[dict]:
    """Unnegated hits, where 'negated' means negated IN THE SAME CLAUSE."""
    out = []
    for clause in CLAUSE_SPLIT.split(" ".join(text.split())):
        low = clause.lower()
        if any(n in low for n in S26.NEGATIONS):
            continue
        for claim, pats in patterns.items():
            for pat in pats:
                for m in re.finditer(pat, low):
                    if low[: m.start()].endswith(NEGATING_PREFIXES):
                        continue
                    out.append({"claim": claim, "match": clause[m.start():m.end()],
                                "clause": " ".join(clause.split())[:200]})
    return out


def scan(text: str, patterns: dict[str, list[str]]) -> list[dict]:
    """Prose scan. Kept as the single entry point so every caller gets the tighter rule."""
    return prose_scan(text, patterns)


def window_scan(text: str, patterns: dict[str, list[str]]) -> list[dict]:
    """Stage 26's own rule, unchanged -- used to show what the looser rule would have allowed."""
    out = []
    for claim, pats in patterns.items():
        for h in S26._scan_text(text, pats):
            if not h["negated"]:
                out.append({"claim": claim, "match": h["match"], "context": h["context"]})
    return out


# =============================================================================================== #
# CL-A — the evidence must verify before a sentence is written (plan §1)
# =============================================================================================== #
def verify_evidence() -> dict:
    t0 = time.perf_counter()
    handoff = _j(EVIDENCE_HANDOFF)
    manifest = _j(EVIDENCE_MANIFEST)
    live = EL.verify_against(manifest, ROOT)
    expected = handoff["lock_digest"]

    plan_text = PLAN.read_text(encoding="utf-8")
    checks = {
        "the evidence lock verdict is GEN1_EVIDENCE_LOCKED":
            handoff["verdict"] == "GEN1_EVIDENCE_LOCKED",
        f"all {manifest['n_artifacts']} locked artifacts still verify": live["clean"],
        "the lock digest is the one this plan was written against":
            manifest["lock_digest"] == expected and expected in plan_text,
        "the frozen ship-plan digest still holds":
            EL.canonical_lf_sha256(SHIP_PLAN) == SHIP_PLAN_DIGEST,
    }
    return {"stage": "CL-A", "lock_digest": manifest["lock_digest"],
            "live_verification": live, "checks": checks, "all_passed": all(checks.values()),
            "runtime_seconds": round(time.perf_counter() - t0, 3)}


# =============================================================================================== #
# CL-B / CL-C — the allowed claims, and the forbidden nine unchanged (plan §2, §3)
# =============================================================================================== #
def build_claims() -> dict:
    t0 = time.perf_counter()
    src = _j(CLAIM_INPUT)
    v25 = _j(RESULTS / "stage25" / "stage25_verdict.json")
    lock = _j(RESULTS / "evidence_lock" / "GEN1_EVIDENCE_LOCK.json")
    h = lock["headline_numbers"]

    allowed = {
        "primary": {
            "text": src["allowed"]["primary"],
            "source": "GEN1_CLAIM_LOCK_INPUT.json, verbatim (§3.1 of the frozen ship plan)",
            "evidence": ["results/stage22_wm989_clones.csv",
                         "results/stage23_wm989_interaction_oof.csv",
                         "results/stage24/stage24_oof_for_stage25.csv"],
            "numbers": "none at abstract level; the supporting figures are the frozen "
                       "out-of-fold results",
            "qualifiers": ["system", "vocabulary", "outcome", "evaluation", "replication"],
        },
        "ranking": {
            "text": src["allowed"]["ranking"],
            "source": f"selected by {src['allowed']['ranking_selected_by']}",
            "evidence": ["results/stage25/stage25_verdict.json",
                         "results/stage24/stage24_oof_for_stage25.csv"],
            "numbers": {"delta_RANK": h["delta_RANK"], "ci95": h["bootstrap_ci95"],
                        "permutation": f"{h['n_null_ge_observed']} of {h['n_perm']} draws "
                                       f"reached the observed value",
                        "p_perm": "p < 0.001 (0 of 1,000); never a point estimate",
                        "eligible_clones": h["eligible_clones"]},
            "qualifiers": ["system", "vocabulary", "outcome", "evaluation", "replication"],
        },
        "supporting_role_A": {
            "text": src["allowed"]["supporting_role_A"],
            "source": "GEN1_CLAIM_LOCK_INPUT.json, verbatim (§3.3 of the frozen ship plan)",
            "evidence": ["results/stage23_2h/stage23_2h_verdict.json",
                         "results/stage23_2h/stage23_2h_power_audit.json"],
            "numbers": "confirmation gate 18.3 FAILED at 0.64; audited to ~0.45",
            "qualifiers": ["replication"],
            "must_travel_with": "the word SUPPORTING. Rewind does not confirm anything, and its "
                                "own confirmation gate failed.",
        },
    }

    # ---- CL-C: the nine, parsed from the frozen plan rather than from a copy ------------------ #
    ship = SHIP_PLAN.read_text(encoding="utf-8")
    section = ship.split("## 3.5 Claims forbidden in Generation 1")[1].split("```")[1]
    in_plan = [ln.strip() for ln in section.strip().splitlines()[1:] if ln.strip()]

    # ---- the narrowing check: nothing may be widened ------------------------------------------ #
    widened = [k for k, v in allowed.items()
               if v["text"] not in (src["allowed"]["primary"], src["allowed"]["ranking"],
                                    src["allowed"]["supporting_role_A"])]

    checks = {
        "every allowed claim is verbatim from the evidence lock's input": not widened,
        "the nine forbidden claims match the frozen plan exactly": src["forbidden"] == in_plan,
        "there are exactly nine, and none was dropped": len(src["forbidden"]) == 9,
        "every allowed claim names the evidence that supports it":
            all(v["evidence"] for v in allowed.values()),
        "every allowed claim carries mandatory scope qualifiers":
            all(v["qualifiers"] for v in allowed.values()),
        "the ranking claim matches the Stage-25 verdict":
            src["allowed"]["ranking_selected_by"] == v25["verdict"],
        "replication is recorded as Generation 2":
            src["replication"]["status"] == "GENERATION 2",
    }
    return {"stage": "CL-B/CL-C", "allowed": allowed, "forbidden": src["forbidden"],
            "forbidden_source": "§3.5 of the frozen ship plan, parsed from the plan",
            "qualifiers": QUALIFIERS,
            "widened": widened,
            "checks": checks, "all_passed": all(checks.values()),
            "runtime_seconds": round(time.perf_counter() - t0, 3)}


# =============================================================================================== #
# CL-D — fire the adversarial corpus at the Stage-26 scanner (plan §4)
# =============================================================================================== #
def adversarial_corpus() -> dict:
    t0 = time.perf_counter()
    canary = S26._canary()
    base, full = S26.FORBIDDEN_CLAIMS, combined_patterns()

    results, missed, missed_by_base = [], [], []
    for group, sentence, permitted in ADVERSARIAL_CLAIMS:
        hits = scan(sentence, full)
        if not scan(sentence, base):
            missed_by_base.append({"group": group, "sentence": sentence})
        if not hits:
            missed.append({"group": group, "sentence": sentence})
        results.append({"group": group, "forbidden_sentence": sentence,
                        "caught": bool(hits), "triggered": sorted({h["claim"] for h in hits}),
                        "caught_by_stage26_alone": bool(scan(sentence, base)),
                        "nearest_permitted_phrasing": permitted})

    # the permitted neighbours must themselves come back clean -- otherwise the table teaches a
    # writer to replace one refused sentence with another
    dirty = [{"permitted": r["nearest_permitted_phrasing"], **h}
             for r in results for h in scan(r["nearest_permitted_phrasing"], full)]

    # A stricter instrument that is only ever pointed at new text is not an instrument. Turn the
    # extended patterns back on the surfaces Stage 26 already passed.
    #
    # Those surfaces are scanned with Stage 26's WINDOW rule, not the prose rule, and that is a
    # judgement worth stating rather than burying. Clause scoping is right for prose and wrong for
    # structured text: it splits a JSON key from the value that negates it, so
    # `"calibrated_probability": "NEVER emitted in Generation 1"` reads as a forbidden claim. The
    # patterns are the extended ones either way; only the negation SCOPE differs, matched to the
    # kind of text being read. Both readings are recorded below.
    resurvey, resurvey_clause = {}, {}
    for f in S26._surfaces():
        if not f.exists():
            continue
        text = f.read_text(encoding="utf-8")
        rel = f.relative_to(ROOT).as_posix()
        if hits := window_scan(text, full):
            resurvey[rel] = hits
        if hits := prose_scan(text, full):
            resurvey_clause[rel] = hits

    # the extension must not have broken the negation rule
    neg_ok = all(not scan(f"This tool cannot {s[0].lower()}{s[1:]}", full)
                 for _g, s, _p in ADVERSARIAL_CLAIMS)

    checks = {
        "the scanner detects all nine claims stated plainly (Stage-26 canary)":
            canary["detects_every_forbidden_claim"],
        "the scanner does not fire on a negation (Stage-26 canary)":
            canary["does_not_fire_on_a_negation"],
        "the prose extension still does not fire on a negation": neg_ok,
        "every adversarial abstract sentence is caught": not missed,
        "every nearest-permitted phrasing comes back clean": not dirty,
        "the extended patterns find nothing on the already-locked shipped surfaces": not resurvey,
    }
    return write_json(ADVERSARIAL_JSON, {
        "stage": "CL-D",
        "instrument": "the Stage-26 scope-lock scanner plus a prose extension; the Stage-26 "
                      "module itself is a locked artifact and was not modified",
        "n_sentences": len(ADVERSARIAL_CLAIMS),
        "n_caught": sum(1 for r in results if r["caught"]),
        "missed": missed,
        "missed_by_stage26_patterns_alone": missed_by_base,
        "why_the_extension_exists":
            "Stage 26's patterns were calibrated on the tool's own surfaces -- code, a model card, "
            "a JSON schema. Fired at prose they let three of fifteen forbidden sentences through. "
            "The tool's documentation would never say 'cancer cells'; an abstract easily would.",
        "prose_patterns_added": PROSE_PATTERNS,
        "permitted_phrasings_that_failed_their_own_scan": dirty,
        "resurvey_of_locked_surfaces": resurvey,
        "resurvey_under_clause_scoping": resurvey_clause,
        "resurvey_note":
            "The locked surfaces are gated under Stage 26's window rule with the EXTENDED "
            "patterns. Clause scoping is reported alongside but not gated on them: it is right "
            "for prose and wrong for structured text, where it splits a JSON key from the value "
            "that negates it. Every clause-scoped hit below was read and is a false positive of "
            "that kind.",
        "stage26_canary": canary,
        "near_miss_table": results,
        "checks": checks, "all_passed": all(checks.values()),
        "runtime_seconds": round(time.perf_counter() - t0, 3),
    })


# =============================================================================================== #
# CL-E — one worked abstract, scanned by the same instrument (plan §5)
# =============================================================================================== #
def worked_abstract(claims: dict) -> dict:
    t0 = time.perf_counter()
    h = _j(RESULTS / "evidence_lock" / "GEN1_EVIDENCE_LOCK.json")["headline_numbers"]
    text = (
        "In one BRAF-V600E melanoma cell line (WM989, 1,401 lineage-traced clones), pretreatment "
        "gene expression carries treatment-specific information about a later observed clonal "
        "detection outcome, beyond treatment identity and captured pretreatment clone abundance, "
        "under clone-held-out evaluation frozen before any result existed. Under a preregistered "
        "test, a frozen state-by-treatment interaction model improves clone-specific ordering of "
        f"the six observed experimental conditions over a non-interactive additive model "
        f"(delta {h['delta_RANK']:+.4f} in equal-clone-weighted within-clone AUROC, 95% CI "
        f"[{h['bootstrap_ci95'][0]:+.4f}, {h['bootstrap_ci95'][1]:+.4f}]; no draw of 1,000 "
        "full-refit permutations reached the observed value, p < 0.001). The outcome is an "
        "observed post-treatment clone-detection proxy and is not death, sensitivity, resistance "
        "or clinical response. The six conditions are the entire supported vocabulary; the model "
        "makes no claim about unseen treatments, other cell lines, or patients, and emits no "
        "calibrated probability. Independent biological replication was not performed and remains "
        "future work."
    )
    hits = scan(text, combined_patterns())

    qualifiers_present = {
        "system": "WM989" in text and "melanoma" in text,
        "vocabulary": "six observed experimental conditions" in text,
        "outcome": "not death" in text,
        "evaluation": "clone-held-out" in text,
        "replication": "replication was not performed" in text,
    }
    checks = {
        "the worked abstract passes the same scanner": not hits,
        "every mandatory qualifier is present": all(qualifiers_present.values()),
        "it never says treatment where it means condition":
            "treatments" not in text.replace("unseen treatments", ""),
    }
    return {"stage": "CL-E", "abstract": text, "unnegated_hits": hits,
            "qualifiers_present": qualifiers_present,
            "status": "a demonstration of the ceiling, not a mandated abstract",
            "checks": checks, "all_passed": all(checks.values()),
            "runtime_seconds": round(time.perf_counter() - t0, 3)}


# =============================================================================================== #
def _claims_document(a: dict, c: dict, d: dict, e: dict, verdict: str) -> str:
    allowed = c["allowed"]
    near = "\n\n".join(
        f"  FORBIDDEN   {r['forbidden_sentence']}\n  PERMITTED   {r['nearest_permitted_phrasing']}"
        for r in d["near_miss_table"])
    return f"""# CellFate-Rx Generation 1 — CLAIM LOCK

```text
  {verdict}

  evidence lock digest   {a['lock_digest']}
  allowed claims         {len(allowed)}
  forbidden claims       {len(c['forbidden'])}
  adversarial sentences  {d['n_caught']} of {d['n_sentences']} caught
```

Every previous stage produced a number. This one fixes **sentences** — because that is where a
project of this kind actually fails: not in the statistics, but in the abstract, where "six
observed experimental conditions in one melanoma line" becomes "treatments in cancer".

## The ceiling may be lowered. It may not be raised.

---

## Allowed claims

### Primary
> {allowed['primary']['text']}

### Ranking — {allowed['ranking']['source']}
> {allowed['ranking']['text']}

```text
  delta_RANK  {allowed['ranking']['numbers']['delta_RANK']:+.6f}
  CI95        [{allowed['ranking']['numbers']['ci95'][0]:+.6f}, {allowed['ranking']['numbers']['ci95'][1]:+.6f}]
  null        {allowed['ranking']['numbers']['permutation']}
  p           {allowed['ranking']['numbers']['p_perm']}
  clones      {allowed['ranking']['numbers']['eligible_clones']} eligible
```

### Supporting — Role A
> {allowed['supporting_role_A']['text']}

**This claim must travel with** {allowed['supporting_role_A']['must_travel_with']}
Its own evidence: {allowed['supporting_role_A']['numbers']}.

---

## Mandatory qualifiers

A claim quoted without its qualifier is a claim this lock did not grant.

```text
  system       {QUALIFIERS['system']}
  vocabulary   {QUALIFIERS['vocabulary']}
  outcome      {QUALIFIERS['outcome']}
  evaluation   {QUALIFIERS['evaluation']}
  replication  {QUALIFIERS['replication']}
```

---

## Forbidden claims — the nine, unchanged

```text
{chr(10).join('  ' + str(i + 1) + '  NEVER  ' + s for i, s in enumerate(c['forbidden']))}
```

Parsed from §3.5 of the frozen ship plan, not from a copy. This lock may add to the list. It may
not subtract from it or reword an entry.

---

## Where the boundary is

A writer told only "no" learns nothing. Each forbidden sentence below was fired at the same scanner
Stage 26 used on the shipped tool, and every one was caught. Each is paired with the nearest
phrasing that is permitted — and every permitted phrasing was itself scanned clean, so this table
never teaches one refused sentence to be swapped for another.

```text
{near}
```

---

## A permitted abstract

Assembled only from locked claims and their mandatory qualifiers, and scanned by the same
instrument. {e['status'].capitalize()}.

> {e['abstract']}

---

## What locking a claim does not do

It grants nothing. It fixes the ceiling of what may be said about evidence that was locked
separately. No claim-lock outcome reopens an earlier stage, changes a recorded number, or
authorizes new data, a new condition or a new model.
"""


def run_all() -> dict:
    t0 = time.perf_counter()
    a = verify_evidence()
    if not a["all_passed"]:
        out = write_json(CLAIMS_JSON, {
            "stage": "GEN-1 CLAIM LOCK", "verdict": "GEN1_CLAIM_LOCK_REFUSED",
            "refused_at": "CL-A", "detail": a,
            "note": "The evidence lock did not verify. Writing claims against evidence that has "
                    "shifted is the failure the previous stage existed to prevent."})
        return out

    c = build_claims()
    d = adversarial_corpus()
    e = worked_abstract(c)

    sub = {"CL-A_evidence_verifies": a["all_passed"], "CL-B_C_claims": c["all_passed"],
           "CL-D_adversarial": d["all_passed"], "CL-E_worked_abstract": e["all_passed"]}
    verdict = "GEN1_CLAIMS_LOCKED" if all(sub.values()) else "GEN1_CLAIM_LOCK_REFUSED"

    out = write_json(CLAIMS_JSON, {
        "stage": "GEN-1 CLAIM LOCK",
        "verdict": verdict,
        "substages": sub,
        "failing": [k for k, v in sub.items() if not v],
        "evidence_lock_digest": a["lock_digest"],
        "plan": PLAN.relative_to(ROOT).as_posix(),
        "plan_canonical_lf_sha256": EL.canonical_lf_sha256(PLAN),
        "ship_plan_digest": SHIP_PLAN_DIGEST,
        "allowed": c["allowed"],
        "forbidden": c["forbidden"],
        "qualifiers": QUALIFIERS,
        "adversarial": {"n": d["n_sentences"], "caught": d["n_caught"], "missed": d["missed"]},
        "worked_abstract": e,
        "evidence_verification": a,
        "grants_nothing": "Locking a claim grants nothing. It fixes the ceiling of what may be "
                          "said. The ceiling may be lowered later; it may not be raised.",
        "next": "MANUSCRIPT + REPRODUCIBILITY PACKAGE",
        "runtime_seconds": round(time.perf_counter() - t0, 3),
    })
    CLAIMS_MD.write_text(_claims_document(a, c, d, e, verdict), encoding="utf-8")

    digest, per_file = claim_digest()
    write_json(DIGEST_JSON, {
        "claim_digest": digest,
        "digest_definition": "SHA-256 over 'path  content-sha256' lines, sorted by path, "
                             "LF-joined. JSON members are normalised -- volatile timing fields "
                             "removed, keys sorted -- so the digest reflects content rather than "
                             "when the stage last ran; other text uses canonical-LF.",
        "volatile_keys_excluded": sorted(VOLATILE_KEYS),
        "covers": per_file,
        "why_it_lives_outside_the_verdict": "a digest stored inside the file it covers would hash "
                                            "itself",
        "verify": "python experiments/run_gen1_claim_lock.py --verify",
    })

    out["claim_digest"] = digest
    write_json(HANDOFF_JSON, {
        "claim_digest": digest,
        "from_stage": "GEN-1 CLAIM LOCK",
        "to_stage": "MANUSCRIPT + REPRODUCIBILITY PACKAGE",
        "verdict": verdict,
        "evidence_lock_digest": a["lock_digest"],
        "claims_document": CLAIMS_MD.relative_to(ROOT).as_posix(),
        "manuscript_must": [
            "re-verify the evidence lock before submission: "
            "python experiments/run_gen1_evidence_lock.py --verify",
            "carry each allowed claim with its mandatory qualifiers; a claim quoted without "
            "them is a claim this lock did not grant",
            "report p_perm as p < 0.001 (0 of 1,000), never as a point estimate",
            "state that independent biological replication is Generation 2, not a Gen-1 gate",
            "scan any new abstract-level sentence against the same nine forbidden claims",
            "bind to BOTH digests: the evidence digest for what the claims are made of, and "
            "the claim digest for what may be said about it"],
        "the_ceiling_may_be_lowered_not_raised": True,
        "no_claim_lock_outcome_reopens_an_earlier_stage": True,
    })
    return out


def run_verify() -> dict:
    """Re-hash the claim files and refuse if one moved -- or if the stage itself refused.

    A verifier that only re-hashes bytes answers 'has this been tampered with?' and NOT 'did
    the stage that produced this pass?'. Those came apart in practice: the manuscript stage
    recorded GEN1_MANUSCRIPT_REFUSED, its covered files hashed exactly as recorded -- because a
    refused run still writes them -- and --verify returned clean, so CI stayed green over a
    refused package. The recorded verdict is therefore part of what is verified.
    """
    if not DIGEST_JSON.exists():
        raise SystemExit("no claim digest: run --stage all first")
    recorded = _j(DIGEST_JSON)
    now, per = claim_digest()
    moved = [k for k, v in per.items() if recorded["covers"].get(k) != v]
    bytes_intact = not moved and now == recorded["claim_digest"]
    stage = _j(CLAIMS_JSON).get("verdict") if CLAIMS_JSON.exists() else None
    passed = stage == "GEN1_CLAIMS_LOCKED"
    return {"clean": bytes_intact and passed, "moved": moved, "claim_digest": now,
            "recorded_digest": recorded["claim_digest"],
            "stage_verdict": stage, "stage_passed": passed,
            "verdict": ("CLAIMS_MOVED" if not bytes_intact
                        else "CLAIMS_INTACT" if passed else "CLAIMS_STAGE_REFUSED")}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Generation-1 claim lock")
    ap.add_argument("--stage", choices=["all"], default=None)
    ap.add_argument("--verify", action="store_true",
                    help="re-hash the claim-lock files and refuse if one has moved")
    a = ap.parse_args(argv)
    if a.verify:
        r = run_verify()
        print(json.dumps(r, indent=2))
        return 0 if r["clean"] else 2
    if a.stage != "all":
        ap.error("pass --stage all or --verify")
    r = run_all()
    payload = {k: r[k] for k in r if k in
               ("stage", "verdict", "substages", "failing", "evidence_lock_digest",
                "claim_digest", "adversarial", "refused_at", "next")}
    if r["verdict"] != "GEN1_CLAIMS_LOCKED":
        # a refused run must not advertise the next stage: it announced "Generation 1 is
        # complete" while its own verdict was REFUSED.
        payload.pop("next", None)
    print(json.dumps(payload, indent=2, default=str))
    return 0 if r["verdict"] == "GEN1_CLAIMS_LOCKED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
