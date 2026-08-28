"""Contracts for the Generation-1 claim lock.

Every earlier stage produced a number. This one produces sentences, which is where a project of
this kind actually fails -- in the abstract, where "six observed experimental conditions in one
melanoma line" becomes "treatments in cancer".

The load-bearing contract is `test_every_forbidden_abstract_sentence_is_caught`. Listing allowed
sentences proves nothing about an abstract assembled from permitted fragments; firing tempting
forbidden sentences at the scanner does. That corpus already earned its place: three of its fifteen
entries walked straight through the unmodified Stage-26 patterns.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = RESULTS / "claim_lock"
SRC = ROOT / "experiments" / "run_gen1_claim_lock.py"
PLAN = ROOT / "plans" / "(newer)practical plans" / "GEN1_CLAIM_LOCK_V1.md"
SHIP_PLAN = (ROOT / "plans" / "(newer)practical plans"
             / "STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md")

CLAIMS = OUT / "GEN1_CLAIM_LOCK.json"
CLAIMS_MD = OUT / "GEN1_CLAIMS.md"
ADVERSARIAL = OUT / "claim_lock_adversarial.json"
HANDOFF = RESULTS / "gen1_handoff_to_manuscript.json"

ran = pytest.mark.skipif(not CLAIMS.exists(), reason="the claim lock has not been run")


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def mod():
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    import run_gen1_claim_lock as CL
    return CL


# ============================================================================================== #
# The corpus is the point
# ============================================================================================== #
def test_the_adversarial_corpus_is_the_plan_s(mod):
    """Declared before it was run, so it cannot be trimmed after seeing what the scanner misses.

    Every sentence, verbatim -- not three sampled fragments. An earlier version of this test
    checked three and passed while seven of the fifteen appeared in the plan only as paraphrases.
    """
    plan = PLAN.read_text(encoding="utf-8")
    groups = {g for g, _s, _p in mod.ADVERSARIAL_CLAIMS}
    assert groups == {"generalisation", "cross_system", "clinical", "causal", "calibration",
                      "replication", "uniformity", "role_a", "single_cell"}
    assert len(mod.ADVERSARIAL_CLAIMS) == 15

    for group, sentence, _p in mod.ADVERSARIAL_CLAIMS:
        assert f'"{sentence}"' in plan, f"{sentence!r} is run but not declared verbatim"
        assert group.upper() in plan, group

    sentences = [s for _g, s, _p in mod.ADVERSARIAL_CLAIMS]
    assert len(sentences) == len(set(sentences))
    # the three that exposed the Stage-26 gap must still be in the corpus
    for fragment in ("cancer cells", "confirmed in a second system", "scores a single cell"):
        assert any(fragment in s for s in sentences), fragment


@ran
def test_every_forbidden_abstract_sentence_is_caught():
    a = _json(ADVERSARIAL)
    assert a["missed"] == []
    assert a["n_caught"] == a["n_sentences"]
    for row in a["near_miss_table"]:
        assert row["caught"] is True, row["forbidden_sentence"]
        assert row["triggered"], row["forbidden_sentence"]


@ran
def test_the_near_miss_table_never_swaps_one_refusal_for_another():
    """Every permitted neighbour is itself scanned clean, or the table teaches bad substitutions."""
    a = _json(ADVERSARIAL)
    assert a["permitted_phrasings_that_failed_their_own_scan"] == []
    for row in a["near_miss_table"]:
        assert row["nearest_permitted_phrasing"], row["forbidden_sentence"]


def test_a_forbidden_sentence_is_rejected_and_its_neighbour_is_not(mod):
    """Directly, not by reading a JSON that says so."""
    full = mod.combined_patterns()
    assert mod.scan("We predict outcomes in cancer cells.", full)
    assert not mod.scan(
        "We predict a clone-detection outcome in one BRAF-V600E melanoma cell line, WM989.", full)


# ============================================================================================== #
# The prose extension
# ============================================================================================== #
@ran
def test_the_stage_26_scanner_alone_was_not_enough_for_prose():
    """Recorded rather than quietly fixed: three of fifteen went through the unmodified patterns."""
    a = _json(ADVERSARIAL)
    missed = [m["sentence"] for m in a["missed_by_stage26_patterns_alone"]]
    assert len(missed) == 3
    assert any("cancer cells" in s for s in missed)
    assert any("second system" in s for s in missed)
    assert any("single cell" in s for s in missed)


def test_the_extension_adds_and_never_substitutes(mod):
    """The nine claims and the twelve negation tokens are Stage 26's, unchanged."""
    import run_stage26_scope_lock as S26
    full = mod.combined_patterns()
    assert set(full) == set(S26.FORBIDDEN_CLAIMS)
    for claim, pats in S26.FORBIDDEN_CLAIMS.items():
        assert pats == full[claim][:len(pats)], f"{claim}: Stage-26 patterns must survive intact"
    assert len(S26.NEGATIONS) == 12


def test_the_locked_stage_26_module_is_not_edited_by_this_stage():
    """It is a locked evidence artifact; extending it in place would break the evidence lock."""
    manifest = _json(RESULTS / "evidence_lock" / "GEN1_EVIDENCE_MANIFEST.json")
    assert "experiments/run_stage26_scope_lock.py" in manifest["artifacts"]
    src = SRC.read_text(encoding="utf-8")
    assert "ADDED here, never substituted" in src
    assert "is not touched" in src


@ran
def test_the_stricter_patterns_were_turned_back_on_the_old_surfaces():
    """A stricter instrument pointed only at new text is not an instrument."""
    a = _json(ADVERSARIAL)
    assert "resurvey_of_locked_surfaces" in a
    assert a["resurvey_of_locked_surfaces"] == {}, \
        "the extended patterns found a forbidden claim in an already-locked surface"


@ran
def test_the_extension_still_refuses_to_fire_on_a_negation():
    a = _json(ADVERSARIAL)
    assert a["checks"]["the prose extension still does not fire on a negation"] is True
    assert a["stage26_canary"]["does_not_fire_on_a_negation"] is True


# ============================================================================================== #
# Negation is scoped to the clause, not to a window
# ============================================================================================== #
def test_a_negation_in_another_clause_does_not_excuse_a_forbidden_claim(mod):
    """The window rule passes all three of these. Prose is full of legitimate negations, and
    proximity cannot tell which clause they govern."""
    full = mod.combined_patterns()
    for sentence in (
        "The model is not calibrated for abundance, and outputs a calibrated probability of death.",
        "We make no claim about dosing; the tool identifies the best treatment for each clone.",
        "This was not replicated internally, but was independently replicated in an external "
        "cohort.",
    ):
        assert mod.window_scan(sentence, full) == [], "these are the cases the window rule misses"
        assert mod.scan(sentence, full), f"clause scoping must catch: {sentence}"


def test_a_negation_in_the_same_clause_still_excuses(mod):
    full = mod.combined_patterns()
    for sentence in ("This is not a clinical tool.",
                     "The model makes no claim about unseen treatments or other cell lines.",
                     "No independent biological replication was performed."):
        assert mod.scan(sentence, full) == [], sentence


def test_a_line_wrap_is_not_a_clause_boundary(mod):
    """Treating \\n as a boundary orphaned "not a / clinical recommendation" in the shipped
    predictor's docstring and reported a negated sentence as a forbidden claim."""
    full = mod.combined_patterns()
    wrapped = "The score is not a calibrated probability, not a measure of death, and not a\n" \
              "clinical recommendation."
    assert mod.scan(wrapped, full) == []


def test_a_word_may_negate_itself(mod):
    """`uncalibrated` contains `calibrated`, and the Stage-26 pattern has no word boundary."""
    full = mod.combined_patterns()
    assert mod.scan("The model outputs an uncalibrated score", full) == []
    # but the prefix must be contiguous -- this must still be caught
    assert mod.scan("We run calibrated probability estimates for each clone", full)


# ============================================================================================== #
# The claims must have an identity of their own
# ============================================================================================== #
@ran
def test_the_claim_lock_has_its_own_digest():
    """The evidence lock hashes 54 artifacts and none of them is this stage's output."""
    d = _json(OUT / "GEN1_CLAIM_DIGEST.json")
    assert len(d["claim_digest"]) == 64
    assert set(d["covers"]) == {
        "plans/(newer)practical plans/GEN1_CLAIM_LOCK_V1.md",
        "experiments/run_gen1_claim_lock.py",
        "tests/test_gen1_claim_lock.py",
        "results/claim_lock/GEN1_CLAIMS.md",
        "results/claim_lock/GEN1_CLAIM_LOCK.json"}
    assert "MISSING" not in d["covers"].values()
    assert _json(HANDOFF)["claim_digest"] == d["claim_digest"]


@ran
def test_the_claim_digest_is_reproducible_and_lives_outside_what_it_covers(mod):
    import hashlib
    d = _json(OUT / "GEN1_CLAIM_DIGEST.json")
    canonical = "\n".join(f"{k}  {d['covers'][k]}" for k in sorted(d["covers"]))
    assert hashlib.sha256(canonical.encode("utf-8")).hexdigest() == d["claim_digest"]
    # the verdict JSON must NOT contain the digest that covers it
    assert "claim_digest" not in CLAIMS.read_text(encoding="utf-8"), \
        "a digest stored inside the file it covers would hash itself"


@ran
def test_the_manuscript_binds_to_both_digests():
    h = _json(HANDOFF)
    assert h["evidence_lock_digest"] and h["claim_digest"]
    assert "bind to BOTH digests" in " ".join(h["manuscript_must"])


# ============================================================================================== #
# Evidence, claims, qualifiers
# ============================================================================================== #
@ran
def test_the_claim_lock_refuses_to_write_against_unverified_evidence():
    a = _json(CLAIMS)["evidence_verification"]
    assert a["live_verification"]["clean"] is True
    assert a["checks"]["all 54 locked artifacts still verify"] is True
    assert _json(CLAIMS)["evidence_lock_digest"] == \
        _json(RESULTS / "evidence_lock" / "GEN1_EVIDENCE_MANIFEST.json")["lock_digest"]


@ran
def test_no_allowed_claim_was_widened():
    c = _json(CLAIMS)
    src = _json(RESULTS / "evidence_lock" / "GEN1_CLAIM_LOCK_INPUT.json")["allowed"]
    permitted = {src["primary"], src["ranking"], src["supporting_role_A"]}
    for name, claim in c["allowed"].items():
        assert claim["text"] in permitted, f"{name} is not verbatim from the evidence lock input"


@ran
def test_every_allowed_claim_names_its_evidence_and_its_qualifiers():
    """A claim with no evidence pointer is not a claim."""
    for name, claim in _json(CLAIMS)["allowed"].items():
        assert claim["evidence"], name
        assert claim["qualifiers"], name
        assert claim["source"], name


@ran
def test_the_supporting_claim_must_travel_with_the_word_supporting():
    s = _json(CLAIMS)["allowed"]["supporting_role_A"]
    assert "SUPPORTING" in s["must_travel_with"]
    assert "does not confirm" in s["must_travel_with"]
    assert "18.3 FAILED" in s["numbers"]


@ran
def test_the_nine_forbidden_claims_are_unchanged():
    c = _json(CLAIMS)
    assert len(c["forbidden"]) == 9
    ship = SHIP_PLAN.read_text(encoding="utf-8")
    section = ship.split("## 3.5 Claims forbidden in Generation 1")[1].split("```")[1]
    in_plan = [ln.strip() for ln in section.strip().splitlines()[1:] if ln.strip()]
    assert c["forbidden"] == in_plan


@ran
def test_the_mandatory_qualifiers_say_what_was_actually_measured():
    q = _json(CLAIMS)["qualifiers"]
    assert set(q) == {"system", "vocabulary", "outcome", "evaluation", "replication"}
    assert "WM989" in q["system"] and "1,401" in q["system"]
    assert all(t in q["vocabulary"] for t in ("Acid", "Cisplatin", "CoCl2", "Dabrafenib",
                                              "Doxorubicin", "Trametinib"))
    assert "not death" in q["outcome"] and "clinical response" in q["outcome"]
    assert "clone-held-out" in q["evaluation"]
    assert q["replication"].startswith("NONE")


# ============================================================================================== #
# The worked abstract
# ============================================================================================== #
@ran
def test_the_worked_abstract_passes_the_same_instrument():
    e = _json(CLAIMS)["worked_abstract"]
    assert e["unnegated_hits"] == []
    assert all(e["qualifiers_present"].values())
    assert "demonstration" in e["status"]


@ran
def test_the_worked_abstract_reports_p_as_a_floor():
    text = _json(CLAIMS)["worked_abstract"]["abstract"]
    assert "p < 0.001" in text
    assert "0.000999" not in text, "the permutation floor is not a point estimate"
    assert "replication was not performed" in text


# ============================================================================================== #
# Verdict and handoff
# ============================================================================================== #
@ran
def test_the_verdict_is_one_of_exactly_two_values():
    v = _json(CLAIMS)
    assert v["verdict"] in ("GEN1_CLAIMS_LOCKED", "GEN1_CLAIM_LOCK_REFUSED")
    if v["verdict"] == "GEN1_CLAIMS_LOCKED":
        assert all(v["substages"].values()) and v["failing"] == []
    else:
        assert v.get("failing") or v.get("refused_at")


@ran
def test_locking_a_claim_grants_nothing():
    v = _json(CLAIMS)
    assert "grants nothing" in v["grants_nothing"]
    assert "may be lowered later; it may not be raised" in v["grants_nothing"]
    assert _json(HANDOFF)["the_ceiling_may_be_lowered_not_raised"] is True
    assert _json(HANDOFF)["no_claim_lock_outcome_reopens_an_earlier_stage"] is True


@ran
def test_the_handoff_binds_the_manuscript_to_the_lock():
    h = _json(HANDOFF)
    assert h["to_stage"] == "MANUSCRIPT + REPRODUCIBILITY PACKAGE"
    assert h["evidence_lock_digest"] == _json(CLAIMS)["evidence_lock_digest"]
    joined = " ".join(h["manuscript_must"])
    assert "--verify" in joined
    assert "mandatory qualifiers" in joined
    assert "p < 0.001" in joined
    assert "Generation 2" in joined


@ran
def test_the_claims_document_shows_the_boundary_not_just_the_ceiling():
    md = CLAIMS_MD.read_text(encoding="utf-8")
    assert "FORBIDDEN" in md and "PERMITTED" in md
    assert "NEVER" in md
    assert md.count("NEVER") >= 9
    assert "may be lowered" in md.lower()
