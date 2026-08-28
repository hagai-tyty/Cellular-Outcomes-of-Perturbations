"""Contracts for the Generation-1 manuscript and reproducibility package.

This is the weakest position an instrument in this project has been in: the prose being checked was
written by the same process doing the checking. The only defence is a checker that refuses
mechanically and has been shown to do so, which is why `test_the_checker_refuses_a_broken_copy` is
the load-bearing contract here and everything else is downstream of it.

The second-sharpest is `test_no_number_in_the_manuscript_is_unsourced`. A manuscript is a place
where a figure gets retyped, and a retyped figure is how a paper ends up disagreeing with its own
data.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = RESULTS / "manuscript"
SRC = ROOT / "experiments" / "run_gen1_manuscript.py"
PLAN = ROOT / "plans" / "(newer)practical plans" / "GEN1_MANUSCRIPT_PACKAGE_V1.md"

MANUSCRIPT = OUT / "MANUSCRIPT.md"
REPRO = OUT / "REPRODUCIBILITY.md"
COMPLIANCE = OUT / "manuscript_compliance.json"
CONTROLS = OUT / "manuscript_controls.json"
VERDICT = OUT / "GEN1_MANUSCRIPT.json"
DIGEST = OUT / "GEN1_PACKAGE_DIGEST.json"
HANDOFF = RESULTS / "gen1_handoff_to_manuscript.json"

ran = pytest.mark.skipif(not VERDICT.exists(), reason="the manuscript stage has not been run")


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def mod():
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    import run_gen1_manuscript as MS
    return MS


# ============================================================================================== #
# The checker has to be able to refuse
# ============================================================================================== #
@ran
def test_the_checker_refuses_a_broken_copy():
    """Four ways to break the manuscript; all four must be caught, and the real one must pass."""
    c = _json(CONTROLS)["controls"]
    assert c["a planted forbidden claim is caught"] is True
    assert c["a dropped qualifier is caught"] is True
    assert c["p quoted as a point estimate is caught"] is True
    assert c["a changed number is caught"] is True
    assert c["the real manuscript passes all four"] is True


def test_a_planted_forbidden_claim_is_caught_end_to_end(mod):
    """Directly, not by reading a JSON that says so."""
    import run_gen1_claim_lock as CL
    full = CL.combined_patterns()
    clean = MANUSCRIPT.read_text(encoding="utf-8")
    assert CL.scan(clean, full) == []
    for planted in ("The model generalises to new treatments.",
                    "This supports clinical decision-making.",
                    "Independently replicated in an external cohort."):
        assert CL.scan(clean + "\n\n" + planted + "\n", full), planted


def test_the_controls_never_modify_the_manuscript(mod):
    src = SRC.read_text(encoding="utf-8")
    assert "in-memory copies" in src
    assert "never modified" in src
    # the file on disk must be untouched by a control run
    before = MANUSCRIPT.read_bytes()
    mod.negative_controls()
    assert MANUSCRIPT.read_bytes() == before


# ============================================================================================== #
# Numbers
# ============================================================================================== #
@ran
def test_no_number_in_the_manuscript_is_unsourced():
    c = _json(COMPLIANCE)
    assert c["untraceable_numbers"] == []
    assert c["numbers_checked"] >= 17


@ran
def test_the_headline_numbers_are_the_locked_ones():
    text = MANUSCRIPT.read_text(encoding="utf-8")
    h = _json(RESULTS / "evidence_lock" / "GEN1_EVIDENCE_LOCK.json")["headline_numbers"]
    assert f"{h['delta_RANK']:+.6f}" in text
    assert f"{h['bootstrap_ci95'][0]:+.6f}" in text and f"{h['bootstrap_ci95'][1]:+.6f}" in text
    assert f"{h['R_W5']:.6f}" in text and f"{h['R_W4']:.6f}" in text
    assert str(h["eligible_clones"]) in text


@ran
def test_p_is_a_floor_and_never_a_point_estimate():
    text = MANUSCRIPT.read_text(encoding="utf-8")
    assert "p < 0.001" in text
    assert "0.000999" not in text
    assert "floor of a 1,000-draw permutation test" in text


def test_a_line_wrap_cannot_hide_a_number(mod):
    """A wrap between "472 were" and "never detected" defeated this check on its first run."""
    src = SRC.read_text(encoding="utf-8")
    assert 'flat = " ".join(text.split())' in src
    assert "line wrap" in src


# ============================================================================================== #
# Claims and qualifiers
# ============================================================================================== #
@ran
def test_no_forbidden_claim_appears_unnegated():
    c = _json(COMPLIANCE)
    assert c["forbidden_hits"] == []
    assert c["package_forbidden_hits"] == []
    assert "clause-scoped" in c["instrument"]


@ran
def test_every_mandatory_qualifier_is_present():
    c = _json(COMPLIANCE)
    assert c["missing_qualifiers"] == []
    assert c["abstract_missing_qualifiers"] == [], \
        "the abstract travels alone and must carry system, vocabulary and outcome by itself"


@ran
def test_the_manuscript_separates_limitations_from_forbidden_claims():
    """What the result cannot support, and what may not be said, are different things."""
    text = MANUSCRIPT.read_text(encoding="utf-8")
    assert "## Limitations" in text
    assert "## What this does not show" in text
    assert text.index("## Limitations") < text.index("## What this does not show")
    forbidden = text.split("## What this does not show")[1]
    assert forbidden.count("NEVER") >= 9


@ran
def test_the_uncomfortable_limitations_survive_into_the_manuscript():
    """A manuscript that drops the inconvenient ones is the failure mode here."""
    text = MANUSCRIPT.read_text(encoding="utf-8")
    assert "3.45x" in text, "abundance still dominates state"
    assert "Cisplatin is negligible" in text and "Doxorubicin is negative" in text
    assert "18.3 FAILED" in text and "0.45" in text, "Role A's failed gate and audited power"
    assert "not death" in text
    assert "No independent biological replication" in text


@ran
def test_every_required_section_is_present():
    c = _json(COMPLIANCE)
    assert c["missing_sections"] == []


# ============================================================================================== #
# The package
# ============================================================================================== #
@ran
def test_every_documented_command_exists_and_accepts_its_flag():
    """A package that documents a flag the code does not accept is worse than no package."""
    v = _json(VERDICT)
    assert v["substages"]["MS-E_package"] is True
    text = REPRO.read_text(encoding="utf-8")
    for script, stage in re.findall(r"python (experiments/[\w./]+\.py)(?: --stage (\S+))?", text):
        assert (ROOT / script).exists(), script
        if stage:
            src = (ROOT / script).read_text(encoding="utf-8")
            m = re.search(r"choices=\[([^\]]*)\]", src)
            if m:
                assert stage in re.findall(r'"([^"]+)"', m.group(1)), f"{script} --stage {stage}"


@ran
def test_the_package_names_what_it_does_not_contain():
    text = REPRO.read_text(encoding="utf-8")
    assert "44 MB" in text and "does NOT contain it" in text
    assert "--stage 24c" in text
    assert "GSE279162" in text and "GSE227151" in text
    assert "Naming a gap is not closing it" in text


@ran
def test_the_package_states_the_long_runtime_honestly():
    text = REPRO.read_text(encoding="utf-8")
    assert "10.7 h" in text
    assert "not an estimate" in text


@ran
def test_verification_comes_before_everything_else():
    text = REPRO.read_text(encoding="utf-8")
    assert text.index("--verify") < text.index("## 2. Environment")
    assert "run_gen1_evidence_lock.py --verify" in text
    assert "run_gen1_claim_lock.py --verify" in text


# ============================================================================================== #
# Verdict, digests, handoff
# ============================================================================================== #
@ran
def test_both_locks_verified_before_a_sentence_was_checked():
    a = _json(VERDICT)
    h = _json(HANDOFF)
    assert a["evidence_digest"] == h["evidence_lock_digest"]
    assert a["claim_digest"] == h["claim_digest"]
    assert a["substages"]["MS-A_locks_verify"] is True


@ran
def test_both_digests_are_quoted_in_the_manuscript():
    text = MANUSCRIPT.read_text(encoding="utf-8")
    a = _json(VERDICT)
    assert a["evidence_digest"] in text
    assert a["claim_digest"] in text


@ran
def test_the_package_has_its_own_digest(mod):
    """The manuscript layer pins itself, as the claim layer pins itself onto the evidence layer."""
    import hashlib
    d = _json(DIGEST)
    assert len(d["package_digest"]) == 64
    assert "MISSING" not in d["covers"].values()
    assert set(d["covers"]) == set(mod.PACKAGE_FILES)
    canonical = "\n".join(f"{k}  {d['covers'][k]}" for k in sorted(d["covers"]))
    assert hashlib.sha256(canonical.encode("utf-8")).hexdigest() == d["package_digest"]


@ran
def test_the_verdict_is_one_of_exactly_two_values():
    v = _json(VERDICT)
    assert v["verdict"] in ("GEN1_MANUSCRIPT_READY", "GEN1_MANUSCRIPT_REFUSED")
    if v["verdict"] == "GEN1_MANUSCRIPT_READY":
        assert all(v["substages"].values()) and v["failing"] == []
    else:
        assert v.get("failing") or v.get("refused_at")


@ran
def test_ready_does_not_claim_the_science_is_good():
    v = _json(VERDICT)
    note = v["what_ready_does_not_mean"]
    assert "NOT a judgement that the science is good" in note
    assert "a reviewer will agree" in note
