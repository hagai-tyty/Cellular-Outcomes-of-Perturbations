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

# Amendment V1.2 added one negative control per structure check.
STRUCTURE_CONTROLS = (
    "sections out of order are caught",
    "a missing Declarations heading is caught",
    "an unstructured abstract is caught",
    "an over-long abstract is caught",
    "too many keywords are caught",
)


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
    """Every way to break the manuscript must be caught, and the real one must pass."""
    c = _json(CONTROLS)["controls"]
    assert c["a planted forbidden claim is caught"] is True
    assert c["a dropped qualifier is caught"] is True
    assert c["p quoted as a point estimate is caught"] is True
    assert c["a changed number is caught"] is True
    for control in STRUCTURE_CONTROLS:
        assert c[control] is True, control
    assert c["overstated freeze wording is caught"] is True
    assert c["an archive check that could run the installed copy is caught"] is True
    assert c["the real package document passes the archive check"] is True
    assert c["the real manuscript passes every control"] is True


def test_a_planted_forbidden_claim_is_caught_end_to_end(mod):
    """Directly, not by reading a JSON that says so.

    Scans through the same reference partition the gate uses. Scanning the raw document instead
    made this contract test a different document than the gate does, and the cited paper titles —
    third-party text containing the word `cancer` — failed it. Third time this project has learned
    that a check must read exactly what the thing it is checking reads.
    """
    import run_gen1_claim_lock as CL
    full = CL.combined_patterns()

    def scan(doc: str) -> list:
        prose, _titles, _problems = mod._partition_references(doc)
        return CL.scan(prose, full)

    clean = MANUSCRIPT.read_text(encoding="utf-8")
    assert scan(clean) == []
    for planted in ("The model generalises to new treatments.",
                    "This supports clinical decision-making.",
                    "Independently replicated in an external cohort."):
        assert scan(clean + "\n\n" + planted + "\n"), planted


def test_the_controls_never_modify_the_manuscript(mod):
    """Running the controls must leave BOTH the manuscript and this stage's own outputs alone.

    `write=False` matters: an earlier version called the writer, so `pytest` rewrote
    manuscript_controls.json with a fresh runtime and left the working tree dirty. A test that
    modifies a committed results artifact is a side effect, not a check.
    """
    src = SRC.read_text(encoding="utf-8")
    assert "in-memory copies" in src
    assert "never modified" in src

    before = MANUSCRIPT.read_bytes()
    controls_before = CONTROLS.read_bytes()
    r = mod.negative_controls(write=False)
    assert r["all_passed"] is True
    assert MANUSCRIPT.read_bytes() == before
    assert CONTROLS.read_bytes() == controls_before, "the suite must not dirty the working tree"


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
    # Amendment V1.2 made both subsections of Discussion. They stay separate, and in this order.
    discussion = text.index("\n## Discussion\n")
    limitations = text.index("\n### Limitations\n")
    forbidden_at = text.index("\n### What this does not show\n")
    generation_2 = text.index("\n### Generation 2\n")
    conclusions = text.index("\n## Conclusions\n")
    assert discussion < limitations < forbidden_at < generation_2 < conclusions
    forbidden = text[forbidden_at:generation_2]
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


# ================================================================================================ #
# The verifiers must check the verdict, not just the bytes
# ================================================================================================ #
@pytest.mark.parametrize("module_name,record_attr,refused", [
    ("run_gen1_evidence_lock", "LOCK_JSON", "GEN1_EVIDENCE_REFUSED"),
    ("run_gen1_claim_lock", "CLAIMS_JSON", "GEN1_CLAIM_LOCK_REFUSED"),
    ("run_gen1_manuscript", "VERDICT_JSON", "GEN1_MANUSCRIPT_REFUSED"),
])
def test_verify_refuses_when_the_stage_itself_refused(tmp_path, module_name, record_attr, refused):
    """Byte-intact is not the same as passed, and the two came apart in practice.

    The manuscript stage recorded GEN1_MANUSCRIPT_REFUSED; its covered files hashed exactly as
    recorded, because a refused run still writes them; and --verify returned clean. CI runs only
    the three --verify commands, so it stayed green over a refused package. Each verifier must
    therefore read the recorded verdict, and this proves it does by planting a refusal.
    """
    import json
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    mod = __import__(module_name)

    real = getattr(mod, record_attr)
    if not real.is_file():
        pytest.skip(f"{module_name} has not been run")

    planted = tmp_path / real.name
    record = json.loads(real.read_text(encoding="utf-8"))
    record["verdict"] = refused
    planted.write_text(json.dumps(record, indent=2), encoding="utf-8")

    try:
        setattr(mod, record_attr, planted)
        r = mod.run_verify()
    finally:
        setattr(mod, record_attr, real)

    assert r["clean"] is False, (
        f"{module_name} --verify called a refused stage clean; CI would stay green over it")
    assert r["verdict"].endswith("STAGE_REFUSED"), (
        f"{module_name} reported {r['verdict']}, which does not say the stage refused")


def test_a_passing_stage_still_verifies_clean():
    """The guard above must not be satisfiable by a verifier that always refuses."""
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    for name in ("run_gen1_evidence_lock", "run_gen1_claim_lock", "run_gen1_manuscript"):
        mod = __import__(name)
        r = mod.run_verify()
        assert r["clean"] is True, f"{name} --verify is not clean on the committed tree: {r}"


# ================================================================================================ #
# The cascade tool, and the documents it is responsible for
# ================================================================================================ #
def test_every_quoted_digest_is_current():
    """A document that quotes a digest must quote the live one.

    `GEN1_MANUSCRIPT_PACKAGE_V1.md` pinned both digests in its Entry block and nothing checked it,
    so it sat stale through several cascades while every lock verified clean -- the package digest
    hashes that file, but hashing it only proves it has not changed, not that what it says is true.
    """
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    import cascade_gen1 as C

    evidence, claim, _ = C.digests()
    stale = []
    for p in C.EVIDENCE_ONLY:
        if evidence not in p.read_text(encoding="utf-8"):
            stale.append(f"{p.name}: evidence digest")
    for p in C.BOTH:
        text = p.read_text(encoding="utf-8")
        if evidence not in text:
            stale.append(f"{p.name}: evidence digest")
        if claim not in text:
            stale.append(f"{p.name}: claim digest")
    assert not stale, ("documents quoting a superseded digest: " + "; ".join(stale) +
                       " -- run python experiments/cascade_gen1.py")


def test_no_document_quotes_a_digest_without_being_registered():
    """A new document that starts quoting digests must be added to the cascade tool.

    Otherwise it is maintained by hand, which is how every stale pin in this project happened. The
    RECORDs are excluded: they quote superseded digests on purpose, as history. So are the
    files a stage generates, which quote a digest because they are derived from one.
    """
    import subprocess
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    import cascade_gen1 as C

    registered = {p.resolve() for p in C.EVIDENCE_ONLY + C.BOTH} | {
        p.resolve() for p in C.GENERATED}
    tracked = subprocess.run(["git", "ls-files", "*.md"], cwd=ROOT,
                             capture_output=True, text=True).stdout.split()
    unregistered = []
    for rel in tracked:
        if "RECORDs/" in rel:
            continue
        p = ROOT / rel
        if not p.is_file() or p.resolve() in registered:
            continue
        text = p.read_text(encoding="utf-8", errors="replace")
        if C.EVIDENCE_PIN.search(text) or C.CLAIM_PIN.search(text):
            unregistered.append(rel)
    assert not unregistered, (
        "these quote a lock digest but are not re-pinned by experiments/cascade_gen1.py: "
        + "; ".join(unregistered))


# ================================================================================================ #
# Amendment V1.2 -- the BMC Research-article structure
# ================================================================================================ #
def test_the_manuscript_has_the_research_article_structure(mod):
    """The shape Amendment V1.2 requires, checked on the real document rather than through a JSON."""
    text = MANUSCRIPT.read_text(encoding="utf-8")
    problems = mod.structure_problems(text)
    assert not any(problems.values()), problems
    assert re.findall(r"^## (.+?)[ \t]*$", text, re.M) == mod.TOP_LEVEL_ORDER


def test_every_structure_check_refuses_a_broken_copy(mod):
    """Directly, not by reading a JSON that says so."""
    controls = mod.negative_controls(write=False)["controls"]
    for control in STRUCTURE_CONTROLS + ("the real manuscript passes every control",):
        assert controls[control] is True, control


def test_every_v1_section_survives_the_restructure():
    """V1.2 renamed and moved sections. It may not have dropped one."""
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for heading in ("### Relation to prior work", "### Data", "### The tool", "### Limitations",
                    "### What this does not show", "### Generation 2",
                    "### Availability of data and materials"):
        assert f"\n{heading}\n" in text, heading


# ================================================================================================ #
# Amendment V1.3 -- when the protocol was frozen, said precisely
# ================================================================================================ #
def test_no_release_document_overstates_when_the_protocol_was_frozen(mod):
    """The ranking test was frozen before any ranking statistic was computed, after earlier predictive
    analyses of the same data -- not before any result existed. Checked on every release document."""
    for p in (MANUSCRIPT, REPRO, ROOT / "README.md", OUT / "SUBMISSION.md", ROOT / ".zenodo.json",
              ROOT / "CITATION.cff"):
        assert mod._overstated_freeze(p.read_text(encoding="utf-8")) == [], p.name


def test_the_precise_freeze_wording_is_present():
    flat = " ".join(MANUSCRIPT.read_text(encoding="utf-8").split())
    assert "before any ranking statistic was computed" in flat
    assert "earlier predictive analyses of the same data" in flat
    assert "fixed before any model was fitted" in flat


def test_overstated_wording_is_caught_across_a_wrap_or_a_blockquote(mod):
    assert mod._overstated_freeze("frozen before any\nresult existed")
    assert mod._overstated_freeze("> frozen before any result\n> existed")
    assert mod._overstated_freeze("fixed before the numbers existed")
    assert not mod._overstated_freeze("fixed before any ranking statistic was computed")
    assert not mod._overstated_freeze("fixed before any of these numbers existed")


# ================================================================================================ #
# Amendment V1.3 -- a checkout is not the archive
# ================================================================================================ #
def test_the_package_distinguishes_a_checkout_from_the_archive(mod):
    text = REPRO.read_text(encoding="utf-8")
    assert mod.archive_problems(text) == []
    assert "44 MB" in text and "does NOT contain it" in text


def test_the_archive_is_verified_with_its_own_code_not_an_installed_copy(mod):
    """Inside an unpacked archive, a plain `python -m cellfate.gen1_cli` imported the checkout's
    editable install and exited 0 -- a pass on the wrong code."""
    text = REPRO.read_text(encoding="utf-8")
    section = text.split(mod.ARCHIVE_SECTION, 1)[1].split("\n## ", 1)[0]
    for command in mod.ARCHIVE_COMMANDS:
        assert command in section, command
    broken = text.replace("PYTHONPATH=src python -m cellfate.gen1_cli", "python -m cellfate.gen1_cli")
    assert mod.archive_problems(broken)


# ---- the reference block, as a list (Amendment V1.6) --------------------------------------------
def test_the_reference_list_exempts_one_title_per_entry(mod):
    """The manuscript's own block is a numbered list now, so that a journal reads what the scanner
    reads. The exemption must stay line-precise: one title per entry, nothing else."""
    text = MANUSCRIPT.read_text(encoding="utf-8")
    block = text.split("## References", 1)[1].split("\n---", 1)[0]
    assert "```" not in block, "the block is a list, not a fenced listing"
    openers = [ln for ln in block.splitlines() if re.match(r"^\d+\.\s", ln)]
    prose, titles, problems = mod._partition_references(text)
    assert problems == []
    assert len(titles) == len(openers) >= 6
    for t in titles:
        assert t in block and not re.match(r"^\d+\.\s", t)


def test_a_sentence_planted_inside_the_reference_block_is_still_scanned(mod):
    """The loophole that started this: exempting too much absorbs a planted claim."""
    import run_gen1_claim_lock as CL
    full = CL.combined_patterns()
    text = MANUSCRIPT.read_text(encoding="utf-8")
    prose, _titles, _problems = mod._partition_references(text)
    assert CL.scan(prose, full) == []
    planted = "    The model generalises to new treatments.\n"
    for anchor in ("\n2. ", "\n5. "):
        doctored = text.replace(anchor, "\n" + planted + anchor[1:], 1)
        p, _t, _pr = mod._partition_references(doctored)
        assert CL.scan(p, full), f"a claim planted before {anchor.strip()} escaped the scan"


def test_the_reference_parser_still_reads_the_fenced_form(mod):
    """Documents written before V1.6 keep their fenced block; both forms parse the same way."""
    fenced = ("## References\n\n```text\n[1] Doe J.\n    A title about cancer drug resistance.\n"
              "    Journal. 2020;1:1. doi:10.0000/x\n```\n\n---\n")
    prose, titles, problems = mod._partition_references(fenced)
    assert titles == ["A title about cancer drug resistance."] and problems == []
    assert "A title about cancer" not in prose


def test_an_entry_without_an_identifier_is_reported(mod):
    listed = ("## References\n\n1. Doe J.\n    A title.\n    Journal. 2020;1:1.\n\n---\n")
    _prose, _titles, problems = mod._partition_references(listed)
    assert any("without a resolvable identifier" in p for p in problems)
