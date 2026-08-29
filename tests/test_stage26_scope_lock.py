"""Contracts for Stage 26 — the known-treatment-only scope lock.

Stage 25 came back positive, and the thing that goes wrong after a positive result is that the
claim grows: six observed conditions become "treatments", one melanoma line becomes "cancer", a
detection proxy becomes "response". Stage 26 is the stage that makes that growth mechanically
impossible, so these contracts exist to check that the lock is a lock and not a sentence in a
document.

The sharpest ones are the negative controls. `test_the_claim_scanner_can_actually_fire` is there
because a scan reporting zero violations is worthless if the scanner cannot detect a violation, and
`test_the_reference_leak_hazard_is_real_and_guarded` is there because `Acid` is the reference level
-- so an unknown condition reaching the dummy encoder would silently come back as the Acid score
under another name.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = RESULTS / "stage26"
SRC = ROOT / "experiments" / "run_stage26_scope_lock.py"
PLAN = (ROOT / "plans" / "(newer)practical plans"
        / "STAGE_26_KNOWN_TREATMENT_SCOPE_LOCK_V1.md")
MODEL_CARD = RESULTS / "stage24" / "tool" / "MODEL_CARD.md"

A_JSON = OUT / "stage26a_vocabulary_closure.json"
B_JSON = OUT / "stage26b_claim_surface.json"
C_JSON = OUT / "stage26c_propagation.json"
D_JSON = OUT / "stage26d_no_rescue.json"
SCOPE_MD = OUT / "GEN1_SCOPE_LIMIT.md"
VERDICT = OUT / "stage26_verdict.json"
HANDOFF = RESULTS / "stage26_handoff_to_evidence_lock.json"

ran = pytest.mark.skipif(not VERDICT.exists(), reason="Stage 26 has not been run")


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def mod():
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    import run_stage26_scope_lock as S26
    return S26


# ============================================================================================== #
# The plan and the module agree on what is being locked
# ============================================================================================== #
def test_the_plan_exists_and_names_its_parent_digest():
    text = PLAN.read_text(encoding="utf-8")
    assert "8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48" in text
    assert "KNOWN_TREATMENT_ONLY_SCOPED_LIMIT" in text
    assert "STAGE_26_SCOPE_HOLE_FOUND" in text


def test_the_module_may_not_fit_anything(mod):
    """§0.1. Stage 26 has no authority to fit, and the check that proves it must be gating."""
    src = SRC.read_text(encoding="utf-8")
    assert "MAY NOT  fit anything" in PLAN.read_text(encoding="utf-8")
    assert "Stage 26 fits nothing" in src


def test_the_adversarial_corpus_is_the_plan_s(mod):
    """The corpus was declared in the plan before it was run, so it cannot be trimmed after.

    An earlier version of this test checked four sentinel strings and passed while the plan
    declared 50 entries and the module ran 56. Spot checks do not detect drift; counts do.
    """
    import re
    plan = PLAN.read_text(encoding="utf-8")
    declared = {m.group(1).lower(): int(m.group(2))
                for m in re.finditer(r"^([A-Z]+) \((\d+)\)", plan, re.M)}
    sizes = {k: len(v) for k, v in mod.ADVERSARIAL.items()}

    assert declared, "plan §2.1 must declare a count per group"
    assert declared == sizes == mod.EXPECTED_GROUP_SIZES
    assert sum(sizes.values()) == mod.EXPECTED_ADVERSARIAL_TOTAL == 56

    flat = [s for group in mod.ADVERSARIAL.values() for s in group]
    assert len(flat) == len(set(flat)), "a duplicate would inflate the refusal count"
    # every printable ASCII entry must appear verbatim in the plan; the unicode confusables and
    # the empty/whitespace entries are described there in words, since they cannot be shown
    for s in flat:
        if s.strip() and s.isascii() and "\t" not in s and "\n" not in s:
            assert s in plan, f"{s!r} is run but not declared in the plan"


def test_the_negation_list_is_exactly_the_twelve_the_plan_declares(mod):
    """A negation list longer than the declared one is a looser gate than the one written down."""
    import re
    plan = PLAN.read_text(encoding="utf-8")
    block = plan.split("tokens must occur within 160 characters")[1].split("```")[1]
    block = block.split("\n", 1)[1]          # drop the ```text fence tag
    # quoted tokens keep their inner spaces ("no "); bare tokens are whitespace separated
    declared = {q or w for q, w in re.findall(r'"([^"]*)"|(\S+)', block)}
    assert declared == set(mod.NEGATIONS)
    assert len(mod.NEGATIONS) == 12
    assert "no " in mod.NEGATIONS, "bare 'no' matches inside not/none/know/cannot"
    assert "no" not in mod.NEGATIONS


def test_the_vocabulary_is_six_and_the_reference_is_acid(mod):
    assert mod.CONDITIONS == ("Acid", "Cisplatin", "CoCl2", "Dabrafenib", "Doxorubicin",
                              "Trametinib")
    assert mod.REFERENCE == "Acid"


def test_all_nine_forbidden_claims_are_gated(mod):
    """§3.5 of the frozen ship plan lists nine. All nine gate; none is demoted to review-only."""
    assert len(mod.FORBIDDEN_CLAIMS) == 9
    assert [k.split("_", 1)[0] for k in mod.FORBIDDEN_CLAIMS] == [str(i) for i in range(1, 10)]
    assert all(v for v in mod.FORBIDDEN_CLAIMS.values()), "every claim needs patterns"


def test_the_scanner_excludes_itself_but_not_the_shipped_surfaces(mod):
    surfaces = {Path(p).name for p in map(str, mod._surfaces())}
    assert {"gen1_predictor.py", "gen1_cli.py", "MODEL_CARD.md", "io_schema.json",
            "GEN1_SCOPE_LIMIT.md"} <= surfaces
    assert "run_stage26_scope_lock.py" not in surfaces, "a scanner that scans itself false-fails"


# ============================================================================================== #
# The negative controls — the instrument has to be able to fail
# ============================================================================================== #
def test_the_claim_scanner_can_actually_fire(mod):
    """Zero violations from an instrument that cannot detect one is the most comfortable lie."""
    c = mod._canary()
    assert c["claims_probed"] == 9
    assert c["missed"] == [], f"the scanner cannot detect {c['missed']}"
    assert c["negated_twin_false_alarms"] == []


def test_a_planted_claim_is_caught(mod):
    """Directly: an unnegated forbidden sentence produces an unnegated hit."""
    hits = mod._scan_text("This model generalizes to any treatment.",
                          mod.FORBIDDEN_CLAIMS["1_unseen_treatment_generalization"])
    assert any(not h["negated"] for h in hits)
    safe = mod._scan_text("This model never generalizes to any treatment outside the six.",
                          mod.FORBIDDEN_CLAIMS["1_unseen_treatment_generalization"])
    assert all(h["negated"] for h in safe)


@ran
def test_the_reference_leak_hazard_is_real_and_guarded():
    """Acid is the reference level, so an unknown string reaching the encoder returns its score."""
    r = _json(A_JSON)["reference_leak_test"]
    assert r["raw_encoder_dummy_sum"] == 0.0
    assert r["unknown_via_raw_encoder"] == r["acid_score"], "the hazard must be demonstrated"
    assert r["hazard_is_real"] is True
    assert r["guard_holds"] is True, "predict() must block what the encoder would let through"


# ============================================================================================== #
# 26A — vocabulary closure
# ============================================================================================== #
@ran
def test_every_adversarial_string_was_refused():
    a = _json(A_JSON)
    assert a["n_adversarial_strings"] >= 50
    assert a["n_refused"] == a["n_adversarial_strings"]
    assert a["leaked"] == []
    for group in a["refusals_by_group"].values():
        for row in group:
            assert row["support_status"] == "UNSUPPORTED_TREATMENT"
            assert row["score"] is None


@ran
def test_the_drugs_a_real_user_would_try_are_refused():
    rows = {r["input"]: r for r in _json(A_JSON)["refusals_by_group"]["pharmacology"]}
    # Vemurafenib is THE drug for this BRAF-V600E line; Carboplatin is one substitution from a
    # condition that IS supported. These are the two a helpful-seeming tool would get wrong.
    for drug in ("Vemurafenib", "Carboplatin", "Cobimetinib", "Pembrolizumab"):
        assert rows[drug]["refused"] is True
        assert rows[drug]["score"] is None


@ran
def test_the_vocabulary_is_closed_by_geometry_not_only_by_a_list():
    s = _json(A_JSON)["structural_closure"]
    assert s["design_columns"] == s["expected"] == 309
    assert s["K"] + s["nuisance"] + s["dummies"] + s["interaction"] == 309
    assert s["dummies"] == 5 and len(s["vocabulary"]) == 6


@ran
def test_batching_does_not_move_a_score_or_flip_an_ordering():
    """§7.1 R2 bounds a cell at 1e-12; R3 requires within-clone ordering to be untouched."""
    m = _json(A_JSON)["mixed_request"]
    assert m["routing_correct"] is True
    assert m["max_batch_vs_solo_difference"] <= 1e-12
    assert m["within_clone_ordering_identical_across_batch_sizes"] is True
    assert "BLAS" in m["named_cause_23_5_7_1_R4"], "§7.1 R4: name the cause, not 'floating point'"


# ============================================================================================== #
# 26B / 26C / 26D
# ============================================================================================== #
@ran
def test_no_forbidden_claim_appears_unnegated_on_a_shipped_surface():
    b = _json(B_JSON)
    assert b["violations"] == []
    assert all(v.get("present") for v in b["per_file"].values())
    assert b["canary"]["detects_every_forbidden_claim"] is True


@ran
def test_the_limit_reaches_the_caller_on_the_failing_paths_too():
    c = _json(C_JSON)
    assert set(c["paths"]) >= {"supported", "unsupported_treatment", "missing_nuisance",
                               "bad_feature_schema"}
    for name, p in c["paths"].items():
        assert p["all_carry_known_limitations"], name
        assert p["status_matches"], name
        assert p["no_calibrated_probability_key"], name


@ran
def test_the_stage_25_verdict_unlocks_a_claim_not_a_computation():
    r = _json(C_JSON)["ranking"]
    assert r["with_verdict_status"] == "SUPPORTED"
    assert r["without_verdict_status"] == "NOT_SUPPORTED"
    assert r["with_verdict_exposes_order"] is True
    assert r["without_verdict_withholds_order"] is True
    assert r["scores_identical_either_way"] is True


@ran
def test_the_cli_never_reports_success_when_it_refused():
    cli = _json(C_JSON)["cli"]
    assert cli["all_known"]["exit_code"] == 0
    assert cli["one_unknown"]["exit_code"] == 2
    assert cli["no_nuisance"]["exit_code"] == 2
    assert cli["unreadable"]["exit_code"] == 3
    assert cli["unreadable"]["stderr_has_status"] is True
    assert "UNSUPPORTED_TREATMENT" in cli["one_unknown"]["statuses"]
    # the CLI is the surface a user actually touches; a refusal without limitations tells them
    # less than a score does
    assert all(cli["every_printed_row_carries_known_limitations"].values())


@ran
def test_stage_26_fitted_nothing_and_moved_no_frozen_artifact():
    d = _json(D_JSON)
    assert d["fitting_tokens_found"] == []
    assert d["artifact_npz_unchanged"] is True
    assert d["ship_plan_digest_holds"] is True
    for k, v in d["hash_recheck"].items():
        if k == "MODEL_CARD.md":
            continue          # the one file §6.2 permits appending to
        assert v["unchanged"] is True, k


# ============================================================================================== #
# 26E — the record
# ============================================================================================== #
@ran
def test_the_verdict_is_one_of_exactly_two_values():
    v = _json(VERDICT)
    assert v["verdict"] in ("KNOWN_TREATMENT_ONLY_SCOPED_LIMIT", "STAGE_26_SCOPE_HOLE_FOUND")
    if v["verdict"] == "KNOWN_TREATMENT_ONLY_SCOPED_LIMIT":
        assert all(v["substages"].values()) and v["failing"] == []
    else:
        assert v["failing"], "a scope hole must name what failed"


@ran
def test_stage_26_grants_no_claim():
    v = _json(VERDICT)
    assert "grants no new claim" in v["grants_no_claim"]
    assert "reopens no earlier stage" in v["grants_no_claim"]
    assert v["parent_plan_digest_holds"] is True


@ran
def test_the_model_card_change_is_append_only_against_the_24f_hash():
    card = _json(VERDICT)["model_card_update"]
    assert card["append_only_proof"] is True
    assert card["base_is_byte_identical_to_the_24F_card"] is True
    assert card["base_sha256"] == card["frozen_24F_sha256"]
    assert card["rerun_is_byte_idempotent"] is True
    assert card["bytes_after"] > card["bytes_frozen_24F"]
    text = MODEL_CARD.read_text(encoding="utf-8")
    assert text.count(card["delimiter"]) == 1, "a re-run must not stack a second section"


@ran
def test_the_scope_document_states_the_boundary_it_is_supposed_to_state():
    md = SCOPE_MD.read_text(encoding="utf-8")
    assert "KNOWN_TREATMENT_ONLY_SCOPED_LIMIT" in md
    for c in ("Acid", "Cisplatin", "CoCl2", "Dabrafenib", "Doxorubicin", "Trametinib"):
        assert c in md
    assert "UNSUPPORTED_TREATMENT" in md
    assert "Vemurafenib" in md
    # the outcome must not be relabelled as death or clinical response
    assert "not proven death" in md
    assert "cannot produce a valid `B`" in md
    # every one of the nine prohibitions is written as its own negation
    assert md.count("NEVER") >= 9


@ran
def test_the_frozen_hashes_are_verified_after_the_run_not_only_before():
    """§5. Checking only at the start proves nothing about what the run then did."""
    post = _json(VERDICT)["frozen_hashes_after_the_run"]
    assert set(post) >= {"io_schema.json", "example_clones.csv", "stage24_oof_for_stage25.csv",
                         "stage24_w5_artifact.json", "stage24_w5_artifact.npz",
                         "stage25_verdict.json"}
    assert all(post.values()), [k for k, v in post.items() if not v]


@ran
def test_every_substage_came_from_the_same_module():
    """§5.1. A verdict merged from substages of different module versions means nothing."""
    s = _json(VERDICT)["substage_module_stamps"]
    assert s["all_equal"] is True
    running = s["running_module"]
    for k in ("26A", "26B", "26C", "26D"):
        assert s[k] == running, f"{k} was produced by a different version of the executor"


@ran
def test_every_evidence_lock_input_is_a_real_path():
    """The evidence lock hashes these. A prose string like 'artifact.npz + .json' is a bug.

    One path is legitimately absent from a fresh clone: the 44 MB model artifact is gitignored and
    rebuilt by `--stage 24c` before anything else runs. This test used to assert every path exists,
    which passed on a machine that had already built it and failed in CI on every commit. The
    exemption is DERIVED from the evidence manifest's own `git_ignored` list rather than hardcoded,
    and the manifest separately gates that list down to exactly that one file -- so a second
    unbuildable path still fails here.
    """
    v = _json(VERDICT)
    assert v["evidence_missing"] == []
    assert v["evidence_paths_verified_present"] is True

    manifest = _json(RESULTS / "evidence_lock" / "GEN1_EVIDENCE_MANIFEST.json")
    rebuildable = set(manifest["git_ignored"])
    assert rebuildable == {"results/stage24/stage24_w5_artifact.npz"}

    absent = []
    for group, paths in _json(HANDOFF)["evidence_to_lock"].items():
        assert isinstance(paths, list), f"{group} must be a list of paths, not a sentence"
        absent += [p for p in paths if not (ROOT / p).exists()]
    assert set(absent) <= rebuildable, \
        f"absent and not rebuildable: {sorted(set(absent) - rebuildable)}"


@ran
def test_the_handoff_names_what_the_evidence_lock_must_freeze():
    h = _json(HANDOFF)
    assert h["to_stage"] == "GEN-1 EVIDENCE LOCK"
    assert h["closed_vocabulary"] == ["Acid", "Cisplatin", "CoCl2", "Dabrafenib", "Doxorubicin",
                                      "Trametinib"]
    assert h["reference_condition"] == "Acid"
    assert len(h["forbidden_claims"]) == 9
    assert set(h["evidence_to_lock"]) == {"benchmark", "out_of_fold_predictions", "tool",
                                          "ranking_verdict", "limitations"}
    assert h["no_stage_26_outcome_reopens_an_earlier_stage"] is True
