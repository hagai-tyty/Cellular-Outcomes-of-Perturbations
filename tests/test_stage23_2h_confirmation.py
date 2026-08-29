"""Contracts for confirmation protocol V5 and Stage 23.2H.

Two jobs.

First, pin V5 itself: what it retains from V4, what it deliberately supersedes, and — most
importantly — that every anti-gaming clause survived the supersession. V5 loosens exactly one gate
(the >= 140 floor becomes a measured power gate) and it would be very easy for that to quietly drag
the rest of the firewall with it. These contracts exist so that it cannot.

Second, pin the execution: that the reconstructions reproduce the authors' own serialized objects,
that the frozen cohort is what V5 says it is, that the power gate was recorded BEFORE the
confirmatory statistics, and that the verdict follows mechanically from the gates rather than from
anything anyone hoped for.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PLANS = ROOT / "plans" / "(newer)practical plans"
OUT = ROOT / "results" / "stage23_2h"

V5 = PLANS / "STAGE_23_2_ROLE_A_CONFIRMATION_V5.md"
V4 = PLANS / "arcive" / "STAGE_23_2_ROLE_A_CONFIRMATION_V4.md"
ADDENDUM1 = PLANS / "STAGE_23_2_V5_ADDENDUM_1_POWER_CURVE.md"

BENCH = OUT / "stage23_2h_benchmark.json"
BENCH_BUG = OUT / "stage23_2h_benchmark_authorbug.json"
REPR = OUT / "stage23_2h_representation.json"
POWER = OUT / "stage23_2h_power.json"
CONFIRM = OUT / "stage23_2h_confirmation.json"
VERDICT = OUT / "stage23_2h_verdict.json"

has_v5 = pytest.mark.skipif(not V5.exists(), reason="V5 not written")
ran_a = pytest.mark.skipif(not BENCH.exists(), reason="23.2H-A has not been run")
ran_bug = pytest.mark.skipif(not BENCH_BUG.exists(), reason="the sensitivity arm has not been run")
ran_b = pytest.mark.skipif(not REPR.exists(), reason="23.2H-B has not been run")
ran_c = pytest.mark.skipif(not POWER.exists(), reason="23.2H-C has not been run")
ran_d = pytest.mark.skipif(not CONFIRM.exists(), reason="23.2H-D has not been run")
ran_e = pytest.mark.skipif(not VERDICT.exists(), reason="23.2H-E has not been run")


def _v5() -> str:
    return V5.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    """Collapse whitespace so a contract pins wording, not the line-wrap it happens to have."""
    return " ".join(text.split())


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


# ============================================================================================== #
# V5 as a document
# ============================================================================================== #
@has_v5
def test_v4_is_preserved_unedited_in_the_archive():
    """V5 supersedes V4; it does not rewrite it. The old verdict must stay readable."""
    assert V4.exists(), "V4 must remain in arcive/"
    assert not (PLANS / "STAGE_23_2_ROLE_A_CONFIRMATION_V4.md").exists(), "V4 is still live"
    for rel in (f"plans/(newer)practical plans/{V4.name}",
                f"plans/(newer)practical plans/arcive/{V4.name}"):
        got = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=ROOT,
                             capture_output=True).stdout
        if got:
            canon = got.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
            assert canon == V4.read_bytes().replace(b"\r\n", b"\n"), "V4 was modified"
            break


@has_v5
def test_v5_names_exactly_what_it_changes():
    v5 = _v5()
    assert "# 13. Record of what V5 changes" in v5
    changes = v5.split("# 13. Record of what V5 changes")[1]
    for item in ("§6", "§7", "§9", "§4", "§8"):
        assert item in changes, item
    assert "Everything else is V4" in changes


@has_v5
def test_v5_carries_forward_the_v4_clauses_it_does_not_touch():
    v5 = _v5()
    carried = v5.split("# 2. Carried forward from V4, unchanged")[1].split("# 3.")[0]
    for clause in ("§11", "§12", "§13", "§14", "§15.1", "§15.2", "§15.4", "§16",
                   "§18.1", "§18.4", "§18.5", "§18.6"):
        assert clause in carried, f"V5 dropped {clause}"
    assert "depth_complete_nuisance_control" in v5, "the corrected hypothesis must be restated"


@has_v5
def test_the_false_premise_is_named_rather_than_quietly_dropped():
    """23.2G's verdict was wrong. V5 has to say so, and say why, not just move on."""
    v5 = _v5()
    assert "QUALIFYING_SET_EMPTY_FROM_FROZEN_SEARCH_SPACE" in v5
    assert "stands unedited" in v5
    assert "never deposited in GEO" in v5
    assert "stage_23_2G_step1_REOPENED_NEW_EVIDENCE.md" in v5


# ---- the one loosened gate, and its guard rails ------------------------------------------------#
@has_v5
def test_the_140_floor_is_replaced_by_a_measured_gate_not_removed():
    v5 = _v5()
    assert "power >= 0.80" in v5
    assert "UNDER ITS OWN REALIZED GEOMETRY" in v5
    assert "The measurement instrument is unchanged" in v5
    # it must be able to fail, and the failure branch must be spelled out
    assert "power <  0.80    §18.3 FAILS" in v5
    assert "**The gate may fail, and no step below is conditioned on its passing.**" in v5
    assert "**not** evidence against the hypothesis" in v5
    assert "ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE" in v5


@has_v5
def test_the_power_gate_cannot_be_reached_by_bending_the_cohort():
    v5 = _v5()
    forbidden = v5.split("## 9.4 What may still not be done to satisfy §9")[1].split("```")[1]
    for item in ("admitting the sorted samples S2/S3",
                 "switching a replicate to the other author linkage list",
                 "changing N, the tie rule or the ranking statistic",
                 "adding a replicate that failed §11 qualification",
                 "counting R1 positives",
                 "different seeds and keeping the better number",
                 "reporting SUPPORTED on a design whose measured power is below 0.80"):
        assert item in forbidden, item


@has_v5
def test_the_power_measurement_must_precede_the_confirmatory_statistics():
    v5 = _v5()
    assert "record power BEFORE any confirmatory statistic exists" in v5
    assert "23.2H-C may not be run after 23.2H-D" in v5


# ---- the rule conflict, resolved without opening a loophole ------------------------------------#
@has_v5
def test_source_faithfulness_replaces_r1s_parameters_and_says_why():
    v5 = _v5()
    assert "R1's top-100-with-ties is R1's source rule" in v5
    assert "not a universal constant" in v5
    assert ("freeze R1's **instance** of the principle in the place where the **principle** "
            "belonged") in _flat(v5)


@has_v5
def test_every_v4_anti_gaming_prohibition_survives_the_rule_change():
    v5 = _v5()
    forbidden = v5.split("**Anti-gaming constraints, retained in full.**")[1].split("```")[1]
    for item in ("changing N, the tie rule, the ranking statistic or the join",
                 "substituting one replicate's rule into another replicate",
                 "pooling outcome libraries across biological replicates before selection",
                 "re-ranking barcodes across units",
                 "because of how it moves positives or performance",
                 "admitting the sorted samples S2/S3"):
        assert item in forbidden, item


@has_v5
def test_the_rules_were_fixed_by_reproducing_author_objects_not_chosen():
    v5 = _v5()
    assert ("reproducing the authors' serialized objects exactly, before any confirmatory"
            ) in _flat(v5)
    assert "They cannot have been tuned to a result because no result existed." in _flat(v5)


@has_v5
def test_v5_argues_benchmark_compatibility_rather_than_assuming_it():
    v5 = _v5()
    compat = v5.split("## 6.1 Why this does not trigger")[1].split("## 6.2")[0]
    assert "R1's outcome rule            unchanged" in compat
    assert "the historical Stage-23 Role-A FAIL   permanent" in compat
    assert "It adds a cohort; it does not edit one." in _flat(compat)
    assert "never by mutating a Stage-22 file" in compat


@has_v5
def test_the_heterogeneous_endpoint_is_declared_as_a_limitation():
    v5 = _v5()
    lim = v5.split("## 6.2 The declared limitation this creates")[1].split("# 7.")[0]
    assert "three different source-defined operationalisations" in lim
    assert "may not be described as a uniform endpoint" in _flat(lim)
    assert "Prevalence heterogeneity" in lim


# ---- the author bug ----------------------------------------------------------------------------#
@has_v5
def test_the_spike_in_bug_decision_is_frozen_in_both_directions():
    v5 = _v5()
    assert "PRIMARY      corrected coefficients" in v5
    assert "SENSITIVITY  author coefficients" in v5
    assert ("may **not** be swapped for the sensitivity arm after results are seen, in either "
            "direction") in _flat(v5)
    assert "Both are run and both are reported" in v5


@has_v5
def test_the_eligibility_exclusion_is_argued_on_population_grounds():
    v5 = _v5()
    el = v5.split("# 4. Frozen eligibility")[1].split("# 5.")[0]
    assert "changes the **pre-state population**" in el
    assert "They may not be admitted later" in _flat(el)
    assert "GSM7092520" in el and "GSM7092521" in el


@has_v5
def test_the_fold_stratification_does_not_repeat_the_v3_self_contradiction():
    """V3 required a y-stratified table frozen before y was inspected. V5 must not restate that."""
    v5 = _v5()
    f = v5.split("## 8.1 Fold construction")[1].split("## 8.2")[0]
    assert "23511" in f
    assert "outcomes are reconstructed first" in f
    assert "no model has been fitted at either point" in f


@has_v5
def test_permutation_strata_carry_the_replicate_term():
    v5 = _v5()
    assert 'stratum = f"{size}|{n_lanes}|{replicate}"' in v5
    assert "prevents a profile from crossing a replicate boundary" in _flat(v5)
    # the null must also contain whatever replicate-identifying signal expression carries
    assert "is present in the **null** as well as in the observed statistic" in _flat(v5)


@has_v5
def test_replicate_identity_is_a_blocking_nuisance_not_a_predictor():
    """The two replicates differ 7-fold in prevalence; without blocking, batch signal reads as
    state signal."""
    v5 = _v5()
    b = v5.split("## 8.2 Replicate identity as a blocking nuisance covariate")[1].split("## 8.3")[0]
    assert "is_biological_replicate_3" in b
    assert "may not be a predictor of interest" in _flat(b)
    assert "may not be\ninteracted with `X`" in b or "may not be interacted with `X`" in _flat(b)
    assert "1.4%" in b and "10.1%" in b
    assert "The blocking term puts that information in the **baseline**" in _flat(b)


@has_v5
def test_the_stage_27_firewall_is_restated():
    v5 = _v5()
    assert "are **not** available as the untouched Stage-27 replication set" in _flat(v5)


@has_v5
def test_v5_adds_no_new_roadmap_exit():
    """The roadmap defines four exits. A protocol may not invent a fifth to escape a failed gate."""
    exits = set(re.findall(r"ROLE_A_[A-Z_]+", _v5()))
    assert exits <= {"ROLE_A_CONFIRMATORY_SUPPORTED", "ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE",
                     "ROLE_A_REDESIGN_REQUIRED", "ROLE_A_ABANDONED_NO_DEFENSIBLE_SIGNAL"}, exits
    assert "V5 adds none" in _v5()


# ============================================================================================== #
# 23.2H-A — the reconstruction is only credible if it reproduces the authors' own objects
# ============================================================================================== #
@ran_a
def test_r2_rule_reproduces_primedcellsind_exactly():
    v = _json(BENCH)["author_object_validation"]["R2"]
    assert v["lineage_sets_identical"] is True
    assert v["cell_sets_identical"] is True
    assert v["author_lineages"] == v["reconstructed_lineages"] == 26
    assert v["author_cells"] == v["reconstructed_cells"] == 79
    assert v["author_SampleNum_values"] == ["S4", "S5"], "R2's object must be replicate 2 only"
    for lib, c in v["coefficient_agreement"].items():
        assert c["relative_error"] < 1e-12, f"{lib} calibration disagrees with the author object"


@ran_a
def test_r3_author_bug_arm_reproduces_primedcellidlist_exactly():
    v = _json(BENCH)["author_object_validation"]["R3"]
    assert v["author_bug_arm_reproduces_author_object"] is True
    assert [u["author_cells"] for u in v["per_unit"]] == [27, 75, 83]
    for u in v["per_unit"]:
        assert u["identical"] is True, u["unit"]


@ran_a
def test_the_spike_in_bug_is_material_and_only_hits_units_2_and_3():
    """FS_1 is the one unit the author code scales correctly, so it must be untouched."""
    ov = _json(BENCH)["author_object_validation"]["R3"]["corrected_vs_bugged_lineage_overlap"]
    assert ov["FS_1"] == 200, "FS_1 uses the right coefficients in both arms"
    assert ov["FS_2"] < 200 and ov["FS_3"] < 200, "the bug must actually move lineages"


@ran_a
def test_the_frozen_cohort_is_what_v5_declares():
    b = _json(BENCH)
    assert b["arm"] == "PRIMARY"
    assert b["expected_vs_realized_mismatch"] == {}
    assert b["realized"]["2"] == {"cells": 3480, "clones": 1827, "positive_cells": 79,
                                  "positive_clones": 26, "prevalence": 0.014231}
    assert b["realized"]["3"]["clones"] == 483
    assert b["totals"]["clones"] == 2310
    assert b["eligibility"]["excluded"] == ["S2", "S3"]
    assert b["folds"]["seed"] == 23511
    assert b["folds"]["stratified_on"] == "(biological_replicate, y_primed)"


@ran_a
def test_every_outer_fold_carries_both_replicates():
    per_fold = _json(BENCH)["folds"]["per_fold"]
    seen = {}
    for row in per_fold:
        seen.setdefault(row["outer_fold"], set()).add(row["biological_replicate"])
    assert len(seen) == 5
    for f, reps in seen.items():
        assert reps == {"2", "3"}, f"fold {f} is missing a replicate: {reps}"


@ran_a
def test_no_sorted_sample_reached_the_cohort():
    import csv
    with open(OUT / "stage23_2h_confirmation_cells.csv", newline="", encoding="utf-8") as fh:
        samples = {r["SampleNum"] for r in csv.DictReader(fh)}
        gsms = set()
        fh.seek(0)
        gsms = {r["gsm"] for r in csv.DictReader(fh)}
    assert samples == {"S1", "S4", "S5"}, samples
    assert "GSM7092520" not in gsms and "GSM7092521" not in gsms


@ran_bug
def test_the_sensitivity_arm_differs_from_the_primary_by_exactly_the_bug():
    p, s = _json(BENCH), _json(BENCH_BUG)
    assert s["arm"] == "R3_MAX_PAIRED_TOP200_UNION_AUTHORBUG"
    # replicate 2 is untouched by the R3 spike-in decision
    assert p["realized"]["2"] == s["realized"]["2"]
    # the eligible cohort is outcome-free and therefore identical
    for rep in ("2", "3"):
        assert p["realized"][rep]["cells"] == s["realized"][rep]["cells"]
        assert p["realized"][rep]["clones"] == s["realized"][rep]["clones"]
    assert p["totals"]["positive_clones"] != s["totals"]["positive_clones"]


# ============================================================================================== #
# 23.2H-B — the representation
# ============================================================================================== #
@ran_b
def test_representation_is_outcome_free_and_normalised_once():
    r = _json(REPR)
    assert r["outcome_free"] is True
    assert r["normalization"] == ("sum raw counts per clone -> CP10K -> log1p, "
                                 "applied exactly once")
    assert r["matrix"]["clones"] == 2310
    assert r["matrix"]["genes"] == 36601, "WM989's Custom lineage features must never appear"
    assert r["all_positive_total_umi"] is True
    assert r["detected_matches_normalised_nonzero_pattern"] is True


# ============================================================================================== #
# 23.2H-C / D / E — the gates
# ============================================================================================== #
@ran_c
def test_power_was_measured_on_the_real_cohort_geometry():
    p = _json(POWER)
    assert p["cohort_geometry"]["clones"] == 2310
    assert p["target_AUC"] == 0.66
    assert p["threshold"] == 0.80
    assert abs(p["calibration"]["achieved_median_oracle_AUC"] - 0.66) < 0.01
    assert p["null"]["n"] == 200 and p["alternative"]["n"] == 100
    assert p["gate_18_3_measured_power"] is (p["power"] >= 0.80)


@ran_c
@ran_d
def test_the_power_gate_was_recorded_and_was_not_tuned_to_pass():
    """V5 §10 forbids running 23.2H-D before 23.2H-C, so the power gate cannot be adjusted after
    the confirmatory statistic is known.

    This was previously asserted from file modification times. **Git does not preserve mtime.** On
    any fresh checkout both files carry the checkout instant in arbitrary order, so the assertion
    was a coin flip that tested the filesystem it happened to run on rather than the evidence — and
    it failed CI on every single commit for exactly that reason.

    Run order is **not recoverable** from the committed artifacts: neither JSON carries a timestamp,
    and only the confirmation records a git commit. Rather than assert something the repository
    cannot support, this checks what the ordering rule was there to protect — that the power gate is
    recorded, and that it **FAILED**, so it demonstrably was not moved to make the design look
    adequate. The run order itself is documented in `stage_23_2H_RECORD.md`; that is a record of the
    claim, not proof of it, and it is not pretended otherwise here.
    """
    p = _json(POWER)
    assert p["verdict"] == "DESIGN_UNDERPOWERED"
    assert p["gate_18_3_measured_power"] is False
    assert p["power"] < p["threshold"], "a gate that failed cannot have been tuned to pass"
    # and the confirmatory run is pinned to the frozen protocol, so it did not run against a revision
    assert _json(CONFIRM)["protocol"]["file"] == "STAGE_23_2_ROLE_A_CONFIRMATION_V5.md"
    assert len(_json(CONFIRM)["protocol"]["canonical_lf_sha256"]) == 64


@ran_d
def test_the_confirmatory_test_uses_the_frozen_dual_gate():
    c = _json(CONFIRM)
    perm = c["permutation"]
    # The primary arm was pre-committed to 2000 draws before any value beyond 200 was inspected.
    # Pinning the committed number is what makes that commitment checkable.
    assert perm["n_permutations"] == 2000
    assert perm["strata"] == "size{1,2,3+} x n_lanes x biological_replicate"
    assert c["gate_18_4_pooled"] is bool(perm["exceeds_null_p95"] and perm["p_perm"] <= 0.05)
    assert set(c["per_replicate"]) == {"2", "3"}
    assert c["gate_18_5_every_replicate_positive"] is bool(
        all(v["delta_AP"] > 0 for v in c["per_replicate"].values()))


@ran_e
def test_the_verdict_follows_mechanically_from_the_gates():
    v = _json(VERDICT)
    assert set(v["gates"]) == {
        "18.1_two_independent_non_R1_biological_replicates",
        "18.2_source_faithful_reconstruction",
        "18.3_measured_design_power",
        "18.4_pooled_dual_gate",
        "18.5_every_replicate_positive",
        "18.6_no_un_re_gated_material_benchmark_change"}
    if all(v["gates"].values()):
        assert v["exit"] == "ROLE_A_CONFIRMATORY_SUPPORTED" and v["stage_24"] == "MAY_OPEN"
    else:
        assert v["exit"] == "ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE"
        assert v["stage_24"] == "BLOCKED"
        assert v["failing_gates"], "a blocked exit must name which gate failed"


@ran_e
def test_the_standing_limitations_are_carried_into_the_verdict():
    lim = " ".join(_json(VERDICT)["standing_limitations"])
    assert "three different source-defined operationalisations" in lim
    assert "NOT available to Stage 27" in lim


# ============================================================================================== #
# The sharded permutation engine.
#
# A 2000-draw run is long enough that it WILL be interrupted, and the pre-commitment attached to it
# only means something if sharding and resuming cannot change the answer. The smoke stage proves
# that against a real sequential computation; these contracts pin its verdict and the audit trail
# around the decision to extend the run.
# ============================================================================================== #
SMOKE = OUT / "stage23_2h_smoke.json"
N200 = OUT / "stage23_2h_confirmation_n200.json"
N200_BUG = OUT / "stage23_2h_confirmation_authorbug_n200.json"

ran_smoke = pytest.mark.skipif(not SMOKE.exists(), reason="the smoke stage has not been run")
has_n200 = pytest.mark.skipif(not N200.exists(), reason="no 200-permutation predecessor on disk")


@ran_smoke
def test_sharding_and_resume_cannot_change_the_answer(s=None):
    s = _json(SMOKE)
    by = {c["check"]: c for c in s["checks"]}
    assert by["3 shards, run out of order, are BIT-IDENTICAL to sequential"]["pass"] is True
    assert by["interrupt-and-resume reproduces sequential exactly"]["pass"] is True
    assert by["mixed-protocol cache is refused"]["pass"] is True
    assert s["all_passed"] is True


@ran_smoke
def test_the_eta_is_measured_not_guessed():
    s = _json(SMOKE)
    assert s["draws_used"] >= 6
    assert s["measured_seconds_per_draw"] > 0
    assert "n_perm=2000" in s["eta"]
    # the honest caveat about parallel contention must ride along with the estimate
    assert "contend" in s["note"]


@has_n200
def test_the_200_permutation_result_is_preserved_before_any_extension():
    """The pre-commitment is only auditable if the superseded result cannot quietly vanish."""
    d = _json(N200)
    assert d["permutation"]["n_permutations"] == 200
    assert d["arm"] == "PRIMARY"
    # the exact numbers that were on the table when the decision to extend was made
    assert abs(d["pooled"]["delta_AP"] - 0.030504) < 1e-5
    assert abs(d["permutation"]["p_perm"] - 0.029851) < 1e-5
    if N200_BUG.exists():
        b = _json(N200_BUG)
        assert b["permutation"]["n_permutations"] == 200
        assert abs(b["permutation"]["p_perm"] - 0.004975) < 1e-5


@ran_d
def test_an_extended_run_reports_its_own_permutation_count():
    """Whatever n_perm the reported result used, it must say so rather than inherit 200."""
    c = _json(CONFIRM)
    n = c["permutation"]["n_permutations"]
    assert n >= 200
    # p_perm must be consistent with (1 + #{null >= obs}) / (n + 1)
    expected = (1 + c["permutation"]["n_null_ge_observed"]) / (n + 1)
    assert abs(c["permutation"]["p_perm"] - expected) < 1e-9


# ============================================================================================== #
# The power CURVE. This is the single most gameable thing added after a failed gate, so the
# contracts here are about what the curve may NOT do.
# ============================================================================================== #
has_add1 = pytest.mark.skipif(not ADDENDUM1.exists(), reason="addendum 1 not written")


@has_add1
def test_v5_itself_was_not_edited_to_add_the_curve():
    """A frozen protocol may not be edited mid-run. The addendum must pin V5's digest."""
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    from run_stage23_learnability_gate import canonical_text_sha256
    a = ADDENDUM1.read_text(encoding="utf-8")
    assert "V5 itself is **not edited**" in a
    assert "REPORTING ONLY" in a
    pinned = re.search(r"`([0-9a-f]{64})`", a).group(1)
    assert canonical_text_sha256(V5) == pinned, "V5 has moved since the addendum pinned it"


@has_add1
def test_the_power_curve_is_reported_and_never_gating():
    c = ADDENDUM1.read_text(encoding="utf-8")
    assert "at **no other** value" in _flat(c) or "and **at no other\nvalue**" in c
    assert "Gate 18.3 is evaluated at 0.66 only." in c
    assert "It does not open" in _flat(c)


@has_add1
def test_the_curve_cannot_be_used_to_manufacture_a_pass():
    c = ADDENDUM1.read_text(encoding="utf-8")
    forbidden = c.split("§9.4 stands unchanged")[1].split("```")[1]
    for item in ("evaluating gate 18.3 at any point other than 0.66",
                 "reporting SUPPORTED because power at some larger alternative reaches 0.80",
                 "selecting a curve point after seeing which one clears the threshold",
                 "observed-power fallacy"):
        assert item in forbidden, item


@has_add1
def test_the_gate_anchor_is_inherited_not_chosen_now():
    a = ADDENDUM1.read_text(encoding="utf-8")
    assert "0.6628" in a, "the 0.66 anchor must name the historical value it matches"
    assert "V5 does not move it" in _flat(a)


@ran_c
def test_the_recorded_gate_result_is_at_the_pre_registered_point():
    p = _json(POWER)
    assert p["target_AUC"] == 0.66
    assert p.get("is_the_pre_registered_gate_point", True) is True
