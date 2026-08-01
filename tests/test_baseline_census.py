"""STAGE 1.5.2 gates G-a and G-b — the baseline census and donor-age wiring.

Both gates are **record-only**. The single most important test in this file is
`test_delta_age_is_bit_identical_with_and_without_the_census`: the plan's hard guard is that
"ΔAge values must come out **bit-identical** before/after. It records, it does not compute. If any
ΔAge moves, the change is wrong — revert, do not rationalise." That is asserted here, not assumed
from reading the diff.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
# `plan_tests/` holds the per-stage verification gates (verify_1a, verify_stage1_5, smoke_stage1).
# They are scripts, not a package, so the directory goes on the path to import from them.
for _p in (ROOT, ROOT / "src", ROOT / "plan_tests"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from verify_stage1_5 import ChunkControlStat, decide_verdict  # noqa: E402

from cellfate.data.aging import (  # noqa: E402
    LinearClock,
    _control_baseline,
    census_warnings,
    delta_age,
)
from cellfate.data.sources import ReprogrammingSource, _maybe_float  # noqa: E402


def _obs(lines, ctrl, **extra):
    d = {"cell_line": list(lines), "is_control": list(ctrl)}
    d.update({k: list(v) for k, v in extra.items()})
    return pd.DataFrame(d)


# ------------------------------------------------------------------ THE HARD GUARD ---- #
def test_delta_age_is_bit_identical_with_and_without_the_census():
    """G-a records; it must not compute. Any drift here means the change is wrong."""
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(8)]
    clock = LinearClock({g: float(w) for g, w in zip(genes, rng.normal(size=8), strict=True)},
                        intercept=41.5)
    expr = rng.normal(3.0, 1.0, size=(12, 8))
    obs = _obs(["A"] * 6 + ["B"] * 6, [True, False, False, True, False, False] * 2,
               batch=["E1", "E2"] * 6, donor_age=[53.0] * 12)
    d_plain, m_plain, _r_plain = delta_age(clock, expr, genes, obs, source="reprogramming")
    census: dict = {}
    d_census, m_census, _r_census = delta_age(clock, expr, genes, obs,
                                              source="reprogramming", census=census)
    assert np.array_equal(d_plain, d_census)          # bit-identical, not "close"
    assert np.array_equal(m_plain, m_census)
    assert census, "the census must actually have been filled, or the test proves nothing"


def test_control_baseline_values_do_not_depend_on_the_composition_argument():
    values = np.array([10.0, 20.0, 30.0, 40.0])
    lines = np.array(["A", "A", "B", "B"])
    ctrl = np.array([True, False, True, False])
    a = _control_baseline(values, lines, ctrl)
    b = _control_baseline(values, lines, ctrl, census={},
                          composition={"batch": np.array(["x", "y", "x", "y"])})
    assert np.array_equal(a, b)


# ---------------------------------------------------------------------- the census ---- #
def test_census_records_the_unreplicated_case_that_was_silent_before():
    values = np.array([10.0, 20.0, 30.0])
    census: dict = {}
    _control_baseline(values, np.array(["A"] * 3), np.array([True, False, False]), census=census)
    assert census["A"]["n_control"] == 1
    assert census["A"]["unreplicated"] is True
    assert census["A"]["source"] == "controls"


def test_census_records_the_no_control_fallback():
    census: dict = {}
    _control_baseline(np.array([1.0, 2.0]), np.array(["A", "A"]), np.array([False, False]),
                      census=census)
    assert census["A"]["source"] == "self_fallback"
    assert census["A"]["n_control"] == 0


def test_census_records_only_the_batches_the_baseline_itself_came_from():
    """Finding D1 in one assertion: controls all Exp2 while the line spans Exp1 and Exp2."""
    census: dict = {}
    _control_baseline(
        np.array([1.0, 2.0, 3.0, 4.0]), np.array(["A"] * 4), np.array([False, True, False, True]),
        census=census, composition={"batch": np.array(["Exp1", "Exp2", "Exp1", "Exp2"])})
    assert census["A"]["batch"] == ["Exp2"]
    assert census["A"]["batch_in_line"] == ["Exp1", "Exp2"]
    assert census["A"]["n_cells"] == 4 and census["A"]["n_control"] == 2


def test_census_warnings_name_each_problem_exactly_once():
    census = {
        "ok": {"n_control": 5, "n_cells": 20, "source": "controls", "unreplicated": False},
        "solo": {"n_control": 1, "n_cells": 20, "source": "controls", "unreplicated": True},
        "none": {"n_control": 0, "n_cells": 20, "source": "self_fallback", "unreplicated": False},
    }
    w = " | ".join(census_warnings(census))
    assert "solo" in w and "n=1" in w
    assert "none" in w and "NO controls" in w
    assert "ok:" not in w


def test_census_warning_flags_a_single_batch_baseline_under_a_multi_batch_line():
    census = {"A": {"n_control": 2, "n_cells": 20, "source": "controls", "unreplicated": False,
                    "batch": ["Exp2"], "batch_in_line": ["Exp1", "Exp2"]}}
    assert any("cross-batch" in x for x in census_warnings(census))


def test_a_column_that_is_constant_within_a_line_never_warns():
    """`donor_age` is a per-donor constant, so a single baseline value is the ONLY possible
    answer. An earlier version warned on every donor -- noise that trains the reader to
    ignore the warnings that matter."""
    census = {"A": {"n_control": 2, "n_cells": 20, "source": "controls", "unreplicated": False,
                    "donor_age": ["53.0"], "donor_age_in_line": ["53.0"],
                    "batch": ["Exp2"], "batch_in_line": ["Exp1", "Exp2"]}}
    w = census_warnings(census)
    assert any("cross-batch" in x for x in w)
    assert not any("donor_age" in x for x in w)


# --------------------------------------------------- verify_stage1_5's extended census ---- #
def test_the_stage_1_5_pass_rule_is_unchanged_by_the_new_flags():
    """Four runs are recorded against what a Stage 1.5 PASS means. It must not move."""
    stats = [ChunkControlStat("c1", "A", 20, 1, control_batches=("E2",),
                              cell_batches=("E1", "E2")),
             ChunkControlStat("c2", "B", 20, 5)]
    v = decide_verdict(stats)
    assert v["status"] == "PASS"                      # still PASS: no fallback fired
    assert v["baseline_warnings"]                     # but the problems are now visible
    assert v["unreplicated_chunks"][0]["cell_line"] == "A"
    assert v["cross_batch_chunks"][0]["control_batches"] == ["E2"]


def test_a_fallback_chunk_still_fails():
    v = decide_verdict([ChunkControlStat("c1", "A", 20, 0)])
    assert v["status"] == "FAIL"


def test_no_batch_information_is_not_reported_as_verified_single_batch():
    """A source that does not stamp `batch` must not be silently cleared of D1."""
    s = ChunkControlStat("c1", "A", 20, 5)            # no batch info at all
    assert s.cross_batch_baseline is False
    assert decide_verdict([s])["cross_batch_chunks"] == []


# ------------------------------------------------------------------------ G-b wiring ---- #
@pytest.mark.parametrize("raw,expect", [("53", 53.0), (" 0 ", 0.0), ("35.0", 35.0),
                                        ("", None), ("N/A", None), (None, None)])
def test_maybe_float_never_defaults_a_missing_age_to_zero(raw, expect):
    """0.0 is a REAL age here (N2/N3 are neonatal), so a missing value must read as None."""
    assert _maybe_float(raw) == expect


def test_build_chunk_carries_extra_metadata_into_obs():
    raw = ReprogrammingSource.build_chunk(
        "t:A", np.zeros((3, 2)), ["G0", "G1"], "A",
        ["control", "OSKM", "OSKM"], [0.0, 24.0, 48.0],
        extra={"donor_age": [53.0] * 3, "batch": ["Exp2", "Exp1", "Exp2"]})
    assert list(raw.obs["donor_age"]) == [53.0] * 3
    assert list(raw.obs["batch"]) == ["Exp2", "Exp1", "Exp2"]


def test_build_chunk_rejects_a_mis_sized_extra_column():
    """Silently broadcasting or truncating would put the wrong age on a cell."""
    with pytest.raises(Exception, match="extra column"):
        ReprogrammingSource.build_chunk(
            "t:A", np.zeros((3, 2)), ["G0", "G1"], "A",
            ["control", "OSKM", "OSKM"], [0.0, 24.0, 48.0], extra={"donor_age": [53.0]})


def test_donor_age_and_batch_are_metadata_not_model_input():
    """They must never reach the model: the deployed request schema forbids extra fields."""
    from pydantic import ValidationError

    from cellfate.inference.schema import Request
    with pytest.raises(ValidationError, match="donor_age"):
        Request(X_raw=[0.0], u_modality="tf", u_descriptor="OSKM", dose_uM=1.0, time_h=24.0,
                donor_age=53.0)


# ------------------------------------------------- REGRESSIONS (found in review) ---- #
def test_census_keys_must_survive_one_cell_line_spanning_many_chunks():
    """The collision that silently discarded 44 of HFF's 45 chunks.

    `cell_line` is NOT unique across chunks: `verify_stage1_5_results.json` records HFF in 45
    of them. `run()` merged each chunk's census with `baseline_census.update(chunk_census)`,
    keyed on the line alone, so every chunk overwrote the previous one and the manifest kept a
    single record -- for the dataset carrying ~99.8% of the age labels. A baseline problem in
    any chunk but the last became invisible, which is precisely what G-a exists to prevent.

    This asserts the merge policy `run()` uses, on the shape that actually broke it.
    """
    per_chunk = {}
    for i in range(3):
        c: dict = {}
        obs = _obs(["HFF"] * 4, [True, False, False, False])
        clock = LinearClock({"G0": 1.0}, intercept=0.0)
        delta_age(clock, np.arange(4.0).reshape(4, 1), ["G0"], obs,
                  source="reprogramming", census=c)
        per_chunk[f"reprogramming:HFF:b{i}"] = c

    merged_buggy: dict = {}
    for c in per_chunk.values():
        merged_buggy.update(c)
    assert len(merged_buggy) == 1, "precondition: keying on cell_line alone collides"

    merged: dict = {}
    for cid, c in per_chunk.items():
        for line, rec in c.items():
            merged[f"{cid}::{line}"] = {**rec, "chunk_id": cid, "cell_line": line}

    assert len(merged) == 3, "every chunk must keep its own census record"
    assert {r["chunk_id"] for r in merged.values()} == set(per_chunk)
    assert {r["cell_line"] for r in merged.values()} == {"HFF"}


def test_census_warnings_ignores_the_chunk_id_and_cell_line_fields():
    """The namespaced record gained two str fields; neither may be read as a composition list."""
    census = {"reprogramming:HFF:b0::HFF": {
        "n_control": 1, "n_cells": 4, "source": "controls", "unreplicated": True,
        "chunk_id": "reprogramming:HFF:b0", "cell_line": "HFF"}}
    w = census_warnings(census)
    assert len(w) == 1 and "n=1" in w[0]


def test_render_handles_an_errored_chunk_without_crashing():
    """G-a widened the table to six columns; the error branch still appended five.

    `render_table` indexes `row[i] for i in range(len(headers))`, so the short row raised
    IndexError -- crashing the renderer on the one path `scan_build` deliberately survives.
    """
    from verify_stage1_5 import _render

    stats = [ChunkControlStat("c0", "A", 0, 0, error="RuntimeError('boom')"),
             ChunkControlStat("c1", "B", 10, 2)]
    _render(stats, decide_verdict(stats))


# ============================================================================ #
# STAGE 1.5.3 STEP 1 — C-6 (`age_mask_reason`) and C-3 (HFF donor metadata)    #
# ============================================================================ #
def test_sample_rejects_a_reason_on_an_unmasked_row():
    """C-6's invariant: `reason is None` exactly when `age_mask` is True. Everything
    downstream is allowed to rely on it, so it is enforced, not documented."""
    from pydantic import ValidationError

    from cellfate.common.schemas import Sample
    kw = dict(cell_id="c", X=[1.0], u_modality="tf", u_tf_emb=[0.0], dose_time=[0.0, 0.0],
              y_cls=[1.0, 0.0, 0.0], sig_scores=[0.0, 0.0, 0.0], cell_line="A",
              pert_id="p", scaffold_id="s", source="synth")
    with pytest.raises(ValidationError, match="age_mask_reason to be None"):
        Sample(y_age=1.0, age_mask=True, age_mask_reason="cancer_source", **kw)
    # the two legal shapes
    assert Sample(y_age=1.0, age_mask=True, **kw).age_mask_reason is None
    assert Sample(y_age=None, age_mask=False, age_mask_reason="cancer_source",
                  **kw).age_mask_reason == "cancer_source"


def test_assemble_defaults_reasons_to_none_and_mirrors_the_mask():
    from cellfate.common.constants import Modality
    from cellfate.data.assemble import assemble_samples
    n = 3
    common = dict(cell_ids=[f"c{i}" for i in range(n)], x_panel=np.zeros((n, 2)),
                  fingerprints=np.zeros((n, 0), dtype=np.uint8), dose_time=np.zeros((n, 2)),
                  y_cls=np.tile([1.0, 0.0, 0.0], (n, 1)), y_age=np.arange(n, dtype=float),
                  sig_scores=np.zeros((n, 3)), cell_line=["A"] * n, pert_id=["p"] * n,
                  scaffold_id=["s"] * n, source="synth", modality=Modality.TF,
                  tf_emb=np.zeros((n, 8)))
    # no reasons supplied -> all None, and masked rows still validate
    s = assemble_samples(age_mask=np.array([True, True, True]), **common)
    assert [x.age_mask_reason for x in s] == [None] * n
    # reasons supplied -> kept only where the mask is False
    s = assemble_samples(age_mask=np.array([True, False, False]),
                         age_mask_reason=[None, "dataset_policy", "cancer_source"], **common)
    assert [x.age_mask_reason for x in s] == [None, "dataset_policy", "cancer_source"]
    assert [x.y_age for x in s] == [0.0, None, None]


def test_assemble_rejects_a_mis_sized_reason_list():
    """Silently truncating would attach the wrong reason to the wrong cell."""
    from cellfate.common.constants import Modality
    from cellfate.data.assemble import assemble_samples
    with pytest.raises(ValueError, match="age_mask_reason has"):
        assemble_samples(cell_ids=["a", "b"], x_panel=np.zeros((2, 2)),
                         fingerprints=np.zeros((2, 0), dtype=np.uint8),
                         dose_time=np.zeros((2, 2)), y_cls=np.tile([1.0, 0.0, 0.0], (2, 1)),
                         y_age=np.zeros(2), age_mask=np.array([False, False]),
                         age_mask_reason=["only_one"], sig_scores=np.zeros((2, 3)),
                         cell_line=["A"] * 2, pert_id=["p"] * 2, scaffold_id=["s"] * 2,
                         source="synth", modality=Modality.TF, tf_emb=np.zeros((2, 8)))


def test_shard_reader_tolerates_a_shard_written_before_c6():
    """The committed shards in runs/ predate the column. Requiring it would break
    training/dataset.py, evaluation/data.py, inference/service.py and three runners for no
    benefit while the masking policies are off -- see the comment at io.shard_to_numpy."""
    import pyarrow as pa

    from cellfate.common import io
    # one real row, minus the new column -- an EMPTY table would exercise numpy reshape edge
    # cases rather than the tolerance this test is about
    row = {"cell_id": ["c"], "X": [[0.0, 0.0]], "u_modality": ["tf"], "u_chem_fp": [None],
           "u_gene_emb": [None], "u_tf_emb": [[0.0]], "dose_time": [[0.0, 0.0]],
           "y_cls": [[1.0, 0.0, 0.0]], "y_age": [1.0], "age_mask": [True],
           "sig_scores": [[0.0, 0.0, 0.0]], "cell_line": ["A"], "pert_id": ["p"],
           "scaffold_id": ["s"], "source": ["synth"]}
    old_schema = pa.schema([f for f in io.SHARD_SCHEMA if f.name != "age_mask_reason"])
    out = io.shard_to_numpy(pa.table(row, schema=old_schema))
    assert out["age_mask_reason"] == [None]          # tolerated, not a KeyError
    assert out["age_mask"].tolist() == [True]        # and the rest still reads


def test_hff_asserts_a_neonatal_donor_age_with_its_provenance_recorded():
    """C-3. The value is asserted, not parsed, so it must be visible and explained."""
    from cellfate.data.sources import GSE242423SingleCellSource as S
    assert S.DONOR_AGE_YEARS == 0.0
    assert "asserted" in S.DONOR_AGE_PROVENANCE and "not in GEO" in S.DONOR_AGE_PROVENANCE


def test_hff_donor_age_sits_outside_the_shipped_clocks_fitted_range():
    """The whole point of C-3: without it, C-2's range rule cannot fire for 99.7% of the data."""
    import json
    from pathlib import Path

    from cellfate.data.sources import GSE242423SingleCellSource as S
    meta = json.loads((Path(__file__).resolve().parents[1] / "configs" / "clocks" /
                       "fleischer_clock.json").read_text(encoding="utf-8"))["meta"]
    lo, hi = meta["age_range"]
    assert S.DONOR_AGE_YEARS < lo, f"{S.DONOR_AGE_YEARS} should be below the clock's {lo}"
    assert lo == 1.0 and hi == 96.0


# ============================================================================ #
# STAGE 1.5.3 STEP 2 — C-1 (`age_label_policy`)                                #
# ============================================================================ #
def test_default_policy_is_bit_identical_to_the_old_expression():
    """THE GUARD. Rules 2 and 3 are off by default, so nothing may move -- even on cells that
    WOULD qualify under them."""
    from cellfate.data.aging import age_label_policy
    obs = _obs(["A"] * 6, [True] + [False] * 5,
               dataset_id=["hff_sc"] * 6, donor_age=[0.0] * 6)
    mask, reasons = age_label_policy(6, "reprogramming", obs)
    assert np.array_equal(mask, np.full(6, True))
    assert reasons == [None] * 6


def test_a_masked_dataset_is_masked_and_its_neighbour_in_the_same_chunk_is_not():
    """The blocking capability C-1 exists for: HFF and Gill separable inside ONE chunk.
    Before this, `age_mask` keyed on `source` alone and both report "reprogramming"."""
    from cellfate.data.aging import age_label_policy
    obs = _obs(["HFF"] * 2 + ["O1"] * 2, [False] * 4,
               dataset_id=["hff_sc", "hff_sc", "gill_bulk", "gill_bulk"])
    mask, reasons = age_label_policy(4, "reprogramming", obs,
                                     masked_datasets=frozenset({"hff_sc"}))
    assert list(mask) == [False, False, True, True]
    assert reasons == ["dataset_policy", "dataset_policy", None, None]


def test_the_cancer_rule_wins_and_is_reported_first():
    """Order stability: a cell excluded twice reports the first rule, so the string does not
    depend on the order the later rules happen to be written in."""
    from cellfate.data.aging import age_label_policy
    obs = _obs(["A"] * 2, [False] * 2, dataset_id=["hff_sc"] * 2)
    mask, reasons = age_label_policy(2, "tahoe", obs, masked_datasets=frozenset({"hff_sc"}))
    assert not mask.any()
    assert reasons == ["cancer_source"] * 2


def test_an_unknown_donor_age_never_masks():
    """Absence of evidence is recorded, not acted on. NaN must not silently exclude."""
    from cellfate.data.aging import age_label_policy
    obs = _obs(["A"] * 3, [False] * 3, donor_age=[float("nan"), 53.0, 0.0])
    mask, reasons = age_label_policy(3, "reprogramming", obs, clock_age_range=(1.0, 96.0))
    assert list(mask) == [True, True, False]
    assert reasons == [None, None, "donor_out_of_clock_range"]


def test_masked_datasets_without_the_column_raises_rather_than_keeping_labels():
    """BUG 3 from the cross-review: the original spec read `if masked_datasets and
    "dataset_id" in obs.columns`, which silently KEEPS labels meant to be withheld -- the
    unsafe direction, and invisible."""
    from cellfate.data.aging import age_label_policy
    with pytest.raises(KeyError, match="dataset_id"):
        age_label_policy(3, "reprogramming", pd.DataFrame({"cell_line": ["a", "b", "c"]}),
                         masked_datasets=frozenset({"hff_sc"}))


def test_clock_age_range_without_donor_age_raises():
    from cellfate.data.aging import age_label_policy
    with pytest.raises(KeyError, match="donor_age"):
        age_label_policy(3, "reprogramming", pd.DataFrame({"cell_line": ["a", "b", "c"]}),
                         clock_age_range=(1.0, 96.0))


def test_a_missing_column_is_not_an_error_when_the_policy_is_off():
    """The distinction that makes the two raises above safe: OFF means inapplicable, not wrong."""
    from cellfate.data.aging import age_label_policy
    mask, reasons = age_label_policy(3, "reprogramming",
                                     pd.DataFrame({"cell_line": ["a", "b", "c"]}))
    assert mask.all() and reasons == [None] * 3


def test_the_reason_is_none_exactly_where_the_mask_is_true():
    """The invariant `Sample` validation depends on (schemas.py)."""
    from cellfate.data.aging import age_label_policy
    obs = _obs(["A"] * 4, [False] * 4, donor_age=[0.0, 53.0, 200.0, 30.0])
    mask, reasons = age_label_policy(4, "reprogramming", obs, clock_age_range=(1.0, 96.0))
    assert [r is None for r in reasons] == list(mask)


def test_delta_age_returns_reasons_and_stays_bit_identical_at_defaults():
    """C-1 widens delta_age to a 3-tuple. At defaults the first two must be unchanged."""
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(8)]
    clock = LinearClock({g: float(w) for g, w in zip(genes, rng.normal(size=8), strict=True)},
                        intercept=41.5)
    expr = rng.normal(3.0, 1.0, size=(12, 8))
    obs = _obs(["A"] * 6 + ["B"] * 6, [True, False, False, True, False, False] * 2,
               dataset_id=["hff_sc"] * 12, donor_age=[0.0] * 12)
    d, mask, reasons = delta_age(clock, expr, genes, obs, source="reprogramming")
    assert len(reasons) == 12 and set(reasons) == {None}
    assert mask.all()
    # and the arithmetic itself is untouched: ΔAge is age minus the per-line control mean
    age = clock.predict_age(expr, genes)
    want = age - _control_baseline(age, obs["cell_line"].to_numpy(),
                                   obs["is_control"].to_numpy().astype(bool))
    assert np.array_equal(d, want)


def test_cancer_sources_are_still_masked_through_delta_age():
    """The pre-existing rule, exercised end to end rather than only in the pure helper."""
    rng = np.random.default_rng(1)
    genes = ["G0", "G1"]
    clock = LinearClock({"G0": 1.0, "G1": 0.5}, intercept=10.0)
    obs = _obs(["A"] * 4, [True, False, False, False])
    _d, mask, reasons = delta_age(clock, rng.normal(size=(4, 2)), genes, obs, source="tahoe")
    assert not mask.any()
    assert reasons == ["cancer_source"] * 4


# ============================================================================ #
# STAGE 1.5.3 STEP 3 — C-2 (the clock's declared age_range, carried not dropped)#
# ============================================================================ #
def test_linear_clock_carries_the_fitted_age_range_from_its_metadata():
    """C-2. `clock_fit.py` WRITES `meta.age_range` and, until now, nothing read it --
    `grep -rn "age_range" src/` returned the single write and no read."""
    clock = LinearClock.from_json(ROOT / "configs" / "clocks" / "fleischer_clock.json")
    assert clock.age_range == (1.0, 96.0)


def test_a_clock_without_metadata_reports_no_range_rather_than_guessing():
    """A clock whose provenance is unknown must not silently claim a validity range."""
    assert LinearClock({"G0": 1.0}).age_range is None


def test_the_range_is_carried_but_never_enforced_by_the_clock_itself():
    """The clock reports what it was fitted on; the POLICY of what to do about
    extrapolation belongs to the label pipeline, not to the clock."""
    clock = LinearClock({"G0": 1.0}, intercept=0.0, age_range=(1.0, 96.0))
    # a neonatal donor still gets an age predicted -- masking is age_label_policy's job
    assert np.isfinite(clock.predict_age(np.array([[5.0]]), ["G0"])[0])


def test_the_range_rule_is_off_by_default_in_dataconfig():
    """Turning it on MOVES LABELS -- N2 and N3 are donor_age 0, below the clock's 1.0 --
    so it is a pre-registered change, never a default."""
    from cellfate.data.build_dataset import DataConfig
    cfg = DataConfig(out="x", gene_panel="y")
    assert cfg.enforce_clock_age_range is False


def test_the_range_rule_masks_the_neonatal_donors_when_switched_on():
    """The consequence, measured rather than assumed: at [1, 96] the age-0 donors go."""
    from cellfate.data.aging import age_label_policy
    obs = _obs(["N2", "N3", "Y1", "O1"], [False] * 4, donor_age=[0.0, 0.0, 29.0, 53.0])
    mask, reasons = age_label_policy(4, "reprogramming", obs, clock_age_range=(1.0, 96.0))
    assert list(mask) == [False, False, True, True]
    assert reasons[:2] == ["donor_out_of_clock_range"] * 2


def test_delta_age_leaves_the_range_rule_off_unless_a_range_is_passed():
    """The step-3 guard in one assertion: the plumbing exists and is inert by default."""
    genes = ["G0"]
    clock = LinearClock({"G0": 1.0}, intercept=0.0, age_range=(1.0, 96.0))
    obs = _obs(["N2"] * 3, [True, False, False], donor_age=[0.0] * 3)
    expr = np.array([[1.0], [2.0], [3.0]])
    # default: the clock KNOWS its range, and delta_age still does not use it
    _d, mask, reasons = delta_age(clock, expr, genes, obs, source="reprogramming")
    assert mask.all() and reasons == [None] * 3
    # opt in explicitly and the same cells are masked
    _d, mask, reasons = delta_age(clock, expr, genes, obs, source="reprogramming",
                                  clock_age_range=clock.age_range)
    assert not mask.any()
    assert reasons == ["donor_out_of_clock_range"] * 3
