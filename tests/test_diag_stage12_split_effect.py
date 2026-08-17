"""Unit tests for the Stage 12 split-effect measurement.

This script re-derives a split from stored rows, so its whole value rests on the derivation being
faithful. The canary (rebuilt map == stored map) is therefore tested in BOTH directions: that it
passes on the real fold, and that a corrupted map actually aborts the run rather than being
reported around.

The real-fold tests are marked slow-ish (they read 51 shards once, via a session fixture).
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments" / "diag_stage12_split_effect.py"
spec = importlib.util.spec_from_file_location("s12e", SRC)
s12 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s12)

FOLD_EXISTS = (ROOT / s12.FOLD / "manifest.parquet").exists()
needs_fold = pytest.mark.skipif(not FOLD_EXISTS, reason="built fold not present")


@pytest.fixture(scope="module")
def result():
    return s12.run()


# ---- the keys ------------------------------------------------------------------------------- #
def test_the_new_key_is_unique_where_the_old_one_collides():
    man = pd.DataFrame({"cell_id": ["reprogramming:HFF:0"] * 3,
                        "shard_id": ["b0", "b1", "b2"], "row_idx": [0, 0, 0]})
    assert len(set(s12.old_key(man))) == 1
    assert len(set(s12.new_key(man))) == 3


def test_the_new_key_carries_the_chunk_and_the_row():
    man = pd.DataFrame({"cell_id": ["x"], "shard_id": ["reprogramming_HFF_b7"], "row_idx": [42]})
    assert s12.new_key(man) == ["reprogramming_HFF_b7:42"]


# ---- the split derivation ------------------------------------------------------------------- #
def test_a_colliding_key_yields_fewer_map_entries_than_cells():
    keys = [f"reprogramming:HFF:{i % 10}" for i in range(100)]
    smap = s12.build_split(keys, ["HFF"] * 100)
    assert len(smap) == 10


def test_the_held_out_line_goes_entirely_to_test():
    keys = [f"k{i}" for i in range(20)]
    lines = ["HFF"] * 10 + [s12.HOLDOUT_LINE] * 10
    smap = s12.build_split(keys, lines)
    assert all(smap[f"k{i}"] == "test" for i in range(10, 20))
    assert all(smap[f"k{i}"] != "test" for i in range(10))


def test_assign_looks_the_split_up_per_cell_which_is_what_gather_split_does():
    """With a colliding key, every cell sharing an id inherits ONE decision."""
    man = pd.DataFrame({"day": ["D0"] * 4})
    keys = ["a", "a", "b", "b"]
    assert s12.assign(man, keys, {"a": "train", "b": "calib"}) == \
        ["train", "train", "calib", "calib"]


def test_composition_counts_d0_share_per_split():
    man = pd.DataFrame({"day": ["D0", "D0", "later", "later"]})
    comp = s12.composition(man, ["train", "later_split", "train", "later_split"])
    assert comp["train"] == {"n": 2, "d0_n": 1, "d0_share": 0.5}


# ---- the canary, in both directions --------------------------------------------------------- #
@needs_fold
def test_the_reconstruction_reproduces_the_real_split_map_exactly(result):
    """Without this the rest is fiction. 1100 entries for 42,600 cells IS the defect, read back
    off the artefact the build actually wrote."""
    c = result["canary"]
    assert c["identical"] is True
    assert c["stored_entries"] == c["rebuilt_entries"] == 1100
    assert result["ABORTED"] is False


def test_a_corrupted_map_aborts_instead_of_reporting(monkeypatch, tmp_path):
    """The canary must be able to FAIL. A check that cannot fail is not a check -- the verify_1a
    lesson, and the reason diag_target_shift's fake result was caught."""
    man = pd.DataFrame({"cell_id": ["a", "b"], "cell_line": ["HFF", "HFF"],
                        "shard_id": ["s0", "s0"], "row_idx": [0, 1],
                        "day_code": [0.0, 1.0], "day": ["D0", "later"]})
    monkeypatch.setattr(s12, "load_fold", lambda fold=None: man)
    fold = tmp_path / "fake"
    (fold / "splits").mkdir(parents=True)
    (fold / "splits" / "holdout.json").write_text(
        json.dumps({"regime": "holdout", "map": {"a": "train", "b": "WRONG"}}), encoding="utf-8")
    monkeypatch.setattr(s12, "ROOT", tmp_path)
    r = s12.run("fake")
    assert r["ABORTED"] is True
    assert r["canary"]["identical"] is False
    assert "first_mismatches" in r["canary"]


# ---- the measurement ------------------------------------------------------------------------ #
@needs_fold
def test_the_recorded_counts_reproduce(result):
    assert result["n_cells"] == 42600
    assert result["n_distinct_ids_old"] == 1100
    assert result["n_distinct_ids_new"] == 42600
    assert result["d0_cells"] == 4988
    assert result["d0_distinct_ids_old"] == 117


@needs_fold
def test_the_old_calib_val_d0_gap_reproduces_the_stage_12_record(result):
    """Stage 12 recorded calib 9.0% vs val 13.3% from the built fold. Re-derived here."""
    old = result["composition_old"]
    assert old["calib"]["d0_share"] == pytest.approx(0.090, abs=0.001)
    assert old["val"]["d0_share"] == pytest.approx(0.133, abs=0.001)
    assert old["train"]["d0_share"] == pytest.approx(0.119, abs=0.001)


@needs_fold
def test_the_fix_pulls_every_split_toward_the_population_d0_rate(result):
    """THE finding. Under the colliding key the three splits scatter over 4.3 points because the
    assignment had ~117 draws to work with; under the fixed key they converge on the population
    rate, because it has 42,600."""
    pop = result["d0_cells"] / result["n_cells"]
    old = result["composition_old"]
    new = result["composition_new"]
    spread_old = max(s["d0_share"] for s in old.values()) - min(s["d0_share"]
                                                                for s in old.values())
    spread_new = max(s["d0_share"] for s in new.values()) - min(s["d0_share"]
                                                                for s in new.values())
    assert spread_old == pytest.approx(0.043, abs=0.002)
    assert spread_new < 0.006
    assert spread_new < spread_old / 5
    for sp in ("train", "val", "calib"):
        assert abs(new[sp]["d0_share"] - pop) < abs(old[sp]["d0_share"] - pop) or \
            abs(old[sp]["d0_share"] - pop) < 0.002


@needs_fold
def test_a_third_of_all_cells_change_split(result):
    """Not a rounding effect: the split map is materially different, which is why the
    model-metric consequence is worth a rebuild rather than being assumed negligible."""
    frac = result["cells_that_change_split"] / result["n_cells"]
    assert 0.30 < frac < 0.40


@needs_fold
def test_calib_is_the_split_that_gains_d0_and_it_is_the_one_that_matters(result):
    """Conformal intervals are computed on calib, and D0 is the control anchor -- which is why
    calib being the most depleted split was the load-bearing part of the Stage 12 harm."""
    d_calib = (result["composition_new"]["calib"]["d0_share"]
               - result["composition_old"]["calib"]["d0_share"])
    assert d_calib > 0.02


# ---- contract ------------------------------------------------------------------------------- #
def test_the_script_writes_only_its_results_file():
    assert SRC.read_text(encoding="utf-8").count(".write_text(") == 1


def test_the_build_parameters_match_the_runner_that_made_the_folds():
    """seed=0 and (0.8, 0.1, 0.1, 0.0) come from local_runners/build_c7_folds.py. If that runner
    changes, this reconstruction silently stops describing the folds on disk."""
    runner = (ROOT / "local_runners" / "build_c7_folds.py").read_text(encoding="utf-8")
    assert "split_fracs=(0.8, 0.1, 0.1, 0.0)" in runner
    assert "seed=0" in runner
    assert s12.SEED == 0 and s12.FRACS == (0.8, 0.1, 0.1, 0.0)


@needs_fold
def test_results_are_json_serialisable(result):
    json.dumps(result)


@needs_fold
def test_d0_is_identified_by_the_smallest_day_code_not_a_hard_coded_float(result):
    src = SRC.read_text(encoding="utf-8")
    assert 'lo = man["day_code"].min()' in src
    assert np.isfinite(result["d0_cells"])
