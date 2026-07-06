import numpy as np
import pytest

from cellfate.common import io
from cellfate.common.errors import BundleError, GenePanelMismatch, ShardIOError
from cellfate.common.io import ArtifactPaths
from cellfate.common.panel import GenePanel
from cellfate.common.schemas import (
    BundleMeta,
    ConformalParams,
    ManifestRow,
    ResParams,
    TemperatureParams,
)
from conftest import G_TEST, make_sample


def test_sanitize_id():
    assert io.sanitize_id("tahoe:HEPG2:rapamycin") == "tahoe_HEPG2_rapamycin"
    assert io.sanitize_id("a//b\\c") == "a_b_c"


def test_shard_roundtrip(tmp_path):
    paths = ArtifactPaths.of(tmp_path)
    samples = [make_sample(cell_id=f"c{i}", y_age=float(-i)) for i in range(5)]
    shard = paths.shard_file("chunk0")
    io.write_shard(shard, samples)

    table = io.read_shard(shard)
    assert table.num_rows == 5

    arr = io.shard_to_numpy(table)
    assert arr["X"].shape == (5, G_TEST)
    assert arr["dose_time"].shape == (5, 2)
    assert arr["y_cls"].shape == (5, 3)
    assert arr["u_chem_fp"].shape == (5, 2048)
    assert arr["age_mask"].all()
    assert np.allclose(arr["y_age"], [0, -1, -2, -3, -4])


def test_shard_masked_age_becomes_nan(tmp_path):
    paths = ArtifactPaths.of(tmp_path)
    s = make_sample(cell_id="m", age_mask=False, y_age=None)
    io.write_shard(paths.shard_file("c"), [s])
    arr = io.shard_to_numpy(io.read_shard(paths.shard_file("c")))
    assert np.isnan(arr["y_age"]).all()


def test_empty_shard_refused(tmp_path):
    paths = ArtifactPaths.of(tmp_path)
    with pytest.raises(ShardIOError):
        io.write_shard(paths.shard_file("empty"), [])


def test_manifest_fragments_consolidate(tmp_path):
    paths = ArtifactPaths.of(tmp_path)
    for cid in ("a", "b"):
        rows = [ManifestRow.from_sample(make_sample(cell_id=f"{cid}{i}"), shard_id=cid, row_idx=i)
                for i in range(3)]
        io.write_manifest_part(paths, cid, rows)
    table = io.load_manifest(paths)             # consolidates parts -> manifest.parquet
    assert table.num_rows == 6
    assert paths.manifest_file.exists()
    rows = io.manifest_rows(table)
    assert {r.shard_id for r in rows} == {"a", "b"}


def test_splits_roundtrip(tmp_path):
    paths = ArtifactPaths.of(tmp_path)
    mapping = {"c0": "train", "c1": "val", "c2": "test", "c3": "calib"}
    io.write_splits(paths, "scaffold", mapping)
    assert io.load_splits(paths, "scaffold") == mapping


def test_gene_panel_assertion():
    panel = GenePanel(["A", "B", "C"])
    io.assert_gene_panel(panel.hash(), panel)   # ok
    with pytest.raises(GenePanelMismatch):
        io.assert_gene_panel("deadbeef", panel)


def test_config_hash_is_stable_and_order_independent():
    h1 = io.hash_config({"a": 1, "b": {"c": 2}})
    h2 = io.hash_config({"b": {"c": 2}, "a": 1})
    assert h1 == h2 and len(h1) == 12


def test_bundle_meta_roundtrip(tmp_path):
    paths = ArtifactPaths.of(tmp_path)
    meta = BundleMeta(n_members=5, gene_panel_hash="abc123", conformal_levels=[0.9])
    io.save_bundle_meta(paths, meta)
    assert io.load_bundle_meta(paths) == meta


def test_bundle_param_files_roundtrip(tmp_path):
    paths = ArtifactPaths.of(tmp_path)
    io.save_temperature(paths, TemperatureParams(temperature=1.7))
    io.save_conformal(paths, ConformalParams(levels=[0.9], q={"0.9": 1.83}))
    io.save_res_params(paths, ResParams(k=2.0))
    assert io.load_temperature(paths).temperature == 1.7
    assert io.load_conformal(paths).q["0.9"] == 1.83
    assert io.load_res_params(paths).k == 2.0


def test_assert_bundle_complete_detects_missing(tmp_path):
    paths = ArtifactPaths.of(tmp_path)
    with pytest.raises(BundleError):
        io.assert_bundle_complete(paths)        # nothing written yet


def test_atomic_write_is_visible(tmp_path):
    p = tmp_path / "sub" / "x.json"
    io.write_json(p, {"k": 1})
    assert io.read_json(p) == {"k": 1}
