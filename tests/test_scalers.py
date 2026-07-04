import numpy as np
import pytest

from cellfate.common.panel import GenePanel
from cellfate.common.scalers import Scalers


def test_fit_transform_standardises():
    rng = np.random.default_rng(0)
    x = rng.normal(5.0, 2.0, size=(1000, 4)).astype(np.float32)
    dt = rng.normal(0.0, 1.0, size=(1000, 2)).astype(np.float32)
    panel = GenePanel([f"G{i}" for i in range(4)])
    sc = Scalers.fit(x, dt, panel, proliferation_coef=(0.5, -1.0))

    x_scaled = sc.transform_x(x)
    assert np.allclose(x_scaled.mean(0), 0.0, atol=1e-2)
    assert np.allclose(x_scaled.std(0), 1.0, atol=1e-2)
    assert sc.proliferation_coef == (0.5, -1.0)
    assert sc.params.gene_panel_hash == panel.hash()


def test_transform_rejects_wrong_width():
    panel = GenePanel([f"G{i}" for i in range(4)])
    sc = Scalers.fit(np.zeros((10, 4), np.float32), np.zeros((10, 2), np.float32), panel)
    with pytest.raises(ValueError):
        sc.transform_x(np.zeros((3, 5), np.float32))


def test_save_load_roundtrip(tmp_path):
    panel = GenePanel([f"G{i}" for i in range(4)])
    sc = Scalers.fit(np.random.rand(50, 4).astype(np.float32),
                     np.random.rand(50, 2).astype(np.float32), panel)
    f = tmp_path / "scalers.json"
    sc.save(f)
    sc2 = Scalers.load(f)
    x = np.random.rand(2, 4).astype(np.float32)
    assert np.allclose(sc.transform_x(x), sc2.transform_x(x))
