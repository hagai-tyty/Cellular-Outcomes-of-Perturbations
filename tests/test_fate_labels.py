"""Fate-label design: markers held out of the panel (anti-circularity) and the
control-relative / distributional threshold."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cellfate.common import constants as C
from cellfate.data.labels import fate_labels
from cellfate.data.normalize import fit_gene_panel


def test_label_markers_held_out_of_panel():
    # LABEL_HOLDOUT genes present in the data must NOT enter the model panel,
    # even when they are the highest-variance genes (they'd be picked otherwise).
    holdout = list(C.LABEL_HOLDOUT)
    genes = holdout + [f"G{i}" for i in range(60)]
    rng = np.random.default_rng(0)
    expr = rng.normal(size=(80, len(genes)))
    expr[:, : len(holdout)] *= 25.0                      # make markers high-variance
    panel = fit_gene_panel(expr, genes, n_top=40, must_exclude=C.LABEL_HOLDOUT)
    assert set(panel.genes).isdisjoint(C.LABEL_HOLDOUT)  # zero leakage
    assert len(panel.genes) == 40


def test_fate_labels_are_valid_distribution():
    genes = list(C.LABEL_HOLDOUT) + [f"G{i}" for i in range(10)]
    rng = np.random.default_rng(1)
    expr = rng.normal(0.5, 0.1, size=(30, len(genes)))
    obs = pd.DataFrame({"cell_line": ["L"] * 30, "is_control": [True] * 10 + [False] * 20})
    y = fate_labels(expr, genes, obs, tau=1.0)
    assert y.shape == (30, 3)
    assert np.allclose(y.sum(axis=1), 1.0) and (y >= 0).all()


def test_fate_labels_control_relative_directions():
    genes = list(C.LABEL_HOLDOUT) + [f"G{i}" for i in range(5)]
    gidx = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(2)
    n = 40
    expr = rng.normal(0.0, 0.05, size=(n, len(genes)))         # tight control baseline
    obs = pd.DataFrame({"cell_line": ["L"] * n, "is_control": [True] * 20 + [False] * 20})
    # cell 20: gains pluripotency AND loses somatic identity -> identity loss
    for g in C.DEFAULT_SIGNATURES["loss"]:
        expr[20, gidx[g]] += 3.0
    for g in C.DEFAULT_SIGNATURES["safe"]:
        expr[20, gidx[g]] -= 3.0
    # cell 21: elevated apoptosis -> death
    for g in C.DEFAULT_SIGNATURES["death"]:
        expr[21, gidx[g]] += 4.0
    y = fate_labels(expr, genes, obs, tau=1.0)               # columns = (safe, loss, death)
    assert y[:20, 0].mean() > 0.5      # controls are mostly "safe"
    assert y[20].argmax() == 1         # loss
    assert y[21].argmax() == 2         # death
