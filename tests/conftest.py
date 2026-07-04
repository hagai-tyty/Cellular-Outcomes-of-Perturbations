"""Shared pytest fixtures for the foundation tests."""
from __future__ import annotations

import pytest

from cellfate.common import constants as C
from cellfate.common.panel import GenePanel
from cellfate.common.schemas import Modality, Sample

G_TEST = 8  # tiny gene count for fast tests


@pytest.fixture
def panel() -> GenePanel:
    return GenePanel([f"GENE{i}" for i in range(G_TEST)])


def make_sample(**overrides) -> Sample:
    """Construct a valid chemical Sample; override any field via kwargs."""
    base = dict(
        cell_id="tahoe:AAA",
        X=[0.0] * G_TEST,
        u_modality=Modality.CHEM,
        u_chem_fp=[0] * C.N_FINGERPRINT_BITS,
        u_gene_emb=None,
        u_tf_emb=None,
        dose_time=[1.0, 1.0],
        y_cls=[0.8, 0.1, 0.1],
        y_age=-3.0,
        age_mask=True,
        sig_scores=[1.0, 0.0, 0.0],
        cell_line="HEPG2",
        pert_id="rapamycin",
        scaffold_id="c1ccccc1",
        source="tahoe",
    )
    base.update(overrides)
    return Sample(**base)


@pytest.fixture
def sample_factory():
    return make_sample
