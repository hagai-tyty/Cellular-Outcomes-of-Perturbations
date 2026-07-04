"""TF/OSKM perturbation encoder: factor-set encoding, the vocabulary embedding,
model dispatch, and the TF data path."""

from __future__ import annotations

import numpy as np
import torch

from cellfate.common import constants as C
from cellfate.data.perturbation import encode_factor_sets, factors_of
from cellfate.models.encoders import TFEncoder
from cellfate.models.network import CellFateNet


def test_factor_sets_share_columns():
    assert factors_of("OSKM") == ("POU5F1", "SOX2", "KLF4", "MYC")
    assert factors_of("OSK") == ("POU5F1", "SOX2", "KLF4")
    assert factors_of("control") == ()
    v = encode_factor_sets(["OSKM", "OSK", "control"])
    assert v.shape == (3, len(C.TF_VOCAB))
    assert (v[0, :3] == 1).all() and (v[1, :3] == 1).all()   # OSKM & OSK share O/S/K
    assert v[0, 3] == 1 and v[1, 3] == 0                      # MYC in OSKM, not OSK
    assert (v[2] == 0).all()                                  # control -> no factors


def test_encode_factor_sets_scales_by_dose():
    v = encode_factor_sets(["OSKM"], [2.5])
    assert v[0, 0] == 2.5 and v[0, 3] == 2.5


def test_tf_encoder_forward_shape():
    enc = TFEncoder(len(C.TF_VOCAB), n_dt=2, d_u=16, p_drop=0.1)
    fv = torch.from_numpy(encode_factor_sets(["OSKM", "OSK"])).float()
    out = enc(fv, torch.zeros(2, 2))
    assert out.shape == (2, 16)


def test_cellfatenet_tf_uses_tf_encoder_and_roundtrips(tmp_path):
    net = CellFateNet(g=20, pert_kind="tf", d_cell=8, d_u=8, latent_dim=8)
    assert type(net.pert).__name__ == "TFEncoder"
    assert net.arch["pert_kind"] == "tf" and net.arch["n_pert"] == len(C.TF_VOCAB)
    x = torch.randn(3, 20)
    u = torch.from_numpy(encode_factor_sets(["OSKM", "OSK", "control"])).float()
    logits, age, z = net(x, u, torch.zeros(3, 2))
    assert logits.shape == (3, 3)
    p = tmp_path / "m.pt"
    net.save_member(p)
    net2 = CellFateNet.load_member(p)                    # arch round-trips the TF encoder
    assert type(net2.pert).__name__ == "TFEncoder" and net2.arch["pert_kind"] == "tf"


def test_chem_model_is_still_the_default():
    net = CellFateNet(g=20)
    assert type(net.pert).__name__ == "ChemEncoder" and net.arch["pert_kind"] == "chem"


def test_tf_dataset_build_writes_tf_shards(tmp_path):
    from cellfate.common import io
    from cellfate.common.io import ArtifactPaths
    from cellfate.data import DataConfig, QCConfig
    from cellfate.data import run as build_run
    from cellfate.data.sources import DataSource, ReprogrammingSource

    genes = [f"G{i}" for i in range(40)]

    class Mem(DataSource):
        name = "reprogramming"

        def plan(self):
            return [{"id": "reprogramming:L1", "cell_line": "L1"},
                    {"id": "reprogramming:L2", "cell_line": "L2"}]

        def fetch(self, chunk):
            r = np.random.default_rng(hash(chunk["id"]) % 100)
            c = r.poisson(20, (40, 40)).astype(np.float32)
            pert = ["control"] * 20 + ["OSKM"] * 20
            time_h = [0.0] * 20 + [312.0] * 20
            return ReprogrammingSource.build_chunk(chunk["id"], c, genes, chunk["cell_line"],
                                                   pert, time_h)

    root = str(tmp_path / "tf")
    build_run(DataConfig(out=root, gene_panel=root + "/panel.json", n_genes=20, clock="random",
                         modality="tf", qc=QCConfig(min_genes=1, max_mito_frac=1.0),
                         split_fracs=(0.6, 0.2, 0.1, 0.1), split_regimes=("cell_line",),
                         primary_regime="cell_line", seed=0),
              sources=[Mem()])
    arr = io.shard_to_numpy(io.read_shard(sorted(ArtifactPaths.of(root).shards_dir.glob("*.parquet"))[0]))
    assert arr["u_modality"][0] == "tf"
    assert arr["u_tf_emb"].shape[1] == len(C.TF_VOCAB)
    assert arr["u_chem_fp"] is None
