"""Unit tests for the cellfate.data scientific functions (Document 2)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cellfate.common import constants as C
from cellfate.common.panel import GenePanel
from cellfate.common.schemas import ManifestRow
from cellfate.data import (
    LinearClock,
    QCConfig,
    RawChunk,
    apply_qc,
    cell_cycle_score,
    deconfound_age,
    delta_age,
    encode_dose_time,
    encode_fingerprints,
    fit_deconfounder,
    fit_gene_panel,
    hashed_fingerprint,
    make_splits,
    morgan_fingerprint,
    normalize_counts,
    recenter_on_controls,
    resolve_scaffolds,
    score_one,
    signature_scores,
    soft_labels,
    to_panel_matrix,
)
from cellfate.data.sources import OBS_COLUMNS
from cellfate.data.splits import CONTROL_SCAFFOLD


# --------------------------------------------------------------------------- #
# labels                                                                      #
# --------------------------------------------------------------------------- #
def test_soft_labels_sum_to_one_and_nonneg():
    rng = np.random.default_rng(0)
    p = soft_labels(rng.normal(size=(20, 3)), tau=1.0)
    assert p.shape == (20, 3)
    assert np.allclose(p.sum(axis=1), 1.0)
    assert (p >= 0).all()


def test_soft_labels_temperature_softens_but_keeps_argmax():
    s = np.array([[3.0, 0.0, 0.0]])
    sharp = soft_labels(s, tau=0.5)
    soft = soft_labels(s, tau=5.0)
    assert sharp.argmax() == soft.argmax() == 0
    assert soft.std() < sharp.std()  # higher tau -> closer to uniform


def test_soft_labels_rejects_nonpositive_tau():
    with pytest.raises(ValueError):
        soft_labels(np.zeros((2, 3)), tau=0.0)


# --------------------------------------------------------------------------- #
# signatures                                                                  #
# --------------------------------------------------------------------------- #
def test_signature_scores_order_and_planted_signal():
    genes: list[str] = []
    for gs in C.DEFAULT_SIGNATURES.values():
        genes += list(gs)
    seen: set[str] = set()
    genes = [g for g in genes if not (g in seen or seen.add(g))]
    gidx = {g: i for i, g in enumerate(genes)}

    rng = np.random.default_rng(1)
    n_per = 12
    expr = rng.normal(0.5, 0.1, size=(3 * n_per, len(genes)))
    for ci, cls in enumerate(C.CLASSES):  # plant each class's markers high
        rows = slice(ci * n_per, (ci + 1) * n_per)
        for g in C.DEFAULT_SIGNATURES[cls]:
            expr[rows, gidx[g]] += 3.0

    sig = signature_scores(expr, genes)
    assert sig.shape == (3 * n_per, 3)
    for ci in range(3):  # the matching column dominates for each planted block
        rows = slice(ci * n_per, (ci + 1) * n_per)
        block_mean = sig[rows].mean(axis=0)
        assert block_mean.argmax() == ci


def test_score_one_absent_geneset_is_zero():
    out = score_one(np.ones((4, 3)), ["A", "B", "C"], ("X", "Y"))
    assert np.allclose(out, 0.0)


# --------------------------------------------------------------------------- #
# proliferation / deconfounding                                               #
# --------------------------------------------------------------------------- #
def test_deconfounder_recovers_slope_and_removes_correlation():
    rng = np.random.default_rng(2)
    cc = rng.uniform(0, 1, size=600)
    delta = 2.5 * cc + rng.normal(0, 0.3, size=600)
    a, b = fit_deconfounder(delta, cc)
    assert abs(a - 2.5) < 0.3
    dec = deconfound_age(delta, cc, (a, b))
    assert abs(np.corrcoef(dec, cc)[0, 1]) < 0.1  # cell-cycle signal removed


def test_cell_cycle_score_tracks_latent_proliferation():
    cc_genes = list(C.S_GENES) + list(C.G2M_GENES)
    genes = cc_genes + [f"FILL{i}" for i in range(40)]
    rng = np.random.default_rng(3)
    latent = rng.uniform(0, 1, size=200)
    expr = rng.normal(0, 0.1, size=(200, len(genes)))
    for j in range(len(cc_genes)):
        expr[:, j] += 5.0 * latent
    score = cell_cycle_score(expr, genes)
    assert np.corrcoef(score, latent)[0, 1] > 0.7


# --------------------------------------------------------------------------- #
# normalise + panel                                                           #
# --------------------------------------------------------------------------- #
def test_normalize_counts_is_cp10k_log1p():
    counts = np.array([[1, 1, 2], [0, 0, 0], [10, 0, 0]], dtype=float)
    norm = normalize_counts(counts)
    assert norm.dtype == np.float32
    assert np.isclose(np.expm1(norm[0]).sum(), 1e4, rtol=1e-3)
    assert np.allclose(norm[1], 0.0)  # empty cell stays zero (no divide-by-zero)


def test_fit_gene_panel_size_and_must_include():
    genes = [f"G{i}" for i in range(50)]
    rng = np.random.default_rng(4)
    expr = rng.normal(size=(100, 50))
    expr[:, 10] *= 20  # a couple of high-variance genes
    expr[:, 20] *= 15
    panel = fit_gene_panel(expr, genes, n_top=8, must_include=("G49", "G48"))
    assert len(panel) == 8
    assert "G49" in panel.genes and "G48" in panel.genes
    assert len(set(panel.genes)) == len(panel.genes)  # no duplicates


def test_to_panel_matrix_reorders_and_zero_fills_missing():
    genes = ["A", "B", "C"]
    norm = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    panel = GenePanel(["C", "A", "Z"])
    x = to_panel_matrix(norm, genes, panel)
    assert x.shape == (2, 3)
    assert np.allclose(x[:, 0], [3, 6])  # C
    assert np.allclose(x[:, 1], [1, 4])  # A
    assert np.allclose(x[:, 2], [0, 0])  # Z absent -> 0


# --------------------------------------------------------------------------- #
# perturbation encoding                                                       #
# --------------------------------------------------------------------------- #
def test_hashed_fingerprint_deterministic_binary_and_sized():
    fp1 = hashed_fingerprint("CCO")
    fp2 = hashed_fingerprint("CCO")
    assert fp1.shape == (C.N_FINGERPRINT_BITS,)
    assert np.array_equal(fp1, fp2)
    assert set(np.unique(fp1)).issubset({0, 1})
    assert not np.array_equal(fp1, hashed_fingerprint("c1ccccc1"))
    assert hashed_fingerprint("").sum() == 0  # control -> empty


def test_morgan_fingerprint_runs_with_or_without_rdkit():
    fp = morgan_fingerprint("CCO")
    assert fp.shape == (C.N_FINGERPRINT_BITS,)
    assert set(np.unique(fp)).issubset({0, 1})
    assert fp.sum() > 0


def test_encode_fingerprints_shape_and_sharing():
    fps = encode_fingerprints(["CCO", "CCO", "c1ccccc1"])
    assert fps.shape == (3, C.N_FINGERPRINT_BITS)
    assert fps.dtype == np.uint8
    assert np.array_equal(fps[0], fps[1])


def test_encode_dose_time_values_and_floor():
    dt = encode_dose_time([1.0, 0.0, 100.0], [24.0, 0.0, 48.0])
    assert dt.shape == (3, 2)
    assert np.isfinite(dt).all()
    assert np.isclose(dt[0, 0], 0.0)        # log10(1)
    assert np.isclose(dt[2, 0], 2.0)        # log10(100)
    assert dt[1, 0] < -3                    # floored control dose
    assert np.isclose(dt[0, 1], np.log(24.0))


def test_resolve_scaffolds_nonempty_with_preference():
    out = resolve_scaffolds(
        smiles=["", "SYN-x", "SYN-y"],
        pert_ids=["control", "cx", "cy"],
        provided=["CONTROL", "SCAF0", None],
    )
    assert out[0] == "CONTROL"
    assert out[1] == "SCAF0"
    assert out[2] == "cy"  # no provided + no rdkit -> falls back to pert_id
    assert all(s for s in out)


# --------------------------------------------------------------------------- #
# QC                                                                          #
# --------------------------------------------------------------------------- #
def _obs(n: int, **cols) -> pd.DataFrame:
    base = {c: ["x"] * n for c in OBS_COLUMNS}
    base.update(cols)
    return pd.DataFrame(base)


def test_apply_qc_drops_low_gene_and_high_mito_cells():
    genes = ["g0", "g1", "g2", "MT-CO1"]
    counts = np.array([
        [5, 5, 5, 1],    # good: 4 genes, mito 1/16
        [5, 0, 0, 0],    # bad: only 1 gene expressed
        [1, 1, 0, 20],   # bad: mito 20/22 high
    ], dtype=float)
    obs = _obs(3, cell_id=["a", "b", "c"], is_control=[True, False, False])
    chunk = RawChunk("chk", "synth", counts, genes, obs)
    kept = apply_qc(chunk, QCConfig(min_genes=2, max_mito_frac=0.5))
    assert len(kept.obs) == 1
    assert kept.obs["cell_id"].tolist() == ["a"]
    assert kept.counts.shape == (1, 4)


# --------------------------------------------------------------------------- #
# aging clock + delta age                                                     #
# --------------------------------------------------------------------------- #
def test_delta_age_centres_on_controls_for_normal_source():
    panel = GenePanel([f"G{i}" for i in range(5)])
    genes = list(panel.genes)
    clock = LinearClock.random(panel, seed=0)
    rng = np.random.default_rng(5)
    x = rng.normal(size=(30, 5))
    obs = pd.DataFrame({
        "cell_line": ["L0"] * 15 + ["L1"] * 15,
        "is_control": ([True] * 5 + [False] * 10) * 2,
    })
    d, mask = delta_age(clock, x, genes, obs, source="synth")
    assert mask.all()  # non-cancer source -> all ages valid
    line0_ctrl = d[:5]
    assert abs(line0_ctrl.mean()) < 1e-6  # ΔAge centred on each line's controls


def test_delta_age_masks_cancer_sources():
    panel = GenePanel([f"G{i}" for i in range(5)])
    clock = LinearClock.random(panel, seed=0)
    obs = pd.DataFrame({"cell_line": ["L0"] * 6, "is_control": [True] * 3 + [False] * 3})
    _, mask = delta_age(clock, np.zeros((6, 5)), list(panel.genes), obs, source="tahoe")
    assert not mask.any()  # tahoe is in CANCER_SOURCES


def test_deconfound_then_recenter_preserves_control_zero_point():
    # Regression: cell-cycle deconfounding subtracts the regression intercept and so
    # re-centres the whole population to mean 0, which shifts the controls off zero.
    # recenter_on_controls must restore the control-relative zero-point (controls ~ 0)
    # WITHOUT reintroducing the cell-cycle slope that deconfounding removed.
    rng = np.random.default_rng(3)
    n = 600
    line = np.where(rng.random(n) < 0.5, "L0", "L1")
    obs = pd.DataFrame({"cell_line": line, "is_control": rng.random(n) < 0.25})
    is_ctrl = obs["is_control"].to_numpy()
    cc = rng.normal(0, 1, n)
    d = 0.6 * cc + rng.normal(0, 0.3, n)              # one shared cell-cycle slope
    d[~is_ctrl] -= rng.uniform(0, 1.5, (~is_ctrl).sum())  # perturbed arm skews younger
    for L in ("L0", "L1"):
        m = line == L
        d[m] -= d[m & is_ctrl].mean()               # control-relative to start
    assert abs(d[is_ctrl].mean()) < 1e-9

    coef = fit_deconfounder(d, cc)
    d_dc = deconfound_age(d, cc, coef)
    assert abs(d_dc[is_ctrl].mean()) > 0.1          # deconfounding alone shifts controls off zero

    d_fixed = recenter_on_controls(d_dc, obs)
    assert abs(d_fixed[is_ctrl].mean()) < 1e-9      # control zero-point restored (overall)
    for L in ("L0", "L1"):                          # ...and per line
        m = line == L
        assert abs(d_fixed[m & is_ctrl].mean()) < 1e-9
    # the cell-cycle confound stays removed (slope 0.6 -> ~0), not reintroduced by recentering
    assert abs(np.polyfit(cc, d_fixed, 1)[0]) < 0.05


# --------------------------------------------------------------------------- #
# splits                                                                      #
# --------------------------------------------------------------------------- #
def _rows() -> list[ManifestRow]:
    rows: list[ManifestRow] = []
    k = 0
    # 6 real scaffolds across 3 cell lines
    for s in range(6):
        for line in range(3):
            for _ in range(4):
                rows.append(ManifestRow(
                    cell_id=f"c{k}", cell_line=f"L{line}", pert_id=f"p{s}",
                    scaffold_id=f"S{s}", source="synth", age_mask=True,
                    shard_id="sh", row_idx=k))
                k += 1
    # controls
    for line in range(3):
        for _ in range(6):
            rows.append(ManifestRow(
                cell_id=f"c{k}", cell_line=f"L{line}", pert_id="control",
                scaffold_id=CONTROL_SCAFFOLD, source="synth", age_mask=True,
                shard_id="sh", row_idx=k))
            k += 1
    return rows


def test_scaffold_split_is_group_disjoint_with_controls_in_train():
    rows = _rows()
    splits = make_splits(rows, (0.5, 0.2, 0.2, 0.1), ("scaffold",), seed=0)["scaffold"]
    assert set(splits) == {r.cell_id for r in rows}

    by_scaffold: dict[str, set[str]] = {}
    cid_to_scaf = {r.cell_id: (r.scaffold_id or r.pert_id) for r in rows}
    for cid, sp in splits.items():
        by_scaffold.setdefault(cid_to_scaf[cid], set()).add(sp)
    # each scaffold lands in exactly one split (no leakage across train/test)
    for scaf, sps in by_scaffold.items():
        assert len(sps) == 1, f"scaffold {scaf} leaked across splits {sps}"
    # controls pinned to train
    assert by_scaffold[CONTROL_SCAFFOLD] == {"train"}
    # the split actually uses more than one bucket
    assert len(set(splits.values())) > 1


def test_cell_line_split_is_group_disjoint():
    rows = _rows()
    splits = make_splits(rows, (0.34, 0.33, 0.33, 0.0), ("cell_line",), seed=0)["cell_line"]
    by_line: dict[str, set[str]] = {}
    cid_to_line = {r.cell_id: r.cell_line for r in rows}
    for cid, sp in splits.items():
        by_line.setdefault(cid_to_line[cid], set()).add(sp)
    for line, sps in by_line.items():
        assert len(sps) == 1, f"cell line {line} leaked across splits {sps}"


def test_random_split_fills_all_four_splits():
    # random regime: cell-level (for single-donor data); every split populated
    rows = _rows()
    splits = make_splits(rows, (0.7, 0.1, 0.1, 0.1), ("random",), seed=0)["random"]
    assert len(splits) == len(rows)
    assert set(splits.values()) >= {"train", "val", "calib", "test"}

