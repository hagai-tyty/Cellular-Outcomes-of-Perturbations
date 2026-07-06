"""Document 5 (cellfate.evaluation): metric correctness, baseline sanity, gate-logic
boundary flips, and an end-to-end evaluation against a real bundle built by Docs 1-4.

Discipline note mirrored from the spec: the acceptance gates are designed for *real*
data. On trivially-separable synthetic data the linear baselines are near-perfect, so
the model does not beat them and several gates fail by design -- the integration tests
therefore assert that the *machinery* is correct (every metric computable, reports and
well-formed gates written, a multi-class regime yields finite metrics), not that the
synthetic model passes every gate.
"""

from __future__ import annotations

import json
import warnings

import numpy as np
import pytest

from cellfate.common.constants import Split
from cellfate.common.seeding import set_global_seed
from cellfate.evaluation import (
    BASELINE_NAMES,
    Estimator,
    EvalConfig,
    KNNFingerprint,
    MeanBaseline,
    PredictControl,
    check_gates,
    des_pds,
    evaluate,
    make_baselines,
    ranking_metrics,
    validate_against_methylation,
    validate_oskm_holdout,
)
from cellfate.evaluation.data import SplitData
from cellfate.evaluation.metrics import (
    brier,
    coverage,
    ece,
    mean_finite,
    per_class_auroc,
    per_class_prauc,
    precision_at_k,
    regression_metrics,
)


# ===================================================================== #
# Pure-unit metric tests (no bundle needed)                             #
# ===================================================================== #
def test_ece_zero_for_confident_correct():
    rng = np.random.default_rng(0)
    y = rng.integers(0, 3, 500)
    p = np.eye(3)[y]  # probability 1 on the true class -> conf==acc==1 in one bin
    assert ece(y, p) == pytest.approx(0.0, abs=1e-12)


def test_ece_positive_for_overconfident_wrong():
    y = np.zeros(100, dtype=int)
    p = np.tile([0.05, 0.95, 0.0], (100, 1))  # 95% confident, always wrong
    assert ece(y, p) > 0.9


def test_coverage_returns_known_fraction():
    y = np.arange(10, dtype=float)
    lo = np.full(10, -1.0)
    hi = np.full(10, 6.5)  # covers y = 0..6 -> 7 of 10
    assert coverage(y, lo, hi, np.ones(10, bool)) == pytest.approx(0.7)


def test_coverage_respects_mask():
    y = np.array([0.0, 100.0, 0.0])
    lo, hi = np.full(3, -1.0), np.full(3, 1.0)
    mask = np.array([True, False, True])  # the out-of-interval point is masked out
    assert coverage(y, lo, hi, mask) == pytest.approx(1.0)


def test_precision_at_k_perfect_and_reversed():
    res = np.array([5.0, 4.0, 3.0, 2.0, 1.0, 0.0])
    quality = np.array([5.0, 4.0, 3.0, 2.0, 1.0, 0.0])
    assert precision_at_k(res, quality, 3) == pytest.approx(1.0)
    assert precision_at_k(res, quality[::-1], 3) == pytest.approx(0.0)


def test_ranking_spearman_aligned_and_reversed():
    res = np.arange(10, dtype=float)
    shift = -np.arange(10, dtype=float)  # quality = -shift = arange -> aligned with res
    assert ranking_metrics(res, shift)["spearman"] == pytest.approx(1.0)
    assert ranking_metrics(res, -shift)["spearman"] == pytest.approx(-1.0)


def test_ranking_nan_when_degenerate():
    res = np.zeros(10)  # constant RES -> undefined ranking
    assert np.isnan(ranking_metrics(res, np.arange(10.0))["spearman"])


def test_brier_known_value():
    y1h = np.array([[1.0, 0.0, 0.0]])
    p = np.array([[0.7, 0.2, 0.1]])
    # (0.3^2 + 0.2^2 + 0.1^2) = 0.14
    assert brier(y1h, p) == pytest.approx(0.14)


def test_regression_metrics_perfect_and_known():
    y = np.array([1.0, 2.0, 3.0, 4.0])
    m = regression_metrics(y, y, np.ones(4, bool))
    assert m["mae"] == 0 and m["rmse"] == 0 and m["r2"] == pytest.approx(1.0)
    assert m["pearson"] == pytest.approx(1.0) and m["n"] == 4
    yhat = y + 1.0  # constant offset: mae=1, rmse=1, pearson=1, r2=1-(4/5)
    m2 = regression_metrics(y, yhat, np.ones(4, bool))
    assert m2["mae"] == pytest.approx(1.0) and m2["rmse"] == pytest.approx(1.0)
    assert m2["pearson"] == pytest.approx(1.0)


def test_regression_metrics_empty_mask_is_nan():
    m = regression_metrics(np.arange(3.0), np.arange(3.0), np.zeros(3, bool))
    assert m["n"] == 0 and np.isnan(m["mae"]) and np.isnan(m["pearson"])


def test_per_class_auroc_prauc_perfect_separation():
    y = np.array([0, 0, 1, 1, 2, 2])
    p = np.eye(3)[y] * 0.9 + 0.05  # confidently correct
    au = per_class_auroc(y, p)
    pr = per_class_prauc(y, p)
    assert all(au[c] == pytest.approx(1.0) for c in range(3))
    assert all(pr[c] == pytest.approx(1.0) for c in range(3))


def test_per_class_metrics_nan_when_class_absent():
    y = np.array([2, 2, 2, 2])  # only class 2 present
    p = np.tile([0.2, 0.3, 0.5], (4, 1))
    au = per_class_auroc(y, p)
    pr = per_class_prauc(y, p)
    assert np.isnan(au[0]) and np.isnan(au[1]) and np.isnan(au[2])  # 2 is all-positive -> nan
    assert np.isnan(pr[0]) and np.isnan(pr[1])  # absent classes -> nan
    assert pr[2] == pytest.approx(1.0)  # present-everywhere class -> trivially 1.0


def test_des_pds_identical_is_one():
    rng = np.random.default_rng(1)
    ctrl = rng.normal(0, 1, (40, 16))
    state = rng.normal(2, 1, (8, 16))
    out = des_pds(state, state, ctrl)
    assert out["des"] == pytest.approx(1.0) and out["pds"] == pytest.approx(1.0)


def test_des_pds_discrimination_drops_with_noise():
    rng = np.random.default_rng(2)
    ctrl = rng.normal(0, 1, (40, 16))
    true = rng.normal(0, 1, (10, 16))
    pred = rng.normal(0, 1, (10, 16))  # unrelated predictions -> poor discrimination
    assert des_pds(pred, true, ctrl)["pds"] < 0.5


def test_mean_finite_ignores_nan():
    assert mean_finite([1.0, np.nan, 3.0]) == pytest.approx(2.0)
    assert np.isnan(mean_finite([np.nan, np.nan]))


# ===================================================================== #
# External validation (pure functions)                                  #
# ===================================================================== #
def test_methylation_anchor_correlation():
    rng = np.random.default_rng(3)
    meth = rng.normal(0, 3, 200)
    pred = meth + 1.0 + rng.normal(0, 0.1, 200)  # tracks with +1 bias
    out = validate_against_methylation(pred, meth)
    assert out["pearson"] > 0.99 and out["bias"] == pytest.approx(1.0, abs=0.1)


def test_oskm_holdout_detects_rejuvenation_and_decoupling():
    # ΔAge falls monotonically; identity-loss prob stays flat until the last step
    ages = [0.0, -1.0, -2.0, -3.0]
    ploss = [0.02, 0.02, 0.03, 0.40]

    def predict_fn(X, fp, dt):
        i = int(X)  # X encodes the timepoint index in this stub
        return np.array([[1 - ploss[i], ploss[i], 0.0]]), np.array([ages[i]])

    tc = [(i, None, None) for i in range(4)]
    out = validate_oskm_holdout(predict_fn, tc)
    assert out["rejuvenates"] is True
    assert out["age_identity_decoupled"] is True
    assert out["delta_age_trajectory"] == ages


# ===================================================================== #
# Baseline sanity (small constructed SplitData)                         #
# ===================================================================== #
def _toy_split(n=60, seed=0) -> SplitData:
    rng = np.random.default_rng(seed)
    return SplitData(
        X=rng.normal(0, 1, (n, 10)).astype(np.float32),
        fp=(rng.random((n, 2048)) > 0.5).astype(np.float32),
        dose_time=rng.normal(0, 1, (n, 2)).astype(np.float32),
        y_cls=rng.integers(0, 3, n).astype(np.int64),
        y_age=rng.normal(-2, 1, n),
        mask=np.ones(n, bool),
        scaffold_id=np.array(["CONTROL" if i < 15 else f"S{i%4}" for i in range(n)]),
        cell_line=np.array(["L0"] * n),
        cell_id=np.array([f"c{i}" for i in range(n)]),
    )


def test_mean_baseline_reproduces_train_marginals():
    tr = _toy_split()
    est = MeanBaseline().fit(tr)
    p, age = est.predict(tr.X, tr.fp, tr.dose_time)
    freq = np.bincount(tr.y_cls, minlength=3) / tr.n
    assert np.allclose(p[0], freq) and np.allclose(p, p[0])  # constant across cells
    assert np.allclose(age, tr.y_age[tr.mask].mean())


def test_predict_control_zero_dage_and_control_marginals():
    tr = _toy_split()
    est = PredictControl().fit(tr)
    p, age = est.predict(tr.X, tr.fp, tr.dose_time)
    assert np.allclose(age, 0.0)
    ctrl = tr.scaffold_id == "CONTROL"
    exp = np.bincount(tr.y_cls[ctrl], minlength=3) / ctrl.sum()
    assert np.allclose(p[0], exp)


def test_all_baselines_obey_estimator_interface():
    tr = _toy_split()
    for name, est in make_baselines(BASELINE_NAMES).items():
        assert isinstance(est, Estimator), name
        est.fit(tr)
        p, age = est.predict(tr.X, tr.fp, tr.dose_time)
        assert p.shape == (tr.n, 3), name
        assert age.shape == (tr.n,), name
        assert np.allclose(p.sum(1), 1.0, atol=1e-6), name
        assert np.all(p >= -1e-9), name


def test_knn_recovers_neighbour_label_on_separable_fp():
    # two fingerprint clusters, each a distinct class -> kNN recovers the label
    rng = np.random.default_rng(5)
    a = np.zeros((20, 2048))
    a[:, :50] = 1
    b = np.zeros((20, 2048))
    b[:, 50:100] = 1
    fp = np.vstack([a, b]).astype(np.float32)
    y = np.array([0] * 20 + [2] * 20)
    tr = SplitData(X=rng.normal(0, 1, (40, 4)).astype(np.float32), fp=fp,
                   dose_time=np.zeros((40, 2), np.float32), y_cls=y, y_age=np.zeros(40),
                   mask=np.zeros(40, bool), scaffold_id=np.array(["S"] * 40),
                   cell_line=np.array(["L"] * 40), cell_id=np.array([f"c{i}" for i in range(40)]))
    est = KNNFingerprint(k=5).fit(tr)
    p, _ = est.predict(fp, fp, tr.dose_time)
    assert p[:20, 0].mean() > 0.9   # cluster A -> class 0
    assert p[20:, 2].mean() > 0.9   # cluster B -> class 2


# ===================================================================== #
# Gate logic: boundary flips (constructed results dicts)                #
# ===================================================================== #
def _results(model_prauc, model_mae, model_ece, cov, spearman,
             base_prauc=0.2, base_mae=1.0):
    def flat(prauc, mae, ece_=0.5):
        return {"prauc_0": prauc, "prauc_1": prauc, "prauc_2": prauc,
                "reg_mae": mae, "ece": ece_}
    R = {"model": flat(model_prauc, model_mae, model_ece),
         "mean": flat(base_prauc, base_mae),
         "ridge": flat(base_prauc, base_mae),
         "coverage": cov, "ranking": {"spearman": spearman}}
    cfg = EvalConfig(bundle="_", dataset="_", baselines=("mean", "ridge"),
                     max_ece=0.05, level=0.90, cov_tol=0.03, min_spearman=0.3)
    return {"r": R}, cfg


def test_gate_beats_all_baselines_flips():
    res, cfg = _results(0.9, 0.1, 0.01, 0.90, 0.5)
    assert check_gates(res, cfg)["r"]["beats_all_baselines"] is True
    # tie a baseline's PR-AUC to the model -> no strict win -> gate fails
    res2, cfg2 = _results(0.9, 0.1, 0.01, 0.90, 0.5, base_prauc=0.9)
    assert check_gates(res2, cfg2)["r"]["beats_all_baselines"] is False
    # baseline better MAE even with higher PR-AUC -> gate fails
    res3, cfg3 = _results(0.9, 1.0, 0.01, 0.90, 0.5, base_prauc=0.2, base_mae=0.1)
    assert check_gates(res3, cfg3)["r"]["beats_all_baselines"] is False


def test_gate_ece_flips_around_threshold():
    below, cfg = _results(0.9, 0.1, 0.049, 0.90, 0.5)
    above, cfg2 = _results(0.9, 0.1, 0.051, 0.90, 0.5)
    assert check_gates(below, cfg)["r"]["ece_ok"] is True
    assert check_gates(above, cfg2)["r"]["ece_ok"] is False


def test_gate_coverage_flips_around_tolerance():
    inside, cfg = _results(0.9, 0.1, 0.01, 0.88, 0.5)   # |0.88-0.90|=0.02 < 0.03
    outside, cfg2 = _results(0.9, 0.1, 0.01, 0.80, 0.5)  # |0.80-0.90|=0.10 > 0.03
    assert check_gates(inside, cfg)["r"]["coverage_ok"] is True
    assert check_gates(outside, cfg2)["r"]["coverage_ok"] is False


def test_gate_ranking_flips_around_threshold():
    hi, cfg = _results(0.9, 0.1, 0.01, 0.90, 0.31)
    lo, cfg2 = _results(0.9, 0.1, 0.01, 0.90, 0.29)
    assert check_gates(hi, cfg)["r"]["ranking_ok"] is True
    assert check_gates(lo, cfg2)["r"]["ranking_ok"] is False


def test_gate_nan_metrics_do_not_pass():
    res, cfg = _results(float("nan"), 0.1, float("nan"), float("nan"), float("nan"))
    g = check_gates(res, cfg)["r"]
    assert g["ece_ok"] is False and g["coverage_ok"] is False and g["ranking_ok"] is False
    assert g["beats_all_baselines"] is False


# ===================================================================== #
# Integration: evaluate() end-to-end against a real multi-regime bundle #
# ===================================================================== #
@pytest.fixture(scope="module")
def eval_bundle(tmp_path_factory):
    warnings.filterwarnings("ignore")
    from cellfate.data import DataConfig, QCConfig, SyntheticSource
    from cellfate.data import run as build_run
    from cellfate.training.train_model import TrainConfig
    from cellfate.training.train_model import run as train_run

    set_global_seed(0)
    root = tmp_path_factory.mktemp("eval_bundle")
    cfg = DataConfig(
        out=str(root), gene_panel=str(root / "panel.json"), n_genes=64,
        qc=QCConfig(min_genes=5, max_mito_frac=0.5), label_tau=0.5, clock="random",
        deconfound=True, split_fracs=(0.6, 0.2, 0.1, 0.1),
        split_regimes=("scaffold", "cell_line", "both"), primary_regime="cell_line", seed=0,
    )
    # n_scaffold_families coprime with 3 decouples scaffold from outcome class, so
    # held-out splits are multi-class (see Document 5 notes).
    sources = [
        SyntheticSource(name="synth", n_lines=5, n_compounds=21, n_cells_per_condition=4,
                        n_filler_genes=50, n_scaffold_families=7, seed=1),
        SyntheticSource(name="tahoe", n_lines=4, n_compounds=21, n_cells_per_condition=4,
                        n_filler_genes=50, n_scaffold_families=7, seed=2),
    ]
    build_run(cfg, sources=sources)
    train_run(TrainConfig(
        dataset_dir=str(root), out=str(root), regime="cell_line",
        d_cell=48, d_u=48, latent_dim=48, p_drop=0.2, lr=3e-3, epochs=8,
        patience=4, batch_size=256, ensemble_size=2, base_seed=0,
        conformal_levels=(0.90,), device="cpu",
    ))
    return root


def test_evaluate_writes_reports_and_wellformed_gates(eval_bundle):
    out = eval_bundle / "reports"
    cfg = EvalConfig(bundle=str(eval_bundle), dataset=str(eval_bundle), out=str(out),
                     regimes=("cell_line", "scaffold"),
                     baselines=("mean", "ridge", "x_only", "predict_control"))
    gates = evaluate(cfg)

    # every reported regime carries all four boolean gate keys
    assert set(gates) <= {"cell_line", "scaffold"}
    assert "cell_line" in gates
    for _regime, g in gates.items():
        assert set(g) == {"beats_all_baselines", "ece_ok", "coverage_ok", "ranking_ok"}
        assert all(isinstance(v, bool) for v in g.values())

    # reports on disk
    assert (out / "cell_line.json").exists() and (out / "cell_line.md").exists()
    assert (out / "summary.json").exists() and (out / "summary.md").exists()
    summary = json.loads((out / "summary.json").read_text())
    assert summary["external"]["status"] == "not_available"  # no external data configured
    # regimes reported separately -- no single pooled number
    assert "cell_line" in summary["regimes"]


def test_cell_line_regime_is_multiclass_with_finite_metrics(eval_bundle):
    from cellfate.common.io import ArtifactPaths
    from cellfate.evaluation.data import gather_split

    paths = ArtifactPaths.of(eval_bundle)
    test = gather_split(paths, "cell_line", Split.TEST.value)
    assert test.n > 0
    assert len(np.unique(test.y_cls)) == 3  # decoupled scaffolds -> all classes present

    out = eval_bundle / "reports"
    R = json.loads((out / "cell_line.json").read_text())
    # model AUROC/PR-AUC are finite on a genuinely multi-class test split
    assert all(np.isfinite(R["model"][f"auroc_{c}"]) for c in range(3))
    assert all(np.isfinite(R["model"][f"prauc_{c}"]) for c in range(3))


def test_model_beats_trivial_baseline_on_cell_line(eval_bundle):
    # the model must clear the central-tendency baselines even if it cannot beat the
    # near-perfect linear baselines on trivially-separable synthetic data
    out = eval_bundle / "reports"
    R = json.loads((out / "cell_line.json").read_text())
    model_prauc = mean_finite(R["model"][f"prauc_{c}"] for c in range(3))
    mean_prauc = mean_finite(R["mean"][f"prauc_{c}"] for c in range(3))
    ctrl_prauc = mean_finite(R["predict_control"][f"prauc_{c}"] for c in range(3))
    assert model_prauc > mean_prauc + 0.1
    assert model_prauc > ctrl_prauc + 0.1


# --------------------------------------------------------------------------- #
# E-distance (measured effect size for RES ranking, Goal 4.5)                   #
# --------------------------------------------------------------------------- #
def test_energy_distance_zero_for_identical_and_grows_with_separation():
    import numpy as np

    from cellfate.evaluation.metrics import energy_distance
    rng = np.random.default_rng(0)
    a = rng.normal(0, 1, (200, 8))
    assert abs(energy_distance(a, a)) < 1e-6                     # identical -> 0
    d = [energy_distance(a, rng.normal(s, 1, (200, 8))) for s in (0.5, 2.0, 5.0)]
    assert all(d[i] < d[i + 1] for i in range(len(d) - 1))       # grows with shift
    assert all(x >= 0 for x in d)                               # non-negative


def test_edistance_to_control_ranks_groups():
    import numpy as np

    from cellfate.evaluation.metrics import edistance_to_control
    rng = np.random.default_rng(1)
    # control at 0; "near" shifted a little, "far" shifted a lot
    X = np.vstack([rng.normal(0, 1, (150, 10)),
                   rng.normal(1, 1, (150, 10)),
                   rng.normal(5, 1, (150, 10))])
    groups = np.array(["ctrl"] * 150 + ["near"] * 150 + ["far"] * 150)
    ed = edistance_to_control(X, groups, "ctrl", n_pcs=8, max_cells=150)
    assert set(ed) == {"near", "far"}          # control excluded
    assert ed["far"] > ed["near"] > 0          # bigger shift -> bigger E-distance
