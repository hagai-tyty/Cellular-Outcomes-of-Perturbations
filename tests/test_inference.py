"""Document 4 (cellfate.inference): RES invariants, OOD/conformal wiring, the
worked four-case ranking, the MC-dropout fix, and full end-to-end serving against
a real bundle built by Documents 1-3. Adversarial and limit-focused."""

from __future__ import annotations

import warnings

import numpy as np
import pytest
import torch
from pydantic import ValidationError

from cellfate.common import io
from cellfate.common.constants import SCHEMA_VERSION, Split
from cellfate.common.errors import ContractViolation, GenePanelMismatch, SchemaError
from cellfate.common.io import ArtifactPaths, load_splits
from cellfate.common.schemas import ResParams
from cellfate.common.seeding import set_global_seed
from cellfate.inference import (
    Predictor,
    Request,
    Response,
    compute_res,
    compute_res_batch,
    enable_mc_dropout,
    interval,
    intervals,
    predict_one,
    score_requests,
    score_shard,
)
from cellfate.inference.encode import _descriptor_to_fp
from cellfate.inference.res import (
    APPROVED,
    REJECTED_NO_REJUVENATION,
    REJECTED_OOD,
    REJECTED_UNSAFE,
    _sigmoid,
)
from cellfate.models import CellFateNet

P = ResParams()  # defaults: tau_safe=0.85, w=0.03, k=2, kappa=5, z_conf=1, lam=0


# ===================================================================== #
# Pure-unit tests: RES score (no bundle needed)                         #
# ===================================================================== #
def test_res_in_unit_interval_random():
    rng = np.random.default_rng(0)
    for _ in range(5000):
        S = rng.uniform(0, 1)
        P_loss = rng.uniform(0, 1 - S) if rng.random() < 0.5 else rng.uniform(0, 1)
        mu = rng.uniform(-30, 30)
        sig = rng.uniform(0, 10)
        res, status = compute_res(S, min(P_loss, 1.0), mu, sig, True, P)
        assert 0.0 <= res < 1.0, (S, P_loss, mu, sig, res)
        assert status in (APPROVED, REJECTED_UNSAFE, REJECTED_NO_REJUVENATION)
        assert 0.0 <= 10.0 * res < 10.0


def test_res_monotone_nondecreasing_in_safety():
    rng = np.random.default_rng(1)
    for _ in range(2000):
        mu, sig, pl = rng.uniform(-20, -1), rng.uniform(0, 2), rng.uniform(0, 0.3)
        s1 = rng.uniform(0, 0.99)
        s2 = rng.uniform(s1, 1.0)
        r1, _ = compute_res(s1, pl, mu, sig, True, P)
        r2, _ = compute_res(s2, pl, mu, sig, True, P)
        assert r2 + 1e-12 >= r1, (s1, s2, r1, r2)


def test_res_monotone_nonincreasing_in_uncertainty():
    rng = np.random.default_rng(2)
    for _ in range(2000):
        S = rng.uniform(0.8, 1.0)
        mu, pl = rng.uniform(-20, -1), rng.uniform(0, 0.3)
        sg1 = rng.uniform(0, 5)
        sg2 = rng.uniform(sg1, 10)
        r1, _ = compute_res(S, pl, mu, sg1, True, P)
        r2, _ = compute_res(S, pl, mu, sg2, True, P)
        assert r2 <= r1 + 1e-12, (sg1, sg2, r1, r2)


def test_res_zero_when_no_confident_rejuvenation():
    # upper age bound mu + z*sigma >= 0  ->  R_eff = 0  ->  RES = 0
    rng = np.random.default_rng(3)
    for _ in range(1000):
        sig = rng.uniform(0, 5)
        mu = rng.uniform(-P.z_conf * sig, 30)  # guarantees mu + z*sig >= 0
        res, status = compute_res(0.99, 0.01, mu, sig, True, P)
        assert res == 0.0
        if 0.99 >= P.tau_safe - 3 * P.w:
            assert status == REJECTED_NO_REJUVENATION


def test_res_ood_short_circuits_to_zero():
    res, status = compute_res(0.99, 0.0, -20.0, 0.1, in_dist=False, p=P)
    assert res == 0.0 and status == REJECTED_OOD


def test_res_confident_modest_outranks_uncertain_large():
    # the whole point of the v2 score
    modest_confident, _ = compute_res(0.92, 0.05, mu_age=-6.0, sigma_age=0.5, in_dist=True, p=P)
    large_uncertain, _ = compute_res(0.92, 0.05, mu_age=-12.0, sigma_age=8.0, in_dist=True, p=P)
    assert modest_confident > large_uncertain > 0.0


def test_res_batch_matches_scalar():
    rng = np.random.default_rng(4)
    n = 500
    S = rng.uniform(0, 1, n)
    pl = rng.uniform(0, 1, n)
    mu = rng.uniform(-20, 20, n)
    sg = rng.uniform(0, 8, n)
    ind = rng.random(n) > 0.2
    res_b, st_b = compute_res_batch(S, pl, mu, sg, ind, P)
    for j in range(n):
        r, s = compute_res(S[j], pl[j], mu[j], sg[j], bool(ind[j]), P)
        assert abs(r - res_b[j]) < 1e-9
        assert s == st_b[j]


def test_worked_four_case_example():
    """C > A; B unsafe; D no-rejuvenation; E out-of-distribution."""
    A = compute_res(0.92, 0.05, -6.0, 1.0, True, P)     # confident modest
    B = compute_res(0.60, 0.35, -10.0, 1.0, True, P)    # unsafe (S < tau-3w)
    C = compute_res(0.95, 0.03, -15.0, 1.0, True, P)    # confident large
    D = compute_res(0.95, 0.03, -2.0, 3.0, True, P)     # uncertain -> R_eff=0
    E = compute_res(0.95, 0.03, -15.0, 1.0, False, P)   # OOD
    assert A[1] == APPROVED and C[1] == APPROVED
    assert C[0] > A[0] > 0.0
    # UNSAFE is carried by status; RES is suppressed to ~0 by the SMOOTH floor (no cliff)
    assert B[1] == REJECTED_UNSAFE and B[0] < 1e-3
    # NO_REJUVENATION (g=0) and OOD (short-circuit) force RES to exactly 0
    assert D[1] == REJECTED_NO_REJUVENATION and D[0] == 0.0
    assert E[1] == REJECTED_OOD and E[0] == 0.0
    assert round(10 * B[0], 2) == 0.0  # and it presents as 0.0


def test_sigmoid_overflow_safe():
    x = np.array([-1e6, -1e3, -50, 0, 50, 1e3, 1e6])
    y = _sigmoid(x)
    assert np.all(np.isfinite(y)) and np.all((y >= 0) & (y <= 1))
    assert y[0] == pytest.approx(0.0) and y[-1] == pytest.approx(1.0)


# ===================================================================== #
# Pure-unit tests: conformal interval & encoding                        #
# ===================================================================== #
def test_conformal_interval_scalar_and_batch():
    assert interval(-3.0, 0.5) == [-3.5, -2.5]
    got = intervals(np.array([-3.0, 1.0]), 0.5)
    assert np.allclose(got, [[-3.5, -2.5], [0.5, 1.5]])


def test_conformal_interval_ordered_and_centered():
    rng = np.random.default_rng(5)
    for _ in range(1000):
        mu, q = rng.uniform(-20, 20), rng.uniform(0, 5)
        lo, hi = interval(mu, q)
        assert lo <= hi and abs((lo + hi) / 2 - mu) < 1e-9


def test_encode_rejects_non_chem_modality():
    with pytest.raises(ContractViolation, match="only chemical"):
        _descriptor_to_fp("CCO", modality=__import__("cellfate.common.constants", fromlist=["Modality"]).Modality.GENETIC)


def test_encode_rejects_bad_fingerprint_length():
    with pytest.raises(ContractViolation, match="2048"):
        _descriptor_to_fp([0.0, 1.0, 0.0], modality=__import__("cellfate.common.constants", fromlist=["Modality"]).Modality.CHEM)


def test_encode_smiles_and_bitvector_paths():
    from cellfate.common.constants import Modality
    fp_smiles = _descriptor_to_fp("CC(=O)Oc1ccccc1C(=O)O", Modality.CHEM)
    assert fp_smiles.shape == (2048,)
    bits = np.zeros(2048, dtype=np.float32)
    bits[[1, 5, 9]] = 1.0
    fp_bits = _descriptor_to_fp(bits.tolist(), Modality.CHEM)
    assert np.array_equal(fp_bits, bits)


# ===================================================================== #
# Pure-unit test: the MC-dropout fix                                    #
# ===================================================================== #
def test_enable_mc_dropout_only_dropout_stochastic():
    net = CellFateNet(g=32, latent_dim=32, p_drop=0.3)
    enable_mc_dropout(net)
    n_dropout = 0
    for m in net.modules():
        if isinstance(m, torch.nn.Dropout):
            assert m.training is True
            n_dropout += 1
        elif isinstance(m, (torch.nn.LayerNorm, torch.nn.Linear)):
            assert m.training is False, f"{type(m).__name__} should stay in eval"
    assert n_dropout >= 1


def test_response_schema_forbids_extra_fields():
    resp = Response(
        status=APPROVED, rejuvenation_efficacy_score=1.0, p_identity_preserved=0.9,
        p_identity_loss=0.05, p_apoptosis=0.05, delta_age_mean=-3.0,
        delta_age_interval=[-3.5, -2.5], in_distribution=True, epistemic_std=0.4,
        predictive_entropy=0.3,
    )
    assert resp.warning is None
    with pytest.raises(ValidationError):
        Response(status="x", rejuvenation_efficacy_score=1.0, p_identity_preserved=0.9,
                 p_identity_loss=0.05, p_apoptosis=0.05, delta_age_mean=-3.0,
                 delta_age_interval=[-3.5, -2.5], in_distribution=True, epistemic_std=0.4,
                 predictive_entropy=0.3, bogus=1)


# ===================================================================== #
# Integration tests against a real bundle (session-scoped fixture)      #
# ===================================================================== #
@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    warnings.filterwarnings("ignore")
    from cellfate.data import DataConfig, QCConfig, SyntheticSource
    from cellfate.data import run as build_run
    from cellfate.training.train_model import TrainConfig
    from cellfate.training.train_model import run as train_run

    set_global_seed(0)
    root = tmp_path_factory.mktemp("bundle")
    cfg = DataConfig(
        out=str(root), gene_panel=str(root / "panel.json"), n_genes=64,
        qc=QCConfig(min_genes=5, max_mito_frac=0.5), label_tau=0.5, clock="random",
        deconfound=True, split_fracs=(0.6, 0.2, 0.1, 0.1),
        split_regimes=("scaffold",), primary_regime="scaffold", seed=0,
    )
    sources = [
        SyntheticSource(name="synth", n_lines=2, n_compounds=10, n_cells_per_condition=10, n_filler_genes=80, seed=1),
        SyntheticSource(name="tahoe", n_lines=2, n_compounds=6, n_cells_per_condition=8, n_filler_genes=80, seed=2),
    ]
    build_run(cfg, sources=sources)
    tc = TrainConfig(
        dataset_dir=str(root), out=str(root), regime="scaffold",
        d_cell=48, d_u=48, latent_dim=48, p_drop=0.2, lr=3e-3, epochs=8,
        patience=4, batch_size=256, ensemble_size=2, base_seed=0,
        conformal_levels=(0.90,), device="cpu",
    )
    train_run(tc)
    return root


def test_predictor_loads_and_is_consistent(bundle):
    pred = Predictor(bundle)
    assert len(pred.members) == 2
    assert pred.temperature > 0 and pred.q > 0
    assert pred.conformal_level == 0.90
    assert pred.meta.gene_panel_hash == pred.scalers.params.gene_panel_hash


def test_request_path_matches_shard_path(bundle):
    pred = Predictor(bundle)
    paths = ArtifactPaths.of(bundle)
    shard = sorted(paths.shards_dir.glob("*.parquet"))[0]
    arr = io.shard_to_numpy(io.read_shard(shard))
    i = 3
    edt = arr["dose_time"][i]
    req = Request(
        X_raw=arr["X"][i].tolist(), u_modality="chem",
        u_descriptor=arr["u_chem_fp"][i].astype(float).tolist(),
        dose_uM=float(10.0 ** edt[0]), time_h=float(np.exp(edt[1])),
    )
    r_single = predict_one(pred, req)
    r_enc = pred.predict_encoded(arr["X"][i:i+1], arr["u_chem_fp"][i:i+1], arr["dose_time"][i:i+1])[0]
    assert r_single.p_identity_preserved == pytest.approx(round(r_enc["S"], 3))
    assert r_single.delta_age_mean == pytest.approx(round(r_enc["mu_age"], 2))
    assert r_single.in_distribution == r_enc["in_dist"]


def test_score_shard_all_responses_valid(bundle):
    pred = Predictor(bundle)
    paths = ArtifactPaths.of(bundle)
    shard = sorted(paths.shards_dir.glob("*.parquet"))[0]
    responses, cell_ids = score_shard(pred, shard)
    assert len(responses) == len(cell_ids) > 0
    for r in responses:
        assert isinstance(r, Response)
        assert 0.0 <= r.rejuvenation_efficacy_score < 10.0
        lo, hi = r.delta_age_interval
        assert lo <= hi
        assert r.status in (APPROVED, REJECTED_OOD, REJECTED_UNSAFE, REJECTED_NO_REJUVENATION)
        assert (r.warning is not None) == (r.status == REJECTED_OOD)
        # APPROVED implies it cleared the (smooth) safety floor
        if r.status == APPROVED:
            assert r.p_identity_preserved >= P.tau_safe - 3 * P.w - 1e-9


def test_ood_cell_is_rejected(bundle):
    pred = Predictor(bundle)
    paths = ArtifactPaths.of(bundle)
    for shard in sorted(paths.shards_dir.glob("*.parquet")):
        responses, _ = score_shard(pred, shard)
        ood = [r for r in responses if r.status == REJECTED_OOD]
        if ood:
            r = ood[0]
            assert not r.in_distribution and r.rejuvenation_efficacy_score == 0.0
            assert r.warning is not None
            return
    pytest.skip("no OOD cells surfaced in this fixture bundle")


def test_batch_and_single_agree(bundle):
    pred = Predictor(bundle)
    paths = ArtifactPaths.of(bundle)
    arr = io.shard_to_numpy(io.read_shard(sorted(paths.shards_dir.glob("*.parquet"))[0]))
    reqs = []
    for i in range(5):
        edt = arr["dose_time"][i]
        reqs.append(Request(X_raw=arr["X"][i].tolist(), u_modality="chem",
                            u_descriptor=arr["u_chem_fp"][i].astype(float).tolist(),
                            dose_uM=float(10.0 ** edt[0]), time_h=float(np.exp(edt[1]))))
    batch = score_requests(pred, reqs)
    for k, req in enumerate(reqs):
        one = predict_one(pred, req)
        assert one.model_dump() == batch[k].model_dump()


def test_determinism_ensemble_mode(bundle):
    pred = Predictor(bundle)
    paths = ArtifactPaths.of(bundle)
    arr = io.shard_to_numpy(io.read_shard(sorted(paths.shards_dir.glob("*.parquet"))[0]))
    a = pred.predict_encoded(arr["X"][:8], arr["u_chem_fp"][:8], arr["dose_time"][:8])
    b = pred.predict_encoded(arr["X"][:8], arr["u_chem_fp"][:8], arr["dose_time"][:8])
    assert a == b


def test_both_inference_modes_carry_their_own_sigma_scale(bundle):
    """Each mode's sigma_age is a different quantity, so each needs its own factor.

    Neither may be left at 1.0 (raw spread = overconfident, the defect Stage 1 removes), and
    neither may borrow the other's (scales the wrong spread).
    """
    ens = Predictor(bundle, mode="ensemble")
    mc = Predictor(bundle, mode="mc_dropout", T=8)

    assert ens.sigma_scale > 1.0, "ensemble mode is serving a raw, uncalibrated spread"
    assert mc.sigma_scale > 1.0, "mc_dropout mode is serving a raw, uncalibrated spread"
    assert ens.sigma_scale != mc.sigma_scale, (
        "identical factors for two different spreads means one mode borrowed the other's"
    )


def test_predictor_refuses_a_mode_the_bundle_was_never_calibrated_for(bundle, tmp_path):
    """A bundle not calibrated for THIS mode must fail loudly, not serve raw sigma.

    Status comes from `sigma_calibrated_modes`, never from the factor's value: the factor is
    clamped at 1.0, so 1.0 means EITHER "measured, spread already adequate" OR "never
    measured". The second case below is the one that distinction exists for -- inferring from
    the value would refuse a correctly-calibrated bundle.
    """
    import shutil

    from cellfate.common.errors import ConfigError

    copy_root = tmp_path / "bundle_copy"
    shutil.copytree(bundle, copy_root)
    paths = ArtifactPaths.of(copy_root)
    conf = io.load_conformal(paths)

    # (a) mode genuinely absent from the calibrated set -> refuse
    io.save_conformal(paths, conf.model_copy(
        update={"sigma_calibrated_modes": ["ensemble"]}))
    Predictor(copy_root, mode="ensemble")
    with pytest.raises(ConfigError, match="mc_dropout"):
        Predictor(copy_root, mode="mc_dropout")

    # (b) mode calibrated but its factor clamped to 1.0 -> must STILL load
    io.save_conformal(paths, conf.model_copy(
        update={"sigma_calibrated_modes": ["ensemble", "mc_dropout"], "sigma_scale_mc": 1.0}))
    p = Predictor(copy_root, mode="mc_dropout", T=4)
    assert p.sigma_scale == 1.0

    # (c) legacy bundle with no per-mode record -> unchanged behaviour, both modes load
    io.save_conformal(paths, conf.model_copy(update={"sigma_calibrated_modes": []}))
    Predictor(copy_root, mode="ensemble")
    Predictor(copy_root, mode="mc_dropout", T=4)


def test_mc_dropout_is_single_batched_call(bundle):
    """T passes must be ONE tiled forward, never a per-sample loop."""
    pred = Predictor(bundle, mode="mc_dropout", T=40)
    paths = ArtifactPaths.of(bundle)
    arr = io.shard_to_numpy(io.read_shard(sorted(paths.shards_dir.glob("*.parquet"))[0]))
    member = pred.members[0]
    calls = {"n": 0}
    orig = member.forward

    def counting(*a, **k):
        calls["n"] += 1
        return orig(*a, **k)

    member.forward = counting
    try:
        rows = pred.predict_encoded(arr["X"][:4], arr["u_chem_fp"][:4], arr["dose_time"][:4])
    finally:
        member.forward = orig
    # one latent pass + one tiled MC pass = 2 forwards, independent of T (=40)
    assert calls["n"] == 2, calls["n"]
    assert all(r["sigma_age"] >= 0.0 for r in rows)


def test_conformal_interval_covers_calibration_set(bundle):
    """The shipped interval reproduces the calibrated coverage on the calib split."""
    pred = Predictor(bundle)
    paths = ArtifactPaths.of(bundle)
    assign = load_splits(paths, "scaffold")
    Xs, fps, dts, ages = [], [], [], []
    for shard in sorted(paths.shards_dir.glob("*.parquet")):
        arr = io.shard_to_numpy(io.read_shard(shard))
        for j, cid in enumerate(arr["cell_id"]):
            if str(assign.get(cid)) == Split.CALIB.value and bool(arr["age_mask"][j]):
                Xs.append(arr["X"][j])
                fps.append(arr["u_chem_fp"][j])
                dts.append(arr["dose_time"][j])
                ages.append(float(arr["y_age"][j]))
    if len(ages) < 20:
        pytest.skip("too few age-valid calib cells in fixture")
    rows = pred.predict_encoded(np.array(Xs), np.array(fps), np.array(dts))
    mu = np.array([r["mu_age"] for r in rows])
    coverage = float(np.mean(np.abs(mu - np.array(ages)) <= pred.q))
    # by construction the empirical quantile gives >= the target level on calib
    assert coverage >= 0.85, coverage


def test_gene_panel_mismatch_detected(bundle, tmp_path):
    """Corrupting the shipped scaler hash trips the consistency check."""
    import shutil

    from cellfate.common.scalers import Scalers
    clone = tmp_path / "clone"
    shutil.copytree(bundle, clone)
    paths = ArtifactPaths.of(clone)
    sc = Scalers.load(paths.bundle_scalers_file)
    sc.params.gene_panel_hash = "deadbeefdeadbeef"
    sc.save(paths.bundle_scalers_file)
    with pytest.raises(GenePanelMismatch):
        Predictor(clone)


def test_schema_version_mismatch_rejected(bundle, tmp_path):
    """A bundle written under a different SCHEMA_VERSION is rejected fail-fast."""
    import shutil
    clone = tmp_path / "clone_ver"
    shutil.copytree(bundle, clone)
    paths = ArtifactPaths.of(clone)
    meta = io.load_bundle_meta(paths).model_dump()
    meta["schema_version"] = SCHEMA_VERSION + "-stale"
    io.write_json(paths.bundle_meta_file, meta)
    with pytest.raises(SchemaError):
        Predictor(clone)
