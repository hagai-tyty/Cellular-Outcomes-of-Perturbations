"""Training orchestrator (Document 3, S7).

Loads the Document-2 dataset, trains the ensemble, calibrates (temperature +
conformal), fits the OOD reference, and writes a complete, contract-valid
deployment bundle plus a metrics report. ``run`` is plain-Python and testable;
``cli`` is the Hydra entry point.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import torch

from cellfate.common import constants as C
from cellfate.common import io
from cellfate.common.calibration import apply_platt, platt_safe
from cellfate.common.errors import ConfigError
from cellfate.common.io import ArtifactPaths
from cellfate.common.logging import get_logger, log_event
from cellfate.common.scalers import Scalers
from cellfate.common.schemas import ResParams
from cellfate.models import CellFateNet

from .bundle import assemble_bundle
from .calibrate import fit_platt_binary, fit_temperature
from .conformal import coverage, fit_conformal
from .dataset import AM_I, FP_I, YA_I, YC_I, load_split_tensors
from .metrics import ece, soft_nll
from .ood import fit_ood
from .train import (
    ensemble_age,
    ensemble_logits,
    ensemble_probs,
    member_outputs,
    train_ensemble,
)
from .xdonor_calib import (
    SIGMA_SCALE_MODE,
    crossdonor_stats,
    n_train_donors,
    save_xstats,
    sigma_scale_factor,
)

log = get_logger("cellfate.training")


@dataclass
class TrainConfig:
    dataset_dir: str
    out: str | None = None              # bundle root; defaults to dataset_dir
    regime: str = "scaffold"
    # model
    d_cell: int = 256
    d_u: int = 256
    latent_dim: int = 256
    p_drop: float = 0.2
    # optimisation
    lr: float = 3e-4
    wd: float = 1e-5
    epochs: int = 60
    patience: int = 8
    min_delta: float = 1e-4
    clip: float = 1.0
    batch_size: int = 512
    ensemble_size: int = 5
    base_seed: int = 0
    # losses
    focal_gamma: float = 2.0
    class_weight_beta: float = 0.999
    huber_delta: float = 2.0
    # calibration / device / RES
    conformal_levels: tuple[float, ...] = (0.90,)
    device: str = "cpu"
    res: dict = field(default_factory=dict)
    # Stage 1b: fit temperature, q and sigma_scale on inner-LODO cross-donor statistics
    # instead of in-distribution splits. (The OOD reference is NOT among them -- see the
    # comment at its call site.) Setting this False is the documented rollback
    # (STAGE_1 S1b.5) -- every call site keeps its original in-distribution branch.
    xdonor_calibration: bool = True
    # sigma_scale calibrates the ENSEMBLE spread, so it is only valid for this mode.
    inference_mode: str = "ensemble"
    # dropout passes used to calibrate the mc_dropout sigma_scale. Must match Predictor's
    # default T, or that mode's factor is fitted to a spread inference never produces.
    mc_dropout_T: int = 50


def _resolve_device(requested: str) -> str:
    if requested == "cuda" and not torch.cuda.is_available():
        warnings.warn("cuda requested but unavailable; using cpu", stacklevel=2)
        return "cpu"
    return requested


def run(cfg: TrainConfig) -> dict:
    """Train + calibrate + bundle; returns a summary/metrics dict."""
    device = _resolve_device(cfg.device)
    data_paths = ArtifactPaths.of(cfg.dataset_dir)
    out_root = cfg.out or cfg.dataset_dir
    scalers = Scalers.load(data_paths.scalers_file)
    g = len(scalers.params.x_mean)

    train_ds = load_split_tensors(data_paths, scalers, cfg.regime, "train")
    val_ds = load_split_tensors(data_paths, scalers, cfg.regime, "val")
    calib_ds = load_split_tensors(data_paths, scalers, cfg.regime, "calib")
    if len(train_ds) == 0:
        raise ConfigError(f"empty train split for regime {cfg.regime!r}")
    log_event(log, "data.loaded", g=g, n_train=len(train_ds),
              n_val=len(val_ds), n_calib=len(calib_ds))

    # perturbation modality: TF-cocktail vs chemical fingerprint. The width comes
    # from the loaded tensor; the kind from the shards' u_modality. Both go into
    # the model config, so the member arch (and inference) get the right encoder.
    def _detect_pert_kind() -> str:
        for shard in sorted(data_paths.shards_dir.glob("*.parquet")):
            mod = io.shard_to_numpy(io.read_shard(shard))["u_modality"]
            return "tf" if (len(mod) and str(mod[0]) == "tf") else "chem"
        return "chem"

    pert_kind = _detect_pert_kind()
    n_pert = int(train_ds.tensors[FP_I].shape[1])
    model_cfg = dict(g=g, d_cell=cfg.d_cell, d_u=cfg.d_u,
                     latent_dim=cfg.latent_dim, p_drop=cfg.p_drop,
                     pert_kind=pert_kind, n_pert=n_pert)

    def make_model() -> CellFateNet:
        return CellFateNet(**model_cfg)

    members, val_losses = train_ensemble(make_model, train_ds, val_ds, cfg, device)
    res_params = ResParams(**cfg.res)

    # -- cross-donor calibration statistics (inner LODO over training donors) ---------- #
    # Every calibrator below is fitted on THESE rather than on val/calib/train, because
    # those are all in-distribution and deployment is not (MASTER_PLAN S5a).
    xstats = None
    n_donors = n_train_donors(train_ds)
    if not cfg.xdonor_calibration:
        log.warning("xdonor_calibration disabled; calibrating IN-DISTRIBUTION (rollback path)")
    elif n_donors < 2:
        log.warning(
            "only %d donor(s) in the training split; inner-LODO impossible, falling back to "
            "IN-DISTRIBUTION calibration. Intervals and fate probabilities will NOT be "
            "trustworthy out-of-donor.", n_donors,
        )
    else:
        xstats = crossdonor_stats(train_ds, val_ds, make_model, cfg, device,
                                  mc_T=cfg.mc_dropout_T, calib_ds=calib_ds)

    # -- temperature: cross-donor logits, falling back to val (then calib) ------------- #
    # cal_* stay defined regardless: the report contrasts in-distribution against
    # cross-donor calibration, which is exactly what Stage 1 is measuring.
    cal_ds = val_ds if len(val_ds) else calib_ds
    if len(cal_ds):
        cal_logits = ensemble_logits(members, cal_ds, device).numpy()
        cal_target = cal_ds.tensors[YC_I].numpy()
    else:
        cal_logits = cal_target = None

    # -- fate calibration: Platt on P(safe), fitted on ALL available held-out cells ------ #
    # WHAT: `res.py` consumes S and P_loss, STAGE_3 §0.1 needs a risk threshold on P(unsafe),
    # and scorecard grades binary ECE on P(safe). MASTER_PLAN §5a names the defective quantity
    # as "S, P_loss" and records "YES -- Platt halves it" (T8.2); STAGE_1's ≲0.17 bar comes from
    # that measurement. Run 2 optimised multi-class NLL instead and regressed the graded metric
    # on every fold (0.281 -> 0.364). Platt subsumes a temperature (its slope IS one on the
    # binary logit) and adds the intercept a scalar cannot express, so `temperature` stays 1.0
    # rather than stacking two interacting calibrators.
    #
    # HOW MUCH DATA: everything held out -- the calib/val split AND the cross-donor pool.
    #
    # This is NOT a departure from Stage 1's cross-donor principle. That principle says to
    # calibrate on data whose ERROR REGIME matches deployment, and it is decisive for ΔAge:
    # ~4 yr in-distribution against ~14 yr out-of-donor. The premise is NOT met for fate --
    #   * discrimination 0.929-0.940 in-distribution vs 0.96-1.00 out-of-donor (T8.1): no
    #     degradation, if anything the held-out donors are easier;
    #   * a CALIB-FITTED Platt halves out-of-donor ECE on 4 of 5 folds (T8.2): it transfers.
    # STAGE_1's <=0.17 bar is itself derived from that calib-fitted Platt -- T8.2's "ECE raw"
    # and "ECE recal" columns are, cell for cell, scorecard's `fate_ece` and `fate_ece_platt`.
    # Holding the calibrator to a bar measured with a method we refused to use would be
    # incoherent; §1b.2's `fit_temperature(xstats...)` is the line that never matched §2's own
    # expected effect.
    #
    # Restricting to the pool alone would fit 2 parameters on ~103 cells while discarding
    # ~4,490 -- and because cells within a donor share that donor's offset, its EFFECTIVE n is
    # nearer 5 (the donor count) than 103. Run 2 measured the cost: cross-donor temperature came
    # out 30% WORSE than the in-distribution one, CI excluding zero. The pool-only variant is
    # still fitted below and REPORTED, so this stays a measurement rather than an assertion.
    from cellfate.common.schemas import TemperatureParams

    cal_probs = ensemble_probs(members, cal_ds, device).numpy() if len(cal_ds) else None
    fate_parts = [(p[:, C.SAFE_IDX], t[:, C.SAFE_IDX]) for p, t in (
        (cal_probs, cal_target),
        (None if xstats is None else xstats.probs_mean, None if xstats is None else xstats.targets),
    ) if p is not None and len(p)]

    if fate_parts:
        p_all = np.concatenate([p for p, _ in fate_parts])
        y_all = np.concatenate([y for _, y in fate_parts])
        a, b = fit_platt_binary(p_all, y_all)
        temperature = TemperatureParams(temperature=1.0, platt_a=a, platt_b=b)
        n_xd = 0 if xstats is None else int(len(xstats.probs_mean))
        fate_calib_n = {"total": int(p_all.size), "in_dist": int(p_all.size) - n_xd,
                        "xdonor": n_xd}
        log_event(log, "fate.calibrated", calibrator="platt_binary",
                  platt_a=round(a, 5), platt_b=round(b, 5), **fate_calib_n)
    else:
        temperature = TemperatureParams(temperature=1.0)
        log.warning("no held-out cells for fate calibration; shipping uncalibrated")
        fate_calib_n = {"total": 0, "in_dist": 0, "xdonor": 0}

    # -- the STRICT cross-donor variant, fitted and REPORTED but NOT shipped ------------ #
    # Stage 1's principle is that calibration must be fitted on data whose error regime matches
    # deployment. For ΔAge that is measured and decisive (4 yr in-distribution vs 14 yr
    # out-of-donor), so `q` and `sigma_scale` use the cross-donor pool alone. For FATE the
    # premise does not hold: discrimination is 0.929-0.940 in-distribution against 0.96-1.00
    # out-of-donor (T8.1) -- no degradation -- and a calib-fitted Platt halves out-of-donor ECE
    # on 4 of 5 folds (T8.2), i.e. it transfers. STAGE_1's <=0.17 bar was itself derived from
    # that calib-fitted Platt.
    #
    # Rather than assert that, MEASURE it: fit the pool-only calibrator too and record both, so
    # the principle is tested on every run instead of quietly dropped. The pool is ~103 cells
    # from 5 donors, and cells within a donor share that donor's offset, so its EFFECTIVE n is
    # nearer 5 than 103 -- which is what a 2-parameter fit has to work with.
    fate_alt = {}
    if xstats is not None and len(xstats.probs_mean):
        p_xd = xstats.probs_mean[:, C.SAFE_IDX]
        y_xd = xstats.targets[:, C.SAFE_IDX]
        a_xd, b_xd = fit_platt_binary(p_xd, y_xd)
        fate_alt = {
            "xdonor_only_platt_a": a_xd,
            "xdonor_only_platt_b": b_xd,
            "xdonor_only_n": int(p_xd.size),
            "xdonor_only_n_donors": int(xstats.n_donors),
            # both calibrators scored on the SAME held-out pool. The pool-only fit is scored on
            # the data it was fitted to, so its number is optimistic -- stated here so the
            # comparison is not read as fairer than it is.
            "xdonor_only_safe_ece_insample": _binary_ece(
                platt_safe(p_xd, a_xd, b_xd), y_xd),
            "shipped_safe_ece_on_pool": _binary_ece(
                platt_safe(p_xd, temperature.platt_a, temperature.platt_b), y_xd
            ) if temperature.has_platt else None,
        }
        log_event(log, "fate.xdonor_only_diagnostic", calibrator="platt_binary",
                  platt_a=round(a_xd, 5), platt_b=round(b_xd, 5), n=int(p_xd.size))

    # -- conformal + per-mode sigma scales: cross-donor residuals ---------------------- #
    # BOTH inference modes are calibrated, from the same held-out rows. `sigma_age` is the
    # ensemble spread in one and the T-pass dropout spread in the other; a mode left at 1.0
    # would serve raw, overconfident uncertainty -- the exact defect Stage 1 removes.
    if xstats is not None and xstats.has_residuals:
        abs_res = xstats.abs_residuals                                    # CROSS-DONOR
        lvl = float(cfg.conformal_levels[0])
        sigma_scale = sigma_scale_factor(xstats, res_params.z_conf, lvl, mode="ensemble")
        sigma_scale_mc = sigma_scale_factor(xstats, res_params.z_conf, lvl, mode="mc_dropout")
        # Record WHICH modes were measured. A factor of 1.0 means "already adequate" here, not
        # "unknown", and only this list can tell those apart at load time.
        calibrated_modes = [m for m, s in (("ensemble", xstats.sigma_pred),
                                           ("mc_dropout", xstats.sigma_pred_mc)) if s.size]
    else:
        if xstats is not None:
            log.warning("xdonor residuals empty; falling back to calib residuals")
        if len(calib_ds):
            age_pred = ensemble_age(members, calib_ds, device).numpy()
            ya = calib_ds.tensors[YA_I].numpy()
            am = calib_ds.tensors[AM_I].numpy().astype(bool)
            abs_res = np.abs(age_pred[am] - ya[am])
        else:
            abs_res = np.array([])
        sigma_scale = sigma_scale_mc = 1.0
        calibrated_modes = []          # legacy behaviour: no per-mode calibration claimed

    conformal = fit_conformal(abs_res, cfg.conformal_levels, sigma_scale=sigma_scale,
                              sigma_scale_mc=sigma_scale_mc,
                              sigma_scale_mode=SIGMA_SCALE_MODE,
                              calibrated_modes=calibrated_modes)

    # -- OOD reference: DELIBERATELY still the deployed model's train features --------- #
    # STAGE_1 S1b.2 Edit 4 says to fit this on `xstats.feats`. That is not implementable:
    # those features come from five independently-seeded INNER models, whose latent bases are
    # related only by arbitrary rotation/permutation. `OODDetector` compares the DEPLOYED
    # member[0]'s latent z against this stored mean/precision, so pooling incomparable
    # coordinates would make the Mahalanobis distance meaningless -- worse than the defect it
    # is meant to fix, and silently so.
    # The OOD refit is independent of the other three (S3), and S1b.4 already anticipates the
    # outcome: "the detector is uninformative regardless of fitting split -> disable the gate
    # (Stage 3d) rather than chase it". T15's AUC of 0.47 is a property of the representation,
    # not of the fitting split. Left unchanged pending that decision.
    ood = fit_ood(member_outputs(members[0], train_ds, device)[2].numpy())

    config_hash = io.hash_config({"model": model_cfg, "train": vars(cfg)})

    paths = assemble_bundle(
        out_root, members, temperature, conformal, res_params, ood,
        gene_panel_hash=scalers.params.gene_panel_hash,
        scalers_src=str(data_paths.scalers_file),
        model_cfg=model_cfg, train_cfg=_jsonable(vars(cfg)),
        deps_hash=io.deps_hash(), config_hash=config_hash,
    )

    # Persist the cross-donor pool: it costs ~35 min/fold to produce and was discarded once the
    # calibrators were fitted, so every calibration experiment cost another full LOOCV pass.
    if xstats is not None:
        save_xstats(paths.bundle_dir, xstats)

    metrics = _report(members, val_losses, train_ds, val_ds, calib_ds,
                      cal_logits, cal_target, temperature, conformal, abs_res, xstats,
                      cfg.mc_dropout_T, res_params.z_conf, fate_calib_n, fate_alt)
    io.write_json(paths.bundle_dir / C.BUNDLE_METRICS_FILENAME, metrics)
    log_event(log, "bundle.done", bundle=str(paths.bundle_dir),
              n_members=len(members), temperature=round(temperature.temperature, 4))
    return {"bundle": str(paths.bundle_dir), **metrics}


def _jsonable(d: dict) -> dict:
    return {k: (list(v) if isinstance(v, tuple) else v) for k, v in d.items()}


def _binary_ece(p: np.ndarray, y: np.ndarray, bins: int = 10) -> float:
    """Binary ECE on P(safe) -- the SAME definition scorecard.py grades as ``fate_ece``.

    Deliberately duplicated from `scorecard.py:_ece` rather than imported: the scorecard is a
    repo-root script, not part of the package, and the whole point of this metric is that the
    bundle's diagnostic and the graded number agree.
    """
    p = np.asarray(p, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    edges = np.linspace(0.0, 1.0, bins + 1)
    err = 0.0
    for i in range(bins):
        hi = edges[i + 1] if i < bins - 1 else 1.0 + 1e-9
        m = (p >= edges[i]) & (p < hi)
        if m.sum():
            err += (m.sum() / len(p)) * abs(p[m].mean() - y[m].mean())
    return float(err)


def _report(members, val_losses, train_ds, val_ds, calib_ds,
            cal_logits, cal_target, temperature, conformal, abs_res, xstats=None,
            mc_T: int = 50, res_z_conf: float = 1.0,
            fate_calib_n: dict | None = None,
            fate_alt: dict | None = None) -> dict:
    from scipy.special import softmax  # available via the data-stage deps

    out = {
        "n_members": len(members),
        "n_train": len(train_ds),
        "n_val": len(val_ds),
        "n_calib": len(calib_ds),
        "val_loss_mean": float(np.mean(val_losses)),
        "val_loss_per_member": [float(v) for v in val_losses],
        "temperature": temperature.temperature,
        "conformal_q": dict(conformal.q),
        "fate_calib_n": fate_calib_n or {},
        **(fate_alt or {}),
        "platt_a": temperature.platt_a,
        "platt_b": temperature.platt_b,
        "sigma_scale": float(conformal.sigma_scale),          # mode="ensemble"
        "sigma_scale_mc": float(conformal.sigma_scale_mc),    # mode="mc_dropout"
        "sigma_scale_mc_T": int(mc_T),   # the T it was fitted at; see mc_dropout_spread
        "sigma_scale_mode": conformal.sigma_scale_mode,
        # whether the calibrators actually saw cross-donor data. Recorded in the bundle so a
        # fallback is auditable after the fact rather than only visible in a log line.
        "xdonor_calibrated": xstats is not None and xstats.n_donors >= 2,
        "xdonor_n_donors": 0 if xstats is None else xstats.n_donors,
        "xdonor_n_residuals": 0 if xstats is None else int(xstats.abs_residuals.size),
        # per-donor pool composition: `q` is a quantile of these, so a lopsided pool means the
        # quantile describes one donor rather than cross-donor error
        "xdonor_residuals_per_donor": {} if xstats is None else xstats.residuals_per_donor,
        # per-donor error scale vs predicted spread -- decides whether an ADAPTIVE conformal
        # interval can meet the coverage bar that a single global q provably cannot
        "xdonor_donor_scales": {} if xstats is None else xstats.donor_scales,
    }
    # IN-DISTRIBUTION (calib/val). From Stage 1b the temperature is no longer FITTED here, so
    # `nll_after_temp <= nll_before_temp` is NO LONGER GUARANTEED on this split -- fit_temperature
    # only promises "never worse than T=1" on the data it saw. In-distribution NLL getting worse
    # is the expected cost of the trade: the model is under-confident in-distribution (baseline
    # T=0.54, which SHARPENS) and over-confident out-of-donor (needs T>1, which SOFTENS). One
    # scalar cannot serve both. Kept as the contrast, not as an invariant.
    if cal_logits is not None:
        labels = cal_target.argmax(axis=1)
        p_before = softmax(cal_logits, axis=1)
        p_after = softmax(cal_logits / temperature.temperature, axis=1)
        out["nll_before_temp"] = soft_nll(p_before, cal_target)
        out["nll_after_temp"] = soft_nll(p_after, cal_target)
        out["ece_before_temp"] = ece(p_before, labels)
        out["ece_after_temp"] = ece(p_after, labels)
    # OUT-OF-DONOR (the fitting set from Stage 1b). The "never worse" guarantee holds HERE,
    # and this is the contrast Stage 1 exists to move: ECE 0.28 -> ~0.13.
    if xstats is not None and xstats.has_logits:
        xl, xt = xstats.logits, xstats.targets
        xlabels = xt.argmax(axis=1)
        xp_before = softmax(xl, axis=1)
        xp_after = softmax(xl / temperature.temperature, axis=1)
        out["xdonor_nll_before_temp"] = soft_nll(xp_before, xt)
        out["xdonor_nll_after_temp"] = soft_nll(xp_after, xt)
        out["xdonor_ece_before_temp"] = ece(xp_before, xlabels)
        out["xdonor_ece_after_temp"] = ece(xp_after, xlabels)
    # ...and the BINARY P(safe) ECE, which is what scorecard.py actually grades as `fate_ece`.
    # Reported alongside the top-1 multi-class figure above so the bundle's own diagnostic and
    # the scorecard metric can never silently disagree again -- run 2's top-1 ECE improved
    # (0.269 -> 0.217) while the graded binary ECE regressed (0.281 -> 0.364).
    if xstats is not None and len(xstats.probs_mean):
        p_safe = xstats.probs_mean[:, C.SAFE_IDX]
        y_safe = xstats.targets[:, C.SAFE_IDX]
        out["xdonor_safe_ece_before"] = _binary_ece(p_safe, y_safe)
        if temperature.has_platt:
            cal = apply_platt(xstats.probs_mean, temperature.platt_a,
                              temperature.platt_b, C.SAFE_IDX)
            out["xdonor_safe_ece_after"] = _binary_ece(cal[:, C.SAFE_IDX], y_safe)
    # -- how much of the sigma fix actually reaches the cells that matter ---------------- #
    # `sigma_scale` is MULTIPLICATIVE, so it fixes the spread's MAGNITUDE but preserves its
    # SHAPE. A cell where the ensemble happens to agree keeps a near-zero sigma even after a
    # 6x scaling -- and RES consumes sigma via R_eff = max(0, -(mu + z*sigma)), so such a cell
    # is scored as if its ΔAge were near-certain and can be APPROVED on that basis, while its
    # true out-of-donor error is ~q. That is the permissive direction, the dangerous one.
    # MASTER_PLAN §5b-bis anticipated this and offered `R_eff = max(0, -(mu + q))` as the
    # "cleaner" alternative; STAGE_1 specified the rescaling instead. These ratios measure how
    # large the gap is, so the choice can be made on evidence in Stage 4 (Change C).
    if xstats is not None and xstats.sigma_pred.size and abs_res.size:
        lvl0 = conformal.levels[0]
        q0 = conformal.q[str(lvl0)]
        if q0 > 0:
            ratio = (xstats.sigma_pred * conformal.sigma_scale * res_z_conf) / q0
            p10, p50, p90 = (float(np.quantile(ratio, x)) for x in (0.10, 0.50, 0.90))
            out["xdonor_sigma_over_q_p10"] = p10   # « 1 means many cells still claim near-
            out["xdonor_sigma_over_q_p50"] = p50   #   certainty despite honest error ~q
            out["xdonor_sigma_over_q_p90"] = p90
            out["xdonor_sigma_under_half_q_frac"] = float((ratio < 0.5).mean())
    if abs_res.size:
        lvl = conformal.levels[0]
        # NOTE: with cross-donor residuals this is coverage ON THE FITTING SET, so it is a
        # sanity check, not evidence. Honest coverage comes from scorecard.py on held-out folds.
        out["conformal_coverage_calib"] = coverage(abs_res, conformal.q[str(lvl)])
        out["abs_residual_mean"] = float(np.mean(abs_res))
    return out


# --------------------------------------------------------------------------- #
# CLI                                                                         #
# --------------------------------------------------------------------------- #
def _config_from_omegaconf(cfg) -> TrainConfig:
    m, t = cfg.model, cfg.train
    return TrainConfig(
        dataset_dir=cfg.data.out,
        out=cfg.get("bundle_out", None) or cfg.data.out,
        regime=str(cfg.data.splits.get("primary", "scaffold")),
        d_cell=int(m.d_cell), d_u=int(m.d_u), latent_dim=int(m.latent_dim),
        p_drop=float(m.p_drop),
        lr=float(t.lr), wd=float(t.wd), epochs=int(t.epochs), patience=int(t.patience),
        min_delta=float(t.min_delta), clip=float(t.clip), batch_size=int(t.batch_size),
        ensemble_size=int(t.ensemble_size), base_seed=int(t.base_seed),
        focal_gamma=float(t.focal.gamma), class_weight_beta=float(t.focal.class_weight_beta),
        huber_delta=float(t.huber_delta),
        conformal_levels=(float(t.conformal_level),),
        device=str(cfg.get("device", "cpu")),
        res=dict(cfg.get("res", {})),
        xdonor_calibration=bool(t.get("xdonor_calibration", True)),
        inference_mode=str(cfg.get("inference_mode", "ensemble")),
    )


def cli() -> None:  # pragma: no cover - exercised in production via Hydra
    try:
        import hydra
        from omegaconf import DictConfig
    except ImportError as exc:
        raise ConfigError("hydra-core/omegaconf required for the CLI") from exc

    config_dir = str(Path(__file__).resolve().parents[3] / "configs")

    @hydra.main(version_base=None, config_path=config_dir, config_name="config")
    def _main(cfg: DictConfig) -> None:
        summary = run(_config_from_omegaconf(cfg))
        log_event(log, "cli.done", bundle=summary["bundle"], n_members=summary["n_members"])

    _main()


if __name__ == "__main__":
    cli()
