"""
END-TO-END SMOKE TEST for Stage 1, on synthetic data shaped like the REAL failure.

    python smoke_stage1.py

WHY THIS EXISTS. Stage 1 run 1 burned 3.5 h of GPU and produced a void experiment, because
`cell_line` merged one bulk corpus (HFF, 33,613 cells) with six tiny donors (~14 cells each) and
the inner-LODO rotated over the corpus as if it were a donor. The existing test-suite fixtures use
BALANCED synthetic sources (2 lines x equal cells), so none of them could ever have reproduced
that geometry -- the bug was invisible to the tests by construction.

This builds a dataset with the SAME SHAPE as the real one:

    BULK_L0     ~300 cells   <- the corpus, must be SKIPPED by the inner-LODO
    DONOR_L0..5  ~20 cells each  <- the donors, must be the ONLY things rotated over

...then runs build -> train -> calibrate -> bundle -> predict and asserts every property Stage 1
is supposed to guarantee. It runs on CPU in a couple of minutes, so a defect costs minutes here
instead of hours on the data machine.

WHAT IT WOULD HAVE CAUGHT: the bulk-corpus rotation (check 4), a silent fallback to
in-distribution calibration (3), an uncalibrated inference mode (6), and a lopsided residual pool
(5). Check 8 verifies the claim that the mc_dropout passes do not disturb training reproducibility
-- dropout consumes RNG, and that claim was argued rather than measured.

Exit code is 0 only if every check passes.
"""
from __future__ import annotations

import os
import shutil
import sys
import tempfile
import time
import traceback
from pathlib import Path

# Before torch initialises CUDA (mirrors retrain_stage1.py). This runs on CPU anyway, but the
# smoke test should exercise the same environment the real run does.
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

import numpy as np  # noqa: E402

BULK = "BULK_L0"
N_DONORS = 6
RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    RESULTS.append((name, ok, detail))
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    return ok


def build_dataset(root: Path):
    """One dominant corpus + six small donors: the geometry that broke run 1."""
    from cellfate.data import DataConfig, QCConfig, SyntheticSource
    from cellfate.data import run as build

    cfg = DataConfig(
        out=str(root), gene_panel=str(root / "panel.json"), n_genes=64,
        qc=QCConfig(min_genes=1, max_mito_frac=0.9), label_tau=0.5,
        split_fracs=(0.6, 0.2, 0.1, 0.1), primary_regime="scaffold", seed=0,
    )
    return build(cfg, sources=[
        # cells per line = n_cells * (2 controls + 2 doses * n_compounds)
        SyntheticSource(name="bulk", n_lines=1, n_compounds=4,          # ~300 cells
                        n_cells_per_condition=30, n_filler_genes=40, seed=1),
        SyntheticSource(name="donor", n_lines=N_DONORS, n_compounds=4,  # ~20 cells each
                        n_cells_per_condition=2, n_filler_genes=40, seed=2),
    ])


def train(root: Path, seed: int = 0):
    from cellfate.training.train_model import TrainConfig
    from cellfate.training.train_model import run as train_run

    return train_run(TrainConfig(
        dataset_dir=str(root), out=str(root), regime="scaffold",
        d_cell=32, d_u=16, latent_dim=32, p_drop=0.2,        # p_drop>0: mc_dropout needs spread
        lr=3e-3, epochs=6, patience=6, batch_size=128,
        ensemble_size=2, base_seed=seed, conformal_levels=(0.90,), device="cpu",
        xdonor_calibration=True, inference_mode="ensemble", mc_dropout_T=16,
    ))


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass

    print("\nSTAGE 1 END-TO-END SMOKE TEST")
    print("synthetic data shaped like the real failure: one bulk corpus + six small donors")

    from cellfate.common.io import ArtifactPaths
    from cellfate.common.scalers import Scalers
    from cellfate.training.dataset import DONOR_I, DONOR_VOCAB, X_I, load_split_tensors
    from cellfate.training.xdonor_calib import MIN_INNER_TRAIN_FRAC

    tmp = Path(tempfile.mkdtemp(prefix="cellfate_smoke_"))
    t0 = time.time()
    try:
        # ---- 1. BUILD -------------------------------------------------------------- #
        print("\n  [1/8] building synthetic dataset ...")
        summary = build_dataset(tmp)
        paths = ArtifactPaths.of(tmp)
        n_shards = len(list(paths.shards_dir.glob("*.parquet")))
        check("dataset builds", n_shards > 0,
              f"{n_shards} shards, {summary.get('n_samples', '?')} cells")

        # ---- 2. DONOR COLUMN (sub-stage 1a) ---------------------------------------- #
        print("\n  [2/8] donor column ...")
        import torch
        sc = Scalers.load(paths.scalers_file)
        ds = load_split_tensors(paths, sc, "scaffold", "train")
        d = ds.tensors[DONOR_I]
        counts = {k: int((d == v).sum()) for k, v in DONOR_VOCAB.items()}
        counts = {k: n for k, n in counts.items() if n}

        check("7 tensor columns", len(ds.tensors) == 7, f"got {len(ds.tensors)}")
        check("donor codes are integer", d.dtype == torch.long, str(d.dtype))
        check("donor column length matches X", len(d) == len(ds.tensors[X_I]))
        check("empty split also returns 7",
              len(load_split_tensors(paths, sc, "scaffold", "__none__").tensors) == 7)
        check("all 7 cell lines present in train", len(counts) == N_DONORS + 1,
              f"{counts}")

        # ---- 3. THE GEOMETRY IS ACTUALLY LOPSIDED ---------------------------------- #
        n_train = len(ds)
        bulk_n = counts.get(BULK, 0)
        check("bulk corpus dominates train (the run-1 geometry)",
              bulk_n > MIN_INNER_TRAIN_FRAC * n_train,
              f"{BULK}={bulk_n}/{n_train} = {100 * bulk_n / max(n_train, 1):.0f}%")

        # ---- 4. TRAIN + CALIBRATE --------------------------------------------------- #
        print(f"\n  [3/8] train + inner-LODO calibrate ({n_train} train cells) ...")
        m = train(tmp)

        check("cross-donor calibration ran", m.get("xdonor_calibrated") is True)
        # THE CENTRAL CHECK: the corpus must be skipped, the donors must not be
        check("bulk corpus SKIPPED by inner-LODO",
              m.get("xdonor_n_donors") == N_DONORS,
              f"rotated over {m.get('xdonor_n_donors')} donors, expected {N_DONORS} "
              f"({N_DONORS + 1} lines minus the corpus)")

        # ---- 5. POOL COMPOSITION ---------------------------------------------------- #
        print("\n  [4/8] residual pool composition ...")
        per = {int(k): int(v) for k, v in (m.get("xdonor_residuals_per_donor") or {}).items()}
        total = sum(per.values())
        top = max(per.values()) if per else 0
        check("residual pool is not dominated by one donor",
              bool(per) and top <= 0.5 * total,
              f"largest donor = {top}/{total} = {100 * top / max(total, 1):.0f}%")
        bulk_code = DONOR_VOCAB.get(BULK)
        check("the corpus contributed NO residuals", bulk_code not in per,
              f"pool = {per}")

        # ---- 6. BOTH INFERENCE MODES CALIBRATED ------------------------------------- #
        print("\n  [5/8] per-mode sigma calibration ...")
        from cellfate.common import io as _io
        s_ens, s_mc = m.get("sigma_scale"), m.get("sigma_scale_mc")
        modes = list(_io.load_conformal(paths).sigma_calibrated_modes)

        # NOT "> 1.0": the factor is clamped at 1.0, so on well-fit synthetic data an already-
        # adequate spread legitimately yields exactly 1.0. Calibration STATUS is the invariant;
        # the magnitude is data-dependent and belongs in the report, not an assertion.
        check("both modes recorded as calibrated",
              sorted(modes) == ["ensemble", "mc_dropout"], f"{modes}")
        check("ensemble factor is finite and >= 1",
              isinstance(s_ens, float) and np.isfinite(s_ens) and s_ens >= 1.0, f"{s_ens}")
        check("mc_dropout factor is finite and >= 1",
              isinstance(s_mc, float) and np.isfinite(s_mc) and s_mc >= 1.0, f"{s_mc}")
        if s_ens > 1.0 and s_mc > 1.0:
            check("the two factors are distinct (not borrowed)", s_ens != s_mc,
                  f"ensemble={s_ens:.3f} mc={s_mc:.3f}")
        else:
            print(f"      note: a factor clamped to 1.0 (ensemble={s_ens:.3f}, mc={s_mc:.3f}) "
                  "-- spread already adequate on this synthetic fit, not a failure")

        # ---- 6b. THE FATE CALIBRATOR + THE PERSISTED POOL ---------------------------- #
        print("\n  [5b/8] fate calibration and the persisted cross-donor pool ...")
        from cellfate.training.xdonor_calib import XSTATS_FILENAME, load_xstats

        a, b = m.get("platt_a"), m.get("platt_b")
        check("Platt fitted on P(safe), not a multi-class temperature",
              isinstance(a, float) and isinstance(b, float), f"a={a} b={b}")
        check("Platt slope is positive (monotone: cannot REORDER cells, so the rank guards hold)",
              isinstance(a, float) and a > 0, f"a={a}")
        check("temperature left at 1.0 (one calibrator, not two stacked)",
              m.get("temperature") == 1.0, str(m.get("temperature")))
        # the graded quantity is reported alongside the top-1 figure
        pre, post = m.get("xdonor_safe_ece_before"), m.get("xdonor_safe_ece_after")
        check("binary P(safe) ECE is reported (the metric scorecard grades)",
              isinstance(pre, float) and isinstance(post, float), f"{pre} -> {post}")
        if isinstance(pre, float) and isinstance(post, float):
            print(f"      binary P(safe) ECE on the cross-donor pool: {pre:.3f} -> {post:.3f}")

        # the strict cross-donor variant: fitted and reported, never shipped, so Stage 1's
        # principle is TESTED on every run rather than quietly dropped
        n = m.get("fate_calib_n") or {}
        check("fit uses all held-out cells, not just the pool",
              n.get("total", 0) > n.get("xdonor", 0),
              f"total={n.get('total')} (in_dist={n.get('in_dist')} + xdonor={n.get('xdonor')})")
        check("strict cross-donor variant reported as a diagnostic",
              isinstance(m.get("xdonor_only_platt_a"), float),
              f"a={m.get('xdonor_only_platt_a')} on n={m.get('xdonor_only_n')} "
              f"from {m.get('xdonor_only_n_donors')} donors")
        shipped, pool_only = m.get("shipped_safe_ece_on_pool"), m.get(
            "xdonor_only_safe_ece_insample")
        if isinstance(shipped, float) and isinstance(pool_only, float):
            print(f"      on the cross-donor pool: shipped(all data) {shipped:.3f} vs "
                  f"pool-only {pool_only:.3f} -- the latter is IN-SAMPLE, so it flatters itself")

        xs_path = ArtifactPaths.of(tmp).bundle_dir / XSTATS_FILENAME
        check("cross-donor pool persisted", xs_path.exists(), str(xs_path.name))
        if xs_path.exists():
            back = load_xstats(ArtifactPaths.of(tmp).bundle_dir)
            check("persisted pool round-trips with the right row count",
                  back.abs_residuals.size == m.get("xdonor_n_residuals")
                  and back.probs_mean.shape[0] == m.get("xdonor_n_residuals"),
                  f"{back.abs_residuals.size} rows vs {m.get('xdonor_n_residuals')}")

        # ---- 7. PREDICTOR ROUND-TRIP ------------------------------------------------ #
        print("\n  [6/8] Predictor round-trip ...")
        from cellfate.common.errors import ConfigError
        from cellfate.inference import Predictor

        pe = Predictor(tmp, mode="ensemble")
        pm = Predictor(tmp, mode="mc_dropout", T=16)
        check("ensemble Predictor picks its own factor",
              abs(pe.sigma_scale - s_ens) < 1e-9, f"{pe.sigma_scale}")
        check("mc_dropout Predictor picks its own factor",
              abs(pm.sigma_scale - s_mc) < 1e-9, f"{pm.sigma_scale}")

        # a bundle NOT CALIBRATED for this mode must refuse, not serve raw sigma
        from cellfate.common import io
        copy = Path(tempfile.mkdtemp(prefix="cellfate_smoke_copy_"))
        shutil.rmtree(copy)
        shutil.copytree(tmp, copy)
        cp = ArtifactPaths.of(copy)
        io.save_conformal(cp, io.load_conformal(cp).model_copy(
            update={"sigma_calibrated_modes": ["ensemble"]}))
        try:
            Predictor(copy, mode="mc_dropout")
            check("refuses a mode it was never calibrated for", False, "loaded anyway")
        except ConfigError:
            check("refuses a mode it was never calibrated for", True)
        # ...but a mode that WAS calibrated still loads, even if its factor clamped to 1.0
        try:
            Predictor(copy, mode="ensemble")
            check("still serves a mode that WAS calibrated", True)
        except ConfigError as exc:
            check("still serves a mode that WAS calibrated", False, repr(exc)[:80])
        shutil.rmtree(copy, ignore_errors=True)

        # sigma must actually be WIDENED at inference, not just stored
        arr = io.shard_to_numpy(io.read_shard(sorted(paths.shards_dir.glob("*.parquet"))[0]))
        rows = pe.predict_encoded(arr["X"][:8], arr["u_chem_fp"][:8], arr["dose_time"][:8])
        check("predictions are finite",
              all(np.isfinite([r["mu_age"], r["sigma_age"], r["S"]]).all() for r in rows))

        # ---- 8. BACK-COMPAT + TEMPERATURE GUARD ------------------------------------- #
        print("\n  [7/8] back-compat and degenerate-input guards ...")
        legacy = Path(tempfile.mkdtemp(prefix="cellfate_smoke_legacy_"))
        lp = ArtifactPaths.of(legacy)
        lp.bundle_dir.mkdir(parents=True, exist_ok=True)
        io.write_json(lp.bundle_conformal_file, {"levels": [0.9], "q": {"0.9": 8.86}})
        lc = io.load_conformal(lp)
        check("pre-Stage-1b bundles still load",
              lc.q["0.9"] == 8.86 and lc.sigma_scale == 1.0 and lc.sigma_scale_mc == 1.0)
        shutil.rmtree(legacy, ignore_errors=True)

        from cellfate.training import fit_temperature, has_class_variation
        one_class = np.tile([1.0, 0.0, 0.0], (200, 1))
        check("single-class targets detected", not has_class_variation(one_class))
        import warnings as _w
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            t_deg = fit_temperature(np.random.default_rng(0).normal(size=(200, 3)),
                                    one_class).temperature
        check("degenerate temperature is NOT shipped", t_deg == 1.0,
              f"T={t_deg} (a collapse to the 0.01 bound would saturate every probability)")

        # ---- 9. REPRODUCIBILITY (the claim that dropout doesn't disturb training) --- #
        print("\n  [8/8] reproducibility with mc_dropout RNG in the loop ...")
        tmp2 = Path(tempfile.mkdtemp(prefix="cellfate_smoke_rep_"))
        shutil.rmtree(tmp2)
        shutil.copytree(tmp, tmp2)
        for p in (ArtifactPaths.of(tmp2).bundle_dir,):
            shutil.rmtree(p, ignore_errors=True)
        m2 = train(tmp2)
        check("re-training reproduces sigma_scale exactly",
              m2["sigma_scale"] == m["sigma_scale"],
              f"{m['sigma_scale']} vs {m2['sigma_scale']}")
        check("re-training reproduces conformal q exactly",
              m2["conformal_q"] == m["conformal_q"],
              f"{m['conformal_q']} vs {m2['conformal_q']}")
        check("re-training reproduces temperature exactly",
              m2["temperature"] == m["temperature"],
              f"{m['temperature']} vs {m2['temperature']}")
        shutil.rmtree(tmp2, ignore_errors=True)

    except Exception:  # noqa: BLE001
        print("\n   EXCEPTION during the smoke test:\n")
        traceback.print_exc()
        RESULTS.append(("smoke test completed without raising", False, "see traceback"))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # ---- report ------------------------------------------------------------------- #
    failed = [r for r in RESULTS if not r[1]]
    print(f"\n{'=' * 78}")
    print(f"  {len(RESULTS) - len(failed)}/{len(RESULTS)} checks passed "
          f"in {time.time() - t0:.0f}s")
    if failed:
        print("\n  FAILED:")
        for name, _, detail in failed:
            print(f"    - {name}" + (f"  ({detail})" if detail else ""))
        print("\n  Stage 1 is NOT ready to run on the real data. Fix these first --")
        print("  each one costs minutes here and hours on the data machine.")
    else:
        print("\n  All Stage 1 invariants hold on the run-1 geometry.")
        print("  This does NOT prove the calibration is CORRECT -- only that the")
        print("  machinery does what it claims. The bars in the lab notebook decide that.")
    print(f"{'=' * 78}\n")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
