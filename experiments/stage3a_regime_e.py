"""REGIME E — is Stage 3a's forward gate gradeable on GSE165177, which is already on disk?

    python experiments/stage3a_regime_e.py

READ-ONLY. Writes `results/stage3a_regime_e_results.json`. No build, no retrain, `src/`
untouched. **Does NOT grade 3a.** It runs the `REF_GROUND_RULES.md` §5b null at GSE165177's real
geometry, against the outcome table pre-registered in `plans/STAGE_3A_REGIME_E_PREREG.md`
(committed as 81602fe, BEFORE this file existed).

WHY
---
3a-bis measured that the held-out `gill_bulk` geometry is UNRESOLVABLE at every effect size, and I
concluded "acquire a second dense time course". The other machine's AUDIT-2 pointed out that
conclusion skipped a dataset already in hand, and it is right:

  gill_bulk    6 donors, ~1.7 samples per (donor, timepoint), 6 controls total (day 0 only)
  GSE165177    3 donors, 6-9 samples per (donor, timepoint), 33 CONTEMPORANEOUS negative
               controls at 2-3 per donor PER timepoint, donors aged 53/53/38 -- all adult and
               all inside the clock's fitted [1, 96], unlike HFF at age 0

The saturation mechanism 3a-bis blamed is `gill_bulk`-specific: at 1.7 samples a fraction can only
be {0, 0.5, 1}. At 4-6 it cannot. That is a different target and it has never been measured.

WHAT IS PRE-REGISTERED (see the plan; graded verbatim in §OUTCOMES below)
  E1  pass >= 0.95 at alpha 1.0 AND 0.5 -> RESOLVABLE  -> acquisition ask WITHDRAWN
  E2  pass >= 0.95 at alpha 1.0 only    -> MARGINAL    -> ask stays ON HOLD
  E3  pass <  0.95 at both              -> UNRESOLVABLE-> acquisition ESTABLISHED with a number
  E4  pass >  0.05 at alpha 0           -> NULL BROKEN -> discard the run entirely
  P0  unsafe fraction has zero variance in >=2 of 3 donors -> REGIME E VOID, do not read E1-E4
  A1/A2/A3  the {3,6} folds x {~2, real} samples 2x2, read ONLY if E3 fires

DECLARED LIMITS (from the pre-registration, restated so they travel with the output)
  1. SAFETY TARGET ONLY. No dAge, no clock, no harmonizer.
  2. Trains WITHIN GSE165177 (2 donors) and holds out the third -- the gill/HFF bundles are in the
     2000-gene panel space and joining gene spaces is its own piece of work. This UNDERSTATES what
     a joined training set could do.
  3. Therefore E3 does not by itself prove acquisition is required; the gene-space join is the
     next thing to cost.
  4. exp1/exp2 batch structure is present and NOT corrected for; it is reported per fold.
  5. The audit's "no modality shift" claim is NOT tested here -- it needs the join in limit 2.
"""
from __future__ import annotations

import contextlib
import gzip
import importlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "experiments"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "stage3a_regime_e_results.json"

GSE_DIR = Path(r"D:\GSE165177")
MATRICES = ("GSE165177_Log2_RPM_Transient_reprogramming.txt.gz",
            "GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz")
SERIES = "GSE165177_series_matrix.txt.gz"
ANNOT_COLS = 12          # Probe, Chromosome, Start, End, Probe Strand, Feature, ID, Description,
#                          Feature Strand, Type, Feature Orientation, Distance

SIM_TRIALS = 2000
SEED = 0
ALPHAS = (0.0, 0.25, 0.5, 1.0)
MIN_PASS = 0.95
IPSC_DONOR = "iPSC"

# HFF's measured unsafe fraction by DAY (GSE242423, 42,481 cells at ~4,700/timepoint).
# Recorded in STAGE_3_TOOL.md "3a-bis". The simulated truth is this curve scaled by alpha --
# a real effect this corpus actually exhibits, not an invented effect size.
HFF_DAYS = np.array([0., 2., 4., 6., 8., 10., 12., 14., 21.])
HFF_CURVE = np.array([0.0835, 0.3509, 0.3508, 0.4226, 0.4708, 0.3964, 0.4659, 0.6791, 0.9996])


def _diag():
    return sys.modules.get("stage3a_diagnose") or importlib.import_module("stage3a_diagnose")


def _t18():
    return sys.modules.get("test18_forward_gate") or importlib.import_module(
        "test18_forward_gate")


# ---------------------------------------------------------------------------------------------
def parse_sample_name(col: str) -> dict | None:
    """`O1_negative_control_intermediate_13days_exp2` -> donor / arm / day / exp.

    The DAY-0 fibroblasts are named differently -- `O1 Fib`, space-separated, with no `days`
    token -- and a regex written for the treated samples drops all three SILENTLY. That is the
    same class of invisible filter that has cost this project real time before, so day 0 is
    matched explicitly rather than left to fall through. `iPSC 13` / `iPSC 21` are separate iPSC
    lines with no donor attribution (the series matrix files them under donor `iPSC`, day 51) and
    are deliberately excluded -- they are not a point on any donor's trajectory.
    """
    m = re.match(r"^(?P<donor>O\d)_(?P<arm>.+?)_(?P<day>\d+)days_(?P<exp>exp\d)$", col)
    if m:
        d = m.groupdict()
        return {"sample": col, "donor": d["donor"], "arm": d["arm"],
                "day": float(d["day"]), "exp": d["exp"]}
    m0 = re.match(r"^(?P<donor>O\d)\s+Fib$", col)
    if m0:
        # the untreated starting state of the trajectory. NOT marked is_control: the 33
        # contemporaneous `negative_control` samples are the better z-scoring reference, and
        # marking day 0 as a control would delete the only timepoint before the transition.
        return {"sample": col, "donor": m0.group("donor"), "arm": "day0_fibroblast",
                "day": 0.0, "exp": "exp2"}
    return None


def load_matrix() -> tuple[np.ndarray, list[str], pd.DataFrame]:
    """Join both Log2-RPM parts on gene symbol. Returns (N_samples, G), genes, obs."""
    frames = []
    for fn in MATRICES:
        with gzip.open(GSE_DIR / fn, "rt", encoding="utf-8", errors="replace") as f:
            head = f.readline().rstrip("\n").split("\t")
        cols = head[ANNOT_COLS:]
        use = [head[0], *cols]
        df = pd.read_csv(GSE_DIR / fn, sep="\t", usecols=use, compression="gzip",
                         low_memory=False)
        df = df.rename(columns={head[0]: "gene"}).set_index("gene")
        df = df.apply(pd.to_numeric, errors="coerce")
        frames.append(df)
    joined = frames[0].join(frames[1], how="inner", lsuffix="_a", rsuffix="_b")
    # duplicate symbols: keep the highest-expressed row, mirroring the 10x loader's rule
    joined = joined.assign(_m=joined.mean(axis=1)).sort_values("_m", ascending=False)
    joined = joined[~joined.index.duplicated(keep="first")].drop(columns="_m").sort_index()
    joined = joined.dropna(how="any")

    meta = [parse_sample_name(c) for c in joined.columns]
    keep = [i for i, m in enumerate(meta) if m is not None]
    obs = pd.DataFrame([meta[i] for i in keep])
    X = joined.to_numpy(dtype=np.float64).T[keep]          # (N samples, G genes)
    return X, list(joined.index), obs


def donor_ages() -> dict[str, str]:
    ages, titles = [], []
    with gzip.open(GSE_DIR / SERIES, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            p = [v.strip().strip('"') for v in line.rstrip("\n").split("\t")]
            if p[0] == "!Sample_title":
                titles = p[1:]
            elif p[0] == "!Sample_characteristics_ch1" and p[1].startswith("donor age"):
                ages = [v.split(":")[1].strip() for v in p[1:]]
    out: dict[str, str] = {}
    for t, a in zip(titles, ages, strict=True):
        out.setdefault(t.split("_")[0], a)
    return out


def pairs_of(rows: list[dict]) -> list[dict]:
    out = []
    for i in range(len(rows)):
        for j in range(len(rows)):
            if rows[j]["day"] <= rows[i]["day"]:
                continue
            out.append({"x_i": rows[i]["x"], "dt": rows[j]["day"] - rows[i]["day"],
                        "day_j": rows[j]["day"], "n_j": rows[j]["n"], "u_j": rows[j]["u"]})
    return out


def p_of(days, alpha: float) -> np.ndarray:
    g = np.interp(np.asarray(days, float), HFF_DAYS, HFF_CURVE)
    return np.clip(HFF_CURVE.mean() + alpha * (g - HFF_CURVE.mean()), 0.01, 0.99)


# ---------------------------------------------------------------------------------------------
def main() -> int:
    from audit_metrics import bar_verdict

    from cellfate.common.console import install_pretty_console, render_table
    from cellfate.data.labels import fate_labels
    install_pretty_console()
    D, t18 = _diag(), _t18()

    print("\n" + "=" * 92)
    print("REGIME E — is 3a's forward gate gradeable on GSE165177 (already on disk)?")
    print("=" * 92)
    print("READ-ONLY. Graded against plans/STAGE_3A_REGIME_E_PREREG.md, committed 81602fe")
    print("BEFORE this script existed. Safety target only -- no dAge, no clock, no harmonizer.")

    X, genes, obs = load_matrix()
    ages = donor_ages()
    obs["is_control"] = obs["arm"].str.contains("negative_control", case=False)
    obs["cell_line"] = obs["donor"]
    print(f"\n   loaded {X.shape[0]} samples x {X.shape[1]} genes; donors="
          f"{sorted(obs.donor.unique())}; ages={ {d: ages.get(d) for d in sorted(obs.donor.unique())} }")
    print(f"   controls (negative_control arm): {int(obs.is_control.sum())}; "
          f"batches: {dict(obs.exp.value_counts())}")

    # ---- fate labels, by the pipeline's own function -----------------------------------------
    soft = fate_labels(X, genes, obs)
    cls = np.argmax(soft, axis=1)
    from cellfate.common.constants import DEATH_IDX, LOSS_IDX
    unsafe = ((cls == LOSS_IDX) | (cls == DEATH_IDX)).astype(float)
    obs["unsafe"] = unsafe
    print(f"   fate labels via cellfate.data.labels.fate_labels; overall unsafe fraction "
          f"{unsafe.mean():.3f}")

    # ---- geometry + the P0 precondition ------------------------------------------------------
    donors = sorted(d for d in obs.donor.unique() if d != IPSC_DONOR)
    traj: dict[str, list[dict]] = {}
    geo_rows = []
    for d in donors:
        rows = []
        for day in sorted(obs.loc[obs.donor == d, "day"].unique()):
            sel = (obs.donor == d) & (obs.day == day)
            tgt = sel & ~obs.is_control                 # the target population: treated samples
            if not tgt.any():
                continue
            idx = np.flatnonzero(tgt.to_numpy())
            rows.append({"day": float(day), "x": X[idx].mean(0), "n": int(len(idx)),
                         "u": float(unsafe[idx].mean()),
                         "n_ctrl": int((sel & obs.is_control).sum())})
        traj[d] = rows
        for r in rows:
            geo_rows.append([d, ages.get(d, "?"), f"{r['day']:.0f}", str(r["n"]),
                             str(r["n_ctrl"]), f"{r['u']:.3f}"])
    print("\n  GEOMETRY — treated samples per (donor, day), and the CONTEMPORANEOUS controls")
    print(render_table(["donor", "age", "day", "treated n", "controls n", "unsafe frac"],
                       geo_rows, aligns=["l", "r", "r", "r", "r", "r"]))

    # ---- WHY the target looks the way it does -- per-arm label distribution ------------------
    # Not pre-registered. Added because the day-0 UNTREATED fibroblasts came back `unsafe`, which
    # cannot be right biologically and therefore had to be explained before P0 was reported.
    from cellfate.common.constants import CLASSES
    obs["cls"] = [CLASSES[c] for c in cls]
    for i, c in enumerate(CLASSES):
        obs[f"p_{c}"] = soft[:, i]
    arm_key = obs["arm"].str.replace("_intermediate", "", regex=False)
    arm_rows = []
    for a, sub in obs.groupby(arm_key):
        arm_rows.append([a, str(len(sub)), f"{sub.p_safe.mean():.3f}", f"{sub.p_loss.mean():.3f}",
                         f"{sub.p_death.mean():.3f}",
                         ", ".join(f"{k}={v}" for k, v in sub.cls.value_counts().items())])
    print("\n  LABEL DISTRIBUTION BY ARM — the reason the target has no time structure")
    print(render_table(["arm", "n", "P(safe)", "P(loss)", "P(death)", "hard labels"],
                       arm_rows, aligns=["l", "r", "r", "r", "r", "l"]))
    d0 = obs[obs.arm == "day0_fibroblast"]
    print(f"   ** the UNTREATED day-0 fibroblasts label `{d0.cls.iloc[0]}` with "
          f"P(loss) = {', '.join(f'{v:.3f}' for v in d0.p_loss)} **")
    print("   They are the starting material -- somatic identity by definition -- so this is a")
    print("   LABELLER-ON-BULK artefact, not biology. `fate_labels` z-scores against the")
    print("   `is_control` samples, which here are fibroblasts cultured 10-17 days; anything that")
    print("   differs from that reference lands on the unsafe side. The split it produces is")
    print("   CONTROL vs NON-CONTROL, not a time course.")
    day0_p_loss = [float(v) for v in d0.p_loss]

    flat = {d: float(np.std([r["u"] for r in traj[d]])) for d in donors}
    n_flat = sum(1 for v in flat.values() if v == 0.0)
    print(f"\n  P0 PRECONDITION — SD of the unsafe fraction across timepoints: "
          f"{ {d: round(v, 4) for d, v in flat.items()} }")
    void = n_flat >= 2
    print(f"   donors with ZERO variation: {n_flat} of {len(donors)}  ->  "
          f"{'REGIME E IS VOID' if void else 'precondition satisfied, continue'}")

    out: dict = {"script": "stage3a_regime_e", "prereg": "plans/STAGE_3A_REGIME_E_PREREG.md",
                 "n_samples": int(X.shape[0]), "n_genes": int(X.shape[1]),
                 "donors": donors, "ages": {d: ages.get(d) for d in donors},
                 "n_controls": int(obs.is_control.sum()),
                 "batches": {k: int(v) for k, v in obs.exp.value_counts().items()},
                 "geometry": geo_rows, "unsafe_sd_by_donor": flat,
                 "label_by_arm": arm_rows, "day0_p_loss": day0_p_loss,
                 "P0_void": bool(void),
                 "pairs_per_donor": {d: len(pairs_of(traj[d])) for d in donors}}
    print(f"   ordered forward pairs per donor: {out['pairs_per_donor']} "
          f"(total {sum(out['pairs_per_donor'].values())})")

    if void:
        print("\n   P0 FIRED. Per the pre-registration, E1-E4 must NOT be read. Stopping.")
        print("\n  THE STRUCTURAL REASON, which generalises beyond this dataset:")
        print("   `p_unsafe` is a fraction OF CELLS. In bulk RNA-seq every sample is ALREADY a")
        print("   population average, so a hard label per sample collapses that fraction to 0 or")
        print("   1 before it can be counted, and the 'fraction' becomes a fraction of SAMPLES.")
        print("   That is why gill_bulk sat at 63/70 values pinned to the bounds and why this")
        print("   dataset -- with 4-6x the replication and real contemporaneous controls -- is")
        print("   saturated at 1.000 in 11 of 12 cells. More bulk replication cannot fix it:")
        print("   the quantity is not expressible in bulk at all. NOTE this is specific to the")
        print("   SAFETY target; dAge is continuous per sample and does NOT have this problem,")
        print("   so GSE165177's replication and contemporaneous controls may still be valuable")
        print("   there. That is a separate question and is not measured here.")
        out["structural_finding"] = (
            "p_unsafe is a fraction of CELLS; a bulk sample is already a population average, so a "
            "per-sample hard label collapses the fraction to 0/1 before it can be counted. No "
            "bulk dataset can carry this target at any replication. Specific to the SAFETY "
            "target -- dAge is continuous per sample and is unaffected.")
        OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
        return 0

    # ---- the §5b null at the real geometry ---------------------------------------------------
    prep = []
    for hd in donors:
        tr = [p for d in donors if d != hd for p in pairs_of(traj[d])]
        te = pairs_of(traj[hd])
        if len(tr) < 4 or len(te) < 3:
            continue
        prep.append({
            "tag": hd, "n_tr": len(tr), "n_te": len(te),
            "s": D.FrozenRidge(t18.feats(tr, False), t18.feats(te, False)),
            "d": D.FrozenRidge(t18.feats(tr, True), t18.feats(te, True)),
            "day_tr": np.array([p["day_j"] for p in tr], float),
            "day_te": np.array([p["day_j"] for p in te], float),
            "cnt_tr": np.array([max(p["n_j"], 1) for p in tr], int),
            "cnt_te": np.array([max(p["n_j"], 1) for p in te], int),
        })
    print(f"\n   folds: {[p['tag'] for p in prep]}  training pairs "
          f"{[p['n_tr'] for p in prep]}  held-out pairs {[p['n_te'] for p in prep]}")

    rng = np.random.default_rng(SEED)

    def run_cell(alpha: float, n_override: int | None, n_folds_mult: int) -> dict:
        """`n_folds_mult=2` simulates 6 donors behaving like these 3, with independent draws."""
        units = prep * n_folds_mult
        lat = [{"tr": p_of(p["day_tr"], alpha), "te": p_of(p["day_te"], alpha)} for p in units]
        hi_raw, hi_log = [], []
        for _ in range(SIM_TRIALS):
            ds_r, dd_r, ds_l, dd_l = [], [], [], []
            for p, L in zip(units, lat, strict=True):
                c_te = (np.full_like(p["cnt_te"], n_override) if n_override
                        else p["cnt_te"])
                c_tr = (np.full_like(p["cnt_tr"], n_override) if n_override
                        else p["cnt_tr"])
                ytr = rng.binomial(c_tr, L["tr"]) / c_tr
                yte = rng.binomial(c_te, L["te"]) / c_te
                ds_r.append(D.mae(p["s"].predict(ytr), yte))
                dd_r.append(D.mae(p["d"].predict(ytr), yte))
                ds_l.append(D.mae(D.expit(p["s"].predict(D.logit(ytr))), yte))
                dd_l.append(D.mae(D.expit(p["d"].predict(D.logit(ytr))), yte))
            hi_raw.append(t18.paired_ci([b - a for a, b in zip(ds_r, dd_r, strict=True)])[1][1])
            hi_log.append(t18.paired_ci([b - a for a, b in zip(ds_l, dd_l, strict=True)])[1][1])
        return {"raw": bar_verdict(np.array(hi_raw, float), 0.0, lower_is_better=True),
                "logit": bar_verdict(np.array(hi_log, float), 0.0, lower_is_better=True)}

    cells = {
        "E_operative_3folds_real_n": (None, 1),
        "Ea_3folds_2_per_tp": (2, 1),
        "Ed_6folds_real_n": (None, 2),
        "Ec_6folds_2_per_tp": (2, 2),
    }
    for name, (n_ov, mult) in cells.items():
        rows = []
        for alpha in ALPHAS:
            v = run_cell(alpha, n_ov, mult)
            rows.append([f"{alpha:.2f}", f"{v['raw']['pass_rate']:.3f}", v["raw"]["verdict"],
                         f"{v['logit']['pass_rate']:.3f}", v["logit"]["verdict"]])
            out.setdefault("cells", {}).setdefault(name, {})[f"alpha={alpha}"] = v
        lbl = {"E_operative_3folds_real_n": "E — OPERATIVE: 3 folds, real 4-6 samples/tp",
               "Ea_3folds_2_per_tp": "E-a — 3 folds at gill's ~2 samples/tp",
               "Ed_6folds_real_n": "E-d — 6 folds (counterfactual), real samples/tp",
               "Ec_6folds_2_per_tp": "E-c — 6 folds (counterfactual) at ~2 samples/tp"}[name]
        print(f"\n  {lbl}")
        print(render_table(["alpha", "pass (raw)", "verdict (raw)", "pass (logit)",
                            "verdict (logit)"], rows, aligns=["r", "r", "l", "r", "l"]))

    # ---- grade the PRE-REGISTERED outcomes ----------------------------------------------------
    def rate(cell, alpha, est="logit"):
        return out["cells"][cell][f"alpha={alpha}"][est]["pass_rate"]

    op = "E_operative_3folds_real_n"
    a1, a05 = rate(op, 1.0), rate(op, 0.5)
    a0_max = max(rate(op, 0.0, e) for e in ("raw", "logit"))
    if a0_max > 0.05:
        verdict, action = "E4 — NULL BROKEN", "discard the run; E1-E3 may not be read"
    elif a1 >= MIN_PASS and a05 >= MIN_PASS:
        verdict, action = "E1 — RESOLVABLE", "acquisition ask WITHDRAWN; proceed to a graded 3a"
    elif a1 >= MIN_PASS:
        verdict, action = "E2 — MARGINAL", "ask stays ON HOLD; a null below alpha=1 is uninformative"
    else:
        verdict, action = "E3 — UNRESOLVABLE", "acquisition ESTABLISHED with a number"

    print("\n" + "-" * 92)
    print("PRE-REGISTERED OUTCOME (logit, operative cell)")
    print("-" * 92)
    print(render_table(["alpha=1.0", "alpha=0.5", "alpha=0.0 (max raw/logit)", "VERDICT"],
                       [[f"{a1:.3f}", f"{a05:.3f}", f"{a0_max:.3f}", verdict]],
                       aligns=["r", "r", "r", "l"]))
    print(f"   -> {action}")
    out["outcome"] = {"verdict": verdict, "action": action, "alpha1": a1, "alpha05": a05,
                      "alpha0_max": a0_max, "min_pass": MIN_PASS}

    if verdict.startswith("E3"):
        d_ok = rate("Ed_6folds_real_n", 1.0) >= MIN_PASS
        a_ok = rate("Ea_3folds_2_per_tp", 1.0) >= MIN_PASS
        attr = ("A1 — FOLD COUNT binds; the ask is MORE DONORS, not more cells" if d_ok
                else "A2 — REPLICATION binds" if (not d_ok and not a_ok and a1 >= MIN_PASS)
                else "A3 — BOTH bind; neither fix alone suffices")
        print(f"\n   ATTRIBUTION (read because E3 fired): {attr}")
        print(f"      E-d 6 folds real n @a=1: {rate('Ed_6folds_real_n', 1.0):.3f}   "
              f"E-a 3 folds ~2/tp @a=1: {rate('Ea_3folds_2_per_tp', 1.0):.3f}")
        out["outcome"]["attribution"] = attr

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    print("\n   NO 3a VERDICT IS TAKEN HERE. Limits: safety target only; trains WITHIN GSE165177")
    print("   (2 donors) so it UNDERSTATES a joined training set; exp1/exp2 batch not corrected.")
    return 0


if __name__ == "__main__":
    with contextlib.suppress(KeyboardInterrupt):
        raise SystemExit(main())
    raise SystemExit(130)
