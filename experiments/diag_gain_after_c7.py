"""Does §4.6's "the sparse clock INVERTS on HFF" survive C-7?

    python experiments/diag_gain_after_c7.py "D:/Gill" "D:/GSE242423"

READ-ONLY. Writes `results/diag_gain_after_c7_results.json`. `src/` untouched, no build touched.

WHY
---
`STAGE_1_5_6_SPARSE_CLOCK.md` §4.5 measured the harmonization gain on HFF at **2.152** (dense
clock) and §4.6 at **2.769** (top-100), concluding that *"§1's finding does not transfer to HFF.
It INVERTS"* — sparsifying concentrates the clock onto exactly the genes harmonization amplifies
most, so the sparse clock improves Gill and degrades HFF.

**Both numbers were computed with `N2_Fib_Sendai_Exp2` inside `sigma_gill`.** §5.7 showed that
column is a near-constant vector, and §5.14 showed removing it moves HFF's day-14 ΔAge from
−26.755 to −8.196 — an implied residual gain of ~0.82 rather than 2.152. So the gain figures the
inversion rests on are contaminated by the same defect C-7 removes, and §5.14 recorded that as an
*implication* without measuring it. This measures it.

WHAT IS COMPUTED
----------------
Gain is the ratio the harmonizer's own algebra defines (§4.3):

    d_harmonized / d_direct  =  SUM_g delta_g * (sigma_gill,g / sigma_hff,g) * w_g
                               ---------------------------------------------------
                                        SUM_g delta_g * w_g

with `delta_g` = HFF's day-14 mean minus its day-0 mean (log1p-CP10k), on the pipeline's own gene
space, and `sigma` floored exactly as `harmonize.py:112-113` does.

Four cells: {dense clock, top-100 clock} x {sigma_gill from ALL SIX controls (as 1c/1d had it),
sigma_gill with N2_Fib EXCLUDED (what C-7 ships)}.

THE QUESTION, STATED SO IT CAN FAIL
------------------------------------
    INVERSION SURVIVES   gain(top100) > gain(dense) with the contaminant removed
                         -> §4.6 stands; sparsifying really does degrade HFF
    INVERSION IS AN ARTEFACT  gain(top100) <= gain(dense) once cleaned
                         -> §4.6's conclusion was contamination, and the option table it drives
                            (§4.6's four ways forward) was answering a question that did not exist

Either way this changes NO decision about the clock: step 2's bar fails on skin & blood at EVERY
k (verified over the recorded sweep: best sb sign-agreement 0.676 against a 0.80 bar), and that
comparison never touches harmonization or HFF. This is about whether a recorded FINDING is right.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "diag_gain_after_c7_results.json"

SHIPPED = ROOT / "runs" / "cellfate_multi" / "harmonization.json"
EPS = 1e-8
DAY14, DAY0 = 14.0, 0.0
TOP_K = 100
CONTAMINANT = "N2"          # donor whose _Fib_ control is the degenerate column


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def sigma_floored(mat: np.ndarray) -> np.ndarray:
    """`harmonize.py:111-113` -- per-gene std over controls, floored at its own median."""
    sig = mat.std(axis=0)
    return np.maximum(sig, float(np.median(sig)))


def main() -> int:
    from cellfate.common.console import install_pretty_console, render_table
    from cellfate.data.normalize import normalize_counts
    from cellfate.data.qc import QCConfig, apply_qc
    from cellfate.data.sources import GSE242423SingleCellSource
    install_pretty_console()

    gill_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("D:/Gill")
    hff_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("D:/GSE242423")

    ship = json.loads(SHIPPED.read_text("utf-8"))
    genes_G = list(ship["genes"])
    sig_hff = np.asarray(ship["stats"]["hff_sc"]["sigma"], float)

    clock = json.loads((ROOT / "configs" / "clocks" / "fleischer_clock.json").read_text("utf-8"))
    W = {k: float(v) for k, v in clock["weights"].items()}
    w_dense = np.array([W.get(g, 0.0) for g in genes_G])
    # top-K by |weight| over the FULL clock, then restricted to the harmonizer's space -- which
    # is what "the pipeline as it actually runs" means (§4.6).
    top = sorted(W, key=lambda g: -abs(W[g]))[:TOP_K]
    topset = set(top)
    w_top = np.array([W.get(g, 0.0) if g in topset else 0.0 for g in genes_G])

    print("\n" + "=" * 78)
    print("DOES §4.6's 'the sparse clock INVERTS on HFF' SURVIVE C-7?")
    print("=" * 78)
    print(f"gene space {len(genes_G)} | clock genes present {int((w_dense != 0).sum())} | "
          f"top-{TOP_K} present {int((w_top != 0).sum())}")

    # ---- Gill controls, with and without the contaminant -------------------------------- #
    dv = _load("dv", ROOT / "experiments" / "diag_dage_variants.py")
    samples, genes, lin = dv.load_gill(gill_dir)
    gnorm = np.asarray(normalize_counts(lin), float)
    idx = {g: i for i, g in enumerate(genes)}
    cols = np.array([idx[g] for g in genes_G])
    ctrl = {s.split("_")[0]: gnorm[i, cols] for i, s in enumerate(samples)
            if re.search(r"_d\d+_", s) is None and "_Fib_" in s}
    donors = sorted(ctrl)
    print(f"Gill controls: {donors}")

    sig_all = sigma_floored(np.stack([ctrl[d] for d in donors]))
    sig_clean = sigma_floored(np.stack([ctrl[d] for d in donors if d != CONTAMINANT]))

    # ---- HFF delta ---------------------------------------------------------------------- #
    hli = _load("hli", ROOT / "experiments" / "diag_hff_label_identity.py")
    gfile = next(hff_dir.glob("*genes.tsv.gz"))
    smp = []
    for mtx in sorted(hff_dir.glob("*.matrix.mtx.gz")):
        bc = mtx.with_name(mtx.name.replace(".matrix.mtx.gz", ".barcodes.tsv.gz"))
        if bc.exists():
            smp.append({"matrix": str(mtx), "barcodes": str(bc),
                        "label": mtx.name.split(".")[0].split("_")[-1]})
    src = GSE242423SingleCellSource(smp, str(gfile), min_genes=hli.MIN_GENES,
                                    max_cells_per_sample=hli.MAX_CELLS,
                                    cells_per_run=hli.CELLS_PER_RUN)
    qc = QCConfig(max_mito_frac=0.20, min_genes=hli.MIN_GENES)
    s14 = np.zeros(len(genes_G))
    s0 = np.zeros(len(genes_G))
    c14 = c0 = 0
    for chunk in src.plan():
        raw = apply_qc(src.fetch(chunk), qc)
        if len(raw.obs) == 0:
            continue
        day = raw.obs["time_h"].to_numpy(float) / 24.0
        m14, m0 = np.isclose(day, DAY14), np.isclose(day, DAY0)
        if not (m14.any() or m0.any()):
            continue
        nm = np.asarray(normalize_counts(raw.counts), float)
        ix = {g: i for i, g in enumerate(raw.genes)}
        cc = np.array([ix.get(g, -1) for g in genes_G])
        take = cc >= 0
        if m14.any():
            s14[take] += nm[m14][:, cc[take]].sum(axis=0)
            c14 += int(m14.sum())
        if m0.any():
            s0[take] += nm[m0][:, cc[take]].sum(axis=0)
            c0 += int(m0.sum())
    delta = s14 / c14 - s0 / c0
    print(f"HFF: {c14} day-14 cells, {c0} day-0 cells")

    def gain(w: np.ndarray, sig_gill: np.ndarray) -> tuple[float, float, float]:
        direct = float((delta * w).sum())
        harm = float((delta * (sig_gill / (sig_hff + EPS)) * w).sum())
        return direct, harm, (harm / direct if direct else float("nan"))

    rows, out = [], {}
    for wname, w in (("dense", w_dense), (f"top{TOP_K}", w_top)):
        for sname, sg in (("all 6 controls (as 1c/1d)", sig_all),
                          ("N2_Fib EXCLUDED (C-7)", sig_clean)):
            d, h, g = gain(w, sg)
            rows.append([wname, sname, f"{d:+.3f}", f"{h:+.3f}", f"{g:.3f}"])
            out[f"{wname}|{sname}"] = {"direct": d, "harmonized": h, "gain": g}

    print("\n" + render_table(["clock", "sigma_gill from", "direct ΔAge", "harmonized", "gain"],
                              rows, aligns=["l", "l", "r", "r", "r"]))

    g_dense_c = out["dense|N2_Fib EXCLUDED (C-7)"]["gain"]
    g_top_c = out[f"top{TOP_K}|N2_Fib EXCLUDED (C-7)"]["gain"]
    g_dense_d = out["dense|all 6 controls (as 1c/1d)"]["gain"]
    g_top_d = out[f"top{TOP_K}|all 6 controls (as 1c/1d)"]["gain"]

    survives = bool(g_top_c > g_dense_c)
    print(f"\n  CONTAMINATED : gain(top{TOP_K}) {g_top_d:.3f}  vs  gain(dense) {g_dense_d:.3f}"
          f"   -> {'top100 higher (the recorded inversion)' if g_top_d > g_dense_d else 'no inversion'}")
    print(f"  CLEAN (C-7)  : gain(top{TOP_K}) {g_top_c:.3f}  vs  gain(dense) {g_dense_c:.3f}"
          f"   -> {'top100 STILL higher' if survives else 'inversion GONE'}")
    print(f"\n  VERDICT: §4.6's inversion "
          f"{'SURVIVES C-7' if survives else 'IS AN ARTEFACT of the degenerate control'}")
    print("\n  Either way this changes no decision about the clock: step 2's bar fails on")
    print("  skin & blood at EVERY k (best sb sign-agreement 0.676 vs a 0.80 bar), and that")
    print("  comparison never touches harmonization or HFF.")

    out["verdict"] = "INVERSION_SURVIVES" if survives else "INVERSION_IS_ARTEFACT"
    out["recorded_1c_1d"] = {"dense_gain": 2.152, "top100_gain": 2.769}
    out["n_hff"] = {"day14": c14, "day0": c0}
    out["n_genes"] = {"space": len(genes_G), "clock_present": int((w_dense != 0).sum()),
                      "top_present": int((w_top != 0).sum())}
    OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
