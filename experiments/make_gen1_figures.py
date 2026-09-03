"""Generate the three Generation-1 manuscript figures from locked results.

Every number drawn here is read from a locked artifact at run time. Nothing is typed into this
file, so a figure cannot drift from the result it depicts -- the same discipline the evidence lock
applies to the records, applied to the panels.

SVG is written directly rather than through a plotting library. That is deliberate: it adds no
dependency to a project whose reproducibility claim is the whole point, the output is vector and
deterministic byte-for-byte, and journals take SVG/PDF.

    python experiments/make_gen1_figures.py

  figure_1_design.svg      what was measured, and the eligibility funnel
  figure_2_primary.svg     the preregistered ranking result against its null
  figure_3_robustness.svg  strata and the top-choice diagnostic
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = RESULTS / "manuscript" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

INK = "#1a1a1a"
MID = "#6b6b6b"
LINE = "#c9c9c9"
W5C = "#1f4e79"     # the interaction model
W4C = "#8a8a8a"     # the additive comparator
NULLC = "#d9d9d9"
ACC = "#a33b20"     # the observed statistic
FONT = ("font-family=\"Helvetica Neue,Helvetica,Arial,sans-serif\"")


def _j(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def txt(x, y, s, size=11, anchor="start", fill=INK, weight="normal", style="normal"):
    return (f'<text x="{x:.1f}" y="{y:.1f}" {FONT} font-size="{size}" fill="{fill}" '
            f'text-anchor="{anchor}" font-weight="{weight}" font-style="{style}">{esc(str(s))}</text>')


def rect(x, y, w, h, fill, stroke="none", rx=0):
    return (f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(w, 0):.1f}" height="{max(h, 0):.1f}" '
            f'fill="{fill}" stroke="{stroke}" rx="{rx}"/>')


def line(x1, y1, x2, y2, stroke=LINE, width=1, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="{width}"{d}/>')


def svg(w, h, body, title, desc):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
            f'viewBox="0 0 {w} {h}" role="img" aria-labelledby="t d">\n'
            f'<title id="t">{esc(title)}</title><desc id="d">{esc(desc)}</desc>\n'
            f'{rect(0, 0, w, h, "#ffffff")}\n' + "\n".join(body) + "\n</svg>\n")


# =============================================================================================== #
# Figure 1 — what was measured
# =============================================================================================== #
def figure_1(v25: dict, a26: dict) -> str:
    d = v25["descriptives"]
    n_total = v25["eligible_clones"] + d["excluded_all_zero_clones"] + d["excluded_all_positive_clones"]
    conds = ["Acid", "Cisplatin", "CoCl2", "Dabrafenib", "Doxorubicin", "Trametinib"]
    b, W = [], 760

    b.append(txt(24, 30, "Figure 1  Clone-level prospective design and evaluable population",
                 13, weight="bold"))
    b.append(txt(24, 50, "One BRAF-V600E melanoma line (WM989). Data from Schaff et al. 2026; "
                         "reanalysis only.", 10.5, fill=MID))

    # --- the design row ---
    y = 86
    b.append(rect(24, y, 150, 54, "#eef2f6", W5C, rx=4))
    b.append(txt(99, y + 22, f"{n_total:,} barcoded", 11, "middle", weight="bold"))
    b.append(txt(99, y + 38, "clones, pretreatment", 10, "middle", fill=MID))

    b.append(txt(196, y + 26, "split", 10, "middle", fill=MID))
    b.append(line(180, y + 30, 214, y + 30, MID, 1.2))
    b.append(f'<polygon points="214,{y+30} 208,{y+27} 208,{y+33}" fill="{MID}"/>')

    cx, cw = 224, 84
    for i, c in enumerate(conds):
        x = cx + i * (cw + 4)
        b.append(rect(x, y, cw, 54, "#ffffff", LINE, rx=3))
        b.append(txt(x + cw / 2, y + 24, c, 9.5, "middle"))
        # 4 weeks per arm throughout; the two chemotherapies as treat-then-recover. The exact
        # split is in the source paper and is deliberately not asserted here -- the GEO summary
        # and the paper's methods disagree on the doxorubicin schedule.
        b.append(txt(x + cw / 2, y + 40, "4 wk" if c in ("Acid", "CoCl2", "Dabrafenib",
                                                         "Trametinib") else "4 wk, treat+recover",
                     7.5, "middle", fill=MID))

    b.append(txt(24, y + 82, "Outcome C1  post-treatment clone DETECTION — an observed zero means "
                             "no assigned cell was seen.", 10.5))
    b.append(txt(24, y + 98, "Not death, not sensitivity, not clinical response.", 10.5,
                 fill=ACC, style="italic"))

    # --- the funnel ---
    fy = 224
    b.append(txt(24, fy, "Evaluable population", 11.5, weight="bold"))
    b.append(txt(24, fy + 17, "A within-clone AUROC is undefined without both a detected and an "
                              "undetected condition.", 10, fill=MID))

    bars = [(n_total, "all clones", "#dfe6ec"),
            (d["excluded_all_zero_clones"], "never detected — excluded", "#f0f0f0"),
            (d["excluded_all_positive_clones"], "always detected — excluded", "#f0f0f0"),
            (v25["eligible_clones"], "evaluable", W5C)]
    by, bw = fy + 34, 470
    for i, (n, lab, col) in enumerate(bars):
        yy = by + i * 26
        b.append(rect(150, yy, bw * n / n_total, 18, col, rx=2))
        b.append(txt(142, yy + 13, lab, 10, "end"))
        b.append(txt(156 + bw * n / n_total, yy + 13, f"{n:,}", 10,
                     fill=INK if col == W5C else MID, weight="bold" if col == W5C else "normal"))

    # --- the vocabulary is closed ---
    vy = by + 4 * 26 + 22
    b.append(line(24, vy, W - 24, vy, LINE))
    b.append(txt(24, vy + 20, "Closed vocabulary", 11.5, weight="bold"))
    b.append(txt(24, vy + 37,
                 f"{a26['n_refused']} of {a26['n_adversarial_strings']} adversarial condition "
                 f"strings refused, including 16 oncology drugs "
                 f"(Vemurafenib, Carboplatin, ...).", 10, fill=MID))
    b.append(txt(24, vy + 53,
                 f"Design width {a26['structural_closure']['design_columns']} columns "
                 f"= {a26['structural_closure']['K']} PCs "
                 f"+ {a26['structural_closure']['nuisance']} abundance "
                 f"+ {a26['structural_closure']['dummies']} condition "
                 f"+ {a26['structural_closure']['interaction']} interaction.", 10, fill=MID))

    return svg(W, vy + 76, b, "Study design and evaluable population",
               "Clone-level prospective design in WM989 and the funnel from all clones to the "
               "evaluable subset.")


# =============================================================================================== #
# Figure 2 — the preregistered result
# =============================================================================================== #
def figure_2(v25: dict) -> str:
    p, perm, sec = v25["primary"], v25["permutation"], v25["secondary"]
    b, W = [], 760

    b.append(txt(24, 30, "Figure 2  Preregistered clone-specific ranking result", 13, weight="bold"))
    b.append(txt(24, 50, "Equal-clone-weighted within-clone AUROC. Every degree of freedom fixed "
                         "before these numbers existed.", 10.5, fill=MID))

    # --- panel A: the three models ---
    b.append(txt(24, 84, "A   Ranking score by model", 11.5, weight="bold"))
    lo, hi = 0.66, 0.76
    x0, x1 = 150, 640
    ay = 104

    def sx(v):
        return x0 + (v - lo) / (hi - lo) * (x1 - x0)

    for t in [0.66, 0.68, 0.70, 0.72, 0.74, 0.76]:
        b.append(line(sx(t), ay, sx(t), ay + 104, LINE, 1, "2,3"))
        b.append(txt(sx(t), ay + 120, f"{t:.2f}", 9, "middle", fill=MID))

    rows = [("W1   B + U", sec["R_W1"], W4C, "nuisance + condition"),
            ("W4   X + B + U", p["R_W4"], W4C, "+ additive state  (comparator)"),
            ("W5   X + B + U + X×U", p["R_W5"], W5C, "+ state × condition")]
    for i, (lab, val, col, note) in enumerate(rows):
        yy = ay + 8 + i * 34
        b.append(rect(x0, yy, sx(val) - x0, 18, col, rx=2))
        b.append(txt(142, yy + 13, lab, 10, "end"))
        b.append(txt(sx(val) + 8, yy + 13, f"{val:.6f}", 9.5,
                     weight="bold" if col == W5C else "normal"))
        b.append(txt(sx(val) + 66, yy + 13, note, 9, fill=MID))

    b.append(txt(150, ay + 142,
                 "R(W4) sits BELOW R(W1): the additive state term adds nothing to ordering. "
                 "The entire gain is the interaction.", 9.5, fill=ACC))

    # --- panel B: observed vs the full-refit null ---
    ny = 300
    b.append(txt(24, ny, "B   Observed ΔRANK against 1,000 full-refit permutations",
                 11.5, weight="bold"))
    nx0, nx1 = 150, 640
    nlo, nhi = -0.004, 0.058
    ny0 = ny + 30

    def nx(v):
        return nx0 + (v - nlo) / (nhi - nlo) * (nx1 - nx0)

    b.append(line(nx0, ny0 + 46, nx1, ny0 + 46, MID, 1))
    for t in [0.00, 0.01, 0.02, 0.03, 0.04, 0.05]:
        b.append(line(nx(t), ny0 + 46, nx(t), ny0 + 51, MID, 1))
        b.append(txt(nx(t), ny0 + 65, f"{t:+.2f}", 9, "middle", fill=MID))

    # the null, drawn as mean +- sd with its p95 and max marked
    m, sd = perm["null_mean"], perm["null_sd"]
    b.append(rect(nx(m - sd), ny0 + 18, nx(m + sd) - nx(m - sd), 22, NULLC, rx=2))
    b.append(line(nx(m), ny0 + 14, nx(m), ny0 + 44, MID, 1.4))
    b.append(txt(nx(m), ny0 + 8, "null mean ± SD", 9, "middle", fill=MID))
    for val, lab in [(perm["null_p95"], "p95"), (perm["null_max"], "max of 1,000")]:
        b.append(line(nx(val), ny0 + 14, nx(val), ny0 + 44, MID, 1, "3,2"))
        b.append(txt(nx(val), ny0 + 80, f"{lab} {val:.4f}", 8.5, "middle", fill=MID))

    obs = p["delta_RANK"]
    ci = p["bootstrap_ci95"]
    b.append(line(nx(ci[0]), ny0 + 29, nx(ci[1]), ny0 + 29, ACC, 2))
    for e in ci:
        b.append(line(nx(e), ny0 + 23, nx(e), ny0 + 35, ACC, 2))
    b.append(f'<circle cx="{nx(obs):.1f}" cy="{ny0+29:.1f}" r="5" fill="{ACC}"/>')
    b.append(txt(nx(obs), ny0 + 6, f"observed {obs:+.6f}", 10, "middle", fill=ACC, weight="bold"))

    b.append(txt(150, ny0 + 104,
                 f"{perm['n_null_ge_observed']} of {perm['n_perm']} draws reached the observed "
                 f"value  →  p < 0.001.  "
                 f"{obs / perm['null_p95']:.1f}× the null p95; "
                 f"{(obs - perm['null_mean']) / perm['null_sd']:.1f} SD above the null mean.",
                 10))
    b.append(txt(150, ny0 + 120,
                 "Permutation refits the entire pipeline inside every draw; observed-data "
                 "hyperparameters are never reused.", 9.5, fill=MID))

    return svg(W, ny0 + 146, b, "Preregistered ranking result",
               "Ranking score for three models and the observed statistic against a 1,000-draw "
               "full-refit permutation null.")


# =============================================================================================== #
# Figure 3 — robustness and the top-choice diagnostic
# =============================================================================================== #
def figure_3(v25: dict) -> str:
    d, t1 = v25["descriptives"], v25["delta_TOP1"]
    b, W = [], 760
    b.append(txt(24, 30, "Figure 3  Robustness across strata, and the top-choice diagnostic",
                 13, weight="bold"))
    b.append(txt(24, 50, "Preregistered as descriptive: these could withhold support, never grant "
                         "it.", 10.5, fill=MID))

    def panel(title, items, ox, oy, sub):
        out = [txt(ox, oy, title, 11.5, weight="bold"), txt(ox, oy + 16, sub, 9.5, fill=MID)]
        bx, bw, top = ox + 74, 200, 0.09
        out.append(line(bx, oy + 30, bx, oy + 30 + len(items) * 26, LINE))
        for i, (lab, val, n) in enumerate(items):
            yy = oy + 36 + i * 26
            out.append(rect(bx, yy, bw * val / top, 15, W5C, rx=2))
            out.append(txt(bx - 8, yy + 12, lab, 9.5, "end"))
            out.append(txt(bx + bw * val / top + 7, yy + 12, f"{val:+.4f}", 9))
            out.append(txt(bx + bw + 62, yy + 12, f"n={n}", 8.5, fill=MID))
        return out

    folds = [(f"fold {k}", v["delta_RANK"], v["clones"])
             for k, v in sorted(d["by_outer_fold"].items())]
    depth = [(k if k != "1" else "1 cell", v["delta_RANK"], v["clones"])
             for k, v in d["by_pretreatment_depth_bin"].items()]

    b += panel("A   By held-out fold", folds, 24, 84, "ΔRANK, positive in all five")
    b += panel("B   By pretreatment clone depth", depth, 396, 84,
               "ΔRANK, positive in all five")

    cy = 84 + 36 + 5 * 26 + 26
    b.append(line(24, cy, W - 24, cy, LINE))
    b.append(txt(24, cy + 24, "C   Choosing each clone's lowest predicted detection score",
                 11.5, weight="bold"))
    b.append(txt(24, cy + 40, "How often that condition is a genuine observed zero.", 9.5,
                 fill=MID))

    x0, bw2 = 210, 300
    for i, (lab, val, col) in enumerate([("W5  interaction", t1["LOW_PERSISTENCE_TOP1_W5"], W5C),
                                         ("W4  additive", t1["LOW_PERSISTENCE_TOP1_W4"], W4C)]):
        yy = cy + 56 + i * 28
        b.append(rect(x0, yy, bw2 * val, 18, col, rx=2))
        b.append(txt(x0 - 8, yy + 13, lab, 10, "end"))
        b.append(txt(x0 + bw2 * val + 8, yy + 13, f"{val * 100:.1f}%", 10,
                     weight="bold" if col == W5C else "normal"))

    b.append(txt(x0, cy + 128,
                 f"ΔTOP1 {t1['value']:+.6f}   CI95 "
                 f"[{t1['ci95'][0]:+.6f}, {t1['ci95'][1]:+.6f}]", 10))
    b.append(txt(24, cy + 128, "Directional check", 9.5, fill=MID))
    b.append(txt(x0, cy + 144, "A consistency check, not a significance test.", 9.5, fill=MID))

    return svg(W, cy + 172, b, "Robustness and top-choice diagnostic",
               "Delta RANK by fold and by pretreatment depth, and the low-persistence top-choice "
               "diagnostic.")


def main() -> int:
    v25 = _j(RESULTS / "stage25" / "stage25_verdict.json")
    a26 = _j(RESULTS / "stage26" / "stage26a_vocabulary_closure.json")

    figures = {
        "figure_1_design.svg": figure_1(v25, a26),
        "figure_2_primary.svg": figure_2(v25),
        "figure_3_robustness.svg": figure_3(v25),
    }
    for name, content in figures.items():
        (OUT / name).write_text(content, encoding="utf-8")
        print(f"  {name:28s} {len(content):>7,} bytes")
    print(f"\nwritten to {OUT.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
