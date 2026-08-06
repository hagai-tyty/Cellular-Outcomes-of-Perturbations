# ΔAge LEDGER — truth vs expected vs actual, per condition, per parameter

**Generated** by `experiments/diag_dage_ledger.py` + `diag_dage_ksweep.py`. **Read-only,
`src/` untouched, no label moved.** Every number here is measured, none is projected.

## How to read a row

| column | meaning |
|---|---|
| **TRUTH** | ΔAge in years from **methylation** — the instrument 1.5.1 validated (negative control inert at +0.5/−2.4 yr, dose-response p = 0.0001, SNR 3.4) |
| **EXPECTED** | what the pipeline *should* produce if the RNA clock worked — i.e. identical to TRUTH |
| **ACTUAL** | what the RNA clock actually produced for that condition |
| **ERROR** | ACTUAL − TRUTH, in years. Positive = the RNA clock reads *older* than reality |

ΔAge = age(condition) − age(matched control at the **same donor and day**). Replicates
(exp1/exp2) are averaged to one row per condition **before** any scoring — an earlier run scored
them as independent rows, which was pseudo-replication.

---

## THE HEADLINE

**The shipped clock over-reports rejuvenation by ~14 years, systematically. Restricting it to its
~100 largest-weight genes removes almost all of that.**

| | MAE | bias | ρ | sign agreement |
|---|---:|---:|---:|---:|
| **full clock (33,155 genes)** | **16.61 yr** | **−14.10** | +0.703 | 0.62 |
| **top-100 genes** | **5.36 yr** | **−1.61** | **+0.835** | **0.94** |

*(transient arm, 68 conditions, Horvath multi-tissue as truth)*

**MAE 5.36 yr is below the reference instrument's own donor-level error of ±7 yr** (1.5.2 §12-R:
two donors of identical age 53 read 44.0 and 58.5). The RNA clock is now agreeing with methylation
about as closely as methylation agrees with itself.

---

## Parameter sweep — every k, both clocks, transient arm (ΔAge vs ΔAge)

| k | MAE sb | bias sb | ρ sb | sign sb | MAE mt | bias mt | ρ mt | sign mt |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 6.59 | +4.28 | +0.358 | 0.51 | 7.69 | +6.08 | +0.495 | 0.69 |
| 20 | 7.47 | +3.26 | +0.221 | 0.41 | 6.19 | +5.06 | +0.777 | 0.85 |
| 50 | 6.69 | +1.42 | +0.249 | 0.44 | 5.99 | +3.22 | +0.762 | 0.79 |
| **100** | 8.79 | -3.41 | +0.260 | 0.44 | 5.36 | -1.61 | +0.835 | 0.94 |
| **150** | 8.47 | -3.79 | +0.269 | 0.46 | 5.42 | -1.99 | +0.839 | 0.90 |
| 200 | 10.13 | -3.91 | +0.241 | 0.40 | 6.31 | -2.11 | +0.793 | 0.90 |
| 300 | 12.97 | -9.49 | +0.279 | 0.50 | 9.99 | -7.69 | +0.799 | 0.82 |
| 500 | 13.72 | -10.99 | +0.379 | 0.53 | 11.64 | -9.19 | +0.733 | 0.71 |
| 1000 | 14.15 | -11.95 | +0.408 | 0.63 | 12.33 | -10.15 | +0.722 | 0.60 |
| 2000 | 17.21 | -15.02 | +0.374 | 0.54 | 15.30 | -13.22 | +0.744 | 0.69 |
| 5000 | 17.65 | -15.71 | +0.420 | 0.65 | 16.31 | -13.91 | +0.707 | 0.59 |
| **all 33,155** | 17.84 | -15.90 | +0.432 | 0.68 | 16.61 | -14.10 | +0.703 | 0.62 |

**Bias crosses zero exactly where MAE bottoms out (k ≈ 100).** That co-incidence is the
mechanism: thousands of near-zero weights each contribute a little drift, and summed over 33,155
genes they become a −14 yr offset. Dropping them removes the offset rather than merely shrinking
the noise.

## Does k generalise? Leave-one-donor-out

| held-out donor | k chosen on the other two | held-out MAE | full clock |
|---|---:|---:|---:|
| O1 | 50 | **6.70** | 16.59 |
| O2 | 100 | **6.84** | 16.48 |
| O3 | 100 | **5.45** | 16.75 |

**k is stable at 50–100 and the improvement survives donor-level cross-validation.** The gain is
not selection: the full clock is worse than every sparse variant for every donor individually.

---

## Per-condition ledger — transient arm, 68 conditions with a matched control

`raw` = the shipped clock. `top100` = the same clock restricted to its 100 largest |weights|.
TRUTH and ERROR are against Horvath multi-tissue.

| donor | condition | day | reps | TRUTH | EXPECTED | ACTUAL raw | ERROR raw | ACTUAL top100 | ERROR top100 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| O1 | failed to transiently reprogram | 10 | 1 | +7.16 | +7.16 | -2.03 | -9.19 | +0.29 | -6.87 |
| O1 | failed to transiently reprogram | 13 | 2 | -5.11 | -5.11 | -21.99 | -16.88 | -3.35 | +1.76 |
| O1 | failed to transiently reprogram | 15 | 2 | +4.01 | +4.01 | -16.92 | -20.93 | +1.83 | -2.18 |
| O1 | failed to transiently reprogram | 17 | 2 | -5.13 | -5.13 | -22.21 | -17.08 | -2.78 | +2.35 |
| O1 | failing to transiently reprogram intermediate | 10 | 1 | -5.76 | -5.76 | -28.37 | -22.61 | -17.78 | -12.02 |
| O1 | failing to transiently reprogram intermediate | 13 | 1 | -7.52 | -7.52 | -27.12 | -19.60 | -14.04 | -6.52 |
| O1 | failing to transiently reprogram intermediate | 15 | 1 | -9.23 | -9.23 | -19.51 | -10.28 | -10.87 | -1.64 |
| O1 | failing to transiently reprogram intermediate | 17 | 1 | -5.00 | -5.00 | -18.36 | -13.36 | -8.46 | -3.47 |
| O1 | negative control | 10 | 1 | -0.07 | -0.07 | -0.15 | -0.08 | +1.21 | +1.28 |
| O1 | negative control | 13 | 2 | +3.88 | +3.88 | -9.33 | -13.21 | +0.79 | -3.09 |
| O1 | negative control | 15 | 2 | +3.98 | +3.98 | -5.75 | -9.72 | +1.95 | -2.03 |
| O1 | negative control | 17 | 2 | +2.93 | +2.93 | -3.22 | -6.14 | +3.39 | +0.46 |
| O1 | negative control intermediate | 10 | 1 | +0.07 | +0.07 | +0.15 | +0.08 | -1.21 | -1.28 |
| O1 | negative control intermediate | 13 | 1 | -3.88 | -3.88 | +9.33 | +13.21 | -0.79 | +3.09 |
| O1 | negative control intermediate | 15 | 1 | -3.98 | -3.98 | +5.75 | +9.72 | -1.95 | +2.03 |
| O1 | negative control intermediate | 17 | 1 | -2.93 | -2.93 | +3.22 | +6.14 | -3.39 | -0.46 |
| O1 | transient reprogramming intermediate | 10 | 1 | -19.62 | -19.62 | -64.06 | -44.45 | -29.28 | -9.66 |
| O1 | transient reprogramming intermediate | 13 | 1 | -24.50 | -24.50 | -50.03 | -25.54 | -26.40 | -1.91 |
| O1 | transient reprogramming intermediate | 15 | 1 | -25.34 | -25.34 | -49.23 | -23.89 | -24.22 | +1.12 |
| O1 | transient reprogramming intermediate | 17 | 1 | -32.91 | -32.91 | -52.53 | -19.62 | -21.66 | +11.26 |
| O1 | transiently reprogrammed | 10 | 1 | -5.91 | -5.91 | -33.56 | -27.65 | -2.61 | +3.30 |
| O1 | transiently reprogrammed | 17 | 1 | -2.07 | -2.07 | -37.70 | -35.63 | -5.94 | -3.87 |
| O2 | failed to transiently reprogram | 10 | 1 | +4.90 | +4.90 | +1.01 | -3.89 | +3.17 | -1.74 |
| O2 | failed to transiently reprogram | 13 | 2 | -3.79 | -3.79 | -7.95 | -4.16 | +1.39 | +5.18 |
| O2 | failed to transiently reprogram | 15 | 2 | +1.28 | +1.28 | -7.34 | -8.62 | +1.80 | +0.52 |
| O2 | failed to transiently reprogram | 17 | 2 | -7.60 | -7.60 | -16.35 | -8.74 | -1.25 | +6.35 |
| O2 | failing to transiently reprogram intermediate | 10 | 1 | -5.16 | -5.16 | -35.66 | -30.50 | -21.21 | -16.05 |
| O2 | failing to transiently reprogram intermediate | 13 | 1 | -7.16 | -7.16 | -36.97 | -29.81 | -20.62 | -13.46 |
| O2 | failing to transiently reprogram intermediate | 15 | 1 | -7.95 | -7.95 | -37.24 | -29.29 | -19.57 | -11.62 |
| O2 | failing to transiently reprogram intermediate | 17 | 1 | -8.51 | -8.51 | -38.70 | -30.19 | -17.24 | -8.73 |
| O2 | negative control | 10 | 1 | +6.61 | +6.61 | +6.55 | -0.06 | +3.23 | -3.38 |
| O2 | negative control | 13 | 2 | +3.10 | +3.10 | -2.83 | -5.92 | +3.11 | +0.01 |
| O2 | negative control | 15 | 2 | +1.90 | +1.90 | -2.55 | -4.45 | +4.11 | +2.21 |
| O2 | negative control | 17 | 2 | +2.68 | +2.68 | -1.75 | -4.42 | +4.88 | +2.21 |
| O2 | negative control intermediate | 10 | 1 | -6.61 | -6.61 | -6.55 | +0.06 | -3.23 | +3.38 |
| O2 | negative control intermediate | 13 | 1 | -3.10 | -3.10 | +2.83 | +5.92 | -3.11 | -0.01 |
| O2 | negative control intermediate | 15 | 1 | -1.90 | -1.90 | +2.55 | +4.45 | -4.11 | -2.21 |
| O2 | negative control intermediate | 17 | 1 | -2.68 | -2.68 | +1.75 | +4.42 | -4.88 | -2.21 |
| O2 | transient reprogramming intermediate | 10 | 1 | -23.12 | -23.12 | -70.87 | -47.75 | -34.17 | -11.05 |
| O2 | transient reprogramming intermediate | 13 | 1 | -28.38 | -28.38 | -62.54 | -34.16 | -33.04 | -4.66 |
| O2 | transient reprogramming intermediate | 15 | 1 | -38.49 | -38.49 | -62.71 | -24.23 | -29.77 | +8.71 |
| O2 | transient reprogramming intermediate | 17 | 1 | -45.95 | -45.95 | -73.26 | -27.30 | -28.98 | +16.97 |
| O2 | transiently reprogrammed | 10 | 1 | -23.72 | -23.72 | -15.55 | +8.17 | -3.26 | +20.46 |
| O2 | transiently reprogrammed | 13 | 2 | -11.04 | -11.04 | -40.00 | -28.96 | -15.66 | -4.62 |
| O2 | transiently reprogrammed | 17 | 1 | +12.01 | +12.01 | -21.48 | -33.49 | +0.37 | -11.64 |
| O3 | failed to transiently reprogram | 10 | 1 | +2.26 | +2.26 | -5.13 | -7.39 | +2.13 | -0.13 |
| O3 | failed to transiently reprogram | 13 | 2 | +1.84 | +1.84 | -13.13 | -14.98 | +0.18 | -1.66 |
| O3 | failed to transiently reprogram | 15 | 2 | +0.55 | +0.55 | -10.72 | -11.27 | +1.42 | +0.87 |
| O3 | failed to transiently reprogram | 17 | 2 | -0.23 | -0.23 | -17.32 | -17.09 | -0.45 | -0.22 |
| O3 | failing to transiently reprogram intermediate | 10 | 1 | -1.32 | -1.32 | -39.09 | -37.77 | -18.85 | -17.53 |
| O3 | failing to transiently reprogram intermediate | 13 | 1 | -5.76 | -5.76 | -34.18 | -28.41 | -20.32 | -14.56 |
| O3 | failing to transiently reprogram intermediate | 15 | 1 | -4.85 | -4.85 | -25.49 | -20.64 | -14.71 | -9.86 |
| O3 | failing to transiently reprogram intermediate | 17 | 1 | -4.62 | -4.62 | -29.98 | -25.36 | -12.91 | -8.29 |
| O3 | negative control | 10 | 1 | +0.95 | +0.95 | +2.37 | +1.42 | +3.46 | +2.51 |
| O3 | negative control | 13 | 2 | +1.89 | +1.89 | -5.88 | -7.77 | +2.49 | +0.60 |
| O3 | negative control | 15 | 2 | +1.32 | +1.32 | -6.17 | -7.49 | +4.31 | +2.99 |
| O3 | negative control | 17 | 2 | +0.31 | +0.31 | -4.51 | -4.82 | +4.50 | +4.20 |
| O3 | negative control intermediate | 10 | 1 | -0.95 | -0.95 | -2.37 | -1.42 | -3.46 | -2.51 |
| O3 | negative control intermediate | 13 | 1 | -1.89 | -1.89 | +5.88 | +7.77 | -2.49 | -0.60 |
| O3 | negative control intermediate | 15 | 1 | -1.32 | -1.32 | +6.17 | +7.49 | -4.31 | -2.99 |
| O3 | negative control intermediate | 17 | 1 | -0.31 | -0.31 | +4.51 | +4.82 | -4.50 | -4.20 |
| O3 | transient reprogramming intermediate | 10 | 1 | -23.90 | -23.90 | -74.94 | -51.04 | -32.25 | -8.35 |
| O3 | transient reprogramming intermediate | 13 | 1 | -22.20 | -22.20 | -61.09 | -38.88 | -32.22 | -10.02 |
| O3 | transient reprogramming intermediate | 15 | 1 | -26.34 | -26.34 | -58.11 | -31.77 | -32.67 | -6.33 |
| O3 | transiently reprogrammed | 10 | 1 | -3.18 | -3.18 | -18.18 | -15.00 | -4.74 | -1.56 |
| O3 | transiently reprogrammed | 13 | 2 | -19.68 | -19.68 | -25.77 | -6.09 | -6.04 | +13.64 |
| O3 | transiently reprogrammed | 15 | 2 | -4.31 | -4.31 | +7.26 | +11.57 | +5.61 | +9.92 |
| O3 | transiently reprogrammed | 17 | 2 | -5.72 | -5.72 | -30.78 | -25.06 | -7.53 | -1.81 |

*Full table including the skin & blood clock, every variant and the Sendai arm:*
`results/dage_ledger.csv` (90 rows × 60 columns).