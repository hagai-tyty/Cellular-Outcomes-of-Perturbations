# Pre-registration: three tests, in order

Written **before** any of the three runs. Decision rules are fixed here so that a result cannot be
re-read after the fact. Nothing in `src/` moves; all three are read-only.

---

## Phase 0 — WHICH CLOCK, and why it comes first

The three tests were framed as "redo everything on top100". That framing hides an assumption worth
testing on its own, and testing it changes what the other phases should use.

**What top100 is:** the 100 largest-|weight| genes of the same Fleischer clock
(`diag_dage_variants.py` docstring), applied as `Σ w_i x_i + b`. Not a different instrument — a
truncation of the same one.

**The mechanism to test.** A dense ridge over 33,155 genes fit from 133 samples carries many small
weights fit to noise. In ABSOLUTE age those largely cancel: the fit was optimised for exactly that
target on exactly that cohort. In a DIFFERENCE between two samples they do **not** cancel — the two
noise terms are independent and compound. ΔAge is a difference. So the prediction is:

> **dense wins on absolute age; sparse wins on differences.**

If true, "use top100 everywhere" is wrong — the right rule is *dense for absolute, sparse for
differences*, and each phase below should use whichever it actually needs.

**The test.** GSE113957 has 133 donors with declared ages, so it has **8,778 pairs with a KNOWN age
difference**. That is ground truth for exactly the quantity ΔAge is, at a scale nothing else here
offers, and it needs no methylation.

For each variant (raw, top100, top500, top2000, covnorm, ranknorm) compute:

| statistic | what it measures |
|---|---|
| `MAE_abs` | error against the donor's true age |
| `MAE_diff` | error of `pred_i - pred_j` against `age_i - age_j`, over all pairs |
| `r_diff` | correlation of predicted vs true age difference |
| `MAE_diff / (sqrt2 * MAE_abs)` | **do errors cancel or compound in a difference?** <1 cancel, >1 compound |

**PRE-REGISTERED READING**
- **MECHANISM CONFIRMED** if the variant ranking by `MAE_diff` differs from the ranking by
  `MAE_abs`, AND a sparse variant ranks better on `MAE_diff` than `raw` does.
- **MECHANISM REFUTED** if the two rankings agree — then sparsity is not doing what is claimed, and
  top100's instrument-floor win needs another explanation.

**Caveats stated now.** GSE113957 is the clock's training cohort, so absolute performance there is
optimistic for `raw` by construction — which makes a `raw` loss on differences *stronger* evidence,
not weaker. Pair differences are not independent (each donor appears in 132 pairs), so pair
statistics are descriptive; no p-value is computed from them.

**What Phase 0 decides:** which age readout Phases 1 and 3 use. Not negotiable after seeing them.

---

## Phase 1 — Redo the ΔAge forward and residual tests on the Phase-0 winner

Two earlier results used the dense clock, whose ΔAge is measured at MAE 22.69 against a
methylation floor of 7.30 — worse than predicting zero. They may have been measuring a corrupted
target.

**But they failed for DIFFERENT reasons, and only one of them is a noise story.** This distinction
is recorded before the run because it determines what a change would mean:

| earlier result | why it failed | can a better clock change it? |
|---|---|---|
| early→late, partial **-0.064** after donor age | a **CONFOUND** — both ends track donor age | Only via attenuation. Measurement error biases a correlation toward zero, so a true effect could be masked. But a better clock also measures donor age more cleanly, which does not help. **Weak prior for change.** |
| residual-expression, **3 of 9** robustness | **FRAGILITY** — plausibly noise | Noise is exactly what a better instrument removes. **Genuine prior for change.** |

**PRE-REGISTERED READING**
- Re-run both, unchanged in every other respect, on the Phase-0 winner.
- **CHANGED** if a verdict flips (partial clears its df=3 bar; or robustness reaches 6 of 9).
- **UNCHANGED** otherwise — and then the earlier conclusions stand on a better instrument, which
  makes them stronger rather than merely repeated.
- Report the partial's magnitude either way. A move from -0.064 to, say, +0.3 that still misses the
  bar is informative and must not be reported as "unchanged".

---

## Phase 2 — Transfer, using the project's own Harmonizer

The earlier transfer test used raw features and a crude per-cohort z-score, and skipped
`cellfate.data.harmonize.Harmonizer`, which exists for exactly this. Raw transfer gave MAE 70-118
with negative correlation — a scale failure, not a statement about age being unreadable.

**Design.** `Harmonizer.fit({"gse113957": [...], "gill": [...]}, ref_dataset=...)` on the two
cohorts' CONTROL populations (GSE113957 is entirely untreated fibroblasts; the held-out cohorts are
day-0 fibroblasts, so both sides are controls). Transform both, train on GSE113957, predict the
held-out cohorts' donor ages. Compare raw / zscore / harmonized on identical rows.

**An honest reframing of "transductive".** I previously called the z-score correction
non-deployable because it needs the whole test cohort. That is true of ANY batch alignment,
harmonizer included — you cannot estimate a batch effect from one sample. In practice samples are
processed in batches, so per-batch alignment is a **deployment constraint, not a defect**. What
would be a defect is needing the test *labels*, which none of these do.

**PRE-REGISTERED READING**
- **TRANSFER WORKS** if harmonized MAE <= 20 yr AND Spearman >= 0.6 on the Gill cohort (the only
  one with usable range), for a majority of the alpha grid.
- **DOES NOT** otherwise.
- GSE165177 is reported but **excluded from the verdict**: its donors span 38-53 yr against a ~12 yr
  instrument error, so it has no dynamic range to rank. That exclusion is declared here, before the
  numbers, not after.

---

## Phase 3 — The forward question, on the best design available

Earlier forward work used the 6 Sendai donors, where donor age explains the outcome and cannot be
separated at n=6. **GSE165177 has a better design that was never used for this question:** every
(donor, arm) is its own trajectory over days 10/13/15/17, and within a donor the arms differ in
OUTCOME while donor age is constant by construction.

That is up to **3 donors x 6 arms = 18 trajectories** with donor age held fixed — the confound that
killed the Sendai analysis cannot arise within a donor.

**The forward question:** from a trajectory's EARLY timepoint(s), predict its LATER ΔAge.

**A precondition that must be checked first, and can kill the phase.** The window is d10→d17, only
7 days. If ΔAge does not move materially within that window, there is nothing to predict and the
phase ends there — that will be reported as the result, not worked around.

**PRE-REGISTERED READING**
- Unit is the (donor, arm) trajectory; inference clusters on DONOR (3), not on trajectory (18),
  because arms within a donor share material.
- **SIGNAL** if leave-one-donor-out prediction of late ΔAge from early expression beats a
  permutation null at the 95th percentile for a majority of the alpha grid.
- **NULL** otherwise.
- With 3 donor clusters this is severely underpowered, and that is stated before the run rather
  than discovered after. A NULL here bounds what the design can show; it does not establish absence.

---

## What none of the three can do

None escapes the circularity finding. `top100` is still a linear function of expression, so
predicting same-timepoint ΔAge from expression stays circular whichever clock supplies the target.
Only the forward formulation escapes it, and Phase 3 is the only one of the three that is forward.
