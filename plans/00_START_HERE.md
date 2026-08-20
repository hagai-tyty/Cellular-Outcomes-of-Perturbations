# 00 — START HERE: the operator's guide

**What this is.** The running order for the whole project, from where you are today to a working
tool and a submitted paper. Each stage has its own detailed document; this file tells you **what
to run, what to send back, and what decision follows**.

---

> ### 🟩 CURRENT STATE — updated 2026-08-20. READ THIS FIRST.
>
> **The running order below is the ORIGINAL plan (Stages 1–6 plus the 1.5.x arc). It is left
> unedited, but it stops at 2026-08-04 and does not mention Stages 10–18, which are where the
> project actually went.** Nothing below is wrong about *its own* stages; it is simply incomplete.
>
> **For what is true today, read [`../ARCHITECTURE.md` §11](../ARCHITECTURE.md), not this file.**
>
> **Stages 10–18 — all executed and closed:**
>
> | | | |
> |---|---|---|
> | **10** | pluripotency → MEDIATION, removal withdrawn | `STAGE_10_…` |
> | **11** | ΔAge is mis-scaled, not broken | `STAGE_11_…` |
> | **12** | `cell_id` was not unique — fixed, rebuilt, null on metrics | `STAGE_12_…` |
> | **13** | scorecard judged signed level shift — 12 of 20 verdicts changed | `STAGE_13_…` |
> | **14** | calibrated ΔAge shipped at the reporting boundary | `STAGE_14_…` |
> | **15** | RES is zero because σ ≫ \|µ\| — attributed | `STAGE_15_RES_ZERO_RECORD.md` |
> | **16** | fate calibrator was fitted on the wrong target — fixed and verified | `STAGE_16_…` |
> | **17** | scorecard rewarded over-coverage — fixed | `STAGE_17_…` |
> | **18** | the fate head is partly a clock; 0.917 on 12 within-timepoint pairs | `STAGE_18_FATE_BEYOND_DAY_RECORD.md` |
>
> **What is still open:** Stage 2 (level correction), Stages 3–5 (tool, validation, publication),
> and `STAGE_6_NEW_DATA_REV.md` / `DATA_REQUIREMENT_SECOND_TIMECOURSE.md` — the data acquisition
> that everything else now waits on.
>
> **Current fold set: `_s16`.** Grade every change with
> `python scorecard.py compare c7t_stage16 <tag>`.

---

## What each file in this folder is

**Three kinds of document. Only one kind describes work you perform.**

| File | Kind | When you open it |
|---|---|---|
| **`00_START_HERE.md`** (this file) | **guide** | first, and between every stage — it says what to type |
| **`STAGE_1_CALIBRATION.md`** | **DO** | when implementing Stage 1 |
| **`STAGE_2_LEVEL_CORRECTION.md`** | **DO** | when implementing Stage 2 |
| **`STAGE_3_TOOL.md`** | **DO** | when implementing Stage 3 (4 internal sub-stages) |
| **`STAGE_4_VALIDATION.md`** | **DO** | when validating |
| **`STAGE_5_PUBLICATION.md`** | **DO** | before writing the manuscript |
| **`STAGE_6_NEW_DATA.md`** | **DO** | acquiring/integrating data — runs in PARALLEL, not after |
| **`STAGE_1_5_1_REV_FINAL.md`** | *closed* | 🆕 the ΔAge methylation anchor. **Read §11** — every open item now has an owner |
| **`STAGE_1_5_2_LABEL_ANCHOR.md`** | *closed* | 🆕 **is the RNA clock calibratable? Answer: NO.** Read §16 (the one-page answer), then §17 (the re-audit) |
| **`STAGE_1_5_3_EXECUTE.md`** | **DO** | 🆕 **the `src/` changes 1.5.2 forces.** Read before touching `aging.py`, `res.py` or any ΔAge label |
| `REF_ARCHITECTURE.md` | *reference* | when you want to know **why** the stages are ordered this way, or what each one buys |
| `REF_GROUND_RULES.md` | *reference* | the rules that apply to **every** stage — read once, then consult |
| `REF_DATA_STRATEGY.md` | *reference* | **why** more data is needed (the *doing* is Stage 6) |
| `MASTER_PLAN.md` | *source* | the document the stages decompose; consult when a stage seems to contradict something |

> **The five `STAGE_*` files are the work.** The three `REF_*` files are context — earlier drafts
> numbered them `01`, `02`, `06`, which made them look like stages 1, 2 and 6. They are not.

---

**How to use it.** Work top to bottom. After each command, paste the output back. Every stage has
a pre-registered accept/reject rule, so the decision is already made before you see the numbers.

---

## Where you are right now

> **Status only — updated 2026-07-22.** Nothing below this heading changes the plan: no stage, no
> command, no acceptance bar has been edited. Only the record of what has already been executed.
> Full detail in `experiments/DELTAAGE_LAB_NOTEBOOK.md`, `CHANGES.md` and
> `plans/STAGE_1_DEVIATIONS.md`.

| | Status |
|---|---|
| Diagnostics | **Done** — ~25 tests, Tests 0 through 18 |
| Root cause found | **Yes** — all calibration fitted in-distribution |
| Main fix validated | **Yes** — T16: k=3 reference cells → MAE 14.3 → 7.1 |
| Baseline snapshot | **Taken 2026-07-19** — `scorecard/baseline.json`, all 6 folds ok |
| Test 18 gate | **Run — STOP.** Δt cannot predict the unsafe fraction forward. Also measured: **1.8 cells/timepoint**, so per-timepoint SE (12.9–15.9 yr) exceeds the ~11.35 yr effect on every donor |
| Stage 1a | **Done and passing** — donor column from `cell_line`; `verify_1a.py` PASS, 7 tensors on all splits, exactly 5 usable donors per fold |
| Stage 1b | **Run twice.** Run 1 (07-21) **invalid** — a bulk corpus (HFF) was rotated as a donor and supplied 99.8% of residuals. Run 2 (07-22) **valid** |
| Stage 1 acceptance | **PARTIAL.** `conformal_coverage` 0.401 → **0.889 ACCEPT**; `fate_ece` 0.281 → **0.364 REGRESSION**; all six guards **bit-identical** (0.00e+00). Under §3's independence clause `q` and `sigma_scale` are **adopted**, the temperature refit **rejected** |
| Open | **Change A″** — recalibrate `P(safe)` (the quantity `res.py` and the scorecard actually use) instead of multi-class NLL. Implemented and tested; awaiting a third LOOCV run. Bar unchanged |

---

## 🆕 THE ORDER OF WORK FROM HERE — decided 2026-08-02

> *Additive. The 2026-07-22 status table above is left exactly as written; it predates the entire
> Stage 1.5 → 1.5.3 arc and is kept as the record of that moment, not corrected.*

### Where the project actually stands

| | |
|---|---|
| **Fate head** | ✅ **works** — `fate_roc` 0.983, `fate_prauc` 0.992. Untouched by every ΔAge problem below |
| **ΔAge, relative** | ⚠️ `rank_model_dage` 0.948 — but that is the model reproducing **its own labels'** ordering. Internally consistent, **not externally validated** |
| **ΔAge, absolute** | 🔴 **broken** — `dage_mae_model` 14.29 yr against a clock whose own error is 12.27 (**SNR ≈ 0.90**); `conformal_coverage` 0.401 against a nominal 0.9 |
| **The product (RES)** | 🔴 `res_median` **0.000 on every fold** — `R_eff` needs `mu_age < −30 yr` against a measured ~−11, so it rejects everything, and not because the treatment failed |
| **Why** | the RNA clock is out of domain on reprogramming cells (1.5.1), and **not calibratable** against methylation (1.5.2 M-2a, SPLIT) |
| **What is real** | methylation measures **−24 to −28 yr** with an inert negative control and a dose-response at p = 0.0001. **The biology is there; the RNA instrument cannot see it** |

**ΔAge is not optional and is not being dropped.** The whole lab notebook exists to make it right.
The order below is how.

### The order

| # | do this | why here | cost |
|---|---|---|---|
| **1** | **Finish 1.5.3 — run step 6.** The only item left | Built, gated, bar registered. Its answer redirects everything downstream: *B better* ⇒ Stage 6 must **replace** HFF's labels; *A better* ⇒ it may **supplement**; *INCONCLUSIVE* ⇒ the blocker is donor count | one retrain |
| **2** | **Stage 1.5.4** *(to be written)* — can a model **learn** age from RNA? | M-2a tested whether the **existing Fleischer clock's output** tracks methylation age → no. **Nobody asked whether a model trained on the transcriptome can.** Different question: that clock was fitted on quiescent fibroblasts against *chronological* age. ~90 paired RNA+methylation samples are already on disk | **free** |
| **3** | **Integrate GSE165177** | 95 adult samples, donors **38 and 53**, all inside the clock's `[1, 96]`, all methylation-paired. Takes in-range labels **~45 → ~140**. Already downloaded, referenced by no training config | **free** |
| **4** | ✅ **`STAGE_6_NEW_DATA_REV.md`** — written 2026-08-03 | The original names **2 of the 6 datasets you hold**. The rewrite carries forward its species blocker and `input_ablation` gate unchanged, adds the real inventory, and **sizes the ask in donors** | done |
| **5** | **Execute Stage 6 — acquire** | The actual fix. Sized by #4 | 💰 real |

**Steps 1–3 are free and none of them fixes ΔAge at scale. Only #5 does.** 1–4 exist so that #5 buys
the right thing.

### The number that decides everything, and nobody has computed it

**How many DONORS does the age arm need?** Not samples — **donors**. Every unresolved statistic in
this project is donor-limited, not sample-limited:

* **M3** — is the per-donor offset real or `n=1` baseline noise? *56 % of variance, 95 % CI
  **[9 %, 100 %]*** — unresolvable at 6 donors.
* **Contrast B** (retention) — needs **≈16 pairs**, has **9**.
* **Step 6's own MDE** — `1.049 × SD(per-fold difference)` at 6 folds; Δ\* = 3.57 yr is detectable
  only if that SD is **≤ ~1 yr**.

**GSE165177 triples the labels and adds ONE donor.** That is the trap to avoid in step #4: sizing
Stage 6 on sample count would buy the wrong thing.

### Two facts to carry into step 6

* **`age_window_k = 4` in BOTH arms.** The default is 1, which means OFF — arm B would be starved
  and the confound 5c exists to remove would come back silently.
* **Report the observed per-fold SD and MDE beside the effect.** Which row of the outcome table
  applies *depends* on the MDE. A CI containing 0 with MDE > Δ\* is **INCONCLUSIVE and licenses
  nothing** — in particular it does **not** license discarding HFF's labels.

### Why HFF may not be fixable at all

HFF is a **neonatal foreskin fibroblast line, donor age 0**, outside the clock's fitted range, and
its day-0 baseline reads **84.5 yr**. But the deeper issue is not the instrument:

> **A neonatal cell has almost no chronological age to remove.** Reprogramming resets epigenetic age
> toward embryonic, and a newborn line is already near zero. Even a *perfect* instrument would read
> HFF's true ΔAge as ≈ 0.

Better measurement on HFF buys an accurate ≈ 0. That is why the age arm needs **adult-donor** data,
and why "B better" is the *expected* outcome of step 6 rather than a surprise.

> **Neither 1.5.4 nor the Stage 6 rewrite exists yet.** They are named here as work, not linked as
> documents — the file table above lists only files that exist, deliberately. `STAGE_6_NEW_DATA.md`
> once carried an acceptance gate naming a test that had never been written, and that gate could
> therefore never fail; the same mistake is not repeated here.

---

## 🆕 2026-08-07 — REPRODUCTION STATUS: the pipeline is the same system it was in July

> *Additive. Nothing above is edited. This records what was executed, not a change of plan.*

The other machine's July results were re-run against arm A, the current true-label build. Both
sides are the same two-dataset harmonized build (42605 cells, 51 shards, panel
`783f269a214aa972`), verified from metadata before running.

| what was re-run | result |
|---|---|
| Test 7 / 7.1 / 7.2 — ranking by ΔAge | **EXACT** on 6/6 folds. `model_dAge` 0.948, `ridge_dAge` 0.955, Δ = −0.000, per-fold identical to three decimals |
| G-c step 1 — HFF ΔAge trajectory | **BIT-FOR-BIT** on the matched fold. ρ −0.9047619047619048, slope −1.5255573306808494, day-14 −24.023, all 8 leave-one-timepoint-out folds identical |
| Change A's two invariance guards | **BOTH HELD** — `rank_model_dage` exactly 0.948, `ood_rate` exactly 0.2732 against a pre-registered 0.273 |
| RES | **degenerate, as pre-registered.** Constant (all-zero) on 3/6 folds; `R_eff = 0` for 100% of cells because σ_age is now ~37 yr. This is `sigma_scale` working, not a regression — it confirms the `res_median` **0.000 on every fold** line in the table above and supplies its mechanism |

**So the July record and the current build describe the same system.** Nothing in the notebook
needs re-deriving on that account.

### The one NEW problem this surfaced

**HFF's ΔAge labels are not stable across LOOCV folds.** Same script, same build family, only the
held-out Gill donor differs:

| fold | N2 | N3 | O1 | O2 | Y1 | Y2 |
|---|---|---|---|---|---|---|
| day-14 ΔAge | **−7.35** | −22.12 | −24.02 | −22.89 | −22.05 | −23.87 |

Spread **16.7 yr** against a pre-registered 2.0 yr tolerance; N2 is a **3.1× compression**. HFF is
42481 of 42605 age-labelled cells (99.7%) and is never the held-out line — the withheld donor is
~21 cells, 0.05% of the corpus. The training target for 99.7% of the data should not move 3× when
0.05% is withheld.

**Why this matters for the order of work above.** Step 6 was INCONCLUSIVE because its per-fold SD
was 4.808, giving MDE 5.045 against Δ* = 3.572 — and the section above already names that SD as
the thing gating everything: *"Δ* = 3.57 yr is detectable only if that SD is ≤ ~1 yr."* A 3×
label swing in one of six folds is a **candidate** source of exactly that variance. If it is, this
is a bug to fix rather than donors to buy — and that would be the cheapest thing on this page.
**Stated as a hypothesis. Nothing here measures its contribution, and no recorded result is
withdrawn on the strength of it.**

Suspect: harmonization is refit per fold (gene set varies 5026–5402, N2 fewest) against the small
`gill_bulk` reference. Not established — Y1 has a similar `gill_bulk` profile and does not
collapse.

Full detail: `experiments/DELTAAGE_LAB_NOTEBOOK.md` (last two sections) and `CHANGES.md`.

---

## 🆕 2026-08-07 — C-2 IS NOT FREE. It masks 99.8% of the age labels.

> *Additive. No ordering above is edited. This records a measurement and states its consequence for
> the order; the order itself is a decision, not a finding.*

The standing shortlist has carried, for weeks:

> *"1. Turn on C-2 — free, code written, worth more than 13 donors"*

**Verified against the data and the code — it is not free, and "worth more than 13 donors" is
backwards.**

| fact | source | value |
|---|---|---|
| the clock's declared validity range | `fleischer_clock.json` -> `meta.age_range`, lifted by `aging.py:80` | **[1.0, 96.0]** |
| Gill donor ages | `GSE165176_series_matrix.txt.gz`, `!Sample_characteristics_ch1` | **N2 = 0, N3 = 0**, O1 = 53, O2 = 53, Y1 = 29, Y2 = 35 |
| HFF donor age | `sources.py:557` `DONOR_AGE_YEARS` | **0.0** |
| is HFF's age stamped? | `sources.py:731` — **C-3 shipped** | **yes** |
| what rule 3 does | `aging.py:200` `_exclude((a < lo) | (a > hi), "donor_out_of_clock_range")` | masks on `donor_age` |

**HFF is a neonatal foreskin line at donor age 0, below the clock's lower bound of 1.0, and C-3
already stamps that age. So enabling `enforce_clock_age_range` masks every HFF cell** — 42,481 of
42,605 age-labelled cells in the arm-A build — **plus N2 and N3.**

### The comment that justified "free" excludes the dataset that dominates

`build_dataset.py:123`:

> *"OFF by default because turning it on MOVES LABELS -- N2 and N3 are donor_age 0 and the shipped
> clock's range starts at 1.0, so 30 of the **75 non-HFF** training labels would go."*

That count is correct **and scoped to non-HFF**. It was written before C-3 stamped HFF's age. With
C-3 shipped, the rule reaches HFF, and the true cost is ~99.8% of the corpus rather than 30 labels.
**Nothing here is a bug — the flag is off, the guard held, and the comment was true when written.
What is wrong is the shortlist's cost estimate.**

### The uncomfortable part: C-2 firing on HFF may be CORRECT

Every HFF age label is the Fleischer clock extrapolated **below its own fitted range**. That is a
validity problem, not a tuning problem, and it is a candidate explanation for why HFF's ΔAge has
resisted every fix attempted in 1.5.2 through 1.5.6. **The rule is not wrong for wanting to mask
them.**

But acting on it before in-range donors exist would end ΔAge as a trainable target on the current
corpus. **So C-2 is not an early free win — it is DOWNSTREAM of acquisition.** That inverts its
position in the shortlist, where it sat first and acquisition sat last.

### It also touches §4.7's fold instability, without explaining it

**The fold that collapses — N2, day-14 ΔAge −7.35 against a median of −22.51 — is a donor-age-0
donor.** N3 is also age 0 and does **not** collapse, so age alone is not the mechanism and this is
**not** an explanation. What N2 additionally has is the **fewest admissible genes (5026)**, which is
the variance-floor lever recorded in `STAGE_1_5_6_SPARSE_CLOCK.md` §5.2 A4. Recorded as a lead, not
a finding.

> **[2026-08-08 — BOTH LEADS ANSWERED, and neither was the mechanism.]**
> *Nothing above is changed.* **Donor age 0** is not it — N3 is also 0 and does not collapse.
> **The variance floor** is not it — §5.5 eliminated it at *maximum* leverage (F = 1.398, worse
> than not correcting at all), and the fold whose floor ratio IS anomalous (Y1, 34% low) has
> entirely normal labels.
>
> **It is the SAMPLE.** `N2_Fib_Sendai_Exp2` is a degenerate column in GEO's deposit — mean−min
> 0.0008, log2 range 1.74, library 1.03e+08 against a cohort ~1.5e+06 — and §5.14's step-3c run
> **ATTRIBUTED** the whole 16.67 yr spread to it: removing it from the O1 fold moves HFF's day-14
> from −26.8 to −8.2, **+18.6 yr, 6.23× the largest healthy control drop and the only one with the
> opposite sign.** Stage 1.5.6 is CLOSED (§5.15); the fix is owned by
> `plans/CHANGE_C7_BULK_SAMPLE_INTEGRITY.md`.

### Where the shortlist stands, item by item

| # | item | status |
|---|---|---|
| 1 | Turn on C-2 | 🔴 **not free** — masks 99.8% of labels; belongs **after** acquisition, not before |
| 2 | Integrate GSE165177 | ✅ still free, still correct, already in the 2026-08-07 order |
| 3 | Stratified shuffle | ✅ **DONE** — arm D landed 2026-08-07 (`scorecard/gc2_D_stratshuffle_hff_s0.json`) |
| 4 | Acquire | ✅ = Stage 6, in the order — **and now C-2's prerequisite, not its successor** |

### Proposed amendment to the 2026-08-07 order — a PROPOSAL, not an edit

*Recorded here for decision. The order above is unchanged until that decision is taken.*

1. **Merge 1.5.6 step 3b with §4.6 option 2** — §5.2 U2 shows they are the same variance-floor
   question. One read-only measurement, not two.
2. **Re-instrument 3b as a per-fold reconstruction** rather than the spread ratio `R` (§5.2 A1/A2),
   and replace the ATTRIBUTED branch's leaky remedy (§5.2 A3).
3. Steps 2 and 3 of 1.5.6 are unaffected and remain free to run now.
4. **Insert C-2 explicitly into the order, positioned after Stage 6**, carrying the 99.8% figure —
   so a switch this destructive can never be flipped as a tidy-up.

---

## 🆕 Where things live (repo tidied 2026-08-01)

Root now holds **only** the tools you run constantly and the libraries other code imports.
Everything else moved, and every folder has a `HOW_TO_RUN` explaining what is in it.

| folder | what is in it | how to run |
|---|---|---|
| **root** | `scorecard.py`, `retrain_stage1.py`, `audit_metrics.py`, the `diag_*`/`dump_*` tools | `python scorecard.py …` |
| **`plan_tests/`** | 🆕 the **per-stage verification gates** — `verify_1a.py`, `verify_stage1_5.py`, `smoke_stage1.py` | [`plan_tests/HOW_TO_RUN.md`](../plan_tests/HOW_TO_RUN.md) |
| **`experiments/`** | exploratory + numbered tests — `test18_forward_gate.py`, `test5_ridge_gap.py`, every `diag_*` | `experiments/HOW_TO_RUN.txt` |
| **`local_runners/`** | the pipeline runners | `local_runners/HOW_TO_RUN.txt` |
| **`tests/`** | the `pytest` suite — runs in CI, needs no data | `python -m pytest tests/ -q` |
| **`results/`** | 🆕 every `*_results.json` a script writes | read-only output |
| **`plans/archive/`** | 🆕 superseded drafts, kept as the audit trail | [`plans/archive/README.md`](archive/README.md) |

**Everything still runs from the repo root**, whichever folder the script lives in.

---

## Every command you will run, in order

```powershell
# always first, in every session
D:\.venv-cellfate\Scripts\Activate.ps1
cd D:\cellfate-rx
```

| # | Command | Stage | Purpose |
|---|---|---|---|
| 0 | `python scorecard.py snapshot --tag baseline` | — | **freeze the reference point** |
| 1 | *(implement Stage 1)* then `python scorecard.py snapshot --tag A_xdonor` | 1 | cross-donor calibration |
| 2 | `python scorecard.py compare baseline A_xdonor` | 1 | accept / reject |
| 3 | *(implement Stage 2)* then `python scorecard.py snapshot --tag B_percalib` | 2 | level correction |
| 4 | `python scorecard.py compare A_xdonor B_percalib` | 2 | accept / reject |
| 5 | `python experiments/test18_forward_gate.py` | 3 | **the tool's existential gate** |
| 6 | *(implement Stage 3)* then `python scorecard.py snapshot --tag C_forward` | 3 | the forward tool |
| 7 | `python scorecard.py compare B_percalib C_forward` | 3 | accept / reject |
| 8 | `python validate_stopping.py` ⚠️ | 4 | does it beat a fixed protocol? |
| 9 | `python test19_second_clock.py` ⚠️ | 5 | is the linearity claim clock-dependent? |
| 10 | `python experiments/test_suite.py y1_probe` | 5 | is the fate edge a Y1 artifact? |

> ⚠️ **Not written yet.** Rows 8 and 9 are Stage 4/5 scripts that `STAGE_4_VALIDATION.md` and
> `STAGE_5_PUBLICATION.md` specify but nobody has implemented — writing them is part of those
> stages, not a prerequisite. Every other command in this table points at a file that exists today
> (verified 2026-08-01).

---

## STEP 0 — Do this now

```powershell
python scorecard.py snapshot --tag baseline
```

**Why first:** every later comparison is measured against this. Without it there is no "before,"
and no change can be judged.

**It also closes two open questions for free** — conformal coverage and OOD rate are captured in
the snapshot.

**Send back:** the whole printed table.

**Decision that follows:** if conformal coverage is far off nominal or OOD fires on nearly
everything, those become confirmed upstream causes and Stage 1 becomes even more clearly the
right first move. *(We already expect this from Tests 14 and 15 — the snapshot records it in the
same format everything else will be measured in.)*

---

## STAGE 1 — Calibration → `STAGE_1_CALIBRATION.md`

**What changes:** three lines in `train_model.py` plus one new module. Fits `temperature`, `q`,
`sigma_age` and the OOD reference on **cross-donor** statistics instead of in-distribution ones.

**You run:**
```powershell
python scorecard.py snapshot --tag A_xdonor
python scorecard.py compare baseline A_xdonor
```

**Accept if:** `conformal_coverage` reaches 0.85–0.95 **and** `fate_ece` says ACCEPT with a ≥40%
drop, **and** `fate_prauc`, `fate_roc`, `rank_model_dage`, `dage_mae_model` all say **noise**.

**Expect intervals to get WIDER.** That is correct — the current narrow interval is the defect.
The scorecard marks width "lower is better"; ignore that arrow for this stage only.

**If a guard says REGRESSION:** stop. Stage 1 alters only calibration, so a moving guard means the
implementation touched something it shouldn't. That is a bug, not a trade-off.

---

## STAGE 2 — Level correction → `STAGE_2_LEVEL_CORRECTION.md`

**Decide before coding:** this stage needs **k≈3 cells per new donor with known true ΔAge** —
clock readings on both control and perturbed samples. **Is that experimentally acceptable to you?**
If no, skip Stage 2 entirely; the tool still works, you just cannot report absolute ΔAge.

**You run:**
```powershell
python scorecard.py snapshot --tag B_percalib
python scorecard.py compare A_xdonor B_percalib
```

**Accept if:** `dage_mae_model` says ACCEPT with ≥25% drop (T16 predicts ~50%), **and**
`rank_model_dage` says **noise or ACCEPT — never REGRESSION**, and no fold worsens by >20%.

**The ranking guard is the real test.** A level shift is rank-invariant, so ranking *must not
move*. If it does, the implementation is doing more than shifting.

---

## STAGE 3 — The tool → `STAGE_3_TOOL.md`

**Run the gate first — before any tool code is written:**
```powershell
python experiments\test18_forward_gate.py
```

| Verdict | What you do |
|---|---|
| **GO** | build the tool (Stage 3 continues) |
| **WEAK GO** | build it, with tempered expectations |
| **STOP** | **do not write tool code.** Ship the scoring model; go to Stage 5 |

**A STOP here is a real result, not a failure** — it means this dataset cannot support forward
prediction, which is worth knowing and worth reporting.

**If GO, you run** (after implementation):
```powershell
python scorecard.py snapshot --tag C_forward
python scorecard.py compare B_percalib C_forward
```

**Accept if:** `dt_response` exceeds 2 yr (currently 0.035) and `forward_coverage` lands in
0.85–0.95.

---

## STAGE 4 — Validation → `STAGE_4_VALIDATION.md`

```powershell
python validate_stopping.py
```

**The question that decides whether this is a product:** does following the recommendation beat a
fixed withdrawal day, on ≥4 of 6 held-out donors?

| Result | What you do |
|---|---|
| Wins ≥4/6 | ship as a **recommender** |
| Wins 6/6 | **audit for leakage first** — given how often simple baselines have won here, a clean sweep is more likely a bug |
| Wins <4/6 but calibration holds | ship as a **calibrated readout**, not a recommender |
| Calibration fails | do not ship the uncertainty at all |

---

## STAGE 5 — Publication → `STAGE_5_PUBLICATION.md`

Two gates before writing anything:

```powershell
python test19_second_clock.py      # is "deep ties linear" an artifact of our linear clock?
python experiments/test_suite.py y1_probe      # does the fate claim survive without Y1?
```

**Both can change what you are allowed to claim.** The second clock is the one a reviewer will
raise first: *"you concluded deep models can't beat linear, but your target was generated by a
linear model you fitted yourself."*

---

## Between-stage checklist

Before starting any stage:

- [ ] `(.venv-cellfate)` is in the prompt
- [ ] you are at `D:\cellfate-rx`
- [ ] the previous stage's result is written into the lab notebook — **including failures**
- [ ] exactly **one** change is being made before the next snapshot

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError: pandas` | venv not active | `D:\.venv-cellfate\Scripts\Activate.ps1` |
| "No folds found" | wrong working directory | `cd D:\cellfate-rx` |
| scorecard finds no snapshots | none taken yet | run Step 0 |
| a metric shows `n/a` | that fold lacked class variation | normal; check the fold count |

## Decision log — fill this in as you go

| Stage | Date | Verdict | Accepted? | Notes |
|---|---|---|---|---|
| 0 baseline | 2026-07-19 | — | — | Frozen as `scorecard/baseline.json`. Every number the plans predicted was confirmed: MAE 14.291, rank 0.948/0.955/0.686, ECE 0.281, coverage 0.401 (0.000 on N2/N3), OOD 0.273 |
| **3 gate (run early)** | 2026-07-20 | **STOP** | n/a | Part C tied — Δt cannot predict the unsafe fraction forward. Also **1.8 cells/timepoint**: per-tp SE 12.9–15.9 yr exceeds the 11.35 yr effect on every donor, which breaks the ±3.7–4.6 yr arithmetic in `MASTER_PLAN` §5b-ter. Part B is structurally void (identical swing on all six folds, guaranteed by a linear model in `[x, dt, dt²]`). Stage 3 closed as specified |
| 1 calibration — run 1 | 2026-07-21 | **INVALID** | no | `cell_line` merges the HFF corpus (33,613 cells) with the six donors (~14 each); the inner-LODO rotated over HFF, whose fold trained on 75 cells and supplied **99.8%** of residuals. `verify_1a` printed the warning and graded it PASS — that scoring bug cost 3.5 h |
| 1 calibration — run 2 | 2026-07-22 | **PARTIAL** | q ✅ σ ✅ T ❌ | Guards **bit-identical on all six** (0.00e+00). `conformal_coverage` 0.401→**0.889** ACCEPT ✅. `fate_ece` 0.281→**0.364** **REGRESSION** ❌. §3's independence clause: adopt `q` + `sigma_scale`, reject the temperature refit |
| 1 calibration — run 3 (A″) | pending | | | Platt on `P(safe)` — the quantity `res.py` ships and the scorecard grades. Bar unchanged: ACCEPT + ≥40% drop |
| 2 level correction | not started | | | **Blocked on a wet-lab decision**: k≈3 reference cells per donor with clock readings on control *and* perturbed samples (`STAGE_2` §0). Run 2 strengthened the case — residuals are offset-dominated, and the offset is what makes coverage bimodal |
| 3 tool | **closed** | — | — | Gate returned STOP; `REF_ARCHITECTURE` §8 counts this as one of three honest exits |
| 4 validation | n/a | — | — | V3 (beat a fixed protocol) is moot without a recommender. V1/V2 are Stage 1's bars measured on held-out donors |
| 5 publication gates | not started | | | Second clock + Y1 probe. Both still required before drafting |

> **Full detail lives in `experiments/DELTAAGE_LAB_NOTEBOOK.md`** (per-run results, predictions
> made before each run, and what they cost). Code-level departures from the stage documents are in
> `plans/STAGE_1_DEVIATIONS.md`; every change is in `CHANGES.md`.
