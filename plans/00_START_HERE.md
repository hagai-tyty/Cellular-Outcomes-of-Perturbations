# 00 — START HERE: the operator's guide

**What this is.** The running order for the whole project, from where you are today to a working
tool and a submitted paper. Each stage has its own detailed document; this file tells you **what
to run, what to send back, and what decision follows**.

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

| | Status |
|---|---|
| Diagnostics | **Done** — ~25 tests, Tests 0 through 18 |
| Root cause found | **Yes** — all calibration fitted in-distribution |
| Main fix validated | **Yes** — T16: k=3 reference cells → MAE 14.3 → 7.1 |
| Code changed so far | **None** — every finding is a measurement, not a modification |
| Baseline snapshot | **Not yet taken** ← *this is your next action* |

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
| 5 | `python test18_forward_gate.py` | 3 | **the tool's existential gate** |
| 6 | *(implement Stage 3)* then `python scorecard.py snapshot --tag C_forward` | 3 | the forward tool |
| 7 | `python scorecard.py compare B_percalib C_forward` | 3 | accept / reject |
| 8 | `python validate_stopping.py` | 4 | does it beat a fixed protocol? |
| 9 | `python test19_second_clock.py` | 5 | is the linearity claim clock-dependent? |
| 10 | `python test_suite.py y1_probe` | 5 | is the fate edge a Y1 artifact? |

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
python test18_forward_gate.py
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
python test_suite.py y1_probe      # does the fate claim survive without Y1?
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
| 0 baseline | | — | — | |
| 1 calibration | | | | |
| 2 level correction | | | | |
| 3 gate | | GO / WEAK / STOP | | |
| 3 tool | | | | |
| 4 validation | | | | |
| 5 publication gates | | | | |
