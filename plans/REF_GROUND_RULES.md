# REFERENCE — Ground rules

**These apply to every stage.** Covers `MASTER_PLAN.md` §8, §9, §10b.

---

## 1. The evidence provenance rule

An external review once asserted a "biphasic MPTR trajectory (+30 to +50 years)" as established
fact, citing analyses that **did not exist in this project's record**. The claim later turned out
to be *right* (Test 13: 4/6 folds show a hump) — but it was tested, not assumed, and that is the
point.

> **A claim enters the plan as an established finding only when reproduced from our own data, in
> our own notebook, with a recorded test. Everything else enters as a hypothesis with a test
> attached — never as a premise.**

This applies to input from humans, AIs, and literature alike. It is the same standard that caught
fabricated citations earlier in the project, and that killed five of our own claims:

- ~~"Fate fails out-of-donor"~~ → discrimination holds (0.96–1.00)
- ~~"The perturbation input is constant"~~ → `dose_time` varies richly
- ~~"Recalibration doesn't help RES"~~ → it does, where it structurally can
- ~~"There is a systematic global ΔAge offset"~~ → calib offset is −0.03; the shift is per-donor
- ~~"RES fails because fate probabilities are miscalibrated"~~ → the primary cause is `g(R_eff)=0`

## 2. The change protocol

```powershell
python scorecard.py snapshot --tag baseline    # BEFORE any change
# ... exactly ONE change ...
python scorecard.py snapshot --tag <name>
python scorecard.py compare baseline <name>
```

`compare` reports the **paired per-fold difference and 95% CI** for every metric and returns one of:

| Verdict | Meaning |
|---|---|
| **ACCEPT (better)** | CI excludes 0 in the improving direction → real |
| **noise (CI incl. 0)** | not distinguishable from fold variation → **not** an improvement |
| **REGRESSION** | CI excludes 0 the wrong way → revert |
| **(context)** | not a quality metric (`res_approvals`, `ood_rate`, `n_cells`) |

## 3. The six rules

1. **Snapshot the baseline before Stage 1.** Without it, no later comparison means anything.
2. **One change per snapshot** — except Stage 2, which bundles B + A′ because B changes the error
   scale A′ depends on. That exception is sanctioned in writing; no others are.
3. **Accept only if the TARGET metric says ACCEPT and no GUARD says REGRESSION.** "It looks a bit
   better" is not acceptance.
4. **Calibration is always fitted last**, to the model that will actually ship.
5. **Split by donor** — never by pair, never by cell.
6. **Record every stage in the lab notebook, including failures.** The failures are what make the
   record trustworthy.

## 4. Metrics whose arrow must be read against the plan, not the scorecard

| Metric | Scorecard says | But in Stage 1 |
|---|---|---|
| `conformal_width` | lower is better | **rising is CORRECT** — 8.9 → ~35–43 yr. The current narrow interval is the defect |
| `res_approvals` | (context) | **more is not better** — the model already over-approves vs the oracle (14 vs 11) |

## 5. Thresholds are set before running, never after

Every stage document states its accept bar **before** the change is implemented. Choosing a
threshold after seeing results is how every change comes to look like an improvement.

If a result lands just outside its bar, the honest options are: **accept the failure**, or **run a
new test with a new pre-registered bar** — never retroactively widen the old one.

## 6. When a result surprises you

The default assumption is a **bug in the test**, not a discovery. Precedents from this project:

- Test 7.3's null was **guaranteed by construction** (rank-invariant metric, monotonic transform)
- Test 7.4's precision half was **structurally void** (quantile threshold is monotone-equivariant)
- Test 7.4.1's `(n=0) -> tied` lines were **artifacts** (nothing to compare)
- A 6/6 clean sweep in Stage 4 should trigger a **leakage audit**, not celebration

**Check the test before believing the result.**
