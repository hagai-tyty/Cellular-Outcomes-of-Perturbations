# REFERENCE — Ground rules

**These apply to every stage.** Covers `MASTER_PLAN.md` §8, §9, §10b.

> **🆕 ADDED 2026-07-26 — two rules earned the hard way.** *Additive; no existing rule is modified.*
> The Stage 1.5 arc spent four diagnostics and concluded "the ΔAge target is unvalidated." That
> conclusion was **wrong**, and both causes were method, not data
> (`STAGE_1_5_1_REVISED.md`). These are cheap, and either alone would have prevented it.
>
> **§10 — Every effect claim needs a negative control, and you must report it.**
> Gill's data contained 47 samples of cells that received OSKM and *failed* to reprogram. They cannot
> rejuvenate — yet they read **+36.5 yr older in 11 days**. That is how the artefact was finally
> seen. Before claiming an effect, ask: *what in this dataset should show no effect, and does it?*
> **A claim whose negative control moves with it is not validated.** Corollary: **do not pool
> responders with non-responders** — pooled, the same contrast read +22.4 yr; split, −28.3 yr.
>
> **§11 — State the effect's expected SHAPE before choosing the statistic.**
> E1/E1b used a monotonic Spearman on an effect that *dips and recovers*. A fall-then-rise
> rank-correlates to ≈0, so the test returned `NO_TREND` and `WRONG_DIRECTION` **on data containing a
> significant −28.3 yr effect**. A statistic that cannot detect the shape you expect is not a test of
> it. For a non-monotonic effect, use a **contrast at the predicted window**, and say where that
> window comes from *before* running.

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

## 5b. A pre-set bar must also be RESOLVABLE — check before running, not after

Setting the threshold before the run (§5) is necessary but **not sufficient**. A bar can be
honestly pre-registered and still be untestable: if the estimator is noisy at the geometry it is
graded on, a system that *fully meets the intent* fails the bar anyway, and both passing and
failing it are uninformative. This is not hypothetical — it happened twice in Stage 1, and both
times we found it **after** the run, not before:

- **`fate_ece ≤ 0.169`, graded as the mean of per-fold ECEs (n≈21, 10 bins).** A *perfectly
  calibrated* model scores 0.183 and clears the bar only **26.9%** of the time. The bar was
  measuring the sample size. (Fix: pool across folds → floor 0.091, a correct system passes 99.6%.)
- **`conformal_coverage ∈ [0.85, 0.95]`.** This one *survived* the check — a correctly-90% system
  lands in-band **93%** of the time — but that was confirmed, not assumed.

> **Before a bar is pre-registered, simulate a system that meets its intent EXACTLY at the geometry
> the bar will be graded on, and confirm that system passes at least 95% of the time. If it does
> not, the bar is unresolvable: move the threshold to `usable_bar`, or change the geometry (pool,
> add cells), or drop the criterion — but do it now, not after a run wears the failure.**

The check is one call: `audit_metrics.bar_verdict(null, bar, lower_is_better)` returns
**RESOLVABLE / UNRESOLVABLE** against `MIN_PASS_RATE = 0.95`, where `null` is the metric simulated
under its own intent (`y ~ Bernoulli(p)` for calibration; `hits ~ Bernoulli(level)` for coverage).
Every registered TARGET bar has an entry in `tests/test_bars_resolvable.py`, so adding a bar means
adding its resolvability test — a bar with no such test is not considered pre-registered.

This is the forward form of the audit that caught both Stage 1 problems. The habit is: **audit the
bar before the run, not the run after the bar.**

## 5b-bis. Two things §5b does not yet say — both learned in Stage 1.5.2

*Additive; §5b above is unmodified.*

**(i) A pre-committed FALLBACK is a bar, and must be resolvability-checked too.**
`STAGE_1_5_2` §6 anticipated that its primary criterion might be unresolvable and pre-registered a
fallback — which is good practice, and it registered the fallback **without checking it**. Measured,
**the fallback was itself UNRESOLVABLE at 92.3%**. On the originally-planned geometry the stage would
have had *no valid decisive criterion at all* and would not have known.

> **Every branch of a pre-registration is a bar. Check the ones you hope not to use.**

**(ii) A bar near the instrument's empirical ceiling is not a bar.**
§5b asks whether a system meeting the intent *exactly* passes. It does not ask whether **anything**
can reach the bar on this data. Stage 1.5.2 set ρ ≥ 0.50 against a null simulated at ρ_true = 0.70 —
and then measured that **two clocks of the same modality, on the same samples, agree with each other
at only ρ_partial +0.568** (§12-R), varying **+0.233 to +0.936** by cell state (§17).

> **Measure the ceiling — the agreement achievable between two instruments of the same kind on the
> same data — before setting a bar near it. A bar the reference standard cannot itself clear tests
> the data, not the candidate.**

Cheapest form: score a second instance of the *reference* modality by the identical criterion, on
the identical samples. If it fails too, the bar is the problem.

## 6. When a result surprises you

The default assumption is a **bug in the test**, not a discovery. Precedents from this project:

- Test 7.3's null was **guaranteed by construction** (rank-invariant metric, monotonic transform)
- Test 7.4's precision half was **structurally void** (quantile threshold is monotone-equivariant)
- Test 7.4.1's `(n=0) -> tied` lines were **artifacts** (nothing to compare)
- A 6/6 clean sweep in Stage 4 should trigger a **leakage audit**, not celebration

**Check the test before believing the result.**
