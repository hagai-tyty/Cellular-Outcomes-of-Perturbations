# plans/archive — superseded drafts, kept deliberately

**Nothing here is current. Nothing here should be acted on.** These files are the **audit trail**,
and they are the reason this directory exists rather than a `git rm`.

## Why they were not deleted

`STAGE_1_5_1_REV_FINAL.md` §10.7 records, as a verified check, that *"the five earlier 1.5.1 drafts"*
are **byte-unmodified**. That check only means something if the files are readable next to the final
document — a reviewer has to be able to see what changed between the drafts and the conclusion, and
what was claimed before it was withdrawn. Git history preserves the bytes either way; it does not
put them where a reviewer will look.

They are also **cited by nine other documents**, including `MASTER_PLAN.md`, `REF_GROUND_RULES.md`
and Stages 2–6. `STAGE_1_5_1_REVISED.md` alone is referenced twelve times.

**Verified on the move (2026-08-01):** all five SHA-256 hashes are identical before and after, and
`git mv` was used so history follows the file.

## What is here

| file | what it was | superseded by |
|---|---|---|
| `STAGE_1_5_1_CLOCK_PRECISION.md` | V1 — the first 1.5.1 plan | `STAGE_1_5_1_REV_FINAL.md` |
| `STAGE_1_5_1_NEW_CHANGES.md` | the change list that followed V1's review | ↑ |
| `STAGE_1_5_1_NEW_V2.md` | V2, adjudicating the critique of V1 | ↑ |
| `STAGE_1_5_1_REVISED.md` | the REVISED plan. **Its mechanism was later refuted by measurement** — it attributed the +36.5 yr artefact to loss of fibroblast identity, and non-responders were then shown to have day-0 identity while still reading +33 yr older | ↑ |
| `STAGE_1_5_1_REVISED_REVIEW.md` | the review of REVISED | ↑ |

## Where to look instead

| you want | read |
|---|---|
| the 1.5.1 method, results and open items | `../STAGE_1_5_1_REV_FINAL.md` — start at §11 |
| whether the RNA clock is calibratable | `../STAGE_1_5_2_LABEL_ANCHOR.md` — §16, then §17 |
| the `src/` changes that follow from that | `../STAGE_1_5_3_EXECUTE.md` |
