# STAGE 21C RECORD — Public Prospective Dataset Qualification

## Goal
Qualify public datasets for the prospective CellFate-Rx path under the frozen two-pass search budget.

## Inputs
- Stage 21 V3 qualification plan
- frozen 21A / 21B results
- external paper/GEO qualification research
- six pre-registered Pass-1 candidate families

## Search budget
Pass 1:
- Rewind / GSE227151
- CellTag / GSE99915
- CellTag-multi / GSE216518 + GSE216521
- ReSisTrace / GSE223003
- GSE279162
- GSE253739

Pass 2:
NOT TRIGGERED

Reason:
Both Role A and Role B were resolved in Pass 1.

## Result

Role A:
- GSE227151 Rewind -> QUALIFIED_ROLE_A
- selected primary Role-A dataset

Role B:
- GSE279162 -> QUALIFIED_ROLE_B
- selected primary Role-B dataset
- GSE253739 retained as strong secondary/sequential candidate
- GSE223003 retained as prospective replication candidate

Surrogate-only:
- GSE99915
- GSE216518 / GSE216521

Correction:
- GSE243933 is not the primary Rewind scRNA prospective dataset.
- GSE227151 is the Rewind scRNA anchor used going forward.

## Final verdict
FULL_DATA_PATH

Role A qualified.
Role B qualified.
Pass 2 not required.

## What this proves
Public data exist with the experimental geometry needed to attempt:
X_before + U -> Y_future

The local-data failure does not block the prospective paper path.

## What this does NOT prove
- that the datasets are reconstructable without ambiguity
- exact independent-unit counts
- exact positive/negative clone counts
- model learnability
- X adds beyond U
- interaction X×U exists
- treatment ranking works

Those belong to Stages 21D onward.

## Files selected for Stage 21D
Role A:
D:\GSE227151_Rewind\

Role B:
D:\GSE279162\

## Engineering state
- no model training
- no src/ changes
- no label tuning
- no architecture changes
- raw public datasets remain outside git

## Next action
Stage 21D — Acquisition + Reconstruction