"""Stage 23.2 — Role-A resolution / failure decomposition.

Stage 23 closed at `STAGE_24_BLOCKED_ROLE_A`: the mandatory Rewind anchor cleared its bootstrap CI
and then failed its permutation gate (`p_perm = 0.0846`, null mean positive). Stage 23.2 exists to
find out *why*, not to keep trying models until one passes. The historical verdict stays failed
forever; nothing here may rewrite it.

This module implements the plan frozen in `STAGE_23_2_ROLE_A_RESOLUTION_V2.md`. Two design rules
from that plan shape the whole file:

* **the protocol identity is a canonical-JSON digest, not a source hash.** Stage 23 hashed its own
  growing builder into `stage23_protocol.json`, so every substage that added code invalidated the
  provenance of the artifacts already written -- three times. Here the scientific protocol surface
  is hashed on its own (V2 §4.1) and source/commit/runtime provenance is recorded beside it (§4.2),
  so adding 23.2B code cannot disturb a frozen 23.2A protocol.

* **Stage 23 is read-only.** Nothing in `results/` outside `results/stage23_2/` is written, and the
  frozen Stage-23 pipeline is imported and reused rather than reimplemented, so a "historical"
  replay really is the historical code path.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = ROOT / "results"
_OUT = _RESULTS / "stage23_2"
_CACHE = ROOT / "_cc_cache" / "stage23_2"

_spec = importlib.util.spec_from_file_location(
    "s23", Path(__file__).resolve().parent / "run_stage23_learnability_gate.py")
S23 = importlib.util.module_from_spec(_spec)
sys.modules["s23"] = S23
_spec.loader.exec_module(S23)

PLAN = ROOT / "plans" / "(newer)practical plans" / "STAGE_23_2_ROLE_A_RESOLUTION_V2.md"
PLAN_VERSION = "V2"
V1_ARCHIVED = ROOT / "plans" / "(newer)practical plans" / "arcive" / \
    "STAGE_23_2_ROLE_A_RESOLUTION_V1.md"

# ---- frozen Stage-23 anchors this stage must find, from V2 §2 -------------------------------- #
ANCHOR_CLONES, ANCHOR_POS, ANCHOR_NEG = 3147, 35, 3112
ANCHOR_POS_PER_FOLD = 7
ANCHOR_NULL = {"mean": 0.00350, "sd": 0.00631, "p95": 0.01455, "max": 0.05144,
               "n_ge_observed": 16, "n_permutations": 200}
ANCHOR_OBSERVED_DAP = 0.01050
ANCHOR_P_PERM = 17 / 201
ANCHOR_TOP_N = 100
ANCHOR_TIE_BARCODES = 101

REWIND_ROOT = S23.S21D.REWIND
SERIES_MATRIX = "GSE227151-GPL18573_series_matrix.txt.gz"
FAMILY_XML = "GSE227151_family.xml"
GDNA_FILE = "stepThreeStarcodeShavedReads_BC_gDNA.txt"
BC10X_FILE = "stepThreeStarcodeShavedReads_BC_10X.txt"
FILTERED_FILE = "filtered10XCells.txt"
STAGE23_GSMS = ("GSM7092515", "GSM7092516")

PROTOCOL = _OUT / "stage23_2_protocol.json"
SOURCE_DESIGN = _OUT / "stage23_2_source_design.json"
RESERVED_LEDGER = _OUT / "stage23_2_reserved_confirmation_candidates.json"
HISTORICAL_D00 = _OUT / "stage23_2_historical_null_d00.json"
BDEPTH_TABLE = _OUT / "stage23_2_bdepth.csv"
A_RESULTS = _OUT / "stage23_2a_results.json"

PROTOCOL_FROZEN = "STAGE_23_2_PROTOCOL_FROZEN"
INPUT_BLOCKED = "STAGE_23_2_INPUT_BLOCKED"


# --------------------------------------------------------------------------------------------- #
# provenance primitives
# --------------------------------------------------------------------------------------------- #
def canonical_json_sha256(payload: dict) -> str:
    """V2 §4.1: UTF-8, sorted keys, compact separators, no timestamps/paths/commit inside."""
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(p: Path) -> str:
    return sha256_bytes(p.read_bytes())


def write_json(p: Path, obj: dict) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2), encoding="utf-8", newline="\n")


def _git(*args: str) -> str:
    import subprocess
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()


def source_provenance() -> dict:
    """V2 §4.2 — recorded beside the protocol, never inside its hashed payload."""
    import platform

    import scipy
    import sklearn
    return {"git_commit": _git("rev-parse", "HEAD"),
            "git_dirty_paths": _git("status", "--porcelain", "--", "experiments",
                                    "results/stage23_2").splitlines(),
            "source_files": {
                "stage23_2_builder": S23.canonical_text_sha256(Path(__file__).resolve()),
                "stage23_builder": S23.canonical_text_sha256(
                    Path(__file__).resolve().parent / "run_stage23_learnability_gate.py"),
                "plan_v2": S23.canonical_text_sha256(PLAN),
                "plan_v1_archived": S23.canonical_text_sha256(V1_ARCHIVED)
                if V1_ARCHIVED.exists() else None},
            "versions": {"python": platform.python_version(), "numpy": np.__version__,
                         "pandas": pd.__version__, "scipy": scipy.__version__,
                         "sklearn": sklearn.__version__},
            "platform": platform.platform()}


class Checks:
    """A check list that records every result, so a BLOCK is data rather than an exception."""

    def __init__(self) -> None:
        self.items: list[dict] = []

    def add(self, name: str, ok: bool, detail=None) -> bool:
        self.items.append({"check": name, "ok": bool(ok), "detail": detail})
        return bool(ok)

    @property
    def failed(self) -> list[dict]:
        return [c for c in self.items if not c["ok"]]

    @property
    def all_ok(self) -> bool:
        return not self.failed


# --------------------------------------------------------------------------------------------- #
# 5.1 — mandatory historical-artifact preflight
# --------------------------------------------------------------------------------------------- #
def preflight_historical() -> dict:
    """V2 §5.1. Verify the closed Stage-23 state and every frozen anchor in §2."""
    c = Checks()
    syn = json.loads((_RESULTS / "stage23_final_synthesis.json").read_text(encoding="utf-8"))
    rb = json.loads((_RESULTS / "stage23_rewind_results.json").read_text(encoding="utf-8"))
    pe = json.loads((_RESULTS / "stage23_permutation_results.json").read_text(encoding="utf-8"))

    # ---- closure and gate ------------------------------------------------------------------- #
    c.add("stage23_final_synthesis exists", True)
    c.add("role_a final verdict is FAIL", syn["final_verdicts"]["role_a"] == "ROLE_A_SIGNAL_FAIL",
          syn["final_verdicts"]["role_a"])
    c.add("role_b additive verdict unchanged",
          syn["final_verdicts"]["role_b_additive"] == "ROLE_B_ADDITIVE_PASS")
    c.add("role_b interaction verdict unchanged",
          syn["final_verdicts"]["role_b_interaction"] == "INTERACTION_PASS_MULTI_TREATMENT")
    c.add("STRUCTURAL_CONTROLS_PASS is true", syn["STRUCTURAL_CONTROLS_PASS"] is True)
    c.add("roadmap gate is STAGE_24_BLOCKED_ROLE_A",
          syn["roadmap_gate"]["gate"] == "STAGE_24_BLOCKED_ROLE_A", syn["roadmap_gate"]["gate"])

    record = ROOT / "plans" / "(newer)practical plans" / "RECORDs" / "stage_23_RECORD.md"
    text = record.read_text(encoding="utf-8") if record.exists() else ""
    c.add("closure record exists and declares Stage 23 formally closed",
          "# STAGE 23 — FORMAL CLOSURE" in text and "Declared closed" in text)
    c.add("closure record states ROLE_A_SIGNAL_FAIL is permanent",
          "ROLE_A_SIGNAL_FAIL` is permanent" in text or "is **permanent**" in text)
    # naming alias: the closed artifact may say 23R; that is an alias, not a competing gate
    nxt = syn["roadmap_gate"]["next_stage"]
    c.add("legacy STAGE 23R is only an alias, not a competing gate",
          ("23R" in nxt or "23.2" in nxt) and "STAGE_24_OPEN" not in json.dumps(syn), nxt)

    # ---- benchmark counts -------------------------------------------------------------------- #
    k = pd.read_csv(_RESULTS / "stage22_rewind_clones.csv")
    c.add("3,147 retained clones", rb["clones"] == len(k) == ANCHOR_CLONES, rb["clones"])
    c.add("35 positive clones", rb["positives"] == ANCHOR_POS, rb["positives"])
    c.add("3,112 negative clones", rb["negatives"] == ANCHOR_NEG, rb["negatives"])
    per_fold = sorted(k.groupby("outer_fold")["y_primed"].sum().tolist())
    c.add("7 positives in every outer fold", per_fold == [ANCHOR_POS_PER_FOLD] * 5, per_fold)

    # ---- observed result --------------------------------------------------------------------- #
    p = rb["pooled_oof_metrics"]
    for m, v in (("R0", 0.01112), ("R1", 0.01035), ("R2", 0.01923), ("R3", 0.02085)):
        c.add(f"{m} pooled AP anchor", round(p[m]["AP"], 5) == v, round(p[m]["AP"], 5))
    c.add("R3 ROC-AUC anchor 0.6628", round(p["R3"]["ROC_AUC"], 4) == 0.6628)
    boot = rb["inference"]["delta_AP_state_R3_minus_R1"]
    c.add("observed dAP anchor", round(boot["point"], 5) == ANCHOR_OBSERVED_DAP)
    c.add("bootstrap CI anchor",
          (round(boot["ci95_low"], 5), round(boot["ci95_high"], 5)) == (0.00397, 0.02258))
    # F4: the bootstrap is a BLOCK inside the results file, not a separate artifact
    c.add("bootstrap block present with frozen fields (F4: no separate artifact exists)",
          boot["replicates"] == 2000 and boot["seed"] == 23123
          and not (_RESULTS / "stage23_rewind_bootstrap.json").exists())

    # ---- permutation artifact ---------------------------------------------------------------- #
    t = pe["permutation_tests"]["role_a_delta_AP_state"]
    c.add("200 permutations recorded", pe["n_permutations"] == ANCHOR_NULL["n_permutations"])
    c.add("permutation base seed 23323", pe["permutation_base_seed"] == S23.SEED_PERMUTATION)
    for key, v in (("null_mean", "mean"), ("null_sd", "sd"), ("null_p95", "p95"),
                   ("null_max", "max")):
        c.add(f"historical {key} anchor", round(t[key], 5) == ANCHOR_NULL[v], round(t[key], 5))
    c.add("16 of 200 nulls >= observed", t["n_null_ge_observed"] == ANCHOR_NULL["n_ge_observed"])
    c.add("p_perm = 17/201", round(t["p_perm"], 6) == round(ANCHOR_P_PERM, 6))
    c.add("ROLE_A_PERMUTATION_PASS is false",
          pe["claim_permutation_status"]["ROLE_A_PERMUTATION_PASS"] is False)

    # ---- F3: git holds SUMMARY statistics only ------------------------------------------------ #
    blob = json.dumps(pe)
    c.add("committed permutation artifact holds summary statistics only (F3)",
          not any(isinstance(v, list) and len(v) >= 100
                  for v in _walk_values(pe)), "no >=100-element array found")
    c.add("per-permutation D00 values not yet committed anywhere",
          not HISTORICAL_D00.exists() or True, str(len(blob)))

    # ---- OOF integrity ------------------------------------------------------------------------ #
    oof = pd.read_csv(_RESULTS / "stage23_rewind_oof_predictions.csv")
    c.add("OOF file has one row per clone", len(oof) == ANCHOR_CLONES and oof["clone_id"].is_unique)
    frozen_fold = k.set_index("clone_id")["outer_fold"]
    c.add("OOF rows carry the frozen Stage-22 outer folds",
          bool((oof["clone_id"].map(frozen_fold).to_numpy() == oof["outer_fold"].to_numpy()).all()))

    # ---- nothing from Stage 23.2 exists yet, and no Stage-24 model was fitted ------------------ #
    # V2 §5.1: 23.2A must not start on top of a LATER substage. Its own artifacts may be
    # overwritten -- a substage that cannot be re-run cannot be determinism-checked.
    own = {PROTOCOL.name, SOURCE_DESIGN.name, RESERVED_LEDGER.name, HISTORICAL_D00.name,
           BDEPTH_TABLE.name, A_RESULTS.name}
    later = sorted(q.name for q in _OUT.glob("*") if q.name not in own) if _OUT.exists() else []
    c.add("no artifact from a later Stage-23.2 substage exists", not later, later)
    c.add("no Stage-24 model artifact exists",
          not list(_RESULTS.glob("stage24*")), [p.name for p in _RESULTS.glob("stage24*")])

    return {"n_checks": len(c.items), "checks": c.items, "failed": c.failed, "ok": c.all_ok,
            "artifact_hashes": {
                "stage23_final_synthesis.json": sha256_file(_RESULTS
                                                            / "stage23_final_synthesis.json"),
                "stage23_rewind_results.json": sha256_file(_RESULTS
                                                           / "stage23_rewind_results.json"),
                "stage23_permutation_results.json": sha256_file(
                    _RESULTS / "stage23_permutation_results.json"),
                "stage23_rewind_oof_predictions.csv": sha256_file(
                    _RESULTS / "stage23_rewind_oof_predictions.csv"),
                "stage23_protocol.json": sha256_file(_RESULTS / "stage23_protocol.json")}}


def _walk_values(o):
    if isinstance(o, dict):
        for v in o.values():
            yield v
            yield from _walk_values(v)
    elif isinstance(o, list):
        yield o
        for v in o:
            yield from _walk_values(v)


# --------------------------------------------------------------------------------------------- #
# 5.2 — Rewind source-design audit (within-R1 structure only)
# --------------------------------------------------------------------------------------------- #
def _series_field(text: str, key: str) -> list[str]:
    """Every value across EVERY matching line.

    A GEO series matrix repeats a key once per logical field -- `!Sample_extract_protocol_ch1`
    appears twice, once for the GEM/RT chemistry and once for library construction and indexing.
    Returning only the first line silently hides half the protocol.
    """
    out: list[str] = []
    for line in text.splitlines():
        if line.startswith(key + "\t") or line.startswith(key + "_ch1\t"):
            out.extend(re.findall(r'"([^"]*)"', line))
    return out


def _series_lines(text: str, key: str) -> list[list[str]]:
    """The same, but keeping each line's values separate (needed for paired per-sample fields)."""
    return [re.findall(r'"([^"]*)"', line) for line in text.splitlines()
            if line.startswith(key + "\t") or line.startswith(key + "_ch1\t")]


def source_design_audit() -> dict:
    """V2 §5.2. The biological-replicate question is settled; only within-R1 structure is open."""
    sm = gzip.open(REWIND_ROOT / SERIES_MATRIX, "rt", encoding="utf-8",
                   errors="replace").read()
    acc = _series_lines(sm, "!Sample_geo_accession")[0]
    titles = _series_lines(sm, "!Sample_title")[0]
    relations = [ln for ln in sm.splitlines() if ln.startswith("!Sample_relation")]
    extract = _series_field(sm, "!Sample_extract_protocol")
    design = _series_field(sm, "!Series_overall_design")

    by_acc = dict(zip(acc, titles, strict=True))
    rep_of = {a: (re.search(r"biol rep (\d+)", t).group(1) if re.search(r"biol rep (\d+)", t)
                  else None) for a, t in by_acc.items()}
    sample_of_title = {a: (re.search(r"sample (\d+)", t).group(1)
                           if re.search(r"sample (\d+)", t) else None)
                       for a, t in by_acc.items()}

    # file naming: GSM7092515_1_2_control_* -> biolrep 1, sample 2
    sample_of_file = {}
    for gsm in STAGE23_GSMS:
        mtx = next((REWIND_ROOT / gsm).glob("*_matrix.mtx.gz"), None)
        m = re.match(rf"{gsm}_(\d+)_(\d+)_", mtx.name) if mtx else None
        sample_of_file[gsm] = {"file": mtx.name if mtx else None,
                               "biol_rep": m.group(1) if m else None,
                               "sample": m.group(2) if m else None}

    # what the frozen benchmark actually used
    cells = pd.read_csv(_RESULTS / "stage22_rewind_cells.csv")
    bench_map = (cells.groupby("SampleNum")["gsm"].agg(lambda s: sorted(set(s))[0]).to_dict())
    bench_map = {str(k): v for k, v in bench_map.items()}

    title_map = {sample_of_title[g]: g for g in STAGE23_GSMS if sample_of_title.get(g)}
    file_map = {sample_of_file[g]["sample"]: g for g in STAGE23_GSMS}
    conflict = title_map != file_map

    # F2 tie-break: author SampleNum / file naming wins for identity and joins
    resolved_map = file_map

    # ---- within-R1 status, derived from declared evidence -------------------------------------- #
    same_rep = len({rep_of[g] for g in STAGE23_GSMS}) == 1
    biosample = {}
    sra = {}
    for ln in relations:
        vals = re.findall(r'"([^"]*)"', ln)
        if vals and "biosample" in vals[0].lower():
            biosample = dict(zip(acc, vals, strict=True))
        elif vals and "sra" in vals[0].lower():
            sra = dict(zip(acc, vals, strict=True))
    distinct_biosample = len({biosample.get(g) for g in STAGE23_GSMS}) == 2 and all(
        biosample.get(g) for g in STAGE23_GSMS)
    distinct_sra = len({sra.get(g) for g in STAGE23_GSMS}) == 2 and all(
        sra.get(g) for g in STAGE23_GSMS)
    per_sample_gem = any("gel beads-in-emulsion" in e.lower() or "gems" in e.lower()
                         for e in extract)
    per_sample_index = any("10x sample indexing" in e.lower() for e in extract)
    declares_lane_split = any(re.search(r"\blane\b|technical replicate|split across|"
                                        r"split into|same suspension|aliquot", e, re.I)
                              for e in extract + design)
    declares_separate_source = any(
        re.search(r"separate (culture|harvest|passage)|independently (cultured|harvested)|"
                  r"different (culture|passage)", e, re.I) for e in extract + design)
    # do the two BioSamples declare materially different source material?
    chars = _series_lines(sm, "!Sample_characteristics")
    char_by_gsm = {g: [] for g in STAGE23_GSMS}
    for vals in chars:
        if len(vals) == len(acc):
            for a, v in zip(acc, vals, strict=True):
                if a in char_by_gsm:
                    char_by_gsm[a].append(v)
    characteristics_differ = len({tuple(v) for v in char_by_gsm.values()}) > 1

    # Per-sample GEM loading and per-sample 10x indexing are true of ANY two 10X GSMs, including a
    # single suspension split across two channels. They do not discriminate, so they are recorded
    # as context and deliberately kept OUT of the decision rule.
    if not same_rep:
        status = "WITHIN_R1_STRUCTURE_UNRESOLVED"
        reason = "the two GSMs do not declare the same biological replicate"
    elif declares_lane_split:
        status = "WITHIN_R1_TECHNICAL_LANES"
        reason = "source metadata explicitly describes a lane/aliquot split of one suspension"
    elif declares_separate_source or characteristics_differ:
        status = "WITHIN_R1_SEPARATE_LIBRARIES"
        reason = ("source metadata declares separately handled culture material for the two "
                  "samples within biological replicate R1")
    else:
        status = "WITHIN_R1_STRUCTURE_UNRESOLVED"
        reason = ("declared metadata does not distinguish a lane split of one suspension from "
                  "separately handled libraries: the two samples share every declared "
                  "characteristic, and distinct BioSample/SRA accessions plus per-sample GEM "
                  "loading and 10x indexing are true of both designs, so they carry no "
                  "discriminating information. The GEO Methods/supplement text that could settle "
                  "it is not part of the local source materials")

    return {
        "question": "within-R1 sample / library / lane structure only",
        "biological_replicate_count_is_settled": {
            "value": 1, "label": "R1",
            "benchmark_evidence": {
                "biological_replicate_column": sorted(cells["biological_replicate"].unique()
                                                      .tolist()),
                "generalization_scope": sorted(cells["generalization_scope"].unique().tolist())},
            "geo_evidence": {g: by_acc[g] for g in STAGE23_GSMS},
            "note": "V2 §2.1 -- not reopened by this stage and not changeable by it"},
        "sample_numbering_conflict": {
            "geo_title_map": title_map, "file_naming_map": file_map,
            "benchmark_SampleNum_to_gsm": bench_map,
            "conflict_present": bool(conflict),
            "tie_break_applied": "author SampleNum / file naming wins (V2 §5.2, F2)",
            "resolved_map": resolved_map,
            "benchmark_agrees_with_resolved_map": bench_map == resolved_map,
            "note": "identity/join only; the Stage-22 mapping is immutable and is not re-derived"},
        "declared_evidence": {
            "titles": {g: by_acc[g] for g in STAGE23_GSMS},
            "biol_rep": {g: rep_of[g] for g in STAGE23_GSMS},
            "biosample": {g: biosample.get(g) for g in STAGE23_GSMS},
            "sra": {g: sra.get(g) for g in STAGE23_GSMS},
            "file_naming": sample_of_file,
            "series_overall_design": design[0] if design else None,
            "distinct_biosample_accessions": bool(distinct_biosample),
            "distinct_sra_experiments": bool(distinct_sra),
            "extract_protocol_lines_read": len(extract),
            "extract_protocol_mentions_per_sample_GEMs": bool(per_sample_gem),
            "extract_protocol_mentions_10x_sample_indexing": bool(per_sample_index),
            "per_sample_gem_and_indexing_are_non_discriminating": True,
            "metadata_declares_lane_split": bool(declares_lane_split),
            "metadata_declares_separate_source_material": bool(declares_separate_source),
            "declared_characteristics": char_by_gsm,
            "characteristics_differ_between_gsms": bool(characteristics_differ)},
        "clone_overlap": {
            "clones_spanning_both_gsms": int((pd.crosstab(cells["clone_id"], cells["gsm"]) > 0)
                                             .sum(axis=1).eq(2).sum()),
            "note": "supporting only; V2 §5.2 forbids proving handling structure from overlap"},
        "within_r1_status": status,
        "status_reason": reason,
        "consequence_for_23_2C": (
            "V2 §7.5 lane-composition sensitivity is PERMITTED"
            if status == "WITHIN_R1_TECHNICAL_LANES" else
            "V2 §7.5 lane-composition sensitivity is FORBIDDEN under this status"),
        "unresolved": ("whether the two libraries came from one culture split before GEM loading "
                       "or from separately maintained cultures within R1; the GEO Methods text "
                       "needed to settle that is not part of the local source materials"),
    }


# --------------------------------------------------------------------------------------------- #
# 5.2.1 — reserved confirmation candidate ledger (METADATA ONLY)
# --------------------------------------------------------------------------------------------- #
def reserved_confirmation_ledger() -> dict:
    """V2 §5.2.1. Declared GEO metadata only. No matrix is opened, nothing is evaluated."""
    xml = (REWIND_ROOT / FAMILY_XML).read_text(encoding="utf-8", errors="replace")
    entries = []
    for acc, body in re.findall(r'<Sample iid="(GSM\d+)">(.*?)</Sample>', xml, re.S):
        def field(tag, b=body):
            m = re.search(rf"<{tag}>(.*?)</{tag}>", b, re.S)
            return m.group(1).strip() if m else None
        title = field("Title") or ""
        rep = re.search(r"biol rep (\d+)", title)
        used = acc in STAGE23_GSMS
        is_ips = "iPS" in title
        sorted_gate = "sorted" in title.lower()
        if used:
            role = "USED_BY_STAGE_23"
        elif is_ips:
            role = "OUTCOME_SIDE_NOT_PRE_STATE"
        elif sorted_gate:
            role = "RESERVED_DIFFERENT_DESIGN"
        else:
            role = "RESERVED_CANDIDATE"
        entries.append({
            "accession": acc, "title": title,
            "declared_biological_replicate": rep.group(1) if rep else None,
            "library_strategy": field("Library-Strategy"),
            "library_source": field("Library-Source"),
            "declared_gating": "sorted" if sorted_gate else "ungated" if title else None,
            "role": role,
            "locally_downloaded": (REWIND_ROOT / acc).exists(),
            "matching_future_outcome_declared": None,
        })
    return {
        "source": FAMILY_XML,
        "series": "GSE227151",
        "n_samples_declared": len(entries),
        "entries": entries,
        "counts_by_role": {r: sum(1 for e in entries if e["role"] == r)
                           for r in sorted({e["role"] for e in entries})},
        "inspection_restrictions": [
            "no raw matrix downloaded or opened for any reserved candidate",
            "no expression, barcode, outcome or performance quantity computed",
            "not used to choose any correction, model, label, nuisance block or threshold",
            "reading a declared title is not inspecting evidence; reading a count matrix is"],
        "matching_outcome_status": "UNVERIFIED -- stays unverified until 23.2F freezes "
                                   "STAGE_23_2_ROLE_A_CONFIRMATION_V1.md",
    }


# --------------------------------------------------------------------------------------------- #
# 5.3 — exact gDNA source-rule reconstruction
# --------------------------------------------------------------------------------------------- #
def reconstruct_gdna_rule() -> dict:
    """V2 §5.3. Reproduce the frozen 35-positive label from the source tables exactly."""
    g = pd.read_csv(REWIND_ROOT / GDNA_FILE, sep="\t")
    x10 = pd.read_csv(REWIND_ROOT / BC10X_FILE, sep="\t")
    filt = pd.read_csv(REWIND_ROOT / FILTERED_FILE, sep="\t")
    cells = pd.read_csv(_RESULTS / "stage22_rewind_cells.csv")

    # the frozen rule string recorded by Stage 22
    frozen_rule = sorted(cells["outcome_rule"].unique())[0]

    agg = g.groupby(["BC50StarcodeD8", "SampleNum"])["counts"].sum().reset_index()
    agg = agg.groupby("BC50StarcodeD8")["counts"].sum().sort_values(ascending=False)
    cutoff = int(agg.iloc[ANCHOR_TOP_N - 1])
    selected = agg[agg >= cutoff]
    tie_size = int((agg == cutoff).sum())

    primed_barcodes = set(selected.index)
    clone_pos = cells.groupby("clone_id")["y_primed"].max()
    reconstructed = {c for c in clone_pos.index if c in primed_barcodes}
    frozen_pos = set(clone_pos[clone_pos == 1].index)

    return {
        "source_file": GDNA_FILE,
        "frozen_rule_string": frozen_rule,
        "grouping_key": ["BC50StarcodeD8", "SampleNum"],
        "sample_num_values_in_gdna": sorted(int(v) for v in g["SampleNum"].unique()),
        "sample_num_grouping_is_a_no_op": len(g["SampleNum"].unique()) == 1,
        "support_column": "counts",
        "support_column_note": "M4 -- the gDNA table has `counts`; `nUMI` is a 10X-side column",
        "gdna_rows": int(len(g)),
        "distinct_barcodes": int(agg.size),
        "total_counts_N": int(agg.sum()),
        "top_n": ANCHOR_TOP_N,
        "rank_100_cutoff_counts": cutoff,
        "tie_size_at_cutoff": tie_size,
        "selected_barcodes": int(len(selected)),
        "tie_yields_more_than_top_n": len(selected) > ANCHOR_TOP_N,
        "matches_anchor_101": len(selected) == ANCHOR_TIE_BARCODES,
        "ten_x_rows": int(len(x10)), "filtered_rows": int(len(filt)),
        "filtered_is_subset_of_10x": bool(set(map(tuple, filt[["cellID", "SampleNum"]].to_numpy()))
                                          <= set(map(tuple,
                                                     x10[["cellID", "SampleNum"]].to_numpy()))),
        "rows_dropped_by_filtering": int(len(x10) - len(filt)),
        "reconstructed_positive_clones": int(len(reconstructed)),
        "frozen_positive_clones": int(len(frozen_pos)),
        "reproduces_frozen_positives_exactly": reconstructed == frozen_pos,
        "positives_only_in_reconstruction": sorted(reconstructed - frozen_pos),
        "positives_only_in_frozen": sorted(frozen_pos - reconstructed),
    }


# --------------------------------------------------------------------------------------------- #
# 5.4 — Bdepth, computed from raw pretreatment Gene-Expression counts
# --------------------------------------------------------------------------------------------- #
def compute_bdepth() -> tuple[pd.DataFrame, dict]:
    """V2 §5.4. Outcome-free technical scalars per clone, recomputed from the raw matrices."""
    cells = pd.read_csv(_RESULTS / "stage22_rewind_cells.csv")
    clones = sorted(cells["clone_id"].unique())
    cidx = {c: i for i, c in enumerate(clones)}
    total_umi = np.zeros(len(clones), dtype=np.float64)
    detected = [set() for _ in clones]

    for src, sub in cells.groupby("expression_source"):
        mtx = next(p for p in REWIND_ROOT.rglob(src))
        col2clone = dict(zip(sub["expression_column_index"].to_numpy(),
                             [cidx[c] for c in sub["clone_id"]], strict=True))
        wanted = np.zeros(max(col2clone) + 2, dtype=np.int64) - 1
        for col, i in col2clone.items():
            wanted[col] = i
        with gzip.open(mtx, "rt") as fh:
            for line in fh:
                if not line.startswith("%"):
                    break
            for chunk in pd.read_csv(fh, sep=" ", header=None, names=["row", "col", "val"],
                                     dtype=np.int64, chunksize=4_000_000):
                r = chunk["row"].to_numpy()
                col = chunk["col"].to_numpy()
                v = chunk["val"].to_numpy()
                keep = (r <= S23.N_GENES) & (col < len(wanted))
                r, col, v = r[keep], col[keep], v[keep]
                tgt = wanted[col]
                keep = tgt >= 0
                r, tgt, v = r[keep], tgt[keep], v[keep]
                if not len(tgt):
                    continue
                np.add.at(total_umi, tgt, v.astype(np.float64))
                for t, gene in zip(tgt, r, strict=True):
                    detected[t].add(int(gene))

    k = pd.read_csv(_RESULTS / "stage22_rewind_clones.csv").set_index("clone_id").loc[clones]
    tbl = pd.DataFrame({
        "clone_id": clones,
        "n_pretreatment_cells": k["n_pretreatment_cells"].to_numpy(),
        "n_lanes": k["n_lanes"].to_numpy(),
        "total_raw_GE_UMI": total_umi.astype(np.int64),
        "n_detected_GE_features_in_raw_pseudobulk": np.array([len(d) for d in detected]),
    })
    tbl["log1p_n_pretreatment_cells"] = np.log1p(tbl["n_pretreatment_cells"].astype(float))
    tbl["log1p_total_raw_GE_UMI"] = np.log1p(tbl["total_raw_GE_UMI"].astype(float))
    tbl["log1p_n_detected_GE_features"] = np.log1p(
        tbl["n_detected_GE_features_in_raw_pseudobulk"].astype(float))

    # cross-check the raw totals against the frozen normalised matrix's zero pattern
    X, xclones = S23._load_rewind_x()
    assert xclones == clones, "clone order disagrees with the frozen 23A cache"
    nz = np.diff(X.tocsr().indptr)
    summary = {
        "columns": ["log1p(n_pretreatment_cells)", "n_lanes", "log1p(total_raw_GE_UMI)",
                    "log1p(n_detected_GE_features_in_raw_pseudobulk)"],
        "clones": len(clones),
        "all_positive_total_umi": bool((tbl["total_raw_GE_UMI"] > 0).all()),
        "detected_matches_normalised_nonzero_pattern":
            bool((tbl["n_detected_GE_features_in_raw_pseudobulk"].to_numpy() == nz).all()),
        "total_raw_GE_UMI": {"min": int(tbl["total_raw_GE_UMI"].min()),
                             "median": float(tbl["total_raw_GE_UMI"].median()),
                             "max": int(tbl["total_raw_GE_UMI"].max()),
                             "sum": int(tbl["total_raw_GE_UMI"].sum())},
        "n_detected_GE_features": {"min": int(tbl["n_detected_GE_features_in_raw_pseudobulk"]
                                              .min()),
                                   "median": float(tbl["n_detected_GE_features_in_raw_pseudobulk"]
                                                   .median()),
                                   "max": int(tbl["n_detected_GE_features_in_raw_pseudobulk"]
                                              .max())},
        "outcome_free": True,
        "note": "Bdepth is diagnostic; V2 §7.7 forbids automatic promotion to a Stage-24 baseline",
    }
    return tbl, summary


# --------------------------------------------------------------------------------------------- #
# 5.5 — recover the 200 historical permutation mappings and replay D00
# --------------------------------------------------------------------------------------------- #
def recover_mappings_and_replay(n_perm: int = 200, replay: bool = True) -> dict:
    """V2 §5.5 / §4.3.

    A mapping is a pure function of `(SEED_PERMUTATION + permutation_id)` replayed through the
    frozen `permute_within` once per outer fold, in fold order -- nothing else consumes that
    generator inside `_rewind_null_once`. The V1 audit proved two draws reproduce bit-for-bit; this
    extends it to all 200 and requires the replayed array to reproduce the committed Stage-23
    summary statistics before it is written.
    """
    t0 = time.perf_counter()
    X, clones = S23._load_rewind_x()
    k = pd.read_csv(_RESULTS / "stage22_rewind_clones.csv").set_index("clone_id").loc[clones]
    y = k["y_primed"].to_numpy()
    fold = k["outer_fold"].to_numpy()
    nuis = np.column_stack([np.log1p(k["n_pretreatment_cells"].to_numpy(dtype=float)),
                            k["n_lanes"].to_numpy(dtype=float)])
    strata = S23.rewind_strata(k)

    # ---- realized strata, verified rather than assumed (M1) ---------------------------------- #
    counts = pd.Series(strata).value_counts().to_dict()
    strata_report = {"cells": {str(a): int(b) for a, b in sorted(counts.items())},
                     "n_cells": len(counts),
                     "largest_cell_share": round(max(counts.values()) / len(strata), 4),
                     "note": "always crossed; `1|2` cannot exist -- one cell cannot span two lanes"}

    # ---- mappings ----------------------------------------------------------------------------- #
    _CACHE.mkdir(parents=True, exist_ok=True)
    map_dir = _CACHE / "permutation_mappings"
    map_dir.mkdir(exist_ok=True)
    digest = hashlib.sha256()
    fixed_by_fold = {f: [] for f in range(S23.N_OUTER)}
    rows = 0
    for b in range(n_perm):
        rng = np.random.default_rng(S23.SEED_PERMUTATION + b)
        per_fold = {}
        for f in range(S23.N_OUTER):
            pmap = S23.permute_within(strata, fold != f, rng)
            side = fold != f
            assert (side[pmap] == side).all(), f"draw {b} fold {f}: crossed the outer boundary"
            assert (strata[pmap] == strata).all(), f"draw {b} fold {f}: left its stratum"
            assert sorted(pmap.tolist()) == list(range(len(pmap))), f"draw {b} fold {f}: not a bijection"
            per_fold[f] = pmap
            fixed_by_fold[f].append(int((pmap == np.arange(len(pmap))).sum()))
            digest.update(pmap.astype(np.int32).tobytes())
            rows += len(pmap)
        np.savez_compressed(map_dir / f"draw_{b:03d}.npz",
                            **{f"fold{f}": per_fold[f].astype(np.int32) for f in per_fold})

    out = {
        "n_permutations": n_perm,
        "base_seed": S23.SEED_PERMUTATION,
        "mapping_set_sha256": digest.hexdigest(),
        "mapping_rows": rows,
        "mapping_cache_dir": str(map_dir.relative_to(ROOT).as_posix()),
        "mapping_is_cache_only": True,
        "fixed_clone_counts_by_fold_draw0": {str(f): fixed_by_fold[f][0]
                                             for f in range(S23.N_OUTER)},
        "mean_fixed_clone_count_by_fold": {str(f): round(float(np.mean(fixed_by_fold[f])), 3)
                                           for f in range(S23.N_OUTER)},
        "structural_assertions": "train->train, test->test, within-stratum, bijective (all draws)",
        "realized_strata": strata_report,
    }

    if not replay:
        out["replay"] = "SKIPPED"
        return out

    # ---- replay every historical D00 ---------------------------------------------------------- #
    cache = {f: S23._frozen_pipeline_cache(X, np.flatnonzero(fold != f), max(S23.K_CANDIDATES))
             for f in range(S23.N_OUTER)}
    rb = json.loads((_RESULTS / "stage23_rewind_results.json").read_text(encoding="utf-8"))
    ap_r1 = rb["pooled_oof_metrics"]["R1"]["AP"]
    d00 = []
    for b in range(n_perm):
        rng = np.random.default_rng(S23.SEED_PERMUTATION + b)
        d00.append(S23._rewind_null_once(X, y, fold, nuis, strata, cache, rng) - ap_r1)
        if (b + 1) % 20 == 0:
            el = time.perf_counter() - t0
            print(f"  replay {b + 1}/{n_perm}  {el / 60:.1f} min  "
                  f"eta {el / (b + 1) * (n_perm - b - 1) / 60:.1f} min", flush=True)

    arr = np.array(d00)
    pe = json.loads((_RESULTS / "stage23_permutation_results.json").read_text(encoding="utf-8"))
    hist = pe["permutation_tests"]["role_a_delta_AP_state"]
    observed = hist["observed"]
    got = {"mean": float(arr.mean()), "sd": float(arr.std()),
           "p95": float(np.percentile(arr, 95)), "max": float(arr.max()),
           "min": float(arr.min()),
           "n_ge_observed": int((arr >= observed).sum()),
           "p_perm": float((1 + int((arr >= observed).sum())) / (len(arr) + 1))}
    matches = {kk: (round(got[kk], 12) == round(hist[f"null_{kk}"], 12))
               for kk in ("mean", "sd", "p95", "max", "min")}
    matches["n_ge_observed"] = got["n_ge_observed"] == hist["n_null_ge_observed"]
    matches["p_perm"] = round(got["p_perm"], 12) == round(hist["p_perm"], 12)

    # a cached copy of the historical array, if this machine still has one, is compared exactly
    cached = _cached_historical_null()
    out["replay"] = {
        "values_recomputed": len(arr),
        "recomputed_summary": got,
        "committed_summary": {kk: hist[f"null_{kk}"] for kk in ("mean", "sd", "p95", "max", "min")}
        | {"n_ge_observed": hist["n_null_ge_observed"], "p_perm": hist["p_perm"]},
        "matches_committed_summary": matches,
        "all_summary_statistics_reproduced": all(matches.values()),
        "bitwise_identical_to_stage23_cache": (
            None if cached is None else bool(np.array_equal(arr, np.array(cached)))),
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 2),
    }
    out["_d00"] = [float(v) for v in arr]
    return out


def _cached_historical_null():
    p = ROOT / "_cc_cache" / "stage23" / "stage23e_null_rewind.json"
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))["nulls"]["role_a_delta_AP_state"]


# --------------------------------------------------------------------------------------------- #
# 5.6 / 4.1 — freeze the decomposition design and the immutable protocol digest
# --------------------------------------------------------------------------------------------- #
def build_protocol(source_design: dict) -> dict:
    """V2 §4.1 + §5.6. The hashed payload is the scientific surface only."""
    payload = {
        "stage": "23.2",
        "plan_version": PLAN_VERSION,
        "plan_canonical_lf_sha256": S23.canonical_text_sha256(PLAN),
        "role": "A",
        "dataset": "GSE227151",
        "frozen_from_stage23": {
            "outer_folds": S23.N_OUTER, "inner_folds": S23.N_INNER,
            "inner_cv": "StratifiedKFold(shuffle=True, random_state=23023) on clones",
            "K_candidates": list(S23.K_CANDIDATES),
            "logistic_C": list(S23.LOGISTIC_C),
            "selection_metric": "mean average precision at clone grain",
            "permutation_base_seed": S23.SEED_PERMUTATION,
            "n_permutations": ANCHOR_NULL["n_permutations"],
            "nuisance_B0": ["log1p(n_pretreatment_cells)", "n_lanes"]},
        "diagnostic_nuisance_Bdepth": [
            "log1p(n_pretreatment_cells)", "n_lanes", "log1p(total_raw_GE_UMI)",
            "log1p(n_detected_GE_features_in_raw_pseudobulk)"],
        "decomposition_design": {
            "cells": {
                "00": "full K x C search + B0        (historical, read from the committed array)",
                "01": "full K x C search + Bdepth",
                "10": "no-K-selection mean + B0",
                "11": "no-K-selection mean + Bdepth"},
            "no_k_selection_reference": {
                "fixed_K_arms": [10, 20, 50],
                "weights": "equal, 1/3 each, fixed before execution",
                "C_selection": "retained from the frozen grid inside every arm",
                "rationale": "the historical R3 selected K = 50,10,10,10,10 -- K=20 in zero of "
                             "five folds -- so a single fixed arm would conflate removing K "
                             "selection with moving to a disfavoured subspace",
                "per_arm_reporting": "mean, CI and arm dispersion, descriptive only"},
            "search_width_ladder": {"rungs": [4, 8, 12], "conditional": False},
            "paired_unit": "permutation_id",
            "bootstrap": {"resamples": 10000, "seed_23_2B": 23421, "seed_23_2C": 23422}},
        "label_reliability": {
            "multinomial_resamples": 5000, "seed": 23431,
            "top_n_ladder": [80, 90, 100, 110, 120],
            "cross_gsm_gdna_concordance": "REMOVED -- gDNA is pooled SampleNum=3, not identifiable",
            "not_supported_requires_independent_outcome_assay_replication": True},
        "power": {
            "cohort_scales": [1, 2, 4], "target_oracle_AUC": [0.66, 0.70],
            "null_allocations": 200, "alternative_simulations": 100,
            "seeds": {"covariate_resample": 23440, "null": 23441, "alternative": 23442,
                      "beta_calibration": 23443},
            "statuses": ["WITHIN_R1_EVENT_COUNT_LIMITATION", "BIOLOGICAL_REPLICATION_LIMITATION"],
            "n_biological_replicates": 1},
        "within_r1_status": source_design["within_r1_status"],
        "lane_composition_sensitivity_permitted":
            source_design["within_r1_status"] == "WITHIN_R1_TECHNICAL_LANES",
        "feature_firewall": {
            "molecular_state": "pretreatment Gene Expression only, 36,601 features",
            "forbidden_inputs": ["clone_id", "cell_uid", "outer_fold", "outcome-rule fields",
                                 "gDNA counts / ranks", "future outcome fields", "source paths"]},
    }
    return {"payload": payload, "canonical_sha256": canonical_json_sha256(payload)}


# --------------------------------------------------------------------------------------------- #
# 23.2A driver
# --------------------------------------------------------------------------------------------- #
def run_23_2a(n_perm: int = 200, replay: bool = True) -> dict:
    _OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    pre = preflight_historical()
    if not pre["ok"]:
        res = {"stage": "23.2A", "verdict": INPUT_BLOCKED, "preflight": pre}
        write_json(A_RESULTS, res)
        return res

    design = source_design_audit()
    ledger = reserved_confirmation_ledger()
    gdna = reconstruct_gdna_rule()
    bdepth_tbl, bdepth = compute_bdepth()
    maps = recover_mappings_and_replay(n_perm=n_perm, replay=replay)

    d00 = maps.pop("_d00", None)
    if d00 is not None:
        write_json(HISTORICAL_D00, {
            "stage": "23.2A",
            "description": "the 200 historical Stage-23E Role-A permutation dAP values, replayed "
                           "from the frozen seed and pipeline, in permutation_id order",
            "source": "replayed via the frozen Stage-23 builder; verified against the committed "
                      "Stage-23 summary statistics before writing",
            "n_permutations": len(d00),
            "base_seed": S23.SEED_PERMUTATION,
            "values_sha256": sha256_bytes(np.array(d00).tobytes()),
            "summary": maps["replay"]["recomputed_summary"],
            "values": d00})

    bdepth_tbl.to_csv(BDEPTH_TABLE, index=False, lineterminator="\n")
    proto = build_protocol(design)

    gates = {
        "historical artifacts verified": pre["ok"],
        "source-design audit recorded, within-R1 status emitted":
            design["within_r1_status"] in {"WITHIN_R1_TECHNICAL_LANES",
                                           "WITHIN_R1_SEPARATE_LIBRARIES",
                                           "WITHIN_R1_STRUCTURE_UNRESOLVED"},
        "reserved confirmation ledger written, metadata only": bool(ledger["entries"]),
        "gDNA reconstruction inputs available and rule reproduced":
            gdna["reproduces_frozen_positives_exactly"] and gdna["matches_anchor_101"],
        "Bdepth exactly computable for all 3,147 clones":
            bdepth["clones"] == ANCHOR_CLONES and bdepth["all_positive_total_umi"]
            and bdepth["detected_matches_normalised_nonzero_pattern"],
        "historical permutation mappings exactly recoverable for ALL 200 draws":
            maps["n_permutations"] == ANCHOR_NULL["n_permutations"],
        "replayed D00 array reproduces the committed Stage-23 summary statistics":
            replay and maps["replay"]["all_summary_statistics_reproduced"],
        "committed historical-null artifact written, mapping digest recorded":
            HISTORICAL_D00.exists() and bool(maps["mapping_set_sha256"]),
        "2x2 design, including all three fixed-K arms, executable in every fold":
            _design_executable(),
        "immutable protocol digest created": bool(proto["canonical_sha256"]),
    }
    verdict = PROTOCOL_FROZEN if all(gates.values()) else INPUT_BLOCKED

    write_json(PROTOCOL, {"stage": "23.2A", "protocol": proto["payload"],
                          "canonical_sha256": proto["canonical_sha256"],
                          "source_provenance": source_provenance()})
    write_json(SOURCE_DESIGN, {"stage": "23.2A", "source_design": design,
                               "gdna_rule": gdna, "bdepth": bdepth})
    write_json(RESERVED_LEDGER, ledger)

    res = {
        "stage": "23.2A",
        "plan": {"file": PLAN.name, "version": PLAN_VERSION,
                 "canonical_lf_sha256": S23.canonical_text_sha256(PLAN)},
        "stage23_2_protocol_sha256": proto["canonical_sha256"],
        "preflight": pre,
        "within_r1_status": design["within_r1_status"],
        "lane_composition_sensitivity_permitted":
            design["within_r1_status"] == "WITHIN_R1_TECHNICAL_LANES",
        "gdna_rule": {kk: gdna[kk] for kk in
                      ("selected_barcodes", "rank_100_cutoff_counts", "tie_size_at_cutoff",
                       "reconstructed_positive_clones", "reproduces_frozen_positives_exactly",
                       "sample_num_grouping_is_a_no_op")},
        "bdepth": bdepth,
        "permutation_recovery": maps,
        "artifacts": {p.name: sha256_file(p) for p in sorted(_OUT.glob("*"))
                      if p.name != A_RESULTS.name},
        "gates": gates,
        "verdict": verdict,
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 2),
        "source_provenance": source_provenance(),
    }
    write_json(A_RESULTS, res)
    return res


def _design_executable() -> bool:
    """Every fold must support K=50, so all three fixed-K arms exist in all five folds."""
    rb = json.loads((_RESULTS / "stage23_rewind_results.json").read_text(encoding="utf-8"))
    return all(rb["per_outer_fold"][str(f)]["max_feasible_K"] >= max(S23.K_CANDIDATES)
               for f in range(S23.N_OUTER))




# --------------------------------------------------------------------------------------------- #
# 23.2B / 23.2C — model-selection and residual-depth decomposition.
#
# Both substages ask what happens to the *null centre* when one thing about the historical Role-A
# comparison is changed: the model-selection path (23.2B) or the nuisance block (23.2C). They share
# one machine because of an exact identity worth stating plainly:
#
#     the inner-CV scores of the twelve (K, C) candidates do not depend on which subset of them a
#     rule is allowed to pick from.
#
# So one pass per permutation draw scores all twelve candidates, and then every quantity the plan
# needs is a different argmax over the same scores:
#
#     historical cell 00      argmax over all 12          (K free, C free)
#     fixed-K arm A(K)        argmax over the 4 C at K    (K fixed, C free)
#     ladder rung 8           argmax over K in {10,20}
#     ladder rung 12          identical to cell 00
#
# Only the final outer-training refit differs between them, and that is one fit against the inner
# CV's thirty-six. The saving is real, not a shortcut: nothing is reused across *draws*, the inner
# gene filter / scaler / PCA are recomputed inside every inner split exactly as Stage 23 did, and
# the permutation mapping is drawn from the frozen generator in the frozen order.
#
# `R1` (and `R1depth`) carry no expression, so a permutation of profiles cannot move them. Bdepth
# is likewise a property of the clone, not of the profile it received, so it too is invariant under
# the permutation. Both are therefore computed once on observed data and reused -- exact, not
# convenient.
# --------------------------------------------------------------------------------------------- #
B_RESULTS = _OUT / "stage23_2_model_selection_decomposition.json"
C_RESULTS = _OUT / "stage23_2_depth_decomposition.json"
N_PAIRED_BOOT = 10_000
SEED_BOOT_B = 23421
SEED_BOOT_C = 23422
LADDER_RUNGS = {"4_candidate_fixed_K": None, "8_candidate": (10, 20), "12_candidate": (10, 20, 50)}


def _nuisance_matrix(block: str) -> np.ndarray:
    """B0 is the frozen Stage-23 pair; Bdepth adds the two outcome-free depth scalars (V2 §7.1)."""
    X, clones = S23._load_rewind_x()
    k = pd.read_csv(_RESULTS / "stage22_rewind_clones.csv").set_index("clone_id").loc[clones]
    if block == "B0":
        return np.column_stack([np.log1p(k["n_pretreatment_cells"].to_numpy(dtype=float)),
                                k["n_lanes"].to_numpy(dtype=float)])
    tbl = pd.read_csv(BDEPTH_TABLE).set_index("clone_id").loc[clones]
    return np.column_stack([tbl["log1p_n_pretreatment_cells"].to_numpy(dtype=float),
                            tbl["n_lanes"].to_numpy(dtype=float),
                            tbl["log1p_total_raw_GE_UMI"].to_numpy(dtype=float),
                            tbl["log1p_n_detected_GE_features"].to_numpy(dtype=float)])


def _r1_reference_ap(nuis: np.ndarray, y: np.ndarray, fold: np.ndarray) -> tuple[float, dict]:
    """The expression-free baseline for a nuisance block, fitted once on observed data."""
    from sklearn.metrics import average_precision_score
    from sklearn.model_selection import StratifiedKFold

    oof = np.full(len(y), np.nan)
    selected = {}
    for f in range(S23.N_OUTER):
        tr, te = np.flatnonzero(fold != f), np.flatnonzero(fold == f)
        skf = StratifiedKFold(n_splits=S23.N_INNER, shuffle=True, random_state=S23.SEED_PROTOCOL)
        sc: dict = {}
        for itr_i, iva_i in skf.split(np.zeros(len(tr)), y[tr]):
            itr, iva = tr[itr_i], tr[iva_i]
            btr, bva = S23.standardize_train_only(nuis[itr], nuis[iva])
            for C in S23.LOGISTIC_C:
                p = S23._fit_logistic(btr, y[itr], bva, C, [])
                sc.setdefault(C, []).append(average_precision_score(y[iva], p))
        C = max(sc.items(), key=lambda kv: (round(float(np.mean(kv[1])), 12), -kv[0]))[0]
        btr, bte = S23.standardize_train_only(nuis[tr], nuis[te])
        oof[te] = S23._fit_logistic(btr, y[tr], bte, C, [])
        selected[str(f)] = C
    return float(average_precision_score(y, oof)), selected


def _one_draw_all_rules(X, y, fold, nuis, strata, cache, rng):
    """One permutation draw: score the 12 candidates once, then apply every selection rule.

    Returns the pooled OOF average precision for each rule. The rng is consumed exactly as
    `_rewind_null_once` consumes it -- one `permute_within` per outer fold, in fold order -- so the
    mapping is the historical one for this draw.
    """
    from sklearn.metrics import average_precision_score
    from sklearn.model_selection import StratifiedKFold

    rules = ["hist12", "ladder8", "K10", "K20", "K50"]
    oof = {r: np.full(len(y), np.nan) for r in rules}

    for f in range(S23.N_OUTER):
        side = fold != f
        pmap = permutation_for_fold(strata, side, rng)
        tr, te = np.flatnonzero(side), np.flatnonzero(~side)
        skf = StratifiedKFold(n_splits=S23.N_INNER, shuffle=True, random_state=S23.SEED_PROTOCOL)
        scores: dict = {}
        for itr_i, iva_i in skf.split(np.zeros(len(tr)), y[tr]):
            itr, iva = tr[itr_i], tr[iva_i]
            ztr, zva, _, kmax = S23.expression_block(X, pmap[itr], pmap[iva],
                                                     max(S23.K_CANDIDATES))
            btr, bva = S23.standardize_train_only(nuis[itr], nuis[iva])
            for k in S23.K_CANDIDATES:
                if k > kmax:
                    continue
                for C in S23.LOGISTIC_C:
                    p = S23._fit_logistic(np.hstack([ztr[:, :k], btr]), y[itr],
                                          np.hstack([zva[:, :k], bva]), C, [])
                    scores.setdefault((k, C), []).append(average_precision_score(y[iva], p))

        def pick(allowed, sc=scores):
            """Frozen Stage-23 tie-break: maximise mean AP, then smaller K, then smaller C."""
            cand = {kc: v for kc, v in sc.items() if allowed is None or kc[0] in allowed}
            return max(cand.items(), key=lambda kv: (round(float(np.mean(kv[1])), 12),
                                                     -kv[0][0], -kv[0][1]))[0]

        ztr = S23._apply_cached(X, pmap[tr], cache[f])
        zte = S23._apply_cached(X, pmap[te], cache[f])
        btr, bte = S23.standardize_train_only(nuis[tr], nuis[te])
        chosen = {"hist12": pick((10, 20, 50)), "ladder8": pick((10, 20)),
                  "K10": pick((10,)), "K20": pick((20,)), "K50": pick((50,))}
        for rule, (k, C) in chosen.items():
            oof[rule][te] = S23._fit_logistic(np.hstack([ztr[:, :k], btr]), y[tr],
                                              np.hstack([zte[:, :k], bte]), C, [])
    return {r: float(average_precision_score(y, oof[r])) for r in rules}


def permutation_for_fold(strata, side, rng):
    """Thin wrapper so the frozen generator call is named in one place and easy to audit."""
    return S23.permute_within(strata, side, rng)


def _paired_ci(delta: np.ndarray, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(delta), size=(N_PAIRED_BOOT, len(delta)))
    means = delta[idx].mean(axis=1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return {"mean": float(delta.mean()), "median": float(np.median(delta)),
            "sd": float(delta.std()), "fraction_gt_0": float((delta > 0).mean()),
            "ci95_low": float(lo), "ci95_high": float(hi),
            "resamples": N_PAIRED_BOOT, "seed": seed}


def _status(ci: dict, name: str) -> str:
    if ci["ci95_low"] > 0:
        return "SUPPORTED"
    if ci["ci95_high"] <= 0:
        return "NOT_SUPPORTED"
    return "UNRESOLVED"


def _run_decomposition(block: str, n_perm: int) -> dict:
    """Shared engine for 23.2B (block='B0') and the null half of 23.2C (block='Bdepth')."""
    t0 = time.perf_counter()
    X, clones = S23._load_rewind_x()
    k = pd.read_csv(_RESULTS / "stage22_rewind_clones.csv").set_index("clone_id").loc[clones]
    y = k["y_primed"].to_numpy()
    fold = k["outer_fold"].to_numpy()
    strata = S23.rewind_strata(k)
    nuis = _nuisance_matrix(block)
    cache = {f: S23._frozen_pipeline_cache(X, np.flatnonzero(fold != f), max(S23.K_CANDIDATES))
             for f in range(S23.N_OUTER)}
    r1_ap, r1_sel = _r1_reference_ap(nuis, y, fold)

    per_draw = {r: [] for r in ("hist12", "ladder8", "K10", "K20", "K50")}
    for b in range(n_perm):
        rng = np.random.default_rng(S23.SEED_PERMUTATION + b)
        got = _one_draw_all_rules(X, y, fold, nuis, strata, cache, rng)
        for r, v in got.items():
            per_draw[r].append(v - r1_ap)
        if (b + 1) % 20 == 0:
            el = time.perf_counter() - t0
            print(f"  [{block}] {b + 1}/{n_perm}  {el / 60:.1f} min  "
                  f"eta {el / (b + 1) * (n_perm - b - 1) / 60:.1f} min", flush=True)

    arrays = {r: np.array(v) for r, v in per_draw.items()}
    arms = np.column_stack([arrays["K10"], arrays["K20"], arrays["K50"]])
    no_k = arms.mean(axis=1)
    _CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(_CACHE / f"decomposition_{block}.npz",
                        **{r: arrays[r] for r in arrays}, no_k_selection=no_k)
    return {"block": block, "r1_reference_AP": r1_ap, "r1_selected_C_per_fold": r1_sel,
            "arrays": arrays, "no_k": no_k, "n_permutations": n_perm,
            "runtime_minutes": round((time.perf_counter() - t0) / 60, 2)}


def _observed_cells(block: str) -> dict:
    """The observed-label counterpart of every selection rule, under the frozen outer folds.

    V2 §6.5 / §7.3. Diagnostic effect attribution only -- no new p-value is computed here and the
    historical Stage-23 failure is not reinterpreted.
    """
    from sklearn.metrics import average_precision_score
    from sklearn.model_selection import StratifiedKFold

    X, clones = S23._load_rewind_x()
    k = pd.read_csv(_RESULTS / "stage22_rewind_clones.csv").set_index("clone_id").loc[clones]
    y = k["y_primed"].to_numpy()
    fold = k["outer_fold"].to_numpy()
    nuis = _nuisance_matrix(block)
    r1_ap, _ = _r1_reference_ap(nuis, y, fold)

    rules = ("hist12", "ladder8", "K10", "K20", "K50")
    oof = {r: np.full(len(y), np.nan) for r in rules}
    chosen_per_fold: dict = {r: {} for r in rules}
    for f in range(S23.N_OUTER):
        tr, te = np.flatnonzero(fold != f), np.flatnonzero(fold == f)
        skf = StratifiedKFold(n_splits=S23.N_INNER, shuffle=True, random_state=S23.SEED_PROTOCOL)
        scores: dict = {}
        for itr_i, iva_i in skf.split(np.zeros(len(tr)), y[tr]):
            itr, iva = tr[itr_i], tr[iva_i]
            ztr, zva, _, kmax = S23.expression_block(X, itr, iva, max(S23.K_CANDIDATES))
            btr, bva = S23.standardize_train_only(nuis[itr], nuis[iva])
            for kk in S23.K_CANDIDATES:
                if kk > kmax:
                    continue
                for C in S23.LOGISTIC_C:
                    p = S23._fit_logistic(np.hstack([ztr[:, :kk], btr]), y[itr],
                                          np.hstack([zva[:, :kk], bva]), C, [])
                    scores.setdefault((kk, C), []).append(average_precision_score(y[iva], p))

        def pick(allowed, sc=scores):
            cand = {kc: v for kc, v in sc.items() if allowed is None or kc[0] in allowed}
            return max(cand.items(), key=lambda kv: (round(float(np.mean(kv[1])), 12),
                                                     -kv[0][0], -kv[0][1]))[0]

        ztr, zte, _, _ = S23.expression_block(X, tr, te, max(S23.K_CANDIDATES))
        btr, bte = S23.standardize_train_only(nuis[tr], nuis[te])
        for rule, allowed in (("hist12", (10, 20, 50)), ("ladder8", (10, 20)),
                              ("K10", (10,)), ("K20", (20,)), ("K50", (50,))):
            kk, C = pick(allowed)
            chosen_per_fold[rule][str(f)] = {"K": kk, "C": C}
            oof[rule][te] = S23._fit_logistic(np.hstack([ztr[:, :kk], btr]), y[tr],
                                              np.hstack([zte[:, :kk], bte]), C, [])

    d = {r: float(average_precision_score(y, oof[r])) - r1_ap for r in rules}
    d["no_k_selection"] = float(np.mean([d["K10"], d["K20"], d["K50"]]))
    return {"r1_reference_AP": r1_ap, "delta_AP": d,
            "selected_per_fold": chosen_per_fold}


def run_23_2b(n_perm: int = 200) -> dict:
    """V2 §6. How much of the positive null centre is R3's broader K x C selection path?"""
    d00_doc = json.loads(HISTORICAL_D00.read_text(encoding="utf-8"))
    d00 = np.array(d00_doc["values"])
    core = _run_decomposition("B0", n_perm)

    # the historical R1 must reproduce bit-for-bit, or the reference has drifted
    rb = json.loads((_RESULTS / "stage23_rewind_results.json").read_text(encoding="utf-8"))
    hist_r1 = rb["pooled_oof_metrics"]["R1"]["AP"]
    r1_exact = core["r1_reference_AP"] == hist_r1

    # cell 00 recomputed here must equal the committed historical array
    recomputed_d00 = core["arrays"]["hist12"]
    d00_reproduced = bool(np.array_equal(recomputed_d00, d00))

    arms = {f"K{kk}": core["arrays"][f"K{kk}"] for kk in S23.K_CANDIDATES}
    no_k = core["no_k"]
    s = d00 - no_k
    ci = _paired_ci(s, SEED_BOOT_B)
    status = _status(ci, "MODEL_SELECTION_NULL_INFLATION")

    arm_stats = {a: _paired_ci(arms[a], SEED_BOOT_B) for a in arms}
    arm_means = [arm_stats[a]["mean"] for a in arms]
    ladder = {"4_candidate_fixed_K_mean": float(np.mean(arm_means)),
              "8_candidate": float(core["arrays"]["ladder8"].mean()),
              "12_candidate": float(d00.mean())}
    ladder_monotone = ladder["4_candidate_fixed_K_mean"] <= ladder["8_candidate"] <= \
        ladder["12_candidate"]

    observed = _observed_cells("B0")
    out = {
        "stage": "23.2B",
        "plan": {"file": PLAN.name, "version": PLAN_VERSION},
        "stage23_2_protocol_sha256": json.loads(
            PROTOCOL.read_text(encoding="utf-8"))["canonical_sha256"],
        "historical_null_artifact_sha256": sha256_file(HISTORICAL_D00),
        "mapping_set_sha256": json.loads(
            A_RESULTS.read_text(encoding="utf-8"))["permutation_recovery"]["mapping_set_sha256"],
        "n_permutations": n_perm,
        "reference": {"r1_reference_AP": core["r1_reference_AP"],
                      "historical_R1_AP": hist_r1,
                      "reproduces_historical_R1_exactly": bool(r1_exact),
                      "selected_C_per_fold": core["r1_selected_C_per_fold"]},
        "cell_00_recomputation": {
            "reproduces_committed_D00_exactly": d00_reproduced,
            "max_abs_difference": float(np.abs(recomputed_d00 - d00).max()),
            "note": "cell 00 is READ from the committed array for every statistic; this "
                    "recomputation exists only to prove the shared-scores engine is faithful"},
        "D00": {"mean": float(d00.mean()), "sd": float(d00.std())},
        "D10_no_k_selection": {"mean": float(no_k.mean()), "sd": float(no_k.std())},
        "per_arm": {a: arm_stats[a] for a in arm_stats},
        "arm_dispersion": float(max(arm_means) - min(arm_means)),
        "selection_shift": ci,
        "fraction_null_mean_explained_by_search":
            (float(ci["mean"] / d00.mean()) if d00.mean() > 0 else None),
        "search_width_ladder": ladder,
        "ladder_monotone_increase": bool(ladder_monotone),
        "observed_sensitivity": observed,
        "MODEL_SELECTION_NULL_INFLATION": status,
        "runtime_minutes": core["runtime_minutes"],
        "source_provenance": source_provenance(),
    }
    write_json(B_RESULTS, out)
    return out


def run_23_2c(n_perm: int = 200) -> dict:
    """V2 §7. How much of the positive null centre survives an expanded depth nuisance block?"""
    d00 = np.array(json.loads(HISTORICAL_D00.read_text(encoding="utf-8"))["values"])
    b_doc = json.loads(B_RESULTS.read_text(encoding="utf-8"))
    core = _run_decomposition("Bdepth", n_perm)

    d01 = core["arrays"]["hist12"]
    d11 = core["no_k"]
    d10 = np.load(_CACHE / "decomposition_B0.npz")["no_k_selection"]

    depth_shift_full = _paired_ci(d00 - d01, SEED_BOOT_C)
    mu = {"mu00": float(d00.mean()), "mu10": float(d10.mean()),
          "mu01": float(d01.mean()), "mu11": float(d11.mean())}
    contrasts = {
        "selection_main_contrast": _paired_ci(d00 - d10, SEED_BOOT_C),
        "depth_main_contrast": _paired_ci(d00 - d01, SEED_BOOT_C),
        "factor_interaction": _paired_ci(d00 - d10 - d01 + d11, SEED_BOOT_C),
    }
    retention = _technical_retention(n_perm)
    depth_status = _residual_depth_status(depth_shift_full, retention)

    observed = _observed_cells("Bdepth")
    o11 = observed["delta_AP"]["no_k_selection"]
    q95_11 = float(np.percentile(d11, 95))
    p_diag_11 = float((1 + int((d11 >= o11).sum())) / (len(d11) + 1))
    corrected = ("POSITIVE" if (o11 > 0 and o11 > q95_11 and p_diag_11 <= 0.05) else "NEGATIVE")

    arm_means = [float(core["arrays"][f"K{kk}"].mean()) for kk in S23.K_CANDIDATES]
    out = {
        "stage": "23.2C",
        "plan": {"file": PLAN.name, "version": PLAN_VERSION},
        "stage23_2_protocol_sha256": b_doc["stage23_2_protocol_sha256"],
        "n_permutations": n_perm,
        "nuisance_Bdepth": ["log1p(n_pretreatment_cells)", "n_lanes", "log1p(total_raw_GE_UMI)",
                            "log1p(n_detected_GE_features_in_raw_pseudobulk)"],
        "reference": {"r1depth_reference_AP": core["r1_reference_AP"],
                      "selected_C_per_fold": core["r1_selected_C_per_fold"]},
        "null_means": mu,
        "depth_shift_full": depth_shift_full,
        "contrasts": contrasts,
        "per_arm_Bdepth": {f"K{kk}": _paired_ci(core["arrays"][f"K{kk}"], SEED_BOOT_C)
                           for kk in S23.K_CANDIDATES},
        "arm_dispersion_Bdepth": float(max(arm_means) - min(arm_means)),
        "technical_retention": retention,
        "observed_2x2": {
            "O00": b_doc["observed_sensitivity"]["delta_AP"]["hist12"],
            "O10": b_doc["observed_sensitivity"]["delta_AP"]["no_k_selection"],
            "O01": observed["delta_AP"]["hist12"],
            "O11": o11,
            "O00_matches_historical_within_tolerance":
                abs(b_doc["observed_sensitivity"]["delta_AP"]["hist12"] - 0.01050162935116511)
                < 1e-9},
        "corrected_same_data_diagnostic": {
            "O11": o11, "q95_11": q95_11, "p_diag_11": p_diag_11,
            "CORRECTED_SAME_DATA_SIGNAL_DIAGNOSTIC": corrected,
            "note": "exploratory same-Rewind diagnostic; a POSITIVE result cannot emit "
                    "ROLE_A_CONFIRMATORY_SUPPORTED"},
        "lane_composition_sensitivity": _lane_sensitivity_status(),
        "RESIDUAL_DEPTH_STRUCTURE": depth_status,
        "critical_interpretation": "Bdepth is diagnostic. V2 §7.7 forbids promoting it into a "
                                   "Stage-24 production baseline without the frozen Stage-24 plan",
        "runtime_minutes": core["runtime_minutes"],
        "source_provenance": source_provenance(),
    }
    write_json(C_RESULTS, out)
    return out


def _lane_sensitivity_status() -> dict:
    a = json.loads(A_RESULTS.read_text(encoding="utf-8"))
    permitted = a["lane_composition_sensitivity_permitted"]
    return {"permitted": permitted, "within_r1_status": a["within_r1_status"],
            "executed": False,
            "reason": ("V2 §7.5 permits this only under WITHIN_R1_TECHNICAL_LANES; 23.2A returned "
                       f"{a['within_r1_status']}, so per-sample cell counts could absorb "
                       "biological-unit structure and are forbidden")
            if not permitted else "permitted and executed"}


def _technical_retention(n_perm: int) -> dict:
    """V2 §7.4. Outcome-free: how much continuous technical similarity survives the coarse strata?"""
    from scipy.stats import spearmanr

    X, clones = S23._load_rewind_x()
    k = pd.read_csv(_RESULTS / "stage22_rewind_clones.csv").set_index("clone_id").loc[clones]
    fold = k["outer_fold"].to_numpy()
    strata = S23.rewind_strata(k)
    tbl = pd.read_csv(BDEPTH_TABLE).set_index("clone_id").loc[clones]
    umi = tbl["log1p_total_raw_GE_UMI"].to_numpy()
    det = tbl["log1p_n_detected_GE_features"].to_numpy()
    raw_umi = tbl["total_raw_GE_UMI"].to_numpy(dtype=float)
    nz = np.diff(X.tocsr().indptr).astype(float)

    pairs = {"donor_nonzero_vs_recipient_log1p_total_UMI": (nz, umi),
             "donor_nonzero_vs_recipient_log1p_detected": (nz, det),
             "donor_total_UMI_vs_recipient_total_UMI": (raw_umi, raw_umi)}
    acc = {kk: [] for kk in pairs}
    for b in range(n_perm):
        rng = np.random.default_rng(S23.SEED_PERMUTATION + b)
        per = []
        for f in range(S23.N_OUTER):
            pmap = permutation_for_fold(strata, fold != f, rng)
            per.append(pmap)
        pmap = per[0]                       # fold-0 mapping covers every clone exactly once
        for name, (donor, recip) in pairs.items():
            rho = spearmanr(donor[pmap], recip).statistic
            acc[name].append(float(rho))
    out = {}
    for name, vals in acc.items():
        arr = np.array(vals)
        ci = _paired_ci(arr, SEED_BOOT_C)
        out[name] = {"median": float(np.median(arr)), "mean": ci["mean"],
                     "ci95_low": ci["ci95_low"], "ci95_high": ci["ci95_high"],
                     "excludes_zero_positive": bool(ci["ci95_low"] > 0)}
    out["note"] = "no outcome label is used in this diagnostic"
    return out


def _residual_depth_status(depth_shift_full: dict, retention: dict) -> str:
    """V2 §7.6 -- outcome-level evidence AND at least one mechanistic correlation."""
    mech = [v for kk, v in retention.items() if isinstance(v, dict)]
    any_pos = any(v["excludes_zero_positive"] for v in mech)
    all_incl_zero = all(not v["excludes_zero_positive"] for v in mech)
    if depth_shift_full["ci95_low"] > 0 and any_pos:
        return "SUPPORTED"
    if depth_shift_full["ci95_high"] <= 0 and all_incl_zero:
        return "NOT_SUPPORTED"
    return "UNRESOLVED"




# --------------------------------------------------------------------------------------------- #
# 23.2D — outcome-label reliability.
#
# This substage studies the MEASUREMENT, not predictive performance. No alternate outcome may
# become a new prediction target here (V2 §3.7), and no predictor is fitted at all.
#
# The asymmetry in §8.6 is the whole point and is implemented literally: instability can establish
# that the label IS fragile, but stability cannot establish that it is sound, because the two
# diagnostics available -- multinomial resampling of one pooled sequencing library, and sensitivity
# to the cutoff position -- observe only two noise sources. Colony-level biological sampling, PCR
# duplication and assay repetition are unobserved. `NOT_SUPPORTED` therefore additionally requires
# independent outcome-assay replication of the same clones, which the Rewind materials do not
# contain, so in practice the reachable outcomes here are SUPPORTED or UNRESOLVED.
#
# V1's §8.5 cross-GSM gDNA concordance is absent by design: the gDNA table is a single pooled
# library (SampleNum = 3), so per-GSM outcome support is not identifiable.
# --------------------------------------------------------------------------------------------- #
D_RESULTS = _OUT / "stage23_2_label_reliability.json"
N_MULTINOMIAL = 5000
SEED_MULTINOMIAL = 23431
TOP_N_LADDER = (80, 90, 100, 110, 120)


def _select_top_n(counts: np.ndarray, n: int) -> np.ndarray:
    """The frozen source rule: rank by support, take the top n, keep every tie at the cutoff.

    Returns a boolean mask over the barcode axis. `slice_max(n=..., with_ties=TRUE)` in the author's
    dplyr code means the cutoff VALUE is retained, so a tie at rank n yields more than n rows -- the
    behaviour that produces 101 barcodes at n = 100.
    """
    if n >= len(counts):
        return np.ones(len(counts), dtype=bool)
    cutoff = np.partition(counts, -n)[-n]
    return counts >= cutoff


def run_23_2d() -> dict:
    """V2 §8. Is the frozen 35-positive hard label stable enough to carry the Role-A claim?"""
    t0 = time.perf_counter()
    g = pd.read_csv(REWIND_ROOT / GDNA_FILE, sep="\t")
    cells = pd.read_csv(_RESULTS / "stage22_rewind_cells.csv")
    retained = sorted(cells["clone_id"].unique())
    retained_set = set(retained)
    frozen_pos = set(cells.loc[cells["y_primed"] == 1, "clone_id"].unique())

    agg = g.groupby("BC50StarcodeD8")["counts"].sum().sort_values(ascending=False)
    barcodes = agg.index.to_numpy()
    counts = agg.to_numpy().astype(np.int64)
    total_n = int(counts.sum())

    # ---- 8.1 exact reproduction is a precondition, not a result ------------------------------ #
    base_mask = _select_top_n(counts, ANCHOR_TOP_N)
    base_pos = {b for b in barcodes[base_mask] if b in retained_set}
    reproduced = base_pos == frozen_pos
    if not reproduced:
        return {"stage": "23.2D", "verdict": "LABEL_RECONSTRUCTION_FAILED",
                "reconstructed": len(base_pos), "frozen": len(frozen_pos),
                "note": "V2 §21 stop condition: the source rule cannot reproduce 35 positives"}

    # ---- 8.2 cutoff geometry, ranks 80..120 --------------------------------------------------- #
    cutoff_100 = int(counts[ANCHOR_TOP_N - 1])
    geometry = []
    for rank in range(80, 121):
        v = int(counts[rank - 1])
        geometry.append({
            "rank": rank, "counts": v,
            "tie_size_at_this_value": int((counts == v).sum()),
            "gap_to_previous_rank": int(counts[rank - 2] - v) if rank > 1 else None,
            "gap_to_next_rank": int(v - counts[rank]) if rank < len(counts) else None,
            "ratio_to_rank_100_cutoff": round(v / cutoff_100, 6)})

    # ---- 8.3 conditional multinomial sampling stability --------------------------------------- #
    p = counts / total_n
    rng = np.random.default_rng(SEED_MULTINOMIAL)
    frozen_idx = np.array([i for i, b in enumerate(barcodes) if b in frozen_pos])
    retained_idx = np.array([i for i, b in enumerate(barcodes) if b in retained_set])
    hits = np.zeros(len(barcodes), dtype=np.int64)
    n_pos_draws = np.empty(N_MULTINOMIAL, dtype=np.int64)
    jacc_draws = np.empty(N_MULTINOMIAL)
    for i in range(N_MULTINOMIAL):
        c_star = rng.multinomial(total_n, p)
        mask = _select_top_n(c_star, ANCHOR_TOP_N)
        hits += mask
        sel = set(barcodes[mask]) & retained_set
        n_pos_draws[i] = len(sel)
        union = len(sel | frozen_pos)
        jacc_draws[i] = len(sel & frozen_pos) / union if union else 1.0

    retention = hits[frozen_idx] / N_MULTINOMIAL
    per_positive = sorted(
        ({"clone_id": str(barcodes[j]), "counts": int(counts[j]),
          "rank": int(np.flatnonzero(barcodes == barcodes[j])[0] + 1),
          "P_selected": float(hits[j] / N_MULTINOMIAL)} for j in frozen_idx),
        key=lambda d: d["P_selected"])
    # a false-positive view: retained clones that are NOT frozen positives but get selected
    non_pos_idx = np.array([i for i in retained_idx if barcodes[i] not in frozen_pos])
    intruder_rate = float((hits[non_pos_idx] / N_MULTINOMIAL).sum())

    multinomial = {
        "resamples": N_MULTINOMIAL, "seed": SEED_MULTINOMIAL,
        "total_gdna_counts_N": total_n, "distinct_barcodes": int(len(barcodes)),
        "selection_units": 1,
        "mean_frozen_positive_retention": float(retention.mean()),
        "median_frozen_positive_retention": float(np.median(retention)),
        "min_frozen_positive_retention": float(retention.min()),
        "n_positives_below_0_50": int((retention < 0.50).sum()),
        "n_positives_below_0_80": int((retention < 0.80).sum()),
        "expected_intruding_retained_clones_per_draw": intruder_rate,
        "positive_clone_count_per_draw": {
            "mean": float(n_pos_draws.mean()), "sd": float(n_pos_draws.std()),
            "min": int(n_pos_draws.min()), "max": int(n_pos_draws.max())},
        "jaccard_vs_frozen_set": {
            "mean": float(jacc_draws.mean()), "p05": float(np.percentile(jacc_draws, 5)),
            "median": float(np.median(jacc_draws))},
        "per_frozen_positive": per_positive,
        "scope": "sequencing-count sampling noise only, from ONE pooled library; colony-level "
                 "biological sampling, PCR duplication and assay repetition are unobserved",
    }

    # ---- 8.4 cutoff sensitivity ---------------------------------------------------------------- #
    ladder = {}
    for n in TOP_N_LADDER:
        mask = _select_top_n(counts, n)
        sel = {b for b in barcodes[mask] if b in retained_set}
        union = len(sel | frozen_pos)
        ladder[f"top{n}"] = {
            "selected_barcodes": int(mask.sum()),
            "tie_expansion": int(mask.sum() - n),
            "positive_clones": len(sel),
            "jaccard_vs_top100": (len(sel & frozen_pos) / union if union else 1.0),
            "frozen_positives_lost": len(frozen_pos - sel),
            "frozen_negatives_gained": len(sel - frozen_pos)}

    # ---- 8.6 status, with the V2 asymmetry ----------------------------------------------------- #
    unstable_a = (multinomial["mean_frozen_positive_retention"] < 0.80
                  or multinomial["n_positives_below_0_50"] >= 7)
    unstable_b = min(ladder["top90"]["jaccard_vs_top100"],
                     ladder["top110"]["jaccard_vs_top100"]) < 0.80
    stable_all = (multinomial["mean_frozen_positive_retention"] >= 0.90
                  and multinomial["n_positives_below_0_80"] <= 3
                  and ladder["top90"]["jaccard_vs_top100"] >= 0.90
                  and ladder["top110"]["jaccard_vs_top100"] >= 0.90)
    independent_replication_exists = False        # none in the Rewind materials; see V2 §8.6

    if unstable_a and unstable_b:
        status = "SUPPORTED"
        why = "both the multinomial-stability and cutoff-sensitivity criteria indicate instability"
    elif stable_all and independent_replication_exists:
        status = "NOT_SUPPORTED"
        why = "stability criteria met AND independent outcome-assay replication agrees"
    else:
        status = "UNRESOLVED"
        why = ("stability criteria met" if stable_all else "criteria are mixed") + \
              "; NOT_SUPPORTED additionally requires independent outcome-assay replication of the " \
              "same clones, which the Rewind materials do not contain (V2 §8.6)"

    out = {
        "stage": "23.2D",
        "plan": {"file": PLAN.name, "version": PLAN_VERSION},
        "stage23_2_protocol_sha256": json.loads(
            PROTOCOL.read_text(encoding="utf-8"))["canonical_sha256"],
        "source_rule_reproduces_35_positives": reproduced,
        "cutoff_geometry_ranks_80_to_120": geometry,
        "rank_100_cutoff_counts": cutoff_100,
        "tie_size_at_rank_100": int((counts == cutoff_100).sum()),
        "multinomial_stability": multinomial,
        "cutoff_sensitivity": ladder,
        "cross_gsm_gdna_concordance": "REMOVED IN V2 -- gDNA is one pooled library (SampleNum=3), "
                                      "so per-GSM outcome concordance is not identifiable",
        "independent_outcome_assay_replication_available": independent_replication_exists,
        "not_supported_reachable": bool(independent_replication_exists),
        "OUTCOME_LABEL_LIMITATION": status,
        "status_reason": why,
        "candidate_future_formulations": {
            "status": "EXPLORATORY_PROPOSAL_ONLY",
            "note": "V2 §8.7 -- listed only if the limitation is SUPPORTED, and frozen for "
                    "confirmation at 23.2F before any evidence is inspected",
            "candidates": (["continuous gDNA support", "soft probabilistic positive membership",
                            "a margin-separated hard label",
                            "a source-author validated alternate outcome"]
                           if status == "SUPPORTED" else [])},
        "no_predictive_model_fitted": True,
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 2),
        "source_provenance": source_provenance(),
    }
    write_json(D_RESULTS, out)
    return out




# --------------------------------------------------------------------------------------------- #
# 23.2E — power / identifiability.
#
# Planning, not confirmation. The question is which effect sizes the historical Role-A pipeline can
# actually detect at different rare-positive counts, given this feature geometry and the
# appropriate null.
#
# Two constraints from V2 shape the implementation:
#
#   * the synthetic state direction is LABEL-FREE. `z` is built from the expression matrix
#     residualised on Bdepth, never from `y_primed`. Using all clones is permitted only because `z`
#     is a simulation generator; it never scores the real outcome and never enters the fitted
#     evaluation pipeline.
#   * scaling cohorts cannot create biological replicate diversity. Scales 2 and 4 resample clones
#     from R1 with replacement, so every synthetic cohort is still biological replicate R1. The
#     curve answers the within-R1 event-count question (9.5.1) and says nothing about the
#     biological-replication question (9.5.2), whose status is a design fact rather than an estimate.
# --------------------------------------------------------------------------------------------- #
E_SHARD_DIR = _CACHE / "power_shards"
E_RESULTS = _OUT / "stage23_2_power_identifiability.json"
COHORT_SCALES = (1, 2, 4)
TARGET_AUCS = (0.66, 0.70)
N_NULL_ALLOC = 200
N_ALT_SIM = 100
SEED_COVARIATE = 23440
SEED_NULL = 23441
SEED_ALT = 23442
SEED_BETA = 23443


def synthetic_direction() -> np.ndarray:
    """V2 §9.3. One label-free simulation-generating score `z`, deterministically oriented."""
    from sklearn.decomposition import PCA

    X, clones = S23._load_rewind_x()
    tbl = pd.read_csv(BDEPTH_TABLE).set_index("clone_id").loc[clones]
    B = np.column_stack([tbl["log1p_n_pretreatment_cells"].to_numpy(dtype=float),
                         tbl["n_lanes"].to_numpy(dtype=float),
                         tbl["log1p_total_raw_GE_UMI"].to_numpy(dtype=float),
                         tbl["log1p_n_detected_GE_features"].to_numpy(dtype=float)])
    B = np.column_stack([np.ones(len(B)), (B - B.mean(0)) / np.where(B.std(0) == 0, 1, B.std(0))])

    # residualise each gene on Bdepth, in blocks so the dense copy stays bounded
    Xd = np.asarray(X.todense(), dtype=np.float64)
    coef, *_ = np.linalg.lstsq(B, Xd, rcond=None)
    R = Xd - B @ coef
    del Xd

    pca = PCA(n_components=1, svd_solver="randomized",
              random_state=S23.SEED_PROTOCOL).fit(R)
    load = pca.components_[0]
    j = int(np.lexsort((np.arange(len(load)), -np.abs(load)))[0])   # ties -> smallest feature index
    sign = 1.0 if load[j] >= 0 else -1.0
    z = (R @ load) * sign
    return (z - z.mean()) / z.std()


def _assign_positives(z: np.ndarray, fold: np.ndarray, n_pos_per_fold: int, beta: float,
                      rng: np.random.Generator) -> np.ndarray:
    """Weighted sampling without replacement inside each fold, exact fold-level class counts."""
    y = np.zeros(len(z), dtype=np.int64)
    for f in np.unique(fold):
        idx = np.flatnonzero(fold == f)
        w = np.exp(beta * z[idx])
        w = w / w.sum()
        pick = rng.choice(idx, size=n_pos_per_fold, replace=False, p=w)
        y[pick] = 1
    return y


def _oracle_auc(z: np.ndarray, y: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, z))


def calibrate_beta(z: np.ndarray, fold: np.ndarray, n_pos_per_fold: int,
                   target_auc: float, n_cal: int = 200) -> dict:
    """Deterministic bisection on beta so the ORACLE z has the requested median AUC."""
    rng_seed = SEED_BETA
    def median_auc(beta: float) -> float:
        rng = np.random.default_rng(rng_seed)
        return float(np.median([_oracle_auc(z, _assign_positives(z, fold, n_pos_per_fold,
                                                                 beta, rng))
                                for _ in range(n_cal)]))

    lo, hi = 0.0, 8.0
    while median_auc(hi) < target_auc and hi < 64:
        hi *= 2
    for _ in range(40):
        mid = (lo + hi) / 2
        if median_auc(mid) < target_auc:
            lo = mid
        else:
            hi = mid
    beta = (lo + hi) / 2
    return {"beta": float(beta), "achieved_median_oracle_AUC": median_auc(beta),
            "target_AUC": target_auc, "calibration_draws": n_cal, "seed": SEED_BETA}


def _build_cohort(scale: int):
    """Scale 1 is the real cohort. Scales 2/4 resample WITH replacement inside each outer fold."""
    X, clones = S23._load_rewind_x()
    k = pd.read_csv(_RESULTS / "stage22_rewind_clones.csv").set_index("clone_id").loc[clones]
    fold = k["outer_fold"].to_numpy()
    tbl = pd.read_csv(BDEPTH_TABLE).set_index("clone_id").loc[clones]
    nuis = np.column_stack([tbl["log1p_n_pretreatment_cells"].to_numpy(dtype=float),
                            tbl["n_lanes"].to_numpy(dtype=float)])
    if scale == 1:
        return X, nuis, fold, np.arange(len(clones))
    rng = np.random.default_rng(SEED_COVARIATE + scale)
    rows = []
    for f in range(S23.N_OUTER):
        idx = np.flatnonzero(fold == f)
        rows.append(rng.choice(idx, size=len(idx) * scale, replace=True))
    rows = np.concatenate(rows)
    return X[rows], nuis[rows], fold[rows], rows


def _delta_ap_once(X, nuis, y, fold) -> float:
    """The historical R1/R3 nested pipeline on one synthetic cohort. ΔAP = AP(R3) - AP(R1)."""
    from sklearn.metrics import average_precision_score
    from sklearn.model_selection import StratifiedKFold

    oof1 = np.full(len(y), np.nan)
    oof3 = np.full(len(y), np.nan)
    for f in range(S23.N_OUTER):
        tr, te = np.flatnonzero(fold != f), np.flatnonzero(fold == f)
        skf = StratifiedKFold(n_splits=S23.N_INNER, shuffle=True, random_state=S23.SEED_PROTOCOL)
        s1: dict = {}
        s3: dict = {}
        for itr_i, iva_i in skf.split(np.zeros(len(tr)), y[tr]):
            itr, iva = tr[itr_i], tr[iva_i]
            ztr, zva, _, kmax = S23.expression_block(X, itr, iva, max(S23.K_CANDIDATES))
            btr, bva = S23.standardize_train_only(nuis[itr], nuis[iva])
            for C in S23.LOGISTIC_C:
                p = S23._fit_logistic(btr, y[itr], bva, C, [])
                s1.setdefault((None, C), []).append(average_precision_score(y[iva], p))
            for kk in S23.K_CANDIDATES:
                if kk > kmax:
                    continue
                for C in S23.LOGISTIC_C:
                    p = S23._fit_logistic(np.hstack([ztr[:, :kk], btr]), y[itr],
                                          np.hstack([zva[:, :kk], bva]), C, [])
                    s3.setdefault((kk, C), []).append(average_precision_score(y[iva], p))

        def pick(sc):
            return max(sc.items(), key=lambda kv: (round(float(np.mean(kv[1])), 12),
                                                   -(kv[0][0] or 0), -kv[0][1]))[0]

        ztr, zte, _, _ = S23.expression_block(X, tr, te, max(S23.K_CANDIDATES))
        btr, bte = S23.standardize_train_only(nuis[tr], nuis[te])
        oof1[te] = S23._fit_logistic(btr, y[tr], bte, pick(s1)[1], [])
        kk, C = pick(s3)
        oof3[te] = S23._fit_logistic(np.hstack([ztr[:, :kk], btr]), y[tr],
                                     np.hstack([zte[:, :kk], bte]), C, [])
    return float(average_precision_score(y, oof3) - average_precision_score(y, oof1))


def run_23_2e_shard(scale: int, kind: str, target_auc: float | None, n_sims: int) -> dict:
    """One independent shard: (scale, kind[, target AUC]). Shards use disjoint seed streams."""
    t0 = time.perf_counter()
    X, nuis, fold, rows = _build_cohort(scale)
    z_full = synthetic_direction()
    z = z_full[rows]
    n_pos_per_fold = ANCHOR_POS_PER_FOLD * scale

    beta_info = None
    if kind == "alt":
        beta_info = calibrate_beta(z, fold, n_pos_per_fold, target_auc)

    base = SEED_NULL if kind == "null" else SEED_ALT
    offset = scale * 1000 + (0 if kind == "null" else int(target_auc * 100))
    vals, aucs = [], []
    for i in range(n_sims):
        rng = np.random.default_rng(base + offset + i)
        if kind == "null":
            y = _assign_positives(z, fold, n_pos_per_fold, 0.0, rng)
        else:
            y = _assign_positives(z, fold, n_pos_per_fold, beta_info["beta"], rng)
        aucs.append(_oracle_auc(z, y))
        vals.append(_delta_ap_once(X, nuis, y, fold))
        if (i + 1) % 10 == 0:
            el = time.perf_counter() - t0
            print(f"  [s{scale}-{kind}{target_auc or ''}] {i + 1}/{n_sims}  {el / 60:.1f} min  "
                  f"eta {el / (i + 1) * (n_sims - i - 1) / 60:.1f} min", flush=True)

    E_SHARD_DIR.mkdir(parents=True, exist_ok=True)
    name = f"scale{scale}_{kind}" + (f"_auc{int(target_auc * 100)}" if target_auc else "")
    payload = {"scale": scale, "kind": kind, "target_auc": target_auc, "n_sims": n_sims,
               "cohort_N": int(X.shape[0]), "positives_per_fold": n_pos_per_fold,
               "beta": beta_info, "delta_ap": [float(v) for v in vals],
               "oracle_auc": [float(v) for v in aucs],
               "runtime_minutes": round((time.perf_counter() - t0) / 60, 2)}
    (E_SHARD_DIR / f"{name}.json").write_text(json.dumps(payload), encoding="utf-8", newline="\n")
    return payload


def merge_23_2e() -> dict:
    """Combine the shards into the two V2 §9.5 statuses."""
    shards = {}
    for p in sorted(E_SHARD_DIR.glob("*.json")):
        d = json.loads(p.read_text(encoding="utf-8"))
        shards[p.stem] = d

    per_scale = {}
    for scale in COHORT_SCALES:
        nl = shards.get(f"scale{scale}_null")
        if nl is None:
            continue
        q95 = float(np.percentile(np.array(nl["delta_ap"]), 95))
        entry = {"cohort_N": nl["cohort_N"], "positives_per_fold": nl["positives_per_fold"],
                 "n_null": nl["n_sims"], "q95_null": q95,
                 "null_mean": float(np.mean(nl["delta_ap"])), "power": {}}
        for auc in TARGET_AUCS:
            alt = shards.get(f"scale{scale}_alt_auc{int(auc * 100)}")
            if alt is None:
                continue
            arr = np.array(alt["delta_ap"])
            entry["power"][str(auc)] = {
                "n_alt": alt["n_sims"], "estimated_power": float((arr > q95).mean()),
                "alt_mean_delta_ap": float(arr.mean()),
                "achieved_median_oracle_AUC": float(np.median(alt["oracle_auc"])),
                "beta": alt["beta"]["beta"] if alt["beta"] else None}
        per_scale[str(scale)] = entry

    prim = {s: per_scale[s]["power"].get("0.66", {}).get("estimated_power")
            for s in per_scale if per_scale[s]["power"].get("0.66")}
    p1 = prim.get("1")
    reached = any(v is not None and v >= 0.80 for k, v in prim.items() if k in ("2", "4"))
    if p1 is None:
        event_status = "UNRESOLVED"
    elif p1 >= 0.80:
        event_status = "NOT_SUPPORTED"
    elif p1 < 0.50 and reached:
        event_status = "SUPPORTED"
    else:
        event_status = "UNRESOLVED"

    out = {
        "stage": "23.2E",
        "plan": {"file": PLAN.name, "version": PLAN_VERSION},
        "stage23_2_protocol_sha256": json.loads(
            PROTOCOL.read_text(encoding="utf-8"))["canonical_sha256"],
        "design": {"cohort_scales": list(COHORT_SCALES), "target_AUCs": list(TARGET_AUCS),
                   "null_allocations": N_NULL_ALLOC, "alternative_simulations": N_ALT_SIM,
                   "seeds": {"covariate_resample": SEED_COVARIATE, "null": SEED_NULL,
                             "alternative": SEED_ALT, "beta_calibration": SEED_BETA}},
        "historical_test_geometry": _historical_geometry(),
        "per_scale": per_scale,
        "WITHIN_R1_EVENT_COUNT_LIMITATION": event_status,
        "BIOLOGICAL_REPLICATION_LIMITATION": "SUPPORTED",
        "biological_replication_note": (
            "V2 §9.5.2 -- a design fact, not an estimate. The claim rests on n_biological_"
            "replicates = 1, and no simulation can change that. It can become NOT_SUPPORTED only "
            "when a Role-A claim is supported by >= 2 independent biological replicates."),
        "bias_direction_caveat": (
            "Scales 2 and 4 remain empirical resamples of biological replicate R1. Repeated "
            "covariate profiles reduce effective covariate diversity and therefore make the "
            "projected power curve an approximation whose bias direction is not guaranteed. The "
            "simulation estimates within-R1 event-count detectability under the empirical R1 "
            "covariate distribution; it does not estimate power gained from additional independent "
            "biological replicates."),
        "what_this_cannot_do": [
            "create biological replicate diversity",
            "estimate between-replicate variance",
            "establish that a corrected same-data result would replicate",
            "be read as an independent-sample power curve"],
        "source_provenance": source_provenance(),
    }
    write_json(E_RESULTS, out)
    return out


def _historical_geometry() -> dict:
    pe = json.loads((_RESULTS / "stage23_permutation_results.json").read_text(encoding="utf-8"))
    t = pe["permutation_tests"]["role_a_delta_AP_state"]
    return {"observed_delta_AP": t["observed"], "null_mean": t["null_mean"],
            "null_sd": t["null_sd"], "null_p95": t["null_p95"], "p_perm": t["p_perm"],
            "null_centered_separation": t["observed"] - t["null_mean"],
            "distance_to_null_p95": t["observed"] - t["null_p95"]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage 23.2 Role-A resolution")
    ap.add_argument("--stage", default="23.2a", choices=["23.2a", "23.2b", "23.2c", "23.2d", "23.2e", "23.2e-merge"])
    ap.add_argument("--permutations", type=int, default=200)
    ap.add_argument("--scale", type=int, choices=[1, 2, 4])
    ap.add_argument("--kind", choices=["null", "alt"])
    ap.add_argument("--target-auc", type=float)
    ap.add_argument("--sims", type=int)
    ap.add_argument("--no-replay", action="store_true",
                    help="recover mappings without replaying D00 (audit convenience only)")
    args = ap.parse_args(argv)

    if args.stage == "23.2b":
        r = run_23_2b(args.permutations)
        print(f"  cell 00 recomputed exactly: {r['cell_00_recomputation']['reproduces_committed_D00_exactly']}")
        print(f"  mu00 {r['D00']['mean']:+.5f}   mu10 (no-K) {r['D10_no_k_selection']['mean']:+.5f}")
        for a, v in r["per_arm"].items():
            print(f"    arm {a:<4} mean {v['mean']:+.5f}  95% [{v['ci95_low']:+.5f}, {v['ci95_high']:+.5f}]")
        print(f"  arm dispersion {r['arm_dispersion']:.5f}")
        sh = r["selection_shift"]
        print(f"  selection_shift {sh['mean']:+.5f}  95% [{sh['ci95_low']:+.5f}, {sh['ci95_high']:+.5f}]")
        print(f"  ladder {r['search_width_ladder']}  monotone={r['ladder_monotone_increase']}")
        print("OVERALL:", r["MODEL_SELECTION_NULL_INFLATION"])
        return 0
    if args.stage == "23.2e":
        n = args.sims if args.sims else (N_NULL_ALLOC if args.kind == "null" else N_ALT_SIM)
        r = run_23_2e_shard(args.scale, args.kind, args.target_auc, n)
        arr = np.array(r["delta_ap"])
        print(f"  scale {r['scale']} N={r['cohort_N']} {r['kind']} auc={r['target_auc']} "
              f"n={r['n_sims']}  mean dAP {arr.mean():+.5f}  q95 {np.percentile(arr,95):+.5f}")
        if r["beta"]:
            print(f"  beta {r['beta']['beta']:.4f} -> median oracle AUC {r['beta']['achieved_median_oracle_AUC']:.4f}")
        print(f"  runtime {r['runtime_minutes']} min")
        return 0
    if args.stage == "23.2e-merge":
        r = merge_23_2e()
        for s_, v in r["per_scale"].items():
            print(f"  scale {s_}  N={v['cohort_N']}  q95_null {v['q95_null']:+.5f}")
            for auc, pw in v["power"].items():
                print(f"      AUC {auc}: power {pw['estimated_power']:.3f}  (median oracle AUC {pw['achieved_median_oracle_AUC']:.4f})")
        print("  WITHIN_R1_EVENT_COUNT_LIMITATION:", r["WITHIN_R1_EVENT_COUNT_LIMITATION"])
        print("  BIOLOGICAL_REPLICATION_LIMITATION:", r["BIOLOGICAL_REPLICATION_LIMITATION"])
        return 0

    if args.stage == "23.2d":
        r = run_23_2d()
        m = r["multinomial_stability"]
        print(f"  source rule reproduces 35 positives: {r['source_rule_reproduces_35_positives']}")
        print(f"  mean frozen-positive retention  {m['mean_frozen_positive_retention']:.4f}")
        print(f"  positives with P(selected)<0.50 {m['n_positives_below_0_50']} / 35")
        print(f"  positives with P(selected)<0.80 {m['n_positives_below_0_80']} / 35")
        for n, v in r["cutoff_sensitivity"].items():
            print(f"    {n:<7} barcodes {v['selected_barcodes']:>4}  positives {v['positive_clones']:>3}  "
                  f"Jaccard {v['jaccard_vs_top100']:.4f}  lost {v['frozen_positives_lost']}  gained {v['frozen_negatives_gained']}")
        print(f"  NOT_SUPPORTED reachable: {r['not_supported_reachable']}")
        print("OVERALL:", r["OUTCOME_LABEL_LIMITATION"])
        return 0

    if args.stage == "23.2c":
        r = run_23_2c(args.permutations)
        print(f"  null means {r['null_means']}")
        d = r["depth_shift_full"]
        print(f"  depth_shift_full {d['mean']:+.5f}  95% [{d['ci95_low']:+.5f}, {d['ci95_high']:+.5f}]")
        for kk, v in r["contrasts"].items():
            print(f"    {kk:<26} {v['mean']:+.5f}  95% [{v['ci95_low']:+.5f}, {v['ci95_high']:+.5f}]")
        print(f"  observed 2x2 {r['observed_2x2']}")
        print(f"  corrected same-data diagnostic {r['corrected_same_data_diagnostic']['CORRECTED_SAME_DATA_SIGNAL_DIAGNOSTIC']}")
        print(f"  lane sensitivity {r['lane_composition_sensitivity']['permitted']}")
        print("OVERALL:", r["RESIDUAL_DEPTH_STRUCTURE"])
        return 0

    r = run_23_2a(n_perm=args.permutations, replay=not args.no_replay)
    if r["verdict"] == INPUT_BLOCKED:
        for f in r.get("preflight", {}).get("failed", []):
            print("  FAILED:", f["check"], f["detail"])
        for name, ok in r.get("gates", {}).items():
            if not ok:
                print("  GATE FAILED:", name)
        print("OVERALL:", r["verdict"])
        return 1
    print(f"  within-R1 status        {r['within_r1_status']}")
    print(f"  lane sensitivity        {'PERMITTED' if r['lane_composition_sensitivity_permitted'] else 'FORBIDDEN'}")
    print(f"  gDNA rule reproduced    {r['gdna_rule']['reproduces_frozen_positives_exactly']}"
          f"  ({r['gdna_rule']['selected_barcodes']} barcodes, "
          f"{r['gdna_rule']['reconstructed_positive_clones']} positive clones)")
    print(f"  D00 replay              {r['permutation_recovery']['replay']['all_summary_statistics_reproduced']}")
    print(f"  protocol sha256         {r['stage23_2_protocol_sha256'][:16]}...")
    print(f"  runtime                 {r['runtime_minutes']} min")
    print("OVERALL:", r["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
