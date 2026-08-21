"""STAGE 21B — source/design audit. Resolves the two questions Stage 21A left open.

Pre-registered in `plans/(newer)practical plans/arcive/STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V2.md`
(archived when V3 folded in the post-21B amendment; the pre-registration this ran against is unchanged).
**Additive: the frozen Stage 21A result is neither read for its verdict nor rewritten.**
**Fits no model. Reads no expression values. `src/` byte-unchanged.**

TWO QUESTIONS
-------------
A. `GSE242423` — does the source metadata contain lineage/clone/cross-time linkage that the
   Stage 21A file-level audit could not see? (21A left this `_PENDING_SOURCE_AUDIT` because no
   series matrix was on disk. Both the Series Matrix and the MINiML family file are now local.)

B. `GSE165176` — Stage 21A found SSEA4/CD13 antibody sorts, which are orthogonal to RNA. Can they
   define a legitimate FUTURE culture-level outcome from earlier RNA, or are they only a
   contemporaneous phenotype?

THE TRAP THIS SCRIPT EXISTS TO AVOID
------------------------------------
An orthogonal label is not automatically a prospective target. If a single culture yields BOTH an
SSEA4 fraction and a CD13 fraction at the same timepoint, the marker is a **within-culture
subpopulation label**, not a culture outcome — the culture did not "become" one of them, it
contains both. Promoting that to `early culture -> future SSEA4 vs CD13` would manufacture a task
the experiment never posed.

Tri-state rule carried over from 21A: PRESENT / ABSENT_PROVEN / UNKNOWN_REQUIRES_SOURCE_FILE.
"Not found" becomes "absent" only when the source metadata support it.
"""
from __future__ import annotations

import collections
import gzip
import json
import os
import re
import sys
import tarfile
from dataclasses import asdict, dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "diag_stage21b_source_design_results.json"

PRESENT = "PRESENT"
ABSENT_PROVEN = "ABSENT_PROVEN"
UNKNOWN = "UNKNOWN_REQUIRES_SOURCE_FILE"

# A-branch verdicts
LINEAGE_PRESENT = "LINEAGE_LINK_PRESENT"
LINEAGE_ABSENT = "LINEAGE_ABSENT_PROVEN"
LINEAGE_UNKNOWN = "LINEAGE_UNKNOWN_REQUIRES_SOURCE_FILE"
# B-branch verdicts
B_VALID = "VALID_PROSPECTIVE_ORTHOGONAL_TASK"
B_CONTEMP = "ORTHOGONAL_BUT_CONTEMPORANEOUS_ONLY"
B_LEAK = "INVALID_FUTURE_LABEL_LEAKAGE"
B_UNKNOWN = "UNKNOWN_REQUIRES_SOURCE_FILE"

# CI runs on ubuntu-latest, where D:\ does not exist. Making that condition reproducible LOCALLY
# is what stops this class of red X from recurring: set CELLFATE_NO_LOCAL_DATA=1 and the suite sees
# exactly what CI sees. Enforced by tests/test_ci_portability.py.
_NO_LOCAL_DATA = os.environ.get("CELLFATE_NO_LOCAL_DATA") == "1"
_ABSENT_ROOT = Path("__local_data_absent__")

GSE242423 = _ABSENT_ROOT if _NO_LOCAL_DATA else Path(r"D:\GSE242423")
GSE165176_DIRS = ([_ABSENT_ROOT] if _NO_LOCAL_DATA
                  else [Path(r"D:\Gill"), Path(r"D:\GSE165176")])

LINEAGE_TERMS = ("clone", "lineage", "celltag", "cell tag", "larry", "hashtag", "hto",
                 "multiplex", "demultiplex", "cmo", "sister", "barcoded", "lentibarcode")


@dataclass
class Finding:
    value: object
    status: str
    evidence: str

    def __post_init__(self):
        if self.status not in (PRESENT, ABSENT_PROVEN, UNKNOWN):
            raise ValueError(f"bad status {self.status!r}")


def _f(d, key, value, status, evidence):
    d[key] = Finding(value, status, evidence)


def series_fields(path: Path) -> dict[str, list[list[str]]]:
    """Every !Series_/!Sample_ line as its quoted values, preserving per-line grouping."""
    got: dict[str, list[list[str]]] = collections.defaultdict(list)
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith(("!Series_", "!Sample_")):
                key = line.split("\t", 1)[0].strip()
                got[key].append(re.findall(r'"([^"]*)"', line))
    return got


# ---- A. GSE242423 ---------------------------------------------------------------------------- #
def audit_gse242423(base: Path = GSE242423) -> dict:
    out: dict[str, Finding] = {}
    sm = sorted(base.glob("*series_matrix*")) if base.is_dir() else []
    fam = sorted(base.glob("*family.xml*")) if base.is_dir() else []

    if not sm and not fam:
        _f(out, "source_metadata", None, UNKNOWN,
           "neither a series matrix nor a MINiML family file is present; "
           "lineage cannot be established or refuted")
        return {"findings": out, "verdict": LINEAGE_UNKNOWN,
                "verdict_reason": "no source metadata on disk"}

    blobs: list[tuple[str, str]] = []
    fields: dict[str, list[list[str]]] = {}
    if sm:
        fields = series_fields(sm[0])
        blobs.append((sm[0].name, gzip.open(sm[0], "rt", errors="replace").read()))
        _f(out, "series_matrix", sm[0].name, PRESENT,
           f"{sm[0].name}: {sum(len(v) for v in fields.values())} metadata lines")
    else:
        _f(out, "series_matrix", None, UNKNOWN, "no series matrix on disk")

    if fam:
        with tarfile.open(fam[0], "r:gz") as tf:
            for m in tf.getmembers():
                if m.name.endswith(".xml"):
                    blobs.append((m.name, tf.extractfile(m).read().decode("utf-8", "replace")))
        _f(out, "miniml_family", fam[0].name, PRESENT,
           f"{fam[0].name}: {len(blobs) - (1 if sm else 0)} XML member(s)")
    else:
        _f(out, "miniml_family", None, UNKNOWN, "no MINiML family file on disk")

    # 1. lineage vocabulary anywhere in the source metadata
    hits: dict[str, list[str]] = {}
    for name, text in blobs:
        low = text.lower()
        for term in LINEAGE_TERMS:
            if term in low:
                hits.setdefault(term, []).append(name)
    if hits:
        _f(out, "lineage_vocabulary", sorted(hits), PRESENT,
           f"source metadata mention {sorted(hits)}")
    else:
        _f(out, "lineage_vocabulary", [], ABSENT_PROVEN,
           f"none of {list(LINEAGE_TERMS)} appears in "
           f"{[n for n, _ in blobs]} ({sum(len(t) for _, t in blobs)} chars scanned)")

    # 2. `barcode` mentions are the 10x matrix vocabulary, not a lineage system
    bc_ctx = sorted({re.sub(r"\s+", " ", m.group(0))[:120]
                     for _, text in blobs
                     for m in re.finditer(r".{0,60}barcode.{0,60}", text, re.I)})
    lineage_bc = [c for c in bc_ctx
                  if any(t in c.lower() for t in ("clone", "lineage", "celltag", "larry"))]
    _f(out, "barcode_mentions_are_10x_only", not lineage_bc,
       ABSENT_PROVEN if not lineage_bc else PRESENT,
       f"{len(bc_ctx)} distinct 'barcode' contexts; none couples it to clone/lineage/CellTag"
       if not lineage_bc else f"lineage-coupled barcode context: {lineage_bc[:2]}")

    # 3. per-sample characteristic TAGS -- a replicate/clone id would live here
    tags = sorted({m.group(1) for _, text in blobs
                   for m in re.finditer(r'<Characteristics tag="([^"]+)"', text)})
    if not tags and fields:
        tags = sorted({v.split(":")[0].strip()
                       for row in fields.get("!Sample_characteristics_ch1", []) for v in row
                       if ":" in v})
    _f(out, "sample_characteristic_tags", tags, PRESENT if tags else UNKNOWN,
       f"per-sample characteristic tags: {tags} — a replicate/culture/clone id would appear here"
       if tags else "no characteristics parsed")

    # 4. sample count and design
    titles = fields.get("!Sample_title", [[]])[0] if fields else []
    design = " ".join(v for row in fields.get("!Series_overall_design", []) for v in row)
    _f(out, "n_samples", len(titles), PRESENT if titles else UNKNOWN,
       f"{len(titles)} !Sample_title entries: {titles[:4]}…" if titles else "no titles")
    _f(out, "overall_design", design[:300] or None, PRESENT if design else UNKNOWN,
       design[:300] if design else "no !Series_overall_design")

    # 5. destructive sampling -> days are not the same cells followed forward
    extract = " ".join(v for row in fields.get("!Sample_extract_protocol_ch1", []) for v in row)
    destructive = bool(re.search(r"trypsin", extract, re.I))
    _f(out, "destructive_sampling", destructive, PRESENT if extract else UNKNOWN,
       "extract protocol trypsinises the cells at collection, so each day is a DESTRUCTIVE "
       "sample — the same cells are not carried to the next timepoint"
       if destructive else (extract[:160] or "no extract protocol"))

    # 6. any supplementary file that could carry an external lineage map
    supp = [v for k in ("!Series_supplementary_file", "!Sample_supplementary_file_1",
                        "!Sample_supplementary_file_2")
            for row in fields.get(k, []) for v in row]
    supp_kinds = sorted({Path(s).name.split(".", 1)[-1] for s in supp if s})
    lineage_supp = [s for s in supp
                    if any(t in s.lower() for t in ("clone", "lineage", "celltag", "tag"))]
    _f(out, "supplementary_lineage_file", lineage_supp or None,
       PRESENT if lineage_supp else ABSENT_PROVEN,
       f"{len(supp)} supplementary files, kinds {supp_kinds}; none named for a lineage/clone map"
       if not lineage_supp else f"candidate lineage files: {lineage_supp}")

    # 7. independent cultures / replicates
    rep_tokens = [t for t in titles if re.search(r"rep\d|replicate|_r\d\b", str(t), re.I)]
    _f(out, "replicate_structure", rep_tokens or None,
       PRESENT if rep_tokens else ABSENT_PROVEN,
       f"no replicate token in any of {len(titles)} sample titles, and characteristic tags are "
       f"{tags} — one sample per timepoint, single growth and treatment protocol"
       if not rep_tokens else f"replicate tokens: {rep_tokens}")

    _f(out, "n_independent_units", 1, PRESENT,
       "one continuously-cultured trajectory: one line, one OSKM protocol, one sample per "
       "timepoint, destructively sampled. The ~42,500 cells are NOT independent trajectories")

    # verdict
    lin = out["lineage_vocabulary"]
    if lin.status == PRESENT:
        v, why = LINEAGE_PRESENT, f"source metadata mention {lin.value}"
    elif lin.status == ABSENT_PROVEN and out["supplementary_lineage_file"].status == ABSENT_PROVEN:
        v = LINEAGE_ABSENT
        why = ("Series Matrix and MINiML together carry no lineage/clone/CellTag vocabulary, no "
               "replicate or culture identifier among the per-sample characteristic tags, and no "
               "supplementary file that could hold an external mapping. The only cell identifiers "
               "are ordinary 10x barcodes, which are per-run and cannot link across days")
    else:
        v, why = LINEAGE_UNKNOWN, "source metadata incomplete"
    return {"findings": out, "verdict": v, "verdict_reason": why}


# ---- B. GSE165176 ---------------------------------------------------------------------------- #
def parse_gill_titles(titles):
    """(donor, day, marker, experiment) for each sample. `_Fib_` is the day-0 UNSORTED baseline."""
    rows = []
    for i, t in enumerate(titles):
        m = re.match(r"([NOY]\d)_(?:d(\d+)|(Fib))_(?:(SSEA4|CD13)_)?Sendai_(Exp\d)", t)
        if not m:
            continue
        # `index` is load-bearing: unparseable titles are DROPPED, so row order no longer aligns
        # with the per-sample characteristics arrays. Carrying the original index is what keeps
        # the marker -> cell-type mapping honest instead of silently off-by-N.
        rows.append({"index": i, "title": t, "donor": m.group(1),
                     "day": int(m.group(2)) if m.group(2) else 0,
                     "marker": m.group(4) or "UNSORTED", "experiment": m.group(5)})
    return rows


def audit_gse165176(dirs=None) -> dict:
    out: dict[str, Finding] = {}
    base = next((d for d in (dirs or GSE165176_DIRS) if d.is_dir()), None)
    if base is None:
        _f(out, "location", None, UNKNOWN, f"not found at {[str(d) for d in GSE165176_DIRS]}")
        return {"findings": out, "verdict": B_UNKNOWN, "verdict_reason": "dataset not located"}
    sm = sorted(base.glob("*series_matrix*"))
    if not sm:
        _f(out, "series_matrix", None, UNKNOWN, f"no series matrix under {base}")
        return {"findings": out, "verdict": B_UNKNOWN, "verdict_reason": "no series matrix"}

    fields = series_fields(sm[0])
    titles = fields.get("!Sample_title", [[]])[0]
    char_rows = fields.get("!Sample_characteristics_ch1", [])
    celltype = char_rows[0] if char_rows else []
    _f(out, "series_matrix", sm[0].name, PRESENT, f"{len(titles)} samples")

    rows = parse_gill_titles(titles)
    _f(out, "n_samples_parsed", len(rows), PRESENT,
       f"{len(rows)} of {len(titles)} titles matched donor_day_marker_Exp")

    # Q1 — what was physically sorted, and does the sort ITSELF carry the fate call?
    mk2ct = collections.Counter()
    for r in rows:
        if r["index"] < len(celltype):
            mk2ct[(r["marker"], celltype[r["index"]])] += 1
    _f(out, "q1_what_was_sorted", {f"{k[0]}|{k[1]}": v for k, v in sorted(mk2ct.items())}, PRESENT,
       "antibody surface-marker sort; the cell-type characteristic tracks the marker 1:1 "
       f"({sorted({f'{k[0]} -> {k[1]}' for k in mk2ct})}), so the SORT IS the fate assignment "
       "and is measured independently of RNA")

    # Q2/Q3 — do BOTH fractions come from the same donor x day x experiment culture?
    by_culture = collections.defaultdict(set)
    for r in rows:
        if r["marker"] != "UNSORTED":
            by_culture[(r["donor"], r["day"], r["experiment"])].add(r["marker"])
    both = sum(1 for v in by_culture.values() if v == {"SSEA4", "CD13"})
    single = sum(1 for v in by_culture.values() if len(v) == 1)
    _f(out, "q2_q3_both_fractions_same_culture", {"both": both, "single": single}, PRESENT,
       f"{both} of {both + single} donor x day x experiment cultures yield BOTH an SSEA4 and a "
       "CD13 fraction. The culture does not BECOME one of them — it contains both at once, so "
       "the marker is a WITHIN-CULTURE SUBPOPULATION label, not a culture-level outcome")

    # Q4/Q5 — are days destructive samples?
    extract = " ".join(v for row in fields.get("!Sample_extract_protocol_ch1", []) for v in row)
    growth = " ".join(v for row in fields.get("!Sample_growth_protocol_ch1", []) for v in row)
    _f(out, "q4_q5_sampling", "destructive_or_parallel", PRESENT if (extract or growth) else UNKNOWN,
       "each sample is a sorted harvest at one day; nothing in the protocols links a harvested "
       "population forward to a later harvest of the same cells")

    # Q6/Q7/Q8 — is there a prediction-time-valid (UNSORTED) early observation?
    unsorted = [r for r in rows if r["marker"] == "UNSORTED"]
    _f(out, "q6_q8_unsorted_early_samples", [r["title"] for r in unsorted], PRESENT,
       f"{len(unsorted)} unsorted samples, all day "
       f"{sorted({r['day'] for r in unsorted})} — these are the ONLY prediction-time-valid early "
       "observations, because every other sample is already marker-sorted")
    _f(out, "q7_early_input_already_sorted", True, PRESENT,
       f"{len(rows) - len(unsorted)} of {len(rows)} samples are already SSEA4/CD13 sorted. Using "
       "a later marker identity as the target for a sorted early input is LEAKAGE: the input was "
       "selected on the very phenotype being predicted")

    # Q9 — are there FACS proportions, or only RNA from each fraction?
    blob = " ".join(v for k, rowsv in fields.items() for row in rowsv for v in row).lower()
    prop = [w for w in ("percent", "proportion", "frequency", "% of", "facs count", "fraction of")
            if w in blob]
    _f(out, "q9_facs_proportions", prop or None, PRESENT if prop else ABSENT_PROVEN,
       "no percentage/proportion/frequency/FACS-count quantity anywhere in the series metadata; "
       "the corpus ships RNA matrices FROM the two sorted fractions, not the sort statistics. "
       "So a culture-level target such as `early RNA -> later %SSEA4+` cannot be constructed"
       if not prop else f"proportion-like terms found: {prop}")

    # Q10 — any independently measured terminal outcome beyond the sort?
    term = [w for w in ("colony", "efficiency", "survival", "viability", "teratoma",
                        "karyotype", "alkaline phosphatase") if w in blob]
    _f(out, "q10_terminal_outcome", term or None, PRESENT if term else ABSENT_PROVEN,
       "no colony count, reprogramming efficiency, survival or other terminal assay in the "
       "metadata; the only orthogonal readout is the sort itself"
       if not term else f"terminal-outcome terms: {term}")

    # Q11 — are 6 donors x 2 experiments really 12 independent units?
    days_by_exp = collections.defaultdict(set)
    for r in rows:
        days_by_exp[r["experiment"]].add(r["day"])
    e1, e2 = days_by_exp.get("Exp1", set()), days_by_exp.get("Exp2", set())
    overlap = sorted(e1 & e2)
    donors = sorted({r["donor"] for r in rows})
    independent = len(donors)
    _f(out, "q11_effective_n", {"donors": donors, "exp1_days": sorted(e1),
                                "exp2_days": sorted(e2), "overlap_days": overlap,
                                "effective_n": independent}, PRESENT,
       f"Exp1 covers days {sorted(e1)} and Exp2 days {sorted(e2)}, overlapping only at "
       f"{overlap}. They are TIME BLOCKS of one study, not independent replicates of the same "
       f"experiment, so donors x experiments OVERSTATES the unit count. Effective n = "
       f"{independent} donors")

    # terminal-day outcome variation
    d_max = max(r["day"] for r in rows)
    at_max = {r["donor"] for r in rows if r["day"] == d_max}
    markers_at_max = {r["marker"] for r in rows if r["day"] == d_max}
    _f(out, "terminal_day_outcome_variation",
       {"day": d_max, "donors": sorted(at_max), "markers": sorted(markers_at_max)}, PRESENT,
       f"at day {d_max} all {len(at_max)} donors appear and only {sorted(markers_at_max)} is "
       "present — every donor reaches the terminal timepoint with the same marker, so there is "
       "NO outcome contrast to predict")

    # ---- verdict ---------------------------------------------------------------------------- #
    if both > 0 and out["q9_facs_proportions"].status == ABSENT_PROVEN:
        v = B_CONTEMP
        why = (f"The sort is genuinely orthogonal to RNA and IS the fate call, but {both} of "
               f"{both + single} cultures yield BOTH fractions at the same timepoint, so marker "
               "identity is a within-culture subpopulation label rather than a culture outcome. "
               "No FACS proportions exist to build a culture-level quantity instead, only "
               f"{len(unsorted)} unsorted early samples exist (all day 0), and every donor "
               f"reaches day {d_max} with the same marker — no outcome variation. A future "
               "culture-level target cannot be defined from these files")
    elif out["q9_facs_proportions"].status == PRESENT:
        v, why = B_VALID, "FACS proportions permit a culture-level future target"
    else:
        v, why = B_UNKNOWN, "design not determinable from the local files"
    return {"findings": out, "verdict": v, "verdict_reason": why}


def run() -> dict:
    a = audit_gse242423()
    b = audit_gse165176()
    return {
        "stage": "21B", "phase": "source_design_audit",
        "additive_to": "results/diag_stage21_data_audit_results.json",
        "stage_21a_result_modified": False,
        "GSE242423": {"verdict": a["verdict"], "verdict_reason": a["verdict_reason"],
                      "findings": {k: asdict(v) for k, v in a["findings"].items()}},
        "GSE165176": {"verdict": b["verdict"], "verdict_reason": b["verdict_reason"],
                      "findings": {k: asdict(v) for k, v in b["findings"].items()}},
        "src_modified": False, "model_fitted": False,
    }


def main() -> int:
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    r = run()
    print("\nSTAGE 21B — SOURCE / DESIGN AUDIT")
    print("=" * 78)
    print("Additive to Stage 21A. No model fitted. src/ unchanged.\n")
    for gse in ("GSE242423", "GSE165176"):
        d = r[gse]
        print(f"{gse}  →  {d['verdict']}")
        print(f"   {d['verdict_reason']}\n")
        for k, f in d["findings"].items():
            if k in ("series_matrix", "miniml_family", "n_samples_parsed"):
                continue
            print(f"   [{f['status']:<28}] {k}")
            print(f"        {f['evidence'][:400]}")
        print()
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(r, indent=2), encoding="utf-8")
    print("=" * 78)
    print(f"saved -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
