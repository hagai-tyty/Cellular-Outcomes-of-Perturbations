"""STAGE 21A — local dataset geometry audit for a PROSPECTIVE fate task.

Pre-registered in `plans/(newer)practical plans/arcive/STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V2.md`
(archived when V3 folded in the post-21B amendment; the pre-registration this ran against is unchanged).
**Fits no model. Reads no expression values. `src/` is not imported for anything but paths.**

WHAT THIS DECIDES
-----------------
Whether any dataset already on disk can pose

    X_before + U  ->  Y_future     (Y obtained INDEPENDENTLY of the input RNA)

and, if not, which of the four task levels it can support:

    LEVEL_3 STRICT_LINEAGE     early RNA + clone/lineage link + later independent outcome
    LEVEL_2 CULTURE_FORWARD    early culture RNA + later independent culture outcome
    LEVEL_1 TRAJECTORY_FORWARD early population RNA -> later population distribution (pilot only)
    LEVEL_0 INVALID            no valid early->late mapping

TWO RULES THIS SCRIPT EXISTS TO ENFORCE
---------------------------------------
**1. Every classification carries its evidence.** A verdict of `TRAJECTORY_FORWARD` is useless if
nobody can see why. Each field is a `Finding` recording the value, HOW it was established, and
WHICH file established it. `_evidence` travels into the JSON alongside every verdict.

**2. "Not found" is NOT "proven absent".** If the file that would answer a question was never
downloaded, the answer is `UNKNOWN_REQUIRES_SOURCE_AUDIT` -- never `ABSENT`. A dataset must not be
discarded because a supplementary mapping file is missing from a local mirror. Concretely:
`GSE242423` has no `series_matrix` on disk, so its replicate/arm/lineage structure is UNKNOWN even
though its barcode ENCODING is provably plain 10x. Those are different claims and the audit keeps
them apart.

Consequence for the verdict: a level is only ruled OUT when the evidence proves it out. If a
required field is UNKNOWN the verdict is suffixed `_PENDING_SOURCE_AUDIT`, which is a request for
one more download -- not a rejection.
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
OUT = _RESULTS / "diag_stage21_data_audit_results.json"

# --- status values ------------------------------------------------------------------------- #
PRESENT = "PRESENT"
ABSENT_PROVEN = "ABSENT_PROVEN"
UNKNOWN = "UNKNOWN_REQUIRES_SOURCE_AUDIT"

# --- task levels --------------------------------------------------------------------------- #
LEVEL_3 = "STRICT_LINEAGE"
LEVEL_2 = "CULTURE_FORWARD"
LEVEL_1 = "TRAJECTORY_FORWARD"
LEVEL_0 = "INVALID_PROSPECTIVE"
PENDING = "_PENDING_SOURCE_AUDIT"

# Where each dataset actually lives. GSE165176 is mirrored under `D:\Gill`, NOT `D:\GSE165176`;
# looking only at the canonical path would report a present dataset as missing -- the exact
# "not found is not absent" error this script exists to avoid.
DEFAULT_LOCATIONS = {
    "GSE242423": [r"D:\GSE242423"],
    "GSE165176": [r"D:\Gill", r"D:\GSE165176"],
    "GSE165177": [r"D:\GSE165177"],
    "GSE165178": [r"D:\GSE165178"],
    "GSE165179": [r"D:\GSE165179"],
    "GSE113957": [r"D:\GSE113957"],
    "GSE297234": [r"D:\GSE297234"],
}

# A predictor carrying any of these is outcome-derived and must hard-fail (plan §19.4).
FORBIDDEN_PREDICTORS = ("final_fate", "response_status", "post_treatment_cluster",
                        "survivor_label", "outcome", "reprogrammed")


@dataclass
class Finding:
    """A single audited fact: what, how sure, and on what evidence."""
    value: object
    status: str
    evidence: str

    def __post_init__(self):
        if self.status not in (PRESENT, ABSENT_PROVEN, UNKNOWN):
            raise ValueError(f"bad status {self.status!r}")


@dataclass
class DatasetAudit:
    dataset: str
    path: str | None
    present: bool
    findings: dict = field(default_factory=dict)
    level: str = LEVEL_0
    level_reason: str = ""
    ruled_out: dict = field(default_factory=dict)

    def add(self, key: str, value, status: str, evidence: str) -> None:
        self.findings[key] = Finding(value, status, evidence)

    def get(self, key: str) -> Finding | None:
        return self.findings.get(key)


# ---- low-level readers ---------------------------------------------------------------------- #
def locate(dataset: str, candidates=None) -> Path | None:
    for c in (candidates or DEFAULT_LOCATIONS.get(dataset, [])):
        p = Path(c)
        if p.is_dir():
            return p
    return None


def head_gz(path: Path, n: int = 3) -> list[str]:
    """First `n` lines of a gzipped text file, without reading the whole thing."""
    out = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i >= n:
                break
            out.append(line.rstrip("\n"))
    return out


def series_matrix_fields(path: Path, keys=("!Sample_title", "!Sample_characteristics_ch1",
                                           "!Sample_geo_accession", "!Sample_source_name_ch1")):
    """Pull the sample-level metadata lines out of a GEO series matrix."""
    found: dict[str, list[str]] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            for k in keys:
                if line.startswith(k):
                    vals = re.findall(r'"([^"]*)"', line)
                    found.setdefault(k, []).extend(vals)
    return found


def barcode_encoding(path: Path) -> tuple[str, str]:
    """Classify a barcodes file's ENCODING. This is provable from the file itself.

    A plain 10x barcode (`AAACCCAAGAAACACT-1`) carries no lineage information. A CellTag /
    LARRY-style file would carry an extra field or a non-10x token. Proving the encoding is NOT
    the same as proving the study had no lineage tracing -- that lives in the series matrix.
    """
    lines = head_gz(path, 3)
    if not lines:
        return "EMPTY", "file has no rows"
    first = lines[0]
    if re.fullmatch(r"[ACGT]{14,20}(-\d+)?", first):
        return "PLAIN_10X", f"first barcode {first!r} matches ^[ACGT]{{14,20}}(-N)?$"
    if "\t" in first or "," in first:
        return "MULTI_FIELD", f"first row {first!r} has a delimiter -- may carry extra annotation"
    return "OTHER", f"first row {first!r} is not a plain 10x barcode"


def count_gz_lines(path: Path, cap: int = 5_000_000) -> int:
    n = 0
    with gzip.open(path, "rt", errors="replace") as fh:
        for _ in fh:
            n += 1
            if n >= cap:
                break
    return n


# ---- the classifier ------------------------------------------------------------------------- #
def classify(a: DatasetAudit) -> DatasetAudit:
    """Assign the strongest DEFENSIBLE level, and record what was ruled out and how.

    A level is ruled OUT only when the evidence proves it out. If a field needed to rule it out is
    UNKNOWN, the level stays undetermined and the verdict is suffixed `_PENDING_SOURCE_AUDIT`.
    """
    if not a.present:
        a.level = LEVEL_0
        a.level_reason = "dataset not found at any candidate path"
        return a

    def st(key):
        f = a.get(key)
        return f.status if f else UNKNOWN

    def val(key):
        f = a.get(key)
        return f.value if f else None

    pending = []

    # --- LEVEL 3 needs a clone/lineage link AND an independent outcome ---------------------- #
    if st("lineage_link") == PRESENT and st("independent_outcome") == PRESENT:
        a.level, a.level_reason = LEVEL_3, "lineage link and independent outcome both established"
        return a
    if st("lineage_link") == UNKNOWN:
        pending.append("lineage_link")
        a.ruled_out[LEVEL_3] = "NOT RULED OUT — lineage structure is unknown locally"
    else:
        a.ruled_out[LEVEL_3] = f"ruled out: lineage_link={val('lineage_link')} " \
                               f"({a.get('lineage_link').evidence if a.get('lineage_link') else ''})"

    # --- LEVEL 2 needs >1 independent culture/replicate AND an independent outcome ---------- #
    n_units = val("n_independent_units")
    if st("independent_outcome") == PRESENT and isinstance(n_units, int) and n_units >= 2:
        a.level, a.level_reason = LEVEL_2, f"independent outcome with {n_units} independent units"
        return a
    if st("independent_outcome") == UNKNOWN or st("n_independent_units") == UNKNOWN:
        pending.append("independent_outcome/n_independent_units")
        a.ruled_out[LEVEL_2] = "NOT RULED OUT — outcome independence or unit count unknown locally"
    else:
        a.ruled_out[LEVEL_2] = (f"ruled out: independent_outcome={val('independent_outcome')}, "
                                f"n_independent_units={n_units}")

    # --- LEVEL 1 needs only a timecourse with >=2 distinct timepoints ------------------------ #
    n_tp = val("n_timepoints")
    if isinstance(n_tp, int) and n_tp >= 2:
        a.level = LEVEL_1 + (PENDING if pending else "")
        a.level_reason = (f"{n_tp} distinct timepoints support a population trajectory; "
                          f"no higher level established"
                          + (f"; UNRESOLVED: {', '.join(pending)}" if pending else ""))
        return a

    if st("n_timepoints") == UNKNOWN:
        pending.append("n_timepoints")
    a.level = LEVEL_0 + (PENDING if pending else "")
    a.level_reason = (f"n_timepoints={n_tp} — no early→late mapping"
                      + (f"; UNRESOLVED: {', '.join(pending)}" if pending else ""))
    return a


# ---- per-dataset auditors -------------------------------------------------------------------- #
def audit_10x_timecourse(a: DatasetAudit) -> DatasetAudit:
    """A per-GSM 10x directory: GSE242423 and GSE297234."""
    p = Path(a.path)
    barcodes = sorted(p.glob("*barcodes.tsv.gz"))
    matrices = sorted(p.glob("*matrix.mtx.gz")) + sorted(p.glob("*.h5"))
    gsms = sorted({m.name.split("_")[0] for m in (barcodes + matrices)
                   if m.name.startswith("GSM")})
    a.add("n_gsm", len(gsms), PRESENT, f"{len(gsms)} GSM prefixes in filenames: {gsms[:4]}…")

    labels = sorted({re.sub(r"^GSM\d+_", "", b.name).split(".")[0] for b in barcodes + matrices})
    # Strip a donor prefix before counting timepoints. `GM00731_D0` and `GM23815_D0` are TWO
    # DONORS AT ONE TIMEPOINT, not two timepoints -- counting raw labels turns a D0-only corpus
    # into a fake timecourse and would wrongly license a forward task.
    donor_re = re.compile(r"^(GM\d+|[NOY]\d+)_")
    donors_10x = sorted({m.group(1) for x in labels if (m := donor_re.match(x))})
    states = sorted({donor_re.sub("", x).replace("_filtered_feature_bc_matrix", "")
                     for x in labels})
    a.add("timepoint_labels", labels, PRESENT, f"parsed from filenames: {labels}")
    a.add("donors", donors_10x, PRESENT if donors_10x else ABSENT_PROVEN,
          f"donor prefixes in filenames: {donors_10x}" if donors_10x
          else "no donor prefix in filenames — single line")
    a.add("n_timepoints", len(states), PRESENT,
          f"{len(states)} distinct timepoint/state tokens after stripping donor prefix: "
          f"{states}"
          + (f" (from {len(labels)} labels across {len(donors_10x)} donors)"
             if donors_10x else ""))

    # replicate structure: one GSM per label means no replication IN THE LOCAL FILE SET
    if len(gsms) and len(labels) and len(gsms) == len(labels):
        a.add("replicates_per_timepoint", 1, PRESENT,
              f"{len(gsms)} GSMs for {len(labels)} distinct labels — one sample per timepoint")
    else:
        a.add("replicates_per_timepoint", None, UNKNOWN,
              f"{len(gsms)} GSMs vs {len(labels)} labels — mapping not determinable from names")

    # barcode ENCODING is provable; study-level lineage tracing is not
    if barcodes:
        enc, why = barcode_encoding(barcodes[0])
        a.add("barcode_encoding", enc, PRESENT, f"{barcodes[0].name}: {why}")
        if enc == "PLAIN_10X":
            a.add("clone_id_in_barcode", False, ABSENT_PROVEN,
                  f"{barcodes[0].name} barcodes are plain 10x; no lineage token is encoded")
        else:
            a.add("clone_id_in_barcode", None, UNKNOWN,
                  f"{barcodes[0].name} encoding {enc} — needs inspection")
    else:
        a.add("barcode_encoding", None, UNKNOWN, "no barcodes file found locally")
        a.add("clone_id_in_barcode", None, UNKNOWN, "no barcodes file found locally")

    sm = sorted(p.glob("*series_matrix*"))
    if sm:
        f = series_matrix_fields(sm[0])
        chars = f.get("!Sample_characteristics_ch1", [])
        a.add("series_matrix", sm[0].name, PRESENT, f"{len(chars)} characteristics fields")
        lineage_hit = [c for c in chars if re.search(r"clone|lineage|barcode|celltag|tag",
                                                     c, re.I)]
        if lineage_hit:
            a.add("lineage_link", True, PRESENT, f"series matrix mentions: {lineage_hit[:2]}")
        else:
            a.add("lineage_link", False, ABSENT_PROVEN,
                  f"{sm[0].name} carries {len(chars)} characteristics, none mentioning "
                  "clone/lineage/barcode/CellTag")
    else:
        # THE case rule 2 exists for
        a.add("series_matrix", None, UNKNOWN,
              "no *series_matrix* file on disk — sample-level metadata was never downloaded")
        a.add("lineage_link", None, UNKNOWN,
              "cannot be established: the series matrix that would carry clone/lineage "
              "annotation is not present locally. The plain-10x barcode encoding proves only "
              "that the BARCODE FILE carries no tag, not that the STUDY lacked lineage tracing")

    ortho = find_orthogonal_phenotype(labels, [x.name for x in p.iterdir()])
    if ortho:
        a.add("independent_outcome", True, PRESENT,
              f"surface-marker sort token in filenames: {ortho}")
        a.add("outcome_is_rna_surrogate", False, PRESENT, f"{ortho} is an antibody phenotype")
    else:
        a.add("independent_outcome", False, ABSENT_PROVEN,
              "directory holds only expression matrices/barcodes; no imaging, colony, sorting or "
              "viability file, and no surface-marker token in any filename. Every fate label in "
              "this project is computed from RNA (src/cellfate/data/labels.py::fate_labels)")
        a.add("outcome_is_rna_surrogate", True, PRESENT,
              "y_cls is derived by fate_labels() from marker-program scores on the cell's own "
              "expression vector")

    # independent units: one culture, so the trajectory is serially dependent
    a.add("n_independent_units", 1 if a.get("replicates_per_timepoint") and
          a.get("replicates_per_timepoint").value == 1 else None,
          PRESENT if a.get("replicates_per_timepoint") and
          a.get("replicates_per_timepoint").value == 1 else UNKNOWN,
          "one sample per timepoint from a single line implies ONE culture; timepoints from one "
          "culture are serially dependent, not independent units")
    a.add("cell_count_is_not_n", True, PRESENT,
          "plan §6: effective n counts independent donors/cultures/clones/experiments, not cells")
    return a


def parse_days(titles):
    """Day tokens appear as `d11`, `D10` AND `13days`. Missing one format silently reports a
    timecourse as having no timepoints, which then drives an INVALID verdict on a parsing bug."""
    days = set()
    for x in titles:
        for m in re.finditer(r"(?:^|[_\s])[dD](\d+)(?:[_\s]|$)", str(x)):
            days.add(int(m.group(1)))
        for m in re.finditer(r"(\d+)\s*days?", str(x), re.I):
            days.add(int(m.group(1)))
    return sorted(days)


# Antibody surface-marker sorts. These are measured by STAINING, not by the RNA vector, so a
# sample's marker identity is an ORTHOGONAL phenotype -- the only non-RNA outcome-like signal
# anywhere in the local corpora. Blanket-asserting "no independent outcome" would have closed
# off the most promising local lead.
ORTHOGONAL_MARKERS = ("SSEA4", "SSEA-4", "CD13", "TRA-1-60", "TRA160")


def find_orthogonal_phenotype(titles, chars):
    hits = sorted({m for x in list(titles) + list(chars)
                   for m in ORTHOGONAL_MARKERS if m.lower() in str(x).lower()})
    return hits


def detect_modality(fields, files):
    strat = " ".join(fields.get("!Sample_library_strategy", []))
    title = " ".join(fields.get("!Series_title", []))
    blob = f"{strat} {title} {' '.join(files)}".lower()
    if "rna-seq" in blob or "rnaseq" in blob:
        return "RNA", "library_strategy/title indicates RNA-Seq"
    if any(k in blob for k in ("array", "methyl", "epic", "450k")):
        return "METHYLATION_ARRAY", "title/strategy indicates a methylation array"
    return None, "modality not determinable from series matrix"


def audit_bulk_series(a: DatasetAudit) -> DatasetAudit:
    """A bulk corpus shipped as a matrix + series matrix: GSE165176/7/8/9, GSE113957."""
    p = Path(a.path)
    sm = sorted(p.glob("*series_matrix*"))
    mats = [f for f in p.glob("*.txt.gz") if "series_matrix" not in f.name] +            [f for f in p.glob("*.tsv.gz")]
    a.add("matrix_files", [m.name for m in mats], PRESENT, f"{len(mats)} matrix file(s)")

    if not sm:
        a.add("series_matrix", None, UNKNOWN, "no series matrix on disk")
        for k in ("n_samples", "n_timepoints", "lineage_link", "n_independent_units",
                  "independent_outcome", "modality"):
            a.add(k, None, UNKNOWN, "series matrix absent locally")
        return a

    f = series_matrix_fields(sm[0], keys=("!Sample_title", "!Sample_characteristics_ch1",
                                          "!Sample_geo_accession", "!Sample_source_name_ch1",
                                          "!Sample_library_strategy", "!Series_title"))
    titles = f.get("!Sample_title", [])
    chars = f.get("!Sample_characteristics_ch1", [])
    a.add("series_matrix", sm[0].name, PRESENT,
          f"{len(titles)} sample titles, {len(chars)} characteristics fields")
    a.add("n_samples", len(titles), PRESENT, f"{len(titles)} !Sample_title entries")

    mod, why = detect_modality(f, [m.name for m in mats])
    a.add("modality", mod, PRESENT if mod else UNKNOWN, why)

    days = parse_days(titles)
    a.add("n_timepoints", len(days) if days else None,
          PRESENT if days else UNKNOWN,
          f"day tokens parsed from sample titles: {days}" if days
          else "no day token matched `dN` or `N days` in sample titles")

    donors = sorted({m.group(1) for x in titles
                     if (m := re.match(r"([NOY]\d+|GM\d+)", str(x)))})
    exps = sorted({m.group(0) for x in titles
                   if (m := re.search(r"[eE]xp\d+", str(x)))})
    a.add("donors", donors, PRESENT if donors else UNKNOWN,
          f"{len(donors)} donor prefixes: {donors}" if donors else "no donor prefix pattern")
    a.add("experiment_replicates", exps, PRESENT if exps else ABSENT_PROVEN,
          f"experiment tokens in titles: {exps}" if exps
          else f"no expN token in {len(titles)} titles")
    n_units = (len(donors) * max(len(exps), 1)) if donors else None
    a.add("n_independent_units", n_units, PRESENT if n_units else UNKNOWN,
          f"{len(donors)} donors x {max(len(exps),1)} experiment(s) = {n_units}" if n_units
          else "no donor structure parsed")

    ortho = find_orthogonal_phenotype(titles, chars)
    if ortho:
        a.add("independent_outcome", True, PRESENT,
              f"surface-marker SORT present in sample titles: {ortho}. These are ANTIBODY "
              "phenotypes, measured independently of the RNA vector -- unlike y_cls, which "
              "fate_labels() computes from expression")
        a.add("outcome_is_rna_surrogate", False, PRESENT,
              f"{ortho} sorting is orthogonal to RNA")
        a.add("outcome_is_contemporaneous", True, PRESENT,
              "the sort is performed AT COLLECTION, so it is an orthogonal phenotype at the "
              "sample's own timepoint -- not yet a FUTURE outcome. Building a forward task "
              "requires pairing an early sample with a LATER sorted sample of the same culture")
    else:
        a.add("independent_outcome", False, ABSENT_PROVEN,
              f"no surface-marker or orthogonal-phenotype token in {len(titles)} titles / "
              f"{len(chars)} characteristics; only expression matrices present")
        a.add("outcome_is_rna_surrogate", True, PRESENT,
              "any fate/age label for this corpus is computed from its own expression")

    lineage_hit = [c for c in list(chars) + list(titles)
                   if re.search(r"clone|lineage|barcode|celltag", str(c), re.I)]
    if lineage_hit:
        a.add("lineage_link", True, PRESENT, f"mentions {lineage_hit[:2]}")
    else:
        a.add("lineage_link", False, ABSENT_PROVEN,
              f"{len(chars)} characteristics + {len(titles)} titles, none mentioning "
              "clone/lineage/barcode/CellTag; and bulk RNA has no per-cell identity to link")
    return a


def audit_dataset(name: str, candidates=None) -> DatasetAudit:
    p = locate(name, candidates)
    a = DatasetAudit(dataset=name, path=str(p) if p else None, present=p is not None)
    if not p:
        a.add("location", None, UNKNOWN,
              f"not found at any of {DEFAULT_LOCATIONS.get(name, [])} — "
              "NOT the same as proven absent")
        return classify(a)
    a.add("location", str(p), PRESENT, f"directory found at {p}")
    files = sorted(x.name for x in p.iterdir() if x.is_file())
    a.add("n_files", len(files), PRESENT, f"{len(files)} files")

    if any(x.endswith((".mtx.gz", ".h5")) for x in files):
        a = audit_10x_timecourse(a)
    else:
        a = audit_bulk_series(a)
    return classify(a)


def check_forbidden_predictors(columns) -> list[str]:
    """Plan §19.4 — an outcome-derived predictor is a hard error, not a warning."""
    return [c for c in columns
            if any(bad in str(c).lower() for bad in FORBIDDEN_PREDICTORS)]


def overall_verdict(audits) -> str:
    levels = {a.dataset: a.level for a in audits}
    if any(v.startswith(LEVEL_3) and PENDING not in v for v in levels.values()):
        return "STRICT_LINEAGE_AVAILABLE"
    if any(v.startswith(LEVEL_2) and PENDING not in v for v in levels.values()):
        return "CULTURE_FORWARD_AVAILABLE"
    if any(PENDING in v for v in levels.values()):
        return "TRAJECTORY_FORWARD_ONLY_PENDING_SOURCE_AUDIT"
    if any(v.startswith(LEVEL_1) for v in levels.values()):
        return "TRAJECTORY_FORWARD_ONLY"
    return "NO_VALID_FORWARD_TASK"


def run(datasets=None) -> dict:
    names = datasets or list(DEFAULT_LOCATIONS)
    audits = [audit_dataset(n) for n in names]
    return {
        "stage": 21, "phase": "data_audit",
        "rule_not_found_is_not_absent": True,
        "datasets": [
            {"dataset": a.dataset, "path": a.path, "present": a.present,
             "level": a.level, "level_reason": a.level_reason,
             "ruled_out": a.ruled_out,
             "findings": {k: asdict(v) for k, v in a.findings.items()}}
            for a in audits],
        "selected_dataset": None, "selected_task_level": None, "selected_forward_task": None,
        "verdict": overall_verdict(audits),
        "src_modified": False,
    }


def main() -> int:
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    r = run()
    print("\nSTAGE 21A — LOCAL DATASET GEOMETRY AUDIT")
    print("=" * 78)
    print("No model fitted. No expression values read. src/ unmodified.\n")
    print(f"{'dataset':<12}{'present':>8}  {'level':<42}{'units':>6}")
    print("-" * 78)
    for d in r["datasets"]:
        f = d["findings"]
        u = f.get("n_independent_units", {}).get("value")
        print(f"{d['dataset']:<12}{str(d['present']):>8}  {d['level']:<42}"
              f"{str(u if u is not None else '?'):>6}")

    print("\nWHY — evidence behind each classification")
    print("=" * 78)
    for d in r["datasets"]:
        print(f"\n{d['dataset']}  →  {d['level']}")
        print(f"   reason: {d['level_reason']}")
        for lvl, why in d["ruled_out"].items():
            print(f"   {lvl}: {why}")
        unknown = {k: v for k, v in d["findings"].items() if v["status"] == UNKNOWN}
        if unknown:
            print("   UNKNOWN (needs one more download, NOT proven absent):")
            for k, v in unknown.items():
                print(f"      {k}: {v['evidence']}")

    print("\n" + "=" * 78)
    print(f"VERDICT: {r['verdict']}")
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(r, indent=2), encoding="utf-8")
    print(f"saved -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
