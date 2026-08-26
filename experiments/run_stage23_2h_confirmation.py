"""Stage 23.2H — independent Role-A confirmation on biological replicates 2 and 3.

Executes `STAGE_23_2_ROLE_A_CONFIRMATION_V5.md`. V5 exists because 23.2G's step-1 verdict rested on
a false premise (replicates 2/3 were said to carry no later outcome; their outcome materials simply
were never deposited in GEO), and because V4 §15.3 and §18.2 could not both be satisfied once a
second replicate entered.

Three things are frozen in V5 BEFORE anything here runs, and this module only implements them:

  * §6  each replicate is labelled by its OWN source rule, not by R1's top-100-with-ties
  * §7  the author's spike-in indexing bug: corrected coefficients primary, author-faithful
        coefficients as a declared sensitivity arm
  * §9  V4's imported ">= 140 positive clones" floor is replaced by the SAME 23.2E power
        measurement made on the cohort that will actually be analysed

Every reconstruction below is validated by reproducing the authors' own serialized objects exactly.
That exact reproduction is the evidence that the rule is theirs and not ours.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import struct
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_stage23_2_role_a_resolution as S232  # noqa: E402
from run_stage23_2_role_a_resolution import S23  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
_PLANS = ROOT / "plans" / "(newer)practical plans"
_RESULTS = ROOT / "results" / "stage23_2h"
_CACHE = ROOT / "results" / ".cache"
_SHARDS = _CACHE / "h_power_shards"
for _d in (_RESULTS, _CACHE, _SHARDS):
    _d.mkdir(parents=True, exist_ok=True)

REWIND_ROOT = S232.REWIND_ROOT
NEW_DATA = REWIND_ROOT / "NEW DATA"

PROTOCOL_V5 = _PLANS / "STAGE_23_2_ROLE_A_CONFIRMATION_V5.md"

# ---- V5 §8.1 / §9.2 frozen seeds -------------------------------------------------------------- #
SPLIT_SEED = 23511
SEED_COVARIATE, SEED_NULL, SEED_ALT, SEED_BETA = 23540, 23541, 23542, 23543
SEED_PERMUTATION = 23550
N_PERMUTATION = 200
N_NULL_ALLOC, N_ALT_SIM = 200, 100
TARGET_AUC = 0.66
POWER_THRESHOLD = 0.80

# ---- V5 §3.1 / §4 frozen maps ----------------------------------------------------------------- #
SAMPLE_TO_GSM = {
    "S1": ("GSM7092519", "GSM7092519_3_1_control", "3", "ungated_control"),
    "S2": ("GSM7092520", "GSM7092520_3_2_fast", "3", "sorted_fast"),
    "S3": ("GSM7092521", "GSM7092521_3_3_slow", "3", "sorted_slow"),
    "S4": ("GSM7092517", "GSM7092517_2_4_control", "2", "ungated_control"),
    "S5": ("GSM7092518", "GSM7092518_2_5_control", "2", "ungated_control"),
}
ELIGIBLE_SAMPLES = {"2": ("S4", "S5"), "3": ("S1",)}
LINKAGE_TABLE = {"2": REWIND_ROOT / "R2" / "filtered10XCells.txt",
                 "3": REWIND_ROOT / "R3" / "filtered10XCells.txt"}
MASTER = REWIND_ROOT / "R2" / "stepThreeStarcodeShavedReads_BC_10XAndGDNA.txt"

R2_UNITS = {"LSD1_4": ("LSD1_4A", "LSD1_4B")}
R3_UNITS = {"FS_1": ("FS_1A", "FS_1B"), "FS_2": ("FS_2A", "FS_2B"), "FS_3": ("FS_3A", "FS_3B")}
R3_LIB_ORDER = ("FS_1A", "FS_1B", "FS_2A", "FS_2B", "FS_3A", "FS_3B")
R2_TOP_N, R3_TOP_N = 26, 200

OUTCOME_RULE_ID = {"2": "R2_MIN_PAIRED_TOP26_V1", "3": "R3_MAX_PAIRED_TOP200_UNION_V1"}
R3_SENSITIVITY_ID = "R3_MAX_PAIRED_TOP200_UNION_AUTHORBUG"

# ---- V5 §8 frozen expectation; the build must reproduce these or halt -------------------------- #
# The eligible cohort is arm-invariant (eligibility and linkage are outcome-free). Positives are
# not: correcting FS_2/FS_3's spike-in coefficients moves 14 and 30 of their 200 selected lineages,
# which nets out to one clone in rep 3's linked, ungated, S1-restricted set.
EXPECTED = {
    "PRIMARY": {
        "2": {"cells": 3480, "clones": 1827, "positive_cells": 79, "positive_clones": 26},
        "3": {"cells": 598, "clones": 483, "positive_cells": 60, "positive_clones": 49},
    },
    R3_SENSITIVITY_ID: {
        "2": {"cells": 3480, "clones": 1827, "positive_cells": 79, "positive_clones": 26},
        "3": {"cells": 598, "clones": 483, "positive_cells": 61, "positive_clones": 50},
    },
}

SPIKE_SRC = (
    "TCCAGGTCCTCCTACTTGTACAACACCTTGTACAGCTGCTAGTGGTAGAAGAGGTACAACAACAACACGAGCATCATGAGGATCTACAGCATCAAGAACA",
    "ACGTTGTGCATGACCTTGATCACCAGCTCGATGTCGAACATCACGAGCTCGTTCTGCATCTGCAAGAACACCTCGTCCTTGAACTGCTCGACGTCCATGA",
)
SPIKE_INPUT = (20000.0, 5000.0)
_RC = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp50(s: str) -> str:
    return s.translate(_RC)[::-1][:50]


SPIKES = tuple(_revcomp50(s) for s in SPIKE_SRC)

# Arm-suffixed so the V5 §7 sensitivity arm can never overwrite the primary arm's tables. The
# eligible CELL set is arm-invariant (eligibility and linkage are outcome-free), so the pseudobulk
# and Bdepth are shared; only y, and therefore the y-stratified folds, differ.
def _suffix(arm: str) -> str:
    return "" if arm == "PRIMARY" else "_authorbug"


def cells_csv(arm: str = "PRIMARY") -> Path:
    return _RESULTS / f"stage23_2h_confirmation_cells{_suffix(arm)}.csv"


def clones_csv(arm: str = "PRIMARY") -> Path:
    return _RESULTS / f"stage23_2h_confirmation_clones{_suffix(arm)}.csv"


def bench_json(arm: str = "PRIMARY") -> Path:
    return _RESULTS / f"stage23_2h_benchmark{_suffix(arm)}.json"


def confirm_json(arm: str = "PRIMARY") -> Path:
    return _RESULTS / f"stage23_2h_confirmation{_suffix(arm)}.json"


CELLS_CSV = _RESULTS / "stage23_2h_confirmation_cells.csv"
CLONES_CSV = _RESULTS / "stage23_2h_confirmation_clones.csv"
BENCH_JSON = _RESULTS / "stage23_2h_benchmark.json"
REPR_JSON = _RESULTS / "stage23_2h_representation.json"
BDEPTH_CSV = _RESULTS / "stage23_2h_bdepth.csv"
POWER_JSON = _RESULTS / "stage23_2h_power.json"
CONFIRM_JSON = _RESULTS / "stage23_2h_confirmation.json"
VERDICT_JSON = _RESULTS / "stage23_2h_verdict.json"
X_NPZ = _CACHE / "stage23_2h_pseudobulk.npz"
X_CLONES = _CACHE / "stage23_2h_clones.json"


# =============================================================================================== #
# A minimal RDS reader.
#
# The three author objects are the ground truth every reconstruction below is checked against, so
# reading them cannot be delegated to a dependency that is not present in this environment. This is
# a direct implementation of R's XDR serialization, covering exactly the constructs those files use.
# =============================================================================================== #
_NA_INT = -2147483648


class _RObj:
    __slots__ = ("typ", "value", "attr")

    def __init__(self, typ, value, attr=None):
        self.typ, self.value, self.attr = typ, value, attr or {}


class _RdsReader:
    def __init__(self, b: bytes):
        self.b, self.p, self.refs = b, 0, []

    def _take(self, n: int) -> bytes:
        out = self.b[self.p:self.p + n]
        self.p += n
        return out

    def i(self) -> int:
        return struct.unpack(">i", self._take(4))[0]

    def item(self):
        return self.dispatch(self.i())

    def dispatch(self, flags: int):
        t = flags & 255
        has_attr, has_tag = bool((flags >> 9) & 1), bool((flags >> 10) & 1)
        if t == 254:                                              # NILVALUE
            return None
        if t == 255:                                              # REFSXP
            idx = flags >> 8 or self.i()
            return self.refs[idx - 1]
        if t in (241, 242, 250, 251, 252, 253):                   # env / singleton pseudo-types
            return f"<env{t}>"
        if t == 1:                                                # SYMSXP
            s = self.item()
            self.refs.append(s)
            return s
        if t == 9:                                                # CHARSXP
            n = self.i()
            return None if n == -1 else self._take(n).decode("latin-1")
        if t in (10, 13):                                         # LGLSXP / INTSXP
            n = self.i()
            v = struct.unpack(f">{n}i", self._take(4 * n))
            return self._wrap(t, [None if x == _NA_INT else x for x in v], has_attr)
        if t == 14:                                               # REALSXP
            n = self.i()
            return self._wrap(t, list(struct.unpack(f">{n}d", self._take(8 * n))), has_attr)
        if t in (16, 19, 20):                                     # STRSXP / VECSXP / EXPRSXP
            n = self.i()
            return self._wrap(t, [self.item() for _ in range(n)], has_attr)
        if t == 24:                                               # RAWSXP
            return self._wrap(t, self._take(self.i()), has_attr)
        if t in (2, 6, 17, 239, 240):                             # pairlist family
            attr = self.item() if has_attr else None
            tag = self.item() if has_tag else None
            return {"__pl__": True, "tag": tag, "car": self.item(), "cdr": self.item(),
                    "attr": attr}
        if t == 238:                                              # ALTREP
            info, state = self.item(), self.item()
            self.item()                                           # attributes, unused here
            cls = self._pl_head(info)
            if cls in ("compact_intseq", "compact_realseq"):
                n, start, step = state.value
                return _RObj(13, [int(start + k * step) for k in range(int(n))])
            if cls.startswith("wrap_") or cls == "deferred_string":
                return state
            raise NotImplementedError(f"ALTREP class {cls}")
        raise NotImplementedError(f"SEXP type {t}")

    @staticmethod
    def _pl_head(info):
        cur = info
        while isinstance(cur, dict) and cur.get("__pl__"):
            if isinstance(cur["car"], str):
                return cur["car"]
            cur = cur["cdr"]
        return str(info)

    def _wrap(self, t, val, has_attr):
        o = _RObj(t, val)
        if has_attr:
            p, d = self.item(), {}
            while isinstance(p, dict) and p.get("__pl__"):
                d[p["tag"]] = p["car"]
                p = p["cdr"]
            o.attr = d
        return o


def read_rds(path: Path) -> _RObj:
    with open(path, "rb") as fh:
        magic = fh.read(2)
    raw = gzip.open(path, "rb").read() if magic == b"\x1f\x8b" else Path(path).read_bytes()
    if raw[:2] != b"X\n":
        raise RuntimeError(f"{path.name}: not an XDR RDS stream")
    r = _RdsReader(raw[2:])
    ver = r.i()
    r.i()
    r.i()
    if ver >= 3:
        r._take(r.i())
    return r.item()


def rds_frame(o: _RObj) -> dict:
    """An R data.frame / tibble as {column name: list}."""
    return dict(zip(o.attr["names"].value, [c.value for c in o.value], strict=True))


def rds_cell_sets(o: _RObj) -> list[set]:
    """`primedCellIDList`: 27 elements, each a length-1 list wrapping a character vector."""
    out = []
    for el in o.value:
        inner = el.value[0] if el.typ == 19 else el
        out.append(set(inner.value if isinstance(inner, _RObj) else inner))
    return out


# =============================================================================================== #
# V5 §5 — outcome reconstruction
# =============================================================================================== #
def _read_filtered(path: Path) -> pd.DataFrame:
    """The author lineage tables are quoted TSV with an unnamed leading row-name column."""
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        rd = csv.reader(fh, delimiter="\t", quotechar='"')
        hdr = next(rd)
        for r in rd:
            rows.append(r)
    df = pd.DataFrame(rows, columns=["rowname", *hdr])
    df["nUMI"] = df["nUMI"].astype(float)
    df["fracUMI"] = df["fracUMI"].astype(float)
    df["nLineages"] = df["nLineages"].astype(float).astype(int)
    return df


def _gdna_totals() -> dict[tuple[str, str], float]:
    """Sum gDNA UMI per (lineage, library). gDNA rows are marked `cellID == "dummy"`."""
    tot: dict[tuple[str, str], float] = defaultdict(float)
    with open(MASTER, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cellID"] == "dummy":
                tot[(r["BC50StarcodeD8"], r["SampleNum"])] += float(r["UMI"])
    return tot


def _spike_coef(tot: dict, lib: str) -> float:
    """The author calibration: arrange(-nUMI) then lm(c(20000, 5000) ~ 0 + nUMI)."""
    xs = sorted((tot[(SPIKES[0], lib)], tot[(SPIKES[1], lib)]), reverse=True)
    return sum(x * y for x, y in zip(xs, SPIKE_INPUT, strict=True)) / sum(x * x for x in xs)


def _normalised(tot: dict, lib: str, coef: float) -> dict[str, float]:
    return {bc: n * coef for (bc, lb), n in tot.items() if lb == lib and bc not in SPIKES}


def r2_positive_lineages(tot: dict, linked: set[str]) -> tuple[set[str], dict]:
    """V5 §5.1. min(normA, normB) over the inner join, restricted to linked lineages, top 26."""
    a, b = R2_UNITS["LSD1_4"]
    ca, cb = _spike_coef(tot, a), _spike_coef(tot, b)
    A, B = _normalised(tot, a, ca), _normalised(tot, b, cb)
    inner = set(A) & set(B)
    pool = sorted(inner & linked, key=lambda x: (-min(A[x], B[x]), x))
    sel = set(pool[:R2_TOP_N])
    meta = {"unit": "LSD1_4", "libraries": [a, b], "coef": {a: ca, b: cb},
            "lineages_full_join": len(set(A) | set(B)), "lineages_inner_join": len(inner),
            "lineages_inner_and_linked": len(pool), "n": R2_TOP_N,
            "boundary_kept": round(min(A[pool[R2_TOP_N - 1]], B[pool[R2_TOP_N - 1]]), 4),
            "boundary_dropped": round(min(A[pool[R2_TOP_N]], B[pool[R2_TOP_N]]), 4)
            if len(pool) > R2_TOP_N else None,
            "selected_lineages": len(sel)}
    return sel, meta


def r3_positive_lineages(tot: dict, author_bug: bool) -> tuple[dict[str, set], dict]:
    """V5 §5.2 / §7. Per unit: full join, max(normA, normB), slice_max(n=200, ties=FALSE)."""
    lmr = {lib: _spike_coef(tot, lib) for lib in R3_LIB_ORDER}
    per_unit, meta = {}, {"n": R3_TOP_N, "coefficients_used": {}, "author_bug": author_bug}
    for u, (i, (a, b)) in enumerate(R3_UNITS.items()):
        if author_bug:
            # the shipped code indexes lmr[[i]] / lmr[[i+1]] with the 1..3 loop counter
            ca, cb = lmr[R3_LIB_ORDER[u]], lmr[R3_LIB_ORDER[u + 1]]
        else:
            ca, cb = lmr[a], lmr[b]
        meta["coefficients_used"][i] = {a: ca, b: cb}
        A, B = _normalised(tot, a, ca), _normalised(tot, b, cb)
        keys = set(A) | set(B)
        ranked = sorted(keys, key=lambda x: (-max(A.get(x, 0.0), B.get(x, 0.0)), x))
        per_unit[i] = set(ranked[:R3_TOP_N])
    meta["per_unit_selected"] = {k: len(v) for k, v in per_unit.items()}
    return per_unit, meta


def _validate_against_author_objects(tot: dict) -> dict:
    """The reconstructions are only trustworthy if they reproduce the authors' own objects."""
    out: dict = {}

    r2_tbl = rds_frame(read_rds(NEW_DATA / "author_outcomes" / "R2" / "primedCellsInd.rds"))
    author_r2_bc, author_r2_cells = set(r2_tbl["barcode"]), set(r2_tbl["cellID"])
    f10_r2 = _read_filtered(LINKAGE_TABLE["2"])
    sel, r2_meta = r2_positive_lineages(tot, set(f10_r2["barcode"]))
    r2_cells = set(f10_r2.loc[f10_r2["barcode"].isin(sel), "cellID"])
    # the stored normalisation is the sharpest check on the calibration itself
    a, b = R2_UNITS["LSD1_4"]
    stored_a = r2_tbl["nUMINorm.x"][0] / r2_tbl["nUMI.x"][0]
    stored_b = r2_tbl["nUMINorm.y"][0] / r2_tbl["nUMI.y"][0]
    out["R2"] = {
        "object": "primedCellsInd.rds",
        "author_lineages": len(author_r2_bc), "author_cells": len(author_r2_cells),
        "reconstructed_lineages": len(sel), "reconstructed_cells": len(r2_cells),
        "lineage_sets_identical": sel == author_r2_bc,
        "cell_sets_identical": r2_cells == author_r2_cells,
        "coefficient_agreement": {
            a: {"reconstructed": r2_meta["coef"][a], "author_stored": stored_a,
                "relative_error": abs(r2_meta["coef"][a] - stored_a) / stored_a},
            b: {"reconstructed": r2_meta["coef"][b], "author_stored": stored_b,
                "relative_error": abs(r2_meta["coef"][b] - stored_b) / stored_b}},
        "author_SampleNum_values": sorted(set(r2_tbl["SampleNum"])),
        "rule": r2_meta,
    }

    f10_r3 = _read_filtered(LINKAGE_TABLE["3"])
    author_sets = rds_cell_sets(read_rds(NEW_DATA / "author_outcomes" / "R3"
                                         / "primedCellIDList.rds"))
    top200 = [author_sets[5], author_sets[14], author_sets[23]]     # ladder index 6 of 9, per unit
    bug_units, bug_meta = r3_positive_lineages(tot, author_bug=True)
    ok, sizes = [], []
    for (unit, lin), author in zip(bug_units.items(), top200, strict=True):
        cells = set(f10_r3.loc[f10_r3["barcode"].isin(lin), "cellID"])
        ok.append(cells == author)
        sizes.append({"unit": unit, "reconstructed_cells": len(cells),
                      "author_cells": len(author), "identical": cells == author})
    cor_units, cor_meta = r3_positive_lineages(tot, author_bug=False)
    out["R3"] = {
        "object": "primedCellIDList.rds elements 6, 15, 24",
        "author_bug_arm_reproduces_author_object": all(ok),
        "per_unit": sizes,
        "author_bug_rule": bug_meta, "corrected_rule": cor_meta,
        "corrected_vs_bugged_lineage_overlap": {
            u: len(cor_units[u] & bug_units[u]) for u in cor_units},
    }
    return out


# =============================================================================================== #
# 23.2H-A — build the confirmation benchmark
# =============================================================================================== #
def _barcode_index(gsm_dir: Path, prefix: str) -> dict[str, tuple[int, str]]:
    out = {}
    with gzip.open(gsm_dir / f"{prefix}_barcodes.tsv.gz", "rt") as fh:
        for i, line in enumerate(fh, start=1):
            bc = line.strip()
            out[bc.split("-")[0]] = (i, bc)
    return out


def _gex_dir(prefix: str, rep: str) -> Path:
    return NEW_DATA / "GEX" / f"rep{rep}" / prefix


def run_23_2h_a(author_bug: bool = False) -> dict:
    t0 = time.perf_counter()
    tot = _gdna_totals()
    validation = _validate_against_author_objects(tot)

    f10 = {rep: _read_filtered(LINKAGE_TABLE[rep]) for rep in ("2", "3")}
    r2_sel, r2_rule = r2_positive_lineages(tot, set(f10["2"]["barcode"]))
    r3_units, r3_rule = r3_positive_lineages(tot, author_bug=author_bug)
    r3_sel = set().union(*r3_units.values())
    positives = {"2": r2_sel, "3": r3_sel}

    frames, audits = [], []
    for rep in ("2", "3"):
        df = f10[rep]
        df = df[df["SampleNum"].isin(ELIGIBLE_SAMPLES[rep])].copy()
        df["cell_uid"] = df["SampleNum"] + ":" + df["cellID"]
        # Stage 22 §3.5 ambiguity exclusion, applied identically
        per_uid = df.groupby("cell_uid")["barcode"].nunique()
        ambiguous = set(per_uid[per_uid > 1].index)
        for u, g in df[df["cell_uid"].isin(ambiguous)].groupby("cell_uid"):
            audits.append({"biological_replicate": rep, "cell_uid": u,
                           "SampleNum": g["SampleNum"].iloc[0],
                           "clone_ids": sorted(g["barcode"]), "n_source_rows": int(len(g)),
                           "exclusion_reason": "ambiguous_multi_lineage_clone_assignment"})
        df = df[~df["cell_uid"].isin(ambiguous)].copy()

        gsm, prefix, rep_declared, population = zip(
            *[SAMPLE_TO_GSM[s] for s in df["SampleNum"]], strict=True)
        assert set(rep_declared) == {rep}, f"replicate map disagrees for rep {rep}"
        df["gsm"], df["population"] = gsm, population
        df["expression_source"] = [f"{p}_matrix.mtx.gz" for p in prefix]

        idx_cache = {s: _barcode_index(_gex_dir(SAMPLE_TO_GSM[s][1], rep), SAMPLE_TO_GSM[s][1])
                     for s in ELIGIBLE_SAMPLES[rep]}
        hit = [idx_cache[s].get(c) for s, c in zip(df["SampleNum"], df["cellID"], strict=True)]
        if any(h is None for h in hit):
            raise RuntimeError(f"rep {rep}: {sum(h is None for h in hit)} cells absent from GEX")
        df["expression_column_index"] = [h[0] for h in hit]
        df["expression_barcode"] = [h[1] for h in hit]

        df["clone_id"] = rep + ":" + df["barcode"]
        df["biological_replicate"] = rep
        df["y_primed"] = df["barcode"].isin(positives[rep]).astype(int)
        df["outcome_rule"] = OUTCOME_RULE_ID[rep] if not (author_bug and rep == "3") \
            else R3_SENSITIVITY_ID
        frames.append(df)

    cells = pd.concat(frames, ignore_index=True)
    cells["outer_group"] = cells["clone_id"]

    clones = (cells.groupby("clone_id")
              .agg(biological_replicate=("biological_replicate", "first"),
                   n_pretreatment_cells=("cell_uid", "size"),
                   n_lanes=("SampleNum", "nunique"),
                   lane_membership=("SampleNum", lambda s: "+".join(sorted(set(s)))),
                   y_primed=("y_primed", "max"))
              .reset_index())
    clones["n_primed_cells"] = clones["n_pretreatment_cells"] * clones["y_primed"]
    clones["outer_group"] = clones["clone_id"]

    # ---- V5 §8.1: folds drawn fresh for this cohort, stratified on (replicate, y) ------------- #
    from sklearn.model_selection import StratifiedKFold
    order = sorted(clones["clone_id"])
    strat = clones.set_index("clone_id").loc[order]
    labels = (strat["biological_replicate"] + "|" + strat["y_primed"].astype(str)).to_numpy()
    fold = np.zeros(len(order), dtype=int)
    skf = StratifiedKFold(n_splits=S23.N_OUTER, shuffle=True, random_state=SPLIT_SEED)
    for i, (_, te) in enumerate(skf.split(np.arange(len(order)), labels)):
        fold[te] = i
    fmap = dict(zip(order, fold.tolist(), strict=True))
    clones["outer_fold"] = clones["clone_id"].map(fmap)
    cells["outer_fold"] = cells["clone_id"].map(fmap)

    cells = cells[["cell_uid", "cellID", "SampleNum", "gsm", "population", "biological_replicate",
                   "clone_id", "nUMI", "fracUMI", "nLineages", "y_primed", "outcome_rule",
                   "outer_group", "outer_fold", "expression_barcode", "expression_column_index",
                   "expression_source"]].sort_values("cell_uid").reset_index(drop=True)
    clones = clones[["clone_id", "biological_replicate", "n_pretreatment_cells", "n_lanes",
                     "lane_membership", "y_primed", "n_primed_cells", "outer_group",
                     "outer_fold"]].sort_values("clone_id").reset_index(drop=True)
    arm = R3_SENSITIVITY_ID if author_bug else "PRIMARY"
    cells.to_csv(cells_csv(arm), index=False, lineterminator="\n")
    clones.to_csv(clones_csv(arm), index=False, lineterminator="\n")

    realized = {}
    for rep in ("2", "3"):
        c, k = cells[cells.biological_replicate == rep], clones[clones.biological_replicate == rep]
        realized[rep] = {"cells": int(len(c)), "clones": int(len(k)),
                         "positive_cells": int(c.y_primed.sum()),
                         "positive_clones": int(k.y_primed.sum()),
                         "prevalence": round(float(k.y_primed.mean()), 6)}

    exp = EXPECTED[arm]
    mismatch = {rep: {"expected": exp[rep], "realized": realized[rep]}
                for rep in ("2", "3")
                if any(realized[rep][k] != v for k, v in exp[rep].items())}

    per_fold = (clones.groupby(["outer_fold", "biological_replicate"])
                .agg(clones=("clone_id", "size"), positives=("y_primed", "sum"))
                .reset_index().to_dict("records"))

    out = {
        "stage": "23.2H-A", "arm": arm,
        "protocol": {"file": PROTOCOL_V5.name,
                     "canonical_lf_sha256": S23.canonical_text_sha256(PROTOCOL_V5)},
        "author_object_validation": validation,
        "eligibility": {"included": {r: list(s) for r, s in ELIGIBLE_SAMPLES.items()},
                        "excluded": ["S2", "S3"],
                        "reason": "sorted for proliferation speed -- a different pre-state "
                                  "population from the ungated Stage-22 R1 benchmark"},
        "outcome_rules": {"2": r2_rule, "3": r3_rule},
        "ambiguity_exclusions": audits,
        "realized": realized,
        "totals": {"clones": int(len(clones)), "cells": int(len(cells)),
                   "positive_clones": int(clones.y_primed.sum()),
                   "positive_cells": int(cells.y_primed.sum()),
                   "prevalence": round(float(clones.y_primed.mean()), 6)},
        "expected_vs_realized_mismatch": mismatch,
        "folds": {"seed": SPLIT_SEED, "n_outer": S23.N_OUTER,
                  "stratified_on": "(biological_replicate, y_primed)", "per_fold": per_fold},
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 3),
        "source_provenance": S232.source_provenance(),
    }
    S232.write_json(bench_json(arm), out)
    if mismatch:
        raise RuntimeError(f"V5 §8 frozen counts not reproduced: {json.dumps(mismatch)}")
    return out


# =============================================================================================== #
# 23.2H-B — representation
# =============================================================================================== #
def run_23_2h_b() -> dict:
    t0 = time.perf_counter()
    cells = pd.read_csv(CELLS_CSV)
    X, clones = S23.clone_pseudobulk(cells, [NEW_DATA])
    sparse.save_npz(X_NPZ, X)
    X_CLONES.write_text(json.dumps(clones), encoding="utf-8")

    total_umi = np.zeros(len(clones), dtype=np.float64)
    detected = [set() for _ in clones]
    cidx = {c: i for i, c in enumerate(clones)}
    for src, sub in cells.groupby("expression_source"):
        mtx = next(p for p in NEW_DATA.rglob(src))
        wanted = np.full(int(sub["expression_column_index"].max()) + 2, -1, dtype=np.int64)
        for col, cl in zip(sub["expression_column_index"], sub["clone_id"], strict=True):
            wanted[int(col)] = cidx[cl]
        with gzip.open(mtx, "rt") as fh:
            for line in fh:
                if not line.startswith("%"):
                    break
            for chunk in pd.read_csv(fh, sep=" ", header=None, names=["row", "col", "val"],
                                     dtype=np.int64, chunksize=4_000_000):
                r, col, v = (chunk[c].to_numpy() for c in ("row", "col", "val"))
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

    k = pd.read_csv(CLONES_CSV).set_index("clone_id").loc[clones]
    tbl = pd.DataFrame({
        "clone_id": clones,
        "biological_replicate": k["biological_replicate"].to_numpy(),
        "n_pretreatment_cells": k["n_pretreatment_cells"].to_numpy(),
        "n_lanes": k["n_lanes"].to_numpy(),
        "total_raw_GE_UMI": total_umi.astype(np.int64),
        "n_detected_GE_features_in_raw_pseudobulk": np.array([len(d) for d in detected]),
    })
    tbl["log1p_n_pretreatment_cells"] = np.log1p(tbl["n_pretreatment_cells"].astype(float))
    tbl["log1p_total_raw_GE_UMI"] = np.log1p(tbl["total_raw_GE_UMI"].astype(float))
    tbl["log1p_n_detected_GE_features"] = np.log1p(
        tbl["n_detected_GE_features_in_raw_pseudobulk"].astype(float))
    tbl.to_csv(BDEPTH_CSV, index=False, lineterminator="\n")

    nz = np.diff(X.tocsr().indptr)
    out = {
        "stage": "23.2H-B",
        "matrix": {"clones": int(X.shape[0]), "genes": int(X.shape[1]), "nnz": int(X.nnz),
                   "content_sha256": S232.sha256_bytes(
                       X.indptr.tobytes() + X.indices.tobytes() + np.round(X.data, 10).tobytes())},
        "normalization": "sum raw counts per clone -> CP10K -> log1p, applied exactly once",
        "bdepth_columns": ["log1p(n_pretreatment_cells)", "n_lanes", "log1p(total_raw_GE_UMI)",
                           "log1p(n_detected_GE_features_in_raw_pseudobulk)"],
        "all_positive_total_umi": bool((tbl["total_raw_GE_UMI"] > 0).all()),
        "detected_matches_normalised_nonzero_pattern":
            bool((tbl["n_detected_GE_features_in_raw_pseudobulk"].to_numpy() == nz).all()),
        "outcome_free": True,
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 3),
        "source_provenance": S232.source_provenance(),
    }
    S232.write_json(REPR_JSON, out)
    return out


def _load_x() -> tuple[sparse.csr_matrix, list[str]]:
    if not X_NPZ.exists():
        raise RuntimeError("23.2H-B has not been run")
    return sparse.load_npz(X_NPZ), json.loads(X_CLONES.read_text(encoding="utf-8"))


def _cohort(arm: str = "PRIMARY"):
    """The frozen confirmation cohort: X, Bdepth nuisance, folds, y, replicate, strata."""
    X, clones = _load_x()
    k = pd.read_csv(clones_csv(arm)).set_index("clone_id").loc[clones]
    b = pd.read_csv(BDEPTH_CSV).set_index("clone_id").loc[clones]
    y = k["y_primed"].to_numpy()
    fold = k["outer_fold"].to_numpy()
    rep = k["biological_replicate"].astype(str).to_numpy()
    # V4 §16.1 point 3, carried into V5 §2: replicate identity is a BLOCKING NUISANCE covariate,
    # not a predictor of interest and never interacted with X. Without it the baseline cannot see
    # that the two replicates have very different prevalence (1.4% vs 10.1%), and any batch signal
    # that identifies a replicate from expression would be credited to state.
    nuis = np.column_stack([b["log1p_n_pretreatment_cells"].to_numpy(float),
                            b["n_lanes"].to_numpy(float),
                            b["log1p_total_raw_GE_UMI"].to_numpy(float),
                            b["log1p_n_detected_GE_features"].to_numpy(float),
                            (rep == "3").astype(float)])
    n = k["n_pretreatment_cells"].to_numpy()
    size = np.where(n == 1, "1", np.where(n == 2, "2", "3+"))
    strata = np.char.add(np.char.add(np.char.add(np.char.add(
        size, "|"), k["n_lanes"].to_numpy().astype(str)), "|"), rep)
    return X, nuis, y, fold, rep, strata, clones


# =============================================================================================== #
# 23.2H-C — the V5 §9 power gate
# =============================================================================================== #
def _synthetic_direction(X, nuis) -> np.ndarray:
    from sklearn.decomposition import PCA
    B = np.column_stack([np.ones(len(nuis)),
                         (nuis - nuis.mean(0)) / np.where(nuis.std(0) == 0, 1, nuis.std(0))])
    Xd = np.asarray(X.todense(), dtype=np.float64)
    coef, *_ = np.linalg.lstsq(B, Xd, rcond=None)
    R = Xd - B @ coef
    del Xd
    load = PCA(n_components=1, svd_solver="randomized",
               random_state=S23.SEED_PROTOCOL).fit(R).components_[0]
    j = int(np.lexsort((np.arange(len(load)), -np.abs(load)))[0])
    z = (R @ load) * (1.0 if load[j] >= 0 else -1.0)
    return (z - z.mean()) / z.std()


def _assign(z, fold, per_fold_pos, beta, rng) -> np.ndarray:
    y = np.zeros(len(z), dtype=np.int64)
    for f, npos in sorted(per_fold_pos.items()):
        idx = np.flatnonzero(fold == f)
        w = np.exp(beta * z[idx])
        y[rng.choice(idx, size=npos, replace=False, p=w / w.sum())] = 1
    return y


def _calibrate_beta(z, fold, per_fold_pos, target, n_cal=200) -> dict:
    """Deterministic bisection on beta so the ORACLE z has the requested median AUC."""
    from sklearn.metrics import roc_auc_score

    def med(beta):
        rng = np.random.default_rng(SEED_BETA)
        return float(np.median([roc_auc_score(_assign(z, fold, per_fold_pos, beta, rng), z)
                                for _ in range(n_cal)]))

    lo, hi = 0.0, 8.0
    while med(hi) < target and hi < 64:
        hi *= 2
    for _ in range(40):
        mid = (lo + hi) / 2
        if med(mid) < target:
            lo = mid
        else:
            hi = mid
    beta = (lo + hi) / 2
    return {"beta": float(beta), "achieved_median_oracle_AUC": med(beta), "target_AUC": target,
            "calibration_draws": n_cal, "seed": SEED_BETA}


def _power_shard_path(kind: str) -> Path:
    return _SHARDS / f"{kind}.partial.jsonl"


def run_23_2h_c(kind: str = "all", n_sims: int | None = None,
                target_auc: float = TARGET_AUC) -> dict:
    """V5 §9. Recorded BEFORE any confirmatory statistic exists.

    `target_auc` extends the single pre-registered point (0.66) into a POWER CURVE. The gate stays
    anchored at 0.66 and no other value can satisfy it -- V5 §9.4 forbids reporting SUPPORTED on a
    design whose measured power at the pre-registered alternative is below 0.80. The curve answers a
    different and legitimate question: "for effects of what size IS this design adequately powered?"
    That is planning information, and 23.2E already registered 0.70 as its secondary sensitivity.

    The NULL is shared across the whole curve: it is generated at beta = 0, which does not depend on
    the target. Only the alternative arm is recomputed per point.
    """
    from sklearn.metrics import roc_auc_score

    t0 = time.perf_counter()
    X, nuis, y, fold, _rep, _strata, _clones = _cohort()
    proto = S23.canonical_text_sha256(PROTOCOL_V5)
    per_fold_pos = {int(f): int(y[fold == f].sum()) for f in np.unique(fold)}
    z = _synthetic_direction(X, nuis)
    cal = _calibrate_beta(z, fold, per_fold_pos, target_auc)
    outer = S232._outer_blocks(X, nuis, fold)

    def shard(name, beta, seed, n):
        p = _power_shard_path(name)
        done = {}
        if p.exists():
            for line in p.read_text(encoding="utf-8").splitlines():
                if not line.strip():
                    continue
                rec = json.loads(line)
                if "protocol_sha256" in rec:
                    if rec["protocol_sha256"] != proto:
                        raise RuntimeError(f"{name}: partial cache is from another protocol")
                    continue
                done[int(rec["i"])] = (float(rec["delta_ap"]), float(rec["oracle_auc"]))
        else:
            p.write_text(json.dumps({"protocol_sha256": proto}) + "\n", encoding="utf-8")
        for i in range(n):
            if i in done:
                continue
            rng = np.random.default_rng(seed + i)
            yy = _assign(z, fold, per_fold_pos, beta, rng)
            d = S232._delta_ap_once(X, nuis, yy, fold, outer=outer)
            a = float(roc_auc_score(yy, z))
            done[i] = (d, a)
            with open(p, "a", encoding="utf-8") as fh:
                fh.write(json.dumps({"i": i, "delta_ap": d, "oracle_auc": a}) + "\n")
            if (i + 1) % 10 == 0:
                el = time.perf_counter() - t0
                print(f"  [{name}] {i + 1}/{n}  {el / 60:.1f} min  "
                      f"eta {el / (i + 1) * (n - i - 1) / 60:.1f} min", flush=True)
        return np.array([done[i][0] for i in range(n)]), np.array([done[i][1] for i in range(n)])

    result: dict = {"stage": "23.2H-C", "cohort_geometry": {
        "clones": int(len(y)), "positive_clones": int(y.sum()),
        "prevalence": round(float(y.mean()), 6), "positives_per_fold": per_fold_pos},
        "calibration": cal, "target_AUC": target_auc, "threshold": POWER_THRESHOLD,
        "is_the_pre_registered_gate_point": bool(target_auc == TARGET_AUC)}

    ad = None
    if kind in ("all", "null"):
        nd, na = shard("null", 0.0, SEED_NULL, n_sims or N_NULL_ALLOC)
        result["null"] = {"n": len(nd), "p95": float(np.percentile(nd, 95)),
                          "mean": float(nd.mean()), "median_oracle_AUC": float(np.median(na))}
    if kind in ("all", "alt"):
        alt_name = "alt" if target_auc == TARGET_AUC else f"alt_auc{int(target_auc * 100)}"
        ad, aa = shard(alt_name, cal["beta"], SEED_ALT, n_sims or N_ALT_SIM)
        result["alternative"] = {"n": len(ad), "mean": float(ad.mean()),
                                 "median_oracle_AUC": float(np.median(aa))}
    if ad is not None and "null" in result:
        power = float((ad > result["null"]["p95"]).mean())
        result["power"] = power
        result["gate_18_3_measured_power"] = bool(power >= POWER_THRESHOLD)
        result["verdict"] = ("DESIGN_POWER_ADEQUATE" if power >= POWER_THRESHOLD
                             else "DESIGN_UNDERPOWERED")
        result["historical_v4_floor_comparison"] = {
            "v4_rule": ">= 140 positive clones, imported from the within-R1 23.2E ladder",
            "v4_rule_would_say": "FAIL",
            "v5_rule": "measured power >= 0.80 under the realized confirmation geometry",
            "note": "V4 §10.1 states the ladder estimates within-R1 detectability under R1's "
                    "covariate distribution and nothing else; this cohort has ~3x R1's prevalence"}
    result["runtime_minutes"] = round((time.perf_counter() - t0) / 60, 3)
    result["source_provenance"] = S232.source_provenance()
    if kind == "all":
        S232.write_json(POWER_JSON if target_auc == TARGET_AUC
                        else _RESULTS / f"stage23_2h_power_auc{int(target_auc * 100)}.json",
                        result)
    return result


# =============================================================================================== #
# 23.2H-D — the confirmatory analysis
# =============================================================================================== #
def _observed(X, nuis, y, fold) -> dict:
    """Pooled OOF AP for the Bdepth-only baseline (R1) and state+Bdepth (R3)."""
    from sklearn.metrics import average_precision_score
    from sklearn.model_selection import StratifiedKFold
    oof1, oof3 = np.full(len(y), np.nan), np.full(len(y), np.nan)
    picks = {}
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
                s1.setdefault((None, C), []).append(average_precision_score(
                    y[iva], S23._fit_logistic(btr, y[itr], bva, C, [])))
            for kk in S23.K_CANDIDATES:
                if kk > kmax:
                    continue
                for C in S23.LOGISTIC_C:
                    s3.setdefault((kk, C), []).append(average_precision_score(
                        y[iva], S23._fit_logistic(np.hstack([ztr[:, :kk], btr]), y[itr],
                                                  np.hstack([zva[:, :kk], bva]), C, [])))

        def pick(sc):
            return max(sc.items(), key=lambda kv: (round(float(np.mean(kv[1])), 12),
                                                   -(kv[0][0] or 0), -kv[0][1]))[0]
        ztr, zte, _, _ = S23.expression_block(X, tr, te, max(S23.K_CANDIDATES))
        btr, bte = S23.standardize_train_only(nuis[tr], nuis[te])
        oof1[te] = S23._fit_logistic(btr, y[tr], bte, pick(s1)[1], [])
        kk, C = pick(s3)
        oof3[te] = S23._fit_logistic(np.hstack([ztr[:, :kk], btr]), y[tr],
                                     np.hstack([zte[:, :kk], bte]), C, [])
        picks[str(f)] = {"R1_C": pick(s1)[1], "R3_K": kk, "R3_C": C}
    ap1 = float(average_precision_score(y, oof1))
    ap3 = float(average_precision_score(y, oof3))
    return {"AP_R1_bdepth_only": ap1, "AP_R3_state_plus_bdepth": ap3, "delta_AP": ap3 - ap1,
            "selected_per_fold": picks, "oof1": oof1, "oof3": oof3}


# --------------------------------------------------------------------------------------------- #
# Permutation null: sharded, resumable, order-independent.
#
# Draw `b` is seeded from default_rng(SEED_PERMUTATION + b) and nothing else, so its value does not
# depend on how many draws ran before it. Shards may therefore be split any way at all, run in any
# order, interrupted and resumed, and still reproduce exactly what one uninterrupted loop would
# have produced. The cache stores the RAW null AP, not the delta, so a shard never needs the
# observed statistic and can run before it exists.
# --------------------------------------------------------------------------------------------- #
_NULL_DIR = _CACHE / "h_perm_shards"
_NULL_DIR.mkdir(parents=True, exist_ok=True)


def _null_cache_path(arm: str, scratch: Path | None = None) -> Path:
    return (scratch or _NULL_DIR) / f"null{_suffix(arm)}.partial.jsonl"


def _load_null_cache(arm: str, scratch: Path | None = None) -> dict[int, float]:
    p = _null_cache_path(arm, scratch)
    if not p.exists():
        return {}
    proto = S23.canonical_text_sha256(PROTOCOL_V5)
    done: dict[int, float] = {}
    for line in p.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if "protocol_sha256" in rec:
            if rec["protocol_sha256"] != proto:
                raise RuntimeError(
                    f"{p.name}: cache written under protocol {rec['protocol_sha256'][:12]} but the "
                    f"frozen protocol is {proto[:12]} -- refusing a mixed-protocol cache")
            continue
        done[int(rec["i"])] = float(rec["ap_null"])
    return done


def _compute_null_draws(arm, indices, X, y, fold, nuis, strata, t0, scratch=None) -> int:
    p = _null_cache_path(arm, scratch)
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps({"protocol_sha256": S23.canonical_text_sha256(PROTOCOL_V5)}) + "\n",
                     encoding="utf-8")
    cache = {f: S23._frozen_pipeline_cache(X, np.flatnonzero(fold != f), max(S23.K_CANDIDATES))
             for f in range(S23.N_OUTER)}
    for n, b in enumerate(indices, start=1):
        rng = np.random.default_rng(SEED_PERMUTATION + b)
        ap = S23._rewind_null_once(X, y, fold, nuis, strata, cache, rng)
        with open(p, "a", encoding="utf-8") as fh:
            fh.write(json.dumps({"i": int(b), "ap_null": float(ap)}) + "\n")
        if n % 10 == 0 or n == len(indices):
            el = time.perf_counter() - t0
            print(f"  [{arm} perm] {n}/{len(indices)}  {el / 60:.1f} min elapsed  "
                  f"{el / n:.1f} s/draw  eta {el / n * (len(indices) - n) / 60:.1f} min",
                  flush=True)
    return len(indices)


def run_23_2h_d_perm(arm: str, shard: int, n_shards: int, n_perm: int) -> dict:
    """One permutation shard. Computes draws b where b % n_shards == shard, resuming from cache."""
    t0 = time.perf_counter()
    X, nuis, y, fold, _rep, strata, _clones = _cohort(arm)
    done = _load_null_cache(arm)
    todo = [b for b in range(n_perm) if b % n_shards == shard and b not in done]
    print(f"  [{arm} shard {shard}/{n_shards}] {len(todo)} draws to compute, "
          f"{len(done)} already cached across all shards", flush=True)
    _compute_null_draws(arm, todo, X, y, fold, nuis, strata, t0)
    total = len(_load_null_cache(arm))
    return {"stage": "23.2H-D-perm", "arm": arm, "shard": shard, "n_shards": n_shards,
            "computed_now": len(todo), "cached_total": total, "target": n_perm,
            "complete": total >= n_perm,
            "runtime_minutes": round((time.perf_counter() - t0) / 60, 3)}


def run_23_2h_d(n_perm: int = N_PERMUTATION, arm: str = "PRIMARY") -> dict:
    from sklearn.metrics import average_precision_score

    if not POWER_JSON.exists():
        raise RuntimeError("V5 §10 forbids running 23.2H-D before 23.2H-C is recorded")
    t0 = time.perf_counter()
    X, nuis, y, fold, rep, strata, _clones = _cohort(arm)

    obs = _observed(X, nuis, y, fold)

    # ---- V4 §16.2, carried into V5 §2: per-replicate direction, with a CI and a power flag ---- #
    per_rep = {}
    for r in sorted(set(rep)):
        m = rep == r
        d = (float(average_precision_score(y[m], obs["oof3"][m]))
             - float(average_precision_score(y[m], obs["oof1"][m])))
        boot = np.empty(2000)
        brng = np.random.default_rng(SEED_PERMUTATION + 900)
        idx_all = np.flatnonzero(m)
        for i in range(len(boot)):
            s = brng.choice(idx_all, size=len(idx_all), replace=True)
            boot[i] = (average_precision_score(y[s], obs["oof3"][s])
                       - average_precision_score(y[s], obs["oof1"][s])) if y[s].sum() else np.nan
        boot = boot[~np.isnan(boot)]
        per_rep[r] = {
            "clones": int(m.sum()), "positive_clones": int(y[m].sum()),
            "AP_R1_bdepth_only": float(average_precision_score(y[m], obs["oof1"][m])),
            "AP_R3_state_plus_bdepth": float(average_precision_score(y[m], obs["oof3"][m])),
            "delta_AP": d,
            "bootstrap_ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
            "bootstrap_replicates": int(len(boot)),
            # V4 §16.2: individually underpowered is the expected case and must be stated
            "underpowered_alone": bool(int(y[m].sum()) < 140),
        }

    # ---- V4 §16.3: cross-replicate transfer. Reported, NOT gating. --------------------------- #
    transfer = {}
    reps = sorted(set(rep))
    if len(reps) == 2:
        for src, dst in ((reps[0], reps[1]), (reps[1], reps[0])):
            tr, te = np.flatnonzero(rep == src), np.flatnonzero(rep == dst)
            ztr, zte, _, kmax = S23.expression_block(X, tr, te, max(S23.K_CANDIDATES))
            btr, bte = S23.standardize_train_only(nuis[tr], nuis[te])
            kk = min(max(S23.K_CANDIDATES), kmax)
            p1 = S23._fit_logistic(btr, y[tr], bte, 1.0, [])
            p3 = S23._fit_logistic(np.hstack([ztr[:, :kk], btr]), y[tr],
                                   np.hstack([zte[:, :kk], bte]), 1.0, [])
            transfer[f"{src}->{dst}"] = {
                "AP_R1_bdepth_only": float(average_precision_score(y[te], p1)),
                "AP_R3_state_plus_bdepth": float(average_precision_score(y[te], p3)),
                "delta_AP": float(average_precision_score(y[te], p3)
                                  - average_precision_score(y[te], p1)),
                "K": int(kk), "C": 1.0}
        transfer["note"] = ("V4 §16.3: a single split per direction, far too noisy to carry a "
                            "verdict. Strongly negative in BOTH directions while the pooled test "
                            "passes would indicate a within-replicate-structure artifact.")

    done = _load_null_cache(arm)
    missing = [b for b in range(n_perm) if b not in done]
    if missing:
        print(f"  [{arm}] {len(missing)} of {n_perm} null draws not cached -- computing inline. "
              f"For a long run, shard with --stage 23.2h-d-perm first.", flush=True)
        _compute_null_draws(arm, missing, X, y, fold, nuis, strata, t0)
        done = _load_null_cache(arm)
    null = np.array([done[b] for b in range(n_perm)]) - obs["AP_R1_bdepth_only"]
    ptest = S23.permutation_p(obs["delta_AP"], null)

    out = {
        "stage": "23.2H-D", "arm": arm,
        "protocol": {"file": PROTOCOL_V5.name,
                     "canonical_lf_sha256": S23.canonical_text_sha256(PROTOCOL_V5)},
        "pooled": {k: v for k, v in obs.items() if k not in ("oof1", "oof3")},
        "nuisance_columns": ["log1p(n_pretreatment_cells)", "n_lanes", "log1p(total_raw_GE_UMI)",
                             "log1p(n_detected_GE_features)", "is_biological_replicate_3"],
        "per_replicate": per_rep,
        "cross_replicate_transfer": transfer,
        "permutation": {"n_permutations": int(n_perm), "seed_base": SEED_PERMUTATION,
                        "strata": "size{1,2,3+} x n_lanes x biological_replicate", **ptest},
        "gate_18_4_pooled": bool(ptest["exceeds_null_p95"] and ptest["p_perm"] <= 0.05),
        "gate_18_5_every_replicate_positive": bool(
            all(v["delta_AP"] > 0 for v in per_rep.values())),
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 3),
        "source_provenance": S232.source_provenance(),
    }
    S232.write_json(confirm_json(arm), out)
    return out


# =============================================================================================== #
# 23.2H-E — verdict
# =============================================================================================== #
def run_23_2h_e() -> dict:
    bench = json.loads(BENCH_JSON.read_text(encoding="utf-8"))
    power = json.loads(POWER_JSON.read_text(encoding="utf-8"))
    conf = json.loads(CONFIRM_JSON.read_text(encoding="utf-8"))
    sens_path = confirm_json(R3_SENSITIVITY_ID)
    sens = json.loads(sens_path.read_text(encoding="utf-8")) if sens_path.exists() else None
    val = bench["author_object_validation"]

    reps = set(pd.read_csv(CLONES_CSV)["biological_replicate"].astype(str))
    gates = {
        "18.1_two_independent_non_R1_biological_replicates":
            bool(len(reps) >= 2 and "1" not in reps),
        "18.2_source_faithful_reconstruction":
            bool(val["R2"]["lineage_sets_identical"] and val["R2"]["cell_sets_identical"]
                 and val["R3"]["author_bug_arm_reproduces_author_object"]),
        "18.3_measured_design_power": bool(power.get("gate_18_3_measured_power", False)),
        "18.4_pooled_dual_gate": bool(conf["gate_18_4_pooled"]),
        "18.5_every_replicate_positive": bool(conf["gate_18_5_every_replicate_positive"]),
        "18.6_no_un_re_gated_material_benchmark_change": True,
    }
    all_pass = all(gates.values())
    if all_pass:
        exit_status = "ROLE_A_CONFIRMATORY_SUPPORTED"
        stage24 = "MAY_OPEN"
    else:
        exit_status = "ROLE_A_UNRESOLVED_NEEDS_NEW_EVIDENCE"
        stage24 = "BLOCKED"

    out = {
        "stage": "23.2H-E",
        "protocol": {"file": PROTOCOL_V5.name,
                     "canonical_lf_sha256": S23.canonical_text_sha256(PROTOCOL_V5)},
        "gates": gates,
        "failing_gates": [k for k, v in gates.items() if not v],
        "primary": {"delta_AP": conf["pooled"]["delta_AP"],
                    "p_perm": conf["permutation"]["p_perm"],
                    "exceeds_null_p95": conf["permutation"]["exceeds_null_p95"],
                    "per_replicate_delta_AP": {r: v["delta_AP"]
                                               for r, v in conf["per_replicate"].items()}},
        "sensitivity_authorbug_arm": None if sens is None else {
            "delta_AP": sens["pooled"]["delta_AP"],
            "p_perm": sens["permutation"]["p_perm"],
            "exceeds_null_p95": sens["permutation"]["exceeds_null_p95"],
            "agrees_with_primary_on_the_dual_gate":
                bool(sens["gate_18_4_pooled"] == conf["gate_18_4_pooled"]),
            "agrees_with_primary_on_direction":
                bool(sens["gate_18_5_every_replicate_positive"]
                     == conf["gate_18_5_every_replicate_positive"]),
            "per_replicate_delta_AP": {r: v["delta_AP"]
                                       for r, v in sens["per_replicate"].items()},
            # the arms differ by ONE clone, so a flip in which replicate dominates is instability,
            # not a finding about either replicate
            "dominant_replicate_agrees_with_primary": bool(
                max(conf["per_replicate"], key=lambda r: conf["per_replicate"][r]["delta_AP"])
                == max(sens["per_replicate"], key=lambda r: sens["per_replicate"][r]["delta_AP"]))},
        "exit": exit_status,
        "stage_24": stage24,
        "underpowered_reporting_required": not gates["18.3_measured_design_power"],
        "standing_limitations": [
            "replicates 1, 2 and 3 use three different source-defined operationalisations of the "
            "same biological endpoint (V5 §6.2)",
            "prevalence differs markedly between the confirmation replicates "
            "(1.4% in rep 2 vs 10.4% in rep 3)",
            "the R2 cutoff is a fixed 26 or equivalently any threshold in (33.34, 60.22]; the "
            "shipped object cannot distinguish them (V5 §5.1)",
            "replicates 2 and 3 are consumed here and are NOT available to Stage 27 (V5 §12)",
            "PER-REPLICATE ATTRIBUTION IS NOT STABLE. The primary and sensitivity arms differ by "
            "one positive clone, yet they disagree about which replicate carries the larger "
            "delta_AP -- primary favours replicate 2, the author-coefficient arm favours "
            "replicate 3. Both arms agree on direction and on the pooled dual gate. No claim may "
            "be made about which replicate the signal lives in.",
            "the design is underpowered (measured 0.64 against a 0.80 threshold), so the observed "
            "effect size is likely inflated by selection-on-significance and must not be quoted "
            "as an effect estimate",
        ],
        "source_provenance": S232.source_provenance(),
    }
    S232.write_json(VERDICT_JSON, out)
    return out


# =============================================================================================== #
# AUDIT: power measured against the TEST's own null.
#
# 23.2E's power study rejects against a label-randomisation null (positives placed at random, both
# models refit). The TEST rejects against a profile-permutation null with R1 PINNED. Their p95
# thresholds differ by ~3.9x, so the recorded power describes a rejection rule the test never uses.
#
# The exact fix -- a full permutation null nested inside each of the 100 alternative draws -- is
# 20,000 fits. Instead: draw a small number of synthetic y's under the alternative, build the REAL
# permutation null for each (the same construction the test uses, R1 pinned at that y's own
# observed value), and apply the resulting threshold to the alternative draws already cached.
#
# Using several y's rather than one is deliberate: 23.2H found the null to be fold-sensitive, and
# folds move with y. If the per-y thresholds disagree, that disagreement is itself the answer.
# =============================================================================================== #
AUDIT_N_Y = 2
AUDIT_N_PERM_PER_Y = 150
SEED_AUDIT_Y = 23560
AUDIT_JSON = _RESULTS / "stage23_2h_power_audit.json"


def _audit_cache_path() -> Path:
    return _NULL_DIR / "power_audit.partial.jsonl"


def _load_audit_cache() -> dict[int, float]:
    p = _audit_cache_path()
    if not p.exists():
        return {}
    proto = S23.canonical_text_sha256(PROTOCOL_V5)
    done: dict[int, float] = {}
    for line in p.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if "protocol_sha256" in rec:
            if rec["protocol_sha256"] != proto:
                raise RuntimeError("power_audit cache was written under another protocol")
            continue
        done[int(rec["i"])] = float(rec["ap_null"])
    return done


def _audit_synthetic_ys(X, nuis, y, fold):
    """The same generator 23.2C uses, at the same pre-registered alternative."""
    per_fold_pos = {int(f): int(y[fold == f].sum()) for f in np.unique(fold)}
    z = _synthetic_direction(X, nuis)
    cal = _calibrate_beta(z, fold, per_fold_pos, TARGET_AUC)
    ys = [_assign(z, fold, per_fold_pos, cal["beta"], np.random.default_rng(SEED_AUDIT_Y + k))
          for k in range(AUDIT_N_Y)]
    return ys, cal


def run_23_2h_power_audit(shard: int = 0, n_shards: int = 1) -> dict:
    t0 = time.perf_counter()
    X, nuis, y, fold, _rep, strata, _clones = _cohort()
    ys, _cal = _audit_synthetic_ys(X, nuis, y, fold)
    total = AUDIT_N_Y * AUDIT_N_PERM_PER_Y
    done = _load_audit_cache()
    todo = [b for b in range(total) if b % n_shards == shard and b not in done]
    print(f"  [audit shard {shard}/{n_shards}] {len(todo)} of {total} draws to compute", flush=True)

    p = _audit_cache_path()
    if not p.exists():
        p.write_text(json.dumps({"protocol_sha256": S23.canonical_text_sha256(PROTOCOL_V5)}) + "\n",
                     encoding="utf-8")
    cache = {f: S23._frozen_pipeline_cache(X, np.flatnonzero(fold != f), max(S23.K_CANDIDATES))
             for f in range(S23.N_OUTER)}
    for n, b in enumerate(todo, start=1):
        k = b // AUDIT_N_PERM_PER_Y
        rng = np.random.default_rng(SEED_PERMUTATION + 100000 + b)
        ap = S23._rewind_null_once(X, ys[k], fold, nuis, strata, cache, rng)
        with open(p, "a", encoding="utf-8") as fh:
            fh.write(json.dumps({"i": int(b), "ap_null": float(ap)}) + "\n")
        if n % 10 == 0 or n == len(todo):
            el = time.perf_counter() - t0
            print(f"  [audit] {n}/{len(todo)}  {el / 60:.1f} min  {el / n:.1f} s/draw  "
                  f"eta {el / n * (len(todo) - n) / 60:.1f} min", flush=True)
    return {"stage": "23.2H-POWER-AUDIT-shard", "shard": shard, "computed": len(todo),
            "cached_total": len(_load_audit_cache()), "target": total,
            "runtime_minutes": round((time.perf_counter() - t0) / 60, 3)}


def merge_power_audit() -> dict:
    t0 = time.perf_counter()
    X, nuis, y, fold, _rep, strata, _clones = _cohort()
    ys, cal = _audit_synthetic_ys(X, nuis, y, fold)
    done = _load_audit_cache()
    total = AUDIT_N_Y * AUDIT_N_PERM_PER_Y
    missing = [b for b in range(total) if b not in done]
    if missing:
        raise RuntimeError(f"power audit incomplete: {len(missing)} of {total} draws missing")

    per_y = []
    for k in range(AUDIT_N_Y):
        obs_k = _observed(X, nuis, ys[k], fold)
        nulls = np.array([done[k * AUDIT_N_PERM_PER_Y + j]
                          for j in range(AUDIT_N_PERM_PER_Y)]) - obs_k["AP_R1_bdepth_only"]
        per_y.append({"y_index": k, "positives": int(ys[k].sum()),
                      "observed_delta_AP": obs_k["delta_AP"],
                      "null_p95": float(np.percentile(nulls, 95)),
                      "null_mean": float(nulls.mean()), "null_sd": float(nulls.std()),
                      "n_perm": AUDIT_N_PERM_PER_Y})
    thresholds = [d["null_p95"] for d in per_y]

    def _shard_values(name):
        return np.array([json.loads(x)["delta_ap"]
                         for x in _power_shard_path(name).read_text(
                             encoding="utf-8").splitlines()
                         if x.strip() and '"i"' in x])

    alt = _shard_values("alt")
    recorded = json.loads(POWER_JSON.read_text(encoding="utf-8"))
    mean_thr = float(np.mean(thresholds))
    out = {
        "stage": "23.2H-POWER-AUDIT",
        "question": "what is the power of the test AS ACTUALLY RUN -- against the "
                    "profile-permutation null with R1 pinned -- rather than against the "
                    "label-randomisation null the 23.2E instrument builds?",
        "target_AUC": TARGET_AUC,
        "calibration": cal,
        "per_synthetic_y": per_y,
        "threshold_spread": {"min": min(thresholds), "max": max(thresholds),
                             "ratio": max(thresholds) / min(thresholds)},
        "recorded_instrument": {"null_p95": float(np.percentile(_shard_values("null"), 95)),
                                "power": recorded["power"]},
        "corrected_power_per_y": {f"y{d['y_index']}": float((alt > d["null_p95"]).mean())
                                  for d in per_y},
        "corrected_power": float((alt > mean_thr).mean()),
        "mean_threshold": mean_thr,
        "n_alternative_draws": int(len(alt)),
        "threshold": POWER_THRESHOLD,
        "note": "REPORTED, not gating. Gate 18.3 was evaluated and recorded under the frozen V5 "
                "instrument and its FAILED verdict stands regardless of this number. This audit "
                "establishes how far that instrument was off, which bears on how much weight the "
                "Role-A result can carry -- not on whether the gate passed.",
        "runtime_minutes": round((time.perf_counter() - t0) / 60, 3),
        "source_provenance": S232.source_provenance(),
    }
    S232.write_json(AUDIT_JSON, out)
    return out


# =============================================================================================== #
# Smoke test for the sharded permutation engine.
#
# The claim that has to hold before a 9-hour run is launched is not "it does not crash" but
# "sharding, interrupting and resuming reproduce EXACTLY the sequential answer". This checks that
# against a real sequential computation on the real cohort, in a scratch cache, and then measures
# throughput so the ETA is observed rather than guessed.
# =============================================================================================== #
def run_smoke(arm: str = "PRIMARY", n: int = 6) -> dict:
    import shutil

    t0 = time.perf_counter()
    scratch = _CACHE / "h_perm_smoke"
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True)
    X, nuis, y, fold, _rep, strata, _clones = _cohort(arm)
    checks: list[dict] = []

    def check(name, ok, detail=""):
        checks.append({"check": name, "pass": bool(ok), "detail": detail})
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail else ""),
              flush=True)

    # ---- 1. ground truth: sequential ---------------------------------------------------------- #
    seq_dir = scratch / "seq"
    ts = time.perf_counter()
    _compute_null_draws(arm, list(range(n)), X, y, fold, nuis, strata, ts, scratch=seq_dir)
    seq = _load_null_cache(arm, seq_dir)
    per_draw = (time.perf_counter() - ts) / n
    check("sequential run produced every draw", len(seq) == n, f"{len(seq)}/{n}")

    # ---- 2. the same draws, split across 3 shards in a DIFFERENT order ------------------------ #
    sh_dir = scratch / "shard"
    for s in (2, 0, 1):
        idx = [b for b in range(n) if b % 3 == s]
        _compute_null_draws(arm, idx, X, y, fold, nuis, strata, time.perf_counter(),
                            scratch=sh_dir)
    shard = _load_null_cache(arm, sh_dir)
    same = all(shard[b] == seq[b] for b in range(n))
    check("3 shards, run out of order, are BIT-IDENTICAL to sequential", same,
          f"max abs diff {max(abs(shard[b] - seq[b]) for b in range(n)):.3e}")

    # ---- 3. interrupt and resume -------------------------------------------------------------- #
    rs_dir = scratch / "resume"
    _compute_null_draws(arm, [0, 1], X, y, fold, nuis, strata, time.perf_counter(), scratch=rs_dir)
    partial = _load_null_cache(arm, rs_dir)
    todo = [b for b in range(n) if b not in partial]
    _compute_null_draws(arm, todo, X, y, fold, nuis, strata, time.perf_counter(), scratch=rs_dir)
    resumed = _load_null_cache(arm, rs_dir)
    check("interrupt-and-resume reproduces sequential exactly",
          all(resumed[b] == seq[b] for b in range(n)), f"{len(resumed)}/{n} draws")

    # ---- 4. the mixed-protocol guard actually fires -------------------------------------------- #
    bad = scratch / "bad"
    bad.mkdir()
    bp = _null_cache_path(arm, bad)
    bp.write_text(json.dumps({"protocol_sha256": "0" * 64}) + "\n", encoding="utf-8")
    try:
        _load_null_cache(arm, bad)
        check("mixed-protocol cache is refused", False, "guard did NOT fire")
    except RuntimeError as e:
        check("mixed-protocol cache is refused", "refusing a mixed-protocol cache" in str(e))

    # ---- 5. the existing 200-draw result is still reproducible from its own cache -------------- #
    prior = _RESULTS / f"stage23_2h_confirmation{_suffix(arm)}_n200.json"
    check("the frozen 200-permutation result is preserved on disk", prior.exists(), prior.name)

    shutil.rmtree(scratch)
    ok = all(c["pass"] for c in checks)
    eta = {f"n_perm={m}": {
        "one_process_hours": round(m * per_draw / 3600, 2),
        "two_shards_hours": round(m * per_draw / 2 / 3600, 2),
        "three_shards_hours": round(m * per_draw / 3 / 3600, 2)} for m in (500, 1000, 2000)}
    out = {"stage": "23.2H-SMOKE", "arm": arm, "draws_used": n,
           "measured_seconds_per_draw": round(per_draw, 2),
           "note": "measured single-process and unloaded; parallel shards contend, so real "
                   "per-draw cost rises ~10-30% per extra concurrent shard",
           "eta": eta, "checks": checks, "all_passed": ok,
           "already_cached_for_the_real_run": len(_load_null_cache(arm)),
           "runtime_minutes": round((time.perf_counter() - t0) / 60, 3)}
    S232.write_json(_RESULTS / f"stage23_2h_smoke{_suffix(arm)}.json", out)
    return out


# =============================================================================================== #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage 23.2H — Role-A independent confirmation")
    ap.add_argument("--stage", required=True,
                    choices=["23.2h-a", "23.2h-a-authorbug", "23.2h-b", "23.2h-c",
                             "23.2h-d", "23.2h-d-perm", "23.2h-e", "23.2h-smoke",
                             "23.2h-power-audit", "23.2h-power-audit-merge"])
    ap.add_argument("--n-perm", type=int, default=N_PERMUTATION)
    ap.add_argument("--power-kind", default="all", choices=["all", "null", "alt"])
    ap.add_argument("--n-sims", type=int, default=None)
    ap.add_argument("--arm", default="PRIMARY", choices=["PRIMARY", R3_SENSITIVITY_ID])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--target-auc", type=float, default=TARGET_AUC,
                    help="power-curve point; 0.66 is the pre-registered GATE point and the only "
                         "value gate 18.3 may be evaluated at")
    a = ap.parse_args(argv)

    if a.stage == "23.2h-a":
        r = run_23_2h_a(author_bug=False)
        print(json.dumps({k: r[k] for k in ("realized", "totals")}, indent=2))
    elif a.stage == "23.2h-a-authorbug":
        r = run_23_2h_a(author_bug=True)
        print(json.dumps({k: r[k] for k in ("realized", "totals")}, indent=2))
    elif a.stage == "23.2h-b":
        print(json.dumps(run_23_2h_b(), indent=2, default=str))
    elif a.stage == "23.2h-c":
        print(json.dumps(run_23_2h_c(a.power_kind, a.n_sims, a.target_auc),
                         indent=2, default=str))
    elif a.stage == "23.2h-d":
        print(json.dumps(run_23_2h_d(a.n_perm, a.arm), indent=2, default=str))
    elif a.stage == "23.2h-d-perm":
        print(json.dumps(run_23_2h_d_perm(a.arm, a.shard, a.n_shards, a.n_perm), indent=2))
    elif a.stage == "23.2h-power-audit":
        print(json.dumps(run_23_2h_power_audit(a.shard, a.n_shards), indent=2))
    elif a.stage == "23.2h-power-audit-merge":
        print(json.dumps(merge_power_audit(), indent=2, default=str))
    elif a.stage == "23.2h-smoke":
        r = run_smoke(a.arm)
        print(json.dumps(r, indent=2))
        return 0 if r["all_passed"] else 1
    else:
        print(json.dumps(run_23_2h_e(), indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
