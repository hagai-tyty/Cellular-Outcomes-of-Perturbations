"""Tests for the real-clock loader/fitter and the reprogramming connector —
the pieces that let CellFate-Rx consume real aging clocks and OSKM data."""

from __future__ import annotations

import numpy as np
import pytest

from cellfate.common.errors import ConfigError
from cellfate.common.panel import GenePanel
from cellfate.data.aging import LinearClock
from cellfate.data.build_dataset import DataConfig, build_clock
from cellfate.data.clock_fit import fit_linear_clock
from cellfate.data.sources import ReprogrammingSource

GENES = [f"G{i}" for i in range(40)]


def _panel() -> GenePanel:
    return GenePanel(GENES)


# --------------------------------------------------------------------------- #
# Clock loader: real weights load; missing/unknown fails loud (no silent random)
# --------------------------------------------------------------------------- #
def test_linear_clock_json_roundtrip(tmp_path):
    clock = LinearClock({"G0": 1.5, "G1": -2.0}, intercept=40.0)
    p = tmp_path / "clock.json"
    clock.to_json(p, meta={"source": "unit"})
    back = LinearClock.from_json(p)
    assert back.weights == clock.weights and back.intercept == clock.intercept


def test_from_json_rejects_empty_weights(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text('{"weights": {}, "intercept": 0.0}')
    with pytest.raises(ValueError):
        LinearClock.from_json(p)


def test_build_clock_random_is_explicit():
    clock = build_clock(DataConfig(out="x", gene_panel="x", clock="random", seed=0), _panel())
    assert isinstance(clock, LinearClock)


def test_build_clock_fails_loud_on_missing_weights():
    # a named clock that isn't 'random' and has no file must NOT silently go random
    cfg = DataConfig(out="x", gene_panel="x", clock="buckley", seed=0)
    with pytest.raises(ConfigError):
        build_clock(cfg, _panel())


def test_build_clock_loads_fitted_file(tmp_path):
    LinearClock({"G0": 1.0}, intercept=40.0).to_json(tmp_path / "c.json")
    cfg = DataConfig(out="x", gene_panel="x", clock=str(tmp_path / "c.json"), seed=0)
    clock = build_clock(cfg, _panel())
    assert clock.weights.get("G0") == 1.0


# --------------------------------------------------------------------------- #
# Clock fitter: recovers an age signal and produces a loadable clock
# --------------------------------------------------------------------------- #
def test_fit_linear_clock_recovers_age_signal():
    rng = np.random.default_rng(0)
    n, g = 80, 40
    ages = rng.uniform(1, 96, n)
    counts = rng.poisson(30, (n, g)).astype(float)
    counts[:, :8] += ages[:, None] * 2.0        # first 8 genes scale with age
    clock, metrics = fit_linear_clock(counts, GENES, ages, seed=0)
    assert metrics["cv_pearson"] > 0.5          # learned a real age signal
    assert set(clock.weights) <= set(GENES) and np.isfinite(clock.intercept)


# --------------------------------------------------------------------------- #
# Reprogramming connector: valid RawChunk, age-valid source, OSKM token encoding
# --------------------------------------------------------------------------- #
def test_reprogramming_build_chunk_is_valid_and_age_valid():
    from cellfate.common.constants import CANCER_SOURCES
    assert "reprogramming" not in CANCER_SOURCES   # age is NOT masked (unlike tahoe)

    counts = np.random.default_rng(0).poisson(20, (10, 40)).astype(np.float32)
    pert = ["control"] * 5 + ["OSKM"] * 5
    time_h = [0.0] * 5 + [312.0] * 5
    rc = ReprogrammingSource.build_chunk("reprogramming:FIB", counts, GENES, "FIB",
                                         pert, time_h)
    assert rc.source == "reprogramming" and len(rc.obs) == 10
    assert rc.obs["is_control"].sum() == 5
    # control carries no SMILES token; OSKM carries the factor-set token
    assert rc.obs.loc[rc.obs.is_control, "smiles"].eq("").all()
    assert (rc.obs.loc[~rc.obs.is_control, "smiles"] == "OSKM").all()
    assert sorted(set(rc.obs["time_h"])) == [0.0, 312.0]


def test_reprogramming_factor_tokens_give_distinct_features():
    from cellfate.data.perturbation import encode_fingerprints
    fp_ctrl, fp_oskm, fp_osk = encode_fingerprints(["", "OSKM", "OSK"])
    assert not np.array_equal(fp_ctrl, fp_oskm)
    assert not np.array_equal(fp_oskm, fp_osk)


# --------------------------------------------------------------------------- #
# The clock consumes the FULL profile, not the HVG panel (decoupling).         #
# --------------------------------------------------------------------------- #
def test_clock_consumes_full_profile_not_panel():
    import pandas as pd

    from cellfate.data.aging import delta_age
    # a clock weight on a gene must drive ΔAge even though delta_age is handed the
    # full gene list (no panel) -- guards the fix that feeds the clock `norm`, not x_panel.
    clock = LinearClock({"AGING_GENE": 5.0}, intercept=40.0)
    genes = ["AGING_GENE", "OTHER1", "OTHER2"]
    expr = np.array([[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]])   # cell0 expresses it, cell1 (control) doesn't
    obs = pd.DataFrame({"cell_line": ["L", "L"], "is_control": [False, True]})
    d, mask = delta_age(clock, expr, genes, obs, source="reprogramming")
    assert mask.all()                       # reprogramming is age-valid
    assert abs(d[0] - 10.0) < 1e-9          # 40+5*2=50 vs control 40 -> ΔAge=10 (clock read the gene)


# --------------------------------------------------------------------------- #
# Gill 2022 (GSE165176) connector: Log2-RPM -> counts, day-0 baseline logic.   #
# --------------------------------------------------------------------------- #
def test_gill_connector_parses_log2rpm_and_baseline(tmp_path):
    import gzip

    from cellfate.data.sources import GillReprogrammingSource

    # tiny fake GEO matrix: 12 annotation cols + 3 sample cols, Log2 RPM values
    expr = tmp_path / "expr.tsv.gz"
    header = ["Probe", "Chromosome", "Start", "End", "Probe Strand", "Feature", "ID",
              "Description", "Feature Strand", "Type", "Feature Orientation", "Distance",
              "D1_Fib_Sendai", "D1_d13_SSEA4_Sendai", "D1_d54_SSEA4_Sendai"]
    ann = ["1", "1", "2", "+", "g", "e", "d", "+", "gene", "N", "0"]  # 11 filler annotation fields
    rows = [
        ["COL1A1", *ann, "3.0", "1.0", "0.0"],
        ["FN1",    *ann, "2.0", "2.0", "1.0"],
        ["VIM",    *ann, "1.0", "3.0", "2.0"],
    ]
    with gzip.open(expr, "wt") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    # tiny series matrix: titles + characteristics (day + cell type)
    sm = tmp_path / "series.txt.gz"
    with gzip.open(sm, "wt") as f:
        f.write('!Sample_title\t"D1_Fib_Sendai"\t"D1_d13_SSEA4_Sendai"\t"D1_d54_SSEA4_Sendai"\n')
        f.write('!Sample_characteristics_ch1\t"cell type: Dermal fibroblast"\t"cell type: Reprogramming fibroblast"\t"cell type: iPSC"\n')
        f.write('!Sample_characteristics_ch1\t"days of reprogramming: 0"\t"days of reprogramming: 13"\t"days of reprogramming: 54"\n')

    src = GillReprogrammingSource(str(expr), str(sm))
    chunks = src.plan()
    assert [c["cell_line"] for c in chunks] == ["D1"]
    rc = src.fetch(chunks[0])
    assert rc.source == "reprogramming" and rc.counts.shape == (3, 3)  # 3 samples x 3 genes
    # Log2 RPM -> linear RPM (2^x - 1): COL1A1 sample0 = 2^3 - 1 = 7
    assert abs(float(rc.counts[0, list(rc.genes).index("COL1A1")]) - 7.0) < 1e-4
    # day-0 dermal fibroblast is the control; d13/d54 are OSKM
    assert rc.obs["is_control"].tolist() == [True, False, False]
    assert rc.obs["pert_id"].tolist() == ["control", "OSKM", "OSKM"]
    assert rc.obs["time_h"].tolist() == [0.0, 13 * 24.0, 54 * 24.0]


# --------------------------------------------------------------------------- #
# GSE242423 (Kundaje) human single-cell 10x connector: droplet filter, symbol  #
# dedup, D0-control baseline.                                                   #
# --------------------------------------------------------------------------- #
def test_gse242423_connector_drops_empty_droplets_and_maps_symbols(tmp_path):
    import gzip

    from cellfate.data.sources import GSE242423SingleCellSource

    # tiny 10x reference: 4 Ensembl rows, GENEA duplicated (tests highest-expressed dedup)
    genes = tmp_path / "genes.tsv.gz"
    with gzip.open(genes, "wt") as f:
        for ens, sym in [("ENSG1", "GENEA"), ("ENSG2", "GENEB"),
                         ("ENSG3", "GENEA"), ("ENSG4", "MT-CO1")]:
            f.write(f"{ens}\t{sym}\tGene Expression\n")

    def write_mtx(path, entries, n_genes, n_bc):
        with gzip.open(path, "wt") as f:
            f.write("%%MatrixMarket matrix coordinate integer general\n%\n")
            f.write(f"{n_genes} {n_bc} {len(entries)}\n")
            for g, b, v in entries:
                f.write(f"{g} {b} {v}\n")

    # 3 barcodes: bc1 (3 genes), bc2 (2 genes) are real; bc3 (1 gene) is an empty droplet
    mtx = tmp_path / "D0.matrix.mtx.gz"
    write_mtx(mtx, [(1, 1, 5), (2, 1, 3), (3, 1, 2), (1, 2, 4), (2, 2, 1), (1, 3, 9)], 4, 3)
    bc = tmp_path / "D0.barcodes.tsv.gz"
    with gzip.open(bc, "wt") as f:
        f.write("BC1-1\nBC2-1\nBC3-1\n")

    src = GSE242423SingleCellSource(
        [{"matrix": str(mtx), "barcodes": str(bc), "label": "D0"}],
        str(genes), cell_line="HFF", min_genes=2)     # min_genes=2 drops bc3
    rc = src.fetch(src.plan()[0])
    assert rc.source == "reprogramming"
    assert rc.counts.shape[0] == 2                      # empty droplet (bc3) dropped
    assert set(rc.genes) == {"GENEA", "GENEB", "MT-CO1"}  # duplicate ENSG collapsed to one GENEA
    assert rc.obs["is_control"].all()                   # D0 -> control baseline
    assert (rc.obs["pert_id"] == "control").all()
    # GENEA kept the highest-expressed Ensembl row (ENSG1 total 9 > ENSG3 total 2)
    gi = list(rc.genes).index("GENEA")
    assert rc.counts[0, gi] == 5.0                      # bc1's ENSG1 value


def test_gse242423_infers_day_and_oskm_from_label():
    from cellfate.data.sources import GSE242423SingleCellSource
    src = GSE242423SingleCellSource([], "x", ipsc_day=21.0)
    assert src._day_of({"label": "D0"}) == (0.0, True)
    assert src._day_of({"label": "D8"}) == (8.0, False)
    assert src._day_of({"label": "iPSC"}) == (21.0, False)
