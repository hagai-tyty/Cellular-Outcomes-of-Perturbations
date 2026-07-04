"""Data sources (Document 2, S3-S4).

A *source* knows how to (1) ``plan`` its work as a list of independently
processable :class:`~cellfate.data.chunking.CellChunk` units and (2) ``fetch``
the raw cells for one chunk as a :class:`RawChunk`. The rest of the ETL is
source-agnostic: it only ever sees ``RawChunk`` objects.

``SyntheticSource`` has **no heavy dependencies** and generates structured fake
data, so the entire pipeline can run end-to-end (in tests / CI) without scanpy,
rdkit, Hugging Face, or any download. The real connectors (Tahoe-100M, sci-Plex)
lazy-import their dependencies and are documented but not exercised here.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
import pandas as pd

from cellfate.common import constants as C
from cellfate.common.errors import DataSourceError

from .chunking import CellChunk

# Columns every source must provide in RawChunk.obs.
OBS_COLUMNS = ("cell_id", "cell_line", "pert_id", "smiles", "scaffold_id",
               "dose_uM", "time_h", "is_control")


@dataclass
class RawChunk:
    """Raw, pre-QC cells for one chunk.

    ``counts`` is a dense (N, G_full) non-negative matrix; ``genes`` are the
    matching gene symbols (length G_full); ``obs`` carries per-cell metadata
    (see :data:`OBS_COLUMNS`). Real connectors may hand back sparse matrices;
    convert to dense per chunk (chunks are sized to fit in memory).
    """

    chunk_id: str
    source: str
    counts: np.ndarray
    genes: list[str]
    obs: pd.DataFrame

    def __post_init__(self) -> None:
        n, g = self.counts.shape
        if g != len(self.genes):
            raise DataSourceError(
                f"{self.chunk_id}: counts has {g} genes but {len(self.genes)} names"
            )
        if len(self.obs) != n:
            raise DataSourceError(
                f"{self.chunk_id}: counts has {n} cells but obs has {len(self.obs)}"
            )
        missing = [c for c in OBS_COLUMNS if c not in self.obs.columns]
        if missing:
            raise DataSourceError(f"{self.chunk_id}: obs missing columns {missing}")


class DataSource(ABC):
    """Abstract source. Subclasses implement :meth:`plan` and :meth:`fetch`."""

    name: str

    @abstractmethod
    def plan(self) -> list[CellChunk]:
        """Return the list of chunks this source will produce."""

    @abstractmethod
    def fetch(self, chunk: CellChunk) -> RawChunk:
        """Load and return the raw cells for one planned chunk."""


# --------------------------------------------------------------------------- #
# Synthetic source (no external dependencies; used for tests + CI)            #
# --------------------------------------------------------------------------- #
class SyntheticSource(DataSource):
    """Deterministic synthetic perturbation data with planted biological signal.

    Each chunk is one cell line containing vehicle controls plus several
    compounds. Every compound is assigned an *effect* (``safe`` / ``loss`` /
    ``death``) that up-regulates the matching signature genes, so soft labels are
    non-trivial. A per-cell proliferation level drives the cell-cycle genes and
    is deliberately correlated with expression, so the deconfounder has a real
    artefact to remove. Compounds are grouped into scaffold families so the
    leave-scaffold-out split is meaningful.
    """

    def __init__(
        self,
        name: str = "synth",
        *,
        n_lines: int = 3,
        n_compounds: int = 6,
        n_cells_per_condition: int = 25,
        n_filler_genes: int = 150,
        n_scaffold_families: int = 3,
        doses: tuple[float, ...] = (0.1, 1.0),
        is_cancer: bool = False,
        seed: int = 0,
    ) -> None:
        self.name = name
        self.n_lines = n_lines
        self.n_compounds = n_compounds
        self.n_cells = n_cells_per_condition
        self.n_scaffold_families = n_scaffold_families
        self.doses = doses
        self.is_cancer = is_cancer
        self.seed = seed

        # Gene panel layout: signature genes + cell-cycle genes + a couple of
        # mito genes (for QC) + filler genes. De-duplicated, stable order.
        genes: list[str] = []
        for gene_set in C.DEFAULT_SIGNATURES.values():
            genes += list(gene_set)
        genes += list(C.S_GENES) + list(C.G2M_GENES)
        genes += ["MT-CO1", "MT-ND1"]
        genes += [f"FILL{i}" for i in range(n_filler_genes)]
        seen: set[str] = set()
        self.genes: list[str] = [g for g in genes if not (g in seen or seen.add(g))]
        self._gidx = {g: i for i, g in enumerate(self.genes)}

        # Compound table (shared across lines): effect + scaffold + smiles.
        effects = ("loss", "death", "safe")
        self.compounds = [
            {
                "pert_id": f"{name}_c{c}",
                "smiles": f"SYN-{name}-c{c}",
                "effect": effects[c % len(effects)],
                "scaffold_id": f"{name}_SCAF{c % n_scaffold_families}",
            }
            for c in range(n_compounds)
        ]

    # -- planning ----------------------------------------------------------- #
    def plan(self) -> list[CellChunk]:
        return [
            CellChunk(
                id=f"{self.name}:line{i}",
                uri="synthetic",
                source=self.name,
                cell_line=f"{self.name.upper()}_L{i}",
                pert_ids=[c["pert_id"] for c in self.compounds],
            )
            for i in range(self.n_lines)
        ]

    # -- fetching ----------------------------------------------------------- #
    def fetch(self, chunk: CellChunk) -> RawChunk:
        rng = np.random.default_rng(self.seed + hash(chunk["id"]) % 10_000)
        g = len(self.genes)
        base_log = rng.normal(2.0, 0.4, size=g)  # per-gene baseline log-rate

        rows: list[np.ndarray] = []
        obs_records: list[dict] = []

        def emit(pert: dict | None, dose: float, time_h: float, k: int) -> None:
            is_ctrl = pert is None
            effect = "safe" if is_ctrl else pert["effect"]
            log_rate = base_log.copy()
            # effect up-regulates the matching signature genes
            for gene in C.DEFAULT_SIGNATURES[effect]:
                if gene in self._gidx:
                    log_rate[self._gidx[gene]] += 1.6
            # proliferation level drives cell-cycle genes (and thus the clock)
            cc_level = float(rng.uniform(0.0, 1.0))
            for gene in list(C.S_GENES) + list(C.G2M_GENES):
                if gene in self._gidx:
                    log_rate[self._gidx[gene]] += 1.5 * cc_level
            lib = float(rng.uniform(0.8, 1.3))
            rate = np.exp(log_rate) * lib
            rows.append(rng.poisson(rate).astype(np.float32))
            obs_records.append({
                "cell_id": f"{self.name}:{chunk['cell_line']}:{k}",
                "cell_line": chunk["cell_line"],
                "pert_id": "control" if is_ctrl else pert["pert_id"],
                "smiles": "" if is_ctrl else pert["smiles"],
                "scaffold_id": "CONTROL" if is_ctrl else pert["scaffold_id"],
                "dose_uM": 0.0 if is_ctrl else dose,
                "time_h": time_h,
                "is_control": is_ctrl,
            })

        k = 0
        for _ in range(self.n_cells * 2):  # controls
            emit(None, 0.0, 24.0, k)
            k += 1
        for comp in self.compounds:
            for dose in self.doses:
                for _ in range(self.n_cells):
                    emit(comp, dose, 24.0, k)
                    k += 1

        counts = np.vstack(rows)
        obs = pd.DataFrame.from_records(obs_records)
        return RawChunk(chunk_id=chunk["id"], source=self.name,
                        counts=counts, genes=list(self.genes), obs=obs)


# --------------------------------------------------------------------------- #
# Real connectors (lazy deps; documented, not exercised in the foundation)    #
# --------------------------------------------------------------------------- #
class TahoeSource(DataSource):
    """Tahoe-100M via the Hugging Face ``datasets`` streaming API (Arc Atlas).

    Plans one chunk per (cell line x drug). Cancer cell lines, so the aging
    label is masked downstream (``tahoe`` is in ``CANCER_SOURCES``).
    """

    name = "tahoe"

    def __init__(self, hf_path: str = "tahoebio/Tahoe-100M", cell_lines=None, drugs=None):
        self.hf_path = hf_path
        self.cell_lines = cell_lines
        self.drugs = drugs

    def plan(self) -> list[CellChunk]:  # pragma: no cover - needs network
        raise DataSourceError(
            "TahoeSource.plan requires the 'data' extra (datasets) and network "
            "access to Hugging Face. Configure cell_lines/drugs and implement the "
            "census query per your deployment."
        )

    def fetch(self, chunk: CellChunk) -> RawChunk:  # pragma: no cover - needs network
        from datasets import load_dataset  # noqa: F401  (lazy import)

        raise DataSourceError("TahoeSource.fetch not wired in the foundation build.")


class SciplexSource(DataSource):
    """sci-Plex (GEO) connector skeleton; reads an .h5ad via anndata (lazy)."""

    name = "sciplex"

    def __init__(self, h5ad_path: str | None = None):
        self.h5ad_path = h5ad_path

    def plan(self) -> list[CellChunk]:  # pragma: no cover - needs data file
        raise DataSourceError("SciplexSource requires a local .h5ad path.")

    def fetch(self, chunk: CellChunk) -> RawChunk:  # pragma: no cover - needs data file
        import anndata  # noqa: F401  (lazy import)

        raise DataSourceError("SciplexSource.fetch not wired in the foundation build.")


class ReprogrammingSource(DataSource):
    """Partial / transient reprogramming connector (OSKM/OSK time-courses).

    This is the **age-valid** arm of the project: unlike Tahoe (cancer, age-masked),
    reprogramming cells carry a real biological-age signal, so ``source`` is kept OUT
    of ``CANCER_SOURCES`` and the ΔAge head trains on them. Targets datasets like
    Gill 2022 (GSE165180) or Roux 2022 (GSE197437) after conversion to .h5ad.

    Perturbation encoding: OSKM is a **transcription-factor cocktail, not a molecule**,
    so there is no SMILES. The factor-set identity is carried as a stable token in the
    ``smiles`` field (``factor_as_token=True``) so the existing hashed-fingerprint path
    yields a distinct, deterministic perturbation feature per cocktail (OSKM vs OSK vs
    control). The reprogramming **timepoint** is the key axis and maps to ``time_h``;
    factor **dose** (if the dataset varies it, e.g. dox level) maps to ``dose_uM``.
    A dedicated TF encoder (schema already reserves ``u_tf_emb``) is the v2 upgrade.

    ``plan``/``fetch`` read a local AnnData; supply a column map for your dataset. The
    RawChunk assembly is exposed as :meth:`build_chunk` so it is unit-testable without
    a file.
    """

    name = "reprogramming"

    def __init__(
        self,
        h5ad_path: str | None = None,
        *,
        cell_line_col: str = "cell_line",
        factor_col: str = "factors",
        timepoint_h_col: str = "time_h",
        dose_col: str | None = None,
        control_value: str = "control",
        factor_as_token: bool = True,
    ) -> None:
        self.h5ad_path = h5ad_path
        self.cell_line_col = cell_line_col
        self.factor_col = factor_col
        self.timepoint_h_col = timepoint_h_col
        self.dose_col = dose_col
        self.control_value = control_value
        self.factor_as_token = factor_as_token

    @staticmethod
    def build_chunk(
        chunk_id: str,
        counts: np.ndarray,
        genes: list[str],
        cell_line: str,
        pert_ids: list[str],
        time_h: list[float],
        *,
        doses: list[float] | None = None,
        control_value: str = "control",
        factor_as_token: bool = True,
        source: str = "reprogramming",
    ) -> RawChunk:
        """Assemble a valid RawChunk from per-cell reprogramming metadata.

        ``pert_ids[i]`` is the factor set (e.g. ``"OSKM"``) or ``control_value`` for
        unreprogrammed cells. Encodes the TF cocktail per the class docstring.
        """
        n = counts.shape[0]
        if not (len(pert_ids) == len(time_h) == n):
            raise DataSourceError(f"{chunk_id}: pert_ids/time_h must match {n} cells")
        doses = doses if doses is not None else [1.0] * n
        recs: list[dict] = []
        for i in range(n):
            is_ctrl = pert_ids[i] == control_value
            factor = pert_ids[i]
            # TF cocktail: no SMILES; carry the factor token so the hashed-fingerprint
            # path gives a distinct, deterministic feature per cocktail.
            token = "" if (is_ctrl or not factor_as_token) else factor
            recs.append({
                "cell_id": f"{source}:{cell_line}:{i}",
                "cell_line": cell_line,
                "pert_id": "control" if is_ctrl else factor,
                "smiles": token,
                "scaffold_id": "CONTROL" if is_ctrl else factor,
                "dose_uM": 0.0 if is_ctrl else float(doses[i]),
                "time_h": float(time_h[i]),
                "is_control": is_ctrl,
            })
        return RawChunk(chunk_id=chunk_id, source=source,
                        counts=np.asarray(counts), genes=list(genes),
                        obs=pd.DataFrame.from_records(recs))

    def _adata(self):  # pragma: no cover - needs data file
        if not self.h5ad_path:
            raise DataSourceError("ReprogrammingSource requires a local .h5ad path.")
        import anndata
        return anndata.read_h5ad(self.h5ad_path, backed="r")

    def plan(self) -> list[CellChunk]:  # pragma: no cover - needs data file
        adata = self._adata()
        chunks: list[CellChunk] = []
        for cl in adata.obs[self.cell_line_col].astype(str).unique():
            chunks.append(CellChunk(id=f"{self.name}:{cl}", cell_line=str(cl)))
        return chunks

    def fetch(self, chunk: CellChunk) -> RawChunk:  # pragma: no cover - needs data file
        adata = self._adata()
        sub = adata[adata.obs[self.cell_line_col].astype(str) == chunk["cell_line"]].to_memory()
        counts = np.asarray(sub.X.todense() if hasattr(sub.X, "todense") else sub.X, dtype=np.float32)
        genes = [str(g) for g in sub.var_names]
        pert_ids = sub.obs[self.factor_col].astype(str).tolist()
        time_h = sub.obs[self.timepoint_h_col].astype(float).tolist()
        doses = (sub.obs[self.dose_col].astype(float).tolist() if self.dose_col else None)
        return self.build_chunk(chunk["id"], counts, genes, chunk["cell_line"],
                                pert_ids, time_h, doses=doses,
                                control_value=self.control_value,
                                factor_as_token=self.factor_as_token)


class GillReprogrammingSource(ReprogrammingSource):
    """Gill et al. 2022 (GSE165176) -- bulk transient-reprogramming RNA-seq.

    Fibroblasts undergoing Sendai (OSKM) iPSC reprogramming, sampled across a
    day-0 -> day-54 time course in 6 donors (neonatal/young/old). This is the
    **age-valid** rejuvenation signal the project needs.

    Format handling:
      * The matrix is **Log2 RPM** (already normalised). We invert it to linear
        RPM (``2**x - 1``) and hand that back as ``counts``; because CP10k(RPM)
        == CP10k(counts), the pipeline's log1p-CP10k normalisation lands the data
        in the exact space the aging clock was fit on.
      * Genes are HGNC symbols (the ``Probe`` column) -- matches the clock keys.
      * Baseline (``is_control``) = the **day-0 "Dermal fibroblast"** of each
        donor, so ΔAge is rejuvenation relative to the un-reprogrammed start.
      * One chunk per donor (= cell_line), so the group-aware split leaves whole
        donors out. Bulk: each sample is one row.
    """

    name = "reprogramming"
    _ANNOT_COLS = 12  # Probe..Distance precede the sample columns in the GEO matrix

    def __init__(self, expr_tsv: str, series_matrix: str) -> None:
        super().__init__(factor_as_token=True)
        self.expr_tsv = expr_tsv
        self.series_matrix = series_matrix
        self._genes: list[str] | None = None
        self._rpm = None            # linear RPM, genes x samples (DataFrame)
        self._meta: dict | None = None

    def _parse_series(self) -> dict:
        import gzip
        rows: dict[str, list[str]] = {}
        opener = gzip.open if str(self.series_matrix).endswith(".gz") else open
        with opener(self.series_matrix, "rt") as f:
            for line in f:
                if line.startswith("!Sample_title"):
                    rows["title"] = [x.strip('"') for x in line.rstrip().split("\t")[1:]]
                elif line.startswith("!Sample_characteristics_ch1"):
                    vals = [x.strip('"') for x in line.rstrip().split("\t")[1:]]
                    key = vals[0].split(":")[0].strip()
                    rows.setdefault(f"char::{key}", vals)
        titles = rows["title"]
        def field(v: str) -> str:
            return v.split(":", 1)[1].strip() if ":" in v else v
        days = [field(x) for x in rows.get("char::days of reprogramming", ["0"] * len(titles))]
        ctype = [field(x) for x in rows.get("char::cell type", [""] * len(titles))]
        donor = [t.split("_")[0] for t in titles]
        return {t: {"donor": donor[i], "day": float(days[i]), "ctype": ctype[i]}
                for i, t in enumerate(titles)}

    def _load(self) -> None:
        if self._rpm is not None:
            return
        self._meta = self._parse_series()
        df = pd.read_csv(self.expr_tsv, sep="\t")
        gene_col = df.columns[0]                       # "Probe" == HGNC symbol
        sample_cols = list(df.columns[self._ANNOT_COLS:])
        log2 = df[sample_cols].to_numpy(dtype=np.float64)
        rpm = np.power(2.0, log2) - 1.0                # Log2 RPM -> linear RPM
        rpm[rpm < 0] = 0.0
        out = pd.DataFrame(rpm, columns=sample_cols)
        out["__sym__"] = df[gene_col].astype(str).to_numpy()
        out["__tot__"] = rpm.sum(axis=1)
        out = out.sort_values("__tot__", ascending=False).drop_duplicates("__sym__", keep="first")
        self._genes = out["__sym__"].tolist()
        self._rpm = out.set_index("__sym__")[sample_cols]

    def plan(self) -> list[CellChunk]:
        self._load()
        donors = sorted({self._meta[c]["donor"] for c in self._rpm.columns if c in self._meta})
        return [CellChunk(id=f"{self.name}:{d}", cell_line=str(d)) for d in donors]

    def fetch(self, chunk: CellChunk) -> RawChunk:
        self._load()
        donor = chunk["cell_line"]
        cols = [c for c in self._rpm.columns if c in self._meta and self._meta[c]["donor"] == donor]
        if not cols:
            raise DataSourceError(f"no samples for donor {donor}")
        counts = self._rpm[cols].to_numpy(dtype=np.float32).T          # samples x genes
        pert, time_h = [], []
        for c in cols:
            m = self._meta[c]
            is_ctrl = m["day"] == 0.0 or m["ctype"] == "Dermal fibroblast"
            pert.append("control" if is_ctrl else "OSKM")
            time_h.append(m["day"] * 24.0)
        return self.build_chunk(chunk["id"], counts, list(self._genes), donor,
                                pert, time_h, factor_as_token=True)


class GSE242423SingleCellSource(ReprogrammingSource):
    """GSE242423 (Kundaje lab) -- human fibroblast OSKM reprogramming, 10x scRNA-seq.

    A densely-sampled reprogramming time course (D0 -> D14 -> iPSC) at single-cell
    resolution -- the **volume** the fate/safety head needs. Handling:

      * Raw **Cell Ranger MTX** (every droplet). We drop **empty droplets** at load
        with a ``min_genes`` gate (you cannot densify ~2.5M barcodes); real cells
        survive (~10-20k per timepoint). Doublets are NOT removed here -- documented
        caveat (Option C); a Scrublet pass is the future upgrade.
      * Genes: Ensembl in col 1, **HGNC symbol in col 2** (10x features) -> use the
        symbol (matches the human clock, no ortholog step); duplicate symbols keep the
        highest-expressed row.
      * All cells are the **same fibroblast line** -> **one chunk** (all timepoints),
        so the D0 controls anchor the control-relative ΔAge + fate labels. **D0** is
        the ``control`` baseline; every later day is ``OSKM`` (``time_h = day*24``).
      * ``max_cells_per_sample`` caps cells per timepoint (memory); ``None`` = all.

    ``samples`` is a list of dicts: ``{"matrix": path, "barcodes": path, "label": "D8",
    "day": 8.0, "is_control": False}`` (day/is_control inferred from the label if absent).
    """

    name = "reprogramming"

    def __init__(
        self,
        samples: list[dict],
        genes_file: str,
        *,
        cell_line: str = "HFF",
        min_genes: int = 500,
        max_cells_per_sample: int | None = None,
        ipsc_day: float = 21.0,
        seed: int = 0,
    ) -> None:
        super().__init__(factor_as_token=True)
        self.samples = samples
        self.genes_file = genes_file
        self.cell_line = cell_line
        self.min_genes = min_genes
        self.max_cells_per_sample = max_cells_per_sample
        self.ipsc_day = ipsc_day
        self.seed = seed

    def _symbols(self) -> list[str]:
        import gzip
        opener = gzip.open if str(self.genes_file).endswith(".gz") else open
        syms: list[str] = []
        with opener(self.genes_file, "rt") as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                syms.append(p[1] if len(p) > 1 else p[0])   # col 2 = HGNC symbol
        return syms

    def _day_of(self, s: dict) -> tuple[float, bool]:
        if "day" in s and "is_control" in s:
            return float(s["day"]), bool(s["is_control"])
        label = str(s.get("label", "")).upper()
        if label.startswith("IPSC"):
            return self.ipsc_day, False
        digits = "".join(ch for ch in label if ch.isdigit())
        day = float(digits) if digits else 0.0
        return day, (day == 0.0)

    def _header_rows(self, path: str) -> tuple[int, int, int, int]:
        """Return (n_skip, n_genes, n_barcodes, nnz) for a MatrixMarket file."""
        import gzip
        opener = gzip.open if str(path).endswith(".gz") else open
        n_skip = 0
        with opener(path, "rt") as f:
            for line in f:
                n_skip += 1
                if not line.startswith("%"):
                    n_genes, n_bc, nnz = (int(x) for x in line.split())
                    return n_skip, n_genes, n_bc, nnz
        raise DataSourceError(f"{path}: no MatrixMarket dimension line")

    def _load_sample_filtered(self, matrix_path: str, rng):
        """Stream a raw 10x MTX in chunks, drop empty droplets (min_genes), and
        subsample -- never holding all entries in memory. Returns genes x kept_cells."""
        import pandas as pd
        from scipy.sparse import csc_matrix

        n_skip, n_genes, n_bc, _ = self._header_rows(matrix_path)
        names = ["g", "b", "v"]
        dt = {"g": np.int32, "b": np.int32, "v": np.float32}
        # pass 1: genes-per-barcode
        per_bc = np.zeros(n_bc + 1, dtype=np.int32)
        for ch in pd.read_csv(matrix_path, sep=r"\s+", skiprows=n_skip, header=None,
                              names=names, dtype=dt, chunksize=5_000_000):
            np.add.at(per_bc, ch["b"].to_numpy(), 1)
        keep = np.where(per_bc >= self.min_genes)[0]          # 1-indexed barcode ids
        if self.max_cells_per_sample and len(keep) > self.max_cells_per_sample:
            keep = np.sort(rng.choice(keep, size=self.max_cells_per_sample, replace=False))
        newcol = np.full(n_bc + 1, -1, dtype=np.int64)
        newcol[keep] = np.arange(len(keep))
        # pass 2: collect only kept-barcode entries
        gs, cs, vs = [], [], []
        for ch in pd.read_csv(matrix_path, sep=r"\s+", skiprows=n_skip, header=None,
                              names=names, dtype=dt, chunksize=5_000_000):
            b = ch["b"].to_numpy()
            nc = newcol[b]
            m = nc >= 0
            if m.any():
                gs.append(ch["g"].to_numpy()[m] - 1)
                cs.append(nc[m])
                vs.append(ch["v"].to_numpy()[m])
        rows = np.concatenate(gs) if gs else np.array([], np.int32)
        cols = np.concatenate(cs) if cs else np.array([], np.int64)
        vals = np.concatenate(vs) if vs else np.array([], np.float32)
        M = csc_matrix((vals, (rows, cols)), shape=(n_genes, len(keep)), dtype=np.float32)
        return M, len(keep)

    def plan(self) -> list[CellChunk]:
        # one chunk: the whole fibroblast line, so D0 controls anchor the baseline
        return [CellChunk(id=f"{self.name}:{self.cell_line}", cell_line=self.cell_line)]

    def fetch(self, chunk: CellChunk) -> RawChunk:
        from scipy.sparse import hstack

        rng = np.random.default_rng(self.seed)
        blocks, pert, time_h = [], [], []
        for s in self.samples:
            M, n = self._load_sample_filtered(s["matrix"], rng)   # genes x kept cells
            blocks.append(M)
            day, is_ctrl = self._day_of(s)
            pert.extend(["control" if is_ctrl else "OSKM"] * n)
            time_h.extend([day * 24.0] * n)

        genes_all = self._symbols()
        combined = hstack(blocks).tocsr()                   # genes x all_kept_cells
        if combined.shape[0] != len(genes_all):
            raise DataSourceError(f"matrix has {combined.shape[0]} genes but genes.tsv has {len(genes_all)}")

        # dedup duplicate symbols: keep the highest-expressed gene row per symbol
        totals = np.asarray(combined.sum(axis=1)).ravel()
        seen: set[str] = set()
        keep_rows: list[int] = []
        for i in np.argsort(-totals):
            g = genes_all[i]
            if g not in seen:
                seen.add(g)
                keep_rows.append(int(i))
        keep_rows.sort()
        genes = [genes_all[i] for i in keep_rows]
        counts = combined[keep_rows, :].T.toarray().astype(np.float32)   # cells x unique symbols

        return self.build_chunk(chunk["id"], counts, genes, self.cell_line,
                                pert, time_h, factor_as_token=True)


# Registry used by the orchestrator to build sources from config.
SOURCE_REGISTRY: dict[str, type[DataSource]] = {
    "synthetic": SyntheticSource,
    "tahoe": TahoeSource,
    "sciplex": SciplexSource,
    "reprogramming": ReprogrammingSource,
    "gill": GillReprogrammingSource,
    "gse242423": GSE242423SingleCellSource,
}
