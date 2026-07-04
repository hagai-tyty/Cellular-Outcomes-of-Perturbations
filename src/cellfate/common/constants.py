"""Project-wide constants: the single source of truth for invariants that every
package depends on. Nothing here may be redefined elsewhere.

If you change a value that affects on-disk artefacts (class order, gene count,
fingerprint size, artefact filenames), bump SCHEMA_VERSION -- it is a breaking
change for every existing shard and bundle.
"""

from __future__ import annotations

from enum import StrEnum

# --------------------------------------------------------------------------- #
# Versioning                                                                   #
# --------------------------------------------------------------------------- #
SCHEMA_VERSION: str = "1.0"

# --------------------------------------------------------------------------- #
# Outcome classes  (ORDER IS FIXED -- never reorder)                          #
#   0 = identity preserved (the cell was refreshed but kept its function)      #
#   1 = identity loss      (unwanted reprogramming / dedifferentiation)        #
#   2 = apoptosis          (toxicity / death)                                  #
# --------------------------------------------------------------------------- #
CLASSES: tuple[str, str, str] = ("safe", "loss", "death")
N_CLASSES: int = len(CLASSES)
CLASS_TO_IDX: dict[str, int] = {c: i for i, c in enumerate(CLASSES)}
IDX_TO_CLASS: dict[int, str] = {i: c for i, c in enumerate(CLASSES)}

SAFE_IDX: int = CLASS_TO_IDX["safe"]
LOSS_IDX: int = CLASS_TO_IDX["loss"]
DEATH_IDX: int = CLASS_TO_IDX["death"]


class Modality(StrEnum):
    """Perturbation modality; selects the perturbation encoder (Document 3)."""

    CHEM = "chem"        # small molecule (Morgan fingerprint)
    GENETIC = "genetic"  # CRISPR target gene (network embedding)
    TF = "tf"            # transcription-factor cocktail / OSKM (factor-set embedding)


class Split(StrEnum):
    """Dataset split. ``calib`` is reserved exclusively for the conformal quantile."""

    TRAIN = "train"
    VAL = "val"
    CALIB = "calib"
    TEST = "test"
    DROP = "drop"  # excluded rows (used by the 'both-unseen' regime)


class Regime(StrEnum):
    """Generalisation regimes, reported separately by evaluation (Document 5)."""

    SCAFFOLD = "scaffold"      # leave-drug-out (Bemis-Murcko scaffold)
    CELL_LINE = "cell_line"    # leave-cell-line / cell-type-out
    BOTH = "both"              # unseen scaffold AND unseen line


# --------------------------------------------------------------------------- #
# Feature-space sizes                                                          #
# --------------------------------------------------------------------------- #
DEFAULT_N_GENES: int = 2000        # G: length of the model input vector X (HVGs)
N_FINGERPRINT_BITS: int = 2048     # Morgan fingerprint length
MORGAN_RADIUS: int = 2
N_DOSE_TIME: int = 2               # [log10(dose_uM), log(time_h)]

# --------------------------------------------------------------------------- #
# Data-source policy                                                           #
# --------------------------------------------------------------------------- #
# Sources on which the aging clock is out-of-distribution: the dAge label is
# masked (age_mask = False) on these, while the safety head still trains on them.
CANCER_SOURCES: frozenset[str] = frozenset({"tahoe"})

# --------------------------------------------------------------------------- #
# Canonical artefact filenames / directory names (centralised so no string is  #
# duplicated across packages). See cellfate.common.io.ArtifactPaths.           #
# --------------------------------------------------------------------------- #
SHARDS_DIRNAME: str = "shards"
MANIFEST_PARTS_DIRNAME: str = "manifest_parts"
MANIFEST_FILENAME: str = "manifest.parquet"
SPLITS_DIRNAME: str = "splits"
SCALERS_FILENAME: str = "scalers.json"
PROGRESS_FILENAME: str = "progress_tracker.json"

BUNDLE_DIRNAME: str = "bundle"
BUNDLE_MEMBERS_DIRNAME: str = "members"
BUNDLE_OOD_DIRNAME: str = "ood"
BUNDLE_META_FILENAME: str = "meta.json"
BUNDLE_CONFIG_FILENAME: str = "config.yaml"
BUNDLE_TEMPERATURE_FILENAME: str = "temperature.json"
BUNDLE_CONFORMAL_FILENAME: str = "conformal.json"
BUNDLE_RES_FILENAME: str = "res_params.json"
BUNDLE_METRICS_FILENAME: str = "metrics.json"

REPORTS_DIRNAME: str = "reports"

# --------------------------------------------------------------------------- #
# Cell-cycle gene sets (Tirosh et al. 2016) for proliferation deconfounding.    #
# Shared because both data labelling and any cell-cycle QC use them.            #
# --------------------------------------------------------------------------- #
S_GENES: tuple[str, ...] = (
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP",
    "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
    "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
    "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
)
G2M_GENES: tuple[str, ...] = (
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
    "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
    "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1",
    "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
)

# --------------------------------------------------------------------------- #
# Default identity / safety signature gene sets (cell-type agnostic defaults;   #
# data labelling overrides 'safe' with cell-type-specific markers via config).  #
# --------------------------------------------------------------------------- #
DEFAULT_SIGNATURES: dict[str, tuple[str, ...]] = {
    "safe": ("COL1A1", "THY1", "VIM", "FN1"),                 # somatic identity (placeholder)
    "loss": ("NANOG", "POU5F1", "LIN28A", "SOX2", "ZFP42"),   # pluripotency / dedifferentiation
    "death": ("BAX", "CASP3", "FAS", "BBC3", "PMAIP1"),       # apoptosis
}

# Genes used to COMPUTE the fate label. They are held OUT of the model's input
# panel (see fit_gene_panel must_exclude) so the network cannot read its own
# label off its input -- the safety head must predict fate from other genes +
# the perturbation. The aging clock still sees them (it reads the full profile).
LABEL_HOLDOUT: tuple[str, ...] = tuple(
    dict.fromkeys(g for gs in DEFAULT_SIGNATURES.values() for g in gs)
)

# Transcription-factor vocabulary for the TF/genetic perturbation encoder. A TF
# cocktail is encoded as a multi-hot vector over this vocabulary (scaled by dose),
# so OSKM and OSK share factors and the model learns factor identity -- not an
# opaque token. Order is fixed (it defines the vector layout).
TF_VOCAB: tuple[str, ...] = ("POU5F1", "SOX2", "KLF4", "MYC", "NANOG", "LIN28A", "GLIS1", "SALL4")

# Yamanaka-style cocktail shorthands -> their canonical factor sets.
TF_COCKTAILS: dict[str, tuple[str, ...]] = {
    "OSKM": ("POU5F1", "SOX2", "KLF4", "MYC"),
    "OSK": ("POU5F1", "SOX2", "KLF4"),
    "OSKMNL": ("POU5F1", "SOX2", "KLF4", "MYC", "NANOG", "LIN28A"),
}
# Single-letter aliases used to parse arbitrary cocktail strings.
TF_LETTERS: dict[str, str] = {
    "O": "POU5F1", "S": "SOX2", "K": "KLF4", "M": "MYC", "N": "NANOG", "L": "LIN28A",
}
