"""
Project 10 — Central Configuration
===================================
All paths, thresholds, and constants live here.
Every pipeline script imports from this module.

Path resolution (no hardcoded paths):
  1. If PROJECT10_DATA_DIR is set, use that for data.
  2. Otherwise, use <repo_root>/data/

To run on your machine:
  1. Clone the repo.
  2. Place your OpenSNP files under data/ (see data/README.md).
  3. Run scripts 01 → 07 in order (see run_pipeline.py).
"""

import os
import logging
from pathlib import Path


def _find_repo_root() -> Path:
    """Walk up from this file to find the repo root (contains lib/)."""
    p = Path(__file__).resolve().parent
    # config.py is at repo root
    if (p / "lib").is_dir():
        return p
    # config.py might be inside lib/ or original_scripts/
    if (p.parent / "lib").is_dir():
        return p.parent
    return p


# ═══════════════════════════════════════════════════════════════════
# PATHS
# ═══════════════════════════════════════════════════════════════════

PROJECT_DIR = Path(
    os.environ.get("PROJECT10_PROJECT_DIR", str(_find_repo_root()))
).resolve()

DATA_DIR = Path(
    os.environ.get("PROJECT10_DATA_DIR", str(PROJECT_DIR / "data"))
).resolve()

RAW_GENO_DIR = DATA_DIR / "opensnp_genotypes_Ancestry__413files"

# ── Input files ──
RAW_ANCESTRY_CSV = DATA_DIR / "opensnp_Ancestry.csv"
MANIFEST_CSV = RAW_GENO_DIR / "manifest.csv"

# ── Ancestry mapping ──
MAPPING_CSV = DATA_DIR / "opensnp_Ancestry_unique_mapping_final_v4.csv"
GROUPINGS_CSV = DATA_DIR / "opensnp_Ancestry_final_groupings_v4.csv"

# ── Pipeline intermediates ──
MASTER_MANIFEST_CSV = DATA_DIR / "master_manifest.csv"
AUDIT_RESULTS_CSV = DATA_DIR / "audit_results_revised.csv"
PIPELINE_MANIFEST_CSV = DATA_DIR / "pipeline_manifest_revised.csv"
CONVERSION_READY_CSV = DATA_DIR / "conversion_ready_manifest.csv"
STAGE4_INPUT_CSV = DATA_DIR / "stage4_input_manifest.csv"

# ── Directories ──
PLINK_INDIVIDUAL_DIR = DATA_DIR / "plink_individual"
MERGE_DIR = DATA_DIR / "plink_merged" / "full_merge" / "strict_retry"
FILTERED_DIR = MERGE_DIR / "filtered_inputs"
QC_DIR = DATA_DIR / "plink_qc"
NETWORK_DIR = DATA_DIR / "network_analysis"
DEDUP_DIR = DATA_DIR / "deduplicated_analysis"
SENSITIVITY_DIR = DATA_DIR / "sensitivity_eur_only"
LIFTOVER_DIR = DATA_DIR / "liftover"
FIG_DIR = DATA_DIR / "figures"

# ── Merged dataset ──
MERGED_PREFIX = MERGE_DIR / "full_merged_strict"

# ── QC step prefixes ──
QC_MISSING_PREFIX = QC_DIR / "qc_missingness"
QC_STEP1_PREFIX = QC_DIR / "step1_sample_filtered"
QC_STEP2_PREFIX = QC_DIR / "step2_snp_filtered"
QC_STEP3_PREFIX = QC_DIR / "step3_maf_filtered"
QC_STEP3B_PREFIX = QC_DIR / "step3b_no_ambiguous"
QC_LD_PRUNE_PREFIX = QC_DIR / "step4_ld_prune"
QC_PRUNED_PREFIX = QC_DIR / "step4_pruned_dataset"
QC_IBS_PREFIX = QC_DIR / "step5_ibs"

# ── Dedup prefixes ──
DEDUP_PREFIX = DEDUP_DIR / "step1_pruned_deduplicated"
DEDUP_IBS_PREFIX = DEDUP_DIR / "step2_ibs_deduplicated"

# ── liftOver ──
LIFTOVER_BED_DIR = LIFTOVER_DIR / "bed_inputs"
LIFTOVER_OUT_DIR = LIFTOVER_DIR / "lifted_bed"
LIFTOVER_UNMAPPED_DIR = LIFTOVER_DIR / "unmapped"
CHAIN_FILE = LIFTOVER_DIR / "hg18ToHg19.over.chain.gz"
LIFTOVER_EXEC = "liftOver"


# ═══════════════════════════════════════════════════════════════════
# TOOLS
# ═══════════════════════════════════════════════════════════════════

PLINK_EXEC = os.environ.get("PLINK_EXEC", "plink")
N_THREADS = int(os.environ.get("PROJECT10_THREADS", "4"))


# ═══════════════════════════════════════════════════════════════════
# THRESHOLDS
# ═══════════════════════════════════════════════════════════════════

# ── Merge ──
SNP_PRESENCE_THRESHOLD = 0.90
MERGE_MAX_ATTEMPTS = 3

# ── Conversion ──
MIN_PARSED_SNPS = 10_000

# ── QC ──
SAMPLE_MISS_THRESHOLD = 0.05
SNP_MISS_THRESHOLD = 0.05
MAF_THRESHOLD = 0.01
LD_WINDOW = 50
LD_STEP = 5
LD_R2 = 0.2

# ── Relatedness ──
IBS2_PRUNE_THRESHOLD = 0.1
IBS2_RELATEDNESS_THRESHOLD = 0.4
PIHAT_3RD_DEGREE = 0.125
PIHAT_2ND_DEGREE = 0.25
PIHAT_1ST_DEGREE = 0.4
PIHAT_DUPLICATE = 0.9

# ── Network ──
PIHAT_EDGE_THRESHOLD = 0.05


# ═══════════════════════════════════════════════════════════════════
# CONTROLLED VOCABULARIES
# ═══════════════════════════════════════════════════════════════════

EXPECTED_TIER0 = {"EUR", "Non-EUR", "Admixed", "Founder", "Unknown"}

TIER0_ORDER = ["EUR", "Non-EUR", "Admixed", "Founder", "Unknown"]

TIER0_COLORS = {
    "EUR":     "#4C72B0",
    "Non-EUR": "#DD8452",
    "Admixed": "#55A868",
    "Founder": "#C44E52",
    "Unknown": "#8C8C8C",
}

AUTOSOMAL_CHROMS = {str(i) for i in range(1, 23)}

VALID_BASES = {"A", "C", "G", "T"}

# Expanded sentinel SNPs for more accurate build detection
# More sentinels across different chromosomes increases confidence
# Selected SNPs with large position differences between builds
SENTINELS = {
    # Chromosome 1
    "rs3131972":  {"chr": "1", "GRCh37": 752721, "GRCh36": 742429},
    "rs12562034": {"chr": "1", "GRCh37": 768448, "GRCh36": 758156},
    "rs9442372":  {"chr": "1", "GRCh37": 54490,  "GRCh36": 45498},
    "rs2814778":  {"chr": "1", "GRCh37": 159174683, "GRCh36": 158821376},  # 353kb diff
    "rs3094315":  {"chr": "1", "GRCh37": 742429, "GRCh36": 752721},
    # Chromosome 2
    "rs11903757": {"chr": "2", "GRCh37": 233248202, "GRCh36": 234542577},
    "rs4988235":  {"chr": "2", "GRCh37": 136608646, "GRCh36": 137973396},  # 1.3Mb diff
    "rs13006529": {"chr": "2", "GRCh37": 234542577, "GRCh36": 233248202},
    # Chromosome 3
    "rs1800628":  {"chr": "3", "GRCh37": 187011472, "GRCh36": 187370191},
    "rs11718699": {"chr": "3", "GRCh37": 186502452, "GRCh36": 186861171},
    # Chromosome 4
    "rs1426654":  {"chr": "4", "GRCh37": 100239319, "GRCh36": 100240000},
    "rs2033529":  {"chr": "4", "GRCh37": 38774644, "GRCh36": 38815467},
    # Chromosome 5
    "rs1872575":  {"chr": "5", "GRCh37": 1286516, "GRCh36": 1286516},
    "rs4703842":  {"chr": "5", "GRCh37": 131995964, "GRCh36": 132535964},
    # Chromosome 6
    "rs1800562":  {"chr": "6", "GRCh37": 26093141, "GRCh36": 26199158},
    "rs9268480":  {"chr": "6", "GRCh37": 32610116, "GRCh36": 32578064},
    # Chromosome 7
    "rs4717568":  {"chr": "7", "GRCh37": 116889148, "GRCh36": 117229406},
    "rs702689":   {"chr": "7", "GRCh37": 151162545, "GRCh36": 151502803},
    # Chromosome 8
    "rs1545390":  {"chr": "8", "GRCh37": 8079549, "GRCh36": 8162049},
    "rs13266634": {"chr": "8", "GRCh37": 118184783, "GRCh36": 119424783},
    # Chromosome 9
    "rs1333049":  {"chr": "9", "GRCh37": 22125504, "GRCh36": 22115503},
    "rs10811661": {"chr": "9", "GRCh37": 22134094, "GRCh36": 22124093},
    # Chromosome 10
    "rs10519107": {"chr": "10", "GRCh37": 52771387, "GRCh36": 54531387},
    "rs1800466":  {"chr": "10", "GRCh37": 94781004, "GRCh36": 96541004},
    # Chromosome 11
    "rs174547":   {"chr": "11", "GRCh37": 61597212, "GRCh36": 61569830},
    "rs7115089":  {"chr": "11", "GRCh37": 47509148, "GRCh36": 47481766},
    # Chromosome 12
    "rs1800497":  {"chr": "12", "GRCh37": 6936739, "GRCh36": 7018739},
    "rs2076848":  {"chr": "12", "GRCh37": 121416650, "GRCh36": 124856650},
    # Chromosome 13
    "rs1805007":  {"chr": "13", "GRCh37": 28881233, "GRCh36": 28360981},
    # Chromosome 15
    "rs1042522":  {"chr": "17", "GRCh37": 7579472, "GRCh36": 7519148},
    "rs11655081": {"chr": "15", "GRCh37": 78882925, "GRCh36": 78882925},
    # Chromosome 16
    "rs8050136":  {"chr": "16", "GRCh37": 53782363, "GRCh36": 53716713},
    "rs4785763":  {"chr": "16", "GRCh37": 53774781, "GRCh36": 53709131},
    # Chromosome 17
    "rs4988235":  {"chr": "17", "GRCh37": 48240029, "GRCh36": 46300029},
    # Chromosome 19
    "rs429358":   {"chr": "19", "GRCh37": 45411941, "GRCh36": 45416178},
    "rs7412":     {"chr": "19", "GRCh37": 45412079, "GRCh36": 45416316},
    # Chromosome 22
    "rs2413583":  {"chr": "22", "GRCh37": 43089873, "GRCh36": 41566739},
    "rs5746467":  {"chr": "22", "GRCh37": 43071077, "GRCh36": 41547943},
}

# Display and processing limits
MAX_FIRST_SNP_DISPLAY_LENGTH = 120
MIN_FILE_SIZE_MB = 1.0
AUDIT_SAMPLE_LIMIT = 10_000
MAX_AUDIT_LINES = 1_000_000  # Prevent OOM on malformed files

# Optional: keep intermediate files for debugging
KEEP_INTERMEDIATE_FILES = (
    os.environ.get("PROJECT10_KEEP_INTERMEDIATE", "false").lower() == "true"
)


# ═══════════════════════════════════════════════════════════════════
# LOGGING CONFIGURATION
# ═══════════════════════════════════════════════════════════════════

LOG_LEVEL = os.environ.get("PROJECT10_LOG_LEVEL", "INFO").upper()
LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logging(level: str = None, log_file: str = None) -> None:
    """
    Configure logging for the pipeline.

    Parameters:
        level: Log level (DEBUG, INFO, WARNING, ERROR). Defaults to LOG_LEVEL.
        log_file: Optional file path to write logs to
    """
    level = level or LOG_LEVEL
    handlers = [logging.StreamHandler()]

    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=level,
        format=LOG_FORMAT,
        datefmt=LOG_DATE_FORMAT,
        handlers=handlers,
        force=True,  # Override any existing configuration
    )

    # Reduce noise from matplotlib and other libraries
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)


# ═══════════════════════════════════════════════════════════════════
# DIRECTORY CREATION
# ═══════════════════════════════════════════════════════════════════

def ensure_dirs():
    """Create all output directories. Safe to call multiple times."""
    for d in [
        PLINK_INDIVIDUAL_DIR, MERGE_DIR, FILTERED_DIR,
        QC_DIR, NETWORK_DIR, DEDUP_DIR, SENSITIVITY_DIR,
        LIFTOVER_DIR, LIFTOVER_BED_DIR, LIFTOVER_OUT_DIR,
        LIFTOVER_UNMAPPED_DIR, FIG_DIR,
    ]:
        d.mkdir(parents=True, exist_ok=True)
