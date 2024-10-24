#' Genome Configuration Parameters
CONFIG <- list(
    PATHS = list(
        BASE_DIR = file.path(Sys.getenv("HOME"), "data"),
        GENOME_DIR = "REFGENS",
        DOCUMENTATION_DIR = "documentation"
    ),
    
    PATTERNS = list(
        GENOME = "S288C_refgenome.fna",
        SAMPLE_TABLE = "sample_table",
        CONTROL_FACTOR_PREFIX = "X__cf_"
    ),
    
    CHROMOSOME_MAPPING = list(
        STYLES = c("UCSC", "Roman", "Numeric"),
        ROMAN = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                 "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"),
        PREFIX = "chr"
    ),
    
    REQUIRED_PACKAGES = c(
        "QuasR",
        "GenomicAlignments",
        "Gviz",
        "rtracklayer",
        "ShortRead",
        "tidyverse"
    )
)
CONFIG <- list(
    PATHS = list(
        FEATURE_DIR = file.path(Sys.getenv("HOME"), "data", "feature_files")
    ),
    
    FEATURE_TYPES = list(
        PEAKS = "eaton_peaks",
        ORIGINS = "origins",
        GENES = "genes"
    ),
    
    REQUIRED_CATEGORIES = list(
        BASE = "antibody",
        OPTIONAL = c("condition", "treatment", "timepoint")
    ),
    
    LABEL_CONFIG = list(
        SEPARATOR = "_",
        MAX_LENGTH = 50,
        TRUNCATE_SUFFIX = "..."
    )
)
#' Merge with existing configurations
CONFIG <- list(
    GENOME = list(
        DEFAULT_STRAND = "*",
        DEFAULT_CHROMOSOME = "chr1",
        STYLES = list(
            UCSC = "UCSC",
            ENSEMBL = "ENSEMBL"
        )
    ),
    
    VISUALIZATION = list(
        TRACKS = list(
            COLORS = c(
                PRIMARY = "#1f77b4",
                SECONDARY = "#ff7f0e",
                HIGHLIGHT = "#2ca02c"
            ),
            TYPES = list(
                DATA = "l",
                ANNOTATION = "gene"
            )
        ),
        DEFAULTS = list(
            GENOME = "hg19",
            MIN_HEIGHT = 0,
            MAX_HEIGHT = 100
        )
    )
)
