#' Merge with existing configurations
CONFIG <- list(
    PATHS = list(
        BASE_DIR = Sys.getenv("HOME"),
        SUBDIRS = list(
            ALIGNMENT = "alignment",
            BIGWIG = "bigwig",
            DOCUMENTATION = "documentation"
        )
    ),
    
    PATTERNS = list(
        SAMPLE_TABLE = "sample_table",
        CONTROL_FACTOR_PREFIX = "X__cf_",
        REFERENCE_GENOME = "S288C",
        BAM_SUFFIX = ".bam$"
    ),
    
    CONTROL = list(
        MAX_CONTROLS = 1,
        DEFAULT_INDEX = 1,
        ANTIBODY_COLUMN = "antibody",
        INPUT_VALUE = "Input"
    )
)
