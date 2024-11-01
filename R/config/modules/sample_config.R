#' Add to existing configurations
CONFIG <- list(
    PATHS = list(
        SUBDIRS = list(
            ALIGNMENT = "alignment",
            BIGWIG = "bigwig"
        )
    ),
    
    FILES = list(
        PATTERNS = list(
            BAM = ".bam$",
            REFERENCE = "S288C"
        )
    ),
    
    VALIDATION = list(
        REQUIRED_COLUMNS = c(
            "sample_ID",
            "short_name"
        )
    )
)
