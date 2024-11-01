#' Add to existing configurations
CONFIG <- list(
    PLOT = list(
        COLORS = list(
            DEFAULT = "#fd0036",
            CONTROL = "#808080"
        ),
        TRACK_TYPE = "l",
        DATE_FORMAT = "%Y%m%d%H%M%S"
    ),
    
    PATHS = list(
        PLOTS = "plots",
        BIGWIG = "bigwig"
    ),
    
    PATTERNS = list(
        COMPARISON_PREFIX = "^comp_",
        DEFAULT_BIGWIG = "S288C_log2ratio",
        TIMECOURSE = "comp_timecourse1108"
    ),
    
    CATEGORIES = list(
        DEFAULT = c(
            "strain_source",
            "rescue_allele",
            "mcm_tag",
            "antibody",
            "timepoint_after_release"
        )
    )
)
#' Add to existing configurations
CONFIG <- list(
    PACKAGES = list(
        REQUIRED = c(
            "QuasR",
            "GenomicAlignments",
            "Gviz",
            "rtracklayer",
            "ShortRead",
            "tidyverse",
            "gtools"
        )
    ),
    
    DEFAULTS = list(
        CHROMOSOME = 10,
        GENOME_DIR = "REFGENS",
        GENOME_PATTERN = "S288C_refgenome.fna",
        FEATURE_PATTERN = "eaton_peaks",
        BIGWIG_PATTERN = "_bamcomp.bw"
    ),
    
    SYNC = list(
        COMMAND = "rsync -nav username@domain:~/data/%s/plots/* /local/dir/%s/plots/",
        REMOTE_USER = "username",
        REMOTE_DOMAIN = "domain"
    )
)
#' Add to existing visualization configurations
CONFIG <- list(
    VISUALIZATION = list(
        COLORS = list(
            SAMPLE = "#E41A1C",
            CONTROL = "#377EB8",
            HIGHLIGHT = list(
                FILL = "#FFE3E6",
                BORDER = "#FF0000",
                ALPHA = 0.3
            )
        ),
        TRACK_TYPES = list(
            LINE = "l",
            HIGHLIGHT = "highlight"
        ),
        LIMITS = list(
            Y_MIN = 0,
            Y_MAX = 100000
        ),
        OUTPUT = list(
            FORMAT = "svg",
            DATE_FORMAT = "%Y%m%d%H%M%S"
        )
    ),
    
    PATHS = list(
        SUBDIRS = list(
            PLOTS = "plots",
            BIGWIG = "bigwig"
        )
    )
)
