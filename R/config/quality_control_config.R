#' Add to existing QC configurations
CONFIG <- list(
    FASTQC = list(
        PATTERNS = list(
            DATA_FILE = "fastqc_data",
            MODULE_START = "^>>",
            MODULE_END = ">>END_MODULE",
            HEADER = "^#"
        ),
        
        OUTPUT = list(
            DATE_FORMAT = "%Y-%m-%d-%H-%M-%S",
            SEPARATOR = "\t",
            EXTENSIONS = list(
                MODULE = ".tab",
                SUMMARY = "_summary.tab"
            )
        ),
        
        COLUMNS = list(
            SUMMARY = c("Stat", "Value")
        )
    ),
    
    PATHS = list(
        BASE_DIR = file.path(Sys.getenv("HOME"), "data"),
        QC_SUBDIR = "qualityControl"
    )
)
