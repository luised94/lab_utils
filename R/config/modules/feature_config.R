#' Add to existing configurations
CONFIG <- list(
    FEATURES = list(
        PATHS = list(
            BASE_DIR = file.path(Sys.getenv("HOME"), "data", "feature_files")
        ),
        
        PATTERNS = list(
            EXCLUDE = c(
                "sample-key.tab",
                "\\.rds$",
                "_converted.bed$"
            ),
            VERIFY = c(
                "\\.rds$",
                "_converted.bed$"
            )
        ),
        
        FILE_TYPES = list(
            NUCLEOSOME = "Nucleosome_calls",
            TIMING = "hawkins",
            TRANSCRIPTION = "Rhee",
            FEATURES = "SGD"
        ),
        
        COLUMNS = list(
            NUCLEOSOME = list(
                EXCLUDE = c("Nucleosome ID", "Nucleosome dyad", "Chromosome"),
                POSITION = "Nucleosome dyad"
            ),
            TIMING = list(
                EXCLUDE = c("Chromosome", "Position"),
                WINDOW = 100  # +/- bp around position
            )
        )
    )
)
