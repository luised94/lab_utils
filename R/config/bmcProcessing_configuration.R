CONFIG <- list(
    ID_REGEX = "\\d{5,6}",
    FILE_PATTERNS = list(
        FASTQ = "*.fastq",
        SAMPLE_TABLE = "sample_table"
    ),
    PATHS = list(
        BASE_DATA = file.path(Sys.getenv("HOME"), "data"),
        DOCUMENTATION = "documentation",
        FASTQ = "fastq"
    )
)
