#!/usr/bin/env Rscript

#' Project Configuration Definition
#' @export PROJECT_CONFIG
PROJECT_CONFIG <- list(
    SYSTEM = list(
        VERSION = "1.0.0",
        R_MIN_VERSION = "4.2.0",
        PLATFORM = .Platform$OS.type
    ),
    
    PATHS = list(
        ROOT = normalizePath("~/lab_utils"),
        R_BASE = file.path(normalizePath("~/lab_utils"), "R"),
        FUNCTIONS = "functions",
        SCRIPTS = "scripts",
        CONFIG = "config",
        DATA = normalizePath("~/data"),
        LOGS = normalizePath("~/logs")
    ),
    
    LOAD_SEQUENCE = list(
        CRITICAL = c(
            "logging_utils.R",
            "environment_utils.R",
            "validation_utils.R"
        )
    ),
    FILE_TYPES = list(
            NGS = list(
                EXTENSIONS = c(
                    BAM = "\\.bam$",
                    FASTQ = "\\.(fastq|fq)(\\.gz)?$",
                    BIGWIG = "\\.bw$",
                    BED = "\\.bed$",
                    NARROWPEAK = "\\.narrowPeak$",
                    MOTIF = "\\.(meme|pwm|jaspar)$"
                ),
                REQUIRED_INDEX = c(
                    BAM = "\\.bai$",
                    BIGWIG = "\\.bw\\.tbi$"
                )
            )
    ),
    VALIDATION = list(
         LIMITS = list(
                CHROMOSOME = c(min = 1, max = 16),
                READ_LENGTH = c(min = 20, max = 150),
                PEAK_SCORE = c(min = 0, max = 1000)
            ),
            REQUIRED_DIRS = c(
                "peak",
                "alignment",
                "plots",
                "documentation",
                "fastq",
                "processedFastq"
            )
        )
)


#    PATTERNS = list(
#        R_FILES = "\\.R$",
#        EXCLUDE = c("^\\.", "^_", "test_", "example_")
#    ),
#
#    ENVIRONMENT = list(
#        REQUIRED_VARS = c(
#            "HOME",
#            "R_LIBS_USER"
#        ),
#        RENV_SETTINGS = list(
#            AUTO_ACTIVATE = TRUE,
#            SNAPSHOT_INTERVAL = 86400  # 24 hours
