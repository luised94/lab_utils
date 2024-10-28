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
