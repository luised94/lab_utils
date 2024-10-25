
#' Initialization Configuration
CONFIG <- list(
    INITIALIZATION = list(
        PATHS = list(
            BASE = "R",
            FUNCTIONS = "functions",
            SCRIPTS = "scripts",
            CONFIG = "config",
            TEMPLATES = "templates"
        ),
        
        PATTERNS = list(
            R_FILES = "\\.R$",
            CONFIG_FILES = "^[^.].*\\.R$",
            EXCLUDE = c("^\\.", "^_", "test_", "example_")
        ),
        
        LOAD_ORDER = list(
            PRIORITY = c(
                "logging_utils.R",
                "environment_utils.R",
                "validation_utils.R"
            ),
            OPTIONAL = c(
                "analysis_utils.R",
                "visualization_utils.R"
            )
        )
    ),
    
    ENVIRONMENT = list(
        REQUIRED_VARS = c(
            "HOME",
            "R_LIBS_USER",
            "WINDOWS_USER"
        ),
        
        PATHS = list(
            USER_LIB = "~/R/library",
            RENV = "renv"
        )
    )
)
#' Initialization Configuration
CONFIG <- list(
    PATHS = list(
        BASE = "R",
        FUNCTIONS = "functions",
        SCRIPTS = "scripts",
        CONFIG = "config",
        TEMPLATES = "templates"
    ),
    
    LOAD_ORDER = list(
        PRIORITY = c(
            "logging_utils.R",
            "environment_utils.R",
            "validation_utils.R"
        ),
        OPTIONAL = c(
            "analysis_utils.R",
            "visualization_utils.R"
        )
    ),
    
    PATTERNS = list(
        R_FILES = "\\.R$",
        EXCLUDE = c("^\\.", "^_", "test_", "example_")
    ),
    
    ENVIRONMENT = list(
        REQUIRED_VARS = c(
            "HOME",
            "R_LIBS_USER",
            "WINDOWS_USER"
        )
    )
)
