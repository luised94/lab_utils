# init_config.R
CONFIG <- list(
    PATHS = list(
        BASE = Sys.getenv("R_PROJECT_ROOT", "R"),
        FUNCTIONS = "functions",
        SCRIPTS = "scripts",
        CONFIG = "config",
        TEMPLATES = "templates",
        LOGS = "logs"
    ),
    
    LOAD_ORDER = list(
        PRIORITY = c(
            "logging_utils.R",
            "environment_utils.R",
            "validation_utils.R"
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
            "R_PROJECT_ROOT"
        )
    ),
    
    CLUSTER = list(
        ENABLED = as.logical(Sys.getenv("R_CLUSTER_MODE", "FALSE")),
        MOUNT_POINT = Sys.getenv("R_CLUSTER_MOUNT", "/cluster/data")
    )
)
cat(CONFIG$PATHS$BASE)
