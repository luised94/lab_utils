# R/config/core_config.R
#' Core Configuration Settings
#' @export CORE_CONFIG
CORE_CONFIG <- list(
    # System
    SYSTEM = list(
        VERSION = "1.0.0",
        R_MIN_VERSION = "4.2.0",
        LOCK_TIMEOUT = 30,
        LOCK_RETRIES = 3
    ),
    # Logging
    LOGGING = list(
        LEVELS = c("TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL"),
        DEFAULT_ROOT = "~/logs",
        FORMAT = "[%s] [%s] [%s] %s",  # timestamp, level, context, message
        MAX_SIZE = 10e6  # 10MB
    ),
    # Paths
    PATHS = list(
        ROOT = "~/lab_utils",
        MODULES = "~/lab_utils/R/modules",
        CONFIG = "~/lab_utils/R/config",
        DATA = "~/data"
    ),
    MODULES = list(
        REQUIRED = c("bmc", "ngs", "genome"),
        OPTIONAL = c("viz")
    )
)
