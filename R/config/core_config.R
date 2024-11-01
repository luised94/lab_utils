# R/config/core_config.R

#' Core Configuration Settings
#' @export CORE_CONFIG
CORE_CONFIG <- list(
    SYSTEM = list(
        VERSION = "1.0.0",
        R_MIN_VERSION = "4.2.0"
    ),
    
    PATHS = list(
        ROOT = "~/lab_utils",
        DATA = "~/data",
        LOGS = "~/logs"
    ),
    
    MODULES = list(
        REQUIRED = c("bmc", "ngs", "genome"),
        OPTIONAL = c("viz")
    )
)
