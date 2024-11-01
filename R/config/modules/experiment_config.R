#' Add to existing configurations
CONFIG <- list(
    EXPERIMENT = list(
        PATHS = list(
            DROPBOX = "/mnt/c/Users/%s/Dropbox (MIT)/",
            CONFIG = "sampleGridConfig.R"
        ),
        
        DIRECTORIES = c(
            "peak",
            "fastq",
            "alignment",
            "qualityControl",
            "bigwig",
            "plots",
            "documentation"
        ),
        
        OUTPUT = list(
            ENABLED = FALSE,
            CONFIG_TEMPLATE = "%s_%s.R"
        )
    ),
    
    ENVIRONMENT = list(
        REQUIRED_VARS = c(
            "WINDOWS_USER"
        )
    )
)
