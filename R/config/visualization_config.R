#' Add to existing visualization configurations
CONFIG <- list(
    MODES = list(
        INTERACTIVE = list(
            DEFAULT_DIR = "240808Bel",
            DEFAULT_CHROMOSOME = 10
        )
    ),
    
    SYNC = list(
        COMMAND = "rsync -nav %s:%s/%s/plots/* %s/%s/plots/",
        USER = "username",
        DOMAIN = "domain",
        REMOTE_BASE = "~/data",
        LOCAL_BASE = "/local/dir"
    )
)
