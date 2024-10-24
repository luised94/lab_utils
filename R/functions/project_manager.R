#' Project Management Functions
get_project_root <- function() {
    root <- tryCatch({
        rprojroot::find_root(
            rprojroot::has_file("renv.lock")
        )
    }, error = function(e) {
        getwd()
    })
    
    normalizePath(root)
}

get_script_directories <- function(config) {
    base_dir <- get_project_root()
    
    paths <- config$PATHS
    
    sapply(paths, function(path) {
        file.path(base_dir, path)
    })
}

find_script <- function(script_name,
                       base_dir,
                       pattern = CONFIG$INITIALIZATION$PATTERNS$R_FILES) {
    log_info("Finding script:", script_name)
    
    matches <- list.files(
        base_dir,
        pattern = script_name,
        recursive = TRUE,
        full.names = TRUE
    )
    
    if (length(matches) == 0) {
        log_warning("Script not found:", script_name)
        return(NULL)
    }
    
    matches[1]
}
