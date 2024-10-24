#' Script Information Functions
get_script_info <- function() {
    log_info("Getting script information")
    
    script_path <- get_script_path()
    
    list(
        path = script_path,
        dir = dirname(script_path),
        name = get_script_name(script_path),
        config = get_script_config(script_path)
    )
}

get_script_path <- function() {
    # Try command line args first
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", cmd_args, value = TRUE)
    
    if (length(file_arg) > 0) {
        return(normalizePath(sub("--file=", "", file_arg)))
    }
    
    # Try sys.frames for sourced scripts
    if (!is.null(sys.frames()[[1]]$ofile)) {
        return(normalizePath(sys.frames()[[1]]$ofile))
    }
    
    log_warning("Unable to determine script path")
    return("interactive")
}

get_script_name <- function(script_path) {
    if (script_path == "interactive") {
        return("interactive")
    }
    
    tools::file_path_sans_ext(basename(script_path))
}

get_script_config <- function(script_path) {
    script_name <- get_script_name(script_path)
    
    if (script_name == "interactive") {
        return(NULL)
    }
    
    CONFIG$SCRIPTS[[script_name]]
}
