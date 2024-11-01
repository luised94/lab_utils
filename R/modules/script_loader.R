#' Script Loading Functions
load_project_scripts <- function(config = CONFIG$INITIALIZATION) {
    log_info("Loading project scripts")
    
    # Load in specific order
    load_priority_scripts(config)
    load_remaining_scripts(config)
    
    log_info("Script loading completed")
}

load_priority_scripts <- function(config) {
    log_info("Loading priority scripts")
    
    base_dir <- get_project_root()
    
    for (script in config$LOAD_ORDER$PRIORITY) {
        script_path <- find_script(script, base_dir)
        if (!is.null(script_path)) {
            source_script_safely(script_path)
        }
    }
}

load_remaining_scripts <- function(config) {
    log_info("Loading remaining scripts")
    
    # Get all script directories
    dirs <- get_script_directories(config)
    
    # Load scripts from each directory
    for (dir in dirs) {
        load_directory_scripts(dir, config)
    }
}

source_script_safely <- function(script_path) {
    log_info("Sourcing:", basename(script_path))
    
    tryCatch({
        source(script_path)
        TRUE
    }, error = function(e) {
        log_error("Failed to source:", script_path)
        log_error("Error:", e$message)
        FALSE
    })
}
#' Script Loading Functions
load_directory_scripts <- function(dir_path) {
    log_info("Loading scripts from:", dir_path)
    
    if (!dir.exists(dir_path)) {
        log_warning("Directory not found:", dir_path)
        return(FALSE)
    }
    
    # Get all R scripts
    scripts <- list.files(
        dir_path,
        pattern = CONFIG$PATTERNS$R_FILES,
        full.names = TRUE
    )
    
    # Filter excluded patterns
    scripts <- filter_scripts(scripts)
    
    # Load each script
    for (script in scripts) {
        source_script(script)
    }
    
    TRUE
}

filter_scripts <- function(scripts) {
    # Remove excluded patterns
    for (pattern in CONFIG$PATTERNS$EXCLUDE) {
        scripts <- scripts[!grepl(pattern, basename(scripts))]
    }
    scripts
}

source_script <- function(script_path) {
    log_info("Sourcing:", basename(script_path))
    
    tryCatch({
        source(script_path)
        TRUE
    }, error = function(e) {
        log_error("Failed to source:", script_path)
        log_error("Error:", e$message)
        FALSE
    })
}

load_priority_scripts <- function() {
    log_info("Loading priority scripts")
    
    for (script in CONFIG$LOAD_ORDER$PRIORITY) {
        script_path <- file.path(CONFIG$PATHS$FUNCTIONS, script)
        if (file.exists(script_path)) {
            source_script(script_path)
        }
    }
}
