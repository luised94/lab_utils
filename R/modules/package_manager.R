#' Package Management Functions
load_required_packages <- function(packages = CONFIG$PACKAGES$REQUIRED) {
    log_info("Loading required packages")
    
    suppressPackageStartupMessages({
        results <- sapply(packages, require, character.only = TRUE)
    })
    
    missing <- packages[!results]
    if (length(missing) > 0) {
        log_error("Missing packages:", paste(missing, collapse = ", "))
        stop("Required packages not available")
    }
    
    return(TRUE)
}
#' Package Management Functions
initialize_environment <- function(config = CONFIG) {
    log_info("Initializing R environment")
    
    # Initialize renv
    setup_renv()
    
    # Setup package management
    setup_package_manager()
    
    # Install and load packages
    install_required_packages(config$PACKAGES)
    
    # Snapshot environment
    snapshot_environment()
    
    log_info("Environment initialization complete")
}

setup_renv <- function() {
    log_info("Setting up renv")
    
    if (!requireNamespace("renv", quietly = TRUE)) {
        log_info("Installing renv")
        install.packages("renv")
    }
    
    renv::init()
}

setup_package_manager <- function() {
    log_info("Setting up BiocManager")
    
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        log_info("Installing BiocManager")
        install.packages("BiocManager")
    }
}

install_required_packages <- function(package_lists) {
    log_info("Installing required packages")
    
    # Combine all packages
    all_packages <- unique(unlist(package_lists))
    
    # Install packages
    result <- tryCatch({
        BiocManager::install(all_packages)
        TRUE
    }, error = function(e) {
        log_error("Package installation failed:", e$message)
        FALSE
    })
    
    if (!result) {
        stop("Package installation failed")
    }
}

load_packages <- function(packages) {
    log_info("Loading packages")
    
    results <- sapply(packages, function(pkg) {
        tryCatch({
            library(pkg, character.only = TRUE)
            TRUE
        }, error = function(e) {
            log_warning("Failed to load package:", pkg)
            FALSE
        })
    })
    
    if (!all(results)) {
        failed <- names(results)[!results]
        log_warning("Failed to load packages:", 
                   paste(failed, collapse = ", "))
    }
}

snapshot_environment <- function() {
    log_info("Creating environment snapshot")
    
    tryCatch({
        renv::snapshot()
    }, error = function(e) {
        log_error("Failed to create snapshot:", e$message)
    })
}
