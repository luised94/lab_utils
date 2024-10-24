#' Environment Management Functions
clean_environment <- function() {
    log_info("Cleaning R environment")
    
    # Detach all packages
    detach_all_packages()
    
    # Clear workspace
    clear_workspace()
    
    # Garbage collection
    gc()
}

detach_all_packages <- function() {
    log_info("Detaching all packages")
    
    attached <- paste0("package:",
                      names(sessionInfo()$otherPkgs))
    
    for (pkg in attached) {
        tryCatch({
            detach(pkg, character.only = TRUE,
                  unload = TRUE, force = TRUE)
        }, error = function(e) {
            log_warning("Failed to detach:", pkg)
        })
    }
}

clear_workspace <- function() {
    log_info("Clearing workspace")
    rm(list = ls(all.names = TRUE))
}
