#' Package Installation Functions
install_package_group <- function(packages,
                                manager = "BiocManager") {
    log_info("Installing package group using:", manager)
    
    results <- switch(manager,
        "BiocManager" = install_bioc_packages(packages),
        "renv" = install_renv_packages(packages),
        "github" = install_github_packages(packages),
        stop("Unknown package manager:", manager)
    )
    
    failed <- names(results)[!results]
    if (length(failed) > 0) {
        log_warning("Failed to install packages:",
                   paste(failed, collapse = ", "))
    }
    
    invisible(results)
}

install_bioc_packages <- function(packages) {
    sapply(packages, function(pkg) {
        tryCatch({
            BiocManager::install(pkg, update = FALSE)
            TRUE
        }, error = function(e) {
            log_error("Failed to install:", pkg)
            FALSE
        })
    })
}

install_renv_packages <- function(packages) {
    sapply(packages, function(pkg) {
        tryCatch({
            renv::install(pkg)
            TRUE
        }, error = function(e) {
            log_error("Failed to install:", pkg)
            FALSE
        })
    })
}

install_github_packages <- function(packages) {
    sapply(packages, function(pkg) {
        tryCatch({
            remotes::install_github(pkg)
            TRUE
        }, error = function(e) {
            log_error("Failed to install:", pkg)
            FALSE
        })
    })
}
