#!/usr/bin/env Rscript

source("../functions/package_installer.R")
source("../functions/environment_manager.R")

#' Main installation function
install_all_packages <- function(config = CONFIG) {
    log_info("Starting package installation")
    
    # Clean environment
    clean_environment()
    
    # Install bioinformatics packages
    for (group in names(config$BIOINFORMATICS)) {
        log_info("Installing", group, "packages")
        install_package_group(
            config$BIOINFORMATICS[[group]],
            "BiocManager"
        )
    }
    
    # Install development packages
    for (group in names(config$DEVELOPMENT)) {
        log_info("Installing", group, "packages")
        install_package_group(
            config$DEVELOPMENT[[group]],
            "renv"
        )
    }
    
    # Install GitHub packages
    install_package_group(
        config$GITHUB_PACKAGES,
        "github"
    )
    
    # Create snapshot
    renv::snapshot()
    
    log_info("Package installation completed")
}

if (!interactive()) {
    install_all_packages()
}
