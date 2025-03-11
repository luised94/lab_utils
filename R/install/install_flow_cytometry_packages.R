#!/usr/bin/env Rscript

message("Starting installation script for flow cytometry packages...")

# Define repository root
repository_root <- "~/lab_utils"

# Define required packages
packages_to_install <- c("flowCore")
repository <- "github::"
user <- "RGLab/"
packages <- c("RProtoBufLib", "cytolib", "flowCore")
packages_to_install <- paste0(repository, user, packages)
renv::install(c(packages_to_install, "ggcyto"))

if (!dir.exists(repository_root)) {
  stop("Error: Repository root not found at ", repository_root)
}

# Activate renv
message("Activating renv environment...")

renv_path <- file.path(repository_root, "renv/activate.R")
if (!file.exists(renv_path)) {
    stop("Error: renv environment not found at ", renv_path)
}
source(renv_path)

# Install and verify packages
message("Installing and verifying genomics packages...")
for (package in packages_to_install) {
    if (!requireNamespace(package, quietly = TRUE)) {
      message("Installing package: ", package)
      renv::install(package)
      if (!requireNamespace(package, quietly = TRUE)) {
        stop("Error: Failed to install or load package '", package, "'")
      }
    } else {
      message("Package '", package, "' is already installed.")
    }
}

# Display success message
message("\nAll packages installed and verified successfully!")
