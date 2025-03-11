#!/usr/bin/env Rscript

message("Starting installation script for flow cytometry packages...")

# Define repository root
repository_root <- "~/lab_utils"
if (!dir.exists(repository_root)) {
  stop("Error: Repository root not found at ", repository_root)
}

# Define required packages
# In R console with activated renv
# Because of issues with boost libraries (updated methods are not the same as methods used by flowCore), you must install a particular version of boost library.
#remove.packages(c("BH", "flowCore", "cytolib"))

renv::install("BH", version = "1.75.0-0")
repository <- "github::"
user <- "RGLab/"
packages <- c("RProtoBufLib", "cytolib", "flowCore")
packages_to_install <- paste0(repository, user, packages)

# Activate renv
message("Activating renv environment...")

renv_activation_path <- file.path(repository_root, "renv/activate.R")
if (!file.exists(renv_path)) {
    stop("Error: renv environment not found at ", renv_path)
}
source(renv_activation_path)

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
#renv::snapshot()
