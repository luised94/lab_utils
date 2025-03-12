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
# Currently RProtoBufLib gives error.
packages <- c("RProtoBufLib", "cytolib", "flowCore")
packages_to_install <- paste0(repository, user, packages)

# Activate renv
message("Activating renv environment...")

renv_activation_path <- file.path(repository_root, "renv/activate.R")
if (!file.exists(renv_activation_path)) {
    stop("Error: renv environment not found at ", renv_activation_path)
}
source(renv_activation_path)

# Install and verify packages
message("Installing and verifying genomics packages...")

for (package in packages_to_install) {
    # Check if the package is already installed
    if (requireNamespace(package, quietly = TRUE)) {
        message("Package '", package, "' is already installed.")
        next
    }

    message("Installing package: ", package)

    # Attempt to install the package
    install_result <- tryCatch(
        expr = {
            renv::install(package)
            TRUE # Returns TRUE on success
        },
        error = function(e){
            print(sprintf("Failed to install package %s", package))
            print(e)
            FALSE # Returns FALSE on failure
        }
    )

    if (!install_result) {
        stop("!!!! Failed to install package '", package, "'")
    }

    # After successful installation, ensure it's loaded
    library_name <- gsub(paste0(repository, user), "", x = package, fixed = TRUE)
    message("Loading package '", library_name, "'...")
    if (!library(library_name, character.only = TRUE, logical.return = TRUE)) {
        stop("!!!! Failed to load package '", package, "' after installation.")
    }
    message("Package '", package, "' has been successfully installed and loaded.")
}
# Display success message
message("\nAll packages installed and verified successfully!")
#renv::snapshot()
