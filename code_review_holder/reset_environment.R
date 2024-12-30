#==============================================================================
# Environment Reset and Setup Script
# Purpose: Create clean, reproducible environment for R analysis
# Updated: 2024-11-14
#==============================================================================

# 1. Complete Environment Reset
#------------------------------------------------------------------------------
rm(list = ls(all.names = TRUE))  # Remove all objects, including hidden ones
gc(full = TRUE)                  # Force garbage collection
if (!is.null(dev.list())) dev.off() # Close all open graphics devices

# 2. Reset Random Seed
#------------------------------------------------------------------------------
set.seed(42)  # Set consistent seed for reproducibility

# 3. Reset Graphics Parameters
#------------------------------------------------------------------------------
par(mar = c(5.1, 4.1, 4.1, 2.1))  # Reset to R defaults
options(scipen = 999)              # Disable scientific notation
options(digits = 7)                # Reset numerical precision

# 4. Reset Working Directory (optional - uncomment if needed)
#------------------------------------------------------------------------------
# setwd("/your/project/path")

# 5. Package Management
#------------------------------------------------------------------------------
# Function to install missing packages
install_if_missing <- function(packages) {
    new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
    if (length(new_packages)) install.packages(new_packages)
    invisible()
}

# Core packages for analysis (customize as needed)
required_packages <- c(
    "tidyverse",
    "data.table",
    "Matrix",
    "parallel"
)

# Install missing packages
install_if_missing(required_packages)

# Load packages with suppressed messages
suppressMessages({
    for (pkg in required_packages) {
        library(pkg, character.only = TRUE)
    }
})

# 6. Set Global Options
#------------------------------------------------------------------------------
options(
    stringsAsFactors = FALSE,     # Prevent automatic factor conversion
    warn = 1,                     # Show warnings immediately
    encoding = "UTF-8",           # Set default encoding
    timeout = 3600               # Set timeout for downloads
)

# 7. Custom Error Handler (optional)
#------------------------------------------------------------------------------
options(error = function() {
    cat("Error occurred at:", date(), "\n")
    traceback(2)
})

# 8. Memory Management
#------------------------------------------------------------------------------
memory.limit(size = NA)  # Maximum memory available (Windows only)

# 9. Session Information
#------------------------------------------------------------------------------
cat("\nEnvironment successfully reset at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R Version:", R.version.string, "\n")
cat("Working Directory:", getwd(), "\n\n")
