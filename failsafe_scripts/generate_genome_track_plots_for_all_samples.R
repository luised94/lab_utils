#!/usr/bin/env Rscript
# Constants and Configuration
#-----------------------------------------------------------------------------
PLOT_CONFIG <- list(
    SAMPLES_PER_PAGE = 4,
    DEFAULT_CHROMOSOME = 10,
    TRACK_COLOR = "#fd0036",
    WIDTH = 10,
    HEIGHT = 8
)
EXPERIMENT_ID <- "241007Bel"
TIMESTAMP <- format(Sys.Date(), "%Y%m%d")
ref_genome_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    pattern = "S288C_refgenome.fna",
    full.names = TRUE,
    recursive = TRUE
)[1]
feature_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "feature_files"),
    pattern = "eaton_peaks",
    full.names = TRUE
)[1]
REQUIRED_PACKAGES <- c("rtracklayer", "GenomicRanges", "Gviz", "tidyverse")
# 1. Load required packages with verification
# 2. Set up initial variables and paths
# 3. Load and verify config and functions
# 4. Load and process sample table with new processing steps
# Enforce factor levels from config
# Sort metadata using config column order
# 5. Load reference genome
# 6. Create genome ranges
# 7. Load feature file
# Load and adjust feature chromosome style
# Function for bigwig validation
# Setup chromosome
# Determine the maximum and minimum across all sample tracks. 
# Plot all genome tracks samples four at time using the previously determined maximum and minimum.
