#!/usr/bin/env Rscript
###############################################################################
# Plot bigwig files generated from bam from same sample but differentially processed
################################################################################
# PURPOSE: Plot bigwig files to see if blacklist filtering helps samples with many reads
#          at the ends or in blacklisted regions
# USAGE: ./plot_genome_tracks_of_different_processed_bam_files.R --experiment-id=<experiment-id> <options>
# DEPENDENCIES: GenomicRanges, rtracklayer
# OUTPUT: svg plots with comparisons for blacklist and non-blacklisted files.
# AUTHOR: LEMR
# DATE: 2025-05-09
################################################################################
current_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
# Bootstrap phase
FUNCTION_FILENAMES <- c("logging", "script_control", "file_operations")
for (function_filename in FUNCTION_FILENAMES) {
    function_filepath <- sprintf("~/lab_utils/core_scripts/functions_for_%s.R", function_filename)
    normalized_path <- normalizePath(function_filepath)
    if (!file.exists(normalized_path)) {
        stop(sprintf("[FATAL] File with functions not found: %s", normalized_path))
    }
    source(normalized_path)
}
message("Bootstrap phase completed...")

# Proceed if packages are installed. Can be disable.
REQUIRED_PACKAGES <- c("rtracklayer", "GenomicRanges", "Gviz")
check_required_packages(REQUIRED_PACKAGES, verbose = TRUE, skip_validation = FALSE)
message("Packages confirmed...")

# End message
message("Script completed succesfully...")
