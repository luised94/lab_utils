################################################################################
# Analyze the files produced by the preprocessing test scripts
################################################################################
# PURPOSE: 
# Conclusion: {{fill}}
# USAGE: source("test_analyze_preprocessing_files_using_R.R")
# DEPENDENCIES: GenomicRanges, rtracklayer
# OUTPUT: {{fill}}
# AUTHOR: LEMR
# DATE: 2025-02-25
################################################################################
#!/usr/bin/env Rscript
# Bootstrap phase
function_filenames <- c("logging", "script_control", "file_operations")
for (function_filename in function_filenames) {
    function_filepath <- sprintf("~/lab_utils/core_scripts/functions_for_%s.R", function_filename)
    normalized_path <- normalizePath(function_filepath)
    if (!file.exists(normalized_path)) {
        stop(sprintf("[FATAL] File with functions not found: %s", normalized_path))
    }
    source(normalized_path)
}

# Proceed if packages are installed. Can be disable.
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
check_required_packages(required_packages, verbose = TRUE, skip_validation = args$skip_validation)
