################################################################################
# Parse fastqc files for an experiment.
################################################################################
# PURPOSE:
#
# USAGE:
#   1. Use '/!!' in vim/neovim to jump to required updates
#
# !! ----> REQUIRED UPDATES:
#
# STRUCTURE:
#
# VALIDATION:
#
# DEPENDENCIES:
#
# COMMON ISSUES:
#
# AUTHOR: Luis
# DATE: 2024-12-02
# VERSION: 1.0.0
#
################################################################################
################################################################################
# Configuration and Debug Settings
################################################################################
# !! Review debug configuration
DEBUG_CONFIG <- list(
    enabled = TRUE,
    verbose = TRUE,
    interactive = TRUE,
    dry_run = FALSE,
    files_to_process_idx = 1
)

FASTQC_CONFIG <- list(
    VERSION = "1.0.0",
    FASTQC_DATA_PATTERN = "fastqc_data",
    OUTPUT_SUFFIX = ".tab",
    QC_SUBDIR = "quality_control"
)

#CONFIG$COLUMN_NAMES <- list(
#    SUMMARY = c("Stat", "Value")
#)

# Time formatting configuration
TIME_CONFIG <- list(
    timestamp_format = "%Y%m%d_%H%M%S",  # YYYYMMDD_HHMMSS
    date_format = "%Y%m%d"               # YYYYMMDD
)

# Generate timestamps once at script start
TIMESTAMPS <- list(
    full = format(Sys.time(), TIME_CONFIG$timestamp_format),
    date = format(Sys.Date(), TIME_CONFIG$date_format)
)


# !! Update experiment_id to denote the data directory that will be processed.
experiment_id <- "241007Bel"
base_experiment_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
stopifnot(
    "Base experiment not set." = !is.null(base_experiment_dir),
    "Directory does not exist." = dir.exists(base_experiment_dir)
)

##########
quality_control_dir <- file.path(base_experiment_dir, FASTQC_CONFIG$QC_SUBDIR)
stopifnot(
    "Quality control directory does not exist." = dir.exists(quality_control_dir)
)

fastqc_files <- list.files(
    quality_control_dir,
    pattern = FASTQC_CONFIG$FASTQC_DATA_PATTERN,
    recursive = TRUE,
    full.names = TRUE
)

stopifnot(
    "No quality control files found in experiment directory." = length(fastqc_files) > 0
)

files_to_process_indexes <- if (DEBUG_CONFIG$enabled) {
    DEBUG_CONFIG$files_to_process_idx
} else {
    seq_along(fastqc_files)
}

for (sample_idx in files_to_process_indexes) {
    if (DEBUG_CONFIG$verbose) {
        message("\nProcessing sample ", sample_idx)
    }

    lines <- readLines(fastqc_files[sample_idx])
    modules_starts <- which(grepl(">>", lines))
    modules_ends <- which(grepl(">>END_MODULE", lines))
    modules_starts <- modules_starts[!(modules_starts %in% modules_ends)]
    output_dir <- dirname(fastqc_files[file_idx])
    fastqc_summary <- list()
    for (module_idx in seq_along(modules_starts)) {
        module_lines <- lines[modules_starts[module_idx]:modules_ends[module_idx]]
        module_summary <- gsub(">>", "", module_lines[1])
        fastqc_summary <- append(fastqc_summary, module_summary)
        module_filename <- gsub(" ", "", strsplit(module_summary, "\t")[[1]][1])
        potential_headers <- which(grepl("^#", module_lines[2:length(module_lines)-1]))
        last_potential_header <- potential_headers[length(potential_headers)] 
        for (header_idx in potential_headers) {
            number_of_elements_in_header <- length(strsplit(module_lines[header_idx], "\t")[[1]])
            number_of_elements_in_line <- length(strsplit(module_lines[last_potential_header], "\t")[[1]])
            if(number_of_elements_in_header == number_of_elements_in_line) {
                header <- gsub("#", "", module_lines[header_idx])
                if (DEBUG_CONFIG$dry_run) {
                    data <- read.table(text = module_lines[(header_idx+1):length(module_lines)-1],
                                   header = FALSE,
                                   col.names = strsplit(header, "\t")[[1]],
                                   sep = "\t")
                }
            }

        }
        output_file_name <- paste0(current_time, "_", "fastqc_", module_filename, ".tab")
        output_file_path <- file.path(output_dir, output_file_name)
        write.table(data, file = output_file_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
    fastqc_summary <- read.table(text = unlist(fastqc_summary),
                       header = FALSE,
                       col.names = c("Stat", "Value"),
                       sep = "\t")
    print(head(fastqc_summary))
    output_file_name <- paste0(current_time, "_", "fastqc_", "summary", ".tab")
    output_file_path <- file.path(output_dir, output_file_name)
    write.table(fastqc_summary, file = output_file_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

#base::sprintf("%s_%s_%s%s",
#    get_timestamp(),
#    prefix,
#    module_name,
#    FASTQC_CONFIG$OUTPUT_SUFFIX)
#
