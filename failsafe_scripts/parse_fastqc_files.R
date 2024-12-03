################################################################################
# Parse fastqc files for an experiment.
################################################################################
# PURPOSE:
#   Process FastQC output files and generate parsed tab-delimited summaries
#
# USAGE:
#   1. Update experiment_id to point to correct data directory
#   2. Adjust DEBUG_CONFIG as needed
#   3. Run script
#
# !! ----> REQUIRED UPDATES:
#   - Set experiment_id
#   - Review debug configuration
#
# STRUCTURE:
#   1. Configuration blocks
#   2. Directory validation
#   3. File discovery
#   4. Module parsing
#   5. Data output
#
# VALIDATION:
#   - Directory existence
#   - FastQC file presence
#   - Module structure
#
# DEPENDENCIES:
#   - Base R
#
# COMMON ISSUES:
#   - Missing quality control directory
#   - Malformed FastQC files
#   - Write permission errors
#
# AUTHOR: Luis
# DATE: 2024-12-02
# VERSION: 1.0.0
################################################################################

################################################################################
# Configuration and Debug Settings
################################################################################
DEBUG_CONFIG <- list( # !! UPDATE THIS
    enabled = TRUE,           # Enable debug mode
    verbose = TRUE,           # Print processing details
    interactive = TRUE,       # Allow interactive processing
    dry_run = TRUE,         # Skip file writes
    files_to_process_idx = 1  # Process specific files in debug mode
)

FASTQC_CONFIG <- list(
    VERSION = "0.11.5",                    # Expected FastQC version
    VERSION_PATTERN = "^##FastQC\\s+",     # Pattern to match version line
    HEADER_PATTERN = "^##FastQC",          # Pattern to identify FastQC header
    module_separator = ">>",
    module_end = ">>END_MODULE",
    header_prefix = "#",
    fastqc_pattern = "fastqc_data",
    output_suffix = ".tab",
    qc_subdir = "quality_control"
)

FASTQC_CONFIG$FILE_PATTERN <- list(
    REGEX = "consolidated_([0-9]{5,6})_sequence_fastqc_data\\.txt$",
    EXPECTED_FORMAT = "consolidated_XXXXXX_sequence_fastqc_data.txt"  # For error messages
)

TIME_CONFIG <- list(
    timestamp_format = "%Y%m%d_%H%M%S",
    date_format = "%Y%m%d"
)

# Generate timestamps
TIMESTAMPS <- list(
    full = format(Sys.time(), TIME_CONFIG$timestamp_format),
    date = format(Sys.Date(), TIME_CONFIG$date_format)
)

################################################################################
# Directory Setup and Validation
################################################################################
experiment_id <- "241007Bel"  # !! UPDATE THIS
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
qc_dir <- file.path(base_dir, FASTQC_CONFIG$qc_subdir)

stopifnot(
    "Base directory does not exist" = dir.exists(base_dir),
    "Quality control directory does not exist" = dir.exists(qc_dir)
)

################################################################################
# File Discovery
################################################################################
fastqc_files <- list.files(
    qc_dir,
    pattern = FASTQC_CONFIG$fastqc_pattern,
    recursive = TRUE,
    full.names = TRUE
)

stopifnot(
    "No FastQC files found" = length(fastqc_files) > 0
)

# Check FastQC versions across files
fastqc_versions <- list()
for (file_path in fastqc_files) {
    first_line <- readLines(file_path, n = 1)
    if (grepl(FASTQC_CONFIG$HEADER_PATTERN, first_line)) {
        version <- gsub(FASTQC_CONFIG$VERSION_PATTERN, "", first_line)
        fastqc_versions[[basename(file_path)]] <- version
    } else {
        warning(sprintf("File %s does not start with expected FastQC header", 
                       basename(file_path)))
    }
}

# Check if all versions match
expected_version <- FASTQC_CONFIG$VERSION
version_check <- sapply(fastqc_versions, function(v) v == expected_version)

if (!all(version_check)) {
    different_versions <- fastqc_versions[!version_check]
    warning(sprintf(
        "Found different FastQC versions:\nExpected: %s\nDifferent versions found in:\n%s",
        expected_version,
        paste(sprintf("- %s: %s", 
                     names(different_versions), 
                     unlist(different_versions)), 
              collapse = "\n")
    ))
}

if (DEBUG_CONFIG$verbose) {
    message("FastQC version check complete")
    message(sprintf("Number of files checked: %d", length(fastqc_versions)))
}

# Extract and validate sample IDs
sample_ids <- character(length(fastqc_files))
invalid_format_files <- character(0)

# Extract sample IDs with direct pattern matching
sample_ids <- gsub(
    pattern = FASTQC_CONFIG$FILE_PATTERN$REGEX,
    replacement = "\\1",
    x = basename(fastqc_files)
)

# Validate specific failure modes with descriptive errors
invalid_format_files <- basename(fastqc_files)[
    !grepl("consolidated_.*_sequence_fastqc_data\\.txt$", basename(fastqc_files))
]
invalid_id_files <- basename(fastqc_files)[!grepl("^[0-9]{5,6}$", sample_ids)]
duplicate_ids <- sample_ids[duplicated(sample_ids)]

# Clear validation with specific error messages
stopifnot(
    `Files not matching expected format` = length(invalid_format_files) == 0,
    `Files with invalid sample IDs` = length(invalid_id_files) == 0,
    `Duplicate sample IDs found` = length(duplicate_ids) == 0,
    `All sample IDs must be extracted` = !any(sample_ids == ""),
    `Number of sample IDs must match number of files` = length(fastqc_files) == length(sample_ids)
)

if (DEBUG_CONFIG$verbose) {
    message("Sample ID validation passed:")
    print(data.frame(
        file = basename(fastqc_files),
        sample_id = sample_ids,
        stringsAsFactors = FALSE
    ))
}

################################################################################
# Process Files
################################################################################
files_to_process <- if (DEBUG_CONFIG$enabled) {
    DEBUG_CONFIG$files_to_process_idx
} else {
    seq_along(fastqc_files)
}

for (file_idx in files_to_process) {
    current_file <- fastqc_files[file_idx]
    
    if (DEBUG_CONFIG$verbose) {
        message(sprintf("\nProcessing file %d of %d: %s", 
                       file_idx, 
                       length(files_to_process), 
                       basename(current_file)))
    }
    
    # Read and parse file
    lines <- readLines(current_file)
    module_starts <- which(grepl(FASTQC_CONFIG$module_separator, lines))
    module_ends <- which(grepl(FASTQC_CONFIG$module_end, lines))
    module_starts <- module_starts[!(module_starts %in% module_ends)]
    
    if (DEBUG_CONFIG$verbose) {
        message(sprintf("Found %d modules", length(module_starts)))
    }
    
    output_dir <- dirname(current_file)
    fastqc_summary <- list()
    
    # Process each module
    for (module_idx in seq_along(module_starts)) {
        if (DEBUG_CONFIG$verbose) {
            message(sprintf("\n  Processing module %d of %d", 
                          module_idx, 
                          length(module_starts)))
        }
        
        module_lines <- lines[module_starts[module_idx]:module_ends[module_idx]]
        module_summary <- gsub(FASTQC_CONFIG$module_separator, "", module_lines[1])
        fastqc_summary <- append(fastqc_summary, module_summary)
        
        # Parse module data
        module_name <- gsub(" ", "", strsplit(module_summary, "\t")[[1]][1])
        if (DEBUG_CONFIG$verbose) {
            message(sprintf("    Module name: %s", module_name))
            message(sprintf("    Module summary: %s", module_summary))
            message(sprintf("    Module lines: %d", length(module_lines)))
        }
        
        # Find potential headers
        module_data <- module_lines[2:(length(module_lines)-1)]
        potential_headers <- which(grepl(paste0("^", FASTQC_CONFIG$header_prefix), module_data))
        if (DEBUG_CONFIG$verbose) {
            message(sprintf("    Found %d potential headers", length(potential_headers)))
            if (length(potential_headers) > 0) {
                message("    Headers content:")
                for(h_idx in potential_headers) {
                    message(sprintf("      Line %d: %s", h_idx, module_data[h_idx]))
                }
            }
        }
        
        # Only process if we found headers
        if (length(potential_headers) > 0) {
            last_potential_header <- potential_headers[length(potential_headers)]
            data <- NULL
            
            # Validate headers against data structure
            for (header_idx in potential_headers) {
                header_elements <- strsplit(module_data[header_idx], "\t")[[1]]
                last_line_elements <- strsplit(module_data[last_potential_header], "\t")[[1]]
                
                if (DEBUG_CONFIG$verbose) {
                    message(sprintf("    Checking header at line %d:", header_idx))
                    message(sprintf("      Header elements: %d", length(header_elements)))
                    message(sprintf("      Data line elements: %d", length(last_line_elements)))
                }
                
                if(length(header_elements) == length(last_line_elements)) {
                    if (DEBUG_CONFIG$verbose) {
                        message("    Found matching header structure")
                    }
                    
                    header <- gsub(FASTQC_CONFIG$header_prefix, "", module_data[header_idx])
                    # Adjust indexes to skip the module header line
                    data_start_idx <- header_idx + 1
                    data_end_idx <- length(module_data)
                    
                    data <- read.table(
                        text = module_data[data_start_idx:data_end_idx],
                        header = FALSE,
                        col.names = strsplit(header, "\t")[[1]],
                        sep = "\t"
                    )
                     
                    if (DEBUG_CONFIG$verbose) {
                        message(sprintf("    Parsed data: %d rows, %d columns", 
                                        nrow(data), ncol(data)))
                        message(sprintf("\n    Data preview for module %s:", module_name))
                        print(head(data))
                    }
                } else {
                    message("Skip header.")
                }
            }
        }
            
        output_file <- file.path(
            output_dir,
            sprintf("%s_%s_fastqc_%s%s", 
                    TIMESTAMPS$full,
                    sample_ids[file_idx],  # Add sample ID to filename
                    module_name,
                    FASTQC_CONFIG$output_suffix)
        )
        # Save module data if we successfully parsed it
        if (!is.null(data) && !DEBUG_CONFIG$dry_run) {
             
            write.table(
                data,
                file = output_file,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE
            )
             
            if (DEBUG_CONFIG$verbose) {
                message(sprintf("    Wrote module data to: %s", 
                                basename(output_file)))
            }
        } else {
            message(sprintf("    Would write module data to: %s", 
                    basename(output_file)))
        }
    }
    
    summary_file <- file.path(
        output_dir,
        sprintf("%s_%s_fastqc_summary%s",
                TIMESTAMPS$full,
                sample_ids[file_idx],  # Add sample ID to filename
                FASTQC_CONFIG$output_suffix)
    )

    summary_data <- read.table(
        text = unlist(fastqc_summary),
        header = FALSE,
        col.names = c("Stat", "Value"),
        sep = "\t"
    )
        
    if (DEBUG_CONFIG$verbose) {
        message("\n  Summary data:")
        print(head(summary_data))
    }

    # Save summary for this file
    if (!DEBUG_CONFIG$dry_run) {
        
        write.table(
            summary_data,
            file = summary_file,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
        )
        
        if (DEBUG_CONFIG$verbose) {
            message(sprintf("  Wrote summary to: %s", basename(summary_file)))
        }
    } else {
        if (DEBUG_CONFIG$verbose) {
            message(sprintf("  Would write summary to: %s", basename(summary_file)))
        }
    }
}

message("\nProcessing complete")
