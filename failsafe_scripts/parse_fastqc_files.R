################################################################################
# Parse FastQC files for ChIP-seq quality control
################################################################################
# PURPOSE:
#   Process FastQC output files and generate parsed tab-delimited summaries for
#   ChIP-seq quality control analysis. Extracts module-specific data and creates
#   standardized output files with sample identifiers.
#
# USAGE:
#   1. Update experiment_id to point to correct data directory
#   2. Adjust DEBUG_CONFIG as needed:
#      - single_file_mode: TRUE for testing single files
#      - verbose: TRUE for detailed processing information
#      - dry_run: TRUE to check without writing files
#   3. Run script
#
# !! ----> REQUIRED UPDATES:
#   - Set experiment_id for data directory
#   - Review debug configuration for testing needs
#   - Verify FastQC version requirements
#
# STRUCTURE:
#   1. Configuration and debug settings
#   2. Directory and version validation
#   3. File discovery and sample ID extraction
#   4. Module parsing and data extraction:
#      - Basic Statistics
#      - Per base sequence quality
#      - Per sequence quality scores
#      - Other FastQC modules
#   5. Data validation and output
#
# VALIDATION:
#   - Directory structure and permissions
#   - FastQC file presence and version
#   - Sample ID format and uniqueness
#   - Module structure and content
#   - Data parsing integrity
#
# DEPENDENCIES:
#   - Base R (>= 4.2.0)
#   - FastQC (version 0.11.5)
#   - submit_fastqc.sh
#   - run_fastqc_array.sh
#
# COMMON ISSUES:
#   - Missing or incorrect quality control directory
#   - Malformed FastQC files or unexpected versions
#   - Write permission errors in output directory
#   - Inconsistent sample ID formats
#   - Memory limitations with large datasets
#
# OUTPUT:
#   - Module-specific tab-delimited files
#   - Summary statistics for each sample
#   - Processing logs (when verbose)
#   - Module Name File for future reference and mappings.
# AUTHOR: Luis
# DATE: 2024-12-02
# VERSION: 1.0.0
################################################################################
################################################################################
# Configuration and Debug Settings
################################################################################
DEBUG_CONFIG <- list( # !! UPDATE THIS
    single_file_mode = FALSE,           # Test single file in main logic.
    verbose = TRUE,           # Print processing details
    interactive = TRUE,       # Allow interactive processing
    dry_run = FALSE,         # Skip file writes
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
    qc_subdir = "quality_control",
    existing_version_limit = 1,
    module_names = character(0),
    module_reference_file = file.path(
        Sys.getenv("HOME"),
        "data",
        "fastqc_module_reference.rds"
    )
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
# Load and Validate Experiment Configuration
################################################################################
# Bootstrap phase
bootstrap_path <- normalizePath("~/lab_utils/failsafe_scripts/functions_for_file_operations.R", 
                              mustWork = FALSE)
if (!file.exists(bootstrap_path)) {
    stop(sprintf("[FATAL] Bootstrap file not found: %s", bootstrap_path))
}
source(bootstrap_path)

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
files_to_process <- if (DEBUG_CONFIG$single_file_mode) {
    DEBUG_CONFIG$files_to_process_idx
} else {
    seq_along(fastqc_files)
}

stopifnot(
    `Files to process must be within valid range` = all(files_to_process <= length(fastqc_files)),
    `Files to process cannot be negative` = all(files_to_process > 0)
)

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
    
    stopifnot(
        `File must contain content` = length(lines) > 0,
        `File must contain FastQC modules` = length(module_starts) > 0,
        `Module start and end markers must match` = length(module_starts) == length(module_ends),
        `Module markers must be in correct order` = all(module_starts < module_ends),
        `Module must have valid content` = all(module_starts > 0 & module_ends <= length(lines))
        #`Module must contain data lines` = all(module_ends - module_starts > 2)
    )

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
        FASTQC_CONFIG$module_names <- unique(c(FASTQC_CONFIG$names, module_name))

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

                    if (!is.null(data)) {
                        stopifnot(
                            `Parsed data must have rows` = nrow(data) > 0,
                            `Parsed data must have columns` = ncol(data) > 0,
                            `Column names must not be empty` = all(nchar(colnames(data)) > 0)
                        )
                    }
                     
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
            existing_files <- find_timestamped_files(output_file)
            
            if (length(existing_files) >= FASTQC_CONFIG$existing_version_limit) {
                if (DEBUG_CONFIG$verbose) {
                    cat("[SKIP] Analysis output limit reached. Existing versions:\n")
                    invisible(lapply(existing_files, function(f) cat(sprintf("  %s\n", basename(f)))))
                }
            } else {
                safe_write_file(
                    data = data,
                    path = output_file,
                    write_fn = write.table,
                    verbose = DEBUG_CONFIG$verbose,
                    interactive = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE
                )
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
        
    stopifnot(
        `Summary must contain data` = nrow(summary_data) > 0,
        `Summary must have required columns` = all(c("Stat", "Value") %in% colnames(summary_data))
    )

    if (DEBUG_CONFIG$verbose) {
        message("\n  Summary data:")
        print(head(summary_data))
    }

    # Save summary for this file
    if (!DEBUG_CONFIG$dry_run) {
        existing_files <- find_timestamped_files(summary_file)
            
        if (length(existing_files) >= FASTQC_CONFIG$existing_version_limit || DEBUG_CONFIG$verbose) {
            cat("Found existing versions:\n")
            invisible(lapply(existing_files, function(f) cat(sprintf("  %s\n", basename(f)))))
        } else {
            safe_write_file(
                data = data,
                path = summary_file,
                write_fn = write.table,
                verbose = DEBUG_CONFIG$verbose,
                interactive = FALSE,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE
            )
        }
    } else {
        if (DEBUG_CONFIG$verbose) {
            message(sprintf("  Would write summary to: %s", basename(summary_file)))
        }
    }
}

message("\nProcessing complete")

if (!DEBUG_CONFIG$dry_run) {
    if (!file.exists(FASTQC_CONFIG$module_reference_file)) {
        success <- safe_write_file(
            data = FASTQC_CONFIG$module_names,
            path = FASTQC_CONFIG$module_reference_file,
            write_fn = saveRDS,
            verbose = DEBUG_CONFIG$verbose,
            interactive = FALSE
        )
        
        if (!success) {
            warning("Failed to save FastQC module reference")
        } else {
            # Simple validation
            tryCatch({
                readRDS(FASTQC_CONFIG$module_reference_file)
                if (DEBUG_CONFIG$verbose) {
                    message("FastQC module reference was read successfully")
                }
            }, error = function(e) {
                warning("Failed to validate saved FastQC module reference")
            })
        }
    } else {
        message("FastQC module reference already exists")
    }
}
