################################################################################
# Parse FastQC files for ChIP-seq quality control
# Author: Luis | Date: 2024-12-02 | Version: 1.0.0
################################################################################
# PURPOSE: Generate tab-delimited summaries for ChIP-seq experiments by converting fastqc output.
#
# USAGE:
# 2. Run from cli with experiment-id args.
#
# DEPENDENCIES: FastQC, submit_fastqc.sh
#
# OUTPUTS:
#   - Module-specific tab-delimited files
################################################################################
source(file.path(Sys.getenv("HOME"), "lab_utils", "core_scripts", "functions_for_script_control.R"))

# Parse arguments and validate configurations
args <- optparse::parse_args(commandArgs(trailingOnly = TRUE))
experiment_id <- args[["experiment-id"]]
source(file.path("~/data", experiment_id, "documentation", 
                paste0(experiment_id, "_configuration_experiment_bmc")))
validate_configs(c("RUNTIME_CONFIG", "EXPERIMENT_CONFIG"))
# Load and Validate Experiment Configuration
################################################################################
# Bootstrap phase
bootstrap_path <- normalizePath("~/lab_utils/core_scripts/functions_for_file_operations.R", 
                              mustWork = FALSE)
if (!file.exists(bootstrap_path)) {
    stop(sprintf("[FATAL] Bootstrap file not found: %s", bootstrap_path))
}
source(bootstrap_path)

# Define required dependencies
required_modules <- list(
    list(
        path = "~/lab_utils/failsafe_scripts/functions_for_logging.R",
        description = "BMC Configuration",
        required = TRUE
    )
)

# Validate module structure
stopifnot(
    "modules must have required fields" = all(sapply(required_modules, function(m) {
        all(c("path", "description", "required") %in% names(m))
    }))
)

# Load dependencies with status tracking
load_status <- lapply(required_modules, function(module) {
    if (RUNTIME_CONFIG$debug_verbose) {
        cat(sprintf("\n[LOADING] %s\n", module$description))
    }
    
    success <- safe_source(module$path, verbose = TRUE)
    
    if (!success && module$required) {
        stop(sprintf(
            "[FATAL] Failed to load required module: %s\n  Path: %s",
            module$description, module$path
        ))
    } else if (!success) {
        warning(sprintf(
            "[WARNING] Optional module not loaded: %s\n  Path: %s",
            module$description, module$path
        ))
    }
    
    return(list(
        module = module$description,
        path = module$path,
        loaded = success
    ))
})

# Display loading summary using ASCII
if (RUNTIME_CONFIG$debug_verbose) {
    cat("\n=== Module Loading Summary ===\n")
    invisible(lapply(load_status, function(status) {
        cat(sprintf(
            "[%s] %s\n    Path: %s\n",
            if(status$loaded) "+" else "-",
            status$module,
            status$path
        ))
    }))
}

#-------------------------------------------------------------------------------
# Directory Setup and Validation
#-------------------------------------------------------------------------------
experiment_id <- "241007Bel"  # !! UPDATE THIS
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
qc_dir <- file.path(base_dir, FASTQC_CONFIG$path_qc_dir)

stopifnot(
    "Base directory does not exist" = dir.exists(base_dir),
    "Quality control directory does not exist" = dir.exists(qc_dir)
)

#-------------------------------------------------------------------------------
# File Discovery
#-------------------------------------------------------------------------------
fastqc_files <- list.files(
    qc_dir,
    pattern = FASTQC_CONFIG$file_pattern,
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
    if (grepl(FASTQC_CONFIG$parse_header, first_line)) {
        version <- gsub(FASTQC_CONFIG$VERSION_PATTERN, "", first_line)
        fastqc_versions[[basename(file_path)]] <- version
    } else {
        warning(sprintf("File %s does not start with expected FastQC header", 
                       basename(file_path)))
    }
}

# Check if all versions match
expected_version <- FASTQC_CONFIG$version_required
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

if (RUNTIME_CONFIG$debug_verbose) {
    message("FastQC version check complete")
    message(sprintf("Number of files checked: %d", length(fastqc_versions)))
}

# Extract and validate sample IDs
sample_ids <- character(length(fastqc_files))
invalid_format_files <- character(0)

# Extract sample IDs with direct pattern matching
sample_ids <- gsub(
    pattern = FASTQC_CONFIG$file_pattern,
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

if (RUNTIME_CONFIG$debug_verbose) {
    message("Sample ID validation passed:")
    print(data.frame(
        file = basename(fastqc_files),
        sample_id = sample_ids,
        stringsAsFactors = FALSE
    ))
}

#-------------------------------------------------------------------------------
# Process Files
#-------------------------------------------------------------------------------
files_to_process <- if (RUNTIME_CONFIG$process_single_file) {
    RUNTIME_CONFIG$process_file_index
} else {
    seq_along(fastqc_files)
}

stopifnot(
    `Files to process must be within valid range` = all(files_to_process <= length(fastqc_files)),
    `Files to process cannot be negative` = all(files_to_process > 0)
)

for (file_idx in files_to_process) {
    current_file <- fastqc_files[file_idx]
    
    if (RUNTIME_CONFIG$debug_verbose) {
        message(sprintf("\nProcessing file %d of %d: %s", 
                       file_idx, 
                       length(files_to_process), 
                       basename(current_file)))
    }
    
    # Read and parse file
    lines <- readLines(current_file)
    module_starts <- which(grepl(FASTQC_CONFIG$parse_module_start, lines))
    module_ends <- which(grepl(FASTQC_CONFIG$parse_module_end, lines))
    module_starts <- module_starts[!(module_starts %in% module_ends)]
    
    stopifnot(
        `File must contain content` = length(lines) > 0,
        `File must contain FastQC modules` = length(module_starts) > 0,
        `Module start and end markers must match` = length(module_starts) == length(module_ends),
        `Module markers must be in correct order` = all(module_starts < module_ends),
        `Module must have valid content` = all(module_starts > 0 & module_ends <= length(lines))
        #`Module must contain data lines` = all(module_ends - module_starts > 2)
    )

    if (RUNTIME_CONFIG$debug_verbose) {
        message(sprintf("Found %d modules", length(module_starts)))
    }
    
    output_dir <- dirname(current_file)
    fastqc_summary <- list()
    
    # Process each module
    for (module_idx in seq_along(module_starts)) {
        if (RUNTIME_CONFIG$debug_verbose) {
            message(sprintf("\n  Processing module %d of %d", 
                          module_idx, 
                          length(module_starts)))
        }
        
        module_lines <- lines[module_starts[module_idx]:module_ends[module_idx]]
        module_summary <- gsub(FASTQC_CONFIG$parse_module_start, "", module_lines[1])
        fastqc_summary <- append(fastqc_summary, module_summary)
        
        # Parse module data
        module_name <- gsub(" ", "", strsplit(module_summary, "\t")[[1]][1])
        FASTQC_CONFIG$module_list <- unique(c(FASTQC_CONFIG$names, module_name))

        if (RUNTIME_CONFIG$debug_verbose) {
            message(sprintf("    Module name: %s", module_name))
            message(sprintf("    Module summary: %s", module_summary))
            message(sprintf("    Module lines: %d", length(module_lines)))
        }
        
        # Find potential headers
        module_data <- module_lines[2:(length(module_lines)-1)]
        potential_headers <- which(grepl(paste0("^", FASTQC_CONFIG$parse_prefix), module_data))
        if (RUNTIME_CONFIG$debug_verbose) {
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
                
                if (RUNTIME_CONFIG$debug_verbose) {
                    message(sprintf("    Checking header at line %d:", header_idx))
                    message(sprintf("      Header elements: %d", length(header_elements)))
                    message(sprintf("      Data line elements: %d", length(last_line_elements)))
                }
                
                if(length(header_elements) == length(last_line_elements)) {
                    if (RUNTIME_CONFIG$debug_verbose) {
                        message("    Found matching header structure")
                    }
                    
                    header <- gsub(FASTQC_CONFIG$parse_prefix, "", module_data[header_idx])
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
                     
                    if (RUNTIME_CONFIG$debug_verbose) {
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
                    FASTQC_CONFIG$file_suffix)
        )
        # Save module data if we successfully parsed it
        if (!is.null(data) && !RUNTIME_CONFIG$output_dry_run) {
            existing_files <- find_timestamped_files(output_file)
            
            if (length(existing_files) >= FASTQC_CONFIG$version_max) {
                if (RUNTIME_CONFIG$debug_verbose) {
                    cat("[SKIP] Analysis output limit reached. Existing versions:\n")
                    invisible(lapply(existing_files, function(f) cat(sprintf("  %s\n", basename(f)))))
                }
            } else {
                safe_write_file(
                    data = data,
                    path = output_file,
                    write_fn = write.table,
                    verbose = RUNTIME_CONFIG$debug_verbose,
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
                FASTQC_CONFIG$file_suffix)
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

    if (RUNTIME_CONFIG$debug_verbose) {
        message("\n  Summary data:")
        print(head(summary_data))
    }

    # Save summary for this file
    if (!RUNTIME_CONFIG$output_dry_run) {
        existing_files <- find_timestamped_files(summary_file)
        if (length(existing_files) >= FASTQC_CONFIG$version_max) {
            if (RUNTIME_CONFIG$debug_verbose) {
                cat("[SKIP] Analysis output limit reached. Existing versions:\n")
                invisible(lapply(existing_files, function(f) cat(sprintf("  %s\n", basename(f)))))
            }
        } else {
            safe_write_file(
                data = data,
                path = summary_file,
                write_fn = write.table,
                verbose = RUNTIME_CONFIG$debug_verbose,
                interactive = FALSE,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE
            )
        }
    } else {
        if (RUNTIME_CONFIG$debug_verbose) {
            message(sprintf("  Would write summary to: %s", basename(summary_file)))
        }
    }
}

message("\nProcessing complete")

if (!RUNTIME_CONFIG$output_dry_run) {
    if (!file.exists(FASTQC_CONFIG$path_module_ref)) {
        success <- safe_write_file(
            data = FASTQC_CONFIG$module_list,
            path = FASTQC_CONFIG$path_module_ref,
            write_fn = saveRDS,
            verbose = RUNTIME_CONFIG$debug_verbose,
            interactive = FALSE
        )
        
        if (!success) {
            warning("Failed to save FastQC module reference")
        } else {
            # Simple validation
            tryCatch({
                readRDS(FASTQC_CONFIG$path_module_ref)
                if (RUNTIME_CONFIG$debug_verbose) {
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

if (RUNTIME_CONFIG$debug_verbose) {
    print_config_settings(RUNTIME_CONFIG)
}
