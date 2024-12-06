#!/usr/bin/env Rscript

#' Validation Utilities
#' @description Core validation functions for project
#' @export

#' Check Project Configuration
#' @param throw_error Logical Whether to throw error if config missing
#' @return Logical TRUE if config exists
check_project_config <- function(throw_error = TRUE) {
    config_exists <- exists("PROJECT_CONFIG", envir = .GlobalEnv)
    if (!config_exists && throw_error) {
        stop("PROJECT_CONFIG not loaded. Source project_config.R first.")
    }
    invisible(config_exists)
}

#' Validate Function Arguments
#' @param args List Arguments to validate
#' @param spec List Argument specifications
#' @return List Validated arguments
validate_function_args <- function(args, spec) {
    check_project_config()
    
    for (arg_name in names(spec)) {
        arg_spec <- spec[[arg_name]]
        arg_value <- args[[arg_name]]
        
        # Check required arguments
        if (is.null(arg_value)) {
            if (arg_spec$required) {
                stop(sprintf("Required argument missing: %s", arg_name))
            }
            args[[arg_name]] <- arg_spec$default
            next
        }
        
        # Type checking
        if (!inherits(arg_value, arg_spec$type)) {
            stop(sprintf(
                "Invalid type for %s: expected %s, got %s",
                arg_name, arg_spec$type, class(arg_value)[1]
            ))
        }
        
        # Custom validation
        if (!is.null(arg_spec$validation)) {
            if (!arg_spec$validation(arg_value)) {
                stop(arg_spec$error_message)
            }
        }
    }
    
    invisible(args)
}

#' Validate NGS File
#' @param file_path Character Path to file
#' @param type Character File type
#' @param config List Configuration settings
#' @return Logical TRUE if valid
validate_ngs_file <- function(
    file_path,
    type = c("bam", "fastq", "bigwig", "bed", "narrowpeak", "motif"),
    config = PROJECT_CONFIG
) {
    type <- match.arg(type)
    
    # Basic existence check
    if (!file.exists(file_path)) {
        log_error(sprintf("File not found: %s", file_path))
        return(FALSE)
    }
    
    # Extension check
    pattern <- config$FILE_TYPES$NGS$EXTENSIONS[[toupper(type)]]
    if (!grepl(pattern, file_path)) {
        log_error(sprintf("Invalid file extension for %s: %s", type, file_path))
        return(FALSE)
    }
    
    # Index check if required
    if (type %in% names(config$FILE_TYPES$NGS$REQUIRED_INDEX)) {
        index_pattern <- config$FILE_TYPES$NGS$REQUIRED_INDEX[[toupper(type)]]
        index_file <- sub(pattern, index_pattern, file_path)
        if (!file.exists(index_file)) {
            log_warning(sprintf("Index file missing: %s", index_file))
            return(FALSE)
        }
    }
    
    TRUE
}

#' Validate BAM File
#' @param file_path Character Path to BAM file
#' @return Logical TRUE if valid
validate_bam_file <- function(file_path) {
    # Basic checks
    if (!grepl("\\.bam$", file_path)) {
        log_error("Not a BAM file:", file_path)
        return(FALSE)
    }
    
    # Check index
    if (!file.exists(paste0(file_path, ".bai"))) {
        log_warning("BAM index missing:", file_path)
        return(FALSE)
    }
    
    TRUE
}

#' Example Usage Configuration
VALIDATION_CONFIG <- list(
    TYPES = c(
        "bam", "fastq", "bigwig", 
        "bed", "bedgraph", "narrowPeak"
    ),
    LIMITS = list(
        CHROMOSOME = c(min = 1, max = 16),
        READ_LENGTH = c(min = 20, max = 150)
    )
)
