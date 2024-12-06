#' Safe File Writing Function
#' @param data The data to write (or source path for file.copy)
#' @param path Destination file path
#' @param write_fn The writing function to use. Supported functions:
#'                 - write.csv(x, file, ...)
#'                 - write.table(x, file, ...)
#'                 - file.copy(from, to, ...)
#'                 - saveRDS(object, file, ...)
#' @param verbose Boolean to control logging output
#' @param interactive Boolean to control user prompts for overwriting
#' @param ... Additional arguments passed to write_fn
#' @return Boolean indicating success
safe_write_file <- function(data, path, write_fn, verbose = FALSE, interactive = TRUE, ...) {
    # Convert to absolute path
    tryCatch({
        abs_path <- normalizePath(path, mustWork = FALSE)
    }, error = function(e) {
        cat("[ERROR] Failed to resolve path. Check path validity.\n")
        return(FALSE)
    })

    # Input validation 
    stopifnot(
        "write_fn must be a function" = is.function(write_fn),
        "verbose must be logical" = is.logical(verbose) && length(verbose) == 1,
        "interactive must be logical" = is.logical(interactive) && length(interactive) == 1,
        "parent directory must be writable" = file.access(dirname(abs_path), mode = 2) == 0
    )

    # Create wrapper for supported functions
    # Handles differences in parameter ordering for different file function.
    # Only allows explicitly handled functions.
    wrapped_write <- function(data, path, ...) {
        if (identical(write_fn, write.csv) || identical(write_fn, write.table)) {
            write_fn(x = data, file = path, ...)
        } else if (identical(write_fn, file.copy)) {
            write_fn(from = data, to = path, ...)
        } else if (identical(write_fn, saveRDS)) {
            write_fn(object = data, file = path, ...)
        } else {
            stop(sprintf(
                "Unsupported writing function: %s. Supported functions are: write.csv, write.table, file.copy, saveRDS",
        } else {
            stop(sprintf(
                "Unsupported writing function: %s. Supported functions are: write.csv, write.table, file.copy",
                deparse(substitute(write_fn))
            ))
        }
    }
    if (verbose) {
        cat(sprintf("[INFO] Processing file: %s\n", abs_path))
    }

    do_write <- function() {
        tryCatch({
            wrapped_write(data, abs_path, ...)
            if (verbose) {
                cat(sprintf("[WROTE] File to: %s\n", abs_path))
            }
            return(TRUE)
        }, error = function(e) {
            cat(sprintf("[ERROR] Failed to write file: %s\n", e$message))
            return(FALSE)
        })
    }

    # If file doesn't exist or non-interactive, write immediately
    if (!file.exists(abs_path) || !interactive) {
        return(do_write())
    }

    # Handle existing file in interactive mode
    if (verbose) {
        cat(sprintf("[WARNING] File already exists: %s\n", abs_path))
    }

    user_input <- readline(prompt="File exists. Overwrite? (y/n): ")
    if (tolower(user_input) == "y") {
        return(do_write())
    }

    # Handle non-overwrite cases
    message <- if(tolower(user_input) == "n") {
        "[WARNING SKIP] No overwrite. Did not write: %s\n"
    } else {
        "[WARNING SKIP] Option not recognized. Did not write: %s\n"
    }
    cat(sprintf(message, abs_path))
    return(FALSE)
}

#' Safe Source Function for R Scripts
#' @param file_path Path to the R script to source
#' @param verbose Boolean to control logging output
#' @param local Logical; evaluate sourced code in local or global environment
#' @param chdir Logical; change directory when sourcing
#' @param ... Additional arguments passed to source()
#' @return Boolean indicating success of sourcing operation
safe_source <- function(file_path, verbose = FALSE, local = FALSE, chdir = FALSE, ...) {
    # Input validation
    stopifnot(
        "file_path must be a character string" = is.character(file_path) && length(file_path) == 1,
        "verbose must be logical" = is.logical(verbose) && length(verbose) == 1,
        "local must be logical" = is.logical(local) && length(local) == 1,
        "chdir must be logical" = is.logical(chdir) && length(chdir) == 1,
        "file must exist" = file.exists(file_path),
        "file must be readable" = file.access(file_path, mode = 4) == 0
    )

    # Convert to absolute path
    tryCatch({
        abs_path <- normalizePath(file_path, mustWork = FALSE)
    }, error = function(e) {
        cat("[ERROR] Failed to resolve path. Check path validity.\n")
        return(FALSE)
    })


    # Attempt to source the file
    if (verbose) {
        cat(sprintf("[INFO] Attempting to source: %s\n", abs_path))
    }

    tryCatch({
        source(abs_path, local = local, chdir = chdir, ...)
        if (verbose) {
            cat(sprintf("[SUCCESS] Successfully sourced: %s\n", abs_path))
        }
        return(TRUE)
    }, error = function(e) {
        # Enhanced error messages for common source() errors
        error_msg <- switch(class(e)[1],
            "parseError" = "Syntax error in script",
            "evaluationError" = "Error evaluating script",
            e$message
        )
        cat(sprintf("[ERROR] Failed to source %s: %s\n", abs_path, error_msg))
        return(FALSE)
    }, warning = function(w) {
        if (verbose) {
            cat(sprintf("[WARNING] Warning while sourcing %s: %s\n", abs_path, w$message))
        }
        return(TRUE)
    })
}

#' Find Files with Timestamp Pattern
#' @param base_path Path to check for timestamped versions
#' @param full_paths Boolean to return full paths or just filenames
#' @return Character vector of matching files
find_timestamped_files <- function(base_path, full_paths = TRUE) {
    # Input validation
    stopifnot(
        "base_path must be a character string" = is.character(base_path) && length(base_path) == 1,
        "full_paths must be logical" = is.logical(full_paths) && length(full_paths) == 1
    )

    # Path normalization
    base_path <- normalizePath(base_path, mustWork = FALSE)
    
    # Extract path components
    dir_path <- dirname(base_path)
    base_name <- basename(base_path)
    
    # Timestamp regex pattern
    timestamp_regex <- "\\d{4}\\d{2}\\d{2}_\\d{2}\\d{2}\\d{2}"
    
    # Remove timestamp from base name
    base_pattern <- gsub(timestamp_regex, "", base_name, perl = TRUE)
    
    # Create search pattern
    search_pattern <- sprintf("^.*%s$", base_pattern)
    
    # Find matching files
    matches <- list.files(
        dir_path,
        pattern = search_pattern,
        full.names = full_paths
    )
    
    return(matches)
}

#' Safe Source Function for R Scripts
#' @param file_path Path to the R script to source
#' @param verbose Boolean to control logging output
#' @param local Logical; evaluate sourced code in local or global environment
#' @param chdir Logical; change directory when sourcing
#' @param ... Additional arguments passed to source()
#' @return Boolean indicating success of sourcing operation
safe_source <- function(file_path, verbose = FALSE, local = FALSE, chdir = FALSE, ...) {
    # Input validation
    stopifnot(
        "file_path must be a character string" = is.character(file_path) && length(file_path) == 1,
        "verbose must be logical" = is.logical(verbose) && length(verbose) == 1,
        "local must be logical" = is.logical(local) && length(local) == 1,
        "chdir must be logical" = is.logical(chdir) && length(chdir) == 1,
        "file must exist" = file.exists(file_path),
        "file must be readable" = file.access(file_path, mode = 4) == 0
    )

    # Convert to absolute path
    tryCatch({
        abs_path <- normalizePath(file_path, mustWork = FALSE)
    }, error = function(e) {
        cat("[ERROR] Failed to resolve path. Check path validity.\n")
        return(FALSE)
    })


    # Attempt to source the file
    if (verbose) {
        cat(sprintf("[INFO] Attempting to source: %s\n", abs_path))
    }

    tryCatch({
        source(abs_path, local = local, chdir = chdir, ...)
        if (verbose) {
            cat(sprintf("[SUCCESS] Successfully sourced: %s\n", abs_path))
        }
        return(TRUE)
    }, error = function(e) {
        # Enhanced error messages for common source() errors
        error_msg <- switch(class(e)[1],
            "parseError" = "Syntax error in script",
            "evaluationError" = "Error evaluating script",
            e$message
        )
        cat(sprintf("[ERROR] Failed to source %s: %s\n", abs_path, error_msg))
        return(FALSE)
    }, warning = function(w) {
        if (verbose) {
            cat(sprintf("[WARNING] Warning while sourcing %s: %s\n", abs_path, w$message))
        }
        return(TRUE)
    })
}
