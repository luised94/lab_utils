#' Safe File Writing Function
#' @param data The data to write (or source path for file.copy)
#' @param path Destination file path
#' @param write_fn The writing function to use. Supported functions:
#'                 - write.csv(x, file, ...)
#'                 - write.table(x, file, ...)
#'                 - file.copy(from, to, ...)
#' @param verbose Boolean to control logging output
#' @param ... Additional arguments passed to write_fn
#' @return Boolean indicating success
safe_write_file <- function(data, path, write_fn, verbose = FALSE, ...) {
    # Input validation 
    stopifnot(
        "write_fn must be a function" = is.function(write_fn),
        "verbose must be logical" = is.logical(verbose) && length(verbose) == 1,
        "parent directory must exist and be writable" = dir.exists(dirname(path)) && file.access(dirname(path), mode = 2) == 0
    )

    # Create wrapper for supported functions
    wrapped_write <- function(data, path, ...) {
        if (identical(write_fn, write.csv) || identical(write_fn, write.table)) {
            write_fn(x = data, file = path, ...)
        } else if (identical(write_fn, file.copy)) {
            write_fn(from = data, to = path, ...)
        } else {
            stop(sprintf(
                "Unsupported writing function: %s. Supported functions are: write.csv, write.table, file.copy",
                deparse(substitute(write_fn))
            ))
        }
    }
    if (verbose) {
        cat(sprintf("[INFO] Processing file: %s\n", path))
    }

    do_write <- function() {
        tryCatch({
            wrapped_write(data, path, ...)
            if (verbose) {
                cat(sprintf("[WROTE] File to: %s\n", path))
            }
            return(TRUE)
        }, error = function(e) {
            cat(sprintf("[ERROR] Failed to write file: %s\n", e$message))
            return(FALSE)
        })
    }

    # If file doesn't exist, write immediately
    if (!file.exists(path)) {
        return(do_write())
    }

    # Handle existing file
    if (verbose) {
        cat(sprintf("[WARNING] File already exists: %s\n", path))
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
    cat(sprintf(message, path))
    return(FALSE)
}
