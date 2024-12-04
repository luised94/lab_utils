#' Safe File Writing Function
#' @param data The data to write (or source path for file.copy)
#' @param path Destination file path
#' @param write_fn The writing function to use (write.csv, write.table, file.copy)
#' @param verbose Boolean to control logging output
#' @param ... Additional arguments passed to write_fn
#' @return Boolean indicating success
safe_write_file <- function(data, path, write_fn, verbose = FALSE, ...) {
    # Input Validation 
    stopifnot(
        "write_fn must be a function" = is.function(write_fn),
        "verbose must be logical" = is.logical(verbose) && length(verbose) == 1,
        "parent directory must exist and be writable" = dir.exists(dirname(path)) && file.access(dirname(path), mode = 2) == 0
    )

    if (verbose) {
        cat(sprintf("[INFO] Processing file: %s\n", path))
    }

    do_write <- function() {
        tryCatch({
            write_fn(data, path, ...)
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
