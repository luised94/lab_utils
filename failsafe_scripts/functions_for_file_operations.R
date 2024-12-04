#' Safe File Writing Function
#' @param data The data to write (or source path for file.copy)
#' @param path Destination file path
#' @param write_fn The writing function to use (write.csv, write.table, file.copy)
#' @param verbose Boolean to control logging output
#' @param ... Additional arguments passed to write_fn
#' @return Boolean indicating success
safe_write_file <- function(data, path, write_fn, verbose = FALSE, ...) {
    if (verbose) {
        cat(sprintf("[INFO] Processing file: %s\n", path))
    }
    
    if (file.exists(path)) {
        if (verbose) {
            cat(sprintf("[WARNING] File already exists: %s\n", path))
        }
        user_input <- readline(prompt="File exists. Overwrite? (y/n): ")
        if (tolower(user_input) == "n") {
            cat(sprintf("[WARNING SKIP] No overwrite. Did not write: %s\n", path))
        } else if (tolower(user_input) == "y") {
        # Execute the provided writing function
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
        } else {
            cat(sprintf("[WARNING SKIP] Option not recognized. Did not write: %s\n", path))
        }
    } else {
        # Execute the provided writing function
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
    
}
