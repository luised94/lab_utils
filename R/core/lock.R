# R/core/lock.R

#' Acquire File Lock
#' @param file_path Character Path to file
#' @param timeout Integer Seconds to wait
#' @return Logical TRUE if lock acquired
#' @export
acquire_lock <- function(
    file_path,
    timeout = CORE_CONFIG$SYSTEM$LOCK_TIMEOUT
) {
    lock_file <- paste0(file_path, ".lock")
    start_time <- Sys.time()
    
    while (difftime(Sys.time(), start_time, units = "secs") < timeout) {
        if (!file.exists(lock_file)) {
            # Create lock file
            tryCatch({
                file.create(lock_file)
                writeLines(as.character(Sys.getpid()), lock_file)
                return(TRUE)
            }, error = function(e) NULL)
        }
        Sys.sleep(0.1)
    }
    return(FALSE)
}

#' Release File Lock
#' @param file_path Character Path to file
#' @return Logical TRUE if successful
release_lock <- function(file_path) {
    lock_file <- paste0(file_path, ".lock")
    if (file.exists(lock_file)) {
        unlink(lock_file)
    }
    return(!file.exists(lock_file))
}
