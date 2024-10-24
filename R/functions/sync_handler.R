#' Sync Command Generation Functions
generate_sync_command <- function(directory,
                                config = CONFIG$SYNC) {
    log_info("Generating sync command")
    
    sprintf(
        config$COMMAND,
        directory,
        directory
    )
}
