#!/usr/bin/env Rscript
# Terminal width detection
#get_width <- function() {
#    width <- Sys.getenv("COLUMNS")
#    if (width != "") {
#        as.integer(width)
#    } else {
#        tryCatch({
#            as.integer(system("tput cols", intern = TRUE))
#        }, error = function(e) {
#            80  # Default width if tput is not available
#        })
#    }
#}

# Enhanced progress tracker
#progress_tracker <- function(current, total, description = "", min_interval = 0.5) {
#    # Static variables using closure
#    if (!exists("last_update", environment(progress_tracker))) {
#        environment(progress_tracker)$last_update <- Sys.time() - 1
#        environment(progress_tracker)$start_time <- Sys.time()
#    }
#    
#    # Check update interval
#    current_time <- Sys.time()
#    if (as.numeric(current_time - environment(progress_tracker)$last_update) < min_interval) {
#        return(invisible())
#    }
#    
#    # Calculate terminal width and adjust bar size
#    term_width <- get_width()
#    static_chars <- 30  # Space needed for numbers, percentage, etc.
#    desc_length <- if(nchar(description) > 0) nchar(description) + 3 else 0
#    bar_width <- max(10, min(50, term_width - static_chars - desc_length))
#    
#    # Calculate progress metrics
#    pct <- sprintf("%3.0f%%", current/total * 100)
#    elapsed <- difftime(current_time, environment(progress_tracker)$start_time, units = "secs")
#    
#    # Calculate ETA
#    eta <- if (current > 1) {
#        remaining <- (elapsed/current) * (total - current)
#        if (remaining < 60) {
#            sprintf(" | ETA: %.0fs", remaining)
#        } else if (remaining < 3600) {
#            sprintf(" | ETA: %.1fm", remaining/60)
#        } else {
#            sprintf(" | ETA: %.1fh", remaining/3600)
#        }
#    } else ""
#    
#    # Create progress bar with blocks
#    filled <- round(bar_width * current/total)
#    bar <- paste0(
#        "[",
#        strrep("Û", filled),
#        if(current < total) "±" else "Û",
#        strrep("°", bar_width - filled - 1),
#        "]"
#    )
#    
#    # Construct message
#    msg <- sprintf("\r%s %s %d/%d%s %s", 
#                  bar, pct, current, total, eta,
#                  if(nchar(description) > 0) paste0("| ", description) else "")
#    
#    # Ensure we don't exceed terminal width
#    if (nchar(msg) > term_width) {
#        msg <- substr(msg, 1, term_width)
#    }
#    
#    # Update console
#    cat(msg)
#    if(current == total) {
#        cat("\n")
#        # Optional: add separator at the end
#        # cat(strrep("Ä", min(term_width, 80)), "\n")
#    }
#    
#    # Update last update time
#    environment(progress_tracker)$last_update <- current_time
#    
#    invisible()
#}
# Terminal width detection
get_terminal_width <- function() {
    # Try to get terminal width from environment variable
    terminal_width <- Sys.getenv("COLUMNS")
    
    if (terminal_width != "") {
        return(as.integer(terminal_width))
    }
    
    # Fallback to using `tput cols` command
    tryCatch({
        as.integer(system("tput cols", intern = TRUE))
    }, error = function(error) {
        # Default width if tput is not available or fails
        80L
    })
}

# Enhanced progress tracker
progress_tracker <- function(current, total, description = "", min_update_interval = 0.5) {
    # Static variables to track progress state
    if (!exists("last_update_time", envir = environment(progress_tracker))) {
        assign("last_update_time", Sys.time() - 1, envir = environment(progress_tracker))
        assign("start_time", Sys.time(), envir = environment(progress_tracker))
    }
    
    # Ensure current does not exceed total
    current <- min(current, total)
    
    # Check if enough time has passed since the last update
    current_time <- Sys.time()
    time_since_last_update <- as.numeric(current_time - get("last_update_time", envir = environment(progress_tracker)))
    if (time_since_last_update < min_update_interval && current < total) {
        return(invisible())
    }
    
    # Calculate terminal width and adjust progress bar size
    terminal_width <- get_terminal_width()
    static_elements_width <- 30L  # Space for numbers, percentage, etc.
    description_width <- if (nchar(description) > 0) nchar(description) + 3L else 0L
    progress_bar_width <- max(10L, min(50L, terminal_width - static_elements_width - description_width))
    
    # Calculate progress metrics
    progress_percentage <- sprintf("%3.0f%%", current / total * 100)
    elapsed_time <- difftime(current_time, get("start_time", envir = environment(progress_tracker)), units = "secs")
    
    # Calculate Estimated Time of Arrival (ETA)
    estimated_time_remaining <- if (current > 1) {
        remaining_time <- (elapsed_time / current) * (total - current)
        if (remaining_time < 60) {
            sprintf(" | ETA: %.0fs", remaining_time)
        } else if (remaining_time < 3600) {
            sprintf(" | ETA: %.1fm", remaining_time / 60)
        } else {
            sprintf(" | ETA: %.1fh", remaining_time / 3600)
        }
    } else ""
    
    # Create progress bar with blocks (ensure no negative values)
    filled_blocks <- round(progress_bar_width * current / total)
    empty_blocks <- max(0, progress_bar_width - filled_blocks - 1L)
    
    progress_bar <- paste0(
        "[",
        strrep("#", filled_blocks),  # Use '#' for filled blocks
        if (current < total) "=" else "#",  # Use '=' for the current block
        strrep("-", empty_blocks),  # Use '-' for empty blocks (clamped to 0)
        "]"
    )
    
    # Construct progress message
    progress_message <- sprintf("\r%s %s %d/%d%s %s", 
                                progress_bar, progress_percentage, current, total, estimated_time_remaining,
                                if (nchar(description) > 0) paste0("| ", description) else "")
    
    # Ensure the message does not exceed terminal width
    if (nchar(progress_message) > terminal_width) {
        progress_message <- substr(progress_message, 1L, terminal_width)
    }
    
    # Update console with progress message
    cat(progress_message)
    
    # Finalize progress bar if current >= total
    if (current >= total) {
        cat("\n")  # Move to a new line after completion
    }
    
    # Update last update time
    assign("last_update_time", current_time, envir = environment(progress_tracker))
    
    invisible()
}
# Example usage
total <- 50
for(i in 1:total) {
    progress_tracker(i, total, "Processing items")
    Sys.sleep(0.1)  # Simulate work
}
total_items <- 50
batch_size <- 5

for (i in seq(1, total_items, by = batch_size)) {
    current <- min(i + batch_size - 1, total_items)  # Ensure current <= total
    progress_tracker(current, total_items, description = "Processing items")
    cat("Processing item", i, "with some additional details...\n")
    Sys.sleep(0.5)  # Simulate work being done
}
