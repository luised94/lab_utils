
# Terminal width detection
get_width <- function() {
    width <- Sys.getenv("COLUMNS")
    if (width != "") {
        as.integer(width)
    } else {
        tryCatch({
            as.integer(system("tput cols", intern = TRUE))
        }, error = function(e) {
            80  # Default width if tput is not available
        })
    }
}

# Enhanced progress tracker
progress_tracker <- function(current, total, description = "", min_interval = 0.5) {
    # Static variables using closure
    if (!exists("last_update", environment(progress_tracker))) {
        environment(progress_tracker)$last_update <- Sys.time() - 1
        environment(progress_tracker)$start_time <- Sys.time()
    }
    
    # Check update interval
    current_time <- Sys.time()
    if (as.numeric(current_time - environment(progress_tracker)$last_update) < min_interval) {
        return(invisible())
    }
    
    # Calculate terminal width and adjust bar size
    term_width <- get_width()
    static_chars <- 30  # Space needed for numbers, percentage, etc.
    desc_length <- if(nchar(description) > 0) nchar(description) + 3 else 0
    bar_width <- max(10, min(50, term_width - static_chars - desc_length))
    
    # Calculate progress metrics
    pct <- sprintf("%3.0f%%", current/total * 100)
    elapsed <- difftime(current_time, environment(progress_tracker)$start_time, units = "secs")
    
    # Calculate ETA
    eta <- if (current > 1) {
        remaining <- (elapsed/current) * (total - current)
        if (remaining < 60) {
            sprintf(" | ETA: %.0fs", remaining)
        } else if (remaining < 3600) {
            sprintf(" | ETA: %.1fm", remaining/60)
        } else {
            sprintf(" | ETA: %.1fh", remaining/3600)
        }
    } else ""
    
    # Create progress bar with blocks
    filled <- round(bar_width * current/total)
    bar <- paste0(
        "[",
        strrep("Û", filled),
        if(current < total) "±" else "Û",
        strrep("°", bar_width - filled - 1),
        "]"
    )
    
    # Construct message
    msg <- sprintf("\r%s %s %d/%d%s %s", 
                  bar, pct, current, total, eta,
                  if(nchar(description) > 0) paste0("| ", description) else "")
    
    # Ensure we don't exceed terminal width
    if (nchar(msg) > term_width) {
        msg <- substr(msg, 1, term_width)
    }
    
    # Update console
    cat(msg)
    if(current == total) {
        cat("\n")
        # Optional: add separator at the end
        # cat(strrep("Ä", min(term_width, 80)), "\n")
    }
    
    # Update last update time
    environment(progress_tracker)$last_update <- current_time
    
    invisible()
}
# Example usage
total <- 50
for(i in 1:total) {
    progress_tracker(i, total, "Processing items")
    Sys.sleep(0.1)  # Simulate work
}
