library(futile.logger)
library(stringi)
library(crayon)

is_valid_log_file <- function(log_file) {
    # Check if log_file is a single character string
    is_single_string <- is.character(log_file) && length(log_file) == 1
    # Extract basename without extension
    basename_no_ext <- tools::file_path_sans_ext(basename(log_file))
    # Check if basename exists and is not NA
    has_valid_basename <- !is.na(basename_no_ext)
    # Check if basename contains only alphanumeric characters and underscores
    has_valid_characters <- grepl("^[[:alnum:]_]+$", basename_no_ext)
    # Get the directory of the log file
    log_directory <- dirname(log_file)
    # Check if directory exists
    directory_exists <- dir.exists(log_directory)
    # Check if directory is writable (mode 2 means writable)
    directory_writable <- file.access(log_directory, mode = 2) == 0
    # Combine directory checks
    valid_directory <- directory_exists && directory_writable

    stopifnot(
        "log_file must be a single character string" = is_single_string,
        "Basename required" = has_valid_basename,
        "Basename must contain only letters, numbers and underscore" = has_valid_characters,
        "Directory must exist and be writable" = valid_directory
    )
}

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

create_separator <- function(
    char = "=",
    width = get_width()
) {
  stri_dup(char, width)
}

print_debug_info <- function(
    debug_info,
    indent_char = "  ",
    title_char = "#",
    log_file = NULL
) {
    # Input validation
    if (!is.list(debug_info) || length(debug_info) == 0) {
        flog.error("debug_info must be a non-empty list")
        stop("debug_info must be a non-empty list")
    }

    if (any(sapply(names(debug_info), nchar) == 0)) {
        flog.error("All elements in debug_info must be named")
        stop("All elements in debug_info must be named")
    }

    # Log the start of the debug information print process
    flog.info("Starting to print debug information.")

    # Get terminal width
    terminal_width <- get_width()

    # Prepare output lines
    output_lines <- c()

    # Process each item in the list, excluding "title"
    for (key in names(debug_info)) {
        if (key == "title") next  # Skip processing for the title

        value <- debug_info[[key]]
        #indent_level <- stri_count_fixed(key, "^.")
        indent_level <- nchar(gsub("[^.].*$", "", key))

        indent <- stri_dup(indent_char, indent_level)
        clean_key <- stri_replace_all_regex(key, "^[.]+", "")

        if (is.null(value)) {
            line <- sprintf("%s%s:", indent, clean_key)
        } else {
            line <- sprintf("%s%s: %s", indent, clean_key, as.character(value))
        }

        output_lines <- c(output_lines, line)
    }

    # Calculate longest line length for box drawing
    #max_line_length <- max(nchar(normalized_matrix_rows))
    # Calculate longest line length for box drawing
    max_line_length <- max(nchar(output_lines))

    # Create top and bottom separators based on longest line
    separator_line <- create_separator("=", max_line_length + nchar(title_char) + 2)

    # Apply padding and box effect to each processed line
    padded_output_lines <- lapply(output_lines, function(line) {
        sprintf("%s-", stri_pad_right(line, max_line_length + nchar(title_char) + 1))
    })


    # Combine everything into final output with box effect
    final_output <- c(
        separator_line,
        sprintf("%s %s-", title_char, stri_pad_right(debug_info$title, max_line_length)),
        separator_line,
        padded_output_lines,
        separator_line
    )

    # Convert final_output to character vector for logging
    final_output_character <- as.character(final_output)


    # Log the final output directly to a file without using logging messages
    if (!is.null(log_file) && is_valid_log_file(log_file)) {
        tryCatch({
          writeLines(final_output_character, con = log_file)
        }, error = function(e) {
          stop(paste("Error writing to log file:", e$message))
        })
    }

    ###
    # Color Output (only if ANSI is supported)
    ###
    ansi_supported <- crayon::has_color()

    colorize_line <- function(line, color_function) {
        if (ansi_supported) {
            color_function(line)
        } else {
            line
        }
    }

    # Count lines made up entirely of '=' for color coding in console
    equal_lines_count <- sum(grepl("^=+$", final_output))

    for (i in seq_along(final_output)) {
        line <- final_output[i]

        if (grepl("^=+$", line)) {
            cat(colorize_line(line, blue), "\n")   # Lines composed entirely of '=' in blue
        } else if (grepl("^#", line)) {
            cat(colorize_line(line, red), "\n")   # Title line in red
        } else if (i < equal_lines_count + 1) {
            cat(colorize_line(line, blue), "\n")  # All other lines until the last one in blue
        } else if (i < length(final_output)) {
            cat(colorize_line(line, green), "\n") # Each processed line in green
        } else {
            cat(colorize_line(line, blue), "\n")  # Bottom separator in blue
        }
    }

    # Log completion
    flog.info("Finished printing debug information.")
}

print_config_settings <- function(
    config_list,
    title = NULL,
    verbose = TRUE,
    log_file = NULL
) {

    # Input validation
    stopifnot(
        "config_list must be a list" = is.list(config_list),
        "config_list cannot be empty" = length(config_list) > 0,
        "verbose must be logical" = is.logical(verbose) && length(verbose) == 1
    )

    if (!verbose) {
        return(invisible(config_list))
    }

    # Format title
    if (is.null(title)) {
        title <- deparse(substitute(config_list))
    }

    # Get terminal width
    terminal_width <- get_width()

    # Create a normalized header with '=' padding
    header_format <- "=== %s "
    header <- sprintf(header_format, title)
    # Calculate the remaining width for '=' padding
    remaining_width <- max(0, terminal_width - nchar(header))
    # Append '=' characters to fill the line
    header <- paste0(header, paste(rep("=", remaining_width), collapse = ""))
    # Ensure the header does not exceed the terminal width
    header <- substr(header, 1, terminal_width)
    separator_line <- create_separator("=")
    cat(sprintf("\n%s\n", header))

    # Process and format each setting
    settings_lines <- lapply(names(config_list), function(setting) {
        value <- config_list[[setting]]

        formatted_value <- if (is.logical(value)) {
            if (value) "YES" else "NO"
        } else if (is.null(value)) {
            "NULL"
        } else if (is.list(value)) {
            "LIST"
        } else if (is.vector(value) & length(value) > 0) {
            escaped_newlines <- gsub("\n", "\\\\n", value)
            paste0("[", paste(escaped_newlines, collapse = ", "), "]")
        } else {
            stop(sprintf(
                "Unsupported configuration value type for setting '%s':\n  Type: %s\n  Value: %s",
                setting,
                class(value)[1],
                deparse(value)
            ))
        }

        max_length <- 100
        if (nchar(formatted_value) > max_length) {
            formatted_value <- paste0(substr(formatted_value, 1, max_length - 3), "...")
        }

        sprintf("%s: %s", gsub("_", " ", toupper(setting)), formatted_value)
    })

    max_item_width <- max(nchar(settings_lines)) + 2
    num_columns <- max(1, floor(terminal_width / max_item_width))
    num_rows <- ceiling(length(settings_lines) / num_columns)

    settings_matrix <- matrix("", nrow = num_rows, ncol = num_columns)
    for (i in seq_along(settings_lines)) {
        row_index <- ((i - 1) %% num_rows) + 1
        col_index <- ((i - 1) %/% num_rows) + 1
        settings_matrix[row_index, col_index] <- settings_lines[[i]]
    }

    normalized_matrix_rows <- apply(settings_matrix, 1, function(row) {
        padded_row_items <- sapply(seq_along(row), function(j) {
            sprintf("%-*s", max_item_width, row[j])
        })
        paste(padded_row_items, collapse = " ")
    })


    if (!is.null(log_file) && is_valid_log_file(log_file)) {
        tryCatch({
          writeLines(normalized_matrix_rows, con = log_file)
        }, error = function(e) {
          stop(paste("Error writing to log file:", e$message))
        })
    }

    cat(paste(normalized_matrix_rows, collapse = "\n"), "\n")
    cat(separator_line, "\n")

    invisible(config_list)
}

structured_log_info <- function(
  message,
  step = NULL,
  skip_functions = c("eval", "envir", "source", "withVisible", "structured_log_info"),
  log_level = "INFO"
) {
    # --------------------------------------------------------------------------
    # 1. Input Validation
    # --------------------------------------------------------------------------
    if (!is.character(message) || length(message) != 1) {
        stop("`message` must be a single character string.")
    }
    if (!is.null(step) && (!is.character(step) || length(step) != 1)) {
        stop("`step`, if provided, must be a single character string.")
    }
    if (!is.character(skip_functions)) {
        stop("`skip_functions` must be a character vector of function names to ignore.")
    }
    if (!is.character(log_level) || length(log_level) != 1) {
        stop("`log_level` must be a single character string (e.g., \"INFO\", \"DEBUG\").")
    }
    # --------------------------------------------------------------------------
    # 2. Capture & Prune the Call Stack
    # --------------------------------------------------------------------------
    calls <- sys.calls()
    # Pull out the function names from each call
    call_names <- sapply(calls, function(a_call) {
        func_part <- a_call[[1]]
        as.character(func_part)
    })
    # Reverse the vector so the most recent call is first
    reversed_names <- rev(call_names)
    # Identify the first user-defined call by skipping known base / wrapper calls
    idx_user_call <- which(!reversed_names %in% skip_functions)[1]
    function_context <- NULL
    if (!is.na(idx_user_call)) {
        # Convert index in reversed array to index in original array
        actual_index <- length(reversed_names) - idx_user_call + 1
        function_context <- call_names[actual_index]
    }
    if (is.null(function_context)) {
        # Because we filtered out everything or we never had a user function
        # Distinguish whether it's top-level or truly unknown
        # For instance, test if "structured_log_info" was the only user-level call:
        #   - If so, you're *probably* at script level.
        #   - Otherwise, you can still call it unknown_function or handle it differently.
        # Simple approach: always call it "script_level"
        function_context <- "script_level"
    }
    # --------------------------------------------------------------------------
    # 3. Construct the Log Message
    # --------------------------------------------------------------------------
    timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
    # Build a logfmt-like message for easy grep/sed/awk usage
    # If you prefer JSON, you could build a JSON string here.

    structured_message <- paste0(
        "level=", log_level, " ",
        "time=", timestamp, " ",
        "function=", function_context, " ",
        if (!is.null(step)) paste0("step=", step, " ") else "",
        "message=\"", message, "\""
    )
    # --------------------------------------------------------------------------
    # 4. Log Using the Appropriate Level
    # --------------------------------------------------------------------------
    # We can dispatch based on `log_level` to use flog.debug(), flog.info(), etc.,
    # or we can always use flog.info and rely on the textual level.  
    # Here we match the actual futile.logger function to the parameter.
    switch(
        toupper(log_level),
        "DEBUG" = flog.debug(structured_message),
        "INFO"  = flog.info(structured_message),
        "WARN"  = flog.warn(structured_message),
        "ERROR" = flog.error(structured_message),
        # Default if unrecognized log level
        flog.info(structured_message)
    )
}

log_system_diagnostics <- function() {
    # Get memory info based on OS
    mem_info <- if (.Platform$OS.type == "windows") {
        memory.limit()
    } else {
        tryCatch({
            as.numeric(system("free -g | awk '/^Mem:/{print $2}'", intern = TRUE))
        }, error = function(e) "not_available")
    }

    list(
        title = "System Diagnostics",
        "system.hostname" = Sys.info()["nodename"],
        "r.version" = R.version.string,
        "r.platform" = R.version$platform,
        "os.type" = .Platform$OS.type,
        "current.directory" = getwd(),
        "session.timezone" = Sys.timezone(),
        "memory" = mem_info,
        "cpu.cores" = parallel::detectCores(),
        "user" = Sys.getenv("USER"),
        "disk_free" = tryCatch({
            system("df -h . | tail -1", intern = TRUE)
        }, error = function(e) "not_available"),
        "git.branch" = tryCatch({
            system("git rev-parse --abbrev-ref HEAD", intern = TRUE)
        }, error = function(e) "git_not_available"),
        "git.commit" = tryCatch({
            system("git rev-parse HEAD", intern = TRUE)
        }, error = function(e) "git_not_available")
    )
}

log_session_info <- function() {
    all_objects <- ls(envir = globalenv())
    
    # Better named object summary
    object_summary <- list()
    for (obj_name in all_objects) {
        obj_value <- get(obj_name, envir = globalenv())
        object_summary[[paste0("object.", obj_name)]] <- sprintf(
            "%s (%s)", 
            class(obj_value)[1],
            format(object.size(obj_value), units = "auto")
        )
    }
    
    list(
        title = "Session Summary",
        "total_objects" = length(all_objects),
        "loaded_packages" = names(sessionInfo()$loadedOnly),
        "functions" = all_objects[sapply(all_objects, 
            function(x) is.function(get(x, envir = globalenv())))],
        "session_time" = Sys.time()
    ) |> 
        c(object_summary)  # Append object summaries at same level
}

get_script_name <- function() {
    tryCatch({
        cmd_args <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("--file=", cmd_args, value = TRUE)
        if (length(file_arg) > 0) {
            basename(sub("--file=", "", file_arg))
        } else {
            ifelse(interactive(), "interactive_session", "batch_session")
        }
    }, error = function(e) "unknown_script")
}
