library(optparse)

validate_configs <- function(required_configs) {
    # 1. Input Validation
    stopifnot(
        "required_configs must be a character vector" = 
            is.character(required_configs),
        "required_configs cannot be empty" = 
            length(required_configs) > 0,
        "required_configs cannot contain NA values" = 
            !any(is.na(required_configs)),
        "required_configs cannot contain empty strings" = 
            !any(required_configs == ""),
        "required_configs must be unique" = 
            !any(duplicated(required_configs))
    )

    # 2. Ensure each item ends with "_CONFIG"
    invalid_pattern <- grep(".*_CONFIG$", required_configs, invert = TRUE, value = TRUE)
    if (length(invalid_pattern) > 0) {
        stop(
          "All required_configs must end with '_CONFIG'. Invalid entries: ",
          paste(invalid_pattern, collapse = ", ")
        )
    }
    
    # 3. Verify each one actually exists in the environment
    missing <- required_configs[!sapply(required_configs, exists)]
    if (length(missing) > 0) {
        stop("Missing required configs: ", paste(missing, collapse = ", "))
    }
    
    # If we get here, all checks have passed
    invisible(TRUE)
}

#' Parse common command line arguments for experiment scripts
#' @param description Character string describing the script's purpose
#' @param prog Character string of program name (optional)
#' @return List of validated arguments
#' @examples
#' args <- parse_common_arguments("Align sequencing reads to reference genome")
parse_common_arguments <- function(description = "Script description", prog = NULL) {
    option_list <- list(
        make_option(
            "--experiment-id",
            type = "character",
            help = "Experiment ID (format: YYMMDD'Bel', e.g., 241228Bel)",
            dest = "experiment_id",
            metavar = "ID"
        ),
        make_option(
            "--accept-configuration",
            action = "store_true",
            default = FALSE,
            help = "Accept configuration and continue with execution (default: stop after config display)",
            dest = "accept_configuration"
        )
    )

    # Create parser with optional program name
    opt_parser <- OptionParser(
        option_list = option_list,
        description = description,
        prog = prog
    )

    # Parse arguments with enhanced error handling
    args <- tryCatch(
        parse_args(opt_parser),
        error = function(e) {
            cat("\nERROR: Argument parsing failed\n", file = stderr())
            cat(as.character(e), "\n\n", file = stderr())
            print_help(opt_parser)
            quit(status = 1, save = "no")
        }
    )

    # Validate required argument
    if (is.null(args$experiment_id)) {
        print_help(opt_parser)
        stop("Missing required argument: --experiment-id")
    }

    # Validate experiment ID format
    if (!grepl("^\\d{6}Bel$", args$experiment_id, perl = TRUE)) {
        stop(sprintf(
            "Invalid experiment-id format.\nExpected: YYMMDD'Bel'\nReceived: %s",
            args$experiment_id
        ))
    }

    return(args)
}

#' Display configuration checkpoint status and handle execution flow
#' @param accept_configuration Logical indicating whether to proceed with execution
#' @param experiment_id Character string of experiment ID (for command example)
#' @param script_name Optional character string of script name (defaults to current script)
#' @return Invisible NULL, exits if configuration not accepted
#' @examples
#' handle_configuration_checkpoint(args$accept_configuration, args$experiment_id)
handle_configuration_checkpoint <- function(
    accept_configuration,
    experiment_id,
    script_name = basename(commandArgs(trailingOnly = FALSE)[4])
) {
    # Input validation
    if (!is.logical(accept_configuration) || length(accept_configuration) != 1 || is.na(accept_configuration)) {
        stop("accept_configuration must be a single logical value (TRUE/FALSE)")
    }
    if (!is.character(experiment_id) || length(experiment_id) != 1 || is.na(experiment_id)) {
        stop("experiment_id must be a single character string")
    }
    
    # Common formatting elements
    separator <- strrep("=", 80)
    format_section <- function(title) {
        paste(
            "\n", separator,
            "\n", title,
            "\n", separator, "\n"
        )
    }
    
    if (!accept_configuration) {
        message(paste(
            format_section("CONFIGURATION CHECKPOINT"),
            "\nThis script is currently configured for configuration confirmation only.",
            "\nNo analysis or processing will be performed at this time.",
            "\n\nTo run the full script and perform the analysis, re-run with the",
            "--accept-configuration flag. For example:\n",
            sprintf("\n    Rscript %s --experiment-id=%s --accept-configuration\n",
                   script_name, experiment_id),
            "\n", separator,
            sep = ""
        ))
        quit(status = 0)
    }
    
    message(paste(
        format_section("CONFIGURATION CONFIRMED"),
        "\nContinuing with script execution...",
        "\n", separator,
        sep = ""
    ))
    
    invisible(NULL)
}
