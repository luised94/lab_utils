library(optparse)
source("~/lab_utils/core_scripts/override_configuration.R")

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
#todo: Add a --track-statistics to measure different stats from the main logic to track certain actions for debugging.
parse_common_arguments <- function(description = "Script description", prog = NULL) {
    option_list <- list(
        make_option(
            "--experiment-id",
            type = "character",
            help = "Experiment ID (format: YYMMDD'Bel', e.g., 241228Bel) or comma separated ID1,ID2,ID3",
            dest = "experiment_id",
            metavar = "ID"
        ),
        make_option(
            "--accept-configuration",
            action = "store_true",
            default = FALSE,
            help = "Accept configuration and continue with execution (default: stop after config display)",
            dest = "accept_configuration"
        ),
        make_option(
            "--log-to-file",
            action = "store_true",
            default = FALSE,
            help = "Enable logging to file (default: FALSE)",
            dest = "log_to_file"
        ),
        make_option(
            "--override",
            type = "character",
            help = paste(
                "Override runtime configuration.",
                "Values:", paste(names(OVERRIDE_PRESETS), collapse = ", "),
                "(default: NULL)"
            ),
            dest = "override",
            metavar = "MODE"
        ),
        make_option(
            "--skip-validation",
            action = "store_true",
            default = TRUE,  # Change to TRUE
            help = "Skip validation for presence of required R packages (default: TRUE). Set to FALSE to force validation.",
            dest = "skip_validation"
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

    # Validate all arguments
    if (is.null(args$experiment_id)) {
        print_help(opt_parser)
        stop("Missing required argument: --experiment-id\nUse -h or --help for usage information")
    }

    # Validate argument types
    if (!is.logical(args$log_to_file)) {
        stop("--log-to-file must be logical value")
    }

    if (!is.logical(args$accept_configuration)) {
        stop("--accept-configuration must be logical value")
    }

    if (!is.null(args$override) && !args$override %in% names(OVERRIDE_PRESETS)) {
        stop(sprintf(
            "--override must be one of: %s",
            paste(names(OVERRIDE_PRESETS), collapse = ", ")
        ))
    }


    # Replace the existing validation section with:
    if (grepl(",", args$experiment_id)) {
        # Split and clean the IDs
        experiment_ids <- stri_split_fixed(args$experiment_id, ",")[[1]]
        experiment_ids <- trimws(experiment_ids)  # Remove any whitespace
        experiment_ids <- experiment_ids[experiment_ids != ""]  # Remove empty elements

        # Check for duplicates
        if (length(unique(experiment_ids)) != length(experiment_ids)) {
            stop("Duplicate experiment IDs detected")
        }

        # Validate format of each ID
        invalid_ids <- experiment_ids[!grepl("^\\d{6}Bel$", experiment_ids, perl = TRUE)]
        if (length(invalid_ids) > 0) {
            stop(sprintf(
                "Invalid experiment-id format(s):\n%s\nExpected format: YYMMDD'Bel'",
                paste(invalid_ids, collapse = ", ")
            ))
        }

        args$experiment_id <- experiment_ids
    } else if (!grepl("^\\d{6}Bel$", args$experiment_id, perl = TRUE) && 
               args$experiment_id != "template") {
        stop(sprintf(
            "Invalid experiment-id format.\nExpected: YYMMDD'Bel' or 'template' or '<exp_id1>,<exp_id2>,...'\nReceived: %s",
            args$experiment_id
        ))
    }

    # Set up experiment directory
    # Modify the directory check section
    if (args$experiment_id[1] == "template") {
        args$experiment_dir <- file.path(Sys.getenv("HOME"), "lab_utils", "core_scripts")
        args$is_template <- TRUE
    } else {
        args$experiment_dir <- sapply(args$experiment_id, function(id) {
            file.path(Sys.getenv("HOME"), "data", id)
        })
        args$is_template <- FALSE

        # Check all directories exist
        missing_dirs <- args$experiment_dir[!dir.exists(args$experiment_dir)]
        if (length(missing_dirs) > 0) {
            stop(sprintf(
                "The following experiment directories were not found:\n%s",
                paste(missing_dirs, collapse = "\n")
            ))
        }
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
    separator <- create_separator()
    format_section <- function(title) {
        paste(
            "\n", separator,
            "\n", title,
            "\n", separator, "\n", sep = ""
        )
    }

    if (!accept_configuration) {
        message(paste(
            format_section("CONFIGURATION CHECKPOINT"),
            "\nThis script is currently configured for configuration confirmation only.",
            "\nNo analysis or processing will be performed at this time.",
            "\n\nTo run the full script and perform the analysis, re-run with the",
            " --accept-configuration flag. For example:\n",
            sprintf("\n    Rscript %s --experiment-id=%s --accept-configuration\n",
                   script_name, experiment_id),
            "\n", separator,
            sep = ""
        ))
        quit(status = 0, save = "no")
    }

    message(paste(
        format_section("CONFIGURATION CONFIRMED"),
        "\nContinuing with script execution...",
        "\n", separator,
        sep = ""
    ))

    invisible(NULL)
}

#' Validate and report required package availability
#' @param packages Character vector of required package names
#' @param verbose Logical, whether to display detailed status messages (defaults to TRUE)
#' @return Invisible list of package status
#' @examples
#' validate_packages(c("rtracklayer", "GenomicRanges", "Gviz"))
check_required_packages <- function(
    packages,
    verbose = TRUE,
    skip_validation = FALSE
) {  # verbose now defaults to TRUE
    if (skip_validation) {
        if (verbose) {
            structured_log_info("Package validation skipped. Use --skip-validation=FALSE to force validation.")
        }
        return(invisible(NULL))
    }

    if (!is.character(packages) || length(packages) == 0) {
        stop("packages must be a non-empty character vector")
    }

    status <- list(
        total = length(packages),
        available = 0,
        missing = character(0)
    )

    if (verbose) { # Use verbose instead of quietly
        message("\nValidating required packages:")
        message(paste(rep("-", 50), collapse = ""))
    }

    for (pkg in packages) {
        if (verbose) cat(sprintf("  %-25s : ", pkg))

        if (requireNamespace(pkg, quietly = TRUE)) {
            status$available <- status$available + 1
            if (verbose) cat("OK\n")
        } else {
            status$missing <- c(status$missing, pkg)
            if (verbose) cat("MISSING\n")
        }
    }

    if (length(status$missing) > 0) {
        missing_pkgs_str <- paste(paste0("   ", status$missing), collapse = "\n")
        install_cmd <- paste0("install.packages(c('", paste(status$missing, collapse = "', '"), "'))")
        bioc_install_cmd <- paste0("BiocManager::install(c('", paste(status$missing, collapse = "', '"), "'))")
        stop(sprintf(
            "\nMissing required packages (%d of %d):\n%s\n\nPlease install them using one of the following commands:\n%s%sIf you are on a cluster, you may need to use your cluster's package installation mechanism.\n", # Corrected sprintf
            length(status$missing),
            status$total,
            missing_pkgs_str,
            if (bioc_install_cmd != ""){paste0("Install from Bioconductor (if applicable):\n    ", bioc_install_cmd, "\n")} else {""},
            paste0("Install from CRAN: \n    ", install_cmd, "\n")
        ))
    }

    if (verbose) {
        message(paste(rep("-", 50), collapse = ""))
        message(sprintf("All required packages available (%d/%d)\n", status$available, status$total))
    }

    invisible(status)
}

#' Setup standard experiment directories
#'
#' This function sets up standard directories for an experiment, including required input directories
#' and an output directory. It validates the existence of input directories (as dependencies) and
#' creates the output directory if it doesn't exist.
#'
#' @param experiment_dir Base experiment directory (must exist).
#' @param output_dir_name Name of the output directory (default: "output").
#' @param required_input_dirs Vector of required input directory names (default: c("fastq", "documentation")).
#' @return A named list of directory paths.
#' @examples
#' \dontrun{
#' # Basic usage
#' dirs <- setup_experiment_dirs("my_experiment")
#' print(dirs)
#'
#' # Custom output directory name
#' dirs <- setup_experiment_dirs("my_experiment", output_dir_name = "results")
#' print(dirs)
#' }
setup_experiment_dirs <- function(experiment_dir, output_dir_name = "output", required_input_dirs = c("fastq", "documentation")) {

    # Input Validation
    stopifnot(
        "experiment_dir must be a character string" = is.character(experiment_dir) && length(experiment_dir) == 1,
        "experiment_dir must be an existing directory" = dir.exists(experiment_dir),
        "output_dir_name must be a non-empty character string" = is.character(output_dir_name) && nchar(output_dir_name) > 0
        # "output_dir_name must be a valid path" = { # More complex validation, may be overkill
        #     tryCatch({
        #         normalizePath(file.path(".", output_dir_name)) # Check if it can be normalized
        #         TRUE
        #     }, error = function(e) FALSE)
        # }
    )

    experiment_dir <- normalizePath(experiment_dir)
    dirs <- list()

    # Check for required input directories
    for (dir_name in required_input_dirs) {
        full_path <- file.path(experiment_dir, dir_name)
        if (!dir.exists(full_path)) {
            stop(sprintf("Required experiment subdirectory '%s' not found: %s", dir_name, full_path))
        }
        dirs[[dir_name]] <- full_path
    }

    # Construct and create the output directory
    output_dir <- file.path(experiment_dir, output_dir_name)
    dirs[[output_dir_name]] <- output_dir # More descriptive name
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    return(dirs)
}

apply_runtime_override <- function(config, preset_name, preset_list) {
    # Current input validation
    if (!is.list(config)) {
        stop("config must be a list")
    }
    if (!is.character(preset_name) || length(preset_name) != 1) {
        stop("preset_name must be a single string")
    }
    if (!is.list(preset_list)) {
        stop("preset_list must be a list of presets")
    }

    # Validate preset exists
    if (!preset_name %in% names(preset_list)) {
        stop(sprintf(
            "Invalid preset name.\nExpected one of: %s\nReceived: %s",
            paste(names(preset_list), collapse = ", "),
            preset_name
        ))
    }

    preset <- preset_list[[preset_name]]

    # Validate preset keys
    invalid_keys <- setdiff(names(preset), names(config))
    if (length(invalid_keys) > 0) {
        stop(sprintf(
            "Preset '%s' contains invalid keys: %s\nValid keys: %s",
            preset_name,
            paste(invalid_keys, collapse = ", "),
            paste(names(config), collapse = ", ")
        ))
    }

    # Validate value types
    for (key in names(preset)) {
        if (!identical(class(preset[[key]]), class(config[[key]]))) {
            stop(sprintf(
                "Type mismatch in preset '%s' for key '%s':\nExpected: %s\nReceived: %s",
                preset_name,
                key,
                class(config[[key]]),
                class(preset[[key]])
            ))
        }
    }

    # Return override info
    list(
        mode = preset_name,
        original = config,
        modified = modifyList(config, preset)
    )
}
