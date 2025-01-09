
#' Script Information Functions
get_script_info <- function() {
    log_info("Getting script information")
    
    script_path <- get_script_path()
    
    list(
        path = script_path,
        dir = dirname(script_path),
        name = get_script_name(script_path),
        config = get_script_config(script_path)
    )
}

get_script_path <- function() {
    # Try command line args first
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", cmd_args, value = TRUE)
    
    if (length(file_arg) > 0) {
        return(normalizePath(sub("--file=", "", file_arg)))
    }
    
    # Try sys.frames for sourced scripts
    if (!is.null(sys.frames()[[1]]$ofile)) {
        return(normalizePath(sys.frames()[[1]]$ofile))
    }
    
    message("Unable to determine script path")
    return("interactive")
}

get_script_name <- function(script_path) {
    if (script_path == "interactive") {
        return("interactive")
    }
    
    tools::file_path_sans_ext(basename(script_path))
}

get_script_config <- function(script_path) {
    script_name <- get_script_name(script_path)
    
    if (script_name == "interactive") {
        return(NULL)
    }
    
    CONFIG$SCRIPTS[[script_name]]
}

#' Add to existing configurations
CONFIG <- list(
    PLOT = list(
        COLORS = list(
            DEFAULT = "#fd0036",
            CONTROL = "#808080"
        ),
        TRACK_TYPE = "l",
        DATE_FORMAT = "%Y%m%d%H%M%S"
    ),
    
    PATHS = list(
        PLOTS = "plots",
        BIGWIG = "bigwig"
    ),
    
    PATTERNS = list(
        COMPARISON_PREFIX = "^comp_",
        DEFAULT_BIGWIG = "S288C_log2ratio",
        TIMECOURSE = "comp_timecourse1108"
    ),
    
    CATEGORIES = list(
        DEFAULT = c(
            "strain_source",
            "rescue_allele",
            "mcm_tag",
            "antibody",
            "timepoint_after_release"
        )
    )
)
#' Add to existing configurations
CONFIG <- list(
    PACKAGES = list(
        REQUIRED = c(
            "QuasR",
            "GenomicAlignments",
            "Gviz",
            "rtracklayer",
            "ShortRead",
            "tidyverse",
            "gtools"
        )
    ),
    
    DEFAULTS = list(
        CHROMOSOME = 10,
        GENOME_DIR = "REFGENS",
        GENOME_PATTERN = "S288C_refgenome.fna",
        FEATURE_PATTERN = "eaton_peaks",
        BIGWIG_PATTERN = "_bamcomp.bw"
    ),
    
    SYNC = list(
        COMMAND = "rsync -nav username@domain:~/data/%s/plots/* /local/dir/%s/plots/",
        REMOTE_USER = "username",
        REMOTE_DOMAIN = "domain"
    )
)
#' Add to existing visualization configurations
CONFIG <- list(
    VISUALIZATION = list(
        COLORS = list(
            SAMPLE = "#E41A1C",
            CONTROL = "#377EB8",
            HIGHLIGHT = list(
                FILL = "#FFE3E6",
                BORDER = "#FF0000",
                ALPHA = 0.3
            )
        ),
        TRACK_TYPES = list(
            LINE = "l",
            HIGHLIGHT = "highlight"
        ),
        LIMITS = list(
            Y_MIN = 0,
            Y_MAX = 100000
        ),
        OUTPUT = list(
            FORMAT = "svg",
            DATE_FORMAT = "%Y%m%d%H%M%S"
        )
    ),
    
    PATHS = list(
        SUBDIRS = list(
            PLOTS = "plots",
            BIGWIG = "bigwig"
        )
    )
)

#!/usr/bin/env Rscript

#' Validation Utilities
#' @description Core validation functions for project
#' @export

#' Check Project Configuration
#' @param throw_error Logical Whether to throw error if config missing
#' @return Logical TRUE if config exists
check_project_config <- function(throw_error = TRUE) {
    config_exists <- exists("PROJECT_CONFIG", envir = .GlobalEnv)
    if (!config_exists && throw_error) {
        stop("PROJECT_CONFIG not loaded. Source project_config.R first.")
    }
    invisible(config_exists)
}

#' Validate Function Arguments
#' @param args List Arguments to validate
#' @param spec List Argument specifications
#' @return List Validated arguments
validate_function_args <- function(args, spec) {
    check_project_config()
    
    for (arg_name in names(spec)) {
        arg_spec <- spec[[arg_name]]
        arg_value <- args[[arg_name]]
        
        # Check required arguments
        if (is.null(arg_value)) {
            if (arg_spec$required) {
                stop(sprintf("Required argument missing: %s", arg_name))
            }
            args[[arg_name]] <- arg_spec$default
            next
        }
        
        # Type checking
        if (!inherits(arg_value, arg_spec$type)) {
            stop(sprintf(
                "Invalid type for %s: expected %s, got %s",
                arg_name, arg_spec$type, class(arg_value)[1]
            ))
        }
        
        # Custom validation
        if (!is.null(arg_spec$validation)) {
            if (!arg_spec$validation(arg_value)) {
                stop(arg_spec$error_message)
            }
        }
    }
    
    invisible(args)
}

#' Validate NGS File
#' @param file_path Character Path to file
#' @param type Character File type
#' @param config List Configuration settings
#' @return Logical TRUE if valid
validate_ngs_file <- function(
    file_path,
    type = c("bam", "fastq", "bigwig", "bed", "narrowpeak", "motif"),
    config = PROJECT_CONFIG
) {
    type <- match.arg(type)
    
    # Basic existence check
    if (!file.exists(file_path)) {
        log_error(sprintf("File not found: %s", file_path))
        return(FALSE)
    }
    
    # Extension check
    pattern <- config$FILE_TYPES$NGS$EXTENSIONS[[toupper(type)]]
    if (!grepl(pattern, file_path)) {
        log_error(sprintf("Invalid file extension for %s: %s", type, file_path))
        return(FALSE)
    }
    
    # Index check if required
    if (type %in% names(config$FILE_TYPES$NGS$REQUIRED_INDEX)) {
        index_pattern <- config$FILE_TYPES$NGS$REQUIRED_INDEX[[toupper(type)]]
        index_file <- sub(pattern, index_pattern, file_path)
        if (!file.exists(index_file)) {
            log_warning(sprintf("Index file missing: %s", index_file))
            return(FALSE)
        }
    }
    
    TRUE
}

#' Validate BAM File
#' @param file_path Character Path to BAM file
#' @return Logical TRUE if valid
validate_bam_file <- function(file_path) {
    # Basic checks
    if (!grepl("\\.bam$", file_path)) {
        log_error("Not a BAM file:", file_path)
        return(FALSE)
    }
    
    # Check index
    if (!file.exists(paste0(file_path, ".bai"))) {
        log_warning("BAM index missing:", file_path)
        return(FALSE)
    }
    
    TRUE
}

#' Example Usage Configuration
VALIDATION_CONFIG <- list(
    TYPES = c(
        "bam", "fastq", "bigwig", 
        "bed", "bedgraph", "narrowPeak"
    ),
    LIMITS = list(
        CHROMOSOME = c(min = 1, max = 16),
        READ_LENGTH = c(min = 20, max = 150)
    )
)

#' Package Loading Functions
load_required_packages <- function(packages = CONFIG$REQUIRED_PACKAGES) {
    log_info("Loading required packages")
    
    suppressPackageStartupMessages({
        results <- sapply(packages, require, character.only = TRUE)
    })
    
    if (!all(results)) {
        missing_packages <- packages[!results]
        log_error("Failed to load packages:", paste(missing_packages, collapse = ", "))
        stop("Missing required packages")
    }
    
    return(TRUE)
}

#' Input Validation Functions
validate_input <- function(args) {
    log_info("Validating input arguments")
    
    if (length(args) != 1) {
        log_error("Invalid number of arguments")
        stop("Usage: Rscript script.R <directory_path>")
    }
    
    directory_path <- file.path(CONFIG$PATHS$BASE_DIR, args[1])
    
    if (!dir.exists(directory_path)) {
        log_error("Directory not found:", directory_path)
        stop("Invalid directory")
    }
    
    return(directory_path)
}

#' Track Creation and Management Functions
create_genome_track <- function(chromosome) {
    log_info("Creating genome axis track")
    
    GenomeAxisTrack(
        name = sprintf("Chr %s Axis", chromosome)
    )
}

create_data_track <- function(bigwig_data,
                            name,
                            chromosome,
                            config = CONFIG$PLOT) {
    log_info("Creating data track:", name)
    
    DataTrack(
        bigwig_data,
        type = config$TRACK_TYPE,
        name = name,
        col = config$COLORS$DEFAULT,
        chromosome = chromosome
    )
}

create_control_track <- function(bigwig_data,
                               chromosome,
                               config = CONFIG$PLOT) {
    log_info("Creating control track")
    
    DataTrack(
        bigwig_data,
        type = config$TRACK_TYPE,
        name = "Input",
        col = config$COLORS$CONTROL,
        chromosome = chromosome
    )
}

source("~/lab_utils/R/functions/001_logging.R")
source("~/lab_utils/R/init.R")
library(assertthat)

sort_columns <- function(df, column_sort_order) {
    # Input validation
    if(!setequal(column_sort_order, colnames(df))) {
        log_error("Column must be sorted using all columns to ensure proper order.")
        stop("Modify column_sort_order variable appropriately.")
    }
    # Create a list of sorting criteria
    sort_criteria <- lapply(column_sort_order, function(col) df[[col]])
    # Sort the dataframe
    df[do.call(order, sort_criteria), ]
}

add_sample_names_to_table <- function(df) {
     df$full_name <- apply(df, 1, paste, collapse = "_")
     colnames_sans_fullname <- !grepl("full_name", colnames(df))
     df_sans_fullname <- df[, colnames_sans_fullname]
     df$short_name <- apply(df_sans_fullname, 1, function(row) paste0(substr(row, 1, 1), collapse = ""))
     return(df)
 }

verify_expected_number_of_samples <- function(df, expected_number_of_samples) {
    # Input validation
    assert_that(is.numeric(expected_number_of_samples), "expected_number_of_samples must be a numeric.")
    if(!nrow(df) == expected_number_of_samples) {
        log_error("Combinations grid does not contain the expected_number_of_samples.")
        log_info("Elements of sample_table:\n")
        print(df)
        log_info("Dimensions of sample_table:\n")
        print(dim(df))
        log_info("Breakdown by antibody:\n")
        print(table(df$antibody))
        stop("Update the impossible_settings and experiment_conditions variable.\nEnsure that it matches configuration_settings$expected_number_of_samples.")
    }
}

add_comparisons <- function(df, comparison_list) {
    log_info("Adding columns with comparison values\n")

    # Input validation
    if (!is.data.frame(df)) {
        stop("Input 'df' must be a data.frame", call. = FALSE)
    }
    if (!is.list(comparison_list) || length(comparison_list) == 0) {
        stop("Input 'comparison_list' must be a non-empty list", call. = FALSE)
    }


    ## Check if all comparison expressions are valid
    #invalid_comps <- sapply(comparison_list, function(comp) !is.language(comp))
    #if (any(invalid_comps)) {
    #    stop("Invalid comparison expressions: ", 
    #         paste(names(comparison_list)[invalid_comps], collapse = ", "), 
    #         call. = FALSE)
    #}
    #                                                                                                       
    ## Get all unique column names referenced in comparisons
    #all_cols <- unique(unlist(lapply(comparison_list, all.vars)))
    #missing_cols <- setdiff(all_cols, names(df))
    #if (length(missing_cols) > 0) {
    #    stop("Columns referenced in comparisons but not in dataframe: ", 
    #         paste(missing_cols, collapse = ", "), 
    #         call. = FALSE)
    #}
    #                                                                                                       
    ## Performance optimization: pre-allocate result list
    #results <- vector("list", length(comparison_list))
    #names(results) <- names(comparison_list)
    #                                                                                                       
    ## Evaluate comparisons
    #for (comp_name in names(comparison_list)) {
    #    tryCatch({
    #        log_info(paste("Evaluating comparison:", comp_name))
    #        results[[comp_name]] <- eval(comparison_list[[comp_name]], df)
    #        if (!is.logical(results[[comp_name]])) {
    #            stop("Comparison result must be logical", call. = FALSE)
    #        }
    #        if (length(results[[comp_name]]) != nrow(df)) {
    #            stop("Comparison result length does not match number of rows in dataframe", call. = FALSE)
    #        }
    #    }, error = function(e) {
    #        stop(paste("Error in comparison", comp_name, ":", e$message), call. = FALSE)
    #    })
    #}
    #                                                                                                       
    ## Add results to dataframe
    #df[names(results)] <- results
    # For all comparisons, create column with that name and TRUE/FALSE values for rows.
    for(comp_name in names(comparison_list)) {
        df[[comp_name]] <- eval(comparison_list[[comp_name]], df)

    }
    return(df)
}

# Function to add new comparisons easily
add_new_comparison <- function(df, name, condition) {
    df[[name]] <- eval(condition, df)
    return(df)
}

add_attributes <- function(df, control_factors) {
    cat("Adding attributes to column names\n")
    control_columns <- unlist(unname(control_factors))
    at_least_one_not_in_df_column <- !all(control_columns %in% colnames(df))
    if(at_least_one_not_in_df_column){
        control_column_not_in_df <- which(!(control_columns %in% colnames(df)))
        log_error("One of the control factors is in the dataframe.")
        stop("Verify the columns in categories and control_factors list to ensure you are assigning correctly\n")
    }

    for (factor in names(control_factors)) {
        new_column_name <- paste0("X__cf_", factor)
        df[[new_column_name]] <- paste(control_factors[[factor]], collapse = ",")
    }
    return(df)

}

create_bmc_table <- function(named_samples_table) {
    cat("Making bmc_table from sample table\n")
    bmc_table <- data.frame(SampleName = named_samples_table$full_name,
       Vol..uL = 10,
       Conc = 0,
       Type = "ChIP",
       Genome = "Saccharomyces cerevisiae",
       Notes = ifelse(named_samples_table$antibody == "Input", "Run on fragment analyzer.", "Run on femto pulse."),
       Pool = "A"
    )
    return(bmc_table)
}

table_has_ID_column <- function(sample_table){
    if(!("sample_ID" %in% colnames(sample_table))){
        cat("No sample_ID column found.\n")
        cat("Must determine sample_IDs from fastq files\n")
        return(FALSE)
    } else {
        cat("Table has sample_ID column.\n")
        return(TRUE)
    }
}
modify_and_output_table <- function(sample_table, sample_ID_array, output_file_path) {
    if(nrow(sample_table) != length(sample_ID_array)) {
        cat("Number of rows is different from length of sample_ID_array.\n")
        cat("Verify fastq file names to ensure proper number is being extracted.\n")
        cat(sprintf("Number of rows: %s\n", nrow(sample_table)))
        cat(sprintf("Length of array: %s\n", length(sample_ID_array)))
        stop()
    } else if ("sample_ID" %in% colnames(sample_table)) {
        cat("sample_ID already part of the sample table.\n")
        print(colnames(sample_table))
        stop()
    } else {
        sample_table$sample_ID <- sample_ID_array
        print(head(sample_table))
        write.table(sample_table, output_file_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
        cat("Output modified sample table with sample ID column.\n")
    }
}

#' BAM File Management Functions
find_sample_bam <- function(sample_info,
                           bam_dir,
                           reference_pattern = CONFIG$FILES$PATTERNS$REFERENCE) {
    log_info("Finding BAM file for:", sample_info$short_name)
    
    bam_files <- find_matching_bams(
        sample_info$sample_ID,
        bam_dir
    )
    
    if (length(bam_files) == 0) {
        log_warning("No BAM files found for:", sample_info$sample_ID)
        return(NULL)
    }
    
    # Filter for reference genome
    ref_files <- filter_reference_bams(
        bam_files,
        reference_pattern
    )
    
    if (length(ref_files) == 0) {
        log_warning("No reference BAM files found for:", sample_info$sample_ID)
        return(NULL)
    }
    
    ref_files[1]
}

find_matching_bams <- function(sample_id, bam_dir) {
    pattern <- sprintf(".*%s.*%s", sample_id, CONFIG$FILES$PATTERNS$BAM)
    list.files(bam_dir, pattern = pattern, full.names = TRUE)
}

filter_reference_bams <- function(bam_files,
                                reference_pattern) {
    bam_files[grepl(reference_pattern, bam_files)]
}

validate_bam_file <- function(file_path, sample_info) {
    if (!file.exists(file_path)) {
        log_warning(sprintf(
            "BAM file not found for sample %s (%s)",
            sample_info$short_name,
            sample_info$sample_ID
        ))
        return(FALSE)
    }
    return(TRUE)
}

#' BAM QC Analysis Functions
analyze_bam_qc <- function(directory,
                          config = CONFIG$BAM_QC) {
    log_info("Analyzing BAM QC for:", directory)
    
    # Get flagstat files
    flagstat_files <- find_flagstat_files(directory)
    
    if (length(flagstat_files) == 0) {
        log_warning("No flagstat files found")
        return(NULL)
    }
    
    # Get genome names
    genome_names <- extract_genome_names(flagstat_files)
    
    # Get sample information
    sample_info <- load_sample_info(directory)
    
    # Calculate mapping percentages
    mapping_stats <- calculate_mapping_stats(
        flagstat_files,
        sample_info,
        genome_names
    )
    
    mapping_stats
}

find_flagstat_files <- function(directory) {
    pattern <- CONFIG$BAM_QC$PATTERNS$FLAGSTAT
    qc_dir <- file.path(directory, CONFIG$PATHS$QC_DIR)
    
    list.files(
        qc_dir,
        pattern = pattern,
        recursive = FALSE,
        full.names = TRUE
    )
}

extract_genome_names <- function(files) {
    unique(sapply(files, function(path) {
        parts <- strsplit(basename(path), "_")[[1]]
        parts[length(parts) - 1]
    }))
}

load_sample_info <- function(directory) {
    pattern <- CONFIG$BAM_QC$PATTERNS$SAMPLE_INFO
    doc_dir <- file.path(directory, CONFIG$PATHS$DOC_DIR)
    
    info_file <- list.files(
        doc_dir,
        pattern = pattern,
        recursive = FALSE,
        full.names = TRUE
    )[1]
    
    if (is.na(info_file)) {
        log_error("Sample info file not found")
        return(NULL)
    }
    
    read.table(info_file, sep = ",", header = TRUE)
}

#' BigWig File Processing Functions
find_bigwig_file <- function(directory,
                           sample_id,
                           pattern) {
    log_info("Finding bigwig file for:", sample_id)
    
    matches <- list.files(
        directory,
        pattern = sample_id,
        full.names = TRUE,
        recursive = TRUE
    )
    
    bigwig_files <- matches[grepl(pattern, matches)]
    
    if (length(bigwig_files) == 0) {
        log_warning("No bigwig file found for:", sample_id)
        return(NULL)
    }
    
    if (length(bigwig_files) > 1) {
        log_warning("Multiple bigwig files found, using first")
    }
    
    bigwig_files[1]
}

import_bigwig_data <- function(file_path, region) {
    log_info("Importing bigwig data from:", file_path)
    
    tryCatch({
        import(con = file_path, which = region)
    }, error = function(e) {
        log_error("Failed to import bigwig:", e$message)
        NULL
    })
}

#' Chromosome Name Conversion Functions
convert_chromosome_names <- function(names,
                                   target_style,
                                   config = CONFIG$CHROMOSOMES) {
    log_info("Converting chromosome names to style:", target_style)
    
    vapply(
        names,
        function(chr) convert_single_chromosome(chr, target_style, config),
        character(1)
    )
}

convert_single_chromosome <- function(chr,
                                    target_style,
                                    config = CONFIG$CHROMOSOMES) {
    # Remove existing prefix if present
    chr <- gsub(CONFIG$NAMING$PATTERNS$PREFIX, "", chr, 
                ignore.case = TRUE)
    
    # Handle special chromosomes
    if (chr %in% config$SPECIAL) {
        return(paste0(config$PREFIX, chr))
    }
    
    # Convert based on current format
    if (grepl(CONFIG$NAMING$PATTERNS$ROMAN, chr)) {
        return(handle_roman_chromosome(chr, target_style))
    } else if (grepl(CONFIG$NAMING$PATTERNS$NUMERIC, chr)) {
        return(handle_numeric_chromosome(chr, target_style))
    }
    
    log_warning("Unable to convert chromosome:", chr)
    return(paste0(config$PREFIX, chr))
}

handle_roman_chromosome <- function(chr, target_style) {
    switch(target_style,
        "ENSEMBL" = as.character(arabic(chr)),
        "UCSC" = paste0(CONFIG$CHROMOSOMES$PREFIX, chr),
        "NCBI" = chr,
        stop("Unknown target style:", target_style)
    )
}

handle_numeric_chromosome <- function(chr, target_style) {
    switch(target_style,
        "ENSEMBL" = chr,
        "UCSC" = paste0(CONFIG$CHROMOSOMES$PREFIX, as.roman(as.integer(chr))),
        "NCBI" = as.roman(as.integer(chr)),
        stop("Unknown target style:", target_style)
    )
}

#' Control Sample Management Functions
find_control_sample <- function(sample_row,
                              sample_table,
                              factors,
                              bigwig_dir,
                              pattern) {
    log_info("Finding control sample")
    
    control_index <- determine_matching_control(
        sample_row,
        sample_table,
        factors
    )
    
    if (length(control_index) == 0) {
        log_warning("No matching control found, using fallback")
        return(find_fallback_control(sample_table, bigwig_dir, pattern))
    }
    
    control_data <- get_control_data(
        sample_table[control_index, ],
        bigwig_dir,
        pattern
    )
    
    if (is.null(control_data)) {
        log_warning("Control data not found, using fallback")
        return(find_fallback_control(sample_table, bigwig_dir, pattern))
    }
    
    control_data
}

get_control_data <- function(control_sample,
                           bigwig_dir,
                           pattern) {
    log_info("Getting control data for:", control_sample$sample_ID)
    
    bigwig_file <- find_bigwig_file(
        bigwig_dir,
        control_sample$sample_ID,
        pattern
    )
    
    if (is.null(bigwig_file)) {
        return(NULL)
    }
    
    list(
        id = control_sample$sample_ID,
        name = control_sample$short_name,
        file = bigwig_file
    )
}

#' Enhanced control sample management
find_valid_control <- function(sample_row,
                             sample_table,
                             bam_dir,
                             factors) {
    log_info("Finding valid control sample")
    
    # Find matching control
    control_index <- determine_matching_control(
        sample_row,
        sample_table,
        factors
    )
    
    # Validate and possibly adjust control
    control_index <- validate_control_index(
        control_index,
        sample_table,
        bam_dir
    )
    
    if (control_index > 0) {
        return(get_control_bam(
            sample_table[control_index, ],
            bam_dir
        ))
    }
    
    log_warning("No valid control found")
    return(NULL)
}

validate_control_index <- function(index,
                                 sample_table,
                                 bam_dir) {
    if (length(index) == 0) {
        log_warning("No matching control found, using fallback")
        return(find_fallback_control(sample_table, bam_dir))
    }
    
    if (length(index) > CONFIG$CONTROL$MAX_CONTROLS) {
        log_warning("Multiple controls found, using first")
        index <- index[1]
    }
    
    # Validate BAM file exists
    control_bam <- find_matching_bam(
        sample_table$sample_ID[index],
        bam_dir
    )
    
    if (!is.null(control_bam)) {
        return(index)
    }
    
    log_warning("Control BAM not found, searching for alternative")
    return(find_fallback_control(sample_table, bam_dir))
}

find_fallback_control <- function(sample_table, bam_dir) {
    log_info("Searching for fallback control")
    
    input_samples <- sample_table[
        sample_table[[CONFIG$CONTROL$ANTIBODY_COLUMN]] == 
        CONFIG$CONTROL$INPUT_VALUE,
    ]
    
    for (i in seq_len(nrow(input_samples))) {
        if (!is.null(find_matching_bam(
            input_samples$sample_ID[i],
            bam_dir
        ))) {
            return(which(sample_table$sample_ID == 
                        input_samples$sample_ID[i]))
        }
    }
    
    return(CONFIG$CONTROL$DEFAULT_INDEX)
}

#' Data Conversion Functions
convert_to_granges <- function(data,
                             file_name,
                             config = CONFIG$FEATURES) {
    log_info("Converting data to GRanges")
    
    if (is(data, "GRanges")) {
        return(data)
    }
    
    converter <- determine_converter(file_name)
    converter(data, file_name)
}

determine_converter <- function(file_name) {
    for (type in names(CONFIG$FEATURES$FILE_TYPES)) {
        if (grepl(CONFIG$FEATURES$FILE_TYPES[[type]], file_name)) {
            return(get(paste0("convert_", tolower(type))))
        }
    }
    
    function(data, file_name) {
        log_warning("No specific converter for:", file_name)
        as(data, "GRanges")
    }
}

convert_nucleosome <- function(data, file_name) {
    log_info("Converting nucleosome data")
    
    config <- CONFIG$FEATURES$COLUMNS$NUCLEOSOME
    
    # Create metadata
    metadata <- data[, !(colnames(data) %in% config$EXCLUDE)]
    colnames(metadata) <- clean_column_names(colnames(metadata))
    
    GRanges(
        seqnames = data$`Nucleosome ID`,
        ranges = IRanges(
            start = data[[config$POSITION]],
            end = data[[config$POSITION]]
        ),
        strand = "*",
        chromosome = data$Chromosome,
        metadata
    )
}

convert_timing <- function(data, file_name) {
    log_info("Converting timing data")
    
    config <- CONFIG$FEATURES$COLUMNS$TIMING
    window <- config$WINDOW
    
    # Create origin names
    names <- paste0(data$Chromosome, "_", data$Position)
    
    # Create metadata
    metadata <- data[, !(colnames(data) %in% config$EXCLUDE)]
    colnames(metadata) <- clean_column_names(colnames(metadata))
    
    GRanges(
        seqnames = names,
        ranges = IRanges(
            start = data$Position - window,
            end = data$Position + window
        ),
        strand = "*",
        chromosome = data$Chromosome,
        metadata
    )
}

#' Data Download Functions
download_feature_data <- function(source_config,
                                base_dir,
                                timestamp) {
    log_info("Downloading feature data")
    
    results <- list()
    
    for (source in names(source_config)) {
        result <- download_source_data(
            source,
            source_config[[source]],
            base_dir,
            timestamp
        )
        
        results[[source]] <- result
    }
    
    results
}

download_source_data <- function(source_name,
                               config,
                               base_dir,
                               timestamp) {
    log_info("Downloading:", source_name)
    
    output_file <- generate_output_path(
        base_dir,
        config$FILE,
        timestamp
    )
    
    # Download file
    download_file(config$URL, output_file)
    
    # Process if compressed
    if (is_compressed(output_file)) {
        output_file <- decompress_file(output_file)
    }
    
    # Validate and process
    process_downloaded_file(
        output_file,
        source_name,
        config$TYPE
    )
}

download_file <- function(url, dest_path) {
    log_info("Downloading:", basename(dest_path))
    
    tryCatch({
        curl::curl_download(url, dest_path, mode = "wb")
        log_info("Download complete")
        TRUE
    }, error = function(e) {
        log_error("Download failed:", e$message)
        FALSE
    })
}

#' Data Preparation Functions
prepare_visualization_data <- function(directory,
                                    chromosome = CONFIG$DEFAULTS$CHROMOSOME,
                                    genome_dir = CONFIG$DEFAULTS$GENOME_DIR,
                                    genome_pattern = CONFIG$DEFAULTS$GENOME_PATTERN) {
    log_info("Preparing visualization data")
    
    # Validate directory
    directory_path <- validate_input(directory)
    
    # Load sample table
    log_info("Loading sample table")
    sample_table <- load_sample_table(directory_path)
    
    # Load reference genome
    log_info("Loading reference genome")
    ref_genome <- load_reference_genome(
        genome_dir = genome_dir,
        genome_pattern = genome_pattern
    )
    
    # Create genome range
    log_info("Creating genome range")
    genome_range <- create_chromosome_GRange(ref_genome)
    
    list(
        directory = directory_path,
        samples = sample_table,
        genome = ref_genome,
        range = genome_range
    )
}

prepare_feature_track <- function(chromosome,
                                genome_range,
                                pattern = CONFIG$DEFAULTS$FEATURE_PATTERN) {
    log_info("Preparing feature track")
    
    feature_data <- load_feature_file_GRange(
        chromosome_to_plot = chromosome,
        feature_file_pattern = pattern,
        genomeRange_to_get = genome_range
    )
    
    AnnotationTrack(
        feature_data,
        name = paste("Origin Peaks", "Eaton 2010", sep = "")
    )
}

#!/usr/bin/env Rscript
# R/functions/environment_utils.R

#' Environment Management Utilities
#' @description Centralized environment management system
#' @export

#' Validate System Environment
#' @param config List Configuration settings
#' @param min_memory Numeric Minimum required memory in GB
#' @return List Environment status
validate_system_environment <- function(
    config = PROJECT_CONFIG,
    min_memory = 4
) {
    tryCatch({
        # Check R version
        check_r_version(
            R_system_version = getRversion(),
            R_project_version = package_version(config$SYSTEM$R_MIN_VERSION)
        )
        
        # Check memory
        check_system_memory(
            min_memory = min_memory,
            current_memory = memory.limit()
        )
        
        # Check environment variables
        check_environment_variables(
            required_vars = config$ENVIRONMENT$REQUIRED_VARS,
            current_vars = names(Sys.getenv())
        )
        
        # Log system info
        log_system_info()
        
        invisible(TRUE)
    }, error = function(e) {
        log_error("Environment validation failed:", e$message)
        stop(e)
    })
}

#' Check System Memory
#' @param min_memory Numeric Required memory in GB
#' @param current_memory Numeric Available memory in MB
#' @return Logical TRUE if check passes
check_system_memory <- function(
    min_memory,
    current_memory
) {
    if (current_memory < min_memory * 1024) {
        log_warning(sprintf(
            "Low memory: %d MB available, %d GB required",
            current_memory,
            min_memory
        ))
    }
    invisible(TRUE)
}

#' Check Environment Variables
#' @param required_vars Character vector Required variables
#' @param current_vars Character vector Current variables
#' @return Logical TRUE if check passes
check_environment_variables <- function(
    required_vars,
    current_vars
) {
    missing_vars <- setdiff(required_vars, current_vars)
    if (length(missing_vars) > 0) {
        stop(sprintf(
            "Missing environment variables: %s",
            paste(missing_vars, collapse = ", ")
        ))
    }
    invisible(TRUE)
}

#' Setup Project Paths
#' @param config List Configuration settings
#' @return List Path information
setup_project_paths <- function(
    config = PROJECT_CONFIG
) {
    tryCatch({
        # Validate base paths
        for (path in c(config$PATHS$ROOT, "~/logs")) {
            if (!dir.exists(path)) {
                dir.create(path, recursive = TRUE)
            }
        }
        
        # Return path information
        list(
            root = normalizePath(config$PATHS$ROOT),
            logs = normalizePath("~/logs"),
            data = if (dir.exists("~/data")) normalizePath("~/data") else NULL
        )
    }, error = function(e) {
        log_error("Path setup failed:", e$message)
        stop(e)
    })
}

#' Clean R Environment
#' @param keep Character vector Names to keep
#' @return None
clean_environment <- function(
    keep = c("PROJECT_CONFIG")
) {
    # Detach packages
    attached <- paste0("package:", names(sessionInfo()$otherPkgs))
    for (pkg in attached) {
        try(detach(pkg, character.only = TRUE, unload = TRUE), silent = TRUE)
    }
    
    # Clear workspace except kept objects
    rm(list = setdiff(ls(all.names = TRUE), keep), envir = .GlobalEnv)
    
    # Force garbage collection
    gc()
    
    invisible(TRUE)
}

#' Feature File Processing Functions
load_feature_file_GRange <- function(chromosome_to_plot = 10,
                                   feature_file_pattern = CONFIG$FEATURE_TYPES$PEAKS,
                                   genomeRange_to_get) {
    log_info("Loading feature file:", feature_file_pattern)
    
    # Validate directory
    feature_dir <- CONFIG$PATHS$FEATURE_DIR
    if (!dir.exists(feature_dir)) {
        log_error("Feature directory not found:", feature_dir)
        stop("Missing feature directory")
    }
    
    # Find feature file
    feature_files <- list.files(feature_dir,
                              pattern = feature_file_pattern,
                              full.names = TRUE,
                              recursive = TRUE)
    
    if (length(feature_files) != 1) {
        log_error("Invalid number of feature files:", length(feature_files))
        stop("Feature file error")
    }
    
    # Load and process features
    feature_grange <- tryCatch({
        import.bed(feature_files[1])
    }, error = function(e) {
        log_error("Failed to import feature file:", e$message)
        stop("Import error")
    })
    
    # Process chromosome styles
    feature_style <- determine_chr_style(seqlevels(feature_grange))
    genome_style <- determine_chr_style(seqlevels(genomeRange_to_get))
    
    log_info("Feature style:", feature_style)
    log_info("Genome style:", genome_style)
    
    # Normalize chromosome styles
    feature_grange_subset <- normalize_feature_range(
        feature_grange,
        genomeRange_to_get,
        feature_style,
        genome_style
    )
    
    return(feature_grange_subset)
}

normalize_feature_range <- function(feature_grange,
                                  genome_range,
                                  feature_style,
                                  genome_style) {
    if (feature_style == genome_style) {
        log_info("Styles match, using direct subsetting")
        return(subsetByOverlaps(feature_grange, genome_range))
    }
    
    log_info("Adjusting styles for compatibility")
    
    # Adjust genome range to match feature style
    adjusted_genome_range <- genome_range
    new_seqlevels <- normalize_chr_names(seqlevels(genome_range),
                                       feature_style)
    seqlevels(adjusted_genome_range) <- new_seqlevels
    
    # Subset features
    feature_subset <- subsetByOverlaps(feature_grange,
                                     adjusted_genome_range)
    
    # Convert back to genome style
    final_seqlevels <- normalize_chr_names(seqlevels(feature_subset),
                                         genome_style)
    seqlevels(feature_subset) <- final_seqlevels
    
    return(feature_subset)
}
#' Feature File Processing Functions
process_feature_files <- function(input_dir = CONFIG$FEATURES$PATHS$BASE_DIR) {
    log_info("Processing feature files")
    
    # Validate directory
    validate_directory(input_dir)
    
    # Get file list
    files <- get_feature_files(input_dir)
    
    if (length(files) == 0) {
        log_warning("No files to process")
        return(NULL)
    }
    
    # Process each file
    results <- process_files(files)
    
    # Verify outputs
    verify_outputs(input_dir)
    
    results
}

process_files <- function(files) {
    log_info("Processing", length(files), "files")
    
    lapply(files, function(file) {
        tryCatch({
            process_single_file(file)
        }, error = function(e) {
            log_error("Failed to process:", basename(file))
            log_error("Error:", e$message)
            NULL
        })
    })
}

process_single_file <- function(file_path) {
    log_info("Processing file:", basename(file_path))
    
    # Read data
    data <- read_feature_file(file_path)
    
    # Process data
    processed <- process_feature_data(data, basename(file_path))
    
    # Convert to GRanges
    granges <- convert_to_granges(processed, basename(file_path))
    
    # Output results
    output_results(granges, file_path)
    
    granges
}
#' Feature Data Processing Functions
process_downloaded_file <- function(file_path,
                                  source_name,
                                  type) {
    log_info("Processing file:", basename(file_path))
    
    # Read file
    data <- read_feature_file(file_path)
    
    # Process based on type
    processed <- switch(type,
        "features" = process_sgd_features(data),
        "acs" = process_eaton_acs(data, file_path),
        process_generic_features(data)
    )
    
    validate_processed_data(processed, type)
}

process_sgd_features <- function(data) {
    log_info("Processing SGD features")
    
    if (ncol(data) != length(CONFIG$SGD_COLUMNS)) {
        log_error("Invalid column count in SGD features")
        return(NULL)
    }
    
    names(data) <- CONFIG$SGD_COLUMNS
    data
}

process_eaton_acs <- function(data, file_path) {
    log_info("Processing Eaton ACS data")
    
    # Fix coordinates
    fixed <- fix_coordinates(data)
    
    # Save fixed version if needed
    if (has_invalid_coordinates(data)) {
        save_fixed_coordinates(fixed, file_path)
    }
    
    fixed
}

fix_coordinates <- function(data) {
    log_info("Fixing coordinates")
    
    data %>%
        mutate(
            width = end - start,
            temp = ifelse(start > end, start, end),
            start = ifelse(start > end, end, start),
            end = temp
        ) %>%
        select(-temp, -width)
}

load_sample_table <- function(directory_path) {
    cat("Loading sample_table from", directory_path, "\n")
    documentation_dir_path <- file.path(directory_path, "documentation")
    sample_table_path <- list.files(documentation_dir_path, pattern = "sample_table", full.names = TRUE)
    if(length(sample_table_path) == 0){
        cat(sprintf("No files with pattern sample_table found in %s\n", documentation_dir_path))
        cat("Consult 000_setupExperimentDir to create sample_table.tsv for the project\n")
        stop()
    } else if(length(sample_table_path) > 1){
        cat(sprintf("Multiple files with pattern sample_table found in %s\n", documentation_dir_path))
        cat(sprintf("Files found in %s\n", documentation_dir_path))
        print(sample_table)
        cat("Consult 000_setupExperimentDir and ensure no duplicates are present for the project\n")
        stop()
    }
    sample_table <- read.delim(sample_table_path, header = TRUE, sep ="\t")
    cat(sprintf("Reading %s\n", sample_table_path))
    cat("Head of sample_table\n")
    print(head(sample_table))
    return(sample_table)
}

determine_sample_id <- function(directory_path) {
    fastq_directory_path <- file.path(directory_path, "fastq")
    fastq_file_paths <- list.files(fastq_directory_path, pattern = "*.fastq", full.names = TRUE)
    if(length(fastq_file_paths) == 0) {
        cat(sprintf("No fastq files found in %s\n", fastq_directory_path))
        cat("Consult 000_setupExperimentDir to create sample_table.tsv for the project and download files from BMC\n")
        stop()
    } else {
        cat(sprintf("Found %s files in %s.\n", length(fastq_file_paths), fastq_directory_path))
    }
    fastq_file_names <- basename(fastq_file_paths)
    ID_regex <- "\\d{5,6}"
    fastq_split_string_list <- strsplit(fastq_file_names, "_|-")
    sample_IDs <- lapply(fastq_split_string_list, function(fastq_split_string_list) {
        for(split_string in fastq_split_string_list) {
            if(grepl(ID_regex, split_string)) {
                return(split_string)
            }
        }
    })
    if(!all(unlist(lapply(sample_IDs, length)) == 1)) {
        cat("At least one of the files did not extract exactly one sample ID.\n")
        cat("Files with problems:\n")
        print(fastq_file_names[unlist(lapply(sample_IDs, length)) != 1])
        cat("Verify sample names. Redownload from BMC if necessary.\n")
        cat(sprintf("Regex pattern used %s:\n", ID_regex))
        stop()
    } else {
        sample_IDs <- unlist(sample_IDs)
        cat(sprintf("Found %s sample_IDs.\n", length(sample_IDs)))
        cat(sprintf("First sample_ID: %s\n",sample_IDs[1]))
        cat(sprintf("Last sample_ID: %s\n", sample_IDs[length(sample_IDs)]))
        cat("Returning sample_ID array.\n")
        return(sample_IDs)
    }
}

#' Genome Loading Functions
load_reference_genome <- function(genome_dir = CONFIG$PATHS$GENOME_DIR, 
                                genome_pattern = CONFIG$PATTERNS$GENOME) {
    log_info("Loading reference genome")
    
    directory_path <- file.path(CONFIG$PATHS$BASE_DIR, genome_dir)
    
    if (!dir.exists(directory_path)) {
        log_error("Genome directory not found:", directory_path)
        stop("Missing genome directory")
    }
    
    genome_files <- list.files(directory_path, pattern = genome_pattern, 
                             full.names = TRUE, recursive = TRUE)
    
    if (length(genome_files) != 1) {
        log_error("Invalid number of genome files found:", length(genome_files))
        stop("Genome file error")
    }
    
    genome <- readFasta(genome_files[1])
    
    genome_df <- data.frame(
        chrom = names(as(genome, "DNAStringSet")),
        basePairSize = width(genome)
    ) %>% 
        filter(chrom != "chrM")
    
    return(genome_df)
}

#' Chromosome Name Processing Functions
normalize_chr_names <- function(chr_names, target_style) {
    log_info("Normalizing chromosome names")
    
    if (!target_style %in% CONFIG$CHROMOSOME_MAPPING$STYLES) {
        log_error("Invalid chromosome style:", target_style)
        stop("Invalid style")
    }
    
    chr_names <- gsub(paste0("^", CONFIG$CHROMOSOME_MAPPING$PREFIX), "", 
                     chr_names)
    
    normalized <- switch(
        target_style,
        "UCSC" = paste0(CONFIG$CHROMOSOME_MAPPING$PREFIX, chr_names),
        "Roman" = paste0(CONFIG$CHROMOSOME_MAPPING$PREFIX, 
                        mapvalues(chr_names, 
                                from = as.character(1:16),
                                to = CONFIG$CHROMOSOME_MAPPING$ROMAN)),
        "Numeric" = mapvalues(chr_names,
                            from = CONFIG$CHROMOSOME_MAPPING$ROMAN,
                            to = as.character(1:16))
    )
    
    return(unname(normalized))
}

#' GRanges Conversion Functions
convert_granges_style <- function(gr,
                                target_style,
                                verbose = FALSE) {
    log_info("Converting GRanges object to style:", target_style)
    
    validate_granges(gr)
    
    # Get current seqnames
    current_seqnames <- as.character(seqnames(gr))
    if (verbose) {
        log_info("Original seqnames:",
                paste(unique(current_seqnames), collapse = ", "))
    }
    
    # Convert names
    new_seqnames <- convert_chromosome_names(
        current_seqnames,
        target_style
    )
    
    if (verbose) {
        log_info("Converted seqnames:",
                paste(unique(new_seqnames), collapse = ", "))
        log_info("Changes made:",
                sum(new_seqnames != current_seqnames))
    }
    
    # Update GRanges object
    update_granges_seqnames(gr, new_seqnames)
}

validate_granges <- function(gr) {
    if (!is(gr, "GRanges")) {
        log_error("Invalid input type")
        stop("Input must be a GenomicRanges object")
    }
}

update_granges_seqnames <- function(gr, new_names) {
    seqlevels(gr) <- unique(new_names)
    seqnames(gr) <- new_names
    gr
}

#' Mapping Statistics Functions
calculate_mapping_stats <- function(flagstat_files,
                                  sample_info,
                                  genome_names) {
    log_info("Calculating mapping statistics")
    
    # Initialize results matrix
    results <- initialize_results_matrix(genome_names)
    
    # Process each sample
    sample_ids <- get_sample_ids(sample_info)
    
    for (sample_id in sample_ids) {
        stats <- process_sample_stats(
            sample_id,
            flagstat_files,
            genome_names
        )
        
        results <- rbind(results, stats)
    }
    
    results
}

process_sample_stats <- function(sample_id,
                               flagstat_files,
                               genome_names) {
    log_info("Processing sample:", sample_id)
    
    # Find relevant flagstat files
    sample_files <- find_sample_flagstats(
        sample_id,
        flagstat_files
    )
    
    # Calculate percentages
    stats <- sapply(sample_files, function(file) {
        calculate_mapping_percentage(file)
    })
    
    # Create result row
    result <- as.data.frame(t(stats))
    colnames(result) <- genome_names
    
    result
}

calculate_mapping_percentage <- function(flagstat_file) {
    log_info("Calculating mapping percentage for:", 
             basename(flagstat_file))
    
    # Read flagstat data
    data <- read.table(flagstat_file, sep = "\t")
    
    # Find relevant rows
    patterns <- CONFIG$BAM_QC$METRICS
    is_relevant <- Reduce("|", lapply(patterns, grepl, data[,3]))
    subset_data <- data[is_relevant, ]
    
    # Calculate percentage
    mapped <- as.numeric(subset_data[2,1])
    total <- as.numeric(subset_data[1,1])
    
    (mapped / total) * 100
}

output_tables_in_list <- function(experiment_directory, list_of_tables, OUTPUT_TABLE = FALSE){
        experiment_name <- basename(experiment_directory)
        if (!(typeof(list_of_tables) == "list")){
            stop("Argument must be a list.")
        }
        names_of_tables <- names(list_of_tables)
        for (name_of_table in names_of_tables){
            output_table <- sample_config_output[[name_of_table]]
            cat("============\n")
            print(head(output_table))
            output_file_path <- file.path(experiment_directory, "documentation", paste(experiment_name, "_", name_of_table, ".tsv", sep = ""))
            cat(sprintf("Outputting to %s: \n", output_file_path))
            if(OUTPUT_TABLE) {
                write.table(output_table, file = output_file_path, sep = "\t", row.names = FALSE)
            } else {
                cat("Skip writing table. MODIFY OUTPUT_TABLE value to output.\n")
            }
        }
}

#' Package Installation Functions
install_package_group <- function(packages,
                                manager = "BiocManager") {
    log_info("Installing package group using:", manager)
    
    results <- switch(manager,
        "BiocManager" = install_bioc_packages(packages),
        "renv" = install_renv_packages(packages),
        "github" = install_github_packages(packages),
        stop("Unknown package manager:", manager)
    )
    
    failed <- names(results)[!results]
    if (length(failed) > 0) {
        log_warning("Failed to install packages:",
                   paste(failed, collapse = ", "))
    }
    
    invisible(results)
}

install_bioc_packages <- function(packages) {
    sapply(packages, function(pkg) {
        tryCatch({
            BiocManager::install(pkg, update = FALSE)
            TRUE
        }, error = function(e) {
            log_error("Failed to install:", pkg)
            FALSE
        })
    })
}

install_renv_packages <- function(packages) {
    sapply(packages, function(pkg) {
        tryCatch({
            renv::install(pkg)
            TRUE
        }, error = function(e) {
            log_error("Failed to install:", pkg)
            FALSE
        })
    })
}

install_github_packages <- function(packages) {
    sapply(packages, function(pkg) {
        tryCatch({
            remotes::install_github(pkg)
            TRUE
        }, error = function(e) {
            log_error("Failed to install:", pkg)
            FALSE
        })
    })
}

#' Package Management Functions
load_required_packages <- function(packages = CONFIG$PACKAGES$REQUIRED) {
    log_info("Loading required packages")
    
    suppressPackageStartupMessages({
        results <- sapply(packages, require, character.only = TRUE)
    })
    
    missing <- packages[!results]
    if (length(missing) > 0) {
        log_error("Missing packages:", paste(missing, collapse = ", "))
        stop("Required packages not available")
    }
    
    return(TRUE)
}
#' Package Management Functions
initialize_environment <- function(config = CONFIG) {
    log_info("Initializing R environment")
    
    # Initialize renv
    setup_renv()
    
    # Setup package management
    setup_package_manager()
    
    # Install and load packages
    install_required_packages(config$PACKAGES)
    
    # Snapshot environment
    snapshot_environment()
    
    log_info("Environment initialization complete")
}

setup_renv <- function() {
    log_info("Setting up renv")
    
    if (!requireNamespace("renv", quietly = TRUE)) {
        log_info("Installing renv")
        install.packages("renv")
    }
    
    renv::init()
}

setup_package_manager <- function() {
    log_info("Setting up BiocManager")
    
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        log_info("Installing BiocManager")
        install.packages("BiocManager")
    }
}

install_required_packages <- function(package_lists) {
    log_info("Installing required packages")
    
    # Combine all packages
    all_packages <- unique(unlist(package_lists))
    
    # Install packages
    result <- tryCatch({
        BiocManager::install(all_packages)
        TRUE
    }, error = function(e) {
        log_error("Package installation failed:", e$message)
        FALSE
    })
    
    if (!result) {
        stop("Package installation failed")
    }
}

load_packages <- function(packages) {
    log_info("Loading packages")
    
    results <- sapply(packages, function(pkg) {
        tryCatch({
            library(pkg, character.only = TRUE)
            TRUE
        }, error = function(e) {
            log_warning("Failed to load package:", pkg)
            FALSE
        })
    })
    
    if (!all(results)) {
        failed <- names(results)[!results]
        log_warning("Failed to load packages:", 
                   paste(failed, collapse = ", "))
    }
}

snapshot_environment <- function() {
    log_info("Creating environment snapshot")
    
    tryCatch({
        renv::snapshot()
    }, error = function(e) {
        log_error("Failed to create snapshot:", e$message)
    })
}

#' Plot Generation Functions
generate_track_plot <- function(tracks,
                              title,
                              chromosome,
                              output_file = NULL) {
    log_info("Generating track plot:", title)
    
    if (!is.null(output_file)) {
        svg(output_file)
        on.exit(dev.off())
    }
    
    plotTracks(
        tracks,
        main = title,
        chromosome = chromosome
    )
}

create_output_filename <- function(base_dir,
                                 chromosome,
                                 pattern,
                                 comparison,
                                 time_id = NULL) {
    if (is.null(time_id)) {
        time_id <- format(Sys.time(), CONFIG$PLOT$DATE_FORMAT)
    }
    
    pattern_clean <- gsub("_|\\.bw", "", pattern)
    comparison_clean <- gsub("_", "", comparison)
    
    file.path(
        base_dir,
        sprintf(
            "%s_%s_%s_%s_%s.svg",
            time_id,
            chromosome,
            pattern_clean,
            comparison_clean,
            format(Sys.time(), "%H%M%S")
        )
    )
}
#' Plot Generation Functions
generate_sample_plot <- function(tracks,
                               chromosome,
                               title,
                               output_file = NULL,
                               config = CONFIG$VISUALIZATION) {
    log_info("Generating plot:", title)
    
    if (!is.null(output_file)) {
        svg(output_file)
        on.exit(dev.off())
    }
    
    plotTracks(
        tracks,
        main = title,
        chromosome = chromosome,
        ylim = c(config$LIMITS$Y_MIN, config$LIMITS$Y_MAX)
    )
}

generate_output_filename <- function(base_dir,
                                   sample_name,
                                   chromosome,
                                   config = CONFIG$VISUALIZATION) {
    log_info("Generating output filename")
    
    date_str <- format(Sys.time(), config$OUTPUT$DATE_FORMAT)
    
    file.path(
        base_dir,
        sprintf(
            "%s_%s_%s_WithInputAndHighlights.%s",
            date_str,
            chromosome,
            sample_name,
            config$OUTPUT$FORMAT
        )
    )
}

#' Main Plot Management Function
plot_all_sample_tracks <- function(sample_table,
                                 directory_path,
                                 chromosome_to_plot = 10,
                                 genomeRange_to_get,
                                 control_track,
                                 annotation_track,
                                 highlight_gr) {
    log_info("Starting sample track plotting")
    
    # Setup
    setup_result <- setup_plotting_environment(
        directory_path,
        chromosome_to_plot
    )
    
    # Process samples
    for (sample_index in seq_len(nrow(sample_table))) {
        if (sample_index == 1) {  # Only process first sample for now
            process_sample(
                sample_table[sample_index, ],
                sample_table,
                setup_result,
                control_track,
                annotation_track,
                highlight_gr
            )
        }
    }
    
    log_info("Plotting completed")
}

#' Helper Functions
setup_plotting_environment <- function(directory_path,
                                     chromosome) {
    log_info("Setting up plotting environment")
    
    list(
        plot_dir = file.path(directory_path, CONFIG$PATHS$SUBDIRS$PLOTS),
        bigwig_dir = file.path(directory_path, CONFIG$PATHS$SUBDIRS$BIGWIG),
        chromosome = paste0("chr", as.roman(chromosome)),
        title = sprintf("Complete View of Chrom %s", chromosome)
    )
}

process_sample <- function(sample_row,
                         sample_table,
                         setup,
                         control_track,
                         annotation_track,
                         highlight_gr) {
    log_info("Processing sample:", sample_row$short_name)
    
    # Get sample data
    sample_data <- get_sample_data(
        sample_row,
        setup$bigwig_dir,
        genomeRange_to_get
    )
    
    if (is.null(sample_data)) {
        log_warning("No data for sample:", sample_row$short_name)
        return(NULL)
    }
    
    # Get control data
    control_data <- get_control_data(
        sample_row,
        sample_table,
        setup$bigwig_dir,
        genomeRange_to_get
    )
    
    # Create tracks
    tracks <- create_track_set(
        sample_data,
        control_data,
        annotation_track,
        highlight_gr,
        setup$chromosome
    )
    
    # Generate plot
    output_file <- generate_output_filename(
        setup$plot_dir,
        sample_row$short_name,
        setup$chromosome
    )
    
    generate_sample_plot(
        tracks,
        setup$chromosome,
        setup$title,
        output_file
    )
}

#' Genome Range Management Functions
create_chromosome_range <- function(genome_data,
                                  config = CONFIG$GENOME) {
    log_info("Creating chromosome range")
    
    validate_genome_data(genome_data)
    
    GRanges(
        seqnames = genome_data$chrom,
        ranges = IRanges(
            start = 1,
            end = genome_data$basePairSize
        ),
        strand = config$DEFAULT_STRAND
    )
}

load_control_range <- function(control_dir,
                             identifier,
                             chromosome,
                             genome_range) {
    log_info("Loading control range data")
    
    bigwig_file <- find_control_bigwig(
        control_dir,
        identifier
    )
    
    if (is.null(bigwig_file)) {
        log_error("Control bigwig file not found")
        return(NULL)
    }
    
    import_control_data(
        bigwig_file,
        chromosome,
        genome_range
    )
}

import_control_data <- function(file_path,
                              chromosome,
                              genome_range) {
    log_info("Importing control data")
    
    control_style <- determine_chr_style(
        seqlevels(import(file_path))
    )
    
    chromosome_name <- normalize_chr_names(
        chromosome,
        control_style
    )
    
    subset_range <- genome_range[
        seqnames(genome_range) == chromosome_name
    ]
    
    import(file_path, which = subset_range)
}

#' Sample Labeling Functions
unique_labeling <- function(table,
                          categories_for_label,
                          config = CONFIG$LABEL_CONFIG) {
    log_info("Creating unique labels for samples")
    
    # Validate inputs
    validate_labeling_inputs(table, categories_for_label)
    
    # Ensure required categories
    categories_for_label <- ensure_required_categories(categories_for_label)
    
    # Get unique values
    unique_values <- get_unique_category_values(table, categories_for_label)
    
    # Create labels
    labels <- create_sample_labels(table,
                                 categories_for_label,
                                 unique_values,
                                 config)
    
    return(labels)
}

validate_labeling_inputs <- function(table, categories) {
    if (!is.data.frame(table)) {
        log_error("Invalid input table type")
        stop("Table must be a data frame")
    }
    
    if (!is.character(categories) || length(categories) == 0) {
        log_error("Invalid categories specification")
        stop("Categories must be non-empty character vector")
    }
    
    missing_cats <- setdiff(categories, colnames(table))
    if (length(missing_cats) > 0) {
        log_error("Missing categories:", paste(missing_cats, collapse = ", "))
        stop("Missing categories in table")
    }
}

ensure_required_categories <- function(categories) {
    required <- CONFIG$REQUIRED_CATEGORIES$BASE
    if (!required %in% categories) {
        log_info("Adding required category:", required)
        categories <- c(required, categories)
    }
    return(categories)
}

create_sample_labels <- function(table,
                               categories,
                               unique_values,
                               config) {
    log_info("Constructing sample labels")
    
    labels <- apply(table, 1, function(sample) {
        # Get relevant category values
        values <- sapply(categories, function(cat) {
            if (length(unique_values[[cat]]) > 1 || 
                cat == CONFIG$REQUIRED_CATEGORIES$BASE) {
                return(sample[cat])
            }
            return(NULL)
        })
        
        # Filter and combine values
        values <- values[!sapply(values, is.null)]
        label <- paste(values, collapse = config$SEPARATOR)
        
        # Truncate if necessary
        if (nchar(label) > config$MAX_LENGTH) {
            label <- paste0(substr(label, 1,
                                 config$MAX_LENGTH - nchar(config$TRUNCATE_SUFFIX)),
                          config$TRUNCATE_SUFFIX)
        }
        
        return(label)
    })
    
    return(unlist(labels))
}

#' Sample and Control Matching Functions
find_matching_samples <- function(sample_table,
                                directory_path,
                                config = CONFIG) {
    log_info("Finding matching samples and controls")
    
    # Setup directories
    directories <- setup_directories(directory_path)
    
    # Get matching factors
    factors <- get_factors_to_match(sample_table)
    
    # Process each sample
    results <- process_all_samples(
        sample_table,
        directories,
        factors
    )
    
    return(results)
}

setup_directories <- function(base_path) {
    list(
        bam = file.path(base_path, CONFIG$PATHS$SUBDIRS$ALIGNMENT),
        bigwig = file.path(base_path, CONFIG$PATHS$SUBDIRS$BIGWIG)
    )
}

process_all_samples <- function(sample_table,
                              directories,
                              factors) {
    log_info("Processing all samples")
    
    results <- list()
    
    for (i in seq_len(nrow(sample_table))) {
        result <- process_single_sample(
            sample_table[i, ],
            sample_table,
            directories,
            factors
        )
        
        if (!is.null(result)) {
            results[[i]] <- result
        }
    }
    
    return(results)
}

process_single_sample <- function(sample_row,
                                sample_table,
                                directories,
                                factors) {
    log_info("Processing sample:", sample_row$short_name)
    
    # Find control
    control_info <- find_control_sample(
        sample_row,
        sample_table,
        directories$bam,
        factors
    )
    
    if (is.null(control_info)) {
        log_warning("No valid control found for:", sample_row$short_name)
        return(NULL)
    }
    
    # Find sample BAM
    sample_bam <- find_sample_bam(
        sample_row,
        directories$bam
    )
    
    if (is.null(sample_bam)) {
        log_warning("No BAM file found for:", sample_row$short_name)
        return(NULL)
    }
    
    list(
        sample = sample_bam,
        control = control_info$file
    )
}

#' Sample Table Processing Functions
process_control_factors <- function(sample_table) {
    log_info("Processing control factors")
    
    cf_cols <- grep(CONFIG$PATTERNS$CONTROL_FACTOR_PREFIX, names(sample_table), 
                   value = TRUE)
    
    if (length(cf_cols) == 0) {
        log_error("No control factor columns found")
        stop("Invalid sample table format")
    }
    
    control_factors <- lapply(sample_table[cf_cols], function(x) {
        strsplit(x[1], ",")[[1]]
    })
    
    names(control_factors) <- sub(CONFIG$PATTERNS$CONTROL_FACTOR_PREFIX, "", 
                                cf_cols)
    
    sample_table[cf_cols] <- NULL
    attr(sample_table, "control_factors") <- control_factors
    
    return(sample_table)
}

#' Control Sample Matching Functions
determine_matching_control <- function(sample_row, sample_table, factors_to_match) {
    log_info("Finding matching control sample")
    
    comparison_row <- sample_row[factors_to_match]
    matches <- apply(sample_table[, factors_to_match], 1, function(row) {
        all(row == comparison_row)
    })
    
    control_indices <- which(matches & sample_table$antibody == "Input")
    
    if (length(control_indices) == 0) {
        log_warning("No matching control found")
        return(1)
    }
    
    if (length(control_indices) > 1) {
        log_warning("Multiple controls found, using first")
    }
    
    return(control_indices[1])
}

#' Script Loading Functions
load_project_scripts <- function(config = CONFIG$INITIALIZATION) {
    log_info("Loading project scripts")
    
    # Load in specific order
    load_priority_scripts(config)
    load_remaining_scripts(config)
    
    log_info("Script loading completed")
}

load_priority_scripts <- function(config) {
    log_info("Loading priority scripts")
    
    base_dir <- get_project_root()
    
    for (script in config$LOAD_ORDER$PRIORITY) {
        script_path <- find_script(script, base_dir)
        if (!is.null(script_path)) {
            source_script_safely(script_path)
        }
    }
}

load_remaining_scripts <- function(config) {
    log_info("Loading remaining scripts")
    
    # Get all script directories
    dirs <- get_script_directories(config)
    
    # Load scripts from each directory
    for (dir in dirs) {
        load_directory_scripts(dir, config)
    }
}

source_script_safely <- function(script_path) {
    log_info("Sourcing:", basename(script_path))
    
    tryCatch({
        source(script_path)
        TRUE
    }, error = function(e) {
        log_error("Failed to source:", script_path)
        log_error("Error:", e$message)
        FALSE
    })
}
#' Script Loading Functions
load_directory_scripts <- function(dir_path) {
    log_info("Loading scripts from:", dir_path)
    
    if (!dir.exists(dir_path)) {
        log_warning("Directory not found:", dir_path)
        return(FALSE)
    }
    
    # Get all R scripts
    scripts <- list.files(
        dir_path,
        pattern = CONFIG$PATTERNS$R_FILES,
        full.names = TRUE
    )
    
    # Filter excluded patterns
    scripts <- filter_scripts(scripts)
    
    # Load each script
    for (script in scripts) {
        source_script(script)
    }
    
    TRUE
}

filter_scripts <- function(scripts) {
    # Remove excluded patterns
    for (pattern in CONFIG$PATTERNS$EXCLUDE) {
        scripts <- scripts[!grepl(pattern, basename(scripts))]
    }
    scripts
}

source_script <- function(script_path) {
    log_info("Sourcing:", basename(script_path))
    
    tryCatch({
        source(script_path)
        TRUE
    }, error = function(e) {
        log_error("Failed to source:", script_path)
        log_error("Error:", e$message)
        FALSE
    })
}

load_priority_scripts <- function() {
    log_info("Loading priority scripts")
    
    for (script in CONFIG$LOAD_ORDER$PRIORITY) {
        script_path <- file.path(CONFIG$PATHS$FUNCTIONS, script)
        if (file.exists(script_path)) {
            source_script(script_path)
        }
    }
}

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

#' Track Assembly Functions
create_track_set <- function(sample_data,
                           control_data,
                           annotation_data,
                           highlight_data,
                           chromosome,
                           config = CONFIG$VISUALIZATION) {
    log_info("Creating track set")
    
    # Create base tracks
    tracks <- list(
        create_genome_axis_track(chromosome),
        create_control_track(control_data, chromosome),
        create_sample_track(sample_data, chromosome),
        create_annotation_track(annotation_data, chromosome)
    )
    
    # Add highlights if present
    if (!is.null(highlight_data)) {
        tracks <- create_highlight_track(
            tracks,
            highlight_data,
            chromosome
        )
    }
    
    return(tracks)
}

create_sample_track <- function(data,
                              chromosome,
                              name,
                              config = CONFIG$VISUALIZATION) {
    log_info("Creating sample track:", name)
    
    DataTrack(
        data,
        type = config$TRACK_TYPES$LINE,
        name = name,
        col = config$COLORS$SAMPLE,
        chromosome = chromosome
    )
}

create_control_track <- function(data,
                               chromosome,
                               name,
                               config = CONFIG$VISUALIZATION) {
    log_info("Creating control track:", name)
    
    DataTrack(
        data,
        type = config$TRACK_TYPES$LINE,
        name = name,
        col = config$COLORS$CONTROL,
        chromosome = chromosome
    )
}

#' Track Generation Functions
create_genome_axis_track <- function(chromosome, name = NULL) {
    log_info("Creating genome axis track")
    
    GenomeAxisTrack(
        name = name %||% sprintf("Chr %s Axis", chromosome)
    )
}

create_data_track <- function(data,
                            chromosome,
                            name,
                            config = CONFIG$VISUALIZATION) {
    log_info("Creating data track:", name)
    
    if (!is(data, "GRanges")) {
        data <- convert_to_granges(data, chromosome)
    }
    
    DataTrack(
        data,
        type = config$TRACKS$TYPES$DATA,
        name = name,
        col = config$TRACKS$COLORS$PRIMARY,
        chromosome = chromosome
    )
}

create_highlight_track <- function(track_list,
                                 highlights,
                                 chromosome) {
    log_info("Creating highlight track")
    
    if (!is(highlights, "GRanges")) {
        stop("Highlights must be a GRanges object")
    }
    
    HighlightTrack(
        trackList = track_list,
        start = start(highlights),
        end = end(highlights),
        chromosome = as.character(seqnames(highlights))
    )
}
