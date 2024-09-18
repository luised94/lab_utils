#DESCRIPTION: 
#USAGE:
#TODO: Need some way to aggregate the column names I use from the sampleConfig.R files to keep track of variables I use to keep consistent experiment to experiment.
main <- function() {
    #Load packages
    suppressPackageStartupMessages({
        library(QuasR)
        library(GenomicAlignments)
        library(Gviz)
        library(rtracklayer)
        library(ShortRead)
        library(tidyverse)
    })
    args <- commandArgs(trailingOnly = TRUE)
    directory_path <- validate_input(args)
    #Add process_control_factors, get_factors_to_match
    sample_table <- load_sample_table(directory_path)
    determine_input_for_all_samples(sample_table, directory_path)
    #determineInput
    #chromosome_to_plot = 10
    #options(ucscChromosomeNames=FALSE)
    #refGenome <- load_reference_genome(genome_dir = "REFGENS", genome_pattern = "S288C_refgenome.fna")
    #genomeRange_to_get <- create_chromosome_GRange(refGenome = refGenome)
    #feature_file_pattern = "eaton_peaks"
    #feature_grange <- load_feature_file_GRange(chromosome_to_plot = chromosome_to_plot, feature_file_pattern = feature_file_pattern,
    #                                         genomeRange_to_get = genomeRange_to_get)
    #control_dir <- "EatonBel"
    #control_track <- load_control_grange_data(control_dir = control_dir, chromosome_to_plot = chromosome_to_plot,
    #                                         genomeRange_to_get = genomeRange_to_get)
    # Convert to input, add determine_matching_control, select_control_index
    #control_dir <- "EatonBel"
    #control_track <- load_control_grange_data(control_dir = control_dir, chromosome_to_plot = chromosome_to_plot,
                #                             genomeRange_to_get = genomeRange_to_get)
    #plot_all_sample_tracks(sample_table = sample_table,
    #                       directory_name = directory_path,
    #                       chromosome_to_plot = chromosome_to_plot, 
    #                       genomeRange_to_get = genomeRange_to_get, 
    #                       control_track = control_track, 
    #                       annotation_track = origin_track)
}

validate_input <- function(args) {
    if (length(args) != 1) {
        cat("Error: Invalid number of arguments.\n")
        cat("Usage: Rscript 001_plotAllSampleTracks.R <directory_path>\n")
        cat("Example: Rscript 001_plotAllSampleTracks.R 240819Bel\n")
        stop()
    }
    directory_path <- file.path(Sys.getenv("HOME"), "data", args[1])
    if(!dir.exists(directory_path)) {
        cat(sprintf("Error: Directory %s does not exist.\n", directory_path))
        stop()
    }
    return(directory_path)
}

process_control_factors <- function(sample_table) {
    cat("Process control factors from __cf_ columns\n")
    df <- sample_table
    cf_cols <- grep("X__cf_", names(df), value = TRUE)
    if(length(cf_cols) == 0) {
        cat("No columns containing __cf_ tag found in sample table")
        stop("Verify sample table was produced with updated sampleGridConfig.")
    }
    control_factors <- lapply(df[cf_cols], function(x) strsplit(x[1], ",")[[1]])
    names(control_factors) <- sub("X__cf_", "", cf_cols)
    df[cf_cols] <- NULL
    attr(df, "control_factors") <- control_factors
    return(df)
}

get_factors_to_match <- function(sample_table, attribute_to_get = "control_factors") {
    cat("Grabbing attributes from sample table\n")
    df <- sample_table
    control_factors <- attr(df, attribute_to_get)
    if (is.null(control_factors)) {
        stop("No control factors defined in sample data.\nVerify 003_updateSampleGrid.R")
    }
    all_factors <- unlist(control_factors)
    return(intersect(all_factors, colnames(df)))
}

determine_matching_control <- function(sample_row, sample_table, factors_to_match) {
    cat("Determining control row for sample row.\n")
    df <- sample_table
    comparison_row <- sample_row[factors_to_match]
    rows_with_same_factors <- apply(df[, factors_to_match], 1, function(row) {
        all(row == comparison_row)
    })
    is_input <- df$antibody == "Input"
    index <- as.numeric(unname(which(is_input & rows_with_same_factors)))
    return(index)
}

select_control_index <- function(control_indices, max_controls = 1) {
    cat("Processing control index to ensure one is used.\n")
    if (length(control_indices) == 0) {
        warning("No matching control found")
        cat("Setting control_index to 1\n")
        control_indices <- 1
    }
    if (length(control_indices) > max_controls) {
    warning(paste("Multiple matching controls found, using first", max_controls))
    control_indices[1:max_controls]
    } else if (length(control_indices) == 1){
        return(control_indices)
    }
}

load_sample_table <- function(directory_name) {
    cat("Loading sample_table from", directory_name, "\n")
    documentation_dir_path <- file.path(directory_name, "documentation")
    sample_table_path <- list.files(documentation_dir_path, pattern = "sample_table", full.names = TRUE)
    if(length(sample_table_path) == 0){
        cat(sprintf("No files with pattern sample_table found in %s\n", documentation_dir_path))
        cat("Consult 000_setupExperimentDir to create sample_table.tsv for the project\n")
        stop()
    } else if(length(sample_table_path) > 1){
        cat(sprintf("Multiple files with pattern sample_table found in %s\n", documentation_dir_path))
        cat(sprintf("Files found in %s\n", documentation_dir_path))
        print(sample_table_path)
        cat("Consult 000_setupExperimentDir and ensure no duplicates are present for the project\n")
        stop()
    }
    cat(sprintf("Reading %s\n", sample_table_path))
    sample_table <- read.delim(sample_table_path, header = TRUE, sep ="\t")
    if (!("sample_ID" %in% colnames(sample_table))) {
           cat("Sample table does not contain the sample_ID column.\n")
           cat("sample_ID column is required for analysis.\n")
           cat("See 000_setupExperimentDir.R and 003_updateSampleGrid.R.\n")
           stop()
    }
    sample_table <- process_control_factors(sample_table)
    cat("Head of sample_table\n")
    print(head(sample_table))
    return(sample_table)
}

unique_labeling <- function(table, categories_for_label) {
    # Input validation
    if (!is.data.frame(table)) {
        stop("Input 'table' must be a data frame")
    }
    if (!is.character(categories_for_label) || length(categories_for_label) == 0) {
        stop("Input 'categories_for_label' must be a non-empty character vector")
    }
    
    # Ensure antibody category is always included
    if (!"antibody" %in% categories_for_label) {
    categories_for_label <- c("antibody", categories_for_label)
    }
    
    print(paste("Categories for label:", paste(categories_for_label, collapse = ", ")))
    
    # Check if all categories exist in the table
    missing_categories <- setdiff(categories_for_label, colnames(table))
    if (length(missing_categories) > 0) {
        stop(paste("The following categories are missing from the table:", 
        paste(missing_categories, collapse = ", ")))
    }
    
    # Identify unique values for each category
    unique_values <- lapply(table[categories_for_label], unique)
    print("Unique values for each category:")
    print(unique_values)
    
    # Function to construct label for a single sample
    construct_label <- function(sample) {
    differing_categories <- sapply(categories_for_label, function(cat) {
        if (length(unique_values[[cat]]) > 1 || cat == "antibody") {
            return(sample[cat])
            #return(paste(cat, sample[cat], sep = ": "))
        } else {
            return(NULL)
        }
    })
        differing_categories <- differing_categories[!sapply(differing_categories, is.null)]
        return(paste(differing_categories, collapse = "_"))
    }
    
    # Apply the construct_label function to each sample (row)
    labels <- apply(table, 1, construct_label)
    
    print("Constructed labels:")
    print(labels)
    
    return(unlist(labels))
}
#determineInput <-  function(sample_row){}
determine_input_for_all_samples <- function(sample_table, directory_path){
    bigwig_directory <- file.path(directory_path, "bigwig")
    for (sample in sample_table) {
        print(sample) 
    }
}

if(!interactive()){
    main()
} else {
    suppressPackageStartupMessages({
        library(QuasR)
        library(GenomicAlignments)
        library(Gviz)
        library(rtracklayer)
        library(ShortRead)
        library(tidyverse)
        library(gtools)
    })
    #@update
    directory_path <- "240819Bel"
    cat("Logic of main function\n")
    main_function_logic <- main
    list_of_functions_and_variables <- ls()
    cat("Main function: Validating directory\n")
    directory_path <- validate_input(directory_path)
    cat("Main function: loading sample table\n")
    sample_table <- load_sample_table(directory_path)
    print(directory_path)
    determine_input_for_all_samples(sample_table, directory_path)
}
