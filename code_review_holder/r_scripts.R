
 # May need to install run in command line because of systemfonts, textshaping, and ragg:
# sudo apt-get install libharfbuzz-dev libfribidi-dev libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libboost-all-dev build-essential
renv::install(c("tidyverse", "R.utils", "ggplot2", "BiocManager", "remotes","devtools", "svglite"))
repository <- "github::"
user <- "RGLab/"
packages <- c("RProtoBufLib", "cytolib", "flowCore")
packages_to_install <- paste0(repository, user, packages)
renv::install(c(packages_to_install, "ggcyto"))
renv::install(c("tidyverse", "R.utils", "ggplot2", "BiocManager", "remotes","devtools"))

options(repos = BiocManager::repositories())

library(flowCore)
renv::snapshot()
# Installs R libraries to library_location.
# Inside R source the file. 
# This should be setup in 001_setupR/000_installingR4.2.0.sh
#home_directory <- Sys.getenv("HOME")
#library_location <- paste(home_directory, "R/x86_64-pc-linux-gnu-library/4.2", sep = "/")

install.packeges(c("tidyverse", "R.utils", "ggplot2", "BiocManager", "renv"),
        lib = library_location)
# If there is an error because of the location of the file, you can use .libPaths() to adjust the location that will be used to install. See one line below.
#.libPaths( library_location, .libPaths())

library(BiocManager)
bioconductor_packages_to_install <- c("QuasR", "GenomicAlignments", "Gviz", "rtracklayer", "ShortRead")
BiocManager::install(bioconductor_packages_to_install, lib = library_location)
#Renv has to be installed overall before it can be used to install 
#install.packeges("renv")
renv::init(bioconductor = "3.16")
# May need to install run in command line because of systemfonts, textshaping, and ragg:
# sudo apt-get install libharfbuzz-dev libfribidi-dev libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libboost-all-dev build-essential
renv::install(c("tidyverse", "R.utils", "ggplot2", "BiocManager", "remotes","devtools"))

options(repos = BiocManager::repositories())
#Wasted a bunch of time trying to figure out how to install flowCore dependencies from source. Solved it by installing through github.
devtools::install_github("RGLab/RProtoBufLib")
devtools::install_github("RGLab/cytolib")
devtools::install_github("RGLab/flowCore")
renv::install("ggcyto")

library(flowCore)

renv::snapshot()

#library(BiocManager)
#bioconductor_packages_to_install <- c("QuasR", "GenomicAlignments", "Gviz", "rtracklayer", "ShortRead")
#BiocManager::install(bioconductor_packages_to_install
#		lib = library_location)

#!/usr/bin/env Rscript
# functions moved

source("../functions/package_manager.R")
source("../functions/data_preparer.R")
source("../functions/sync_handler.R")
source("../functions/plot_generator.R")  # From previous analysis

#' Main execution function
run_visualization <- function(directory,
                            chromosome = CONFIG$DEFAULTS$CHROMOSOME,
                            bigwig_pattern = CONFIG$DEFAULTS$BIGWIG_PATTERN) {
    log_info("Starting genome visualization")
    
    # Load required packages
    load_required_packages()
    
    # Prepare data
    data <- prepare_visualization_data(
        directory,
        chromosome
    )
    
    # Prepare feature track
    feature_track <- prepare_feature_track(
        chromosome,
        data$range
    )
    
    # Generate visualization
    plot_all_sample_tracks(
        sample_table = data$samples,
        directory_path = data$directory,
        chromosome_to_plot = chromosome,
        genomeRange_to_get = data$range,
        annotation_track = feature_track,
        highlight_gr = feature_track@range,
        pattern_for_bigwig = bigwig_pattern
    )
    
    log_info("Visualization complete")
    
    # Generate sync command
    sync_cmd <- generate_sync_command(directory)
    log_info("Sync command:", sync_cmd)
}

#' Main entry point
main <- function() {
    if (!interactive()) {
        tryCatch({
            run_visualization("240819Bel")
        }, error = function(e) {
            log_error("Visualization failed:", e$message)
            quit(status = 1)
        })
    } else {
        run_visualization(
            directory = "240819Bel",
            chromosome = 10,
            bigwig_pattern = "_bamcomp.bw"
        )
    }
}

if (!interactive()) {
    main()
}
#!/usr/bin/env Rscript

#' Load all required components
source("../functions/track_generator.R")
source("../functions/feature_processor.R")
source("../functions/sample_processor.R")
source("../functions/visualization_manager.R")

#' Main execution function with mode handling
main <- function(interactive_mode = FALSE) {
    log_info("Starting genome visualization")
    
    if (interactive_mode) {
        run_interactive_mode()
    } else {
        run_batch_mode()
    }
}

#' Batch mode execution
run_batch_mode <- function() {
    tryCatch({
        args <- commandArgs(trailingOnly = TRUE)
        process_visualization(args)
        
        # Output sync command
        log_info("To sync results:")
        log_info("rsync -nav username@domain:~/data/<dir>/plots/* /local/dir/<dir>/plots/")
        
    }, error = function(e) {
        log_error("Visualization failed:", e$message)
        quit(status = 1)
    })
}

#' Interactive mode execution
run_interactive_mode <- function() {
    tryCatch({
        # Default settings for interactive mode
        directory_path <- "240808Bel"
        chromosome <- 10
        
        # Process visualization
        process_visualization(directory_path)
        
    }, error = function(e) {
        log_error("Interactive visualization failed:", e$message)
    })
}

#' Main visualization process
process_visualization <- function(args) {
    # Setup
    setup <- initialize_visualization(args)
    
    # Load data
    data <- load_visualization_data(setup)
    
    # Generate tracks
    tracks <- generate_visualization_tracks(data)
    
    # Create plots
    create_visualization_plots(tracks, setup)
}

#' Initialize visualization environment
initialize_visualization <- function(args) {
    # Validate input
    directory <- validate_input(args)
    
    # Load required packages
    load_required_packages()
    
    # Setup visualization parameters
    list(
        directory = directory,
        chromosome = CONFIG$VISUALIZATION$DEFAULT_CHROMOSOME,
        options = list(
            ucscChromosomeNames = FALSE
        )
    )
}

#' Load required data
load_visualization_data <- function(setup) {
    # Load sample table
    sample_table <- load_sample_table(setup$directory)
    
    # Load reference genome
    ref_genome <- load_reference_genome(
        CONFIG$GENOME$DIR,
        CONFIG$GENOME$PATTERN
    )
    
    # Create genome range
    genome_range <- create_chromosome_range(ref_genome)
    
    # Load features
    features <- load_feature_data(
        setup$directory,
        genome_range
    )
    
    list(
        samples = sample_table,
        genome = ref_genome,
        range = genome_range,
        features = features
    )
}

#' Main entry point
if (!interactive()) {
    main(FALSE)
} else {
    main(TRUE)
}

#!/usr/bin/env Rscript
# functions moved

source("../functions/track_generator.R")
source("../functions/range_handler.R")

#' Generate example visualization
create_example_visualization <- function(chromosome = CONFIG$GENOME$DEFAULT_CHROMOSOME,
                                      start = 1000,
                                      end = 5000) {
    log_info("Creating example visualization")
    
    # Create example data
    highlights <- create_example_highlights(chromosome, start, end)
    tracks <- create_example_tracks(chromosome, start, end)
    
    # Create highlight track
    highlight_track <- create_highlight_track(
        tracks,
        highlights,
        chromosome
    )
    
    # Generate plot
    plotTracks(
        highlight_track,
        from = start,
        to = end,
        chromosome = chromosome
    )
}

#' Example data generation functions
create_example_highlights <- function(chromosome, start, end) {
    GRanges(
        seqnames = chromosome,
        ranges = IRanges(
            start = c(start + 500, start + 2000),
            end = c(start + 1000, start + 2500)
        )
    )
}

create_example_tracks <- function(chromosome, start, end) {
    length <- end - start + 1
    
    list(
        create_genome_axis_track(chromosome),
        create_data_track(runif(length), chromosome, "Random"),
        create_data_track(rnorm(length), chromosome, "Normal"),
        create_annotation_track(chromosome, start, end)
    )
}

#!/usr/bin/env Rscript
# functions moved

source("../functions/feature_processor.R")
source("../functions/data_converter.R")

#' Main feature processing function
main <- function() {
    log_info("Starting feature processing")
    
    # Process features
    results <- process_feature_files()
    
    if (is.null(results)) {
        log_error("Feature processing failed")
        quit(status = 1)
    }
    
    # Output transfer command
    log_info("Processing complete")
    log_info("To transfer files:")
    log_info("scp -r user@server:from_dir to_dir")
}

if (!interactive()) {
    main()
}

#!/usr/bin/env Rscript
# functions moved

source("../functions/track_manager.R")
source("../functions/bigwig_processor.R")
source("../functions/control_handler.R")
source("../functions/plot_generator.R")

#' Main plotting function with improved structure
plot_all_sample_tracks <- function(sample_table,
                                 directory_path,
                                 chromosome_to_plot = 10,
                                 genomeRange_to_get,
                                 annotation_track,
                                 highlight_gr = NULL,
                                 pattern_for_bigwig = CONFIG$PATTERNS$DEFAULT_BIGWIG) {
    # ... (implementation using the above functions)
}

#!/usr/bin/env Rscript
# functions moved

source("../functions/package_installer.R")
source("../functions/environment_manager.R")

#' Main installation function
install_all_packages <- function(config = CONFIG) {
    log_info("Starting package installation")
    
    # Clean environment
    clean_environment()
    
    # Install bioinformatics packages
    for (group in names(config$BIOINFORMATICS)) {
        log_info("Installing", group, "packages")
        install_package_group(
            config$BIOINFORMATICS[[group]],
            "BiocManager"
        )
    }
    
    # Install development packages
    for (group in names(config$DEVELOPMENT)) {
        log_info("Installing", group, "packages")
        install_package_group(
            config$DEVELOPMENT[[group]],
            "renv"
        )
    }
    
    # Install GitHub packages
    install_package_group(
        config$GITHUB_PACKAGES,
        "github"
    )
    
    # Create snapshot
    renv::snapshot()
    
    log_info("Package installation completed")
}

if (!interactive()) {
    install_all_packages()
}
#!/usr/bin/env Rscript
# functions moved

source("../functions/utils.R")
source("../functions/sample_processor.R")
source("../functions/genome_processor.R")

main <- function() {
    log_info("Starting genome track plotting")
    
    # Load required packages
    load_required_packages()
    
    # Process command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    directory_path <- validate_input(args)
    
    # Load and process sample table
    sample_table <- load_sample_table(directory_path)
    
    # Set up genome visualization
    options(ucscChromosomeNames = FALSE)
    refGenome <- load_reference_genome()
    
    # Create genome ranges
    genomeRange <- create_chromosome_GRange(refGenome)
    
    # Additional visualization logic here
    
    log_info("Genome track plotting completed")
}

if (!interactive()) {
    main()
}

#!/usr/bin/env Rscript
# functions moved

source("../functions/sample_matcher.R")
source("../functions/bam_finder.R")

#' Main sample input determination function
determine_input_for_all_samples <- function(sample_table,
                                          directory_path,
                                          reference_pattern = CONFIG$FILES$PATTERNS$REFERENCE) {
    log_info("Starting sample input determination")
    
    # Validate inputs
    validate_inputs(sample_table, directory_path)
    
    # Find matching samples
    results <- find_matching_samples(
        sample_table,
        directory_path
    )
    
    # Output results
    output_results(results)
    
    log_info("Sample input determination completed")
}

validate_inputs <- function(sample_table, directory_path) {
    if (!all(CONFIG$VALIDATION$REQUIRED_COLUMNS %in% colnames(sample_table))) {
        log_error("Missing required columns in sample table")
        stop("Invalid sample table")
    }
    
    if (!dir.exists(directory_path)) {
        log_error("Directory not found:", directory_path)
        stop("Invalid directory")
    }
}

output_results <- function(results) {
    for (result in results) {
        if (!is.null(result)) {
            cat(result$sample, result$control, sep = "\n")
        }
    }
}

#!/usr/bin/env Rscript

source("../functions/input_validator.R")

#' Main plotting function
main <- function() {
    # Get command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    
    # Validate inputs
    tryCatch({
        validate_script_args(
            directory = args[1],
            chromosome = if(length(args) > 1) as.integer(args[2]) else NULL
        )
    }, error = function(e) {
        log_error("Validation failed:", e$message)
        quit(status = 1)
    })
    
    # Continue with script execution...
}

if (!interactive()) {
    main()
}

#!/usr/bin/env Rscript
# functions moved

source("../functions/data_downloader.R")
source("../functions/feature_processor.R")

#' Main feature download function
main <- function() {
    log_info("Starting feature download")
    
    # Setup
    timestamp <- format(Sys.time(), "%y%m%d")
    base_dir <- file.path(Sys.getenv("HOME"), "data", "feature_files")
    
    ensure_directory(base_dir)
    
    # Download and process data
    results <- download_feature_data(
        CONFIG$SOURCES,
        base_dir,
        timestamp
    )
    
    # Validate results
    validate_results(results)
    
    log_info("Feature download complete")
    log_info("Cleanup hint: rm *_eaton_acs.bed if not needed")
}

if (!interactive()) {
    main()
}

#!/usr/bin/env Rscript
# functions moved

source("../functions/file_validator.R")
source("../functions/control_manager.R")

#' Main control determination function
determine_sample_controls <- function(directory,
                                    task_id,
                                    reference_pattern = CONFIG$PATTERNS$REFERENCE_GENOME) {
    log_info("Starting control determination")
    
    # Load and validate sample table
    sample_table <- load_sample_table(directory)
    
    # Get matching factors
    factors <- get_factors_to_match(sample_table)
    
    # Setup directories
    bam_dir <- file.path(directory, CONFIG$PATHS$SUBDIRS$ALIGNMENT)
    
    # Get sample BAM
    sample_bam <- find_matching_bam(
        sample_table$sample_ID[task_id],
        bam_dir
    )
    
    if (is.null(sample_bam)) {
        log_error("Sample BAM not found")
        return(NULL)
    }
    
    # Find control
    control_bam <- find_valid_control(
        sample_table[task_id, ],
        sample_table,
        bam_dir,
        factors
    )
    
    if (is.null(control_bam)) {
        log_error("No valid control found")
        return(NULL)
    }
    
    # Return paths
    cat(sample_bam, control_bam, sep = "\n")
}

#' Main entry point
main <- function() {
    if (!interactive()) {
        args <- commandArgs(trailingOnly = TRUE)
        tryCatch({
            arg_list <- validate_input(args)
            determine_sample_controls(
                arg_list$directory_path,
                arg_list$slurm_array_task_id
            )
        }, error = function(e) {
            log_error("Execution failed:", e$message)
            quit(status = 1)
        })
    }
}

if (!interactive()) {
    main()
}

#!/usr/bin/env Rscript
# functions moved

source("../functions/bam_qc_analyzer.R")
source("../functions/mapping_calculator.R")

#' Main BAM QC analysis function
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    
    if (length(args) == 0) {
        log_error("No directory specified")
        stop("Usage: Rscript analyze_bam_qc.R <directory>")
    }
    
    # Find directory
    directory <- find_data_directory(args[1])
    
    if (is.null(directory)) {
        log_error("Directory not found:", args[1])
        stop("Invalid directory")
    }
    
    # Analyze QC
    results <- analyze_bam_qc(directory)
    
    if (!is.null(results)) {
        print(results)
    }
}

if (!interactive()) {
    main()
}

# functions moved
# revisit
#library(flowCore)
#library(tidyverse)
#library(svglite)
#library(ggplot2)
#library(ggridges)
#
#path_to_fcs_files <- "/mnt/c/Users/Luis/Dropbox (MIT)/Lab/Experiments/Yeast Genetics/2023_10_04 Cdc6 Overexpression/CL4_Cdc6Oe_CarbonSourceShift_01"
##=======
#
##' Extract FITC Width Data from an FCS File
##'
##' This function reads an FCS file and extracts the fluorescence data
##' corresponding to the FITC width of individual cells.
##'
##' @param file_path The path to the FCS file.
##' @return A data frame with one column, FITC_Width, containing the extracted data.
##' @examples
##' extract_fitc_width_from_fcs("path/to/your/file.fcs")
##' @export
#extract_fitc_width_from_fcs <- function(file_path) {
#	    if (!file.exists(file_path)) {
#		            stop("File does not exist: ", file_path)
#    }
#    
#    fcs_data <- tryCatch({
#	            read.FCS(file_path, transformation = FALSE)
#		        }, error = function(e) {
#				        stop("Failed to read FCS file: ", e$message)
#		        })
#        
#        if (!'FITC-Width' %in% colnames(exprs(fcs_data))) {
#		        stop("FITC-Width data not found in the file.")
#	    }
#        
#        data_frame <- data.frame(FITC_Width = exprs(fcs_data[,'FITC-Width']))
#	    return(data_frame)
#}
#
##' Generate a Current Datetime String for File Naming
##'
##' This function returns a string representing the current date and time,
##' formatted for use in file names. The format used is "YYYY-MM-DD-HH-MM-SS",
##' and the system's local time zone is assumed.
##'
##' @return A character string of the current date and time.
##' @examples
##' get_current_datetime_string()
##' @export
#get_current_datetime_string <- function() {
#	          return(format(Sys.time(), "%Y-%m-%d-%H-%M-%S"))
#}
######
#
#windows_user <- list.files("/mnt/c/Users")[grepl(Sys.info()[["user"]], list.files("/mnt/c/Users"), ignore.case = TRUE)]
#path_to_fcs_files <- sprintf("/mnt/c/Users/%s/Dropbox (MIT)/Lab/Experiments/Yeast Genetics/2023_10_04 Cdc6 Overexpression/CL4_Cdc6Oe_CarbonSourceShift_01", windows_user)
#fcs_file_paths <- list.files(path_to_fcs_files, pattern = "\\.fcs$", full.names = TRUE)
##fcs_data <- read.FCS(fcs_file_paths[1])
#
##TODO: Initialize variables. for each one, use length to assign module class. then maye use switch or some other tidyverse to replace values with variables.
##REMINDER: Calculated the types of samples that would account for my FACS samples. Seems like I had 4 genotypes, 2 carbon sources and 7 timepoints.
#genotypes <- c("ORC4", "orc4-R267A", "ORC1", "orc1-K485T")
#carbon_sources <- c("GLU", "GAL")
#time_points <- c(seq(0, 150, by = 30))
#df_sample_info <- as.data.frame(expand.grid(genotype_Orc4 = genotypes, 
#					    carbonSources = carbon_sources, 
#					    timePoints = time_points)) %>% arrange(genotype_Orc4,carbonSources,timePoints)
#df_sample_info$filePath <- gtools::mixedsort(fcs_file_paths)#
#is_WTandGLU <- df_sample_info$genotype_Orc4 == "ORC4" & df_sample_info$carbonSources == "GLU"
#subset_df <- df_sample_info %>% filter(is_WTandGLU)
#for (i in 1:nrow(subset_df)) {
#
##Check the first six to ensure creation went well. 
#for (i in 1:6) {
##	flow_data <- read.FCS(df_sample_info$filePath[i])
##	exprs(flow_data[,'FITC-Width'])
#	formatted_output <- sprintf("%s | %s ", i, basename(subset_df$filePath[i]))
#	print(formatted_output)
#}
#
#subset_df$fitcData <- map(subset_df$filePath, extract_fitc_width_from_fcs)
#combined_data <- bind_rows(subset_df$fitcData, .id = "sample_id")
#
#head(combined_data)
#dim(combined_data)
#str(combined_data)
#summary(combined_data)
##Plot the data on the same plot and use color gradient to see. 
#output_directory <- sprintf("/mnt/c/Users/%s/Dropbox (MIT)/Lab/Experiments/Yeast Genetics/2023_10_04 Cdc6 Overexpression/", windows_user)
#plot_output <- paste0(output_directory, get_current_datetime_string(), "_densityPlot.svg")
##svglite(plot_output, width = 6, height = 4)
#ggplot(combined_data, aes(x = `FITC.Width`, y = after_stat(density), fill = sample_id)) +
#	geom_density(alpha = .2) +
#	xlim(0, median(combined_data$FITC.Width) + 1500) + 
#	labs(title = "Fluorescence Intensity Distribution", x = "Fluorescence Intensity", y = "Density")
##dev.off()
#
#
## Plot a staggered version using an offset of the samples on the same plot. 
#combined_data$offset <- as.numeric(as.factor(combined_data$sample_id)) * 100 # Create the offset by multiplying sample id with 10
#plot_output <- paste0(output_directory, get_current_datetime_string(), "_plot_count.svg")
