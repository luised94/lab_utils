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
