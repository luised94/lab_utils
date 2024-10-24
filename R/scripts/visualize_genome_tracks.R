#!/usr/bin/env Rscript

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
