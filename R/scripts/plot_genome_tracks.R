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
