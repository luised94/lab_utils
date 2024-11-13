#!/usr/bin/env Rscript
# Configuration
#-----------------------------------------------------------------------------
DEBUG_CONFIG <- list(
    enabled = FALSE,           # TRUE for testing single group, FALSE for all
    group = 10,               # Which group to process when in debug mode
    samples_per_group = 4,    # Samples per plot
    save_plots = FALSE,       # Whether to save plots to files
    verbose = TRUE,           # Print debug information
    chromosome = 10,
    interactive = TRUE,
    display_time = 2
)


# Time formatting configuration
TIME_CONFIG <- list(
    timestamp_format = "%Y%m%d_%H%M%S",  # YYYYMMDD_HHMMSS
    date_format = "%Y%m%d"               # YYYYMMDD
)

# Generate timestamps once at script start
TIMESTAMPS <- list(
    full = format(Sys.time(), TIME_CONFIG$timestamp_format),
    date = format(Sys.Date(), TIME_CONFIG$date_format)
)

PLOT_CONFIG <- list(
    width = 10,
    height = 8,
    track_color = "#fd0036",
    placeholder_color = "#cccccc",
    track_name_format = "%s: %s - %s",
    placeholder_suffix = "(No data)",
    background_colors = list(
        title = "white",
        border =  "#E0E0E0"
    )

)

PLOT_CONFIG$title_format <- list(
    main = "%s\nChromosome %s (%d samples)\n%s\nNormalization: %s",
    background = "#F0F0F0"
)

# Load required packages
#-----------------------------------------------------------------------------
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
}

#source("~/lab_utils/failsafe_scripts/all_functions.R")
source("~/lab_utils/failsafe_scripts/bmc_config.R")
# Load metadata and files
#-----------------------------------------------------------------------------
experiment_id <- "241007Bel"
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
plots_dir <- file.path(base_dir, "plots", "genome_tracks", "overview")
metadata_path <- file.path(base_dir, "documentation", 
                          paste0(experiment_id, "_sample_grid.csv"))
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Find fastq files and extract sample IDs
fastq_files <- list.files(
    path = file.path(base_dir, "fastq"),
    pattern = "consolidated_.*_sequence\\.fastq$",
    full.names = FALSE
)

if (length(fastq_files) == 0) {
    stop("No fastq files found in specified directory")
}

# Extract sample IDs from fastq filenames
sample_ids <- gsub(
    pattern = "consolidated_([0-9]{5,6})_sequence\\.fastq",
    replacement = "\\1",
    x = fastq_files
)

# 3. Load and process metadata
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

# 4. Enforce factor levels from config
for (col_name in names(EXPERIMENT_CONFIG$CATEGORIES)) {
    if (col_name %in% colnames(metadata)) {
        metadata[[col_name]] <- factor(
            metadata[[col_name]],
            levels = EXPERIMENT_CONFIG$CATEGORIES[[col_name]],
            ordered = TRUE
        )
    }
}

# 5. Sort metadata using config column order
sorted_metadata <- metadata[do.call(
    order, 
    metadata[EXPERIMENT_CONFIG$COLUMN_ORDER]
), ]

# 6. Add sample IDs to metadata
sorted_metadata$sample_id <- sample_ids

# Function to generate distinct colors
generate_distinct_colors <- function(n) {
    # RColorBrewer provides good distinct colors
    if (n <= 8) {
        RColorBrewer::brewer.pal(max(3, n), "Set2")[1:n]
    } else {
        # For more categories, use rainbow with better spacing
        rainbow(n, s = 0.7, v = 0.9)
    }
}

# Create color mapping for antibodies
unique_antibodies <- unique(sorted_metadata$antibody)
antibody_colors <- generate_distinct_colors(length(unique_antibodies))
names(antibody_colors) <- unique_antibodies

# Update PLOT_CONFIG with dynamic colors
PLOT_CONFIG$track_colors <- list(
    antibody = antibody_colors,
    placeholder = PLOT_CONFIG$placeholder_color  # Maintain consistent placeholder
)

# Create color legend text
legend_text <- sprintf(
    "Track Colors:\n%s\n%s",
    paste("?", names(PLOT_CONFIG$track_colors$antibody), 
          sprintf("(%s)", PLOT_CONFIG$track_colors$antibody), 
          collapse = "\n"),
    sprintf("? No Data (%s)", PLOT_CONFIG$placeholder_color)
)

if (DEBUG_CONFIG$verbose) {
    message("Metadata processing summary:")
    message(sprintf("Found %d fastq files", length(fastq_files)))
    message(sprintf("Processed %d metadata rows", nrow(sorted_metadata)))
    message("Columns with enforced factors:")
    print(names(EXPERIMENT_CONFIG$CATEGORIES))
}

# Load reference genome
ref_genome_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    pattern = "S288C_refgenome.fna",
    full.names = TRUE,
    recursive = TRUE
)[1]
genome_data <- Biostrings::readDNAStringSet(ref_genome_file)


# Create chromosome range
chromosome_to_plot <- DEBUG_CONFIG$chromosome
chromosome_width <- genome_data[chromosome_to_plot]@ranges@width
chromosome_roman <- paste0("chr", utils::as.roman(chromosome_to_plot))

genome_range <- GenomicRanges::GRanges(
    seqnames = chromosome_roman,
    ranges = IRanges::IRanges(start = 1, end = chromosome_width),
    strand = "*"
)

# Find bigwig files
bigwig_files <- list.files(
    file.path(base_dir, "coverage"),
    pattern = "_CPM\\.bw$",
    full.names = TRUE
)
normalization_method <- sub(".*_([^_]+)\\.bw$", "\\1", 
                          basename(bigwig_files[1]))

# Load feature file (annotation)
feature_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "feature_files"),
    pattern = "eaton_peaks",
    full.names = TRUE
)[1]

if (!is.null(feature_file)) {
    features <- rtracklayer::import(feature_file)
    # Convert to chrRoman format
    GenomeInfoDb::seqlevels(features) <- paste0(
        "chr", 
        utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(features)))
    )
}

# Process samples in groups
#-----------------------------------------------------------------------------
# Create sample groups
sample_groups <- split(
    seq_len(nrow(sorted_metadata)),
    ceiling(seq_len(nrow(sorted_metadata)) / DEBUG_CONFIG$samples_per_group)
)

# Determine which groups to process
groups_to_process <- if (DEBUG_CONFIG$enabled) {
    DEBUG_CONFIG$group
} else {
    seq_along(sample_groups)
}

# Calculate global y-limits for all plots (before the plotting loop)
#-----------------------------------------------------------------------------
if (DEBUG_CONFIG$verbose) {
    message("\nCalculating global range for all tracks...")
}

# Initialize vectors to store all min/max values
all_track_values <- c()

# Process each bigwig file
for (bigwig_file in bigwig_files) {
    if (file.exists(bigwig_file)) {
        tryCatch({
            # Import data for specific chromosome
            track_data <- rtracklayer::import(
                bigwig_file,
                which = genome_range
            )
            
            if (length(track_data) > 0) {
                values <- GenomicRanges::values(track_data)$score
                if (length(values) > 0) {
                    all_track_values <- c(all_track_values, values)
                }
            }
        }, error = function(e) {
            if (DEBUG_CONFIG$verbose) {
                message("Skipping ", basename(bigwig_file), ": ", e$message)
            }
        })
    }
}

# Calculate global limits with 10% padding
if (length(all_track_values) > 0) {
    y_min <- min(all_track_values, na.rm = TRUE)
    y_max <- max(all_track_values, na.rm = TRUE)
    y_range <- y_max - y_min
    y_limits <- c(
        y_min - (y_range * 0.1),  # Add 10% padding
        y_max + (y_range * 0.1)
    )
    
    if (DEBUG_CONFIG$verbose) {
        message(sprintf("Global y-limits: [%.2f, %.2f]", y_limits[1], y_limits[2]))
    }
} else {
    y_limits <- NULL
    if (DEBUG_CONFIG$verbose) {
        message("No valid track data found for y-limit calculation")
    }
}

# Process each group
for (group_idx in groups_to_process) {
    if (DEBUG_CONFIG$verbose) {
        message("\nProcessing group ", group_idx)
    }
    
    # Get current group's samples
    current_samples <- sorted_metadata[sample_groups[[group_idx]], ]
    
    # Initialize tracks list with chromosome axis

    tracks <- list(
        Gviz::GenomeAxisTrack(
            name = paste("Chr ", chromosome_to_plot, " Axis", sep = "")
        )
    )
    
    # Add tracks for each sample
    for (i in seq_len(nrow(current_samples))) {
        sample_id <- current_samples$sample_id[i]
        
        # Find matching bigwig file
        bigwig_file <- bigwig_files[grepl(sample_id, bigwig_files)][1]

        track_name <- sprintf(
            PLOT_CONFIG$track_name_format,
            sample_id,
            current_samples$short_name[i],
            current_samples$antibody[i]
        )

        placeholder_name <- sprintf(
            "%s %s",
            track_name,
            PLOT_CONFIG$placeholder_suffix
        )

        if (!is.na(bigwig_file) && file.exists(bigwig_file)) {
            if (DEBUG_CONFIG$verbose) {
                message("Adding track for sample: ", sample_id)
            }
            
            track_data <- rtracklayer::import(bigwig_file)
            
            tracks[[length(tracks) + 1]] <- Gviz::DataTrack(
                track_data,
                name = track_name,
                type = "l",
                col = PLOT_CONFIG$track_color
            )
        } else {
            if (DEBUG_CONFIG$verbose) {
                message("Creating placeholder for sample: ", sample_id)
            }
            
            # Calculate reasonable number of points (e.g., one point per 100bp)
            sampling_rate <- 100  # Adjust this value based on your needs
            num_points <- ceiling(chromosome_width / sampling_rate)
            
            # Create proper placeholder track with genome coordinates
            empty_ranges <- GenomicRanges::GRanges(
                seqnames = chromosome_roman,
                ranges = IRanges::IRanges(
                    # Create evenly spaced points across chromosome
                    start = seq(1, chromosome_width, length.out = num_points),
                    width = 1
                ),
                score = rep(0, num_points),  # One score per position
                strand = "*"
            )
            
            empty_track <- Gviz::DataTrack(
                empty_ranges,
                name = placeholder_name,
                type = "l",
                col = PLOT_CONFIG$placeholder_color,
                chromosome = chromosome_roman
            )
            tracks[[length(tracks) + 1]] <- empty_track
        }
    }

    # Add feature track if available
    if (exists("features")) {
        tracks[[length(tracks) + 1]] <- Gviz::AnnotationTrack(
            features,
            name = "Features"
        )
    }

    plot_title <- sprintf(
        PLOT_CONFIG$title_format$main,
        experiment_id,
        chromosome_to_plot,
        nrow(current_samples),
        TIMESTAMPS$full,
        normalization_method
    )
    
    Gviz::plotTracks(
        trackList = tracks,
        chromosome = chromosome_roman,
        from = 1,
        to = chromosome_width,
        ylim = y_limits,
        main = plot_title,
        cex.main = 1,
        background.title = PLOT_CONFIG$background_colors$title,
        col.border.title = PLOT_CONFIG$background_colors$border,

        col.main = "black",
        main.width = 0.8,
        legend = legend_text,
        cex.legend = 0.8,
        margin = 15,        # Increase margin for readability
        innerMargin = 5    # Space between tracks
    )
    
    # Save plot if needed
    if (DEBUG_CONFIG$save_plots) {
        plot_filename <- sprintf(
                "%s_%s_chr%s_n%d_group%d.svg",
                TIMESTAMPS$full,
                experiment_id,
                chromosome_to_plot,
                nrow(current_samples),
                group_idx
        )
        plot_file <- file.path(
            plots_dir,
            plot_filename
        )
        
        if (DEBUG_CONFIG$verbose) {
            message(sprintf(
                "Saving plot to: %s", 
                basename(plot_file)
            ))
        }
        
        svg(plot_file, width = PLOT_CONFIG$width, height = PLOT_CONFIG$height)
        Gviz::plotTracks(
            trackList = tracks,
            chromosome = chromosome_roman,
            from = 1,
            to = chromosome_width,
            ylim = y_limits,
            main = plot_title,
            cex.main = 1,
            col.main = "black",
            margin = 10,        # Increase margin for readability
            innerMargin = 5,    # Space between tracks
            background.title = "#F0F0F0"  # Light gray background for track titles
        )
        dev.off()
    }
    
    # Interactive viewing options
    if (DEBUG_CONFIG$interactive) {
        user_input <- readline(
            prompt = "Options: [Enter] next plot, 's' skip rest, 'q' quit: "
        )
        if (user_input == "q") break
        if (user_input == "s") DEBUG_CONFIG$save_plots <- FALSE
    } else {
        Sys.sleep(DEBUG_CONFIG$display_time)  # Pause between plots
    }
}

message("Processing complete")
