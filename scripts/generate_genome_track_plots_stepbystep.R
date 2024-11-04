# 1. Load required packages with verification
packages <- c("QuasR", "GenomicAlignments", "Gviz", "rtracklayer", "ShortRead", "tidyverse")
for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        stop(paste("Package", pkg, "not found"))
    }
    message(paste("Loaded package:", pkg))
}

# 2. Set up initial variables and paths
experiment_id <- "241007Bel"
#experiment_id <- "241010Bel"
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
message("Working directory: ", base_dir)

# 3. Load and verify config
source("~/lab_utils/scripts/bmc_config.R")
source("~/lab_utils/scripts/extract_bmcIDmetadata_process.R")
source("~/lab_utils/scripts/genome_core.R")
source("~/lab_utils/scripts/sample_processing.R")
message("Loaded config with ", length(EXPERIMENT_CONFIG$COMPARISONS), " comparisons")

# Add timestamp constant
TIMESTAMP <- format(Sys.Date(), "%Y%m%d")

# 4. Load and process sample table with new processing steps
metadata_file <- file.path(
    base_dir, 
    "documentation", 
    sprintf("%s_processed_grid.csv", experiment_id)
)
message("Processing metadata file: ", metadata_file)

# First enforce factor levels
sample_table <- read.csv(metadata_file, stringsAsFactors = FALSE)
message("Initial sample table rows: ", nrow(sample_table))

# Enforce factor levels from config
sample_table <- enforce_factor_levels(
    data_frame = sample_table,
    categories = EXPERIMENT_CONFIG$CATEGORIES
)
message("Factor levels enforced")

# Sort metadata using config column order
sample_table <- sort_metadata_frame(
    data_frame = sample_table,
    column_order = EXPERIMENT_CONFIG$COLUMN_ORDER
)
message("Sample table sorted by: ", 
        paste(EXPERIMENT_CONFIG$COLUMN_ORDER, collapse = ", "))
if("sample_id" %in% names(sample_table)){
    sample_table$experiment_number <- sort(sample_table$sample_id)
} else {
    sample_table$experiment_number <- sort(sample_table$experiment_number)
}

# 5. Load reference genome
ref_genome_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    pattern = "S288C_refgenome.fna",
    full.names = TRUE,
    recursive = TRUE
)[1]
message("Found reference genome: ", ref_genome_file)
ref_genome <- readFasta(ref_genome_file)
ref_genome_df <- data.frame(
    chrom = names(as(ref_genome, "DNAStringSet")),
    basePairSize = width(ref_genome)
) %>% filter(chrom != "chrM")
message("Processed reference genome with ", nrow(ref_genome_df), " chromosomes")

# 6. Create genome ranges
genome_ranges <- GRanges(
    seqnames = ref_genome_df$chrom,
    ranges = IRanges(start = 1, end = ref_genome_df$basePairSize),
    strand = "*"
)
message("Created genome ranges object")

# 7. Load feature file
feature_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "feature_files"),
    pattern = "eaton_peaks",
    full.names = TRUE
)[1]
message("Found feature file: ", feature_file)
# Load and adjust feature chromosome style
feature_ranges <- rtracklayer::import.bed(feature_file)
feature_style <- genome_detect_chr_style(seqlevels(feature_ranges))
genome_style <- genome_detect_chr_style(seqlevels(genome_ranges))
if (feature_style != genome_style) {
    message("Adjusting feature chromosome style to match genome")
    new_seqlevels <- genome_convert_chr_names(
        seqlevels(feature_ranges),
        genome_style
    )
    seqlevels(feature_ranges) <- new_seqlevels
    feature_style <- genome_detect_chr_style(seqlevels(feature_ranges))
}

feature_track <- AnnotationTrack(
    feature_ranges,
    name = "Origin Peaks (Eaton 2010)"
)
message("Created feature track with style: ", 
        genome_detect_chr_style(seqlevels(feature_ranges)))

# Function for bigwig validation
validate_bigwig <- function(bigwig_path, experiment_number) {
    if (is.na(bigwig_path) || !file.exists(bigwig_path)) {
        warning(sprintf("Bigwig file not found for sample %s", experiment_number))
        return(NULL)
    }
    tryCatch({
        rtracklayer::import(bigwig_path)
        return(bigwig_path)
    }, error = function(e) {
        warning(sprintf("Invalid bigwig file for sample %s: %s", experiment_number, e$message))
        return(NULL)
    })
}

find_fallback_control <- function(sample_table, bigwig_dir, pattern) {
    # Get all Input samples
    input_samples <- sample_table[sample_table$antibody == "Input", ]
    
    for (i in seq_len(nrow(input_samples))) {
        control_bigwig <- list.files(
            bigwig_dir,
            pattern = paste0(input_samples$experiment_number[i], ".*normalized.*\\.bw$"),
            full.names = TRUE
        )[1]
        
        valid_control <- validate_bigwig(control_bigwig, input_samples$experiment_number[i])
        if (!is.null(valid_control)) {
            message("Using fallback control: ", basename(valid_control))
            return(valid_control)
        }
    }
    
    warning("No valid Input controls found in dataset")
    return(NULL)
}

# Setup chromosome
chromosome <- 10
chromosome_roman <- paste0("chr", as.roman(chromosome))
message("Processing chromosome: ", chromosome_roman)

# 8. Process all comparisons
calculate_track_limits <- function(tracks) {
    y_ranges <- lapply(tracks, function(track) {
        if (inherits(track, "DataTrack")) {
            values <- values(track)
            if (length(values) > 0) {
                return(range(values, na.rm = TRUE))
            }
        }
        return(NULL)
    })
    
    # Remove NULL entries and calculate overall range
    y_ranges <- y_ranges[!sapply(y_ranges, is.null)]
    if (length(y_ranges) > 0) {
        y_min <- min(sapply(y_ranges, `[`, 1), na.rm = TRUE)
        y_max <- max(sapply(y_ranges, `[`, 2), na.rm = TRUE)
        
        # Add 10% padding
        y_range <- y_max - y_min
        y_min <- y_min - (y_range * 0.1)
        y_max <- y_max + (y_range * 0.1)
        
        return(c(y_min, y_max))
    }
    return(NULL)
}
for (comp_name in names(EXPERIMENT_CONFIG$COMPARISONS)) {
    message("\nProcessing comparison: ", comp_name)
    # Get samples for this comparison
    comp_samples <- sample_table[eval(
        EXPERIMENT_CONFIG$COMPARISONS[[comp_name]], 
        envir = sample_table
    ), ]
    if (nrow(comp_samples) == 0) {
        warning("No samples found for comparison: ", comp_name)
        next
    }
    message("Found ", nrow(comp_samples), " samples for comparison")
    # Create tracks list (rest of track creation code remains the same)
    tracks <- list(GenomeAxisTrack(name = sprintf("Chr %s Axis", chromosome)))
    # Find control sample first
    control_sample <- sample_table[
        sample_table$antibody == "Input" &
        sample_table$rescue_allele == comp_samples$rescue_allele[1],
    ][1,]
    
    if (!is.null(control_sample)) {
        control_bigwig <- list.files(
            file.path(base_dir, "coverage"),
            pattern = paste0(control_sample$experiment_number, ".*normalized.*\\.bw$"),
            full.names = TRUE
        )[1]
        
        # Try primary control
        valid_control <- validate_bigwig(control_bigwig, control_sample$experiment_number)
        
        # If primary control fails, try fallback
        if (is.null(valid_control)) {
            message("Primary control not found, searching for fallback")
            valid_control <- find_fallback_control(sample_table, file.path(base_dir, "coverage"), pattern)
        }
        
        if (!is.null(valid_control)) {
            tracks[[length(tracks) + 1]] <- DataTrack(
                import(valid_control),
                type = "l",
                name = "Input Control",
                col = "#808080"
            )
        }
    }

    label_categories <- EXPERIMENT_CONFIG$COLUMN_ORDER  # Using the ordered columns from config
    
    # Generate labels for all samples in comparison
    comparison_labels <- unique_labeling(
        comp_samples,
        label_categories
    )

    #
    test_labels <- sample_generate_labels(
        comp_samples,
        EXPERIMENT_CONFIG$COLUMN_ORDER,
        verbose = TRUE
    )
    message("Generated labels for tracks:")
    print(test_labels)
    
    # Process each sample
    for (i in seq_len(nrow(comp_samples))) {
        bigwig_file <- list.files(
            file.path(base_dir, "coverage"),
            pattern = paste0(comp_samples$experiment_number[i], ".*normalized.*\\.bw$"),
            full.names = TRUE
        )[1]
        # Validate sample bigwig
        valid_bigwig <- validate_bigwig(bigwig_file, comp_samples$experiment_number[i])
        if (!is.null(valid_bigwig)) {
            message("Adding sample track from: ", basename(valid_bigwig))
            tracks[[length(tracks) + 1]] <- DataTrack(
                import(valid_bigwig),
                type = "l",
                name = comparison_labels[i],
                col = "#fd0036",
                chromosome = chromosome_roman
            )
        }
    }

    # Add feature track if it exists
    if (!is.null(feature_track)) {
        tracks[[length(tracks) + 1]] <- feature_track
    }

    output_file <- file.path(
        base_dir,
        "plots",
        sprintf(
            "%s_%s_%s_chr%s_all_comparisons.svg",
            TIMESTAMP,
            experiment_id,
            comp_name,
            chromosome
        )
    )
    message("Generating plot: ", basename(output_file))
    # Create plot
    tryCatch({
        y_limits <- calculate_track_limits(tracks)
        if (!is.null(y_limits)) {
            message(sprintf("Setting y-limits to: [%.2f, %.2f]", y_limits[1], y_limits[2]))
            svg(output_file)
            plotTracks(
                tracks,
                main = sprintf("Chr %s - %s", chromosome, sub("comp_", "", comp_name)),
                chromosome = chromosome_roman,
                ylim = y_limits
            )
            dev.off()
            message("Plot created successfully")
        } else {
            svg(output_file)
            plotTracks(
                tracks,
                main = sprintf("Chr %s - %s", chromosome, sub("comp_", "", comp_name)),
                chromosome = chromosome_roman
            )
            dev.off()
            message("Plot created successfully")
        }
    }, error = function(e) {
        warning("Failed to create plot: ", e$message)
    })
}

message("All plots generated.")
# rsync -avz --progress luised94@luria.mit.edu:/home/luised94/data/241007Bel/plots/ $HOME/data/241007Bel/plots/
