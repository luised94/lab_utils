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
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
message("Working directory: ", base_dir)

# 3. Load and verify config
source("~/lab_utils/scripts/bmc_config.R")
source("~/lab_utils/scripts/extract_bmcIDmetadata_process.R")
source("~/lab_utils/scripts/genome_core.R")
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
sample_table$experiment_number <- sort(sample_table$experiment_number)
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
validate_bigwig <- function(bigwig_path, sample_id) {
    if (is.na(bigwig_path) || !file.exists(bigwig_path)) {
        warning(sprintf("Bigwig file not found for sample %s", sample_id))
        return(NULL)
    }
    tryCatch({
        rtracklayer::import(bigwig_path)
        return(bigwig_path)
    }, error = function(e) {
        warning(sprintf("Invalid bigwig file for sample %s: %s", sample_id, e$message))
        return(NULL)
    })
}

# Setup chromosome
chromosome <- 10
chromosome_roman <- paste0("chr", as.roman(chromosome))
message("Processing chromosome: ", chromosome_roman)
# 8. Process all comparisons
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
        # Validate control bigwig
        valid_control <- validate_bigwig(control_bigwig, control_sample$experiment_number)
        if (!is.null(valid_control)) {
            message("Adding control track from: ", basename(valid_control))
            tracks[[length(tracks) + 1]] <- DataTrack(
                import(valid_control),
                type = "l",
                name = "Input Control",
                col = "#808080"
            )
        }
    }
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
                name = paste(
                    comp_samples$antibody[i],
                    comp_samples$rescue_allele[i],
                    comp_samples$time_after_release[i],
                    sep = "_"
                ),
                col = "#fd0036"
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
            "%s_%s_%s_chr%s_test_for_all_comparisons.svg",
            TIMESTAMP,
            experiment_id,
            comp_name,
            chromosome
        )
    )
    message("Generating plot: ", basename(output_file))
    # Create plot
    tryCatch({
        svg(output_file)
        plotTracks(
            tracks,
            main = sprintf("Chr %s - %s", chromosome, sub("comp_", "", comp_name)),
            chromosome = chromosome_roman
        )
        dev.off()
        message("Plot created successfully")
    }, error = function(e) {
        warning("Failed to create plot: ", e$message)
    })
}
message("All plots generated.")
