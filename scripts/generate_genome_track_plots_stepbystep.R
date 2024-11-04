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
message("Loaded config with ", length(EXPERIMENT_CONFIG$COMPARISONS), " comparisons")

# 4. Load sample table
metadata_file <- list.files(
    file.path(base_dir, "documentation"),
    pattern = "sample_table",
    full.names = TRUE
)[1]
message("Found metadata file: ", metadata_file)
sample_table <- read.delim(metadata_file, sep = "\t")
message("Loaded sample table with ", nrow(sample_table), " rows")

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
feature_ranges <- rtracklayer::import.bed(feature_file)
feature_track <- AnnotationTrack(
    feature_ranges,
    name = "Origin Peaks (Eaton 2010)"
)
message("Created feature track")

# 8. Process a single comparison (for testing)
# Pick first comparison
comp_name <- names(EXPERIMENT_CONFIG$COMPARISONS)[1]
message("Testing with comparison: ", comp_name)

# Get samples for this comparison
comp_samples <- sample_table[eval(
    EXPERIMENT_CONFIG$COMPARISONS[[comp_name]], 
    envir = sample_table
), ]
message("Found ", nrow(comp_samples), " samples for comparison")

# 9. Create test plot for first comparison
# Setup chromosome
chromosome <- 10
chromosome_roman <- paste0("chr", as.roman(chromosome))
message("Processing chromosome: ", chromosome_roman)

# Create basic tracks
tracks <- list(
    GenomeAxisTrack(
        name = sprintf("Chr %s Axis", chromosome)
    )
)

# Add control track if available
control_sample <- sample_table[
    sample_table$antibody == "Input" &
    sample_table$rescue_allele == comp_samples$rescue_allele[1],
][1,]

if (!is.null(control_sample)) {
    message("Found control sample: ", control_sample$sample_ID)
    control_bigwig <- list.files(
        file.path(base_dir, "coverage"),
        pattern = paste0(control_sample$sample_ID, ".*normalized.*\\.bw$"),
        full.names = TRUE
    )[1]
    if (!is.na(control_bigwig)) {
        message("Found control bigwig: ", control_bigwig)
        tracks[[length(tracks) + 1]] <- DataTrack(
            import(control_bigwig),
            type = "l",
            name = "Input Control",
            col = "#808080"
        )
    }
}

# 10. Add sample tracks
for (i in seq_len(nrow(comp_samples))) {
    bigwig_file <- list.files(
        file.path(base_dir, "coverage"),
        pattern = paste0(comp_samples$sample_ID[i], ".*normalized.*\\.bw$"),
        full.names = TRUE
    )[1]
    
    if (!is.na(bigwig_file)) {
        message("Processing sample ", i, ": ", bigwig_file)
        tracks[[length(tracks) + 1]] <- DataTrack(
            import(bigwig_file),
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

# 11. Add feature track and generate test plot
tracks[[length(tracks) + 1]] <- feature_track
message("Generated ", length(tracks), " tracks")

# 12. Create test plot
output_file <- file.path(
    base_dir, 
    "plots", 
    sprintf("test_%s_%s.svg", comp_name, chromosome_roman)
)
message("Creating test plot: ", output_file)

svg(output_file)
plotTracks(
    tracks,
    main = sub("comp_", "", comp_name),
    chromosome = chromosome_roman
)
dev.off()
message("Plot created successfully")
