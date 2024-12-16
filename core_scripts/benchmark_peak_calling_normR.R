################################################################################
# Benchmarking peak calling sample processing using control sample
################################################################################
# PURPOSE: Perform peak calling on ChIP-seq data using normR
# USAGE: source("peak_calling_normr.R")
# !! ----> REQUIRED UPDATES: experiment_id, input files paths
# STRUCTURE: Configuration -> Validation -> Processing -> Output
# VALIDATION: Input file existence, format compatibility, memory requirements
# DEPENDENCIES: normR, GenomicRanges, rtracklayer
# COMMON ISSUES: Memory limitations for large BAM files, sparse control samples
# OUTPUT: BED files with called peaks, quality metrics
# AUTHOR: [Your Name]
# DATE: [Current Date]
# VERSION: 1.0.0
################################################################################
################################################################################
# Configuration and Debug Settings
################################################################################
DEBUG_CONFIG <- list(
    single_file_mode = TRUE,
    verbose = TRUE,
    dry_run = TRUE,
    files_to_process_idx = 7
)

NORMR_CONFIG <- list(
    # Core analysis parameters
    bin_size = 1000L,    # Base pairs, adjusted for yeast genome
    min_mapq = 30L,      # Minimum mapping quality for high-quality alignments
    iterations = 10L,    # Number of EM iterations
    processors = 1L,     # Number of processors to use
    paired_end = FALSE,  # Specify single-end data
    
    # Peak calling thresholds
    default_fdr = 10^-6,                    # Default threshold for significant peaks
    min_enrichment_score = 1.5,            # Minimum fold change
    min_reads_per_bin = 1,                 # Minimum reads required per bin
    
    # S. cerevisiae ORC-specific parameters
    expected_peak_range = list(
        min = 100,
        max = 500,
        typical = 250:400
    ),
    
    # Genome specifications
    expected_chromosomes = 16,
    chromosome_prefix = "chr",
    chromosome_adjustment = TRUE,  # Round chromosome sizes to bin size multiples
    
    # File patterns and templates
    bam_pattern = "consolidated_([0-9]{5,6})_sequence_to_S288C_sorted\\.bam$",
    genome_pattern = "S288C_refgenome.fna",
    
    # Output templates
    output_name_template = "%s_peaks_%s_vs_%s_%s.bed",        # timestamp, chip, input, package
    region_file_template = "%s_regions_%s_vs_%s_%s.tsv",
    bedgraph_template = "%s_enrichment_%s_vs_%s_%s.bedGraph",
    
    # Output specifications
    region_columns = c(
        "chromosome", "start", "end",
        "treatment_count", "control_count",
        "enrichment", "qvalue", "peak_class"
    )
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

################################################################################
# Load Required Libraries
################################################################################
required_packages <- c("normr", "GenomicRanges", "rtracklayer", "cosmo", "ggplot2", "Gviz", "stats", "boot")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
}

################################################################################
# Directory Setup
################################################################################
experiment_id <- "241010Bel"  # !! UPDATE THIS
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
bam_dir <- file.path(base_dir, "alignment")
peak_dir <- file.path(base_dir, "peak")
metadata_path <- file.path(base_dir, "documentation",
                          paste0(experiment_id, "_sample_grid.csv"))

config_path <- file.path(base_dir, "documentation", paste0(experiment_id, "_bmc_config.R"))
# Validate directories
stopifnot(
    "Base directory does not exist" = dir.exists(base_dir),
    "BAM directory does not exist" = dir.exists(bam_dir),
    "Metadata path file does not exist" = file.exists(metadata_path)
)

# Create output directory if it doesn't exist
if (!dir.exists(peak_dir)) {
    dir.create(peak_dir, recursive = TRUE, showWarnings = FALSE)
}

################################################################################
# Load and Validate Experiment Configuration and Dependencies
################################################################################
# Bootstrap phase
bootstrap_path <- normalizePath("~/lab_utils/core_scripts/functions_for_file_operations.R",
                              mustWork = FALSE)
if (!file.exists(bootstrap_path)) {
    stop(sprintf("[FATAL] Bootstrap file not found: %s", bootstrap_path))
}
source(bootstrap_path)

# Define required dependencies
required_modules <- list(
    list(
        path = "~/lab_utils/core_scripts/functions_for_logging.R",
        description = "Logging functions",
        required = TRUE
    ),
    list(
        path = "~/lab_utils/core_scripts/functions_for_metadata_processing.R",
        description = "Process metadata grid for downstream analysis.",
        required = TRUE
    ),
    list(
        path = "~/lab_utils/core_scripts/functions_for_peak_calling.R",
        description = "Load functions required for peak calling and benchmarking script",
        required = TRUE
    ),
    list(
        path = config_path,
        description = "Load EXPERIMENT CONFIG object.",
        required = TRUE
    )
)

# Validate module structure
stopifnot(
    "modules must have required fields" = all(sapply(required_modules, function(m) {
        all(c("path", "description", "required") %in% names(m))
    }))
)

# Load dependencies with status tracking
load_status <- lapply(required_modules, function(module) {
    if (DEBUG_CONFIG$verbose) {
        cat(sprintf("\n[LOADING] %s\n", module$description))
    }

    success <- safe_source(module$path, verbose = TRUE)

    if (!success && module$required) {
        stop(sprintf(
            "[FATAL] Failed to load required module: %s\n  Path: %s",
            module$description, module$path
        ))
    } else if (!success) {
        warning(sprintf(
            "[WARNING] Optional module not loaded: %s\n  Path: %s",
            module$description, module$path
        ))
    }

    return(list(
        module = module$description,
        path = module$path,
        loaded = success
    ))
})

# Display loading summary using ASCII
if (DEBUG_CONFIG$verbose) {
    cat("\n=== Module Loading Summary ===\n")
    invisible(lapply(load_status, function(status) {
        cat(sprintf(
            "[%s] %s\n    Path: %s\n",
            if(status$loaded) "+" else "-",
            status$module,
            status$path
        ))
    }))
}

if (DEBUG_CONFIG$verbose) {
    print_config_settings(DEBUG_CONFIG)
}

################################################################################
# File Discovery and Validation
################################################################################
bam_files <- list.files(
    path = bam_dir,
    pattern = NORMR_CONFIG$bam_pattern,
    recursive = TRUE,
    full.names = TRUE
)

# Extract sample IDs from fastq filenames
sample_ids <- gsub(
    pattern = NORMR_CONFIG$bam_pattern,
    replacement = "\\1",
    x = basename(bam_files)
)

# Find bigwig files
bigwig_pattern <- sprintf("_%s\\.bw$", EXPERIMENT_CONFIG$NORMALIZATION$active)
bigwig_files <- list.files(
    file.path(base_dir, "coverage"),
    pattern = bigwig_pattern,
    full.names = TRUE
)
normalization_method <- sub(".*_([^_]+)\\.bw$", "\\1",
                          basename(bigwig_files[1]))

stopifnot(
    "No BAM files found" = length(bam_files) > 0,
    "No bigwig files found" = length(bigwig_files) > 0,
    "No sample_ids found" = length(sample_ids) > 0,
    "Number of sample ids must match number of files" = length(bam_files) == length(sample_ids)
)

################################################################################
# Load and process sample metadata
################################################################################
# Load and process metadata
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

# Enforce factor levels from config
for (col_name in names(EXPERIMENT_CONFIG$CATEGORIES)) {
    if (col_name %in% colnames(metadata)) {
        metadata[[col_name]] <- factor(
            metadata[[col_name]],
            levels = EXPERIMENT_CONFIG$CATEGORIES[[col_name]],
            ordered = TRUE
        )
    }
}

# Sort metadata using config column order
sorted_metadata <- metadata[do.call(
    order,
    metadata[EXPERIMENT_CONFIG$COLUMN_ORDER]
), ]

# Add sample IDs to metadata
sorted_metadata$sample_id <- sample_ids


# Add after processing metadata but before track creation
short_sample_ids <- create_minimal_identifiers(
    sorted_metadata$sample_id,
    verbose = DEBUG_CONFIG$verbose
)

# Create mapping between full and short IDs
sample_id_mapping <- setNames(short_sample_ids, sorted_metadata$sample_id)

################################################################################
# Load reference genome
################################################################################
ref_genome_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    pattern = NORMR_CONFIG$genome_pattern,
    full.names = TRUE,
    recursive = TRUE
)[1]

stopifnot(
    "No reference genome found. One expected." = length(ref_genome_file) == 1
)

genome_data <- Biostrings::readDNAStringSet(ref_genome_file)

# Convert genome data to required format
genome_info <- data.frame(
    chrom = names(genome_data),
    size = genome_data@ranges@width
)

# Validation
stopifnot(
    "genome_info must be a data frame" = is.data.frame(genome_info),
    "genome_info must have exactly 16 rows and 2 columns" = all(dim(genome_info) == c(16, 2)),
    "Column names must be 'chrom' and 'size'" = all(colnames(genome_info) == c("chrom", "size")),
    "Chromosome sizes must be positive" = all(genome_info$size > 0),
    "Chromosome names must start with 'chr'" = all(grepl("^chr", genome_info$chrom)),
    "No missing values allowed" = !any(is.na(genome_info))
)

################################################################################
# Load feature file (annotation)
################################################################################
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

################################################################################
# Process given control file.
################################################################################
# Set the appropriate control for the given experiment.
control_idx <- 7
chromosome_to_plot <- 10
chromosome_width <- genome_data[chromosome_to_plot]@ranges@width
chromosome_roman <- paste0("chr", utils::as.roman(chromosome_to_plot))
#TODO: Add interactive that displays the chosen control to confirm benchmarking.
# Find control sample
control_sample <- find_control_sample(
    experimental_sample = sorted_metadata[control_idx, ],
    metadata = sorted_metadata,
    control_factors = EXPERIMENT_CONFIG$CONTROL_FACTORS
)
# Get sample IDs
chip_id <- sorted_metadata[control_idx, "sample_id"]
input_id <- control_sample$sample_id
track_name <- sprintf(
    "%s: %s - %s",
    sample_id_mapping[chip_id],
    sorted_metadata$short_name[control_idx],
    current_samples$antibody[control_idx]
)
if (DEBUG_CONFIG$verbose) {
    message(sprintf("Processing sample: %s", chip_id))
    message(sprintf("Using control sample: %s", input_id))
    message(sprintf("Track name for visualization: %s", track_name))
}
# Find BAM files using grepl for more robust matching
chip_bam <- bam_files[grepl(paste0("consolidated_", chip_id, "_sequence_to_S288C_sorted\\.bam$"), bam_files)]
input_bam <- bam_files[grepl(paste0("consolidated_", input_id, "_sequence_to_S288C_sorted\\.bam$"), bam_files)]
bigwig_file <- bigwig_files[grepl(chip_id, bigwig_files)][1]
input_bigwig_file <- bigwig_files[grepl(input_bam, bigwig_files)][1]

# Validate BAM files
stopifnot(
    "ChIP BAM file not found" = length(chip_bam) == 1,
    "Input BAM file not found" = length(input_bam) == 1,
    "ChIP BAM file does not exist" = file.exists(chip_bam),
    "Input BAM file does not exist" = file.exists(input_bam),
    "ChIP bigwig file does not exist" = file.exists(bigwig_file)
    "Input bigwig file does not exist" = file.exists(input_bigwig_file)
)
# Generate output filename using short IDs
output_filename <- sprintf(
    NORMR_CONFIG$output_name_template,
    TIMESTAMPS$full,
    sample_id_mapping[chip_id],
    sample_id_mapping[input_id],
    "normr"
)
output_path <- file.path(peak_dir, output_filename)
if (DEBUG_CONFIG$verbose) {
    message(sprintf("Output will be written to: %s", output_path))
}
# Perform peak calling
tryCatch({
    if (DEBUG_CONFIG$verbose) {
        message("\nStarting peak calling...")
        message(sprintf("Genome size: %d bp across %d chromosomes",
                    sum(genome_info$size),
                    nrow(genome_info)))
    }
    if (DEBUG_CONFIG$verbose) {
        message("\nPre-processing count data...")
    }
    # Create count configuration for single-end data
    count_config <- normr::countConfigSingleEnd(
        binsize = NORMR_CONFIG$bin_size,
        mapq = NORMR_CONFIG$min_mapq,
        filteredFlag = 1024,  # Filter duplicates
        shift = 0  # No shift in 3' direction
    )
    if (DEBUG_CONFIG$verbose) {
        message("\nStarting enrichment analysis...")
        message(sprintf("Processing ChIP: %s", basename(chip_bam)))
        message(sprintf("Using Input: %s", basename(input_bam)))
    }
    # Run normR enrichR with count configuration
    enrichment_results <- normr::enrichR(
        treatment = chip_bam,
        control = input_bam,
        genome = genome_info,
        countConfig = count_config,
        iterations = NORMR_CONFIG$iterations,
        procs = NORMR_CONFIG$processors,
        verbose = DEBUG_CONFIG$verbose
    )

    ranges <- normr::getRanges(enrichment_results)
    counts <- normr::getCounts(enrichment_results)
    qvals <- normr::getQvalues(enrichment_results)
    stopifnot("Q-values must be numeric and between 0 and 1" = 
              is.numeric(qvals) && all(qvals >= 0 & qvals <= 1, na.rm = TRUE))
    enrichment_scores <- normr::getEnrichment(enrichment_results)
    peak_classes <- normr::getClasses(enrichment_results)

    # Simple 500 peaks
    top_500_peaks <- ranges[order(qvals)[1:500]]
    
    # Elbow method (using the geometric approach we discussed)
    elbow_peaks <- ranges[qvals <= find_elbow_point(qvals)$elbow_qval]
    
    # Benjamini-Hochberg FDR
    bh_threshold <- max(stats::p.adjust(qvals, method="BH") <= 0.05)
    bh_peaks <- ranges[qvals <= bh_threshold]
    
    # Bootstrapping (using 95% CI upper bound for conservativeness)
    boot_results <- boot::boot(qvals, bootstrap_elbow, R=1000)
    ci <- boot::boot.ci(boot_results, type="perc")
    boot_peaks <- ranges[qvals <= ci$percent[5]]

    # --- Output and Visualization ---
    # Print peak counts
    print(paste("Number of top 500 peaks:", length(top_500_peaks)))
    print(paste("Number of elbow peaks:", length(elbow_peaks)))
    print(paste("Number of BH peaks:", length(bh_peaks)))
    print(paste("Number of bootstrap peaks:", length(boot_peaks)))
    
    # Print Elbow Method Results
    print("Elbow Method Results:")
    print(paste("Elbow q-value:", elbow_results$elbow_qval))
    print(paste("Elbow Peaks:", length(elbow_peaks)))
    
    # Print BH Threshold
    print(paste("Benjamini-Hochberg threshold at 5% FDR:", bh_threshold))
    
    # Print Bootstrap Results
    print(paste("95% CI for elbow q-value:", round(ci$percent[4], 6), "-", round(ci$percent[5], 6)))
    
    # Stability Check
    stability_results <- stability_check(qvals)
    print("Stability check results:")
    print(stability_results)
    
    # Genome Track Plot
    # Use tryCatch to handle potential plotting errors gracefully
    tryCatch({

        track_data <- rtracklayer::import(bigwig_file)
      plotTracks(list(
        Gviz::GenomeAxisTrack(),
        Gviz::DataTrack(track_data, name = track_name, type = "l")
        Gviz::AnnotationTrack(top_500_peaks, name = "Top 500"),
        Gviz::AnnotationTrack(elbow_peaks, name = "Elbow Method"),
        Gviz::AnnotationTrack(bh_peaks, name = "BH FDR"),
        Gviz::AnnotationTrack(boot_peaks, name = "Bootstrap"),
        Gviz::AnnotationTrack(features, name = "Eaton 2010")
      ), chromosome = chromosome_roman)
    }, error = function(e) {
      print(paste("Error plotting tracks:", e$message))
    })
    
    # Overlap Calculation and Output
    overlap_500 <- calculate_overlap(top_500_peaks, features)
    overlap_elbow <- calculate_overlap(elbow_peaks, features)
    overlap_bh <- calculate_overlap(bh_peaks, features)
    overlap_boot <- calculate_overlap(boot_peaks, features)
    
    print(paste("Top 500 overlap:", round(overlap_500$percent_overlap,2), "%"))
    print(paste("Elbow method overlap:", round(overlap_elbow$percent_overlap,2), "%"))
    print(paste("BH FDR overlap:", round(overlap_bh$percent_overlap,2), "%"))
    print(paste("Bootstrap overlap:", round(overlap_boot$percent_overlap,2), "%"))
    
    # Output unique peaks for further analysis if needed
    #Example: write unique peaks from elbow method to bed file
    #if(length(overlap_elbow$unique_peaks)>0){
    #    rtracklayer::export(overlap_elbow$unique_peaks, "unique_elbow_peaks.bed", format = "BED")
    #}

    # Export results if not in dry run mode
    if (!DEBUG_CONFIG$dry_run) {
        normr::exportR(
            obj = enrichment_results,
            filename = output_path,
            type = "bed",
            fdr = NORMR_CONFIG$fdr_threshold
        )
    }
}, error = function(e) {
    message(sprintf("Error processing sample %s: %s",
                    chip_id,
                    e$message))
})
