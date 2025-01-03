################################################################################
# Peak Calling using normR for ChIP-seq Analysis
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
    files_to_process_idx = 1
)

NORMR_CONFIG <- list(
    # Core analysis parameters
    bin_size = 1000L,    # Base pairs, adjusted for yeast genome
    min_mapq = 30L,      # Minimum mapping quality for high-quality alignments
    iterations = 10L,    # Number of EM iterations
    processors = 1L,     # Number of processors to use
    paired_end = FALSE,  # Specify single-end data

    # Peak calling thresholds
    fdr_thresholds = c(0.01, 0.05, 0.1),  # For classification of peak strength
    default_fdr = 0.05,                    # Default threshold for significant peaks
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
required_packages <- c("normr", "GenomicRanges", "rtracklayer")
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

stopifnot(
    "No BAM files found" = length(bam_files) > 0,
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
# Process Files
################################################################################
files_to_process <- if (DEBUG_CONFIG$single_file_mode) {
    DEBUG_CONFIG$files_to_process_idx
} else {
    seq_along(bam_files)
}
stop("this is the strategic breakpoint.")
for (file in files_to_process) {
    # Find control sample
    control_sample <- find_control_sample(
        experimental_sample = sorted_metadata[file, ],
        metadata = sorted_metadata,
        control_factors = EXPERIMENT_CONFIG$CONTROL_FACTORS
    )

    # Get sample IDs
    chip_id <- sorted_metadata[file, "sample_id"]
    input_id <- control_sample$sample_id

    if (DEBUG_CONFIG$verbose) {
        message(sprintf("Processing sample: %s", chip_id))
        message(sprintf("Using control sample: %s", input_id))
    }

    # Find BAM files using grepl for more robust matching
    chip_bam <- bam_files[grepl(paste0("consolidated_", chip_id, "_sequence_to_S288C_sorted\\.bam$"), bam_files)]
    input_bam <- bam_files[grepl(paste0("consolidated_", input_id, "_sequence_to_S288C_sorted\\.bam$"), bam_files)]

    # Validate BAM files
    stopifnot(
        "ChIP BAM file not found" = length(chip_bam) == 1,
        "Input BAM file not found" = length(input_bam) == 1,
        "ChIP BAM file does not exist" = file.exists(chip_bam),
        "Input BAM file does not exist" = file.exists(input_bam)
    )

    # Check if using same file as treatment and control
    is_self_control <- identical(chip_bam, input_bam)
    if (is_self_control && DEBUG_CONFIG$verbose) {
        message("Note: Using same file as treatment and control - expect minimal/no enrichment")
    }

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
        stopifnot(class(ranges) == "GRanges")
        counts <- normr::getCounts(enrichment_results)
        stopifnot(is.list(counts))
        qvals <- normr::getQvalues(enrichment_results)
        stopifnot(is.numeric(qvals) | is.null(qvals))
        enrichment_scores <- normr::getEnrichment(enrichment_results)
        stopifnot(is.numeric(enrichment_scores))
        peak_classes <- normr::getClasses(enrichment_results)
        stopifnot(is.numeric(peak_classes) | is.null(peak_classes))

        # --- 4. Filter significant peaks ---
        significant_peaks <- ranges[!is.na(qvals) & qvals <= NORMR_CONFIG$default_fdr]
        
        if (DEBUG_CONFIG$verbose) {
            message(sprintf("Number of significant peaks: %d", length(significant_peaks)))
            message(sprintf("FDR threshold used: %g", NORMR_CONFIG$fdr_threshold))
            if(length(significant_peaks) > 0) {
                first_peak_index <- which(GenomicRanges::ranges(ranges) == GenomicRanges::ranges(significant_peaks[1]))
                message(sprintf("First significant peak - Treatment count: %d, Control count: %d", 
                                counts$treatment[first_peak_index], counts$control[first_peak_index]))
            }
        }
        
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
}

if (DEBUG_CONFIG$verbose) {
    print_config_settings(DEBUG_CONFIG)
}
