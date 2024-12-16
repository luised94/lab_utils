################################################################################
# BMC Experiment Configuration
################################################################################
#
# PURPOSE:
#   Defines and validates experimental design for BMC ChIP-seq experiments,
#   including sample categories, valid combinations, and comparison groups.
#
# USAGE:
#   1. Use '/!!' in vim/neovim to jump to required updates
#   2. Modify METADATA section with experiment details
#   3. Update CATEGORIES if experimental design changes
#   4. Review INVALID_COMBINATIONS and EXPERIMENTAL_CONDITIONS
#
# !! ----> REQUIRED UPDATES:
# !! EXPERIMENT_CONFIG$METADATA <- list(
# !!     EXPERIMENT_ID = "241010Bel",
# !!     EXPECTED_SAMPLES = 65,
# !!     VERSION = "1.0.0"
# !! )
#
# STRUCTURE:
#   EXPERIMENT_CONFIG/
#   +-- METADATA/
#   |   +-- EXPERIMENT_ID    # Format: YYMMDD'Bel'
#   |   +-- EXPECTED_SAMPLES # Total valid combinations
#   |   +-- VERSION         # Configuration version
#   +-- CATEGORIES/         # Valid values for each factor
#   +-- INVALID_COMBINATIONS/# Excluded experimental combinations
#   +-- EXPERIMENTAL_CONDITIONS/# Valid sample definitions
#   +-- COMPARISONS/        # Analysis groupings
#   +-- CONTROL_FACTORS/    # Control sample definitions
#   +-- COLUMN_ORDER/       # Standard column arrangement
#
# VALIDATION:
#   1. Category Values: Must be character vectors, unique
#   2. Column References: All referenced columns must exist
#   3. Column Order: Must include all category columns
#   4. Sample Count: Must match EXPECTED_SAMPLES
#
# DEPENDENCIES:
#   - R base packages only
#   - functions_for_bmc_config_validation.R for validation functions
#
# COMMON ISSUES:
#   1. Mismatched EXPERIMENT_ID -> Check format YYMMDD'Bel'
#   2. Wrong sample count -> Review INVALID_COMBINATIONS
#   3. Missing categories -> Check CATEGORIES vs COLUMN_ORDER
#
# AUTHOR: Luis
# DATE: 2024-11-27
# VERSION: 2.0.0
#
################################################################################
# !! Update EXPERIMENT_CONFIG if starting a new experiment.
EXPERIMENT_CONFIG <- list(
    METADATA = list(
        EXPERIMENT_ID = "100303Bel",
        EXPECTED_SAMPLES = 5,
        VERSION = "1.0.0"
    ),

    CATEGORIES = list(
        orc_phenotype = c("WT", "orc1-161"),
        cell_cycle = c("alpha", "nocodazole", "async"),
        temperature = c("23", "37"),
        antibody = c("ORC", "Nucleosomes")
    ),

    INVALID_COMBINATIONS = list(
        # Group 1: orc1-161 restrictions
        orc1_161_restrictions = quote(
            orc_phenotype == "orc1-161" & 
            (temperature == "23" |                  # no orc1-161 at 23øC
             antibody == "ORC" |                    # no orc1-161 with ORC
             cell_cycle %in% c("async", "alpha"))      # no orc1-161 in async or alpha
        ),
        
        # Group 2: ORC antibody restrictions
        orc_restrictions = quote(
            antibody == "ORC" & 
            (temperature == "23" |                  # no ORC at 23øC
             cell_cycle %in% c("alpha", "async"))      # no ORC in alpha or async
        ),
        
        # Group 3: Nucleosome and temperature restrictions
        nucleosome_temp_restrictions = quote(
            (antibody == "Nucleosomes" & temperature == "23" & cell_cycle %in% c("alpha", "nocodazole")) | # no Nucleosomes at 23øC in alpha/nocodazole
            (temperature == "37" & cell_cycle == "async")    # no async at 37øC
        )
    ),

    EXPERIMENTAL_CONDITIONS = list(
        is_orc = quote(orc_phenotype == "WT" & cell_cycle == "nocodazole" & antibody == "ORC"),
        is_nucleosomes = quote(antibody == "Nucleosomes"),
        is_wt_37c = quote(orc_phenotype == "WT" & temperature == "37"),
        is_async = quote(orc_phenotype == "WT" & cell_cycle == "async" & temperature == "23")
    ),

    SAMPLE_CLASSIFICATIONS = list(
        is_positive = quote(orc_phenotype == "WT" & antibody == "ORC"),
        is_negative = quote(orc_phenotype == "orc1-161" & antibody == "ORC" & temperature == "37")
    ),

    COMPARISONS = list(
        wt_vs_mutant = quote(orc_phenotype == "WT" | orc_phenotype == "orc1-161"),
        temp_effect = quote(temperature == "23" | temperature == "37"),
        cell_cycle_effect = quote(cell_cycle %in% c("alpha", "nocodazole", "async"))
    ),

    CONTROL_FACTORS = list(
        genotype = c("orc_phenotype")
    ),

    COLUMN_ORDER = c("antibody", "orc_phenotype", "cell_cycle", "temperature"),

    NORMALIZATION = list(
        methods = c("CPM", "BPM", "RPGC", "RPKM"),
        active = "CPM"
    )
)

################################################################################
# Time Configurations
################################################################################
TIME_CONFIG <- list(
    timestamp_format = "%Y%m%d_%H%M%S",  # YYYYMMDD_HHMMSS
    date_format = "%Y%m%d"               # YYYYMMDD
)
TIME_CONFIG$full <- format(Sys.time(), TIME_CONFIG$timestamp_format)
TIME_CONFIG$date <- format(Sys.Date(), TIME_CONFIG$date_format)

################################################################################
# DEBUG CONFIGURATIONS
################################################################################
DEBUG_CONFIG <- list(
    # Runtime control
    enabled = TRUE,
    interactive = FALSE,
    verbose = TRUE,
    validate_config = TRUE,
    
    # Processing scope
    comparison = "comp_1108forNoneAndWT",
    chromosome = 10,
    group = 10,
    samples_per_group = 4,
    
    # Output control
    save_plots = FALSE,
    dry_run = TRUE,
    display_time = 2
)

DEBUG_CONFIG <- list( # !! UPDATE THIS
    single_file_mode = FALSE,           # Test single file in main logic.
    verbose = TRUE,           # Print processing details
    interactive = TRUE,       # Allow interactive processing
    dry_run = FALSE,         # Skip file writes
    files_to_process_idx = 1  # Process specific files in debug mode
)

################################################################################
# Quality control configuration
################################################################################
FASTQC_CONFIG <- list(
    VERSION = "0.11.5",                    # Expected FastQC version
    VERSION_PATTERN = "^##FastQC\\s+",     # Pattern to match version line
    HEADER_PATTERN = "^##FastQC",          # Pattern to identify FastQC header
    module_separator = ">>",
    module_end = ">>END_MODULE",
    header_prefix = "#",
    fastqc_pattern = "fastqc_data",
    output_suffix = ".tab",
    qc_subdir = "quality_control",
    existing_version_limit = 1,
    module_names = character(0),
    module_reference_file = file.path(
        Sys.getenv("HOME"),
        "data",
        "fastqc_module_reference.rds"
    )
)

FASTQC_CONFIG$FILE_PATTERN <- list(
    REGEX = "consolidated_([0-9]{5,6})_sequence_fastqc_data\\.txt$",
    EXPECTED_FORMAT = "consolidated_XXXXXX_sequence_fastqc_data.txt"  # For error messages
)
################################################################################
# GENOME TRACK CONFIG
################################################################################
PLOT_CONFIG <- list(

    # Basic plot settings
    dimensions = list(
        width = 10,
        height = 8
    ),
    
    # Track configuration
    tracks = list(
        # Visual settings for all tracks
        display = list(
            width = 0.9,
            fontface = 1,
            cex = 0.6,
            background = "white",
            fontcolor = "black",
            border_color = "#E0E0E0"
        ),
        
        # Track-specific colors
        colors = list(
            placeholder = "#cccccc",
            input = "#808080"
        ),
        
        # Track name formatting
        names = list(
            format = "%s: %s",
            control_format = "%s: %s - %s",
            placeholder_suffix = "(No data)"
        )
    ),
    
    # Main title configuration
    main_title = list(
        mode = "development",  # or "publication"
        development = list(
            format = paste(
                "%s",
                "Comparison: %s",
                "Chromosome %s (%d samples)",
                "%s",
                "Normalization: %s",
                sep = "\n"
            ),
            cex = 0.7,      # Size specific to main title
            fontface = 2    # Bold for main title
        ),
        publication = list(
            format = "%s: Chr%s (%s)",
            cex = 1,
            fontface = 2
        ),
        format = list(
            max_width = 40,
            max_lines = 5
        )
    )
)

PLOT_CONFIG <- list(
    width = 10,
    height = 8,
    placeholder_color = "#cccccc",
    input_color ="#808080",
    track_name_format = "%s: %s - %s",
    placeholder_suffix = "(No data)",
    title_format = "%s\nChromosome %s (%d samples)\n%s\nNormalization: %s"

)
################################################################################
# Peak Calling Configurations
################################################################################
NORMR_CONFIG <- list(
    # Peak calling parameters
    fdr_thresholds = c(0.01, 0.05, 0.1),
    default_fdr = 0.000001,
    bin_size = 100,  # Base pairs
    min_mapq = 30,   # Minimum mapping quality

    # S. cerevisiae ORC-specific parameters
    expected_peak_range = list(
        min = 100,
        max = 500,
        typical = 250:400
    ),

    # File patterns
    bam_pattern = "consolidated_([0-9]{5,6})_sequence_to_S288C_sorted\\.bam$",
    genome_pattern = "S288C_refgenome.fna",
    # Output formatting
    output_name_template = "%s_peaks_%s_vs_%s_%s.bed",  # timestamp, chip, input, package
    region_file_template = "%s_regions_%s_vs_%s_%s.tsv",
    bedgraph_template = "%s_enrichment_%s_vs_%s_%s.bedGraph",

    # Genome requirements
    expected_chromosomes = 16,
    chromosome_prefix = "chr",

    # Quality thresholds
    min_enrichment_score = 1.5,
    min_read_count = 10,

    # Binning parameters
    min_reads_per_bin = 1,
    paired_end = FALSE,

    # Count configuration
    bin_size = 1000L,    # Reasonable bin size for yeast
    min_mapq = 30L,      # High quality alignments
    iterations = 10L,    # Number of EM iterations
    processors = 1L,     # Number of processors to use

    # Column specifications for region output
    region_columns = c(
        "chromosome", "start", "end",
        "treatment_count", "control_count",
        "enrichment", "qvalue", "peak_class"
    )
)

# Plot Viewer Configuration
#-----------------------------------------------------------------------------
VIEWER_CONFIG <- list(
    base_dir = file.path(Sys.getenv("HOME"), "data"),
    patterns = list(
        svg = "\\.svg$",
        timestamp = "^[0-9]{8}_[0-9]{6}",  # YYYYMMDD_HHMMSS
        experiment = "^[0-9]{6}Bel"
    ),
    device = list(
        width = 10,
        height = 8
    )
)
################################################################################
# Configuration Validation
################################################################################
experiment_id <- EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID
stopifnot(
    "Experiment ID must be a character string" = is.character(experiment_id),
    "Invalid experiment ID format. Expected: YYMMDD'Bel'" = grepl("^\\d{6}Bel$", experiment_id)
)
source("~/lab_utils/core_scripts/functions_for_bmc_config_validation.R")

# !! Update if you want thourough messages during validation.
validation_verbose <- FALSE  # Set to TRUE for detailed validation output

# Validate configuration structure
if (validation_verbose) cat("\nValidating EXPERIMENT_CONFIG structure...\n")

required_sections <- c("METADATA", "CATEGORIES", "INVALID_COMBINATIONS",
                      "EXPERIMENTAL_CONDITIONS", "COMPARISONS",
                      "CONTROL_FACTORS", "COLUMN_ORDER", "NORMALIZATION",
                      "SAMPLE_CLASSIFICATIONS")


missing_sections <- setdiff(required_sections, names(EXPERIMENT_CONFIG))
if (length(missing_sections) > 0) {
    stop(sprintf("Missing required config sections: %s",
                paste(missing_sections, collapse = ", ")))
}

if (validation_verbose) cat("[PASS] All required sections present\n\n")

# Validate configuration structure
stopifnot(
    "Missing required config sections" =
        all(required_sections %in% names(EXPERIMENT_CONFIG))
)

# Validate each section
validate_category_values(
    EXPERIMENT_CONFIG$CATEGORIES,
    verbose = validation_verbose
)

validate_column_references(
    categories = EXPERIMENT_CONFIG$CATEGORIES,
    comparisons = EXPERIMENT_CONFIG$COMPARISONS,
    control_factors = EXPERIMENT_CONFIG$CONTROL_FACTORS,
    conditions = EXPERIMENT_CONFIG$EXPERIMENTAL_CONDITIONS,
    verbose = validation_verbose
)

validate_column_order(
    categories = EXPERIMENT_CONFIG$CATEGORIES,
    column_order = EXPERIMENT_CONFIG$COLUMN_ORDER,
    verbose = validation_verbose
)

cat("\n[VALIDATED] Experiment configuration loaded successfully\n")
