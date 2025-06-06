################################################################################
# BMC Experiment Configuration
################################################################################
# PURPOSE:
#   Defines and validates experimental design for BMC ChIP-seq experiments,
#   including sample categories, valid combinations, and comparison groups.
# USAGE:
#   1. Use '/!!' in vim/neovim to jump to required updates
#   2. Modify METADATA section with experiment details
#   3. Update CATEGORIES if experimental design changes
#   4. Review INVALID_COMBINATIONS and EXPERIMENTAL_CONDITIONS
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
# Comparison Naming Convention:
# comp_[Antibody]_[Variables]_vs_[Baseline]_[Context]_[Modifier]
# - [Antibody]: Official antibody code (e.g., HM1108, V5, ALFA)
# - [Variables]: Main experimental factors being compared
# - [vs_Baseline]: Reference state for comparison (if applicable)
# - [Context]: Additional conditions (e.g., AuxinTreated, TimeCourse)
# - [Modifier]: Special analysis flags or time points (e.g., 2h)
# Examples:
# comp_HM1108_NoRescue_vs_WT
# comp_V5_Auxin_2h
# comp_ALFA_AllAlleles_Auxin_2h
# Note: Terms like "None", "For", "With", and "And" are avoided for clarity.
# Time is abbreviated as "0h", "2h", etc. "NoRescue" is used instead of "None".
#
# VALIDATION:
#   1. Category Values: Must be character vectors, unique
#   2. Column References: All referenced columns must exist
#   3. Column Order: Must include all category columns
#   4. Sample Count: Must match EXPECTED_SAMPLES
# DEPENDENCIES: R base packages only, functions_for_bmc_config_validation.R for validation functions
# COMMON ISSUES:
#   1. Mismatched EXPERIMENT_ID -> Check format YYMMDD'Bel'
#   2. Wrong sample count -> Review INVALID_COMBINATIONS
#   3. Missing categories -> Check CATEGORIES vs COLUMN_ORDER
# AUTHOR: Luis
# DATE: 2024-11-27
# VERSION: 2.1.0
################################################################################
# !! Update EXPERIMENT_CONFIG if starting a new experiment.
EXPERIMENT_CONFIG <- list(
    METADATA = list(
        EXPERIMENT_ID = "100303Bel",
        EXPECTED_SAMPLES = 5,
        VERSION = "1.0.0",
        PROJECT_ID = "project_001"
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
            (temperature == "23" |                  # no orc1-161 at 23�C
             antibody == "ORC" |                    # no orc1-161 with ORC
             cell_cycle %in% c("async", "alpha"))      # no orc1-161 in async or alpha
        ),
        # Group 2: ORC antibody restrictions
        orc_restrictions = quote(
            antibody == "ORC" &
            (temperature == "23" |                  # no ORC at 23�C
             cell_cycle %in% c("alpha", "async"))      # no ORC in alpha or async
        ),
        # Group 2: ORC antibody restrictions
        orc_restrictions = quote(
            antibody == "ORC" & 
            (temperature == "23" |                  # no ORC at 23�C
             cell_cycle %in% c("alpha", "async"))      # no ORC in alpha or async
        ),
        # Group 3: Nucleosome and temperature restrictions
        nucleosome_temp_restrictions = quote(
            (antibody == "Nucleosomes" & temperature == "23" & cell_cycle %in% c("alpha", "nocodazole")) | # no Nucleosomes at 23�C in alpha/nocodazole
            (temperature == "37" & cell_cycle == "async")    # no async at 37�C
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

#########
# setup_bmc_experiment_flags: Control what is displayed during setup_bmc_experiment.R
# Distinct from debug configuration, this is relevant when setting up the experiment only.
#########
show_all_metadata <- TRUE
show_particular_metadata <- TRUE
category_to_show <- "antibody"
value_to_show <- "HM1108"
values_in_category <- unlist(EXPERIMENT_CONFIG$CATEGORIES[category_to_show])

stopifnot(
    "show_all_metadata has to be logical" = is.logical(show_all_metadata),
    "flag must be logical type" = is.logical(show_particular_metadata),
    "category must be in grid" = category_to_show %in% names(EXPERIMENT_CONFIG$CATEGORIES),
    "Value must be in category" = value_to_show %in% values_in_category
)

#todo: need to add guard statement around this section to prevent double loading of config objects when I run multi-experiment scripts. Also use default config with all values. Would need to adjust the parse arguments script as well. Pretty complex refactor. Could just to override_configuration.
################################################################################
# Time Configurations
################################################################################
TIME_CONFIG <- list(
    # Format specifications
    timestamp_format = "%Y%m%d_%H%M%S",    # YYYYMMDD_HHMMSS
    date_format = "%Y%m%d",                # YYYYMMDD

    # Current values
    current_timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
    current_date = format(Sys.Date(), "%Y%m%d")
)

################################################################################
# DEBUG CONFIGURATIONS
################################################################################
RUNTIME_CONFIG <- list(
    # Execution Mode
    #show_debug_output = TRUE,      # (formerly debug_verbose)
    #require_confirmation = FALSE,   # (formerly debug_interactive)
    #validate_extensively = TRUE,    # (formerly debug_validate)

    ## Processing Scope
    #target_comparison = "comp_1108forNoneAndWT",
    #target_chromosome = 10,
    #target_batch = 10,
    #samples_per_batch = 4,

    ## Output Control
    #plot_display_duration = 2

    # Core control flags
    debug_interactive = FALSE,
    debug_verbose = TRUE,
    debug_validate = TRUE,

    # Processing control
    process_single_file = FALSE,
    process_single_comparison = TRUE,
    process_comparison = "comp_1108forNoneAndWT",
    process_chromosome = 14,
    #process_group = 10,
    process_batch = 10,
    process_samples_per_batch = 4,
    #process_samples_per_group = 4,
    process_file_index = 1,
    # Output control
    output_save_plots = FALSE,
    output_dry_run = TRUE,
    output_display_time = 2
)

################################################################################
# Quality control configuration
################################################################################
FASTQC_CONFIG <- list(
    # Version control
    version_required = "0.11.5",
    version_pattern = "^##FastQC\\s+",
    version_max = 1,

    # Parsing patterns
    parse_header = "^##FastQC",
    parse_module_start = ">>",
    parse_module_end = ">>END_MODULE",
    parse_prefix = "#",

    # File handling
    file_pattern = "consolidated_([0-9]{5,6})_sequence_fastqc_data\\.txt$",
    file_format = "consolidated_XXXXXX_sequence_fastqc_data.txt",
    file_suffix = ".tab",
    file_base = "fastqc_data",

    # Paths and references
    path_qc_dir = "quality_control",
    path_module_ref = file.path(Sys.getenv("HOME"), "data", "fastqc_module_reference.rds"),

    # Module configuration
    module_list = character(0)
)

################################################################################
# VISUALIZATION AND DISPLAY CONFIGURATIONS
################################################################################
VIEWER_CONFIG <- list(
    # Paths
    path_base = file.path(Sys.getenv("HOME"), "data"),

    # Patterns
    pattern_svg = "\\.svg$",
    pattern_timestamp = "^[0-9]{8}_[0-9]{6}",  # YYYYMMDD_HHMMSS
    pattern_experiment = "^[0-9]{6}Bel",

    # Display dimensions
    display_width = 10,
    display_height = 8
)

# TODO: Requires overhaul. I like having them inside the list but //
# TODO: the prefix strategy may be too verbose and annoying. //
GENOME_TRACK_CONFIG <- list(
    use_custom_visualization = FALSE,  # Control flag

    # Display dimensions
    display_width = 10,
    display_height = 8,

    # Track Creation
    track_points_default = 1000,
    #track_show_title = TRUE,

    # Track defaults
    track_ylim = c(0, 1000),  # Default y-limits, adjust as needed
    track_sampling_rate = 100,  # Points per base pair for empty tracks

    # Track colors
    color_placeholder = "#cccccc",
    color_input = "#808080",

    # Track naming
    format_sample_track_name = "%s: %s",
    format_control_track_name = "%s: %s - %s",
    format_placeholder_track_name = "%s: %s - %s",
    format_suffix = "(No data)",
    format_genome_axis_track_name = "Chr %s Axis",

    # Labels
    label_always_show = "antibody",
    label_never_show = c("sample_id", "full_name", "short_name", "X__cf_genotype"),
    label_separator = "-",

    # File handling
    file_pattern = "consolidated_.*_sequence\\.fastq$",
    file_sample_id = "consolidated_([0-9]{1,6})_sequence\\.fastq",
    file_sample_id_from_bigwig = "processed_([0-9]{1,6})_sequence_to_S288C_(RPKM|CPM|BPM|RPGC)\\.bw",
    #file_sample_id_from_bigwig = "processed_([0-9]{1,6})_sequence_to_S288C_.*_(RPKM|CPM|BPM|RPGC)\\.bw",
    file_genome_pattern = "S288C_refgenome.fna",
    file_genome_directory = file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    file_feature_directory = file.path(Sys.getenv("HOME"), "data", "feature_files"),
    file_feature_pattern = "eaton_peaks",

    # File Names
    filename_format_group_template = "%s_%s_group%02d_chr%s_%s.svg",
    filename_format_comparison_template = "%s_%s_%s_chr%s_%s.svg",
    title_group_template = paste(
        "%s",               # Title
        "Group: %s",   # Comparison ID
        "Chromosome %s (%d samples)", # Chr info
        "%s",               # Additional info
        "Normalization: %s", # Norm method
        sep = "\n"
    ),
    title_comparison_template = paste(
        "%s",               # Title
        "Comparison: %s",   # Comparison ID
        "Chromosome %s (%d samples)", # Chr info
        "%s",               # Additional info
        "Normalization: %s", # Norm method
        sep = "\n"
    ),
    # Development mode title
    #title_dev_mode = "development",  # Enum: "development" | "publication"
    #title_dev_style = 2,    # Bold
    ## Publication mode title
    #title_pub_template = "%s: Chr%s (%s)",
    #title_pub_size = 1,
    #title_pub_style = 2,    # Bold
    ## Title constraints
    #title_max_width = 40,
    #title_max_lines = 5,
    # Interactive mode
    #interactive_prompt = "Options: [Enter] next plot, 's' skip rest, 'q' quit: ",

    track_defaults_sample = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 1.2,
        background.title = "white",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.6,
        fontface = 1,
        title.width = 1.2
    ),

    track_defaults_placeholder = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 0.8,
        background.title = "white",
        background.panel = "#f5f5f5",    # light gray to indicate "empty"
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9,
        alpha = 0.5,
        grid = FALSE
        #ylim = c(0, 1)                  # fixed range for empty tracks
    ),
    track_defaults_control = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 0.8,
        background.title = "white",
        background.panel = "white",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9
        #alpha = 0.8
    ),
    track_defaults_feature = list(
        showaxis = FALSE,
        showtitle = TRUE,
        size = 0.5,
        background.title = "white",
        background.panel = "#8b7355",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9,
        fill = "#8b4513",
        col = "#8b4513"
    ),

    # New Plot Defaults
    plot_defaults = list(
        margin = 15,
        innerMargin = 5,
        spacing = 10,
        extend.left = 0,
        extend.right = 0,
        col.axis = "black",
        cex.axis = 0.8,
        cex.main = 0.8,
        fontface.main = 2,
        background.panel = "transparent"
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
    #conditions = EXPERIMENT_CONFIG$EXPERIMENTAL_CONDITIONS,
    verbose = validation_verbose
)

validate_column_order(
    categories = EXPERIMENT_CONFIG$CATEGORIES,
    column_order = EXPERIMENT_CONFIG$COLUMN_ORDER,
    verbose = validation_verbose
)

cat("\n[VALIDATED] Experiment configuration loaded successfully\n")
