################################################################################
# BMC ChIP-seq Configuration Setup
# Author: Luis | Date: 2024-11-27 | Version: 2.1.0 Experiment: 250715Bel
################################################################################
# PURPOSE: 
#   Write down the categories and quote expressions to obtain metadata for an experiment
#   Assign experiment id according to BMC, set comparisons, and desired order.
#
# USAGE:
# 1. Update experiment_id (format: YYMMDDBel, e.g., "241122Bel")
# 2. Run setup_bmc_experiment.R with the experiment-id as the parameter.
#
# Comparison Naming Convention:
# comp_[Antibody]_[Variables]_vs_[Baseline]_[Context]_[Modifier]
#
# DEPENDENCIES: NONE
#
# OUTPUTS:
#   No outputs produced by script. Use via setup_bmc_experiment.R
# DESCRIPTION:
#   Describe the experiment.
#   - 
################################################################################
# !! Update EXPERIMENT_CONFIG if starting a new experiment.
# 1) Update description and metadata.
# 2) Set the categories of the experiment. This will produce all of the combinations.
# 3) Use quote expressions and logical conditions to write what should not be in the metadata out of all the combinations.
# 4) Update column order to sort the metadata after filtering.
EXPERIMENT_CONFIG <- list(
  METADATA = list(
    # YYMMDDBel
    EXPERIMENT_ID = "250821Bel",
    EXPECTED_SAMPLES = 5,
    VERSION = "1.0.0",
    PROJECT_ID = "project_002"
  ),

  CATEGORIES = list(
    category1 = c("Value1", "Value2"),
    orc_phenotype = c("WT", "orc1-161"),
    cell_cycle = c("alpha", "nocodazole", "async"),
    temperature = c("23", "37"),
    antibody = c("ORC", "Nucleosomes")
  ),

  CONTROL_FACTORS = list(
    genotype = c("orc_phenotype")
  ),

  COLUMN_ORDER = c(
    "antibody", "orc_phenotype",
    "cell_cycle", "temperature"
  ),

  INVALID_COMBINATIONS = list(
    # Group 1: orc1-161 restrictions
    orc1_161_restrictions = quote(
      orc_phenotype == "orc1-161" & 
      (temperature == "23" |          # no orc1-161 at 23øC
       antibody == "ORC" |          # no orc1-161 with ORC
       cell_cycle %in% c("async", "alpha"))    # no orc1-161 in async or alpha
    ),
    # Group 2: ORC antibody restrictions
    orc_restrictions = quote(
      antibody == "ORC" &
      (temperature == "23" |          # no ORC at 23øC
       cell_cycle %in% c("alpha", "async"))    # no ORC in alpha or async
    ),
    # Group 2: ORC antibody restrictions
    orc_restrictions = quote(
      antibody == "ORC" & 
      (temperature == "23" |          # no ORC at 23øC
       cell_cycle %in% c("alpha", "async"))    # no ORC in alpha or async
    ),
    # Group 3: Nucleosome and temperature restrictions
    nucleosome_temp_restrictions = quote(
      (antibody == "Nucleosomes" & temperature == "23" & cell_cycle %in% c("alpha", "nocodazole")) | # no Nucleosomes at 23øC in alpha/nocodazole
      (temperature == "37" & cell_cycle == "async")  # no async at 37øC
    )
  ),

  # Set the options for which bigwig file is used to create genome track plots.
  # They look mostly the same.
  # @TODO Remove this section of code?
  NORMALIZATION = list(
    methods = c("CPM", "BPM", "RPGC", "RPKM"),
    active = "CPM"
  )

  #EXPERIMENTAL_CONDITIONS = list(
  #  is_orc = quote(orc_phenotype == "WT" & cell_cycle == "nocodazole" & antibody == "ORC"),
  #  is_nucleosomes = quote(antibody == "Nucleosomes"),
  #  is_wt_37c = quote(orc_phenotype == "WT" & temperature == "37"),
  #  is_async = quote(orc_phenotype == "WT" & cell_cycle == "async" & temperature == "23")
  #),

  #SAMPLE_CLASSIFICATIONS = list(
  #  is_positive = quote(orc_phenotype == "WT" & antibody == "ORC"),
  #  is_negative = quote(orc_phenotype == "orc1-161" & antibody == "ORC" & temperature == "37")
  #),

  #COMPARISONS = list(
  #  wt_vs_mutant = quote(orc_phenotype == "WT" | orc_phenotype == "orc1-161"),
  #  temp_effect = quote(temperature == "23" | temperature == "37"),
  #  cell_cycle_effect = quote(cell_cycle %in% c("alpha", "nocodazole", "async"))
  #),

)

#########
# setup_bmc_experiment_flags: Control what is displayed during setup_bmc_experiment.R
# Distinct from debug configuration, this is relevant when setting up the experiment only.
#########
#show_all_metadata <- TRUE
#show_particular_metadata <- TRUE
#category_to_show <- "antibody"
#value_to_show <- "HM1108"
#values_in_category <- unlist(EXPERIMENT_CONFIG$CATEGORIES[category_to_show])
#
#stopifnot(
#  "show_all_metadata has to be logical" = is.logical(show_all_metadata),
#  "flag must be logical type" = is.logical(show_particular_metadata),
#  "category must be in grid" = category_to_show %in% names(EXPERIMENT_CONFIG$CATEGORIES),
#  "Value must be in category" = value_to_show %in% values_in_category
#)

#todo: need to add guard statement around this section to prevent double loading of config objects when I run multi-experiment scripts. Also use default config with all values. Would need to adjust the parse arguments script as well. Pretty complex refactor. Could just to override_configuration.
#-------------------------------------------------------------------------------
# Time Configurations
#-------------------------------------------------------------------------------
TIME_CONFIG <- list(
  # Format specifications
  timestamp_format = "%Y%m%d_%H%M%S",  # YYYYMMDD_HHMMSS
  date_format = "%Y%m%d",        # YYYYMMDD

  # Current values
  current_timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
  current_date = format(Sys.Date(), "%Y%m%d")
)

#-------------------------------------------------------------------------------
# DEBUG CONFIGURATIONS
#-------------------------------------------------------------------------------
RUNTIME_CONFIG <- list(
  # Execution Mode
  #show_debug_output = TRUE,    # (formerly debug_verbose)
  #require_confirmation = FALSE,   # (formerly debug_interactive)
  #validate_extensively = TRUE,  # (formerly debug_validate)

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

#-------------------------------------------------------------------------------
# Quality control configuration
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# VISUALIZATION AND DISPLAY CONFIGURATIONS
#-------------------------------------------------------------------------------
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
    "%s",         # Title
    "Group: %s",   # Comparison ID
    "Chromosome %s (%d samples)", # Chr info
    "%s",         # Additional info
    "Normalization: %s", # Norm method
    sep = "\n"
  ),
  title_comparison_template = paste(
    "%s",         # Title
    "Comparison: %s",   # Comparison ID
    "Chromosome %s (%d samples)", # Chr info
    "%s",         # Additional info
    "Normalization: %s", # Norm method
    sep = "\n"
  ),
  # Development mode title
  #title_dev_mode = "development",  # Enum: "development" | "publication"
  #title_dev_style = 2,  # Bold
  ## Publication mode title
  #title_pub_template = "%s: Chr%s (%s)",
  #title_pub_size = 1,
  #title_pub_style = 2,  # Bold
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
    background.panel = "#f5f5f5",  # light gray to indicate "empty"
    fontcolor.title = "black",
    col.border.title = "#e0e0e0",
    cex.title = 0.7,
    fontface = 1,
    title.width = 0.9,
    alpha = 0.5,
    grid = FALSE
    #ylim = c(0, 1)          # fixed range for empty tracks
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

#-------------------------------------------------------------------------------
# Configuration Validation
#-------------------------------------------------------------------------------
experiment_id <- EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID
category_names <- names(EXPERIMENT_CONFIG$CATEGORIES)
stopifnot(
  "Experiment ID must be a character string" =
    is.character(experiment_id),
  "Invalid experiment ID format. Expected: YYMMDD'Bel'" =
    grepl("^\\d{6}Bel$", experiment_id)
)
#source("~/lab_utils/core_scripts/functions_for_configuration_bmc_validation.R")

# !! Update if you want thourough messages during validation.
#validation_verbose <- FALSE  # Set to TRUE for detailed validation output

# Validate configuration structure
#if (validation_verbose) cat("\nValidating EXPERIMENT_CONFIG structure...\n")

required_sections <- c(
  "METADATA", "CATEGORIES", "INVALID_COMBINATIONS",
   "COMPARISONS", "CONTROL_FACTORS",
  "COLUMN_ORDER", "NORMALIZATION"
  #"EXPERIMENTAL_CONDITIONS", "SAMPLE_CLASSIFICATIONS"
)


missing_sections <- setdiff(required_sections, names(EXPERIMENT_CONFIG))
if (length(missing_sections) > 0) {
  stop(sprintf("Missing required config sections: %s",
        paste(missing_sections, collapse = ", ")))
}

#if (validation_verbose) cat("[PASS] All required sections present\n\n")

# Validate configuration structure
stopifnot(
  "Missing required config sections" =
    all(required_sections %in% names(EXPERIMENT_CONFIG))
)

# @TODO Gotta figure out where to move this. Should it just be in the setup script?
sapply(category_names,
  function(category_name) {
    category_values <- EXPERIMENT_CONFIG$CATEGORIES[[category_name]]
    is_not_character <- !is.character(category_values)
    is_duplicated <- duplicated(category_values)
    if (any(is_not_character)) {
      stop(sprintf("Some category values in '%s' category are not character: %s",
        category_name,
        paste(category_values[is_not_character], collapse = ", ")
      ))
    }
    if (any(is_duplicated)) {
      stop(sprintf("Some category values in '%s' category are duplicated: %s",
        category_name,
        paste(category_values[is_duplicated], collapse = ", ")
      ))
    }
  }
)

categories_and_column_order_are_not_identical <- !identical(
  sort(category_names),
  sort(EXPERIMENT_CONFIG$COLUMN_ORDER)
)
if (categories_and_column_order_are_not_identical) {
  stop("Column order must include all category columns.")
}

lapply(names(EXPERIMENT_CONFIG$INVALID_COMBINATIONS),
  function(invalid_combination_condition) {
    combination_condition <- EXPERIMENT_CONFIG$INVALID_COMBINATIONS[[invalid_combination_condition]]
    invalid_category <- setdiff(
      all.vars(combination_condition),
      category_names
    )
    if (length(invalid_category) > 0) {
      stop(sprintf(
        "Invalid columns in INVALID_COMBINATIONS:\n%s",
          paste(invalid_category, collapse = ", ")
      ))

    }
  }
)

lapply(names(EXPERIMENT_CONFIG$CONTROL_FACTORS),
  function(control_factor_names) {
    control_factor_categories <- EXPERIMENT_CONFIG$CONTROL_FACTORS[[control_factor_names]]
    invalid_category <- setdiff(
      control_factor_categories,
      category_names
    )
    if (length(invalid_category) > 0) {
      stop(sprintf(
        "Invalid columns in CONTROL_FACTORS '%s':\n%s",
        control_factor_names, paste(invalid_category, collapse = ", ")
      ))

    }
  }
)

# Validate each section
#validate_category_values(
#  EXPERIMENT_CONFIG$CATEGORIES,
#  verbose = validation_verbose
#)

#validate_column_references(
#  categories = EXPERIMENT_CONFIG$CATEGORIES,
#  comparisons = EXPERIMENT_CONFIG$COMPARISONS,
#  control_factors = EXPERIMENT_CONFIG$CONTROL_FACTORS,
#  #conditions = EXPERIMENT_CONFIG$EXPERIMENTAL_CONDITIONS,
#  verbose = validation_verbose
#)
#
#validate_column_order(
#  categories = EXPERIMENT_CONFIG$CATEGORIES,
#  column_order = EXPERIMENT_CONFIG$COLUMN_ORDER,
#  verbose = validation_verbose
#)

#rm("experiment_id", "category_names")
cat("\n[VALIDATED] Experiment configuration loaded successfully\n")
