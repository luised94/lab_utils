################################################################################
# BMC ChIP-seq Configuration Setup
# Author: Luis | Date: 2025-10-27 | Version: 3.0.0 | Experiment: 250715Bel
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
#   No outputs produced by file or script.
#   Scripts (such as setup_bmc_experiment.R) source this script.
#   Setups environment with the values.
#
# DESCRIPTION:
#   Describe the experiment.
################################################################################
# !! Update EXPERIMENT_CONFIG if starting a new experiment.
# 1) Update description and metadata.
# 2) Set the categories of the experiment. This will produce all of the combinations.
# 3) Use quote expressions and logical conditions to write what should not be in the metadata out of all the combinations.
# 4) Update column order to sort the metadata after filtering.
# 5) Specify particular genome track plots to create using comparison expressions to filter dataset.
EXPERIMENT_CONFIG <- list(
  METADATA = list(
    # YYMMDDBel
    EXPERIMENT_ID = "250821Bel",
    EXPECTED_SAMPLES_SETUP = 5,
    EXPECTED_SAMPLES_POST = 5,
    VERSION = "1.0.0",
    PROJECT_ID = "project_002"
  ),

  # Always include antibody category for chip experiments.
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

  COMPARISONS = list(
  ),

  # Set the options for which bigwig file is used to create genome track plots.
  # They look mostly the same.
  # @TODO Remove this section of code?
  NORMALIZATION = list(
    methods = c("CPM", "BPM", "RPGC", "RPKM"),
    active = "CPM"
  ),

  # 
  TARGET_COMPARISON_COLUMNS = c(
    "rescue_allele", "timepoints"
  ),

  COLUMNS_TO_EXCLUDE_FROM_REPLICATE_DETERMINATION = c(
    "bigwig_file_paths",
    "experiment_id",
    "full_name",
    "repeats",
    "sample_ids",
    "sample_type",
    "short_name"
  ),

  METADATA_COLUMNS_TO_EXCLUDE = c(
      "sample_type", "sample_ids",
      "bigwig_file_paths", "full_name",
      "short_name"
  ),

  # Optional: specify a row reordering correction if a sample was misplaced during submission
  # Use NULL (or omit) when no correction is needed.
  # Example: move sample from row 12 to position 24 (1-based indexing)
  ROW_ORDER_CORRECTION = list(
    from_row = 12,
    to_row   = 24
  ),

  # Optional: specify a sample that was dropped during prep and is missing from core data
  # Use NULL when no dropout occurred.
  # Identify the sample using a logical condition on metadata columns.
  # If you add, do not forget to adjust the EXPECTED_NUMBER_OF_SAMPLES above.
  SAMPLE_DROPOUT_CONDITION = quote(
    suppressor_allele == "4PS" & 
    antibody == "ORC" & 
    cell_cycle_treatment == "ALPHA" & 
    rescue_allele == "4R"
  )

  ## ROW_ORDER_CORRECTION = quote(
  #  suppressor_allele == "4PS" & antibody == "ORC" & 
  #  cell_cycle_treatment == "ALPHA" & rescue_allele == "4R"
  #)

  #SAMPLE_CLASSIFICATIONS = list(
  #  is_positive = quote(orc_phenotype == "WT" & antibody == "ORC"),
  #  is_negative = quote(orc_phenotype == "orc1-161" & antibody == "ORC" & temperature == "37")
  #),

)
# @TODO: Add once genome tracks is fixed.
#metadata <- subset(metadata, valid_idx)
# Set to NULL
# Should probably move to configuration
#columns_to_show <- c(
#  "rescue_allele", "suppressor_allele",
#  "cell_cycle_arrest", "repeats",
#  "antibody", "experiment_id"
#  )
## Rename to focus list or something
#message("Sourcing subset expressions for reproducible examples...")
##EXPECTED_NUMBER_OF_REPRODUCIBLE_SAMPLES <-
#SUBSET_REPRODUCIBLE_SAMPLES <- list(
#    inputs = quote(
#      antibody == "Input"
#    ),
#    orc_repeats_in_noco = quote(
#      antibody == "HM1108" & cell_cycle_arrest == "NOCO" &
#      !(repeats == "2" & experiment_id == "250324Bel")
#    ),
#    orc_repeats_in_alpha = quote(
#      antibody == "HM1108" & cell_cycle_arrest == "ALPHA" &
#      ((repeats == "1" & experiment_id == "250324Bel") |
#      (repeats == "2" & experiment_id == "250207Bel"))
#    ),
#    mcm_repeats_in_alpha = quote(
#      repeats == "2" & antibody == "UM174" & cell_cycle_arrest == "ALPHA"
#    ),
#    mcm_repeats_in_noco = quote(
#      antibody == "UM174" & cell_cycle_arrest == "NOCO" &
#      !(repeats == "1" & antibody == "UM174" & suppressor_allele %in% c("5EK", "6EK"))
#    ),
#    cdc6_repeats = quote(
#      suppressor_allele == "Cdc6OE"
#    )
#)
