#CATEGORIES = list(
#    rescue_allele = c("WT", "4R"),
#    suppressor_allele = c("NONE", "1EK","3PL","4PS","5EK","6EK","TGE","Cdc6OE"),
#    cell_cycle_arrest = c("ALPHA", "NOCO"),
#    antibody = c("Input", "HM1108", "UM174"),
#    repeats = c("1", "2")
#)
#experiment_id <- c("250207Bel", "250324Bel")
#if (length(EXPERIMENT_CONFIG$INVALID_COMBINATIONS) > 0) {
#    invalid_idx <- Reduce(
#        `|`,
#        lapply(EXPERIMENT_CONFIG$INVALID_COMBINATIONS, eval, envir = metadata)
#    )
#    metadata <- subset(metadata, !invalid_idx)
#}

# Apply experimental conditions
#valid_idx <- Reduce(
#    `|`,
#    lapply(EXPERIMENT_CONFIG$EXPERIMENTAL_CONDITIONS, eval, envir = metadata)
#)
#metadata <- subset(metadata, valid_idx)
SUBSET_REPRODUCIBLE_SAMPLES <- list(
    inputs = quote(
      antibody == "Input"
    ),
    orc_repeats_in_noco = quote(
      antibody == "ORC" & cell_cycle_arrest == "NOCO" &
      !(repeats != "2" & experiment_id == "250324Bel")
    ),
    orc_repeats_in_alpha = quote(
      antibody == "ORC" & cell_cycle_arrest == "ALPHA" &
      ((repeats == "1" & experiment_id == "250324Bel") |
      (repeats == "2" & experiment_id == "250207Bel"))

    ),
    mcm_repeats_in_alpha = quote(
      repeats == "2" & antibody == "UM174" & cell_cycle_arrest == "ALPHA"
    ),
    mcm_repeats_in_noco = quote(
      !(repeats == "1" & antibody == "UM174" & suppressor_allele %in% c("5EK", "6EK"))
    ),
    cdc6_repeats = quote(
      suppressor_allele == "Cdc6OE"
    )
)
