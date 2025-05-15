
if(interactive()) {
  message("Interactive job... sourcing configuration file.")
  script_configuration_path <- "~/lab_utils/core_scripts/script_configuration.R"
  stopifnot(
    "Script configuration file does not exist. Please copy the template." =
    file.exists(script_configuration_path)
  )
  source(script_configuration_path)
  message("Configuration file sourced...")
}

# Ensure the variables expected in the script were //
# defined in the configuration file. //
variables_to_check <- c(
  "EXPERIMENT_IDS",
  "EXPERIMENT_DIR",
  "CHROMOSOMES_TO_PLOT",
  "OUTPUT_FORMAT",
  "ACCEPT_CONFIGURATION",
  "SKIP_PACKAGE_CHECKS"
)
sapply(variables_to_check, function(variable){
  if (!exists(variable)) {
    stop(sprintf(
      fmt = "Variable %s not defined.\nAdjust configuration file.",
      variable
    ))
  }
})
message("All variables defined in the configuration file...")


#} else {
# add the args code?
#}
################################################################################
# Verify Required Libraries
################################################################################
# Add the packages that are used in the script.
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
if (!is.character(required_packages) || length(required_packages) == 0) {
  stop("required_packages must be a non-empty character vector")
}

if (!SKIP_PACKAGE_CHECKS) {
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        fmt = "Package '%s' is missing.\nPlease install using renv or base R.",
        pkg
      ))
    }
  }
  SKIP_PACKAGE_CHECKS <- TRUE
}

message("All required packages available...")
################################################################################
# Setup experiment-specific configuration path
################################################################################
number_of_experiments <- length(EXPERIMENT_DIR)
config_paths <- vector("character", length = number_of_experiments)
metadata_paths <- vector("character", length = number_of_experiments)
for (experiment_idx in seq_len(number_of_experiments)) {
  config_paths[experiment_idx] <- file.path(
    EXPERIMENT_DIR[experiment_idx], "documentation",
    paste0(EXPERIMENT_IDS[experiment_idx], "_bmc_config.R")
  )
  metadata_paths[experiment_idx] <- file.path(
    EXPERIMENT_DIR[experiment_idx], "documentation",
    paste0(EXPERIMENT_IDS[experiment_idx], "_sample_grid.csv")
  )
}
sapply(c(config_paths, metadata_paths), function(file_path){
  if (!file.exists(file_path)) {
    stop(sprintf(
      fmt = "File %s does not exist.\nRun setup_bmc_experiment.R script.",
      file_path
    ))
  }

})
