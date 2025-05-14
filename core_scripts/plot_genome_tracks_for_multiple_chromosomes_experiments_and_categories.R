
if(interactive()) {
  message("Interactive job... sourcing configuration file.")
  source("~/lab_utils/core_scripts/script_configuration.R")
  message("Configuration file sourced...")
}
#variables_to_check <- list(
#  "experiment_dir"
#)
#} else {
# add the args code?
#}
################################################################################
# Verify Required Libraries
################################################################################
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
if (!is.character(required_packages) || length(required_packages) == 0) {
  stop("required_packages must be a non-empty character vector")
}

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      fmt = "Package '%s' is missing.\nPlease install using renv or base R.",
      pkg
    ))
  }
}

message("All required packages available...")
