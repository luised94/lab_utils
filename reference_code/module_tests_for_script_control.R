#!/usr/bin/env Rscript
source(file.path(Sys.getenv("HOME"), "lab_utils", "core_scripts", "functions_for_script_control.R"))

# Parse arguments and validate configurations
args <- parse_args(commandArgs(trailingOnly = TRUE))
experiment_id <- args[["experiment-id"]]
source(file.path("~/data", experiment_id, "documentation", 
                paste0(experiment_id, "_bmc_config.R")))
validate_configs(c("RUNTIME_CONFIG", "EXPERIMENT_CONFIG"))

print_config_settings(DEBUG_CONFIG)
