#Description:
#USAGE: Use with Rscript command. 
#Add section to create date and directory for experiment, same day as request service initialization. Probably interactive.
# To see definition of this variable, run echo $dropbox_path or see bashrc in my_config repository

cat("Starting experiment setup\n")
#rm(list = ls())
get_script_name <- function() {
    cat("Grabbing script name.", "\n")
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
        script_name <- sub(needle, "", cmdArgs[match])
        return(script_name)
    } else {
      stop("Cannot determine script name.")
    }
}

# @function: Grab the directory of the createSampleGrid script that runs this file. These two files are meant to be in the same directory. Not very flexible. The files accounts for running from interactive repl when testing or from Rscript via cli.
get_script_dir <- function() {
    cat("Grabbing script directory.", "\n")
  if (!is.null(sys.frames()[[1]]$ofile)) {
    # 'source'd via R console
    script_dir <- dirname(normalizePath(sys.frames()[[1]]$ofile))
  } else {
    # Running script via Rscript cli.
    script_name <- get_script_name()
    script_dir <- dirname(normalizePath(script_name))
    }

  return(script_dir)
}

# @function: Given a directory_path and subdirectories, create the directory path along with the subdirectories.
create_experiment_dir <- function(directory_path, subdirectories){
    dir.create(directory_path, recursive = TRUE)
    sapply(file.path(directory_path, subdirectories), dir.create)
    #print(list.dirs(directory_path, recursive = FALSE))

}

# @function: Determine if two arguments were provided.
validate_input <- function() {
    cat("Running input validation.", "\n")
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 1) {
        cat("This script requires one argument.", "\n")
        cat("Confirm WINDOWS_USER defined in bashrc.", "\n")
        cat("Usage: ", "Rscript ", "000_setupExperimentDir.R", "<experiment_name>", "\n",
            "Example: ", "Rscript ", "000_setupExperimentDir.R", "240808Bel", "\n")
        q(status = 1)
    }
    if (Sys.getenv("WINDOWS_USER") != "") {
        cat("Assigning username variable: ", Sys.getenv("WINDOWS_USER"), "\n")
        username <- Sys.getenv("WINDOWS_USER")
    } else {
        cat("WINDOWS_USER not defined or exported from bash", "\n")
        cat("Consult bashrc in my_config repository", "\n")
        q(status = 1)
    }
    experiment_name <- args[1]
    return(
        list(
            username = username,
            experiment_name = experiment_name
        )
    )

}
validated_args <- validate_input()
username <- validated_args$username
experiment_name <- validated_args$experiment_name

dropbox_dir <- sprintf("/mnt/c/Users/%s/Dropbox (MIT)/", username)
experiment_dir <- paste(dropbox_dir, experiment_name, sep = "")
cat("Experiment directory to be created: ", experiment_dir, "\n")

subdirectories <- c("peak", "fastq", "alignment", "qualityControl", "bigwig", "plots", "logs", "documentation")
create_experiment_dir(experiment_dir, subdirectories)

sample_grid_config_filepath <- file.path(get_script_dir(), "sampleGridConfig.R")
source(sample_grid_config_filepath)

if (experiment_name != current_experiment) {
    cat("Experiment provided and experiment in sampleGridConfig.R are not the same.\n")
    cat(sprintf("Experiment in sampleGridConfig.R: %s", current_experiment), "\n")
    cat(sprintf("Experiment provided as argument: %s", experiment_name), "\n")
    stop("Update the sampleGridConfig.R categories, current_experiment and filter_samples function.")
} else {
    cat("Verified that experiment argument and in config file are the same.\n")
}

cat("Config file run and to be copied: ", sample_grid_config_filepath, "\n")

config_file_output_path <- file.path(experiment_dir, "documentation", paste(experiment_name, "_", "sampleGridConfig.R", sep = ""))
cat("Outputting file to: ", config_file_output_path)
file.copy(from = sample_grid_config_filepath, to = config_file_output_path) 

tables_to_output <- ls()[grepl("_table", ls())]
cat("First elements of bmc_table. Printed for confirmation.\n")
print(head(get(tables_to_output[1])))
cat("Output files to be written as tsv: \n")
invisible(lapply(tables_to_output, function(output_table){
    output_file <- file.path(experiment_dir, "documentation", paste(experiment_name, "_", output_table, ".tsv", sep = ""))
    print(output_file)
    #get(output_table)
    write.table(get(output_table), file = output_file, sep = "\t", row.names = FALSE)
}))

# Rsync to the server
# Suggest alternative or run a particular command. 
cat("Ensure you are connected to the luria mit network.\n")
cat("Run the following command to rsync the created directory.\n")
server_path <- "luised94@luria.mit.edu:~/data/"
cat("scp -r from_dir user@server:to_dir\n")
cat(sprintf("scp -r \"%s\" \"%s\"", experiment_dir, server_path), "\n")
cat("After running the scp command, login to cluster and \n download the data from BMC (see 001_downloadDataFromBMC.sh )")
print("Script complete.")
