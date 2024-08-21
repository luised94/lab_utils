rm(list = ls())
#Description:
#USAGE: Use with Rscript command. 
#Add section to create date and directory for experiment, same day as request service initialization. Probably interactive.
# To see definition of this variable, run echo $dropbox_path or see bashrc in my_config repository

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
    if (length(args) < 2) {
        cat("This script requires two arguments.", "\n")
        cat("Usage: ", "Rscript ", "000_setupExperimentDir.R", "<username>", "<experiment_name>", "\n",
            "Example: ", "Rscript ", "000_setupExperimentDir.R", "${windows_user}", "240808Bel", "\n")
        q(status = 1)
    }
    username <- args[1]
    experiment_name <- args[2]
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
cat("Experiment directory to be created: ", experiment_dir)

subdirectories <- c("peak", "fastq", "alignment", "qualityControl", "bigwig", "plots", "logs", "documentation")
create_experiment_dir(experiment_dir, subdirectories)

sample_grid_config_filepath <- file.path(get_script_dir(), "sampleGridConfig.R")
source(sample_grid_config_filepath)
cat("Config file run and to be copied: ", sample_grid_config_filepath)

config_file_output_path <- file.path(experiment_dir, "documentation", paste(experiment_name, "_", "sampleGridConfig.R", sep = ""))
cat("Outputting file to: ", config_file_output_path)
file.copy(from = sample_grid_config_filepath, to = config_file_output_path) 

tables_to_output <- ls()[grepl("_table", ls())]
print(head(get(tables_to_output[1])))
lapply(tables_to_output, function(output_table){
    output_file <- file.path(experiment_dir, "documentation", paste(experiment_name, "_", output_table, ".tsv", sep = ""))
    print(output_file)
    #get(output_table)
    write.table(get(output_table), file = output_file, sep = "\t", row.names = FALSE)
})

# Rsync to the server
# Suggest alternative or run a particular command. 
cat("Ensure you are connected to the luria mit network.")
cat("Run the following command to rsync the created directory.")
server_path <- "luised94@luria.mit.edu:~/data/"
cat(sprintf("scp -r %s %s", experiment_dir, server_path))
print("Script complete.")
