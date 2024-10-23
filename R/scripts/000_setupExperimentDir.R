#Description:
#USAGE: Use with Rscript command. 
#Add section to create date and directory for experiment, same day as request service initialization. Probably interactive.
# To see definition of this variable, run echo $dropbox_path or see bashrc in my_config repository
source("~/lab_utils/R/init.R")
cat("Starting experiment setup\n")
OUTPUT_TO_FILE <- FALSE
main <- function() {
    validated_args <- validate_input()
    username <- validated_args$username
    experiment_name <- validated_args$experiment_name

    dropbox_dir <- sprintf("/mnt/c/Users/%s/Dropbox (MIT)/", username)
    experiment_dir <- paste(dropbox_dir, experiment_name, sep = "")
    cat("Experiment directory to be created: ", experiment_dir, "\n")

    sample_grid_config_filepath <- file.path(get_script_dir(), "sampleGridConfig.R")
    source(sample_grid_config_filepath) # Initializes sample_config_output and current_experiment

    if (experiment_name != current_experiment) {
        cat("Experiment provided and experiment in sampleGridConfig.R are not the same.\n")
        cat(sprintf("Experiment in sampleGridConfig.R: %s", current_experiment), "\n")
        cat(sprintf("Experiment provided as argument: %s", experiment_name), "\n")
        stop("Update the sampleGridConfig.R categories, current_experiment and filter_samples function.")
    } else {
        cat("Verified that experiment argument and in config file are the same.\n")
    }

    subdirectories <- c("peak", "fastq", "alignment", "qualityControl", "bigwig", "plots", "logs", "documentation")
    create_experiment_dir(experiment_dir, subdirectories)

    cat("Config file run and to be copied: ", sample_grid_config_filepath, "\n")
    config_file_output_path <- file.path(experiment_dir, "documentation", paste(experiment_name, "_", "sampleGridConfig.R", sep = ""))
    cat("Outputting file to:", config_file_output_path,"\n" )
    if (OUTPUT_TO_FILE){
        file.copy(from = sample_grid_config_filepath, to = config_file_output_path) 
    } else {
        cat(sprintf("Skip writing %s to file. Modify OUTPUT_TO_FILE value.\n", sample_grid_config_filepath))
    }
    #print("Files currently loaded")


    output_tables_in_list(experiment_directory = experiment_dir, sample_config_output, OUTPUT_TABLE = OUTPUT_TO_FILE )

    # Rsync to the server
    # Suggest alternative or run a particular command. 
   scp_command_reminder() 

}
# @function: Given a directory_path and subdirectories, create the directory path along with the subdirectories.

# @function: Determine if arguments were provided.
validate_input <- function() {
    cat("Running input validation.", "\n")
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 1) {
        cat("This script requires one argument.", "\n")
        cat("Confirm WINDOWS_USER defined in bashrc.", "\n")
        cat("Usage: ", "Rscript ", "000_setupExperimentDir.R", "<experiment_name>", "\n",
            "Example: ", "Rscript ", "000_setupExperimentDir.R", "240808Bel", "\n")
        stop("Provide one argument corresponding to experiment ID from BMC.")
    }
    if (Sys.getenv("WINDOWS_USER") != "") {
        cat("Assigning username variable: ", Sys.getenv("WINDOWS_USER"), "\n")
        username <- Sys.getenv("WINDOWS_USER")
    } else {
        cat("WINDOWS_USER not defined or exported from bash", "\n")
        cat("Consult bashrc in my_config repository", "\n")
        stop()
    }
    experiment_name <- args[1]
    return(
        list(
            username = username,
            experiment_name = experiment_name
        )
    )
}

if(!interactive()) {
    main()
}
