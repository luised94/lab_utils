rm(list = ls())
#Description:
#USAGE: Use with Rscript command. 
#Add section to create date and directory for experiment, same day as request service initialization. Probably interactive.
# To see definition of this variable, run echo $dropbox_path or see bashrc in my_config repository

# @function: Grab the directory of the createSampleGrid script that runs this file. These two files are meant to be in the same directory. Not very flexible. The files accounts for running from interactive repl when testing or from Rscript via cli.
get_script_dir <- function() {
  if (!is.null(sys.frames()[[1]]$ofile)) {
    # 'source'd via R console
    script_dir <- dirname(normalizePath(sys.frames()[[1]]$ofile))
  } else {
    # Rscript
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
      script_dir <- dirname(normalizePath(sub(needle, "", cmdArgs[match])))
    } else {
      stop("Cannot determine script directory")
    }
  }
  return(script_dir)
}

# @function: Given a directory_path and subdirectories, create the directory path along with the subdirectories.
create_experiment_dir <- function(directory_path, subdirectories){
    dir.create(directory_path, recursive = TRUE)
    sapply(file.path(directory_path, subdirectories), dir.create)
    print(list.dirs(directory_path, recursive = FALSE))

}
args <- commandArgs(trailingOnly = TRUE)

dropbox_dir <- sprintf("/mnt/c/Users/%s/Dropbox (MIT)/", args[2])
experiment_dir <- paste(dropbox_dir, args[1], sep = "")
print(experiment_dir)

subdirectories <- c("peak", "fastq", "alignment", "qualityControl", "bigwig", "plots", "logs", "documentation")
#create_experiment_dir(experiment_dir, subdirectories)
source(file.path(get_script_dir(), "sampleGridConfig.R"))
print(file.path(get_script_dir(), "sampleGridConfig.R"))

print(file.path(experiment_dir, "documentation", "sampleGridConfig.R"))
config_file_output_path <- file.copy(file.path(experiment_dir, "documentation", paste(args[1], "_", )))
tables_to_output <- ls()[grepl("_table", ls())]
print(get(tables_to_output[1]))
lapply(tables_to_output, function(output_table){
    output_file <- file.path(experiment_dir, "documentation", paste(args[1], "_", output_table, ".tsv", sep = ""))
    print(output_file)
    #get(output_table)
#    write.table(get(output_table), file = output_file, sep = "\t", row.names = FALSE)
})
