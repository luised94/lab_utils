
dropbox_dir <- "/mnt/c/Users/Luis/Dropbox (MIT)/"
args <- commandArgs(trailingOnly = TRUE)
experiment_dir <- paste(dropbox_dir, args[1], sep = "")
print(experiment_dir)
dir.create(experiment_dir)
dir_creation_command <- paste("mkdir -p '", experiment_dir, "'/{peak,fastq,alignment,qualityControl,bigwig,plots,logs,documentation}", sep = "")
print(dir_creation_command)
cat("Command:", dir_creation_command, "\n")
#system(dir_creation_command)
