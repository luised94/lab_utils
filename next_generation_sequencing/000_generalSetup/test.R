
dropbox_dir <- "/mnt/c/Users/Luis/Dropbox (MIT)/"
args <- commandArgs(trailingOnly = TRUE)
experiment_dir <- paste(dropbox_dir, args[1], sep = "")
print(experiment_dir)
dir.create(experiment_dir)
"mkdir -p"
system("mkdir -p""$dir"/{peak,fastq,alignment,qualityControl,bigwig,plots,logs,documentation}
