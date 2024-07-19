#Load packages required by script
package_to_check <- c("curl", "tidyverse","tools","R.utils")
package_was_loaded <- unlist(suppressPackageStartupMessages(lapply(package_to_check, library, character.only = TRUE, logical.return=TRUE, quietly = TRUE)))
if (length(lapply(package_to_check[!package_was_loaded], function(x){
  paste(x, "Package did not install")
})) == 0) print("All packages loaded.")

#Use subversion to get feature files from Rossi et al.
system('wsl svn export https://github.com/CEGRcode/2021-Rossi_Nature.git/trunk/02_References_and_Features_Files')

#Get folder that has features in it. 
feature_folder <- list.dirs(recursive = FALSE)[grepl("Features", list.dirs(recursive = FALSE))]

#Download origin timing data from https://www.sciencedirect.com/science/article/pii/S2211124713005834?via%3Dihub
hawkins_timing_url <- "https://ars.els-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx"
hawkins_timing <- paste(feature_folder, "hawkins-origins-timing.xlsx", sep = "/")
if (!file.exists(hawkins_timing)){
  curl_download(hawkins_timing_url, 
                hawkins_timing)
  
  print(paste(hawkins_timing, "file downloaded"))
} else {
  print(paste(hawkins_timing, "file exists"))
}

#Download called peaks for ORC in nocodazole from https://pubmed.ncbi.nlm.nih.gov/20351051/
eaton_orc_bed_url <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM424nnn/GSM424494/suppl/GSM424494_wt_G2_orc_chip_combined.bed.gz"
eaton_orc_bed <- paste(feature_folder, "G2_orc_chip.bed.gz", sep = "/")
if (!file.exists(eaton_orc_bed)){
  curl_download(eaton_orc_bed_url, 
                eaton_orc_bed)
  gunzip(eaton_orc_bed)
  print(paste(eaton_orc_bed, "file downloaded and extracted"))
} else {
  print(paste(eaton_orc_bed, "file exists"))
}

#Cant use system function because it requires password
# system('wsl rsync --stats -nv --update /mnt/c/Users/Luis/Projects/working-on-a-cluster/02_References_and_Features_Files/* luised94@luria.mit.edu:/home/luised94/data/feature-files/')
