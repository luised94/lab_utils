
cat(Sys.time())

.libPaths( c( "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2", .libPaths() ) )

.libPaths()

library(tidyverse, lib.loc = "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2")
library(stringr, lib.loc = "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2")

#Create folders if not present. Added conditional to run project startup script before
#assigning variables. Accounts for system OS (cluster or local).
# if (!dir.exists("./data")){
#   if(grepl("Windows", osVersion)){
#     renv::run("./scripts/2-ngs-project-startup.r")
#   } else if (grepl("Linux", osVersion)){
#     renv::run("../R-scripts/2-ngs-project-startup.r")
#   }
# }
# source("./scripts/2-ngs-project-startup.R")
cat(date())
#Folder to get directories for. 
data_folder <- "./data"
dir_names <- stringr::str_replace(list.dirs(data_folder, recursive = FALSE), pattern = data_folder, replace = "")
data_folder;dir_names

#Create the variables that will be assign the directory locations
#Some processing has to be done based on how they were created. Remove ./data/, -files and replace and - with _
dir_variables <- paste(paste(stringr::str_replace(stringr::str_replace(list.dirs(data_folder, recursive = FALSE), pattern = "./data/", replace = ""), pattern = "-files", replace=""), "_directory", sep = ""))
dir_variables <- stringr::str_replace(dir_variables, pattern = "-", "_")                       
#Create the data frame to iterate through
dir_df <- data.frame(names = dir_names, variables = dir_variables)

#Initialize file extensions in same order as dir_variables based on initial script
file_extension <- c(".bam", ".bw", ".fastq", ".html", ".fasta", ".bed", ".fastq")           

if (!all(unlist(lapply(dir_variables, exists)))) {
  
  dir_list <- c()
  for (i in 1:length(dir_df$names)){
    
    folder_variable <- paste(data_folder, dir_df$names[i], sep = "")
    assign(dir_df$variables[i], folder_variable)
    dir_list <- append(dir_list, folder_variable)
  }
  
  print("Directory variables created.")
  
} else {
  
  print("All variables assigned.")
}

print(dir_variables)
print(dir_df)
print(dir_list)

getwd()

cat(Sys.time())
print("assign-directory-variables.r complete")


# mapply(function(names, variables, dir_list = c(), folder){
# 
#   folder_variable <- paste(folder, names, sep = "")
#   assign(variables, folder_variable)
#   variables
#   dir_list <- append(dir_list, folder_variable)
# 
# }, names = dir_df$names, variables = dir_df$variables, folder = data_folder)

# sink("sessionInfo.txt")
# sessionInfo()
# sink()