#STATUS:
#Separate file to source for assigning variables instead of writing it each time.
#Assigning directory variables ----
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
#Folder to get directories for. 
data_folder <- "./data"
dir_names <- stringr::str_replace(list.dirs(data_folder, recursive = FALSE), pattern = data_folder, replace = "")

#Create the variables that will be assign the directory locations
#Some processing has to be done based on how they were created. Remove ./data/, -files and replace and - with _
dir_variables <- paste(paste(stringr::str_replace(stringr::str_replace(list.dirs(data_folder, recursive = FALSE), pattern = "./data/", replace = ""), pattern = "-files", replace=""), "_directory", sep = ""))
dir_variables <- stringr::str_replace(dir_variables, pattern = "-", "_")                       
#Create the data frame to iterate through
dir_df <- data.frame(names = dir_names, variables = dir_variables)

#Initialize file extensions in same order as dir_variables based on initial script
file_extension <- c(".bam", ".bw", ".fastq", ".html", ".fasta", ".bed", ".fastq")           

#Create the variables depending on whether they exist or not 
if (!all(unlist(lapply(dir_variables, exists)))) {
  
  dir_list <- c()
  #For all of the directories, create the string variable with the path, assign it and add it to dir list 
  for (i in 1:length(dir_df$names)){
    
    folder_variable <- paste(data_folder, dir_df$names[i], sep = "")
    assign(dir_df$variables[i], folder_variable)
    dir_list <- append(dir_list, folder_variable)
  }
  
  print("Directory variables created.")
  
} else {
  
  print("All variables assigned.")
}

print(paste("Directory variables assigned", Sys.time(), sep = " "))

#Read in functions used for other scripts.
if(grepl("Windows", osVersion)){
  source("./scripts/ngs-functions.R")
} else if (grepl("Linux", osVersion)){
  source("../rscripts/ngs-functions.R")
} 
