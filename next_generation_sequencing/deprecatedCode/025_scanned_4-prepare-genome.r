
#Assigning directory variables ----
Sys.time()

.libPaths( c( "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2", .libPaths() ) )

.libPaths()

library(tidyverse, lib.loc = "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2")
library(stringr, lib.loc = "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2")

#Create folders if not present. Added conditional to run project startup script before
# #assigning variables. Accounts for system OS (cluster or local).
# if (!dir.exists("./data")){
#   if(grepl("Windows", osVersion)){
#     renv::run("./scripts/2-ngs-project-startup.r")
#   } else if (grepl("Linux", osVersion)){
#     renv::run("../R-scripts/2-ngs-project-startup.r")
#   }
# }
# source("./scripts/2-ngs-project-startup.R")
Sys.time()
#Folder to get directories for. 
data_folder <- "./data"
dir_names <- stringr::str_replace(list.dirs(data_folder, recursive = FALSE), pattern = data_folder, replace = "")
print(data_folder);print(dir_names)

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
#Initialize directory variables 
#Runs project startup script if folders dont exist. File location is based on operating system.
# if (!dir.exists("./data")){
#   if(grepl("Windows", osVersion)){
#     renv::run("./scripts/3-assign-directory-variables.r")
#   } else if (grepl("Linux", osVersion)){
#     renv::run("../R-scripts/3-assign-directory-variables.r")
#   }
# }

#Preparing genome ----

genome_file_path <- paste(genome_directory, "sacCer3.fasta", sep = "/")
index_dir <- paste(genome_directory, "sacCer3-index", sep = "/")
print(genome_file_path)
print(index_dir)


if (!file.exists(genome_file_path)){
  
  library(BSgenome.Scerevisiae.UCSC.sacCer3)
  sacCer3 <- BSgenome.Scerevisiae.UCSC.sacCer3
  genome_file_path <- paste(genome_directory, "sacCer3.fasta", sep = "/")
  index_dir <- paste(genome_directory, paste(metadata(sacCer3)$genome, "-index", sep = ""), sep = "/")
  #Export to local folder
  export(sacCer3, genome_file_path, format = "fasta")
  
  #Create directory for index
  dir.create(index_dir)
  print(paste(index_dir, "directory created", sep = " "))
  
  #Using bowtie2 to build index for sacCer3
  #Have to download python and perl
  #https://www.bioconductor.org/packages/release/bioc/manuals/Rbowtie2/man/Rbowtie2.pdf
  #Works after installing python and perl. See https://www.biostars.org/p/9507828/#9508926 and https://github.com/wzthu/Rbowtie2 
  library(Rbowtie2)
  bt2_index_path <- paste(index_dir, "sacCer3", sep = "/") 
  
  #Make index 
  cmdout <- bowtie2_build(references = tools::file_path_as_absolute(genome_file_path), 
                          bt2Index = file.path(tools::file_path_as_absolute(index_dir), "sacCer3"),
                          "--threads 4 --verbose", overwrite=FALSE)
  
    
  
  #Remove at the end since you cant align to BSgenome object 
  rm(sacCer3)
  
} else if (!dir.exists(index_dir)){
  
  print(paste(genome_file_path, "file exists but index folder is not present", sep = " "))
  dir.create(index_dir)
  print(paste(index_dir, "directory created", sep = " "))
  
  library(Rbowtie2)
  bt2_index_path <- paste(index_dir, "sacCer3", sep = "/")
  cmdout <- bowtie2_build(references = tools::file_path_as_absolute(genome_file_path), 
                          bt2Index = file.path(tools::file_path_as_absolute(index_dir), "sacCer3"),
                          "--threads 4 --verbose", overwrite=FALSE)
  
  print(paste("Indexed ", genome_file_path, sep = ""))
  
} else if (!(length(list.files(index_dir) == 6))){
  print(paste(index_dir, "directory exists but there arent six files in it", sep = " "))
  
  library(Rbowtie2)
  bt2_index_path <- paste(index_dir, "sacCer3", sep = "/")
  cmdout <- bowtie2_build(references = tools::file_path_as_absolute(genome_file_path), 
                          bt2Index = file.path(tools::file_path_as_absolute(index_dir), "sacCer3"),
                          "--threads 4 --verbose", overwrite=FALSE)
  print(paste("Indexed ", genome_file_path, sep = ""))
        
} else {
  
  print(paste(genome_file_path, "file exists", sep = " "))
  print(paste(index_dir, "directory exists", sep = " "))
  print(paste("Exactly six files present in", index_dir, sep = ""))
  
}

print("Verified if reference genome is ready")

Sys.time()
print("prepare-genome.r complete")


