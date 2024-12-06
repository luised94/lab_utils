#STATUS:
# source("./scripts/4-prepare-genome.r")
library(ShortRead)
library(QuasR)
library(R.utils)
require(tidyverse)

#Assigning directory variables ----
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

print(dir_variables)
print(dir_df)
print(dir_list)

getwd()

Sys.time()
print("assign-directory-variables.r complete")

#Get toy data if user wants it and if fastq-files folder is empty ----
#download_data <- readline(prompt = "Would you like to download toy data? (T or F): ")

# if (as.logical(download_data) &
#     (length(list.files(fastq_directory)) == 0)) {
#   toy_data <-
#     paste(fastq_directory, "eaton_toy_data.fastq.gz", sep = "/")
#   library(curl)
#   print("Getting some toy data")
#   curl_download(
#     "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR034/SRR034475/SRR034475.fastq.gz",
#     toy_data
#   )
#   gunzip(toy_data)
#   print(paste(toy_data, "file downloaded and extracted"))
# 
#   print("Splitting file.")
#   toy_data <- str_replace(toy_data, pattern = ".gz", "")
# 
#   stream <-
#     open(FastqStreamer(toy_data, n = 500000, readerBlockSize = 1000))
# 
#   counter <- 1
#   repeat {
#     tmp_fq <- yield(stream)
#     destination_file <-
#       paste(
#         str_replace(toy_data, pattern = ".fastq$", replacement = ""),
#         "_chunk",
#         counter,
#         ".fastq",
#         sep = ""
#       )
#     #If input is 0, stop the function
#     if (length(tmp_fq) == 0) {
#       print("Split file")
#       break
#     }
#     writeFastq(tmp_fq,
#                destination_file,
#                "a",
#                compress = FALSE ,
#                qualityType = "FastqQuality")
#     counter <- counter + 1
#   }
#   close(stream)
#   file.remove(toy_data)
# 
# } else {
#   print("Following files are in the directory:")
#   list.files(fastq_directory)
# }

#Get the index directory
index_dir <- list.dirs(genome_directory)[grepl("index", list.dirs(genome_directory))]

#Get all fastq files in fastq_directory obtained from assign-directory-variable.R ----
fastq_file_list <- list.files(fastq_directory, pattern = ".fastq$", full.names = TRUE, recursive = FALSE)[!grepl("unmapped", list.files(fastq_directory, pattern = ".fastq$", full.names = TRUE, recursive = FALSE))]

#Create variable for genome
genome_file_path <- list.files(path = genome_directory, pattern = ".fasta", full.names = TRUE)

#Obtain the sample name for each fastq file
sample_name <- str_replace(basename(fastq_file_list), pattern = paste(".", tools::file_ext(fastq_file_list), sep = ""), replacement = "")   

#Create the dataframe that has all paths required to preprocessing and alignment. 
sample_paths_df <- mapply(function(dir_, ext_){
  #Create the name for each file depending on the folder it will be outputted to. 
  if (grepl("processed", dir_)){
    paste(dir_, paste("processed_", basename(sample_name), ext_, sep = ""), sep = "/")
  } else if (!grepl("genome|fastqc", dir_)){
    paste(dir_, paste(basename(sample_name), ext_, sep = ""), sep = "/")
  } else {
  }
}, dir_list, file_extension) %>% discard(is.null) %>% data.frame(sample_names = sample_name, original_file = fastq_file_list, .)
dim(sample_paths_df)
#Add the variables that can be used to assign different data types 
sample_paths_df <- sample_paths_df %>% mutate(fq_var = paste(sample_names, "_fq", sep = ""),
                                              bw_var = paste(sample_names, "_bw", sep = ""),
                                              track_var = paste(sample_names, "_track", sep = ""))

#Load in myfilter and trim function and my align,sort and index function.
source("/home/luised94/data/rscripts/ngs-functions.r")

#Try it for one.
# myFilterAndTrim(sample_paths_df$original_file[1], destination_file = sample_paths_df$..data.processed.fastq.files[1])

mapply(myFilterAndTrim, filelocation = sample_paths_df$original_file, destination_file = sample_paths_df$..data.processed.fastq.files)

#Set the alignment parameters for Bowtie2
alignment_params <- "-q --mp 4 --met-stderr"

#Try it for one.
# library(Rbowtie2)
# bowtie2_samtools(bt2Index = file.path(tools::file_path_as_absolute(index_dir), "sacCer3"),
#                  output = str_replace(sample_paths_df$..data.alignment.files[1], pattern = ".bam", replacement = ""),
#                  outputType = str_replace(".bam", pattern = ".", replacement = ""),
#                  seq1 = sample_paths_df$..data.processed.fastq.files[1],
#                  seq2 = NULL,
#                  bamFile = NULL,
#                  alignment_params)
# 
# sortBam(sample_paths_df$..data.alignment.files[1], destination = paste(tools::file_path_sans_ext(sample_paths_df$..data.alignment.files[1]), "_sorted", sep = ''))
# indexBam(paste(tools::file_path_sans_ext(sample_paths_df$..data.alignment.files[1]), "_sorted.bam", sep = ''))

#Try my function for one file 
alignSortandIndex(index_dir = index_dir, index_prefix = "sacCer3", fastq_file = sample_paths_df$..data.processed.fastq.files[1],
                  alignment_file = sample_paths_df$..data.alignment.files[1], alignment_params = alignment_params)

mapply(alignSortandIndex, index_dir=index_dir, index_prefix="sacCer3", fastq_file = sample_paths_df$..data.processed.fastq.files, 
       alignment_file = sample_paths_df$..data.alignment.files, alignment_params = alignment_params)

expected_files_present <- list.files(unique(dirname(sample_paths_df$..data.alignment.files)), pattern = "_sorted.bam$", full.names = TRUE) == str_replace(sample_paths_df$..data.alignment.files, pattern = ".bam", "_sorted.bam")

  if(all(expected_files_present)){
    print("All expected sorted alignment files present.")
    list.files(unique(dirname(sample_paths_df$..data.alignment.files)), pattern = "_sorted.bam$", full.names = TRUE)
  } else {
    paste(sample_paths_df$..data.alignment.files[!expected_files_present], "not present", sep = " ")
    paste(sample_paths_df$..data.alignment.files[expected_files_present], "created", sep = " ")
  }

 
#Will attempt to parallelize. Have to delete the files created from previous lines
# library(snow)
# z=vector('list',4)
# z=1:4
# system.time(lapply(z,function(x) Sys.sleep(1)))
# cl<-makeCluster(detectCores()-1,type="SOCK")
# system.time(clusterApply(cl, z,function(x) Sys.sleep(1)))
# stopCluster(cl)



# library(doParallel)
# library(parallel)
# 
# no_cores <- 3
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)
# original_file_list <- sample_paths_df$original_file
# processed_file_list <- sample_paths_df$..data.processed.fastq.files
# export_to_cluster <- c('myFilterAndTrim', 'processed_file_list')
# clusterExport(cl,export_to_cluster)
# 
# parLapply(cl, original_file_list, function(x){
#   myFilterAndTrim(filelocation = x, destination_file = processed_file_list)
# })
# 
# # mapply(myFilterAndTrim, filelocation = sample_paths_df$original_file, destination_file = sample_paths_df$..data.processed.fastq.files)
# stopCluster(cl)

download_data <- as.logical(readline(prompt = "Were you running toy data? (T or F): "))
remove_toy_data <- as.logical(readline(prompt = "Would you like to delete the toy-data generated files? (T or F): "))
also_fastq <- readline(prompt = "Would you also like to delete download toy data? (T or F) Keep to continue testing pipeline.")

if (download_data & remove_toy_data){
  unlink(paste(processed_fastq_directory, "/*", sep = ""))
  unlink(paste(alignment_directory, "/*", sep = ""))
  
}






