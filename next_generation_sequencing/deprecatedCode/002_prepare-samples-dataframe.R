#STATUS:
#Gets relative paths for the directories in data folder and assigns them to variables
#Based on the operating system since the relative paths are different. 
if(grepl("Windows", osVersion)){
  source("./scripts/assign-directory-variables.R")
} else if (grepl("Linux", osVersion)){
  source("../rscripts/assign-directory-variables.R")
}

#Prepare the dataframe with the input/output files ----

#Get the genome index directory
index_dir <- list.dirs(genome_directory)[grepl("index", list.dirs(genome_directory))]

#Get all fastq files in fastq_directory obtained from assign-directory-variable.R 
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
  } 
}, dir_list, file_extension) %>% discard(is.null) %>% data.frame(sample_names = sample_name, original_file = fastq_file_list, .)
#dim(sample_paths_df)

#Add columns for variables that can be used to assign different data types 
sample_paths_df <- sample_paths_df %>% mutate(fq_var = paste(sample_names, "_fq", sep = ""),
                                              bw_var = paste(sample_names, "_bw", sep = ""),
                                              track_var = paste(sample_names, "_track", sep = ""))
#Add the columns for sorted bam files and their indexes 
sample_paths_df <- sample_paths_df %>% mutate(sorted_bam = paste(tools::file_path_sans_ext(..data.alignment.files), "_sorted.bam", sep = ''),
                                              bam_index = paste(tools::file_path_sans_ext(..data.alignment.files), "_sorted.bam.bai", sep = ''))

print(paste("Dataframe with paths and variables to all files created", Sys.time(), sep = " "))
