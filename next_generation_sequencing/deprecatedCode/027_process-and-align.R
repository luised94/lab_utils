#STATUS: REMOVE.
source("./scripts/assign-directory-variables.R")
#Prepare the dataframe with the input/output files ----

#Get the index directory
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
dim(sample_paths_df)
#Add the variables that can be used to assign different data types 
sample_paths_df <- sample_paths_df %>% mutate(fq_var = paste(sample_names, "_fq", sep = ""),
                                              bw_var = paste(sample_names, "_bw", sep = ""),
                                              track_var = paste(sample_names, "_track", sep = ""))

#Load in myfilter and trim function and my align,sort and index function
source("/home/luised94/data/rscripts/ngs-functions.r")

#Load the packages required by my functions ----
package_to_check <- c("Rbowtie2", "ShortRead")
package_was_loaded <- unlist(lapply(package_to_check, library, character.only = TRUE, logical.return=TRUE)) 
if (length(lapply(package_to_check[!package_was_loaded], function(x){
  paste(x, "Package did not install ")
})) == 0) print("All packages loaded.")


#Run my functions on one sample to sbatch submission ----
# myFilterAndTrim(sample_paths_df$original_file[1], destination_file = sample_paths_df$..data.processed.fastq.files[1])
#Run Rbowtie2 to align and index my file 
# alignSortandIndex(index_dir = index_dir, index_prefix = "sacCer3", fastq_file = sample_paths_df$..data.processed.fastq.files[1],
#                   alignment_file = sample_paths_df$..data.alignment.files[1], alignment_params = alignment_params)
 

#Align files ----

#Set alignment parameters for Rbowtie2
alignment_params <- "-q --mp 4"

#Use mapply to filter,align,sort and index all of the files by referencing the different columns of the data frame. 
mapply(myFilterAndTrim, filelocation = sample_paths_df$original_file, destination_file = sample_paths_df$..data.processed.fastq.files)

mapply(alignSortandIndex, index_dir=index_dir, index_prefix="sacCer3", fastq_file = sample_paths_df$..data.processed.fastq.files, 
       alignment_file = sample_paths_df$..data.alignment.files, alignment_params = alignment_params)

#Determine if all files were created 
expected_files_present <- list.files(unique(dirname(sample_paths_df$..data.alignment.files)), pattern = "_sorted.bam$", full.names = TRUE) == str_replace(sample_paths_df$..data.alignment.files, pattern = ".bam", "_sorted.bam")

if(all(expected_files_present)){
  print("All expected sorted alignment files present.")
  list.files(unique(dirname(sample_paths_df$..data.alignment.files)), pattern = "_sorted.bam$", full.names = TRUE)
} else {
  paste(sample_paths_df$..data.alignment.files[!expected_files_present], "not present", sep = " ")
  paste(sample_paths_df$..data.alignment.files[expected_files_present], "created", sep = " ")
}

