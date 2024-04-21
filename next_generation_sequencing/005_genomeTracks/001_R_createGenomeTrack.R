nstall.packages(c("tidyverse", "stringr","R.utils"), lib.loc = "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2")

bPaths( c( "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2", .libPaths() ) )

.libPaths()

library(tidyverse, lib.loc = "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2")
library(stringr, lib.loc = "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2")

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

_folder <- "./data"
dir_names <- stringr::str_replace(list.dirs(data_folder, recursive = FALSE), pattern = data_folder, replace = "")

# TODO Dont think I need to do all this variable assignment and I aslo included some a lot of if statements to prevent rerunning some code which got pretty crazy. Just really need to run code once. 
# Create environment directory structures for different analysis. 
# Lots of unused code in comment as well. 

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


























#Load packages using through a list of strings and suppress the messages, return a TRUE if loading was succesful
package_list <- c("QuasR", "GenomicAlignments", "Gviz", "rtracklayer", "ShortRead")
package_was_loaded <- unlist(
    suppressPackageStartupMessages(
        lapply(
            package_list, 
	    library, 
    	    character.only = TRUE, 
            logical.return=TRUE, 
	    quietly = TRUE
    )
  )
)

#LOG
#Determine which packages were not loaded, print all packages loaded or which packages were not loaded. 
packages_not_loaded <- package_list[!package_was_loaded]
if (length(packages_not_loaded) == 0) {
    print("All packages loaded.")
} else {
  lapply(
    packages_not_loaded,
    function(x) { 
      message <- paste(x, "Package did not install") 
      print(message)
    }
  )
}



#Create variable for genome ----
# TODO Search prepare-samples-dataframe.R to see how variables were initialized

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

eparate file to source for assigning variables instead of writing it each time.
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

# BASH_SECTION
Rcript -e 'source("/home/luised94/data/rscripts/3-assign-directory-variables.r")'


genome_file_pathgg- paste(genome_directory[1], "sacCer3.fasta", sep = "/")
sacCer3 <- readFasta(genome_file_path)
# names(as(sacCer3, "DNAStringSet"))
# width(sacCer3)
sacCer3_df <- data.frame(chrom = names(as(sacCer3, "DNAStringSet")), size = width(sacCer3)) %>% filter(chrom != "chrM")
rm(sacCer3)
# TODO Go through ngs-functions to extract 
#Output bigwig files by chunking through chromosomes ----
#Uses mapply to go through bam files and their bw connections. Uses function bamReadPosAndQwidthByChromToRle mostly taken from the introduction to Rsamtools manual
#https://www.bioconductor.org/packages/release/bioc/vignettes/Rsamtools/inst/doc/Rsamtools-Overview.pdf

mapply(function(bam_files, bw_connection) {
  print(paste("Working on:", bam_files, sep = " "))
  cvg <-
    RleList(
      mapply(
        FUN = bamReadPosAndQwidthByChromToRle,
        chrom_name = sacCer3_df$chrom,
        bam_file = bam_files,
        chromosome_size = sacCer3_df$size,
        SIMPLIFY = FALSE
      )
    )
  print(paste("Exporting:", bw_connection, sep = " "))
  export.bw(con = bw_connection,
            object = cvg,
            format = "bw")

}, bam_files = sample_paths_df$sorted_bam, bw_connection = sample_paths_df$..data.bw.files)

print(paste("Bigwig Files created. Time Elasped (in mins)", difftime(start_time, Sys.time(), units = "mins"), sep = " "))

# Info was read using read.delim
tq_info <- read.delim("../rscripts/fastq_info.txt", header = TRUE) %>% arrange(sname)
}




utputfiles=(WT-G2-ORC-rep1.fastq.gz WT-G2-ORC-rep2.fastq.gz)
websites=(ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR034/SRR034475/SRR034475.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR034/SRR034476/SRR034476.fastq.gz)

COUNTER=0
for file in ${outputfiles[@]};
do
  echo $file
  echo $COUNTER
  echo ${websites[COUNTER]}
  wget --output-document=$file ${websites[COUNTER]}
  cat $file >> nnNnH.fastq.gz
  rm $file
  let COUNTER++
done

# #Check that file is right size
# ls -l --block-size=M
#Do this after files have been renamed
gunzip nnNnH.fastq.gz
#cat 221024Bel_CHIP.txt | sed s/"    "/""/g | sed s/":"/""/g | sed s/"  "/" "/g | sed s/" "/"\t"/g










#### BASH END
