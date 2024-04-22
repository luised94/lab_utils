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






svg(paste("./reports/plots-files/", underscoreDate(), "_Quick_comparison_of_eaton_and_V5_at_Noc", ".svg", sep = ""))
plotTracks(all_tracks, main = "Complete View of Chromosome 14", chromosome = "chrXIV")
dev.off()


reate GRanges object to read in a particular chromosome
chrXIV.gr <- GRanges(seqnames=c(sacCer3_df$chrom[14]), ranges = IRanges(start = 100000, end = sacCer3_df$size[14]), strand = "*")

#Create genome axis to display
all_tracks <- list(GenomeAxisTrack())


ll_tracks <-append(all_tracks, value = mapply(function(bw_con, bw_variable, track_variable, track_name){
  assign(bw_variable,  import(con = bw_con, which = chrXIV.gr))
  assign(track_variable, DataTrack(get(bw_variable), type = "l", name = track_name, chromosome = "chrXIV"))
}, bw_con = subset_df$..data.bw.files, bw_variable = subset_df$bw_var, track_variable =  subset_df$track_var, track_name = subset_df$sample_names)
)

fastq_ids <- fastq_ids %>% mutate(sname = str_c(str_sub(complement,1,1), str_sub(suppressor,1,1),str_sub(cellcycle,1,1),str_sub(auxin,1,1),str_sub(antibody,1,1)))

#Condition that determines if sample is
is_input <- fastq_ids$antibody == 'input'

#Conditions that determine if sample is negative control
is_negative<- (fastq_ids$antibody == 'V5' & fastq_ids$complement == 'none') |
  (fastq_ids$antibody == 'Myc' & fastq_ids$auxin == 'yes') |
  (fastq_ids$antibody == 'UM174' & fastq_ids$cellcycle == 'Noc') |
  (fastq_ids$antibody == 'UM174' & fastq_ids$complement == 'none' & fastq_ids$auxin == 'yes')

#Creating a factor vector depending on the conditions.
fastq_ids$Sample_type <-  factor(case_when(is_input ~ 'input',
                 is_negative ~ 'negative',
                 TRUE ~ 'experiment'))
factor(case_when(as.numeric(rownames(fastq_ids)) <= 24 ~ 'A',
                 as.numeric(rownames(fastq_ids)) > 24 ~ 'B'))

fastq_ids <- fastq_ids %>% mutate(Pool = factor(case_when(as.numeric(rownames(fastq_ids)) <= 24 ~ 'A',
                                          as.numeric(rownames(fastq_ids)) > 24 ~ 'B')))

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

eset all_tracks variable just in case. Add GenomeAxisTrack at top
  all_tracks <- list(GenomeAxisTrack(name=paste("Chr ", chr_num, " Axis", sep = "")))

  #Assign tracks using mapply and sample dataframe in order of subset dataframe. May add rm to get rid of bw_variable if it is too large
  #Append to all tracks the result of mapply. Assigns bigwig variable subset by chromosome and creates track variable
  all_tracks <-append(all_tracks, value = mapply(function(bw_con, bw_variable, track_variable, track_name){
    # print(bw_con)
    assign(bw_variable,  import(con = bw_con, which = get(sacCer3_df$gr_var[chr_num])))
    assign(track_variable, DataTrack(get(bw_variable), type = "l", name = track_name, ylim = c(0, 2000)))
  }, bw_con = subset_df$..data.bw.files, bw_variable = subset_df$bw_var, track_variable =  subset_df$track_var, track_name = subset_df$sample_names)
  )

  #Add origin track at the end so that it shows up below the tracks.
  all_tracks <- append(all_tracks, get(sacCer3_df$origin_track_var[chr_num])


#Define files that I want to load
annofile_extension <- c("*.bed", "*.gff3", "*.xls", "*.tab")

#Create list, find files in feature_folder that contain annotations of interest
feature_files <- c()
for (extension in 1:length(annofile_extension)){
  feature_files<- c(feature_files,
                    list.files(feature_folder, pattern = annofile_extension[extension]))
}






  if (grepl("xls|xlsx", tools::file_ext(filepath))){
    assign(df_names[i], readxl::read_excel(filepath))
    if(df_names[i] == "timing"){
      assign(df_names[i], get(df_names[i]) %>% as.data.frame() %>% dplyr::select(1:7) %>% filter(!is.na(Chromosome)))
    }

  } else if(tools::file_ext(filepath) == "bed"){
    assign(df_names[i], readBed(filepath))
    print(seqlevelsStyle(eval(parse(text = df_names[i]))))

  } else if(tools::file_ext(filepath) == "gff3"){
    assign(df_names[i], toGRanges(makeTxDbFromGFF(filepath), format = 'gene'))
    print(seqlevelsStyle(eval(parse(text = df_names[i]))))

  } else if(tools::file_ext(filepath) == "tab"){
    assign(df_names[i], read.delim(filepath, header = FALSE))
    if(df_names[i] == "SGD_features"){
      assign(df_names[i], get(df_names[i]) %>% as.data.frame() %>% dplyr::select(c(9,10,11,12,2,4,5)))
    }
  }
}

#Prepare the sacCer3 dataframe with the information to assign necessary variables for reading in track/bigwig data ----
sacCer3_df <- sacCer3_df %>% mutate(gr_var = paste(chrom, "_gr", sep = ""),
                                    bw_var = paste(chrom, "_bw", sep = ""),
                                    track_var = paste(chrom, "_track", sep = ""),
                                    origin_gr_var = paste("origin_", chrom, "_track", sep = ""),
                                    chr_df= paste(chrom, "_df",sep = ""),
                                    origin_track_var = paste("origin_", chrom, "_track",sep = ""))


#Get the files with the info to repeat the plots already made.
library(stringr)
plot_input_files <- list.files(plot_folder, pattern = ".txt", full.names = TRUE)[grepl(file_id, list.files(plot_folder, pattern = ".txt"))]
plot_input <- lapply(plot_input_files, function(file_names){
  as.data.frame(read.table(file_names, header = TRUE, sep = "\t"))
})
plot_input <- bind_rows(plot_input)
plot_input <- plot_input %>% distinct(plot_names, .keep_all = TRUE)

#Process the set_of_file to turn into column of lists.
set_of_files_lists <- str_split(str_replace_all(plot_input$set_of_file, pattern = "c|\\(|\\)|,", replacement = ""), pattern = " ")
plot_input$set_of_files <- I(set_of_files_lists)

plots_to_generate <- plot_input
plots_to_generate$chromosome_for_plot <- rep(chromosome_number, length(plots_to_generate$set_of_files))

fastq_ids <- read.table("../rscripts/fastq_info.txt") #or add an extra ../ depending on current working dir. Most scripts will be set relative to 221024Bel_CHIP
knitr::kable(fastq_ids)
#Obtain files that will be renamed
fastq_files_to_rename <- list.files(path ='.', pattern = "*.fastq$", recursive = TRUE, full.names = TRUE)[!grepl("unmapped", list.files(path ='.', pattern = "*.fastq$", recursive = TRUE))]
print(fastq_files_to_rename)
#Subset the dataframe. I think it was recycling the TRUE values
bmc_ids <- fastq_ids[!grepl("nnNnH", fastq_ids$sname),]
prefix <- './data/fastq-files/'
grepl()


for (i in 1:length(fastq_files_to_rename)) {
  correct_name <- bmc_ids$sname[as.logical(na.omit(unlist(lapply(bmc_ids$BMC_ID2, grepl, x = fastq_files_to_rename[i]))))]


MAX <- -Inf
for (track in 1:length(all_tracks)) {
  if(class(all_tracks[[track]]) != "GenomeAxisTrack"){
    if(max(all_tracks[track][[1]]@data) > MAX) MAX <- max(all_tracks[track][[1]]@data)
  }
}
plotTracks(all_tracks, main = "Complete View of Chromosome 14", chromosome = "chrXIV", ylim = c(0, MAX * 1.20))

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

# BASH_SECTION
 #Rcript -e 'source("/home/luised94/data/rscripts/3-assign-directory-variables.r")'
