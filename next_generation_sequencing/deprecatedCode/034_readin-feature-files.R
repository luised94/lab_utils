#STATUS:

#Gets relative paths for the directories in data folder and assigns them to variables
if(grepl("Windows", osVersion)){
  source("./scripts/assign-directory-variables.R")
} else if (grepl("Linux", osVersion)){
  source("../rscripts/assign-directory-variables.R")
}

#Load packages required by script ----
package_to_check <- c("tidyverse","tools","genomation","ChIPpeakAnno","GenomicFeatures","ShortRead")
loadPackages(package_to_check)

# Loading feature files in working directory ----

#Get folder that contains feature files. Identified by having Feature in name
feature_folder <- list.dirs(recursive = FALSE)[grepl("Features", list.dirs(recursive = FALSE))]
if(grepl("Windows", osVersion)){
  feature_folder <- list.dirs(recursive = FALSE)[grepl("Features", list.dirs(recursive = FALSE))]
} else if (grepl("Linux", osVersion)){
  feature_folder <- '../feature-files'
}
#/home/luised94/data/feature-files or ../feature-files when on cluster 

#Define files that I want to load 
annofile_extension <- c("*.bed", "*.gff3", "*.xls", "*.tab")

#Create list, find files in feature_folder that contain annotations of interest
feature_files <- c()
for (extension in 1:length(annofile_extension)){
  feature_files<- c(feature_files, 
                    list.files(feature_folder, pattern = annofile_extension[extension]))
}

#Check file extensions
unique(tools::file_ext(feature_files))

#Create variable for genome to normalize chromosome name ----
genome_file_path <- paste(genome_directory[1], "sacCer3.fasta", sep = "/")
sacCer3 <- readFasta(genome_file_path)
sacCer3_df <- data.frame(chrom = names(as(sacCer3, "DNAStringSet")), size = width(sacCer3)) %>% filter(chrom != "chrM")
rm(sacCer3)



#For all of the dataframes to be created, print the name of the file to be imported ----
#Names for the data frames to be created 
# df_names <- c("ChExMix","MEME","XUTs","CUTs","ORF","Nucleosome",
#               "timing", "Rhee" , "SGD_features", "SUTs", "G2_orc")

for(i in 1:length(df_names)){
  #Print file used
  print(feature_files[grepl(df_names[i], feature_files)])
  #Create filepath to the file by identifying filename that contains name to be used for dataframe
  #Use grepl to determine if df_name is in feature_files and index feature_files with it
  filepath <- file.path(feature_folder, feature_files[grepl(df_names[i], feature_files)])
  
  
  #if structure to load in file depending on file extension
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
#Construct a GRange object from timing dataframe
#GRanges(seqnames = timing$Chromosome, ranges = IRanges(start = timing$Position-100, end = timing$Position+100), strand = "*", timing$T1/2)
