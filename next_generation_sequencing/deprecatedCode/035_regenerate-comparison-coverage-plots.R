#STATUS: REMOVE.
#Run with: ----
#Login to luria via cmd line using ssh  -A -Y account@luria.mit.edu using your username and password.
#Add R module using module add r/4.2.0
#$sbatch sbatch-regenerate-comparison-plots.sh 7 "chr7_plot_input"
#See ./regenerate-comparison-coverage-plots.R to get parameters used for plot. 

#Get time at start of script ----
start_date <- Sys.time()
start_time <- as.character(stringr::str_replace_all(start_date, pattern = ":| |-", replacement="_"))
print(paste("Script start:", start_date, sep = " "))

#Prepares directory variables and the dataframe that contains paths, connections and variables to the fastq samples ----
#This script sources assign-directory-variables.R, which in turn also sources ngs-functions.R
if(grepl("Windows", osVersion)){
  source("./scripts/prepare-samples-dataframe.R")
} else if (grepl("Linux", osVersion)){
  source("../rscripts/prepare-samples-dataframe.R")
}

package_to_check <- c("QuasR", "GenomicAlignments","Gviz","rtracklayer","ShortRead")
loadPackages(package_to_check) #Load packages and output message when done or if any package didnt load. Takes list of characters.

#Read in feature files according to their file type. ----
#Variables to be visualized with track data. 
#Also creates the sacCer3_df variable with the variables for assigining track data.
#Set df_names variables with any of the below options.
# df_names <- c("ChExMix","MEME","XUTs","CUTs","ORF","Nucleosome",
#               "timing", "Rhee" , "SGD_features", "SUTs", "G2_orc")
df_names <- c("timing")

if(grepl("Windows", osVersion)){
  source("./scripts/readin-feature-files.R")
} else if (grepl("Linux", osVersion)){
  source("../rscripts/readin-feature-files.R")
}

#Read in tab-delimited text file with info for each sample. Made in ./scripts/fastq-bmc-info-dataframe.R -----
if(grepl("Windows", osVersion)){
  fastq_info <- read.delim("./scripts/fastq_info.txt", header = TRUE) %>% arrange(sname)
} else if (grepl("Linux", osVersion)){
  fastq_info <- read.delim("../rscripts/fastq_info.txt", header = TRUE) %>% arrange(sname)
  fastq_info <- bind_cols(fastq_info, sample_paths_df)
}


#Turn character columns into factors to reorder according to the levels variable.
#Also use for plotting on Gviz. 
#Columns to leave as characters. Then select the columns to be turned into factors
columns_to_exclude <- c("BMC", "Index", "V2", "V3","sname", "Sample_names", "..data","original","track","fq", "bw","sorted")
character_to_factor <- fastq_info %>% select_if(is.character) %>% dplyr::select(-contains(columns_to_exclude)) %>% colnames(.)

#Use to reorder the columns. Does not generalize!!! 
order_of_columns <- list(c(1,3,2), c(3,1,2,4), c(1:2),c(1:2),c(2,1,3,5,4),c(3,4,2,1),c(1:3))

invisible(mapply(function(char_to_factor, order_for_columns){
  #Get the reordered levels. Determined by watching the original order. 
  reorder_levels <- levels(factor(as.factor(fastq_info[, char_to_factor])))[order_for_columns]
  fastq_info[, char_to_factor] <<- factor(as.factor(fastq_info[, char_to_factor]), levels = reorder_levels)
}, char_to_factor = character_to_factor, order_for_columns = order_of_columns, SIMPLIFY = FALSE))

if(fastq_info %>% select_if(is.factor) %>% colnames(.) %>% length() == length(character_to_factor)){
  print("Columns refactored")
}

#Handle script arguments ---- 
#If running in windows, assign chromosome number and file id manually. 
#If running on linux system ( i. e. the Luria cluster) assign based on arguments passed to bash script
if(grepl("Windows", osVersion)){
  chromosome_number <- 7
  #String that will be used to locate the files with inputs for plots to repeat.
  #String has to be found in any of the txt files in plot_folder. Example is "chr7_plot_input"
  file_id <- args[2]
  
} else if (grepl("Linux", osVersion)){
  args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0) {
    stop("Pass chromosome number to script as a number.", call.=FALSE)
  } else if (length(args)==1) {
    # default output file
    stop("Pass string to find a plot text file to repeat.", call.=FALSE)
  }
  print(paste("Args used in script are:", args[1], args[2], sep = " "))
  chromosome_number <- as.numeric(args[1])
  
  #String that will be used to locate the files with inputs for plots to repeat.
  #String has to be found in any of the txt files in plot_folder. Example is "chr7_plot_input"
  file_id <- args[2]
}

#Use the dataframe to assign the tracks for the chromosome GRanges and tracks ----
#See ./ngs-functions.R for assignChrTracks details.

assignChrTracks(sacCer3_df, chromosome_number, timing, 200)

#Get the info to redo the plots ---- 
#If you want to read the plot file back in, it has to be processed first and then it should be equivalent. 
if(grepl("Windows", osVersion)){
  plot_folder <- './reports-and-figures/plots-files'
} else if (grepl("Linux", osVersion)){
  plot_folder <- './reports/plots-files'
}


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

#Create the plots ----
name_of_files <- c()
mapply(function(files, names, chr_num){
  # print(files);print(names)
  # print(str_replace_all(names, pattern = "_", replacement = " "))
  #Subset data frame according to files and arrange by sample columns to display according to that order in Gviz
  subset_df <- fastq_info %>% filter(sname %in% files) %>% arrange(antibody, complement, suppressor, cellcycle, auxin)
  
  #check that subset occured correctly. File length should be same as subset row number.
  if(length(subset_df$sname) != length(files)){
    print("Subsetting didnt occur properly.")
    break
  }
  
  #Reset all_tracks variable just in case. Add GenomeAxisTrack at top
  all_tracks <- list(GenomeAxisTrack(name=paste("Chr ", chr_num, " Axis", sep = "")))
  
  #Assign tracks using mapply and sample dataframe in order of subset dataframe. May add rm to get rid of bw_variable if it is too large
  #Append to all tracks the result of mapply. Assigns bigwig variable subset by chromosome and creates track variable
  all_tracks <-append(all_tracks, value = mapply(function(bw_con, bw_variable, track_variable, track_name){
    # print(bw_con)
    assign(bw_variable,  import(con = bw_con, which = get(sacCer3_df$gr_var[chr_num])))
    assign(track_variable, DataTrack(get(bw_variable), type = "l", name = track_name))
  }, bw_con = subset_df$..data.bw.files, bw_variable = subset_df$bw_var, track_variable =  subset_df$track_var, track_name = subset_df$sample_names)
  )
  MAX <- -Inf
  for (track in 1:length(all_tracks)) {
    if(class(all_tracks[[track]]) != "GenomeAxisTrack"){
      if(max(all_tracks[track][[1]]@data) > MAX) MAX <- max(all_tracks[track][[1]]@data)
    }
  }
  #Add origin track at the end so that it shows up below the tracks.
  all_tracks <- append(all_tracks, get(sacCer3_df$origin_track_var[chr_num]))
  
  #Create the plot.
  # print(paste("Creating ", paste("./reports/plots-files/", paste(start_time,"_",sep = ""), "chr", chr_num, "_", names, ".svg", sep = ""), sep = ""))
  svg(paste("./reports/plots-files/", paste(start_time,"_",sep = ""), "chr", chr_num, "_", names, ".svg", sep = ""))
  plotTracks(all_tracks, main = paste(str_replace_all(names, pattern = "_",replacement = " "),  sep = " "), ylim = c(0, MAX * 1.20))
  dev.off()
  print(paste(paste("./reports/plots-files/", paste(start_time,"_",sep = ""), "chr", chr_num, "_", names, ".svg", sep = ""), "Plot created.", sep = " "))
  name_of_files <<- append(name_of_files, paste("./reports/plots-files/", paste(start_time,"_",sep = ""), "chr", chr_num, "_", names, ".svg", sep = ""))
}, files = plots_to_generate$set_of_files, names = plots_to_generate$plot_names, chr_num = chromosome_number)

#Output a table with the values for each variable in the 221024Bel_CHIP experiment ---- 
#Add name_of_files, chr name and start_time to plots_to_generate dataframe
plots_to_generate$file_name <- name_of_files
plots_to_generate$plot_time <- rep(start_time, length(plots_to_generate$set_of_files))

plot_input_file <- paste("./reports/plots-files/", paste(start_time,"_",sep = ""), "chr", chromosome_number, "_", "plot_input.txt", sep = "")
write.table(plots_to_generate, file = plot_input_file, row.names = FALSE, quote = FALSE, sep = "\t")

print(paste("Plots created. Time Elasped (in mins)", difftime(start_date, Sys.time(), units = "mins"), sep = " "))
sessionInfo()
