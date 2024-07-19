#Run with: ----
#$sbatch sbatch-generate-comparison-plots.sh 7 
#See ./regenerate-comparison-coverage-plots.R to get parameters used for plot. 

#Get time at start of script ----
start_date <- Sys.time()
start_time <- as.character(stringr::str_replace_all(start_date, pattern = ":| |-", replacement="_"))
print(paste("Script start:", start_date, sep = " "))
#Handle script arguments ----
if(grepl("Windows", osVersion)){
  chromosome_number <- 7
} else if (grepl("Linux", osVersion)){
  chromosome_number <- as.numeric(args[1])
  args = commandArgs(trailingOnly=TRUE)
  if (length(args)==0) {
    stop("Pass the chromosome number to batch script", call.=FALSE)
  } 
  print(paste("Args used in script are:", args[1], sep = " "))
}

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
#Variables to be visualized with track data. #Also creates the sacCer3_df variable with the variables for assigining track data.
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
columns_to_exclude <- c("BMC", "Index", "V2", "V3","sname", "Sample_names", "..data","original","track","fq", "bw","track","sorted")
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

#Use the dataframe to assign the tracks for the chromosome GRanges and tracks ----
#See ./ngs-functions.R for assignChrTracks details.

assignChrTracks(sacCer3_df, chromosome_number, timing, 200)

#Create the dataframe with the different comparisons to make and the name to assign to the plot ----
#Dataframe has two columns, one with the files that will be accessed and the other with the name of the plot.
#Use mapply to iterate through the dataframe. 
#Files that will be used to subset dataframe
#Removed nnNyU from effect of Auxin on MCM because it was giving trouble and wanted to check other samples. 
files_to_visualize = list(c("nnani","R1ani","R4ani","Rnani","RTani","Wnani"))
#List of names for the plots that will be generated.
name_of_plots = c("Input_Samples")
#Variable to construct the track name

#Chromosome visualized in plots
chromosome_for_plot = rep(chromosome_number, length(files_to_visualize))
#Create the dataframe. Use I() function to treat list of lists as is. 
plots_to_generate <- data.frame(set_of_files = I(files_to_visualize), plot_names = name_of_plots,
                                plot_chromosome = chromosome_for_plot)

#Setting the output folder
if(grepl("Windows", osVersion)){
  output_folder = "./reports-and-figures/plots-files/"
} else if (grepl("Linux", osVersion)){
  output_folder = "../reports/plots-files/"
}

#Create the plots ----
name_of_files <- c()
mapply(function(files, names, chr_num){
  # print(files);print(names)
  # print(str_replace_all(names, pattern = "_", replacement = " "))
  #Subset data frame according to files and arrange by antibody to display according to that order in Gviz
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
    assign(track_variable, DataTrack(get(bw_variable), type = "l", name = track_name, ylim = c(0, 2000)))
  }, bw_con = subset_df$..data.bw.files, bw_variable = subset_df$bw_var, track_variable =  subset_df$track_var, track_name = subset_df$sample_names)
  )

  #Add origin track at the end so that it shows up below the tracks.
  all_tracks <- append(all_tracks, get(sacCer3_df$origin_track_var[chr_num]))

  #Create the plot.
  # print(paste("Creating ", paste(output_folder, paste(start_time,"_",sep = ""), "chr", chr_num, "_", names, ".svg", sep = ""), sep = ""))
  svg(paste(output_folder, paste(start_time,"_",sep = ""), "chr", chr_num, "_", names, ".svg", sep = ""))
  plotTracks(all_tracks, main = paste(str_replace_all(names, pattern = "_",replacement = " "),  sep = " "))
  dev.off()
  print(paste(paste(output_folder, paste(start_time,"_",sep = ""), "chr", chr_num, "_", names, ".svg", sep = ""), "Plot created.", sep = " "))
  
  name_of_files <<- append(name_of_files, paste(output_folder, paste(start_time,"_",sep = ""), "chr", chr_num, "_", names, ".svg", sep = ""))
}, files = plots_to_generate$set_of_files, names = plots_to_generate$plot_names, chr_num = chromosome_number)
#Output a table with the values for each variable in the 221024Bel_CHIP experiment. ----
#Add name_of_files, chr name and start_time to plots_to_generate dataframe
plots_to_generate$file_name <- name_of_files
plots_to_generate$plot_time <- rep(start_time, length(files_to_visualize))

plot_input_file <- paste(output_folder, paste(start_time,"_",sep = ""), "chr", chromosome_number, "_", "plot_input.txt", sep = "")
write.table(plots_to_generate, file = plot_input_file, row.names = FALSE, quote = FALSE, sep = "\t")

print(paste("Plots created. Time Elasped (in mins)", difftime(start_date, Sys.time(), units = "mins"), sep = " "))
sessionInfo()

#Output a table with the unique values of the treatment column for Steve ----
#Install kableExtra and run webshot::install_phamtomjs()
# knitr::kable(data.frame(
#   Variables = treatments,
#   Values = str_replace_all(
#     as.character(unique_values),
#     pattern = "c|\\(|\\)|,|\"",
#     replacement = ""
#   )
# ),
# format = "html") %>% kable_classic(full_width = F, position = "center") %>% save_kable(file = "table_with_variable_values_for_221024Bel_CHIP.png", zoom = 1.5)

# 
# html_table_width <- function(kable_output, width){
#   width_html <- paste0(paste0('<col width="', width, '">'), collapse = "\n")
#   sub("<table>", paste0("<table>\n", width_html), kable_output)
# }
# fastq_info[, sapply(fastq_info, class) == 'character']
# fastq_info %>% select_if(colnames(.) %in% treatments) %>% apply(., 2, table)
# fastq_info$sample_names == fastq_info$sname
# levels(fastq_info$antibody)
# unique_values <- lapply(character_to_factor, FUN = function(column) {
# unique(fastq_info[, names(fastq_info)[grepl(column, names(fastq_info))]])
# })
# fastq_info  %>% mutate_at(character_to_factor, as.factor) 
