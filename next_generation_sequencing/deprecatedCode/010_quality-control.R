#STATUS: REMOVE.
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

package_to_check <- c("QuasR", "GenomicAlignments","Gviz","rtracklayer","ShortRead","Rsamtools","Rqc")
loadPackages(package_to_check) #Load packages and output message when done or if any package didnt load. Takes list of characters.

#Read in feature file to get sacCer3_df then delete timing variables
df_names <- c("timing")

if(grepl("Windows", osVersion)){
  source("./scripts/readin-feature-files.R")
} else if (grepl("Linux", osVersion)){
  source("../rscripts/readin-feature-files.R")
}
rm("timing")

#Read in tab-delimited text file with info for each sample. Made in ./scripts/fastq-bmc-info-dataframe.R -----
if(grepl("Windows", osVersion)){
  fastq_info <- read.delim("./scripts/fastq_info.txt", header = TRUE) %>% arrange(sname)
} else if (grepl("Linux", osVersion)){
  fastq_info <- read.delim("../rscripts/fastq_info.txt", header = TRUE) %>% arrange(sname)
}
if (grepl("Linux", osVersion)) {
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

#Use Rsamtools package to count mapped, unmapped and secondary alignments ----
unmappedparam <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = TRUE))
mappedparam <- ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE))
secondaryAlignmentparam <-ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = TRUE))
totalparam <- ScanBamParam()
paramlist <- c(unmappedparam,mappedparam,secondaryAlignmentparam,totalparam)
# scanBamWhat()[c(5:7,13)]

bam_stats <- as.data.frame(t(as.data.frame(lapply(sample_paths_df$sorted_bam, function(bam_file){
  unlist(lapply(paramlist, function(param){
    countBam(bam_file, param = param)$records
  }))
}))), row.names = FALSE) %>% rename("V1" = "Unmapped", "V2" = "Mapped", "V3" = "Secondary", "V4" = "Total") %>% mutate(across(1:2, ~ .x/Total, .names = "percent.{.col}"))

#Use ShortRead to get statistics for fastq files ----
countFastq(sample_paths_df$..data.processed.fastq.files)$records

stream <- open(FastqStreamer(sample_paths_df$..data.fastq.files[1]))
on.exit(close(stream))
tmp_fq <- yield(stream)

#If input is 0, stop the function
if (length(tmp_fq) == 0){
  print("Processed all")
  break
}

length(tmp_fq)
length(clean(tmp_fq))
length(tmp_fq[!srduplicated(tmp_fq)])
rowSums(as(quality(tmp_fq), "matrix"))/width(tmp_fq)

width(trimEnds(tmp_fq, a ="5"))
##Calculate meanq_reads to obtain reads with mean quality larger than 20
meanq_reads <- rowSums(as(quality(tmp_fq), "matrix"), na.rm = TRUE)/width(tmp_fq)
#meanq_reads[is.na(meanq_reads)] <- 0]
#Had to include a filter for NA values. Not sure why some of them have it. 
tmp_fq <- tmp_fq[meanq_reads[!is.na(meanq_reads)] >= 20]


repeat {
  ## input chunk
  tmp_fq <- yield(stream)
  
  #If input is 0, stop the function
  if (length(tmp_fq) == 0){
    print("Processed all")
    break
  }
  ##Trim reads that contain 'N'
  tmp_fq <- clean(tmp_fq)
  ##Filter duplicated reads 
  tmp_fq <- tmp_fq[!srduplicated(tmp_fq)]
  ##Trim reads
  tmp_fq <- trimEnds(tmp_fq, a ="5")
  ##Calculate meanq_reads to obtain reads with mean quality larger than 20
  meanq_reads <- rowSums(as(quality(tmp_fq), "matrix"), na.rm = TRUE)/width(tmp_fq)
  #meanq_reads[is.na(meanq_reads)] <- 0]
  #Had to include a filter for NA values. Not sure why some of them have it. 
  tmp_fq <- tmp_fq[meanq_reads[!is.na(meanq_reads)] >= 20]
  
  ## append to destination
  ##Do not compress, output quality as FastqQuality 
  # writeFastq(tmp_fq, destination_file, "a", compress = FALSE , qualityType = "FastqQuality")
} 

closeAllConnections()
print("Script ended. Time Elapsed:", difftime(start_date, Sys.time(), units = "mins"))
  #Run on cluster command ----
  #$sbatch sbatch-quality-control.sh
  # Attempt to use rqc package failed too much. Thought it was due to memory issues but unclear. Wed Nov 16 19:23:37 2022 ------------------------------
  # See slurm-7275251.out to see last attempt.   
  #Run rqc on processed fastq files by several groups 
  #Use lapply and mapply to run rqcQA on each sample with group info, combine results and output report for each group
  
  #Used commented lines to test on my toy data ----
  # sample_paths_df$pools <- unlist(c(rep("A", 5), rep("B", 4)))
  # sample_paths_df$groups <- unlist(c(rep("potato", 4), rep("corn", 5)))
  # groups_for_rqc <- c("pools","groups")
  # 
  # #Strings that contain part of column to subset by 
  # groups_for_rqc <- c("Sample","Pool")
  # 
  # #For each of the columns in groups_for_rqc, grab the column name it specifies
  # column_names_for_rqc <- unlist(lapply(groups_for_rqc, function(x) {
  #   names(fastq_info)[grepl(x, names(fastq_info))]
  # }))
  # 
  # #Apply rqcQA function to each sample and its group in processed.fastq.files for each group that we want. Use mapply inside lapply
  # lapply(column_names_for_rqc, function(column_name, fastq_files, dataframe_with_groups){
  #   #Subset dataframe using column name
  #   group_info <- dataframe_with_groups[, column_name]
  #   print(paste("Rqc by", column_name, sep = " "))
  #   
  #   #rqcQA on each fastq file and its group value. To get group value, subset dataframe with column name 
  #   qcresults <- mapply(function(fastq_files_for_rqc, rqc_groups) {
  #     rqcQA(x = fastq_files_for_rqc, sample = FALSE, group = rqc_groups, workers = 1)
  #   }, fastq_files_for_rqc = fastq_files, rqc_groups = group_info, SIMPLIFY = FALSE)
  #   
  #   rqcReport(qcresults, outdir = './reports/reports-files', file = paste("rqc_report_by_",column_name, sep = ""))
  #   
  # }, fastq_files = sample_paths_df$..data.processed.fastq.files, dataframe_with_groups = fastq_info)
  
  
  # %>% mutate(total = rowSums(across(1:3), na.rm = T))
