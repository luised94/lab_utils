#Exporting and visualizing coverage data 
library(tidyr)
library(dplyr)
library(ChIPpeakAnno)
library(tools)
library(genomation)
library(normr)
library(rtracklayer)
library(QuasR)
library(seqinr)
library(Gviz)
library(Rsamtools)
library(GenomicAlignments)
#Getting directories 
processed_fastq_folder <- list.dirs(path = "./data")[grepl("processed", list.dirs(path = "./data"))]
alignment_target_directory <- list.dirs(path = "./data")[grepl("alignment", list.dirs(path = "./data"))]
genome_directory <- list.dirs(path = "./data")[grepl("genome-files$", list.dirs(path = "./data"))]
peak_directory <- list.dirs(path = "./data")[grepl("peak-files$", list.dirs(path = "./data"))]
bw_directory<- list.dirs(path = "./data")[grepl("bw", list.dirs(path = "./data"))]

#Reading in bam_summary file.
bam_summary <- read.table(paste(alignment_target_directory, "bam-summary.txt", sep = "/"))
#Create column with location 
bam_summary <- bam_summary %>% mutate(sorted_file_location = paste(row.names(bam_summary), "_sorted.bam", sep="")) 
bam_summary <- bam_summary %>% mutate(bw_location = paste(bw_directory, stringr::str_replace(basename(sorted_file_location), ".bam", ".bw"), sep = "/"))

#Generate bigwig files for all bam files 
for (i in length(bam_summary$sorted_file_location)){
  
}


lapply(bam_summary$sorted_file_location, function(x){
 alns <- readGAlignments(x)
 covs <- coverage(alns)
 # paste(bw_directory, stringr::str_replace(basename(x), ".bam", ".bw"), sep = "/")
 export.bw(con = paste(bw_directory, stringr::str_replace(basename(x), ".bam", ".bw"), sep = "/"),
                     object = covs,
                     format = "bw")
})


#Generate bigwig file 
# param <- ScanBamParam()
# alns <- readGAlignments(input_control)
# alns[seqlevels(alns) == "chrVII"]
# covs <- coverage(alns, weight = ) #or coverage(input_control)
# covs$chrVII
# slice(covs$chrVII)

#Export to bigwig 
# export.bw(con = paste(bw_directory, stringr::str_replace(basename(input_control), ".bam", ".bw"), sep = "/"),
#           object = covs, 
#           format = "bw") 
