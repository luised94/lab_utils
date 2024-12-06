#STATUS: REMOVE.
#Get time at start of script ----
start_time <- Sys.time()
print("Script start")

#Prepares directory variables and the dataframe that contains paths, connections and variables to the fastq samples 
if(grepl("Windows", osVersion)){
  source("./scripts/prepare-samples-dataframe.R")
} else if (grepl("Linux", osVersion)){
  source("../rscripts/prepare-samples-dataframe.R")
}

#Load the packages required by the script ----
package_to_check <- c("QuasR", "GenomicAlignments","Gviz","rtracklayer","ShortRead")
loadPackages(package_to_check)

#Create variable for genome ----
genome_file_path <- paste(genome_directory[1], "sacCer3.fasta", sep = "/")
sacCer3 <- readFasta(genome_file_path)
# names(as(sacCer3, "DNAStringSet"))
# width(sacCer3)
sacCer3_df <- data.frame(chrom = names(as(sacCer3, "DNAStringSet")), size = width(sacCer3)) %>% filter(chrom != "chrM")
rm(sacCer3)


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

paste(as.character(file.size(sample_paths_df$..data.bw.files[1])/10^3), "Kb file size for first file.", sep = "")

#Read in tab-delimited file with info on fastq_files 
#Read in tab-delimited text file with info for each sample. Made in C:/Users/Luis/Projects/working-on-a-cluster/scripts/fastq-bmc-info-dataframe.R -----
if(grepl("Windows", osVersion)){
  fastq_info <- read.delim("./scripts/fastq_info.txt", header = TRUE) %>% arrange(sname)
} else if (grepl("Linux", osVersion)){
  fastq_info <- read.delim("../rscripts/fastq_info.txt", header = TRUE) %>% arrange(sname)
  fastq_info <- bind_cols(fastq_info, sample_paths_df)
}

#Create GRanges object to read in a particular chromosome
chrXIV.gr <- GRanges(seqnames=c(sacCer3_df$chrom[14]), ranges = IRanges(start = 100000, end = sacCer3_df$size[14]), strand = "*")

#Create genome axis to display 
all_tracks <- list(GenomeAxisTrack())

#Assign Tracks of interest
#Can use index or file name to subset the dataframe Wnani,WnayV,nnNnH
# files_to_visualize <- c("WnNyV", "nnNnH","Wnani")
# subset_df <- fastq_info %>% filter(sname %in% files_to_visualize)
subset_df <- sample_paths_df[1:3,]
#Assign tracks using mapply and sample dataframe. May add rm to get rid of bw_variable if it is too large
#Does not return bw_variables

all_tracks <-append(all_tracks, value = mapply(function(bw_con, bw_variable, track_variable, track_name){
  assign(bw_variable,  import(con = bw_con, which = chrXIV.gr))
  assign(track_variable, DataTrack(get(bw_variable), type = "l", name = track_name, chromosome = "chrXIV"))
}, bw_con = subset_df$..data.bw.files, bw_variable = subset_df$bw_var, track_variable =  subset_df$track_var, track_name = subset_df$sample_names)
)
MAX <- -Inf
for (track in 1:length(all_tracks)) {
  if(class(all_tracks[[track]]) != "GenomeAxisTrack"){
    if(max(all_tracks[track][[1]]@data) > MAX) MAX <- max(all_tracks[track][[1]]@data)
  }
}
plotTracks(all_tracks, main = "Complete View of Chromosome 14", chromosome = "chrXIV", ylim = c(0, MAX * 1.20))
#format(object.size(get(subset_df$track_var[1]), units = "Kb"))

#Plot Tracks 
#plotTracks(list(tracks[[1]],tracks[[2]],tracks[[3]]), main = "Eaton view of Chromosome 14", chromosome = "chrXIV")

#mapply output can be used directly in plotTracks. Save to svg. 
svg(paste("./reports/plots-files/", underscoreDate(), "_Quick_comparison_of_eaton_and_V5_at_Noc", ".svg", sep = ""))
plotTracks(all_tracks, main = "Complete View of Chromosome 14", chromosome = "chrXIV")
dev.off()

print(paste("Script complete. Time Elasped (in mins)", difftime(start_time, Sys.time(), units = "mins"), sep = " "))

#Extra code and chunks to confirm how functions works ----
# #Create the bigwig coverage files for all the bam files. Uses dataframe that holds paths to all the different files
#This works for toy data
# mapply(function(bam_files, bw_connection){
#   alns <- readGAlignments(bam_files)
#   covs <- coverage(alns)
#   export.bw(con = bw_connection,
#             object = covs,
#             format = "bw")
# }, bam_files = sample_paths_df$sorted_bam, bw_connection = sample_paths_df$..data.bw.files)

#Export function from rtracklayer package overwrites by default
# export.bw(con = sample_paths_df$..data.bw.files[1],
#           object = RleList(cvg1),
#           format = "bw")


# #Establishes connection 
# BamFile(sample_paths_df$sorted_bam[1])
# #Reads in the alignment file according to param argument, pass BamFile object.
# scanBam(sample_paths_df$sorted_bam[1])

#Confirm that levels where same for bam files and sacCer3 
# levels(aln[[1]]$rname)[!grepl("chrM", levels(aln[[1]]$rname))] == sacCer3_df$chrom
