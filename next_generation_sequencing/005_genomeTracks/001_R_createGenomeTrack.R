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

genome_file_path <- paste(genome_directory[1], "sacCer3.fasta", sep = "/")
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

fastq_info <- bind_cols(fastq_info, sample_paths_df)

reate GRanges object to read in a particular chromosome
chrXIV.gr <- GRanges(seqnames=c(sacCer3_df$chrom[14]), ranges = IRanges(start = 100000, end = sacCer3_df$size[14]), strand = "*")

#Create genome axis to display
all_tracks <- list(GenomeAxisTrack())

#Assign Tracks of interest
#Can use index or file name to subset the dataframe Wnani,WnayV,nnNnH
files_to_visualize <- c("WnNyV", "nnNnH","Wnani")
subset_df <- fastq_info %>% filter(sname %in% files_to_visualize)



#Assign tracks using mapply and sample dataframe. May add rm to get rid of bw_variable if it is too large
#Does not return bw_variables
all_tracks <-append(all_tracks, value = mapply(function(bw_con, bw_variable, track_variable, track_name){
  assign(bw_variable,  import(con = bw_con, which = chrXIV.gr))
  assign(track_variable, DataTrack(get(bw_variable), type = "l", name = track_name, chromosome = "chrXIV"))
}, bw_con = subset_df$..data.bw.files, bw_variable = subset_df$bw_var, track_variable =  subset_df$track_var, track_name = subset_df$sample_names)
)

Extra code and chunks to confirm how functions works ----
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

























