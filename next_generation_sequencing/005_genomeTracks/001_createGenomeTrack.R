#DESCRIPTION: 
#USAGE:
##TODO: Define the comparisons and plots to be generated in my sampleConfig.R template. 
##TODO: Find the best way to have the same levels and factors when I read in the sampleGridConfig.R file. This is defined by the categories list variable. Can grab that during 002_loadSampleGrid and use it to equalize. 
##TODO: Use 002_loadSampleGrid to open the sample_table given a directory. Use conditional statement to determine the ID. Could potentially use system function to call find and print with awk statement. R would require list.files(), strsplit, and grabbing regular expression for 5 digits.
#Load packages using through a list of strings and suppress the messages, return a TRUE if loading was succesful
##TODO: Need to make sure the chromosome IDs are formatted properly.
##TODO: need to automate the creation of the experiments to plot. 
##TODO: Need to ensure that the names I use in the columns are compatible with my short name convention for subsetting df
main <- function() {
    suppressPackageStartupMessage({
        library(QuasR)
        library(GenomicAlignments)
        library(Gviz)
        library(rtracklayer)
        library(ShortRead)
        library(tidyverse)
    })
    args <- commandArgs(trailingOnly = TRUE)
    directory_path <- validate_input(args)
    sample_table <- load_sample_table(directory_path)
    refGenome <- load_reference_genome()
    # Has to be run every time if you are using chromosome ID to get tracks from the bigwig file.
    #options(ucscChromosomeNames=FALSE)
}
#validate_input(args) {
#}
#load_sample_table <- function(directory_name) {
#    cat("Loading sample_table from", directory_name, "\n")
#    documentation_dir_path <- file.path(directory_name, "documentation")
#    sample_table_path <- list.files(documentation_dir_path, pattern = "sample_table", full.names = TRUE)
#    if(length(sample_table_path) == 0){
#        cat(sprintf("No files with pattern sample_table found in %s\n", documentation_dir_path))
#        cat("Consult 000_setupExperimentDir to create sample_table.tsv for the project\n")
#        stop()
#    } else if(length(sample_table_path) > 1){
#        cat(sprintf("Multiple files with pattern sample_table found in %s\n", documentation_dir_path))
#        cat(sprintf("Files found in %s\n", documentation_dir_path))
#        print(sample_table)
#        cat("Consult 000_setupExperimentDir and ensure no duplicates are present for the project\n")
#        stop()
#    }
#    sample_table <- read.delim(sample_table_path, header = TRUE, sep ="\t")
#    if (!("sample_ID" %in% colnames(sample_table)) {
#           cat("Sample table does not contain the sample_ID column.\n")
#           cat("sample_ID column is required for analysis.\n")
#           cat("See 000_setupExperimentDir.R and 002_verifySampleGrid.R.\n")
#           stop()
#    }
#    cat(sprintf("Reading %s\n", sample_table_path))
#    cat("Head of sample_table\n")
#    print(head(sample_table))
#    return(sample_table)
#}
#load_reference_genome(genome_dir = "REFGENS, genome_pattern = "S288C_refgenome.fna") {
#   directory_of_refgenomes <- file.path(Sys.getenv("HOME"), "data", genome_dir)
#   if(!dir.exists(directory_of_refgenomes)) {
#       stop("Directory with reference genomes doesnt exist.\n")
#   }
#   genome_file_path <- list.files(directory_of_refgenomes, pattern = genome_pattern, full.names = TRUE)
#   if(!file.exists(genome_file_path)) {
#       stop("Reference genome doesnt exist.\n")
#   }
#   refGenome <- readFasta(genome_file_path)
#   refGenome <- data.frame(chrom = names(as(df_sacCer_refGenome, "DNAStringSet")), 
#                   basePairSize = width(df_sacCer_refGenome)) %>% filter(chrom != "chrM")
#   return(refGenome)    
#}
#
#create_chromosome_GRange <- function(refGenome, chromosome_to_plot){
##Create GRanges object to read in a particular chromosome
#chromosome_to_plot <- 10
#genomeRange_to_get <- GRanges(seqnames=c(df_sacCer_refGenome$chrom[chromosome_to_plot]), 
#        ranges = IRanges(start = 1, 
#        end = df_sacCer_refGenome$basePairSize[chromosome_to_plot]), 
#        strand = "*")
#}
#load_feature_file_GRange <- function(chromosome_to_plot, feature_file_pattern = "eaton_acs")
## Create list with samples to plot. 
#experiments_to_plot <- list(
#    c("lnnnMI", "onnnMA", "lnnnMA", "lnnyMA"),
#    c("onnnMI", "onnnG7", "on7nG7", "onnnM7", "on7nM7")
#)
#
#descriptive_names_for_plots <- c(
#    "AlFA_comparison",
#    "MCM7_MCMpoly"
#)
##TODO: Need to implement in sampleGrid and experiment design.
##experiments_to_plot <- list(
##    c("lnnnMI", "onnnMA", "lnnnMA", "lnnyMA"),
##    c("lnnnMI", "onnnMH", "lnnnMH", "lnnyMH"), 
##    c("lnnnMI", "lwnnG7", "lwnnM7"),
##    c("onnnMI", "onnnG7", "on7nG7", "onnnM7", "on7nM7"), 
##    c("onnnMI", "onnnG1", "on7nG1", "onnnM1", "on7nM1"), 
##    c("onnnMI", "onnnGC", "on7nGC", "onnnMC", "on7nMC"), 
##    c("lnnnMI", "lwnnG1", "lw2nG1", "lw2nM1"), 
##    c("lnnnMI", "lwnnGC", "lw2nGC", "lw2nMC"), 
##    c("lnnnMI", "lwnnG7", "lw2nG7", "lw2nM7")
##) 
#
##descriptive_names_for_plots <- c(
##    "AlFA_comparison",
##    "ORC_comparison", 
##    "Cell_Cycle_MCMpoly",
##    "MCM7tag_MCMpoly", 
##    "MCM7tag_HA11", 
##    "MCM7tag_HAC", 
##    "MCM2tag_HA11", 
##    "MCM2tag_HAC", 
##    "MCM2tag_MCMpoly"
##)
#
#
##Create variables to name plot
#main_title_of_plot_track <- paste("Complete View of Chrom", as.character(chromosome_to_plot, 
#                  sep = " "))
#date_plot_created <- stringr::str_replace_all(Sys.time(), pattern = ":| |-", replacement="")  
#
#bigwig_directory <- paste(working_directory, "bigwig", sep = "/")
#plot_output_dir <- paste(working_directory, "plots", sep = "/")
#for (experiment_index in 1:length(experiments_to_plot)) {
#    all_tracks_to_plot <- list(GenomeAxisTrack(name = paste("Chr ", chromosome_to_plot, " Axis", sep = "") ) )
#    #Subset the dataframe using strings    
#    #TODO: Need some way to aggregate the column names I use from the sampleConfig.R files to keep track of variables I use to keep consistent experiment to experiment.
#    #Files are already ordered from the sampleGridConfig.R file. No need to arrange.
#    df_sample_info_subset <- df_sample_info %>% filter(short_name %in% unlist(experiments_to_plot[experiment_index])) #%>% arrange(antibody, rescue_allele, mcm_tag, cell_cycle, auxin_treatment)
#    #print("Head output of subset dataframe")
#    #print(head(df_sample_info_subset))
#    #Create all of the tracks and append them to the all_tracks_to_plot variable
#    for (sample_index in 1:nrow(df_sample_info_subset)) {
#        sample_ID_pattern <- df_sample_info_subset$sample_ID[df_sample_info_subset$short_name == unlist(experiments_to_plot[experiment_index])[sample_index]]
#        initial_matches <- list.files(bigwig_directory, pattern = as.character(sample_ID_pattern), full.names = TRUE, recursive = TRUE)     
#        path_to_bigwig <- initial_matches[grepl("S288C", initial_matches)]
#        bigwig_to_plot <- import(con = path_to_bigwig, which = genomeRange_to_get)
#        sample_short_name <- df_sample_info_subset$short_name[df_sample_info_subset$short_name == unlist(experiments_to_plot[experiment_index])[sample_index]]
#        track_to_plot <- DataTrack(bigwig_to_plot, type = "l", name = sample_short_name, chromosome = df_sacCer_refGenome$chrom[chromosome_to_plot])
#        all_tracks_to_plot <- append(all_tracks_to_plot, track_to_plot)                                                                          
#
#        cat(sprintf("Experiment_to_plot sample: %s", unlist(experiments_to_plot[experiment_index])[sample_index]), "\n")
#        cat(sprintf("Sample Short_name: %s", sample_short_name), "\n")
#        cat(sprintf("Sample ID_pattern: %s", sample_ID_pattern), "\n")
#        cat("Using experiments_to_plot subset","\n")
#        cat(sprintf("Bigwig file to access: %s", path_to_bigwig), "\n")
#    }
##TODO figure out if there is a way to normalize the samples 
#    # Determine the scale of the plotTracks plot by getting the max 
#    MAX <- -Inf
#    for (track in 1:length(all_tracks_to_plot)) {
#        print("Name of track")
#        print(names(all_tracks_to_plot[[track]]))
#      if(class(all_tracks_to_plot[[track]]) != "GenomeAxisTrack"){
#        if(max(all_tracks_to_plot[track][[1]]@data) > MAX) MAX <- max(all_tracks_to_plot[track][[1]]@data)
#      }
#    }
#   #Generate the plot
#    output_plot_name <- paste(plot_output_dir, "/", date_plot_created, "_", chromosome_to_plot, "_", descriptive_names_for_plots[experiment_index], ".svg", sep = "")
#    print("Output plot name: ")
#    print(output_plot_name)
#    #svg(output_plot_name)
#    #plotTracks(all_tracks_to_plot, main = main_title_of_plot_track, chromosome = df_sacCer_refGenome$chrom[chromosome_to_plot], ylim = c(0, MAX))
#    #dev.off()
#}
#for (sample_index in 1:nrow(df_sample_info)) {
#    all_tracks_to_plot <- list(GenomeAxisTrack(name = paste("Chr ", chromosome_to_plot, " Axis", sep = "") ) )
#    sample_ID_pattern <- df_sample_info$sample_ID[sample_index]
#    initial_matches <- list.files(bigwig_directory, pattern = as.character(sample_ID_pattern), full.names = TRUE, recursive = TRUE)     
#    path_to_bigwig <- initial_matches[grepl("S288C", initial_matches)]
#    print("Name of the bigwig path")
#    print(path_to_bigwig)
#    if (length(path_to_bigwig) > 0){
#        bigwig_to_plot <- import(con = path_to_bigwig, which = genomeRange_to_get)
#        sample_short_name <- df_sample_info$short_name[sample_index]
#        track_to_plot <- DataTrack(bigwig_to_plot, type = "l", name = sample_short_name, chromosome = df_sacCer_refGenome$chrom[chromosome_to_plot])
#        all_tracks_to_plot <- append(all_tracks_to_plot, track_to_plot)                                                                          
#
#        #Generate the plot
#        print("Name of the plot to be generated")
#        output_plot_name <- paste(plot_output_dir, "/", date_plot_created, "_", chromosome_to_plot, "_", df_sample_info$short_name[sample_index], ".svg", sep = "")
#        print(output_plot_name)
#        svg(output_plot_name)
#        plotTracks(all_tracks_to_plot, main = main_title_of_plot_track, chromosome = df_sacCer_refGenome$chrom[chromosome_to_plot], ylim = c(0, 100000))
#        dev.off()
#    }
#
#}
#print(utils::head(df_sample_info, 15)) 
#print(utils::tail(df_sample_info, 15)) 
