#DESCRIPTION: 
#USAGE:
#TODO: Need some way to aggregate the column names I use from the sampleConfig.R files to keep track of variables I use to keep consistent experiment to experiment.
#TODO figure out if there is a way to normalize the samples 
##TODO: Define the comparisons and plots to be generated in my sampleConfig.R template. 
##TODO: Find the best way to have the same levels and factors when I read in the sampleGridConfig.R file. This is defined by the categories list variable. Can grab that during 002_loadSampleGrid and use it to equalize. 
##TODO: Use 002_loadSampleGrid to open the sample_table given a directory. Use conditional statement to determine the ID. Could potentially use system function to call find and print with awk statement. R would require list.files(), strsplit, and grabbing regular expression for 5 digits.
#Load packages using through a list of strings and suppress the messages, return a TRUE if loading was succesful
##TODO: Need to make sure the chromosome IDs are formatted properly.
##TODO: need to automate the creation of the experiments to plot. 
##TODO: Need to ensure that the names I use in the columns are compatible with my short name convention for subsetting df
main <- function() {
    #Load packages
    suppressPackageStartupMessages({
        library(QuasR)
        library(GenomicAlignments)
        library(Gviz)
        library(rtracklayer)
        library(ShortRead)
        library(tidyverse)
    })
    #Validate the arguments
    args <- commandArgs(trailingOnly = TRUE)
    directory_path <- validate_input(args)

    #Load the reference genome and proces into dataframe with chromosome and size columns.
    refGenome <- load_reference_genome()

    #Create the grange object to load grange from bigwig files for features, eaton control, and sample
    genomeRange_to_get <- create_chromosome_GRange(referenceGenome = refGenome, chromosome_to_plot = 10)

    #Load the origin AnnotationTrack from eaton acs using genomeRange_to_get
    origin_track <- load_origin_annotation_track(feature_file_path = "240830Bel_eaton_peaks.bed", genomeRange_to_get = genomeRange_to_get)

    #Load the eaton track data to overlay with sample data.
    eaton_track <- load_control_track_data(eaton_path = "EatonBel", genomeRange_to_get = genomeRange_to_get)
    #Load the sample table.
    sample_table <- load_sample_table(directory_path)

    #Plot all tracks from the sample_table using the origin track for annotation and the eaton track to create an overlayed track
    plotAllSamples(sample_table, annotationtrack = origin_track, overlaytrack = eaton_track)

    # Has to be run every time if you are using chromosome ID to get tracks from the bigwig file.
    #options(ucscChromosomeNames=FALSE)
}

validate_input(args) {
    if (length(args) != 1) {
        cat("Error: Invalid number of arguments.\n")
        cat("Usage: Rscript 001_plotAllSampleTracks.R <directory_path>\n")
        cat("Example: Rscript 001_plotAllSampleTracks.R 240819Bel\n")
        stop()
    }
    directory_path <- file.path(Sys.getenv("HOME"), "data", args[1])
    if(!dir.exists(directory_path)) {
        cat(sprintf("Error: Directory %s does not exist.\n", directory_path))
        stop()
    }
    return(directory_path)
}
load_reference_genome(genome_dir = "REFGENS", genome_pattern = "S288C_refgenome.fna") {
    cat("Loading reference genome\n")
    directory_of_refgenomes <- file.path(Sys.getenv("HOME"), "data", genome_dir)
    if(!dir.exists(directory_of_refgenomes)) {
        stop("Directory with reference genomes doesnt exist.\n")
    }
    genome_file_path <- list.files(directory_of_refgenomes, pattern = genome_pattern, full.names = TRUE, recursive = TRUE)
    if (length(genome_file_path) > 1) {
        cat(sprintf("More than one file matched genome pattern %s", genome_pattern))
        print(genome_file_path)
        stop()
    }
    if(!file.exists(genome_file_path)) {
        stop("Reference genome doesnt exist.\n")
    }
    refGenome <- readFasta(genome_file_path)
    refGenome <- data.frame(chrom = names(as(refGenome, "DNAStringSet")), 
                    basePairSize = width(refGenome)) %>% filter(chrom != "chrM")
    return(refGenome)
}

load_sample_table <- function(directory_name) {
    cat("Loading sample_table from", directory_name, "\n")
    documentation_dir_path <- file.path(directory_name, "documentation")
    sample_table_path <- list.files(documentation_dir_path, pattern = "sample_table", full.names = TRUE)
    if(length(sample_table_path) == 0){
        cat(sprintf("No files with pattern sample_table found in %s\n", documentation_dir_path))
        cat("Consult 000_setupExperimentDir to create sample_table.tsv for the project\n")
        stop()
    } else if(length(sample_table_path) > 1){
        cat(sprintf("Multiple files with pattern sample_table found in %s\n", documentation_dir_path))
        cat(sprintf("Files found in %s\n", documentation_dir_path))
        print(sample_table)
        cat("Consult 000_setupExperimentDir and ensure no duplicates are present for the project\n")
        stop()
    }
    sample_table <- read.delim(sample_table_path, header = TRUE, sep ="\t")
    if (!("sample_ID" %in% colnames(sample_table)) {
           cat("Sample table does not contain the sample_ID column.\n")
           cat("sample_ID column is required for analysis.\n")
           cat("See 000_setupExperimentDir.R and 003_updateSampleGrid.R.\n")
           stop()
    }
    cat(sprintf("Reading %s\n", sample_table_path))
    cat("Head of sample_table\n")
    print(head(sample_table))
    return(sample_table)
}

load_feature_file <- function(chromosome_to_plot, feature_file_name = "240830_eaton_peaks.bed", chromosome_grange) {
create_chromosome_GRange <- function(refGenome, chromosome_to_plot){
refGenome <- readFasta(genome_file_path)
refGenome <- data.frame(chrom = names(as(refGenome, "DNAStringSet")), basePairSize = width(refGenome)) %>% filter(chrom != "chrM")
genomeRange_to_get <- GRanges(seqnames=10,ranges = IRanges(start = 1,end = refGenome$basePairSize[10]),strand = "*")
origin_grange <- import.bed("240830_eaton_peaks.bed", which = genomeRange_to_get)
origin_track <- AnnotationTrack(origin_grange, name = "Origin Peaks (Eaton et al 2010)")
#Create GRanges object to read in a particular chromosome
chromosome_to_plot <- 10
genomeRange_to_get <- GRanges(seqnames=c(df_sacCer_refGenome$chrom[chromosome_to_plot]), 
        ranges = IRanges(start = 1, 
        end = df_sacCer_refGenome$basePairSize[chromosome_to_plot]), 
        strand = "*")
}
load_feature_file_GRange <- function(chromosome_to_plot, feature_file_pattern = "eaton_acs")
main_title_of_plot_track <- paste("Complete View of Chrom", as.character(chromosome_to_plot, 
                  sep = " "))
date_plot_created <- stringr::str_replace_all(Sys.time(), pattern = ":| |-", replacement="")  
bigwig_directory <- paste(working_directory, "bigwig", sep = "/")
for (sample_index in 1:nrow(df_sample_info)) {
    all_tracks_to_plot <- list(GenomeAxisTrack(name = paste("Chr ", chromosome_to_plot, " Axis", sep = "") ) )
    sample_ID_pattern <- df_sample_info$sample_ID[sample_index]
    initial_matches <- list.files(bigwig_directory, pattern = as.character(sample_ID_pattern), full.names = TRUE, recursive = TRUE)     
    path_to_bigwig <- initial_matches[grepl("S288C", initial_matches)]
    print("Name of the bigwig path")
    print(path_to_bigwig)
    if (length(path_to_bigwig) > 0){
        bigwig_to_plot <- import(con = path_to_bigwig, which = genomeRange_to_get)
        sample_short_name <- df_sample_info$short_name[sample_index]
        track_to_plot <- DataTrack(bigwig_to_plot, type = "l", name = sample_short_name, chromosome = df_sacCer_refGenome$chrom[chromosome_to_plot])
        all_tracks_to_plot <- append(all_tracks_to_plot, track_to_plot)                                                                          

        #Generate the plot
        print("Name of the plot to be generated")
        output_plot_name <- paste(plot_output_dir, "/", date_plot_created, "_", chromosome_to_plot, "_", df_sample_info$short_name[sample_index], ".svg", sep = "")
        print(output_plot_name)
        svg(output_plot_name)
        plotTracks(all_tracks_to_plot, main = main_title_of_plot_track, chromosome = df_sacCer_refGenome$chrom[chromosome_to_plot], ylim = c(0, 100000))
        dev.off()
    }

}
