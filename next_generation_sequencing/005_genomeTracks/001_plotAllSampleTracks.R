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
    args <- commandArgs(trailingOnly = TRUE)
    directory_path <- validate_input(args)
    sample_table <- load_sample_table(directory_path)
    chromosome_to_plot = 10
    options(ucscChromosomeNames=FALSE)
    refGenome <- load_reference_genome(genome_dir = "REFGENS", genome_pattern = "S288C_refgenome.fna")

    genomeRange_to_get <- create_chromosome_GRange(refGenome = refGenome, chromosome_to_plot = chromosome_to_plot)

    feature_file_pattern = "eaton_peaks"
    feature_grange <- load_feature_file_GRange(chromosome_to_plot = chromosome_to_plot, feature_file_pattern = feature_file_pattern,
                                             genomeRange_to_get = genomeRange_to_get)

    control_dir <- "EatonBel"
    control_track <- load_control_grange_data(control_dir = control_dir, chromosome_to_plot = chromosome_to_plot,
                                             genomeRange_to_get = genomeRange_to_get)

    plot_all_sample_tracks(sample_table = sample_table,
                           directory_name = directory_path,
                           chromosome_to_plot = chromosome_to_plot, 
                           genomeRange_to_get = genomeRange_to_get, 
                           control_track = control_track, 
                           annotation_track = origin_track)

}

validate_input <- function(args) {
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
        print(sample_table_path)
        cat("Consult 000_setupExperimentDir and ensure no duplicates are present for the project\n")
        stop()
    }
    sample_table <- read.delim(sample_table_path, header = TRUE, sep ="\t")
    if (!("sample_ID" %in% colnames(sample_table))) {
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

load_reference_genome <- function(genome_dir = "REFGENS", genome_pattern = "S288C_refgenome.fna") {
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
    cat("Head of refGenome.\n")
    print(head(refGenome))
    return(refGenome)
}

#Create GRanges object to read in a particular chromosome
create_chromosome_GRange <- function(refGenome, chromosome_to_plot = 10) {
    cat("Creating chromosome GRange for loading feature, samples, etc\n")
    genomeRange_to_get <- GRanges(seqnames = refGenome$chrom[chromosome_to_plot],
                                  ranges = IRanges(start = 1, 
                                                   end = refGenome$basePairSize[chromosome_to_plot]),
                                  strand = "*")
    cat("Head of Genome Range for loading other files.\n")
    print(head(genomeRange_to_get))
    return(genomeRange_to_get)
}
# Chromosome mapping functions
chr_to_roman <- c(
  "1" = "I", "2" = "II", "3" = "III", "4" = "IV", "5" = "V", "6" = "VI", "7" = "VII", "8" = "VIII",
  "9" = "IX", "10" = "X", "11" = "XI", "12" = "XII", "13" = "XIII", "14" = "XIV", "15" = "XV", "16" = "XVI"
)

roman_to_chr <- setNames(names(chr_to_roman), chr_to_roman)

normalize_chr_names <- function(chr_names, target_style) {
  chr_names <- gsub("^chr", "", chr_names)
  if (target_style == "UCSC") {
    return(paste0("chr", chr_names))
  } else if (target_style == "Roman") {
    return(sapply(chr_names, function(x) paste0("chr", ifelse(x %in% names(chr_to_roman), chr_to_roman[x], x))))
  } else if (target_style == "Numeric") {
    return(sapply(chr_names, function(x) ifelse(x %in% chr_to_roman, roman_to_chr[x], x)))
  }
}

determine_chr_style <- function(chr_names) {
  if (all(grepl("^chr[0-9]+$", chr_names))) return("UCSC")
  if (all(grepl("^chr[IVX]+$", chr_names))) return("Roman")
  if (all(grepl("^[0-9]+$", chr_names))) return("Numeric")
  return("Unknown")
}

load_feature_file_GRange <- function(chromosome_to_plot = 10, feature_file_pattern = "eaton_peaks", genomeRange_to_get) {
  cat(sprintf("Loading %s feature file.\n", feature_file_pattern))
  
  # Input validation
  feature_file_dir <- file.path(Sys.getenv("HOME"), "data", "feature_files")
  if(!dir.exists(feature_file_dir)) {
    stop(sprintf("Directory %s does not exist.", feature_file_dir))
  }
  feature_file_path <- list.files(feature_file_dir, pattern = feature_file_pattern, full.names = TRUE, recursive = TRUE)
  if(length(feature_file_path) != 1) {
    stop(sprintf("Error finding feature file. Found %d files: %s", length(feature_file_path), paste(feature_file_path, collapse = ", ")))
  }
  # Load feature file and determine its style
  feature_grange <- import.bed(feature_file_path)
  feature_style <- determine_chr_style(seqlevels(feature_grange))
  cat("Feature file chromosome style:", feature_style, "\n")
  # Determine genomeRange style
  genome_style <- determine_chr_style(seqlevels(genomeRange_to_get))
  cat("Genome range chromosome style:", genome_style, "\n")
  if (feature_style == genome_style) {
    # Styles match, use genomeRange_to_get as is
    cat("Styles match. Using provided genome range.\n")
    feature_grange <- subsetByOverlaps(feature_grange, genomeRange_to_get)
  } else {
    # Styles don't match, adjust genomeRange_to_get
    cat("Styles don't match. Adjusting genome range to match feature file.\n")
    adjusted_genomeRange <- genomeRange_to_get
    seqlevels(adjusted_genomeRange) <- normalize_chr_names(seqlevels(genomeRange_to_get), feature_style)
    seqnames(adjusted_genomeRange) <- normalize_chr_names(as.character(seqnames(genomeRange_to_get)), feature_style)
    
    # Subset feature_grange using adjusted genomeRange
    feature_grange <- subsetByOverlaps(feature_grange, adjusted_genomeRange)
    
    # Convert back to original genome style
    #seqlevels(feature_grange_subset) <- normalize_chr_names(seqlevels(feature_grange_subset), genome_style)
    #seqnames(feature_grange_subset) <- normalize_chr_names(as.character(seqnames(feature_grange_subset)), genome_style)
  }
  
  return(feature_grange)
}


#load_feature_file_GRange <- function(chromosome_to_plot = 10, feature_file_pattern = "eaton_peaks", genomeRange_to_get) {
#    cat(sprintf("Loading %s feature file.\n", feature_file_pattern))
#    feature_file_dir <- file.path(Sys.getenv("HOME"), "data", "feature_files")
#    if(!dir.exists(feature_file_dir)) {
#        cat(sprintf("Directory %s does not exist.\n", feature_file_dir))
#        stop()
#    }
#    feature_file_path <- list.files(feature_file_dir, pattern = feature_file_pattern, full.names = TRUE, recursive = TRUE)
#    if(length(feature_file_path) != 1) {
#        cat(sprintf("Error finding feature file. A single file wasn't identified."))
#        print(feature_file_path)
#        stop()
#    }
#    cat("Seq levels of feature file and style.\n")
#    #seqlevels(genomeRange_to_get) <- as.character(as.numeric(roman2int(mapSeqlevels(as.character(seqnames(genomeRange_to_get)), style = "NCBI"))))
#    feature_grange <- import.bed(feature_file_path)
#    print(seqlevels(feature_grange))
#    print(seqlevelsStyle(feature_grange))
#    #seqnames(genomeRange_to_get) <- chromosome_to_plot
#    #feature_grange <- import.bed(feature_file_path, which = genomeRange_to_get)
#
#    #seqnames(feature_grange) <- paste("chr", as.roman(chromosome_to_plot), sep = "")
#    annotation_track <- AnnotationTrack(feature_grange, name = "Origin Peaks (Eaton2010)")
#    cat("Annotation Track for feature file.\n")
#    print(head(annotation_track))
#    return(annotation_track)
#}

load_control_grange_data <- function(control_dir = "EatonBel", chromosome_to_plot = 10, genomeRange_to_get) {
    cat("Loading control track data from", control_dir, "\n")
    bigwig_dir_path <- file.path(Sys.getenv("HOME"), "data", control_dir, "bigwig")
    bigwig_file_path <- list.files(bigwig_dir_path, pattern = "S288C", full.names = TRUE, recursive = TRUE)
    if(!file.exists(bigwig_file_path)) {
        cat(sprintf("File %s doesnt exist.\n", bigwig_file_path))
        stop()
    }
    #control_grange <- import.bed(bigwig_file_path, which = genomeRange_to_get)
    control_grange <- import(bigwig_file_path, which = genomeRange_to_get)
    return(control_track)
}

plot_all_sample_tracks <- function(sample_table, directory_name, chromosome_to_plot = 10, genomeRange_to_get, control_track, annotation_track) {
    main_title_of_plot_track <- paste("Complete View of Chrom", as.character(chromosome_to_plot), sep = " ")
    date_plot_created <- stringr::str_replace_all(Sys.time(), pattern = ":| |-", replacement="")  
    cat("Plotting all sample tracks.\n")
    plot_output_dir <- file.path(Sys.getenv("HOME"), "data", directory_name, "plots")
    for (sample_index in 1:nrow(sample_table)) {
        gtrack <- GenomeAxisTrack(name = paste("Chr ", chromosome_to_plot, " Axis", sep = ""))
        sample_ID_pattern <- sample_table$sample_ID[sample_index]
        initial_matches <- list.files(bigwig_directory, pattern = as.character(sample_ID_pattern), full.names = TRUE, recursive = TRUE)
        path_to_bigwig <- initial_matches[grepl("S288C", initial_matches)]
        print("Name of the bigwig path")
        print(path_to_bigwig)
        chromosome_as_chr_roman <- paste("chr", as.roman(chromosome_to_plot), sep = "")
        if (length(path_to_bigwig) > 0){
            bigwig_to_plot <- import(con = path_to_bigwig, which = genomeRange_to_get)
            cat("Bigwig plot output\n")
            head(bigwig_to_plot)
            sample_short_name <- sample_table$short_name[sample_index]
            track_to_plot <- DataTrack(bigwig_to_plot, type = "l", name = sample_short_name, chromosome = chromosome_to_plot)
            overlay <- OverlayTrack(trackList = list(control_track, track_to_plot))
            cat("Overlay object\n")
            head(overlay)
            #Generate the plot
            print("Name of the plot to be generated")
            output_plot_name <- paste(plot_output_dir, "/", date_plot_created, "_", chromosome_to_plot, "_", sample_short_name, "_WithEaton", ".svg", sep = "")
            print(output_plot_name)
#            svg(output_plot_name)
#            plotTracks(list(gtrack, overlay), main = main_title_of_plot_track, chromosome = chromosome_as_chr_roman, ylim = c(0, 100000))
#            dev.off()
        }
    }

}

if(!interactive()){
    main()
} else {
    suppressPackageStartupMessages({
        library(QuasR)
        library(GenomicAlignments)
        library(Gviz)
        library(rtracklayer)
        library(ShortRead)
        library(tidyverse)
        library(gtools)
    })
    directory_path <- "240819Bel"
    cat("Logic of main function\n")
    main_function_logic <- main
    list_of_functions_and_variables <- ls()
    directory_path <- validate_input(directory_path)
    chromosome_to_plot = 10
    #chromosome_to_plot = c(10, paste0("chr", as.roman(10)))
    sample_table <- load_sample_table(directory_path)
    refGenome <- load_reference_genome(genome_dir = "REFGENS", genome_pattern = "S288C_refgenome.fna")
    genomeRange_to_get <- create_chromosome_GRange(refGenome = refGenome, chromosome_to_plot = chromosome_to_plot)
    feature_file_pattern = "eaton_peaks"
    feature_grange <- load_feature_file_GRange(chromosome_to_plot = chromosome_to_plot, feature_file_pattern = feature_file_pattern,genomeRange_to_get = genomeRange_to_get)
    feature_track <- AnnotationTrack(feature_grange, name = paste("Origin Peaks","Eaton 2010", sep = ""))
  
    #seqlevelsStyle, genomeStyles, extractSeqlevels, mapSeqlevels
    control_dir <- "EatonBel"
    control_grange <- load_control_grange_data(control_dir = control_dir, chromosome_to_plot = chromosome_to_plot,genomeRange_to_get = genomeRange_to_get)

    control_track <- DataTrack(control_grange, name = "Eaton 2010")
    #plot_all_sample_tracks(sample_table = sample_table,directory_name = directory_path,chromosome_to_plot = chromosome_to_plot,genomeRange_to_get = genomeRange_to_get,control_track = control_track,annotation_track = origin_track)
}
    #print(seqlevels(genomeRange_to_get))
    #print(as.character(as.numeric(roman2int(mapSeqlevels(as.character(seqnames(genomeRange_to_get)), style = "NCBI")))))
    #seqlevels(genomeRange_to_get) <- as.character(as.numeric(roman2int(mapSeqlevels(as.character(seqnames(genomeRange_to_get)), style = "NCBI"))))
    #print(seqlevels(genomeRange_to_get))
