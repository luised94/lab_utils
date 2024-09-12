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
    #Add process_control_factors, get_factors_to_match
    sample_table <- load_sample_table(directory_path)
    chromosome_to_plot = 10
    options(ucscChromosomeNames=FALSE)
    refGenome <- load_reference_genome(genome_dir = "REFGENS", genome_pattern = "S288C_refgenome.fna")

    genomeRange_to_get <- create_chromosome_GRange(refGenome = refGenome)

    #feature_file_pattern = "eaton_peaks"
    #feature_grange <- load_feature_file_GRange(chromosome_to_plot = chromosome_to_plot, feature_file_pattern = feature_file_pattern,
    #                                         genomeRange_to_get = genomeRange_to_get)

    #control_dir <- "EatonBel"
    #control_track <- load_control_grange_data(control_dir = control_dir, chromosome_to_plot = chromosome_to_plot,
    #                                         genomeRange_to_get = genomeRange_to_get)

    # Convert to input, add determine_matching_control, select_control_index
    #control_dir <- "EatonBel"
    #control_track <- load_control_grange_data(control_dir = control_dir, chromosome_to_plot = chromosome_to_plot,
                #                             genomeRange_to_get = genomeRange_to_get)

    #plot_all_sample_tracks(sample_table = sample_table,
    #                       directory_name = directory_path,
    #                       chromosome_to_plot = chromosome_to_plot, 
    #                       genomeRange_to_get = genomeRange_to_get, 
    #                       control_track = control_track, 
    #                       annotation_track = origin_track)

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

process_control_factors <- function(sample_table) {
    cat("Process control factors from __cf_ columns\n")
    df <- sample_table
    cf_cols <- grep("X__cf_", names(df), value = TRUE)
    if(length(cf_cols) == 0) {
        cat("No columns containing __cf_ tag found in sample table")
        stop("Verify sample table was produced with updated sampleGridConfig.")
    }
    control_factors <- lapply(df[cf_cols], function(x) strsplit(x[1], ",")[[1]])
    names(control_factors) <- sub("X__cf_", "", cf_cols)
    df[cf_cols] <- NULL
    attr(df, "control_factors") <- control_factors
    return(df)
}

get_factors_to_match <- function(sample_table, attribute_to_get = "control_factors") {
    cat("Grabbing attributes from sample table\n")
    df <- sample_table
    control_factors <- attr(df, attribute_to_get)
    if (is.null(control_factors)) {
        stop("No control factors defined in sample data.\nVerify 003_updateSampleGrid.R")
    }
    all_factors <- unlist(control_factors)
    return(intersect(all_factors, colnames(df)))
}

determine_matching_control <- function(sample_row, sample_table, factors_to_match) {
    cat("Determining control row for sample row.\n")
    df <- sample_table
    comparison_row <- sample_row[factors_to_match]
    rows_with_same_factors <- apply(df[, factors_to_match], 1, function(row) {
        all(row == comparison_row)
    })
    is_input <- df$antibody == "Input"
    index <- as.numeric(unname(which(is_input & rows_with_same_factors)))
    return(index)
}

select_control_index <- function(control_indices, max_controls = 1) {
    cat("Processing control index to ensure one is used.\n")
    if (length(control_indices) == 0) {
        warning("No matching control found")
        cat("Setting control_index to 1\n")
        control_indices <- 1
    }
    if (length(control_indices) > max_controls) {
    warning(paste("Multiple matching controls found, using first", max_controls))
    control_indices[1:max_controls]
    } else if (length(control_indices) == 1){
        return(control_indices)
    }
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
    cat(sprintf("Reading %s\n", sample_table_path))
    sample_table <- read.delim(sample_table_path, header = TRUE, sep ="\t")
    if (!("sample_ID" %in% colnames(sample_table))) {
           cat("Sample table does not contain the sample_ID column.\n")
           cat("sample_ID column is required for analysis.\n")
           cat("See 000_setupExperimentDir.R and 003_updateSampleGrid.R.\n")
           stop()
    }
    sample_table <- process_control_factors(sample_table)
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
create_chromosome_GRange <- function(refGenome) {
    cat("Creating chromosome GRange for loading feature, samples, etc\n")
    genomeRange_to_get <- GRanges(seqnames = refGenome$chrom,
                                  ranges = IRanges(start = 1, 
                                                   end = refGenome$basePairSize),
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
  normalized_chr_name <- switch(target_style,
    "UCSC" = paste0("chr", chr_names),
    "Roman" = sapply(chr_names, function(x) paste0("chr", ifelse(x %in% names(chr_to_roman), chr_to_roman[x], x))),
    "Numeric" = sapply(chr_names, function(x) ifelse(x %in% chr_to_roman, roman_to_chr[x], x)),
    stop("Unknown target style")
    )
    cat("Structure of normalized_chr_name\n")
    print(str(normalized_chr_name))
    return(unname(normalized_chr_name)) 
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
    feature_grange_subset <- subsetByOverlaps(feature_grange, genomeRange_to_get)
  } else {
    # Styles don't match, adjust genomeRange_to_get
    cat("Styles don't match. Adjusting genome range to match feature file.\n")
    adjusted_genomeRange <- genomeRange_to_get
    new_seqlevels <- normalize_chr_names(seqlevels(genomeRange_to_get), feature_style)
    seqlevels(adjusted_genomeRange) <- new_seqlevels

    feature_grange_subset <- subsetByOverlaps(feature_grange, adjusted_genomeRange)

    new_seqlevels <- normalize_chr_names(seqlevels(feature_grange_subset), genome_style)
    seqlevels(feature_grange_subset) <- new_seqlevels
    cat(sprintf("Confirming Feature GRange file style: %s\n", determine_chr_style(seqlevels(feature_grange_subset))))

  }
  return(feature_grange_subset)
}

load_control_grange_data <- function(control_dir, file_identifier, chromosome_to_plot = 10, genomeRange_to_get) {
    cat("Loading control track data from", control_dir, "\n")
    bigwig_dir_path <- file.path(Sys.getenv("HOME"), "data", control_dir, "bigwig")
    bigwig_file_paths <- list.files(bigwig_dir_path, pattern = file_identifier, full.names = TRUE, recursive = TRUE) 
    S288C_bigwigs <- grepl("S288C", bigwig_file_paths)
    bigwig_file_path <- bigwig_file_paths[S288C_bigwigs]
    if(!file.exists(bigwig_file_path)) {
        cat(sprintf("File %s doesnt exist.\n", bigwig_file_path))
        stop()
    }
    control_style <- determine_chr_style(seqlevels(import(bigwig_file_path)))
    chromosome_to_subset <- normalize_chr_names(chromosome_to_plot, control_style)
    subset_genome_range <- genomeRange_to_get[seqnames(genomeRange_to_get) == chromosome_to_subset]
    control_grange <- import(bigwig_file_path, which = genomeRange_to_get)
    return(control_grange)
}

plot_all_sample_tracks <- function(sample_table, directory_path, chromosome_to_plot = 10, genomeRange_to_get, control_track, annotation_track, highlight_gr) {
    main_title_of_plot_track <- paste("Complete View of Chrom", as.character(chromosome_to_plot), sep = " ")
    date_plot_created <- stringr::str_replace_all(Sys.time(), pattern = ":| |-", replacement="")  
    factors_to_match <- get_factors_to_match(sample_table)
    cat("Factors in attributes of sample_table\n")
    print(factors_to_match)
    plot_output_dir <- file.path(directory_path, "plots")
    bigwig_dir <- file.path(directory_path, "bigwig")
    cat("Plotting all sample tracks.\n")
    gtrack <- GenomeAxisTrack(name = paste("Chr ", chromosome_to_plot, " Axis", sep = ""))
    for (sample_index in 1:nrow(sample_table)) {
        if(sample_index == 1){
            cat("===============\n")
            sample_ID_pattern <- sample_table$sample_ID[sample_index]
            initial_matches <- list.files(bigwig_dir, pattern = as.character(sample_ID_pattern), full.names = TRUE, recursive = TRUE)
            path_to_bigwig <- initial_matches[grepl("S288C", initial_matches)]
            print("Name of the bigwig path")
            print(path_to_bigwig)
            chromosome_as_chr_roman <- paste("chr", as.roman(chromosome_to_plot), sep = "")
            if (length(path_to_bigwig) == 0){
                cat(sprintf("No bigwig found for sample_ID: %s\n", sample_ID_pattern))
                cat("Results of initial matches\n")
                print(initial_matches)
            }
            if (length(path_to_bigwig) > 0){
                control_index <- determine_matching_control(sample_row = sample_table[sample_index, ], sample_table, factors_to_match = factors_to_match)
                if(length(control_index) == 0) {
                    cat("No control index found\n")
                    cat("Printing sample row\n")
                    print(sample_table[sample_index, ])
                }
                control_index <- select_control_index(control_indices = control_index, max_controls = 1)
                control_ID_pattern <- sample_table$sample_ID[control_index]
                control_sample_name <-sample_table$short_name[control_index] 
                control_initial_matches <- list.files(bigwig_dir, pattern = as.character(control_ID_pattern), full.names = TRUE, recursive = TRUE)
                control_path_to_bigwig <- control_initial_matches[grepl("S288C", control_initial_matches)]
                if(length(control_path_to_bigwig) == 0){
                    cat("Appropriate control bigwig not found. Setting to first sample.\n")
                    control_ID_pattern <- sample_table$sample_ID[1]
                    control_sample_name <-sample_table$short_name[1] 
                    control_initial_matches <- list.files(bigwig_dir, pattern = as.character(control_ID_pattern), full.names = TRUE, recursive = TRUE)
                    control_path_to_bigwig <- control_initial_matches[grepl("S288C", control_initial_matches)]
                }
                print("Control ID pattern")
                print(control_ID_pattern)
                print("Name of the control bigwig path")
                print(control_path_to_bigwig)
                bigwig_to_plot <- import(con = path_to_bigwig, which = genomeRange_to_get)
                cat("Bigwig plot output\n")
                head(bigwig_to_plot)
                sample_short_name <- sample_table$short_name[sample_index]
                track_to_plot <- DataTrack(bigwig_to_plot, type = "l", name = sample_short_name, col = "#E41A1C", chromosome = chromosome_to_plot)
                sample_control_bigwig_to_plot <- import(con = control_path_to_bigwig, which = genomeRange_to_get)
                sample_control_track_to_plot <- DataTrack(sample_control_bigwig_to_plot, type = "l", name = control_sample_name, col = "#377EB8", chromosome = chromosome_to_plot)
                all_tracks <- list(gtrack, sample_control_track_to_plot, track_to_plot, control_track, annotation_track)
                output_plot_name <- paste(plot_output_dir, "/", date_plot_created, "_", chromosome_to_plot, "_", sample_short_name, "_", "WithInputAndEaton", ".svg", sep = "")
                print("Name of the plot to be generated")
                print(output_plot_name)
                #svg(output_plot_name)
                #plotTracks(all_tracks, main = main_title_of_plot_track, chromosome = chromosome_as_chr_roman, ylim = c(0, 100000))
                #dev.off()
                output_plot_name <- paste(plot_output_dir, "/", date_plot_created, "_", chromosome_to_plot, "_", sample_short_name, "_", "WithInputAndHighlights", ".svg", sep = "")
                print("Verifying GdObject")
                for (i in seq_along(all_tracks)) {
                  if (!inherits(all_tracks[[i]], "GdObject")) {
                      cat("Track", i, "is not a valid Gviz track object\n")
                          print(class(all_tracks[[i]]))
                    } else {
                        cat("Track", i, "inherits GdObject\n")
                        print(class(all_tracks[[i]]))
                        }
                }
                cat("===============\n")
                #tryCatch({
                sample_style <- determine_chr_style(seqlevels(track_to_plot))
                chromosome_to_subset <- normalize_chr_names(chromosome_to_plot, sample_style)
                subset_highlight_gr <- highlight_gr[seqnames(highlight_gr) == chromosome_to_subset]
                print("Subset highlight_gr")
                cat("===============\n")
                print(start(subset_highlight_gr))
                cat("===============\n")
                print(end(subset_highlight_gr))
                cat("===============\n")
                print(as.character(seqnames(subset_highlight_gr)))
                cat("===============\n")
                highlight_track <- HighlightTrack(tracklist = all_tracks,
                                        start = start(subset_highlight_gr),
                                        end = end(subset_highlight_gr),
                                        chromosome = as.character(seqnames(subset_highlight_gr)),
                                        fill = c("#FFE3E6", "#E6FFE3"),
                                        col = c("#FF0000", "#00FF00"),
                                        alpha = 0.3)
                #        }, error = function(e) {
                #                cat("Error creating HighlightTrack:", conditionMessage(e), "\n")
                #                })
                #tryCatch({
                print("Name of the plot to be generated")
                output_plot_name <- paste(plot_output_dir, "/", date_plot_created, "_", chromosome_to_plot, "_", sample_short_name, "_", "WithInputAndHighlights", ".svg", sep = "")
                print(output_plot_name)
                svg(output_plot_name)
                plotTracks(highlight_track, 
                #plotTracks(all_tracks, 
                            main = main_title_of_plot_track,
                            chromosome = chromosome_as_chr_roman,
                            ylim = c(0, 100000))
                dev.off()
               #         }, error = function(e) {
               #             cat("Error creating plotTracks:", conditionMessage(e), "\n")
               #             })
              }
        }
    }
    #cat("===============\n")
    #print(str(highlight_gr))
    #cat("===============\n")
    #print(start(highlight_gr))
    #cat("===============\n")
    #print(end(highlight_gr))
    #cat("===============\n")
    #print(as.character(seqnames(highlight_gr)))
    #cat("===============\n")
    #print(str(highlight_track))
    #cat("===============\n")
    #print(range(highlight_track))
    cat("Reached end of for loop\n")
    cat("Finished plotting samples\n")

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
    directory_path <- "240808Bel"
    cat("Logic of main function\n")
    main_function_logic <- main
    list_of_functions_and_variables <- ls()
    cat("Main function: Validating directory\n")
    directory_path <- validate_input(directory_path)
    chromosome_to_plot = 10
    #chromosome_to_plot = c(10, paste0("chr", as.roman(10)))
    cat("Main function: loading sample table\n")
    sample_table <- load_sample_table(directory_path)
    print(attr(sample_table, "control_factors"))
    cat("Main function: loading reference genome\n")
    refGenome <- load_reference_genome(genome_dir = "REFGENS", genome_pattern = "S288C_refgenome.fna")
    cat("Main function: Creating genomeRange_to_get\n")
    genomeRange_to_get <- create_chromosome_GRange(refGenome = refGenome)
    # Load peaks feature file for highlight track and annotation. 
    feature_file_pattern = "eaton_peaks"
    cat("Main function: load feature grange\n")
    feature_grange <- load_feature_file_GRange(chromosome_to_plot = chromosome_to_plot, feature_file_pattern = feature_file_pattern,genomeRange_to_get = genomeRange_to_get)
    feature_track <- AnnotationTrack(feature_grange, name = paste("Origin Peaks","Eaton 2010", sep = ""))
    print(head(feature_grange))
    # Load the Eaton Bel Track data as second comparison. 
    cat("Main function: loading control grange\n")
    control_dir <- "EatonBel"
    file_identifier <- "nnNnH"
    control_grange <- load_control_grange_data(control_dir = control_dir, file_identifier = file_identifier, chromosome_to_plot = chromosome_to_plot,genomeRange_to_get = genomeRange_to_get)
    control_track <- DataTrack(control_grange, type = "l",  name = "Eaton 2010", col = "#377EB8")
    print(head(control_grange))
    # Plot samples, determine the input control for each sample. No need to modify the files provided then. Just the logic.
    plot_all_sample_tracks(sample_table = sample_table,directory_path = directory_path,chromosome_to_plot = chromosome_to_plot,genomeRange_to_get = genomeRange_to_get,control_track = control_track,annotation_track = feature_track, highlight_gr = feature_grange)
}
