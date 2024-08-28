#DESCRIPTION: Convert features files from Rossi 2021, Hawkins 2014, and Eaton 2010 to a bed and rds format for plotting tracks and factor analysis.
#USAGE: Rscript ./004_processFeaturesToRdsAndBed.R 
# Load libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(GenomicRanges)
    library(rtracklayer)
    library(readxl)
    library(ChIPpeakAnno)
    library(GenomicFeatures)
})

main <- function(input_dir) {
    cat("Entering main function. All libraries loaded.\n")
    feature_file_dir <- paste0(Sys.getenv("HOME"), "/data/feature_files")
    validate_input(feature_file_dir)
    pattern_to_exclude  <- "sample-key.tab|\\.rds$|\\_converted.bed$"
    files_to_convert <- get_file_list(feature_file_dir, pattern_to_exclude = pattern_to_exclude)
    for (file_path in files_to_convert) {
        file_basename <- basename(file_path)
        data <- read_file(file_path)
        processed_data <- process_data(data, file_basename)
        grange_data <- convert_to_granges(processed_data, file_basename)
        output_processed_data(grange_data, file_basename, feature_file_dir)
    }
    pattern_to_verify <- "\\.rds$|\\_converted.bed$"
    verify_output(feature_file_dir, pattern_in_file = pattern_to_verify)
    cat("Transfer the files to dropbox for inspection\n")
    cat("scp -r user@server:from_dir to_dir\n")
}

validate_input <- function(input_dir) {
    cat(sprintf("Ensuring %s directory exists\n", input_dir))
    if (!dir.exists(input_dir)) {
        cat(sprintf("Directory %s does not exist.\n", input_dir))
        stop()
    }
}

get_file_list <- function(input_dir, pattern_to_exclude = "sample-key.tab|\\.rds$|\\_converted.bed$") {
    cat("Getting file list...\n")
    all_files <- list.files(path = input_dir)
    file_to_exclude <- grepl(pattern_to_exclude, all_files)
    files <- list.files(path = input_dir, full.names = TRUE)[!file_to_exclude]
    cat(sprintf("Number of files to process: %s\n", length(files)))
    return(files)
}

read_file <- function(file_path) {
    cat(sprintf("Reading file: %s\n", basename(file_path)))
    file_extension <- tools::file_ext(file_path)
    file_readers <- list(
        xls = readxl::read_excel,
        xlsx = readxl::read_excel,
        bed = rtracklayer::import.bed,
        gff3 = function(file_path) import(file_path, format = 'gff3'),
        tab = function(file_path) read.delim(file_path, header = FALSE, sep = "\t"),
        rds = readRDS
    )

    if (!(file_extension %in% names(file_readers))) {
        stop(sprintf("Unsupported file type: %s", file_extension))
    }
    data <- file_readers[[file_extension]](file_path)
    return(data)
}

process_data <- function(data, file_name) {
    cat(sprintf("Processing data for: %s\n", file_name))
    if(grepl("timing", file_name)) {
        data <- data %>%
                as.data.frame() %>%
                dplyr::select(1:7) %>%
                filter(!is.na(Chromosome))
        return(data)
    } else if (grepl("SGD_features", file_name)) {
        data <- data %>%
                as.data.frame() %>%
                dplyr::select(c(9:12, 2, 4, 5))
        print(head(data))
        return(data)
    } else {
        return(data)
    }
}

convert_to_granges <- function(data, file_basename) {
    cat("Converting to Granges...\n")
    if(!is(data, "GRanges")) {
        if (grepl("Nucleosome_calls", file_basename)) {
            cat("Processing Nucleosome calls xlsx file\n")
            columns_to_exclude <- c("Nucleosome ID", "Nucleosome dyad", "Chromosome")
            metadata_dataframe <- data.frame(data[, !(colnames(data) %in% columns_to_exclude)])
            colnames(metadata_dataframe) <- gsub(" |\\.", "_", colnames(metadata_dataframe))
            data <- GRanges(seqnames = data$`Nucleosome ID`, ranges = IRanges(start = data$`Nucleosome dyad`, end = data$`Nucleosome dyad`), strand = "*", chromosome = data$Chromosome, metadata_dataframe)
            print(head(data))
            return(data)
        } else if (grepl("hawkins", file_basename)) {
            cat("Processing hawkins timing xlsx file\n")
            name_of_origins <- paste0(data$Chromosome, "_", data$Position)
            columns_to_exclude <- c()
            metadata_dataframe <- data.frame(data[, !(colnames(data) %in% columns_to_exclude)])
            colnames(metadata_dataframe) <- gsub(" |\\.", "_", colnames(metadata_dataframe))
            data <- GRanges(seqnames = name_of_origins, ranges = IRanges(start = data$Position-100,end = data$Position+100), strand = "*", chromosome = data$Chromosome, metadata_dataframe)
            cat("Finished processing hawkins xlsx into grange.\n")
            print(head(data))
            return(data)
        } else if (grepl("Rhee", file_basename)) {
            cat("Processing Rhee transcription data file\n")
            data <- data[!apply(data, 1, function(row) all(is.na(row))), ]
            columns_to_exclude <- c("chrom", "TATA_coor", "strand", colnames(data)[grepl("color =class;", colnames(data))])
            metadata_dataframe <- data.frame(data[, !(colnames(data) %in% columns_to_exclude)])
            colnames(metadata_dataframe) <- gsub(" |\\.", "_", colnames(metadata_dataframe))
            data$strand <- ifelse(data$strand == "W", "+", ifelse(data$strand == "C", "-", ifelse(data$strand, "*")))
            data <- GRanges(seqnames = data$chrom, ranges = IRanges(start = data$TATA_coor, end = data$TATA_coor), strand = data$strand, chromosome = data$chrom, metadata_dataframe)
            cat("Finished processign Rhee data file into grange")
            print(head(data))
            return(data)
        } else if (grepl("SGD", file_basename)) {
            data <- data[!apply(data, 1, function(row) all(is.na(row))), ]
            columns_to_exclude <- c("chrom", "TATA_coor", "strand", colnames(data)[grepl("color =class;", colnames(data))])
            metadata_dataframe <- data.frame(data[, !(colnames(data) %in% columns_to_exclude)])
            colnames(metadata_dataframe) <- gsub(" |\\.", "_", colnames(metadata_dataframe))
            data$strand <- ifelse(data$strand == "W", "+", ifelse(data$strand == "C", "-", ifelse(data$strand, "*")))
            data <- GRanges(seqnames = data$chrom, ranges = IRanges(start = data$TATA_coor, end = data$TATA_coor), strand = data$strand, chromosome = data$chrom, metadata_dataframe)
            cat("Finished processign Rhee data file into grange")
            print(head(data))
            return(data)
        } else {
            tryCatch({
                as(data, "GRanges")
            }, error = function(e) {
                warning("Could not convert to GRanges. Returning as-is.")
                print(head(data))
                return(data)
            })
            print(head(data))
       }
    } else {
        print(head(data))
        return(data)
    }
}

output_processed_data <- function(data, file_name, output_dir) {
    cat(sprintf("Saving output for: %s to %s\n", file_name, output_dir))
    rds_output_file <- file.path(output_dir,
        paste0(tools::file_path_sans_ext(basename(file_name)), ".rds"))
   # saveRDS(data, file = rds_output_file)
    bed_output_file <- file.path(output_dir,
        paste0(tools::file_path_sans_ext(basename(file_name)), "_converted.bed"))
    #tryCatch(
    #    rtracklayer::export.bed(data, con = bed_output_file),
    #    error = function(e) stop(paste("Error exporting BED file:", e$message))
    #)
    cat(sprintf("Saved to: %s\n", bed_output_file))
    cat(sprintf("Saved to: %s\n", rds_output_file))
}

verify_output <- function(output_dir, pattern_in_file = "_converted\\.bed") {
    cat("Verifying generated output\n")
    files_to_verify <- list.files(output_dir, pattern = pattern_in_file)
    cat(sprintf("Number of files to verify: %s\n", length(files_to_verify)))
    #file_readers <- list(
    #        bed = rtracklayer::import,
    #        rds = readRDS
    #)

    #file_converter <- list(
    #        rds = function(file_path) makeGRangesFromDataFrame(file_path, keep.extra.columns = TRUE)
    #)
    if (length(files_to_verify) == 0) {
        cat("No files to verify")
        return(0)
    }
    for (file_path in files_to_verify) {
         data <- rtracklayer::import.bed(file_path)
         cat(sprintf("Successfully read: %s\n", basename(file_path)))
        if (is(data, "GRanges")) {
            cat(sprintf("File %s is a valid GRanges object\n", basename(file_path)))
            cat("Sample data:\n")
            print(head(data))
        } else {
            cat(sprintf("Warning: File %s is not a GRanges object\n", basename(file_path)))
            print(head(data))
        }
    }
}

main()
    #        file_extension <- tools::file_ext(file_path)
    #        tryCatch({
    #            data <- file_readers[[file_extension]](file_path)
    #            cat(sprintf("Successfully read: %s\n", basename(file_path)))
    #            if (is(data, "GRanges")) {
    #                cat(sprintf("File %s is a valid GRanges object\n", basename(file_path)))
    #                cat("Sample data:\n")
    #                print(head(data))
    #            } else {
    #                cat(sprintf("Warning: File %s is not a GRanges object\n", basename(file_path)))
    #                data <- file_converter[[file_extension]](data)
    #                if(is(data, "GRanges")) {
    #                    cat(sprintf("File %s is a valid GRanges object\n", basename(file_path)))
    #                    print(head(data))
    #                }
    #            }
    #        }, error = function(e) {
    #            cat(sprintf("Error reading file %s: %s\n", basename(file_path), e$message))
    #        })
    #        cat("\n")

#for (file_path in feature_files) {
#    file_extension <- tools::file_ext(file_path)
#
#    if(!(file_ext %in% names(file_readers))) {
#        stop(paste("Unsupported file type:", file_ext))
#    }
#
#    data <- file_readers[[file_extension]](file_path)
#
#    if (grepl("timing", basename(file_path)) {
#        data <- data %>% 
#            as.data.frame() %>% 
#            dplyr::select(1:7) %>% 
#            filter(!is.na(Chromosome))
#    } else if (grepl("SGD_features", basename(file_path)) {
#        data <- data %>% 
#            as.data.frame() %>% 
#            dplyr::select(c(9:12,2,4,5))
#    }
#
#    if (!is(data, "GRanges")) {
#        data <- try(as(data, "GRanges"), silent = TRUE)
#        if (inherits(data, "try-error")) {
#            warning(paste("Could not convert", basename(file_path), "to Granges. Returning as-is."))
#        }
#    }
#
#}
##
#df_names <- c("ChExMix","MEME","XUTs","CUTs","ORF","Nucleosome", "timing", "Rhee" , "SGD_features", "SUTs", "G2_orc")
#
#for(i in 1:length(df_names)){
#  print(feature_files[grepl(df_names[i], feature_files)])
#  filepath <- file.path(feature_folder, feature_files[grepl(df_names[i], feature_files)])
#  if (grepl("xls|xlsx", tools::file_ext(filepath))){
#    assign(df_names[i], readxl::read_excel(filepath))
#    if(df_names[i] == "timing"){
#      assign(df_names[i], get(df_names[i]) %>% as.data.frame() %>% dplyr::select(1:7) %>% filter(!is.na(Chromosome)))
#    }
#  } else if(tools::file_ext(filepath) == "bed"){
#    assign(df_names[i], readBed(filepath))
#    print(seqlevelsStyle(eval(parse(text = df_names[i]))))
#  } else if(tools::file_ext(filepath) == "gff3"){
#    assign(df_names[i], toGRanges(makeTxDbFromGFF(filepath), format = 'gene'))
#    print(seqlevelsStyle(eval(parse(text = df_names[i]))))
#  } else if(tools::file_ext(filepath) == "tab"){
#    assign(df_names[i], read.delim(filepath, header = FALSE))
#    if(df_names[i] == "SGD_features"){
#      assign(df_names[i], get(df_names[i]) %>% as.data.frame() %>% dplyr::select(c(9,10,11,12,2,4,5)))
#    }
#  }
#}
#
#AnnotationTrack(variable_name, name = paste(chromosome, "origins", sep = " "))
#assign(df_names[i], readxl::read_excel(filepath))
#if(df_names[i] == "timing"){
#assign(df_names[i], get(df_names[i]) %>% as.data.frame() %>% dplyr::select(1:7) %>% filter(!is.na(Chromosome)))}
#origin_GRange[seqnames(origin_GRange) == chromosome_to_plot]
#origin_track <- AnnotationTrack(origin_GRange[seqnames(origin_GRange) == chromosome_to_plot], name = "Origins")
#seqnames(origin_track) <- df_sacCer_refGenome$chrom_ID[chromosome_to_plot]
#timing <- readxl::read_excel(list.files(feature_file_directory, pattern = "timing", full.names =TRUE)
#) %>% as.data.frame %>% dplyr::select(1:7) %>% filter(!is.na(Chromosome))

#assign(genome_df$origin_track_var[index_], AnnotationTrack(get(genome_df$origin_gr_var[index_]), name=paste(sacCer3_df$chrom[index_], "Origins", sep = " ")), envir = .GlobalEnv)
