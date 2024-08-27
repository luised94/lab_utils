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
    feature_file_dir <- paste0(Sys.getenv("HOME"), "/data/feature_files")
    validate_input(feature_file_dir)
    files_to_convert <- get_file_list(feature_file_dir)
    cat(sprintf("Number of files to convert: %s", length(files_to_convert)), "\n")
    for (file_path in files_to_convert) {
        file_basename <- basename(file_path)
        data <- read_file(file_path)
        processed_data <- process_data(data, file_basename)
        grange_data <- convert_to_granges(processed_data)        
        output_processed_data(grange_data, file_basename, feature_file_dir)
    }
    cat("Transfer the files to dropbox for inspection\n")
    cat("scp -r user@server:from_dir to_dir\n")
}

validate_input <- function(input_dir) {
    cat("Ensuring directory exists", "\n")
    if (!dir.exists(input_dir)) {
        cat(sprintf("Directory %s does not exist.", input_dir))
        stop()
    }
}

get_file_list <- function(input_dir) {
    cat("Getting file list...\n")
    file_to_exclude <- list.files(path = input_dir) == "sample-key.tab"
    files <- list.files(path = input_dir, full.names = TRUE)[!file_to_exclude]
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
        tab = read.delim
    )

    if (!(file_extension %in% names(file_readers))) {
        stop(sprintf("Unsupported file type: %s", file_extension))
    }

    return(file_readers[[file_extension]](file_path))
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
        return(data)
    } else {
        return(data)
    }
}

convert_to_granges <- function(data) {
    cat("Converting to Granges...\n")
    if(!is(data, "GRanges")) {
        tryCatch({
            as(data, "GRanges")
        }, error = function(e) {
            warning("Could not convert to GRanges. Returning as-is.")
            return(data)
        })
    } else {
        return(data)
    }
}

output_processed_data <- function(data, file_name, output_dir) {
    cat(sprintf("Saving output for: %s to %s\n", file_name, output_dir))
    rds_output_file <- file.path(output_dir,
        paste0(tools::file_path_sans_ext(basename(file_name)), ".rds"))
    #saveRDS(data, file = rds_output_file)
    bed_output_file <- file.path(output_dir,
        paste0(tools::file_path_sans_ext(basename(file_name)), "_converted.bed"))
    #tryCatch(
    #    rtracklayer::export.bed(data, con = bed_output_file),
    #    error = function(e) stop(paste("Error exporting BED file:", e$message))
    #)
    cat(sprintf("Saved to: %s\n", bed_output_file))
    cat(sprintf("Saved to: %s\n", rds_output_file))
}

verify_output <- function(output_dir, pattern_in_file) {
    files_to_verify <- list.files(output_dir, pattern = pattern_in_file)
    cat(sprintf("Number of files to verify: %s\n", length(files_to_verify))

}
main()
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
#origin_GRange <- GRanges(seqnames = timing$Chromosome, ranges = IRanges(start = timing$Position-100,end = timing$Position+100), strand = "*", chromosome = df_sacCer_refGenome$chrom_ID[chromosome_to_plot], timing$T1/2)
#origin_GRange[seqnames(origin_GRange) == chromosome_to_plot]
#origin_track <- AnnotationTrack(origin_GRange[seqnames(origin_GRange) == chromosome_to_plot], name = "Origins")
#seqnames(origin_track) <- df_sacCer_refGenome$chrom_ID[chromosome_to_plot]
#timing <- readxl::read_excel(list.files(feature_file_directory, pattern = "timing", full.names =TRUE)
#) %>% as.data.frame %>% dplyr::select(1:7) %>% filter(!is.na(Chromosome))
#assign(genome_df$origin_gr_var[index_], GRanges(seqnames = genome_df$chrom[index_], ranges = IRanges(start = start, end = end), strand = "*"), envir = .GlobalEnv)
#assign(genome_df$origin_track_var[index_], AnnotationTrack(get(genome_df$origin_gr_var[index_]), name=paste(sacCer3_df$chrom[index_], "Origins", sep = " ")), envir = .GlobalEnv)
