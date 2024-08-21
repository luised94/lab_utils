#DESCRIPTION: Convert features files downloaded from Rossi 2021, Hawkins 2014 and Eaton 2010 to a uniform format for later use in plotting tracks and factor analysis.
#USAGE: ./004_processFeaturesToRdsAndBed.R or run each line from an interactive console if you want to verify the results before committing.
# Load libraries
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(readxl)
library(ChIPpeakAnno)
library(GenomicFeatures)

feature_file_dir <- paste0(Sys.getenv("HOME"), "/data/feature_files")
if (!dir.exists(feature_file_dir)) {
    cat(sprintf("Directory %s does not exist.", feature_file_dir))
    stop()
}
file_to_exclude <- list.files(path = feature_file_dir) == "sample-key.tab"
feature_files <- list.files(path = feature_file_dir, full.names = TRUE)[!file_to_exclude]

file_readers <- list(
    xls = readxl::read_excel,
    xlsx = readxl::read_excel,
    bed = rtracklayer::import.bed,
    gff3 = function(file_path) toGRanges(makeTxDbFromGFF(file_path), format = 'gene'),
    tab = read.delim
)

for (file_path in feature_files) {
    file_extension <- tools::file_ext(file_path)

    if(!(file_ext %in% names(file_readers))) {
        stop(paste("Unsupported file type:", file_ext))
    }

    data <- file_readers[[file_extension]](file_path)

    if (grepl("timing", basename(file_path)) {
        data <- data %>% 
            as.data.frame() %>% 
            dplyr::select(1:7) %>% 
            filter(!is.na(Chromosome))
    } else if (grepl("SGD_features", basename(file_path)) {
        data <- data %>% 
            as.data.frame() %>% 
            dplyr::select(c(9:12,2,4,5))
    }

    if (!is(data, "GRanges")) {
        data <- try(as(data, "GRanges"), silent = TRUE)
        if (inherits(data, "try-error")) {
            warning(paste("Could not convert", basename(file_path), "to Granges. Returning as-is."))
        }
    }

}
#
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
