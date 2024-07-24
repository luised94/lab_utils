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
