#STATUS:
#DESCRIPTION: Read in and analyse the quality control files for Bam.
#USAGE: 003_checkQcBam.R <dir>
#Run with: Rscript lab_utils/next_generation_sequencing/004_bamProcessing/003_R_node_checkQcBam.R dirname > output.log 2>&1
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 ) { print("No argument provided. Provide a directory to process") ; q() }

get_current_datetime_string <- function() {
                  return(format(Sys.time(), "%Y-%m-%d-%H-%M-%S"))
}

find_the_directory <- function(dirname){
    data_directory <- paste0(Sys.getenv("HOME"), "/data")
    contains_dirname <- grepl(dirname, list.dirs(data_directory, recursive = FALSE))
    directory_path <- list.dirs(data_directory, recursive = FALSE)[contains_dirname]
    if (length(directory_path) > 0 && dir.exists(directory_path)) { return(directory_path) }
    else { print(paste0("Directory ", dirname, " not found")) ; return(NULL) }
}

directory_to_process <- find_the_directory(args[1])
if (is.null(directory_to_process)) { print("Directory to process is null. Exiting script.") ; q() }

message <- sprintf("Directory to process is %s.", directory_to_process)
print(message)
setwd(directory_to_process)

flagstat_file_path <- list.files("./qualityControl", pattern = "bamFlagstat", recursive = FALSE, full.names = TRUE)
print(flagstat_file_path[1])
print(length(unlist(strsplit(flagstat_file_path[1], "_"))))
genome_names <- unique(sapply(flagstat_file_path, function(path) {
parts <- unlist(strsplit(path, "_"))
amount_of_underscores <- length(parts)
parts[length(parts)-1]
}))
#print(genome_names)
print(paste0("Number of files to process ", length(flagstat_file_path)))

#print(sapply(bam_flagstat_filepath, function(path) {
#parts <- unlist(strsplit(path, "_"))
#parts[3]  # Assuming the genome name is always the fourth element
#}))

#print(sample_info)
#print(sample_info$short_name)
#print(seq_along(sample_info$short_name))
pattern_for_subset <- c("total.*reads", "^mapped$")

sample_info <- read.table(list.files("./documentation", pattern = "sample_info", recursive = FALSE, full.names = TRUE), sep = ",", header = TRUE)
#percent_mapped_df <- data.frame(matrix(ncol = length(genome_names), nrow = nrow(sample_info)))

percent_mapped_df <- data.frame(matrix(ncol = length(genome_names), nrow = 0))
colnames(percent_mapped_df) <- genome_names 
#print(percent_mapped_df)
#print(colnames(sample_info))
if(any(is.na(sample_info$sample_ID))) {
    for (sample in sample_info$short_name) {
    #    print(sample)
        is_flagstat <- grepl("bamFlagstat", list.files("./qualityControl", pattern = as.character(sample), recursive = FALSE, full.names = TRUE))
        bam_flagstat_filepath <- list.files("./qualityControl", pattern = as.character(sample), recursive = FALSE, full.names = TRUE)[is_flagstat]
    #    print(bam_flagstat_filepath)
        data <- list()
        for (flagstat_filepath in bam_flagstat_filepath) {
            flagstat_table <- read.table(flagstat_filepath, sep = "\t")
    #        print(lapply(pattern_for_subset, function(x) { grepl(x, flagstat_table[,3]) }))
                
            is_row_to_subset <- Reduce("|", lapply(pattern_for_subset, function(pattern) { grepl(pattern, flagstat_table[,3]) }))
            subset_table <- flagstat_table[is_row_to_subset,]
    #        print(subset_table[2,1]) ; print(subset_table[1,1])
            percent_mapped <- (as.numeric(subset_table[2,1]) / as.numeric(subset_table[1,1])) * 100
            data <- append(data, percent_mapped)
    #        print(percent_mapped)
    #        print(head(flagstat_table))
        }
#        print(as.data.frame(t(unlist(data))))
        percent_mapped_df_to_append <- as.data.frame(t(unlist(data)))
        colnames(percent_mapped_df_to_append) <- genome_names
#        print(percent_mapped_df_to_append)
        percent_mapped_df <- rbind(percent_mapped_df, percent_mapped_df_to_append)
    #    print(percent_mapped_df)
    }
} else {
    for (sample in sample_info$sample_ID) {
    #   print(sample)
        is_flagstat <- grepl("bamFlagstat", list.files("./qualityControl", pattern = as.character(sample), recursive = FALSE, full.names = TRUE))
        bam_flagstat_filepath <- list.files("./qualityControl", pattern = as.character(sample), recursive = FALSE, full.names = TRUE)[is_flagstat]
    #   print(bam_flagstat_filepath)
        data <- list()
        for (flagstat_filepath in bam_flagstat_filepath) {
            flagstat_table <- read.table(flagstat_filepath, sep = "\t")
    #       print(lapply(pattern_for_subset, function(x) { grepl(x, flagstat_table[,3]) }))

            is_row_to_subset <- Reduce("|", lapply(pattern_for_subset, function(pattern) { grepl(pattern, flagstat_table[,3]) }))
            subset_table <- flagstat_table[is_row_to_subset,]
    #       print(subset_table[2,1]) ; print(subset_table[1,1])
            percent_mapped <- (as.numeric(subset_table[2,1]) / as.numeric(subset_table[1,1])) * 100
            data <- append(data, percent_mapped)
    #       print(percent_mapped)
    #       print(head(flagstat_table))
        }
#        print(as.data.frame(t(unlist(data))))
        percent_mapped_df_to_append <- as.data.frame(t(unlist(data)))
        colnames(percent_mapped_df_to_append) <- genome_names
#        print(percent_mapped_df_to_append)
        percent_mapped_df <- rbind(percent_mapped_df, percent_mapped_df_to_append)
    #   print(percent_mapped_df)
    }
}

print(percent_mapped_df[2:nrow(percent_mapped_df),])
