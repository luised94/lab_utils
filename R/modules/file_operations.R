load_sample_table <- function(directory_path) {
    cat("Loading sample_table from", directory_path, "\n")
    documentation_dir_path <- file.path(directory_path, "documentation")
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
    cat(sprintf("Reading %s\n", sample_table_path))
    cat("Head of sample_table\n")
    print(head(sample_table))
    return(sample_table)
}

determine_sample_id <- function(directory_path) {
    fastq_directory_path <- file.path(directory_path, "fastq")
    fastq_file_paths <- list.files(fastq_directory_path, pattern = "*.fastq", full.names = TRUE)
    if(length(fastq_file_paths) == 0) {
        cat(sprintf("No fastq files found in %s\n", fastq_directory_path))
        cat("Consult 000_setupExperimentDir to create sample_table.tsv for the project and download files from BMC\n")
        stop()
    } else {
        cat(sprintf("Found %s files in %s.\n", length(fastq_file_paths), fastq_directory_path))
    }
    fastq_file_names <- basename(fastq_file_paths)
    ID_regex <- "\\d{5,6}"
    fastq_split_string_list <- strsplit(fastq_file_names, "_|-")
    sample_IDs <- lapply(fastq_split_string_list, function(fastq_split_string_list) {
        for(split_string in fastq_split_string_list) {
            if(grepl(ID_regex, split_string)) {
                return(split_string)
            }
        }
    })
    if(!all(unlist(lapply(sample_IDs, length)) == 1)) {
        cat("At least one of the files did not extract exactly one sample ID.\n")
        cat("Files with problems:\n")
        print(fastq_file_names[unlist(lapply(sample_IDs, length)) != 1])
        cat("Verify sample names. Redownload from BMC if necessary.\n")
        cat(sprintf("Regex pattern used %s:\n", ID_regex))
        stop()
    } else {
        sample_IDs <- unlist(sample_IDs)
        cat(sprintf("Found %s sample_IDs.\n", length(sample_IDs)))
        cat(sprintf("First sample_ID: %s\n",sample_IDs[1]))
        cat(sprintf("Last sample_ID: %s\n", sample_IDs[length(sample_IDs)]))
        cat("Returning sample_ID array.\n")
        return(sample_IDs)
    }
}
