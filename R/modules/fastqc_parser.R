#' FastQC Data Parsing Functions
parse_fastqc_files <- function(directory,
                             config = CONFIG$FASTQC) {
    log_info("Parsing FastQC files in:", directory)
    
    # Find FastQC files
    files <- find_fastqc_files(directory)
    
    if (length(files) == 0) {
        log_warning("No FastQC files found")
        return(NULL)
    }
    
    log_info("Found", length(files), "files to process")
    
    # Process each file
    lapply(files, process_fastqc_file, config)
}

find_fastqc_files <- function(directory) {
    qc_dir <- file.path(directory, CONFIG$PATHS$QC_SUBDIR)
    
    list.files(
        qc_dir,
        pattern = CONFIG$FASTQC$PATTERNS$DATA_FILE,
        recursive = TRUE,
        full.names = TRUE
    )
}

process_fastqc_file <- function(file_path,
                              config = CONFIG$FASTQC) {
    log_info("Processing file:", file_path)
    
    # Read file content
    lines <- readLines(file_path)
    
    # Find module boundaries
    modules <- find_module_boundaries(lines, config)
    
    # Process each module
    results <- process_modules(lines, modules, file_path)
    
    # Write summary
    write_fastqc_summary(results$summary, file_path)
}

find_module_boundaries <- function(lines, config) {
    starts <- which(grepl(config$PATTERNS$MODULE_START, lines))
    ends <- which(grepl(config$PATTERNS$MODULE_END, lines))
    
    # Remove end markers from starts
    list(
        starts = starts[!(starts %in% ends)],
        ends = ends
    )
}

process_modules <- function(lines,
                          boundaries,
                          file_path) {
    log_info("Processing modules")
    
    summaries <- list()
    
    for (i in seq_along(boundaries$starts)) {
        module_lines <- lines[
            boundaries$starts[i]:boundaries$ends[i]
        ]
        
        result <- process_single_module(
            module_lines,
            file_path
        )
        
        summaries <- append(summaries, result$summary)
    }
    
    list(
        summary = create_summary_table(summaries)
    )
}

process_single_module <- function(lines, file_path) {
    # Extract module name
    module_name <- extract_module_name(lines[1])
    
    # Find and process data section
    data <- extract_module_data(lines[-c(1, length(lines))])
    
    # Write module data
    if (!is.null(data)) {
        write_module_data(data, module_name, file_path)
    }
    
    list(
        summary = module_name,
        data = data
    )
}
