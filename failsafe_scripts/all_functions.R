#' @title Validate Comparison Input Parameters
#' @description Validates input parameters for comparison analysis
#' @param metadata data.frame Experiment metadata
#' @param comparisons list List of comparison expressions
#' @return logical TRUE if validation passes
#' @throws ValidationError if inputs don't meet requirements
#' @examples
#' validate_comparison_inputs(mtcars, list(high_mpg = quote(mpg > 20)))
#' @importFrom methods is
validate_comparison_inputs <- function(metadata, comparisons) {
    validation_result <- tryCatch({
        stopifnot(
            "metadata must be a data.frame" = is.data.frame(metadata),
            "comparisons must be a list" = is.list(comparisons),
            "comparisons must be named" = !is.null(names(comparisons)),
            "metadata must have rows" = nrow(metadata) > 0,
            "comparisons must not be empty" = length(comparisons) > 0
        )
        
        all_expressions <- vapply(comparisons, is.language, logical(1))
        if (!all(all_expressions)) {
            stop("All comparisons must be language expressions")
        }
        
        list(
            valid = TRUE,
            error = NULL
        )
    }, error = function(e) {
        list(
            valid = FALSE,
            error = e$message
        )
    })
    
    return(validation_result)
}

#' @title Execute Single Comparison
#' @description Executes a single comparison expression against metadata
#' @param metadata data.frame Experiment metadata
#' @param comparison language Quoted comparison expression
#' @return logical vector Rows matching comparison
#' @throws EvalError if comparison evaluation fails
#' @precondition metadata must contain all variables referenced in comparison
#' @postcondition return value length equals nrow(metadata)
#' @examples
#' execute_comparison(mtcars, quote(mpg > 20))
#' @seealso validate_comparison_inputs
execute_comparison <- function(metadata, comparison) {
    result <- tryCatch({
        stopifnot(
            "metadata must be a data.frame" = is.data.frame(metadata),
            "comparison must be a language object" = is.language(comparison)
        )
        
        match_result <- eval(comparison, metadata, parent.frame())
        stopifnot(
            "comparison must return logical vector" = is.logical(match_result),
            "comparison result length must match metadata rows" = length(match_result) == nrow(metadata)
        )
        
        list(
            success = TRUE,
            matches = match_result,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            matches = logical(nrow(metadata)),
            error = e$message
        )
    })
    
    return(result)
}

#' @title Create Sample Summary
#' @description Generates a formatted summary string for a sample row
#' @param sample_row data.frame Single row of sample data
#' @param required_cols character Required column names
#' @return character Formatted summary string
#' @throws MissingColumnError if required columns are absent
#' @examples
#' create_sample_summary(mtcars[1,], c("mpg", "cyl"))
#' @seealso execute_comparison
create_sample_summary <- function(sample_row, required_cols) {
    summary_result <- tryCatch({
        stopifnot(
            "sample_row must be a data.frame" = is.data.frame(sample_row),
            "required_cols must be character vector" = is.character(required_cols),
            "sample_row must have exactly one row" = nrow(sample_row) == 1,
            "required columns must exist" = all(required_cols %in% names(sample_row))
        )
        
        values <- vapply(required_cols, function(col) {
            as.character(sample_row[[col]])
        }, character(1))
        
        summary <- paste(names(values), values, sep = ": ", collapse = ", ")
        
        list(
            success = TRUE,
            summary = summary,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            summary = "",
            error = e$message
        )
    })
    
    return(summary_result)
}

#' @title Process Comparison Results
#' @description Creates structured results from comparison matches
#' @param metadata data.frame Original metadata
#' @param match_rows logical Vector of matching rows
#' @param comparison_name character Name of comparison
#' @return list Structured comparison results
#' @postcondition return value contains original row indices
#' @examples
#' process_comparison_results(mtcars, mtcars$mpg > 20, "high_mpg")
#' @seealso create_sample_summary
process_comparison_results <- function(metadata, match_rows, comparison_name) {
    result <- tryCatch({
        stopifnot(
            "metadata must be a data.frame" = is.data.frame(metadata),
            "match_rows must be logical vector" = is.logical(match_rows),
            "comparison_name must be character" = is.character(comparison_name),
            "match_rows length must match metadata rows" = length(match_rows) == nrow(metadata)
        )
        
        matching_indices <- which(match_rows)
        matching_data <- metadata[match_rows, , drop = FALSE]
        
        list(
            success = TRUE,
            results = list(
                name = comparison_name,
                matching_indices = matching_indices,
                matching_data = matching_data,
                match_count = length(matching_indices)
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            results = NULL,
            error = e$message
        )
    })
    
    return(result)
}
#' @title Analyze Comparisons
#' @description Main function orchestrating comparison analysis workflow
#' @param metadata data.frame Experiment metadata
#' @param comparisons list Named list of comparison expressions
#' @return list Results for each comparison
#' @importFrom purrr map safely
#' @examples
#' analyze_comparisons(
#'   mtcars,
#'   list(
#'     high_mpg = quote(mpg > 20),
#'     high_cyl = quote(cyl > 6)
#'   )
#' )
#' @seealso
#'   validate_comparison_inputs,
#'   execute_comparison,
#'   process_comparison_results
analyze_comparisons <- function(metadata, comparisons) {
    # Validate inputs
    validation <- validate_comparison_inputs(metadata, comparisons)
    if (!validation$valid) {
        return(list(
            success = FALSE,
            results = NULL,
            error = validation$error
        ))
    }
    
    # Process each comparison
    results <- lapply(names(comparisons), function(comp_name) {
        comparison <- comparisons[[comp_name]]
        
        # Execute comparison
        exec_result <- execute_comparison(metadata, comparison)
        if (!exec_result$success) {
            return(list(
                name = comp_name,
                success = FALSE,
                error = exec_result$error
            ))
        }
        
        # Process results
        proc_result <- process_comparison_results(
            metadata,
            exec_result$matches,
            comp_name
        )
        
        if (!proc_result$success) {
            return(list(
                name = comp_name,
                success = FALSE,
                error = proc_result$error
            ))
        }
        
        proc_result$results
    })
    
    names(results) <- names(comparisons)
    
    return(list(
        success = TRUE,
        results = results,
        error = NULL
    ))
}

#' @title Validate Category Names
#' @description Validates category names against table columns
#' @param table data.frame The input data frame
#' @param categories character Vector of category names
#' @return list Validation result with status and validated categories
#' @throws CategoryValidationError if categories are invalid
#' @precondition table must have at least one column
#' @postcondition returned categories contain 'antibody'
#' @examples
#' validate_category_names(data.frame(antibody = 1, type = 2), c("type"))
validate_category_names <- function(table, categories) {
    result <- tryCatch({
        stopifnot(
            "table must be a data.frame" = is.data.frame(table),
            "categories must be character vector" = is.character(categories),
            "table must have columns" = ncol(table) > 0,
            "categories cannot be empty" = length(categories) > 0
        )
        
        missing_cats <- setdiff(categories, colnames(table))
        if (length(missing_cats) > 0) {
            stop(sprintf("Missing categories: %s", 
                        paste(missing_cats, collapse = ", ")))
        }
        
        list(
            success = TRUE,
            data = categories,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}
#' @title Ensure Required Categories
#' @description Ensures required categories are present
#' @param categories character Vector of category names
#' @param required character Vector of required category names
#' @return character Updated category vector
#' @precondition categories must not be empty
#' @examples
#' ensure_required_categories(c("type"), "antibody")
#' @seealso validate_category_names
ensure_required_categories <- function(categories, required) {
    result <- tryCatch({
        stopifnot(
            "categories must be character vector" = is.character(categories),
            "required must be character vector" = is.character(required)
        )
        
        unique_cats <- unique(c(required, categories))
        
        list(
            success = TRUE,
            data = unique_cats,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Extract Category Values
#' @description Extracts unique values for each category
#' @param table data.frame Input data frame
#' @param categories character Validated category names
#' @return list Named list of unique values per category
#' @throws CategoryAccessError if categories cannot be accessed
#' @examples
#' extract_category_values(data.frame(antibody = c(1,1,2)), "antibody")
#' @seealso validate_category_names
extract_category_values <- function(table, categories) {
    result <- tryCatch({
        stopifnot(
            "table must be a data.frame" = is.data.frame(table),
            "categories must be character vector" = is.character(categories),
            "categories must exist in table" = all(categories %in% colnames(table))
        )
        
        unique_values <- lapply(table[categories], unique)
        
        list(
            success = TRUE,
            data = unique_values,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Filter Relevant Categories
#' @description Filters categories based on value uniqueness
#' @param sample data.frame Single sample row
#' @param categories character Category names
#' @param unique_values list Unique values per category
#' @return character Vector of relevant category values
#' @examples
#' filter_relevant_categories(
#'   data.frame(antibody = 1, type = "A"),
#'   c("antibody", "type"),
#'   list(antibody = 1:2, type = "A")
#' )
filter_relevant_categories <- function(sample, categories, unique_values) {
    result <- tryCatch({
        stopifnot(
            "sample must be a data.frame" = is.data.frame(sample),
            "categories must be character vector" = is.character(categories),
            "unique_values must be a list" = is.list(unique_values),
            "sample must have one row" = nrow(sample) == 1
        )
        
        relevant_values <- vapply(categories, function(cat) {
            if (length(unique_values[[cat]]) > 1 || cat == "antibody") {
                return(as.character(sample[[cat]]))
            }
            return(NA_character_)
        }, character(1))
        
        relevant_values <- relevant_values[!is.na(relevant_values)]
        
        list(
            success = TRUE,
            data = relevant_values,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}
#' @title Format Sample Label
#' @description Creates formatted label from category values
#' @param category_values character Vector of category values
#' @param separator character Label component separator
#' @return character Formatted sample label
#' @examples
#' format_sample_label(c("AB1", "TypeA"), "_")
#' @seealso filter_relevant_categories
format_sample_label <- function(category_values, separator = "_") {
    result <- tryCatch({
        stopifnot(
            "category_values must be character vector" = is.character(category_values),
            "separator must be character" = is.character(separator),
            "separator must be length 1" = length(separator) == 1
        )
        
        label <- paste(category_values, collapse = separator)
        
        list(
            success = TRUE,
            data = label,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}
#' @title Create Sample Labels
#' @description Generates labels for all samples
#' @param table data.frame Input data frame
#' @param categories character Category names
#' @param unique_values list Unique values per category
#' @return character Vector of sample labels
#' @examples
#' create_sample_labels(
#'   data.frame(antibody = c(1,2), type = c("A","B")),
#'   c("antibody", "type")
#' )
#' @seealso
#'   filter_relevant_categories,
#'   format_sample_label
create_sample_labels <- function(table, categories, unique_values) {
    result <- tryCatch({
        stopifnot(
            "table must be a data.frame" = is.data.frame(table),
            "categories must be character vector" = is.character(categories),
            "unique_values must be a list" = is.list(unique_values)
        )
        
        labels <- vapply(seq_len(nrow(table)), function(i) {
            sample_row <- table[i, , drop = FALSE]
            relevant_cats <- filter_relevant_categories(
                sample_row, categories, unique_values
            )
            if (!relevant_cats$success) {
                return(NA_character_)
            }
            
            label_result <- format_sample_label(relevant_cats$data)
            if (!label_result$success) {
                return(NA_character_)
            }
            
            return(label_result$data)
        }, character(1))
        
        list(
            success = TRUE,
            data = labels,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Generate Sample Labels
#' @description Main function orchestrating label generation
#' @param table data.frame Input data frame
#' @param categories character Desired category names
#' @param options list Optional configuration parameters
#' @return list Results containing labels and metadata
#' @throws LabelGenerationError if label creation fails
#' @examples
#' generate_sample_labels(
#'   data.frame(antibody = 1:2, type = c("A","B")),
#'   c("type"),
#'   list(separator = "_", verbose = TRUE)
#' )
#' @seealso
#'   validate_category_names,
#'   create_sample_labels
generate_sample_labels <- function(table, categories, 
                                 options = list(separator = "_", 
                                              verbose = FALSE)) {
    # Validate inputs
    valid_cats <- validate_category_names(table, categories)
    if (!valid_cats$success) {
        return(valid_cats)
    }
    
    # Ensure required categories
    updated_cats <- ensure_required_categories(valid_cats$data, "antibody")
    if (!updated_cats$success) {
        return(updated_cats)
    }
    
    # Extract unique values
    unique_vals <- extract_category_values(table, updated_cats$data)
    if (!unique_vals$success) {
        return(unique_vals)
    }
    
    # Generate labels
    labels <- create_sample_labels(table, updated_cats$data, 
                                 unique_vals$data)
    
    if (options$verbose && labels$success) {
        message("Generated ", length(labels$data), " labels")
        message("Categories used: ", 
                paste(updated_cats$data, collapse = ", "))
    }
    
    return(labels)
}

#' @title Validate Directory Structure
#' @description Validates experiment directory structure and required files
#' @param directory_path character Path to experiment directory
#' @return list Validation result {success, data, error}
#' @throws DirectoryValidationError if structure invalid
#' @precondition directory must exist
#' @postcondition confirms presence of fastq and documentation subdirectories
#' @examples
#' directory_structure_validate("path/to/experiment")
#' @seealso experiment_files_scan
directory_structure_validate <- function(directory_path) {
    result <- tryCatch({
        stopifnot(
            "directory_path must be character" = is.character(directory_path),
            "directory_path must have length 1" = length(directory_path) == 1,
            "directory_path must not be empty" = nzchar(directory_path)
        )
        
        required_dirs <- c("fastq", "documentation")
        dir_exists <- dir.exists(file.path(directory_path, required_dirs))
        
        if (!all(dir_exists)) {
            missing_dirs <- required_dirs[!dir_exists]
            stop(sprintf("Missing directories: %s", 
                        paste(missing_dirs, collapse = ", ")))
        }
        
        list(
            success = TRUE,
            data = list(
                path = directory_path,
                subdirectories = required_dirs
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Scan Experiment Files
#' @description Scans directory for experiment files matching pattern
#' @param directory_path character Path to experiment directory
#' @param file_pattern character Regex pattern for matching files
#' @return list File scan result {success, data, error}
#' @throws FileScanError if file access fails
#' @importFrom base list.files
#' @examples
#' experiment_files_scan("path/to/experiment", "consolidated_.*_sequence\\.fastq$")
#' @seealso directory_structure_validate
experiment_files_scan <- function(directory_path, file_pattern) {
    result <- tryCatch({
        stopifnot(
            "directory_path must be character" = is.character(directory_path),
            "file_pattern must be character" = is.character(file_pattern),
            "file_pattern must have length 1" = length(file_pattern) == 1
        )
        
        files <- list.files(
            path = file.path(directory_path, "fastq"),
            pattern = file_pattern,
            full.names = FALSE
        )
        
        if (length(files) == 0) {
            stop("No matching files found")
        }
        
        list(
            success = TRUE,
            data = files,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Extract Experiment Identifiers
#' @description Extracts experiment IDs from file names
#' @param file_names character Vector of file names
#' @param pattern_extract character Regex pattern for ID extraction
#' @return list Extraction result {success, data, error}
#' @throws ExtractionError if pattern matching fails
#' @examples
#' experiment_identifiers_extract(files, "consolidated_([0-9]{5,6})_sequence\\.fastq")
#' @seealso experiment_files_scan
experiment_identifiers_extract <- function(file_names, pattern_extract) {
    result <- tryCatch({
        stopifnot(
            "file_names must be character vector" = is.character(file_names),
            "pattern_extract must be character" = is.character(pattern_extract),
            "pattern_extract must have length 1" = length(pattern_extract) == 1
        )
        
        identifiers <- gsub(pattern_extract, "\\1", file_names)
        
        if (any(identifiers == file_names)) {
            stop("Pattern extraction failed for some files")
        }
        
        list(
            success = TRUE,
            data = identifiers,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Validate Metadata Schema
#' @description Validates metadata file against expected schema
#' @param metadata_frame data.frame Metadata table
#' @param expected_columns character Required column names
#' @return list Validation result {success, data, error}
#' @throws SchemaValidationError if columns missing
#' @examples
#' metadata_schema_validate(metadata, c("sample_id", "condition"))
metadata_schema_validate <- function(metadata_frame, expected_columns) {
    result <- tryCatch({
        stopifnot(
            "metadata_frame must be data.frame" = is.data.frame(metadata_frame),
            "expected_columns must be character" = is.character(expected_columns)
        )
        
        missing_cols <- setdiff(expected_columns, colnames(metadata_frame))
        if (length(missing_cols) > 0) {
            stop(sprintf("Missing columns: %s", 
                        paste(missing_cols, collapse = ", ")))
        }
        
        list(
            success = TRUE,
            data = colnames(metadata_frame),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Validate Factor Categories
#' @description Validates factor levels against allowed categories
#' @param column_data character Vector of column values
#' @param allowed_levels character Vector of valid levels
#' @return list Validation result {success, data, error}
#' @throws CategoryValidationError if invalid levels found
#' @examples
#' factor_categories_validate(data$treatment, c("control", "treated"))
#' @seealso factor_categories_enforce
factor_categories_validate <- function(column_data, allowed_levels) {
    result <- tryCatch({
        stopifnot(
            "column_data must be atomic vector" = is.atomic(column_data),
            "allowed_levels must be character" = is.character(allowed_levels)
        )
        
        invalid_levels <- setdiff(unique(column_data), allowed_levels)
        if (length(invalid_levels) > 0) {
            stop(sprintf("Invalid levels found: %s", 
                        paste(invalid_levels, collapse = ", ")))
        }
        
        list(
            success = TRUE,
            data = allowed_levels,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Enforce Factor Categories
#' @description Converts columns to factors with specified levels
#' @param metadata_frame data.frame Metadata table
#' @param category_definitions list Named list of allowed levels
#' @return list Factor enforcement result {success, data, error}
#' @throws FactorEnforcementError if conversion fails
#' @examples
#' factor_categories_enforce(metadata, list(treatment = c("control", "treated")))
#' @seealso factor_categories_validate
factor_categories_enforce <- function(metadata_frame, category_definitions) {
    result <- tryCatch({
        stopifnot(
            "metadata_frame must be data.frame" = is.data.frame(metadata_frame),
            "category_definitions must be list" = is.list(category_definitions)
        )
        
        for (col_name in names(category_definitions)) {
            validation <- factor_categories_validate(
                metadata_frame[[col_name]],
                category_definitions[[col_name]]
            )
            
            if (!validation$success) {
                stop(sprintf("Validation failed for column %s: %s", 
                            col_name, validation$error))
            }
            
            metadata_frame[[col_name]] <- factor(
                metadata_frame[[col_name]],
                levels = category_definitions[[col_name]],
                ordered = TRUE
            )
        }
        
        list(
            success = TRUE,
            data = metadata_frame,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Sort Metadata Frame
#' @description Sorts metadata by specified column order
#' @param metadata_frame data.frame Metadata table
#' @param sort_columns character Vector of column names
#' @return list Sorting result {success, data, error}
#' @throws SortError if columns missing
#' @examples
#' metadata_frame_sort(metadata, c("sample_id", "condition"))
metadata_frame_sort <- function(metadata_frame, sort_columns) {
    result <- tryCatch({
        stopifnot(
            "metadata_frame must be data.frame" = is.data.frame(metadata_frame),
            "sort_columns must be character" = is.character(sort_columns)
        )
        
        schema_valid <- metadata_schema_validate(metadata_frame, sort_columns)
        if (!schema_valid$success) {
            stop(schema_valid$error)
        }
        
        sorted_frame <- metadata_frame[do.call(order, 
                                             metadata_frame[sort_columns]), ]
        
        list(
            success = TRUE,
            data = sorted_frame,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Process Experiment Metadata
#' @description Main function orchestrating metadata processing
#' @param directory_path character Path to experiment directory
#' @param configuration list Processing configuration
#' @param output_options list Output configuration
#' @return list Processing result {success, data, error}
#' @throws MetadataProcessingError if processing fails
#' @importFrom utils read.csv write.csv
#' @examples
#' experiment_metadata_process(
#'   "path/to/experiment",
#'   list(categories = list(...), column_order = c(...)),
#'   list(output_file = TRUE, output_path = "...")
#' )
#' @seealso
#'   directory_structure_validate,
#'   experiment_files_scan,
#'   factor_categories_enforce,
#'   metadata_frame_sort
experiment_metadata_process <- function(directory_path, configuration, 
                                      output_options) {
    result <- tryCatch({
        # Validate directory structure
        dir_valid <- directory_structure_validate(directory_path)
        if (!dir_valid$success) return(dir_valid)
        
        # Scan for files
        files <- experiment_files_scan(
            directory_path,
            "consolidated_.*_sequence\\.fastq$"
        )
        if (!files$success) return(files)
        
        # Extract identifiers
        ids <- experiment_identifiers_extract(
            files$data,
            "consolidated_([0-9]{5,6})_sequence\\.fastq"
        )
        if (!ids$success) return(ids)
        
        # Read metadata
        metadata_path <- file.path(directory_path, "documentation", 
                                 paste0(basename(directory_path), 
                                      "_sample_grid.csv"))
        metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
        
        # Process metadata
        processed <- factor_categories_enforce(metadata, 
                                            configuration$categories)
        if (!processed$success) return(processed)
        
        # Sort metadata
        sorted <- metadata_frame_sort(processed$data, 
                                    configuration$column_order)
        if (!sorted$success) return(sorted)
        
        # Add identifiers
        sorted$data$sample_id <- ids$data
        
        # Save if requested
        if (output_options$output_file) {
            write.csv(sorted$data, output_options$output_path, 
                     row.names = FALSE)
        }
        
        list(
            success = TRUE,
            data = sorted$data,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Validate Metadata File Path
#' @description Validates existence and structure of metadata file path
#' @param directory_path character Path to experiment directory
#' @param file_name_pattern character Pattern for metadata filename
#' @return list Validation result {success, data, error}
#' @throws PathValidationError if file not found
#' @precondition directory must exist
#' @postcondition returns valid file path
#' @examples
#' metadata_path_validate("experiment/dir", "%s_processed_grid.csv")
#' @seealso metadata_file_read
metadata_path_validate <- function(directory_path, file_name_pattern) {
    result <- tryCatch({
        stopifnot(
            "directory_path must be character" = is.character(directory_path),
            "file_name_pattern must be character" = is.character(file_name_pattern),
            "directory_path must not be empty" = nzchar(directory_path),
            "file_name_pattern must not be empty" = nzchar(file_name_pattern)
        )
        
        file_path <- file.path(
            directory_path,
            "documentation",
            sprintf(file_name_pattern, basename(directory_path))
        )
        
        if (!file.exists(file_path)) {
            stop(sprintf("File not found: %s", file_path))
        }
        
        list(
            success = TRUE,
            data = file_path,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Read Metadata File
#' @description Reads and validates metadata CSV file
#' @param file_path character Full path to metadata file
#' @param read_options list CSV reading options
#' @return list Read result {success, data, error}
#' @throws FileReadError if read operation fails
#' @importFrom utils read.csv
#' @examples
#' metadata_file_read("path/to/metadata.csv", list(stringsAsFactors = FALSE))
#' @seealso metadata_path_validate
metadata_file_read <- function(file_path, read_options) {
    result <- tryCatch({
        stopifnot(
            "file_path must be character" = is.character(file_path),
            "read_options must be list" = is.list(read_options)
        )
        
        read_args <- c(list(file = file_path), read_options)
        metadata <- do.call(read.csv, read_args)
        
        if (nrow(metadata) == 0) {
            stop("Metadata file is empty")
        }
        
        list(
            success = TRUE,
            data = metadata,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Validate Comparison Expression
#' @description Validates syntax and structure of comparison expression
#' @param comparison_expression language Quoted comparison expression
#' @param metadata_columns character Available column names
#' @return list Validation result {success, data, error}
#' @throws ExpressionValidationError if expression invalid
#' @examples
#' comparison_expression_validate(quote(treatment == "control"), c("treatment"))
#' @seealso comparison_expression_evaluate
comparison_expression_validate <- function(comparison_expression, metadata_columns) {
    result <- tryCatch({
        stopifnot(
            "comparison_expression must be language" = is.language(comparison_expression),
            "metadata_columns must be character" = is.character(metadata_columns)
        )
        
        expr_vars <- all.vars(comparison_expression)
        missing_vars <- setdiff(expr_vars, metadata_columns)
        
        if (length(missing_vars) > 0) {
            stop(sprintf("Unknown variables in expression: %s", 
                        paste(missing_vars, collapse = ", ")))
        }
        
        list(
            success = TRUE,
            data = comparison_expression,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Evaluate Comparison Expression
#' @description Evaluates comparison expression against metadata
#' @param comparison_expression language Validated comparison expression
#' @param metadata_frame data.frame Metadata to evaluate against
#' @return list Evaluation result {success, data, error}
#' @throws ExpressionEvaluationError if evaluation fails
#' @examples
#' comparison_expression_evaluate(quote(treatment == "control"), metadata)
#' @seealso comparison_expression_validate
comparison_expression_evaluate <- function(comparison_expression, metadata_frame) {
    result <- tryCatch({
        stopifnot(
            "comparison_expression must be language" = is.language(comparison_expression),
            "metadata_frame must be data.frame" = is.data.frame(metadata_frame)
        )
        
        matches <- eval(comparison_expression, metadata_frame)
        
        if (!is.logical(matches)) {
            stop("Expression must evaluate to logical vector")
        }
        
        if (length(matches) != nrow(metadata_frame)) {
            stop("Expression result length mismatch")
        }
        
        list(
            success = TRUE,
            data = matches,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Format Sample Details
#' @description Creates formatted string of sample information
#' @param sample_row data.frame Single sample metadata row
#' @param format_template character Template for sample formatting
#' @return list Formatting result {success, data, error}
#' @throws FormattingError if required fields missing
#' @examples
#' sample_details_format(sample, "%s_%s_%s")
#' @seealso comparison_results_format
sample_details_format <- function(sample_row, format_template) {
    result <- tryCatch({
        stopifnot(
            "sample_row must be data.frame" = is.data.frame(sample_row),
            "format_template must be character" = is.character(format_template),
            "sample_row must have one row" = nrow(sample_row) == 1
        )
        
        required_fields <- str_count(format_template, "%s")
        if (ncol(sample_row) < required_fields) {
            stop("Insufficient columns for format template")
        }
        
        formatted <- do.call(
            sprintf,
            c(list(format_template), as.list(sample_row[1, 1:required_fields]))
        )
        
        list(
            success = TRUE,
            data = formatted,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Format Comparison Results
#' @description Formats complete comparison results
#' @param matched_samples data.frame Matching sample data
#' @param comparison_name character Name of comparison
#' @return list Formatting result {success, data, error}
#' @throws ResultFormattingError if formatting fails
#' @examples
#' comparison_results_format(matches, "control_samples")
#' @seealso sample_details_format
comparison_results_format <- function(matched_samples, comparison_name) {
    result <- tryCatch({
        stopifnot(
            "matched_samples must be data.frame" = is.data.frame(matched_samples),
            "comparison_name must be character" = is.character(comparison_name)
        )
        
        formatted_samples <- lapply(seq_len(nrow(matched_samples)), function(i) {
            sample_details_format(
                matched_samples[i, , drop = FALSE],
                "Sample %d: %s_%s_%s_%s (Exp: %s)"
            )
        })
        
        success <- all(vapply(formatted_samples, `[[`, logical(1), "success"))
        if (!success) {
            errors <- vapply(formatted_samples, `[[`, character(1), "error")
            stop(paste(unique(errors), collapse = "; "))
        }
        
        list(
            success = TRUE,
            data = list(
                name = comparison_name,
                count = nrow(matched_samples),
                samples = vapply(formatted_samples, `[[`, character(1), "data")
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Process Metadata Comparisons
#' @description Main function orchestrating metadata comparison analysis
#' @param metadata_frame data.frame Processed metadata
#' @param comparison_definitions list Named list of comparison expressions
#' @param output_options list Output formatting options
#' @return list Analysis results {success, data, error}
#' @throws ComparisonProcessingError if analysis fails
#' @examples
#' metadata_comparisons_process(
#'   metadata,
#'   list(control = quote(treatment == "control")),
#'   list(verbose = TRUE)
#' )
#' @seealso
#'   comparison_expression_validate,
#'   comparison_expression_evaluate,
#'   comparison_results_format
metadata_comparisons_process <- function(metadata_frame, comparison_definitions, 
                                       output_options) {
    result <- tryCatch({
        stopifnot(
            "metadata_frame must be data.frame" = is.data.frame(metadata_frame),
            "comparison_definitions must be list" = is.list(comparison_definitions),
            "output_options must be list" = is.list(output_options)
        )
        
        results <- lapply(names(comparison_definitions), function(comp_name) {
            # Validate expression
            valid_expr <- comparison_expression_validate(
                comparison_definitions[[comp_name]],
                colnames(metadata_frame)
            )
            if (!valid_expr$success) {
                return(valid_expr)
            }
            
            # Evaluate expression
            matches <- comparison_expression_evaluate(
                valid_expr$data,
                metadata_frame
            )
            if (!matches$success) {
                return(matches)
            }
            
            # Format results
            matched_data <- metadata_frame[matches$data, , drop = FALSE]
            formatted <- comparison_results_format(matched_data, comp_name)
            
            if (output_options$verbose && formatted$success) {
                message(sprintf("\nComparison: %s", comp_name))
                message(sprintf("Found %d matches", formatted$data$count))
                if (formatted$data$count > 0) {
                    message("Samples:")
                    message(paste(formatted$data$samples, collapse = "\n"))
                }
            }
            
            formatted
        })
        
        names(results) <- names(comparison_definitions)
        
        list(
            success = TRUE,
            data = results,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Initialize Chromosome Style Mappings
#' @description Creates bidirectional mappings between chromosome naming styles
#' @return list Mapping configuration {success, data, error}
#' @throws ConfigurationError if mapping initialization fails
#' @examples
#' chromosome_mappings_initialize()
#' @seealso chromosome_style_validate
chromosome_mappings_initialize <- function() {
    result <- tryCatch({
        numeric_to_roman <- c(
            "1" = "I", "2" = "II", "3" = "III", "4" = "IV",
            "5" = "V", "6" = "VI", "7" = "VII", "8" = "VIII",
            "9" = "IX", "10" = "X", "11" = "XI", "12" = "XII",
            "13" = "XIII", "14" = "XIV", "15" = "XV", "16" = "XVI"
        )
        
        roman_to_numeric <- setNames(
            names(numeric_to_roman),
            numeric_to_roman
        )
        
        list(
            success = TRUE,
            data = list(
                numeric_to_roman = numeric_to_roman,
                roman_to_numeric = roman_to_numeric
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Validate Chromosome Style Configuration
#' @description Validates chromosome style mapping configuration
#' @param mapping_config list Chromosome mapping configuration
#' @return list Validation result {success, data, error}
#' @throws ValidationError if configuration invalid
#' @precondition mapping must contain required styles
#' @postcondition confirms bidirectional mapping integrity
#' @examples
#' chromosome_style_validate(mapping_config)
#' @seealso chromosome_mappings_initialize
chromosome_style_validate <- function(mapping_config) {
    result <- tryCatch({
        stopifnot(
            "mapping_config must be a list" = is.list(mapping_config),
            "mapping_config must contain numeric_to_roman" = !is.null(mapping_config$numeric_to_roman),
            "mapping_config must contain roman_to_numeric" = !is.null(mapping_config$roman_to_numeric)
        )
        
        # Validate bidirectional mapping
        forward_names <- names(mapping_config$numeric_to_roman)
        reverse_names <- mapping_config$roman_to_numeric[mapping_config$numeric_to_roman]
        
        if (!identical(forward_names, reverse_names)) {
            stop("Bidirectional mapping integrity violated")
        }
        
        list(
            success = TRUE,
            data = mapping_config,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Detect Chromosome Naming Style
#' @description Identifies chromosome naming convention from input
#' @param chromosome_names character Vector of chromosome names
#' @param style_patterns list Regular expressions for style detection
#' @return list Style detection result {success, data, error}
#' @throws DetectionError if style cannot be determined
#' @precondition chromosome names must be non-empty
#' @examples
#' chromosome_style_detect(c("chr1", "chr2"), style_patterns)
#' @seealso chromosome_style_validate
chromosome_style_detect <- function(chromosome_names, style_patterns) {
    result <- tryCatch({
        stopifnot(
            "chromosome_names must be character" = is.character(chromosome_names),
            "chromosome_names cannot be empty" = length(chromosome_names) > 0,
            "style_patterns must be list" = is.list(style_patterns)
        )
        
        style_matches <- vapply(style_patterns, function(pattern) {
            all(grepl(pattern, chromosome_names))
        }, logical(1))
        
        detected_style <- names(style_matches)[style_matches][1]
        if (is.na(detected_style)) detected_style <- "Unknown"
        
        list(
            success = TRUE,
            data = detected_style,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Clean Chromosome Names
#' @description Standardizes chromosome name format
#' @param chromosome_names character Vector of chromosome names
#' @return list Cleaning result {success, data, error}
#' @throws CleaningError if standardization fails
#' @examples
#' chromosome_names_clean(c("chr1", "chrII"))
#' @seealso chromosome_style_detect
chromosome_names_clean <- function(chromosome_names) {
    result <- tryCatch({
        stopifnot(
            "chromosome_names must be character" = is.character(chromosome_names)
        )
        
        clean_names <- gsub("^chr", "", chromosome_names)
        
        list(
            success = TRUE,
            data = clean_names,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Convert to UCSC Style
#' @description Converts chromosome names to UCSC format
#' @param chromosome_names character Cleaned chromosome names
#' @param mapping_config list Chromosome mapping configuration
#' @return list Conversion result {success, data, error}
#' @throws ConversionError if UCSC conversion fails
#' @examples
#' chromosome_names_to_ucsc(c("1", "2"), mapping_config)
#' @seealso chromosome_names_clean
chromosome_names_to_ucsc <- function(chromosome_names, mapping_config) {
    result <- tryCatch({
        stopifnot(
            "chromosome_names must be character" = is.character(chromosome_names),
            "mapping_config must be list" = is.list(mapping_config)
        )
        
        ucsc_names <- paste0("chr", chromosome_names)
        
        list(
            success = TRUE,
            data = ucsc_names,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Convert to Roman Style
#' @description Converts chromosome names to Roman numeral format
#' @param chromosome_names character Cleaned chromosome names
#' @param mapping_config list Chromosome mapping configuration
#' @return list Conversion result {success, data, error}
#' @throws ConversionError if Roman conversion fails
#' @examples
#' chromosome_names_to_roman(c("1", "2"), mapping_config)
#' @seealso chromosome_names_clean
chromosome_names_to_roman <- function(chromosome_names, mapping_config) {
    result <- tryCatch({
        stopifnot(
            "chromosome_names must be character" = is.character(chromosome_names),
            "mapping_config must be list" = is.list(mapping_config)
        )
        
        roman_names <- vapply(chromosome_names, function(x) {
            if (x %in% names(mapping_config$numeric_to_roman)) {
                mapping_config$numeric_to_roman[x]
            } else {
                x
            }
        }, character(1))
        
        list(
            success = TRUE,
            data = roman_names,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Convert to Numeric Style
#' @description Converts chromosome names to numeric format
#' @param chromosome_names character Cleaned chromosome names
#' @param mapping_config list Chromosome mapping configuration
#' @return list Conversion result {success, data, error}
#' @throws ConversionError if numeric conversion fails
#' @examples
#' chromosome_names_to_numeric(c("I", "II"), mapping_config)
#' @seealso chromosome_names_clean
chromosome_names_to_numeric <- function(chromosome_names, mapping_config) {
    result <- tryCatch({
        stopifnot(
            "chromosome_names must be character" = is.character(chromosome_names),
            "mapping_config must be list" = is.list(mapping_config)
        )
        
        numeric_names <- vapply(chromosome_names, function(x) {
            if (x %in% names(mapping_config$roman_to_numeric)) {
                mapping_config$roman_to_numeric[x]
            } else {
                x
            }
        }, character(1))
        
        list(
            success = TRUE,
            data = numeric_names,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Convert Chromosome Names
#' @description Main function orchestrating chromosome name conversion
#' @param chromosome_names character Vector of chromosome names
#' @param target_style character Desired naming convention
#' @param options list Conversion options and configuration
#' @return list Conversion result {success, data, error}
#' @throws ConversionError if process fails
#' @examples
#' chromosome_names_convert(
#'   c("chr1", "chr2"),
#'   "Roman",
#'   list(validate = TRUE)
#' )
#' @seealso
#'   chromosome_style_detect,
#'   chromosome_names_clean,
#'   chromosome_names_to_ucsc,
#'   chromosome_names_to_roman,
#'   chromosome_names_to_numeric
chromosome_names_convert <- function(chromosome_names, target_style, options = list()) {
    result <- tryCatch({
        stopifnot(
            "chromosome_names must be character" = is.character(chromosome_names),
            "target_style must be character" = is.character(target_style),
            "target_style must be valid" = target_style %in% c("UCSC", "Roman", "Numeric"),
            "options must be list" = is.list(options)
        )
        
        # Initialize mappings
        mappings <- chromosome_mappings_initialize()
        if (!mappings$success) return(mappings)
        
        # Clean names
        clean_result <- chromosome_names_clean(chromosome_names)
        if (!clean_result$success) return(clean_result)
        
        # Convert based on target style
        converted <- switch(target_style,
            "UCSC" = chromosome_names_to_ucsc(clean_result$data, mappings$data),
            "Roman" = chromosome_names_to_roman(clean_result$data, mappings$data),
            "Numeric" = chromosome_names_to_numeric(clean_result$data, mappings$data)
        )
        
        if (!converted$success) return(converted)
        
        list(
            success = TRUE,
            data = converted$data,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Validate Bigwig File Existence
#' @description Validates existence of bigwig file path
#' @param file_path character Path to bigwig file
#' @return list Validation result {success, data, error}
#' @throws FileValidationError if file not found
#' @precondition file path must be non-empty
#' @examples
#' bigwig_file_exists("path/to/file.bw")
#' @seealso bigwig_file_validate
bigwig_file_exists <- function(file_path) {
    result <- tryCatch({
        stopifnot(
            "file_path must be character" = is.character(file_path),
            "file_path must have length 1" = length(file_path) == 1,
            "file_path must not be empty" = nzchar(file_path)
        )
        
        if (!file.exists(file_path)) {
            stop("File does not exist")
        }
        
        list(
            success = TRUE,
            data = file_path,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}
#' @title Validate Bigwig File Format
#' @description Validates bigwig file format integrity
#' @param file_path character Path to validated bigwig file
#' @param experiment_identifier character Experiment identifier
#' @return list Validation result {success, data, error}
#' @throws FormatValidationError if file format invalid
#' @importFrom rtracklayer import
#' @examples
#' bigwig_file_validate("path/to/file.bw", "exp001")
#' @seealso bigwig_file_exists
bigwig_file_validate <- function(file_path, experiment_identifier) {
    result <- tryCatch({
        stopifnot(
            "file_path must be character" = is.character(file_path),
            "experiment_identifier must be character" = is.character(experiment_identifier)
        )
        
        # Validate file existence first
        exists_result <- bigwig_file_exists(file_path)
        if (!exists_result$success) return(exists_result)
        
        # Try importing bigwig file
        track_data <- rtracklayer::import(file_path)
        
        list(
            success = TRUE,
            data = list(
                path = file_path,
                track = track_data,
                experiment = experiment_identifier
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = sprintf("Invalid bigwig format for %s: %s", 
                          experiment_identifier, e$message)
        )
    })
    
    return(result)
}

#' @title Filter Input Control Samples
#' @description Filters sample table for input controls
#' @param sample_table data.frame Sample metadata table
#' @param control_criteria list Control sample criteria
#' @return list Filtering result {success, data, error}
#' @throws FilterError if no controls found
#' @examples
#' input_controls_filter(samples, list(antibody = "Input"))
#' @seealso bigwig_controls_find
input_controls_filter <- function(sample_table, control_criteria) {
    result <- tryCatch({
        stopifnot(
            "sample_table must be data.frame" = is.data.frame(sample_table),
            "control_criteria must be list" = is.list(control_criteria)
        )
        
        control_samples <- sample_table
        for (field in names(control_criteria)) {
            control_samples <- control_samples[
                control_samples[[field]] == control_criteria[[field]], 
            ]
        }
        
        if (nrow(control_samples) == 0) {
            stop("No control samples found matching criteria")
        }
        
        list(
            success = TRUE,
            data = control_samples,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Find Control Bigwig Files
#' @description Locates bigwig files for control samples
#' @param control_samples data.frame Filtered control samples
#' @param directory_path character Bigwig directory path
#' @param file_pattern character Pattern for matching files
#' @return list File search result {success, data, error}
#' @throws FileSearchError if files not found
#' @examples
#' bigwig_controls_find(controls, "data/bigwig", ".*normalized.*\\.bw$")
#' @seealso input_controls_filter
bigwig_controls_find <- function(control_samples, directory_path, file_pattern) {
    result <- tryCatch({
        stopifnot(
            "control_samples must be data.frame" = is.data.frame(control_samples),
            "directory_path must be character" = is.character(directory_path),
            "file_pattern must be character" = is.character(file_pattern)
        )
        
        control_files <- lapply(control_samples$experiment_number, function(exp_id) {
            pattern <- sprintf(file_pattern, exp_id)
            files <- list.files(directory_path, pattern = pattern, 
                              full.names = TRUE)
            if (length(files) > 0) files[1] else NULL
        })
        
        valid_controls <- !sapply(control_files, is.null)
        if (!any(valid_controls)) {
            stop("No control files found")
        }
        
        list(
            success = TRUE,
            data = list(
                files = unlist(control_files[valid_controls]),
                samples = control_samples[valid_controls, ]
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Extract Track Values
#' @description Extracts numeric values from data track
#' @param track_object GenomicRanges Track object
#' @return list Extraction result {success, data, error}
#' @throws ExtractionError if values cannot be extracted
#' @importFrom GenomicRanges values
#' @examples
#' track_values_extract(data_track)
#' @seealso track_limits_calculate
track_values_extract <- function(track_object) {
    result <- tryCatch({
        stopifnot(
            "track_object must be GenomicRanges" = inherits(track_object, "GenomicRanges")
        )
        
        track_values <- GenomicRanges::values(track_object)
        if (length(track_values) == 0) {
            stop("No values found in track")
        }
        
        list(
            success = TRUE,
            data = track_values,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Calculate Track Value Range
#' @description Calculates range of track values with padding
#' @param track_values numeric Vector of track values
#' @param padding_percentage numeric Padding percentage
#' @return list Range calculation result {success, data, error}
#' @throws RangeError if calculation fails
#' @examples
#' track_range_calculate(values, 0.1)
#' @seealso track_values_extract
track_range_calculate <- function(track_values, padding_percentage = 0.1) {
    result <- tryCatch({
        stopifnot(
            "track_values must be numeric" = is.numeric(track_values),
            "padding_percentage must be numeric" = is.numeric(padding_percentage),
            "padding_percentage must be between 0 and 1" = 
                padding_percentage >= 0 && padding_percentage <= 1
        )
        
        value_range <- range(track_values, na.rm = TRUE)
        range_size <- diff(value_range)
        padding <- range_size * padding_percentage
        
        list(
            success = TRUE,
            data = list(
                min = value_range[1] - padding,
                max = value_range[2] + padding
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Calculate Track Visualization Limits
#' @description Main function calculating track visualization limits
#' @param track_list list List of track objects
#' @param options list Calculation options
#' @return list Limit calculation result {success, data, error}
#' @throws LimitCalculationError if process fails
#' @examples
#' track_limits_calculate(
#'   tracks,
#'   list(padding = 0.1, min_range = TRUE)
#' )
#' @seealso
#'   track_values_extract,
#'   track_range_calculate
track_limits_calculate <- function(track_list, options = list(padding = 0.1)) {
    result <- tryCatch({
        stopifnot(
            "track_list must be list" = is.list(track_list),
            "options must be list" = is.list(options)
        )
        
        # Extract values from all tracks
        track_values <- lapply(track_list, track_values_extract)
        successful_extracts <- vapply(track_values, `[[`, logical(1), "success")
        
        if (!any(successful_extracts)) {
            stop("No valid track values found")
        }
        
        # Combine all values
        all_values <- unlist(lapply(track_values[successful_extracts], 
                                  function(x) x$data))
        
        # Calculate range
        range_result <- track_range_calculate(all_values, options$padding)
        
        list(
            success = TRUE,
            data = range_result$data,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Find Fallback Control Track
#' @description Main function for finding fallback control tracks
#' @param sample_table data.frame Sample metadata
#' @param bigwig_directory character Directory path
#' @param search_pattern character File search pattern
#' @return list Control search result {success, data, error}
#' @throws ControlSearchError if no valid control found
#' @examples
#' control_track_find(
#'   samples,
#'   "data/bigwig",
#'   ".*normalized.*\\.bw$"
#' )
#' @seealso
#'   input_controls_filter,
#'   bigwig_controls_find,
#'   bigwig_file_validate
control_track_find <- function(sample_table, bigwig_directory, search_pattern) {
    result <- tryCatch({
        # Filter control samples
        controls <- input_controls_filter(sample_table, 
                                       list(antibody = "Input"))
        if (!controls$success) return(controls)
        
        # Find control files
        control_files <- bigwig_controls_find(controls$data, 
                                            bigwig_directory, 
                                            search_pattern)
        if (!control_files$success) return(control_files)
        
        # Validate each control file
        valid_controls <- lapply(seq_len(nrow(control_files$data$samples)), 
                               function(i) {
            bigwig_file_validate(
                control_files$data$files[i],
                control_files$data$samples$experiment_number[i]
            )
        })
        
        successful_controls <- vapply(valid_controls, `[[`, logical(1), "success")
        if (!any(successful_controls)) {
            stop("No valid control files found")
        }
        
        first_valid <- valid_controls[[which(successful_controls)[1]]]
        
        list(
            success = TRUE,
            data = first_valid$data,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Filter Comparison Samples
#' @description Filters samples based on comparison expression
#' @param sample_table data.frame Sample metadata table
#' @param comparison_expression language Quoted comparison expression
#' @return list Filtering result {success, data, error}
#' @throws FilterError if evaluation fails
#' @examples
#' comparison_samples_filter(samples, quote(treatment == "control"))
#' @seealso comparison_tracks_create
comparison_samples_filter <- function(sample_table, comparison_expression) {
    result <- tryCatch({
        stopifnot(
            "sample_table must be data.frame" = is.data.frame(sample_table),
            "comparison_expression must be language" = is.language(comparison_expression)
        )
        
        matched_samples <- eval(comparison_expression, envir = sample_table)
        filtered_samples <- sample_table[matched_samples, , drop = FALSE]
        
        if (nrow(filtered_samples) == 0) {
            stop("No samples match comparison criteria")
        }
        
        list(
            success = TRUE,
            data = filtered_samples,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Locate Control Sample
#' @description Finds matching control sample for comparison
#' @param sample_table data.frame Sample metadata
#' @param comparison_sample data.frame Single comparison sample
#' @param control_criteria list Control matching criteria
#' @return list Control sample result {success, data, error}
#' @throws ControlError if no matching control found
#' @examples
#' control_sample_locate(samples, comparison, list(antibody = "Input"))
#' @seealso control_track_create
control_sample_locate <- function(sample_table, comparison_sample, control_criteria) {
    result <- tryCatch({
        stopifnot(
            "sample_table must be data.frame" = is.data.frame(sample_table),
            "comparison_sample must be data.frame" = is.data.frame(comparison_sample),
            "control_criteria must be list" = is.list(control_criteria),
            "comparison_sample must have one row" = nrow(comparison_sample) == 1
        )
        
        # Combine base criteria with comparison-specific matching
        full_criteria <- c(
            control_criteria,
            list(rescue_allele = comparison_sample$rescue_allele[1])
        )
        
        # Apply all criteria
        control_matches <- Reduce(function(acc, criterion_name) {
            acc & sample_table[[criterion_name]] == full_criteria[[criterion_name]]
        }, names(full_criteria), init = TRUE)
        
        control_sample <- sample_table[control_matches, , drop = FALSE][1, ]
        
        if (nrow(control_sample) == 0) {
            stop("No matching control sample found")
        }
        
        list(
            success = TRUE,
            data = control_sample,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Create Control Track
#' @description Creates visualization track for control sample
#' @param control_sample data.frame Control sample data
#' @param bigwig_directory character Directory containing bigwig files
#' @param track_options list Track visualization options
#' @return list Track creation result {success, data, error}
#' @throws TrackError if track creation fails
#' @importFrom rtracklayer import
#' @examples
#' control_track_create(control, "data/bigwig", list(color = "#808080"))
#' @seealso comparison_tracks_create
control_track_create <- function(control_sample, bigwig_directory, track_options) {
    result <- tryCatch({
        stopifnot(
            "control_sample must be data.frame" = is.data.frame(control_sample),
            "bigwig_directory must be character" = is.character(bigwig_directory),
            "track_options must be list" = is.list(track_options),
            "control_sample must have one row" = nrow(control_sample) == 1
        )
        
        # Locate bigwig file
        bigwig_pattern <- sprintf(
            "%s.*normalized.*\\.bw$",
            control_sample$experiment_number
        )
        bigwig_file <- list.files(
            bigwig_directory,
            pattern = bigwig_pattern,
            full.names = TRUE
        )[1]
        
        if (is.na(bigwig_file)) {
            stop("Control bigwig file not found")
        }
        
        # Create track
        track_data <- rtracklayer::import(bigwig_file)
        track <- DataTrack(
            track_data,
            type = "l",
            name = paste("Input", control_sample$rescue_allele, sep = "_"),
            col = track_options$color %||% "#808080"
        )
        
        list(
            success = TRUE,
            data = list(
                track = track,
                source = bigwig_file
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Create Sample Track
#' @description Creates visualization track for comparison sample
#' @param sample_data data.frame Single sample data
#' @param bigwig_directory character Directory containing bigwig files
#' @param track_options list Track visualization options
#' @return list Track creation result {success, data, error}
#' @throws TrackError if track creation fails
#' @examples
#' sample_track_create(sample, "data/bigwig", list(color = "#fd0036"))
#' @seealso comparison_tracks_create
sample_track_create <- function(sample_data, bigwig_directory, track_options) {
    result <- tryCatch({
        stopifnot(
            "sample_data must be data.frame" = is.data.frame(sample_data),
            "bigwig_directory must be character" = is.character(bigwig_directory),
            "track_options must be list" = is.list(track_options),
            "sample_data must have one row" = nrow(sample_data) == 1
        )
        
        # Locate bigwig file
        bigwig_pattern <- sprintf(
            "%s.*normalized.*\\.bw$",
            sample_data$experiment_number
        )
        bigwig_file <- list.files(
            bigwig_directory,
            pattern = bigwig_pattern,
            full.names = TRUE
        )[1]
        
        if (is.na(bigwig_file)) {
            stop("Sample bigwig file not found")
        }
        
        # Create track
        track_data <- rtracklayer::import(bigwig_file)
        track <- DataTrack(
            track_data,
            type = "l",
            name = track_options$label,
            col = track_options$color %||% "#fd0036",
            chromosome = track_options$chromosome
        )
        
        list(
            success = TRUE,
            data = list(
                track = track,
                source = bigwig_file
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Create Feature Track
#' @description Creates genome feature visualization track
#' @param feature_data GRanges Feature genomic ranges
#' @param track_options list Track visualization options
#' @return list Track creation result {success, data, error}
#' @throws TrackError if track creation fails
#' @examples
#' feature_track_create(features, list(type = "annotation"))
#' @seealso comparison_tracks_create
feature_track_create <- function(feature_data, track_options) {
    result <- tryCatch({
        stopifnot(
            "feature_data must be GRanges" = inherits(feature_data, "GRanges"),
            "track_options must be list" = is.list(track_options)
        )
        
        track <- Gviz::AnnotationTrack(
            feature_data,
            name = if (is.null(track_options$name)) "Features" else track_options$name, # Use to replace %||%
            chromosome = track_options$chromosome,
            genome = track_options$genome
        )
        
        list(
            success = TRUE,
            data = list(
                track = track,
                features = feature_data
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Create Plot Output Path
#' @description Generates standardized plot output file path
#' @param base_directory character Base output directory
#' @param plot_parameters list Plot identification parameters
#' @return list Path generation result {success, data, error}
#' @throws PathError if path creation fails
#' @examples
#' plot_path_create("output", list(timestamp = "20240108", id = "exp1"))
#' @seealso comparison_plot_create
plot_path_create <- function(base_directory, plot_parameters) {
    result <- tryCatch({
        stopifnot(
            "base_directory must be character" = is.character(base_directory),
            "plot_parameters must be list" = is.list(plot_parameters),
            "required parameters missing" = all(c("timestamp", "experiment_id", 
                                               "comparison_name", "chromosome") %in% 
                                             names(plot_parameters))
        )
        
        file_name <- sprintf(
            "%s_%s_%s_chr%s_all_comparisons.svg",
            plot_parameters$timestamp,
            plot_parameters$experiment_id,
            plot_parameters$comparison_name,
            plot_parameters$chromosome
        )
        
        output_path <- file.path(base_directory, "plots", file_name)
        
        # Ensure directory exists
        dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
        
        list(
            success = TRUE,
            data = output_path,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Create Comparison Plot
#' @description Creates visualization plot for comparison
#' @param track_list list List of visualization tracks
#' @param plot_options list Plot configuration options
#' @param output_path character Plot output file path
#' @return list Plot creation result {success, data, error}
#' @throws PlotError if plot creation fails
#' @examples
#' comparison_plot_create(tracks, list(title = "Comparison"), "plot.svg")
#' @seealso
#'   comparison_tracks_create,
#'   plot_path_create
comparison_plot_create <- function(track_list, plot_options, output_path) {
    result <- tryCatch({
        stopifnot(
            "track_list must be list" = is.list(track_list),
            "plot_options must be list" = is.list(plot_options),
            "output_path must be character" = is.character(output_path)
        )
        
        # Calculate track limits if requested
        y_limits <- if (plot_options$calculate_limits) {
            track_values <- lapply(track_list, function(track) {
                if (inherits(track, "DataTrack")) {
                    values(track)
                } else {
                    NULL
                }
            })
            track_values <- unlist(track_values[!sapply(track_values, is.null)])
            if (length(track_values) > 0) {
                range_raw <- range(track_values, na.rm = TRUE)
                range_size <- diff(range_raw)
                c(range_raw[1] - range_size * 0.1,
                  range_raw[2] + range_size * 0.1)
            } else {
                NULL
            }
        } else {
            NULL
        }
        
        # Create plot
        svg(output_path)
        on.exit(dev.off(), add = TRUE)
        
        plotTracks(
            track_list,
            main = plot_options$title,
            chromosome = plot_options$chromosome,
            ylim = y_limits
        )
        
        list(
            success = TRUE,
            data = list(
                path = output_path,
                limits = y_limits
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Process Comparison Analysis
#' @description Main function orchestrating comparison analysis
#' @param sample_table data.frame Sample metadata
#' @param comparison_config list Comparison configuration
#' @param output_options list Output configuration
#' @return list Analysis result {success, data, error}
#' @throws ComparisonError if analysis fails
#' @examples
#' comparison_analysis_process(
#'   samples,
#'   list(name = "test", expression = quote(treatment == "control")),
#'   list(directory = "output", format = "svg")
#' )
#' @seealso
#'   comparison_samples_filter,
#'   control_track_create,
#'   sample_track_create,
#'   comparison_plot_create
comparison_analysis_process <- function(sample_table, comparison_config, 
                                      output_options) {
    result <- tryCatch({
        stopifnot(
            "sample_table must be data.frame" = is.data.frame(sample_table),
            "comparison_config must be list" = is.list(comparison_config),
            "output_options must be list" = is.list(output_options)
        )
        
        # Filter comparison samples
        samples_result <- comparison_samples_filter(
            sample_table, 
            comparison_config$expression
        )
        if (!samples_result$success) return(samples_result)
        
        # Initialize tracks list
        tracks <- list(GenomeAxisTrack(
            name = sprintf("Chr %s Axis", output_options$chromosome)
        ))
        
        # Create control track
        control_result <- control_sample_locate(
            sample_table,
            samples_result$data[1, ],
            list(antibody = "Input")
        )
        if (control_result$success) {
            control_track <- control_track_create(
                control_result$data,
                output_options$bigwig_directory,
                list(color = "#808080")
            )
            if (control_track$success) {
                tracks[[length(tracks) + 1]] <- control_track$data$track
            }
        }
        
        # Create sample tracks
        sample_tracks <- lapply(seq_len(nrow(samples_result$data)), function(i) {
            sample_track_create(
                samples_result$data[i, ],
                output_options$bigwig_directory,
                list(
                    label = comparison_config$labels[i],
                    color = "#fd0036",
                    chromosome = output_options$chromosome
                )
            )
        })
        
        successful_tracks <- vapply(sample_tracks, `[[`, logical(1), "success")
        tracks <- c(
            tracks,
            lapply(sample_tracks[successful_tracks], function(x) x$data$track)
        )
        
        # Create plot
        plot_result <- comparison_plot_create(
            tracks,
            list(
                title = comparison_config$name,
                chromosome = output_options$chromosome,
                calculate_limits = TRUE
            ),
            output_options$output_path
        )
        
        list(
            success = TRUE,
            data = list(
                tracks = tracks,
                plot = plot_result$data
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Validate Required Package Dependencies
#' @description Validates presence of required R packages
#' @param package_list character Vector of package names
#' @return list Validation result {success, data, error}
#' @throws DependencyError if packages missing
#' @examples
#' packages_required_validate(c("rtracklayer", "GenomicRanges"))
#' @seealso experiment_environment_validate
packages_required_validate <- function(package_list) {
    result <- tryCatch({
        stopifnot(
            "package_list must be character vector" = is.character(package_list),
            "package_list cannot be empty" = length(package_list) > 0
        )
        
        missing_packages <- vapply(package_list, function(pkg) {
            !requireNamespace(pkg, quietly = TRUE)
        }, logical(1))
        
        if (any(missing_packages)) {
            stop(sprintf(
                "Missing required packages: %s",
                paste(names(missing_packages)[missing_packages], collapse = ", ")
            ))
        }
        
        list(
            success = TRUE,
            data = package_list,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Validate Experiment Environment
#' @description Validates experiment directory structure and requirements
#' @param experiment_identifier character Experiment ID
#' @param environment_requirements list Required directory structure
#' @return list Validation result {success, data, error}
#' @throws EnvironmentError if validation fails
#' @examples
#' experiment_environment_validate("exp001", list(directories = c("coverage")))
#' @seealso experiment_files_validate
experiment_environment_validate <- function(experiment_identifier, 
                                         environment_requirements) {
    result <- tryCatch({
        stopifnot(
            "experiment_identifier must be character" = is.character(experiment_identifier),
            "environment_requirements must be list" = is.list(environment_requirements),
            "directories must be specified" = "directories" %in% 
                names(environment_requirements)
        )
        
        base_path <- file.path(Sys.getenv("HOME"), "data", experiment_identifier)
        
        if (!dir.exists(base_path)) {
            stop("Experiment directory not found")
        }
        
        missing_dirs <- vapply(environment_requirements$directories, function(dir) {
            !dir.exists(file.path(base_path, dir))
        }, logical(1))
        
        if (any(missing_dirs)) {
            stop(sprintf(
                "Missing required directories: %s",
                paste(names(missing_dirs)[missing_dirs], collapse = ", ")
            ))
        }
        
        list(
            success = TRUE,
            data = list(
                base_path = base_path,
                validated_directories = environment_requirements$directories
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}
#' @title Validate Experiment Files
#' @description Validates presence and format of experiment files
#' @param experiment_directory character Base directory path
#' @param file_requirements list Required file patterns
#' @return list Validation result {success, data, error}
#' @throws FileValidationError if files invalid
#' @examples
#' experiment_files_validate("path/to/exp", list(bigwig = "\\.bw$"))
#' @seealso experiment_environment_validate
experiment_files_validate <- function(experiment_directory, file_requirements) {
    result <- tryCatch({
        stopifnot(
            "experiment_directory must be character" = is.character(experiment_directory),
            "file_requirements must be list" = is.list(file_requirements)
        )
        
        validation_results <- lapply(names(file_requirements), function(file_type) {
            pattern <- file_requirements[[file_type]]
            files <- list.files(
                experiment_directory,
                pattern = pattern,
                recursive = TRUE
            )
            
            if (length(files) == 0) {
                return(sprintf("No %s files found", file_type))
            }
            return(NULL)
        })
        
        errors <- unlist(validation_results[!sapply(validation_results, is.null)])
        
        if (length(errors) > 0) {
            stop(paste(errors, collapse = "; "))
        }
        
        list(
            success = TRUE,
            data = file_requirements,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Create Genomic Range Object
#' @description Creates standardized genomic range for visualization
#' @param chromosome_number numeric Chromosome number
#' @param range_parameters list Range configuration parameters
#' @return list Range creation result {success, data, error}
#' @throws RangeCreationError if range invalid
#' @importFrom GenomicRanges GRanges
#' @examples
#' genomic_range_create(10, list(start = 1, end = 1e6))
#' @seealso track_group_create
genomic_range_create <- function(chromosome_number, range_parameters) {
    result <- tryCatch({
        stopifnot(
            "chromosome_number must be numeric" = is.numeric(chromosome_number),
            "range_parameters must be list" = is.list(range_parameters),
            "chromosome_number must be between 1 and 16" = 
                chromosome_number >= 1 && chromosome_number <= 16,
            "start and end must be specified" = all(c("start", "end") %in% 
                                                     names(range_parameters))
        )
        
        chromosome_roman <- paste0("chr", utils::as.roman(chromosome_number))
        
        range <- GenomicRanges::GRanges(
            seqnames = chromosome_roman,
            ranges = IRanges::IRanges(
                start = range_parameters$start,
                end = range_parameters$end
            ),
            strand = "*"
        )
        
        list(
            success = TRUE,
            data = range,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Calculate Track Data Range
#' @description Calculates value range from track data
#' @param track_data GenomicRanges Track data object
#' @param range_options list Range calculation options
#' @return list Range calculation result {success, data, error}
#' @throws RangeCalculationError if calculation fails
#' @examples
#' track_range_calculate(track_data, list(padding = 0.1))
#' @seealso track_group_range_calculate
track_range_calculate <- function(track_data, range_options) {
    result <- tryCatch({
        stopifnot(
            "track_data must be GenomicRanges" = inherits(track_data, "GenomicRanges"),
            "range_options must be list" = is.list(range_options),
            "padding percentage must be specified" = "padding" %in% 
                names(range_options)
        )
        
        track_values <- GenomicRanges::values(track_data)
        if (length(track_values) == 0) {
            stop("No values found in track data")
        }
        
        value_range <- range(track_values, na.rm = TRUE)
        range_size <- diff(value_range)
        padding <- range_size * range_options$padding
        
        list(
            success = TRUE,
            data = list(
                min = value_range[1] - padding,
                max = value_range[2] + padding
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}
#' @title Calculate Track Group Range
#' @description Calculates combined range for track group
#' @param track_list list List of track objects
#' @param range_options list Range calculation options
#' @return list Range calculation result {success, data, error}
#' @throws GroupRangeError if calculation fails
#' @examples
#' track_group_range_calculate(tracks, list(padding = 0.1))
#' @seealso track_range_calculate
track_group_range_calculate <- function(track_list, range_options) {
    result <- tryCatch({
        stopifnot(
            "track_list must be list" = is.list(track_list),
            "range_options must be list" = is.list(range_options)
        )
        
        # Calculate ranges for each track
        track_ranges <- lapply(track_list, function(track) {
            if (inherits(track, "DataTrack")) {
                track_range_calculate(track, range_options)
            } else {
                NULL
            }
        })
        
        valid_ranges <- track_ranges[vapply(track_ranges, 
                                          function(x) !is.null(x) && x$success, 
                                          logical(1))]
        
        if (length(valid_ranges) == 0) {
            stop("No valid ranges found in track group")
        }
        
        # Combine ranges
        global_min <- min(vapply(valid_ranges, function(x) x$data$min, numeric(1)))
        global_max <- max(vapply(valid_ranges, function(x) x$data$max, numeric(1)))
        
        list(
            success = TRUE,
            data = list(
                min = global_min,
                max = global_max
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Create Placeholder Track
#' @description Creates an empty visualization track for missing bigwig data
#' @param sample_name character Name of the sample
#' @param track_options list Track visualization options
#' @return list Track creation result {success, data, error}
track_placeholder_create <- function(sample_name, track_options) {
    result <- tryCatch({
        stopifnot(
            "sample_name must be character" = is.character(sample_name),
            "track_options must be list" = is.list(track_options)
        )
        
        # Create empty GRanges object
        empty_ranges <- GenomicRanges::GRanges(
            seqnames = track_options$chromosome,
            ranges = IRanges::IRanges(start = 1, end = 1),
            score = 0
        )
        
        # Create empty track with clear visual indication
        empty_track <- Gviz::DataTrack(
            empty_ranges,
            type = "l",
            name = paste(sample_name, "(No Data)"),
            col = track_options$placeholder_color %||% "#cccccc",
            chromosome = track_options$chromosome,
            background.title = "#f0f0f0"  # Light gray background for placeholder
        )
        
        list(
            success = TRUE,
            data = list(
                track = empty_track,
                source = "placeholder",
                metadata = list(
                    is_placeholder = TRUE,
                    original_name = sample_name
                )
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Create Single Track
#' @description Creates visualization track for single sample
#' @param sample_data list Sample metadata and configuration
#' @param track_options list Track visualization options
#' @return list Track creation result {success, data, error}
#' @throws TrackCreationError if creation fails
#' @importFrom Gviz DataTrack
#' @examples
#' track_single_create(sample_data, list(color = "#fd0036"))
#' @seealso track_group_create

track_single_create <- function(sample_data, track_options) {
    result <- tryCatch({
        stopifnot(
            "sample_data must be track_configuration" = inherits(sample_data, "track_configuration"),
            "track_options must be list" = is.list(track_options),
            "required track configuration fields missing" = 
                all(c("bigwig_file", "name") %in% names(sample_data))
        )
        
        # Extract validated paths and names
        bigwig_path <- sample_data$bigwig_file
        sample_name <- sample_data$name
        
        # Validate file existence
        if (is.na(bigwig_path) || !file.exists(bigwig_path)) {
            if (track_options$verbose) {
                message(sprintf("Creating placeholder track for sample: %s", sample_name))
            }
            return(track_placeholder_create(
                sample_name = sample_name,
                track_options = track_options
            ))
        }
        
        # Validate file size
        if (file.size(bigwig_path) == 0) {
            if (track_options$verbose) {
                message(sprintf("Empty bigwig file for sample: %s", sample_name))
            }
            return(track_placeholder_create(
                sample_name = sample_name,
                track_options = track_options
            ))
        }
        
        if (track_options$verbose) {
            message(sprintf("Creating track from: %s", basename(bigwig_path)))
        }
        
        # Import track data
        track_data <- rtracklayer::import(bigwig_path)
        
        # Create visualization track
        track <- Gviz::DataTrack(
            track_data,
            type = track_options$type %||% "l",
            name = sample_name,
            col = track_options$color %||% "#fd0036",
            chromosome = track_options$chromosome
        )
        
        list(
            success = TRUE,
            data = list(
                track = track,
                source = bigwig_path,
                metadata = sample_data$metadata  # Preserve any additional metadata
            ),
            error = NULL
        )
    }, error = function(e) {
        if (track_options$verbose) {
            message(sprintf("Error creating track: %s", e$message))
        }
        # On any error, create placeholder track
        return(track_placeholder_create(
            sample_name = sample_data$name,
            track_options = track_options
        ))
    })
    
    return(result)
}

#' @title Create Track Group
#' @description Creates group of visualization tracks
#' @param sample_list list List of sample data
#' @param group_options list Group visualization options
#' @return list Track group creation result {success, data, error}
#' @throws GroupCreationError if creation fails
#' @examples
#' track_group_create(samples, list(samples_per_page = 4))
#' @seealso track_single_create
track_group_create <- function(sample_list, group_options) {
    result <- tryCatch({
        stopifnot(
            "sample_list must be list" = is.list(sample_list),
            "group_options must be list" = is.list(group_options)
        )
        
        # Create axis track
        tracks <- list(
            Gviz::GenomeAxisTrack(
                name = sprintf("Chr %s", group_options$chromosome)
            )
        )
        
        # Create individual tracks (now always succeeding with placeholders)
        sample_tracks <- lapply(sample_list, function(sample) {
            track_single_create(sample, group_options)
        })
        
        # All tracks should succeed now (either real or placeholder)
        tracks <- c(
            tracks,
            lapply(sample_tracks, function(x) x$data$track)
        )
        
        list(
            success = TRUE,
            data = list(
                tracks = tracks,
                count = length(tracks)
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Generate Plot Output Path
#' @description Generates standardized plot file path
#' @param base_directory character Base output directory
#' @param plot_parameters list Plot identification parameters
#' @return list Path generation result {success, data, error}
#' @throws PathGenerationError if generation fails
#' @examples
#' plot_path_generate("output", list(chromosome = 10, group = 1))
#' @seealso plot_tracks_create
plot_path_generate <- function(base_directory, plot_parameters) {
    result <- tryCatch({
        stopifnot(
            "base_directory must be character" = is.character(base_directory),
            "plot_parameters must be list" = is.list(plot_parameters),
            "required parameters missing" = all(c("chromosome", "group") %in% 
                                               names(plot_parameters))
        )
        
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        filename <- sprintf(
            "%s_overview_chr%s_group%d.svg",
            timestamp,
            plot_parameters$chromosome,
            plot_parameters$group
        )
        
        output_path <- file.path(base_directory, "plots", filename)
        dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
        
        list(
            success = TRUE,
            data = output_path,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}
#' @title Create Track Visualization
#' @description Creates track visualization plot
#' @param track_group list Group of visualization tracks
#' @param plot_options list Plot configuration options
#' @return list Plot creation result {success, data, error}
#' @throws PlotCreationError if creation fails
#' @examples
#' plot_tracks_create(tracks, list(width = 10, height = 8))
#' @seealso track_group_create
plot_tracks_create <- function(track_group, plot_options) {
    result <- tryCatch({
        stopifnot(
            "track_group must be list" = is.list(track_group),
            "plot_options must be list" = is.list(plot_options),
            "required options missing" = all(c("width", "height", "output_path") %in% 
                                            names(plot_options))
        )
        
        # Calculate track limits if requested
        y_limits <- if (plot_options$calculate_limits) {
            range_result <- track_group_range_calculate(
                track_group$tracks,
                list(padding = 0.1)
            )
            if (range_result$success) {
                c(range_result$data$min, range_result$data$max)
            } else {
                NULL
            }
        } else {
            NULL
        }
        
        # Create plot
        svg(plot_options$output_path, 
            width = plot_options$width, 
            height = plot_options$height)
        on.exit(dev.off(), add = TRUE)
        
        Gviz::plotTracks(
            track_group$tracks,
            chromosome = plot_options$chromosome,
            ylim = y_limits
        )
        
        list(
            success = TRUE,
            data = list(
                path = plot_options$output_path,
                limits = y_limits
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Calculate Global Range Across Bigwig Files
#' @description Calculates combined value range across multiple bigwig files
#' @param bigwig_files character Vector of bigwig file paths
#' @param genome_range GRanges Genomic range to consider
#' @return list Range calculation result {success, data, error}
#' @throws RangeCalculationError if calculation fails
#' @examples
#' get_global_range(c("sample1.bw", "sample2.bw"), genome_range)
#' @seealso track_range_calculate
get_global_range <- function(bigwig_files, genome_range) {
    result <- tryCatch({
        stopifnot(
            "bigwig_files must be character vector" = is.character(bigwig_files),
            "bigwig_files cannot be empty" = length(bigwig_files) > 0,
            "genome_range must be GRanges" = inherits(genome_range, "GRanges")
        )
        
        # Calculate range for each file
        ranges <- lapply(bigwig_files, function(bw_file) {
            tryCatch({
                # Import bigwig data for specified range
                track_data <- rtracklayer::import(bw_file, which = genome_range)
                if (length(track_data) == 0) {
                    return(NULL)
                }
                
                # Extract values and calculate range
                values <- GenomicRanges::values(track_data)$score
                if (length(values) == 0) {
                    return(NULL)
                }
                
                range(values, na.rm = TRUE)
            }, error = function(e) {
                warning(sprintf("Failed to process %s: %s", 
                              basename(bw_file), e$message))
                return(NULL)
            })
        })
        
        # Remove NULL results and combine ranges
        valid_ranges <- do.call(rbind, ranges[!sapply(ranges, is.null)])
        
        if (nrow(valid_ranges) == 0) {
            stop("No valid ranges found in any bigwig file")
        }
        
        # Calculate global min and max
        global_range <- c(
            min(valid_ranges[, 1], na.rm = TRUE),
            max(valid_ranges[, 2], na.rm = TRUE)
        )
        
        list(
            success = TRUE,
            data = global_range,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}


create_bigwig_sample_mapping <- function(sample_table, bigwig_files) {
    result <- tryCatch({
        stopifnot(
            "sample_table must be data.frame" = is.data.frame(sample_table),
            "bigwig_files must be character" = is.character(bigwig_files),
            "sample_id column must exist" = "sample_id" %in% colnames(sample_table)
        )
        
        # Create mapping between sample IDs and bigwig files
        bigwig_mapping <- sapply(
            X = sample_table$sample_id,
            FUN = function(sample_id) {
                # Remove file extension for matching
                clean_sample_id <- sub("\\.bw$", "", sample_id)
                matching_file <- bigwig_files[grepl(clean_sample_id, bigwig_files)]
                if (length(matching_file) > 0) matching_file[1] else NA_character_
            },
            USE.NAMES = TRUE
        )
        
        list(
            success = TRUE,
            data = bigwig_mapping,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}


# Define a formal data structure contract
#' @title Create Standard Track Configuration
#' @description Creates standardized track configuration structure
track_configuration_create <- function(sample_data) {
    result <- tryCatch({
        stopifnot(
            "sample_data must be list" = is.list(sample_data),
            "required fields missing" = all(c("bigwig_file", "name") %in% names(sample_data))
        )
        
        # Create standardized structure
        track_config <- list(
            bigwig_file = sample_data$bigwig_file,
            name = sample_data$name,
            metadata = list()  # Optional additional metadata
        )
        
        class(track_config) <- "track_configuration"
        
        list(
            success = TRUE,
            data = track_config,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}


create_sample_track_configs <- function(group_samples, bigwig_mapping) {
    result <- tryCatch({
        configs <- lapply(seq_len(nrow(group_samples)), function(i) {
            current_sample_id <- group_samples$sample_id[i]
            config_result <- track_configuration_create(list(
                bigwig_file = bigwig_mapping[current_sample_id],
                name = current_sample_id
            ))
            if (!config_result$success) stop(config_result$error)
            config_result$data
        })
        
        list(
            success = TRUE,
            data = configs,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}
