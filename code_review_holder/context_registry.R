# context_registry.R

FUNCTION_REGISTRY <- list(
    data_processing = c(
        "metadata_path_validate",
        "metadata_file_read",
        "metadata_schema_validate",
        "metadata_comparisons_process",
        "comparison_samples_filter",
        "comparison_expression_validate",
        "comparison_expression_evaluate",
        "comparison_results_format"
    ),
    genomic_data = c(
        "chromosome_mappings_initialize",
        "chromosome_style_validate",
        "chromosome_names_clean",
        "chromosome_names_to_ucsc",
        "chromosome_names_to_roman",
        "chromosome_names_to_numeric",
        "chromosome_names_convert"
    ),
    visualization = c(
        "bigwig_file_exists",
        "bigwig_file_validate",
        "track_values_extract",
        "track_range_calculate",
        "track_group_range_calculate",
        "track_single_create",
        "track_group_create",
        "plot_path_generate",
        "plot_tracks_create",
        "get_global_range"
    )
)

DATA_CONTRACTS <- list(
    result_structure = list(
        success = "logical(1)",
        data = "list()",
        error = "character(1)"
    ),
    track_data = list(
        track = "GenomicRanges",
        source = "character(1)",
        metadata = "list()"
    ),
    plot_options = list(
        width = "numeric(1)",
        height = "numeric(1)",
        color = "character(1)",
        chromosome = "character(1)"
    )
)

SHARED_CONSTANTS <- list(
    PLOT_CONFIG = list(
        SAMPLES_PER_PAGE = 4,
        DEFAULT_CHROMOSOME = 10,
        TRACK_COLOR = "#fd0036",
        WIDTH = 10,
        HEIGHT = 8
    ),
    REQUIRED_PACKAGES = c(
        "rtracklayer",
        "GenomicRanges",
        "Gviz",
        "tidyverse"
    )
)

FUNCTION_DEPENDENCIES <- list(
    metadata_comparisons_process = c(
        "comparison_samples_filter",
        "comparison_expression_validate",
        "comparison_expression_evaluate"
    ),
    track_group_create = c(
        "track_single_create",
        "bigwig_file_validate",
        "track_values_extract"
    ),
    plot_tracks_create = c(
        "track_group_range_calculate",
        "track_range_calculate"
    ),
    get_global_range = c(
        "track_range_calculate"
    )
)

# Function to validate function existence
validate_function_existence <- function() {
    all_functions <- unlist(FUNCTION_REGISTRY)
    missing_functions <- all_functions[!sapply(all_functions, exists)]
    if (length(missing_functions) > 0) {
        stop("Missing functions: ", paste(missing_functions, collapse = ", "))
    }
    return(TRUE)
}

# Function to check dependencies
#validate_dependencies <- function() {
#    missing_deps <- lapply(names(FUNCTION_DEPENDENCIES), function(func) {
#        deps <- FUNCTION_DEPENDENCIES[[func]]
#        missing <- deps[!sapply(deps, exists)]
#        if (length(missing) > 0) {
#            return(list(function = func, missing = missing))
#        }
#        return(NULL)
#    })
#    missing_deps <- missing_deps[!sapply(missing_deps, is.null)]
#    if (length(missing_deps) > 0) {
#        stop("Unmet dependencies found")
#    }
#    return(TRUE)
#}
