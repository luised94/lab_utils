# Define analysis-focused override presets
OVERRIDE_PRESETS <- list(
    inspect_pipeline = list(  # For checking pipeline behavior
        debug_interactive = TRUE,
        debug_verbose = TRUE,
        debug_validate = TRUE,
        process_single_file = TRUE,
        output_dry_run = FALSE
    ),
    full_analysis = list(  # For complete ChIP-seq analysis
        debug_interactive = FALSE,
        debug_verbose = FALSE,
        debug_validate = FALSE,
        process_single_file = FALSE,
        output_dry_run = FALSE
    ),
    test_pipeline = list(  # For quick pipeline validation
        debug_verbose = TRUE,
        debug_validate = TRUE,
        process_single_file = TRUE,
        process_samples_per_batch = 2,
        output_display_time = 5
    ),
    full_inspect_pipeline = list(  # For quick pipeline validation
        debug_verbose = TRUE,
        debug_validate = TRUE,
        process_single_file = FALSE,
        process_samples_per_batch = 2,
        output_dry_run = TRUE
    ),
    check_single_sample = list(  # For investigating specific sample issues
        debug_verbose = TRUE,
        output_dry_run = TRUE,
        process_single_file = TRUE,
        process_file_index = 58
    )
)
