# Define analysis-focused override presets
OVERRIDE_PRESETS <- list(
    full_script_with_debug = list(  # For checking pipeline behavior
        debug_verbose = TRUE,
        debug_validate = TRUE,
        process_single_file = TRUE,
        process_single_comparison = TRUE,
        output_dry_run = FALSE
    ),
    full_analysis = list(  # For complete ChIP-seq analysis
        debug_verbose = FALSE,
        debug_validate = FALSE,
        process_single_file = FALSE,
        output_dry_run = FALSE
    ),
    smoke_test = list(  # For quick pipeline validation
        debug_verbose = TRUE,
        debug_validate = TRUE,
        process_single_file = TRUE,
        process_samples_per_batch = 2,
        output_display_time = 5
    ),
    dry_run_mode = list(  # For quick pipeline validation
        debug_verbose = TRUE,
        debug_validate = TRUE,
        process_single_file = FALSE,
        process_single_comparison = FALSE,
        output_dry_run = TRUE
    ),
    dry_run_single_sample = list(  # For investigating specific sample issues
        debug_verbose = TRUE,
        output_dry_run = TRUE,
        process_single_file = TRUE,
        process_comparison = "comp_1108forNoneAndWT",
        process_single_comparison = FALSE,
        process_file_index = 58
    ),
    output_single_sample = list(  # For investigating specific sample issues
        debug_verbose = TRUE,
        output_dry_run = FALSE,
        process_single_file = TRUE,
        process_comparison = "comp_1108forNoneAndWT",
        process_single_comparison = FALSE,
        process_file_index = 58
    )
)
