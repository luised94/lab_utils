# Test validation
test_validation <- validate_comparison_inputs(
    mtcars,
    list(high_mpg = quote(mpg > 20))
)
print(test_validation$valid) # Should be TRUE

# Test comparison execution
test_exec <- execute_comparison(mtcars, quote(mpg > 20))
print(sum(test_exec$matches)) # Should show count of high MPG cars

# Test summary creation
test_summary <- create_sample_summary(
    mtcars[1, ],
    c("mpg", "cyl")
)
print(test_summary$summary)

# Test full analysis
test_analysis <- analyze_comparisons(
    mtcars,
    list(
        high_mpg = quote(mpg > 20),
        high_cyl = quote(cyl > 6)
    )
)
print(sapply(test_analysis$results, function(x) x$match_count))


# Test data
test_df <- data.frame(
    antibody = c("AB1", "AB2"),
    type = c("T1", "T2"),
    treatment = c("TR1", "TR1"),
    stringsAsFactors = FALSE
)

# Test validation
test_valid <- validate_category_names(test_df, c("type", "treatment"))
print(test_valid$success) # Should be TRUE

# Test category extraction
test_extract <- extract_category_values(test_df, c("antibody", "type"))
print(test_extract$data) # Should show unique values

# Test full label generation
test_labels <- generate_sample_labels(
    test_df,
    c("type", "treatment"),
    list(verbose = TRUE)
)
print(test_labels$data) # Should show formatted labels

# Test directory validation
test_dir_valid <- directory_structure_validate("test/path")

# Test file scanning
test_files <- experiment_files_scan("test/path", 
                                  "consolidated_.*_sequence\\.fastq$")

# Test identifier extraction
test_ids <- experiment_identifiers_extract(
    c("consolidated_12345_sequence.fastq"),
    "consolidated_([0-9]{5,6})_sequence\\.fastq"
)

# Test metadata processing
test_metadata <- experiment_metadata_process(
    "test/path",
    list(
        categories = list(treatment = c("control", "treated")),
        column_order = c("sample_id", "treatment")
    ),
    list(
        output_file = FALSE,
        output_path = NULL
    )
)

# Test mapping initialization
test_mappings <- chromosome_mappings_initialize()

# Test style detection
test_detect <- chromosome_style_detect(
    c("chr1", "chr2"),
    list(
        UCSC = "^chr[0-9]+$",
        Roman = "^chr[IVX]+$",
        Numeric = "^[0-9]+$"
    )
)

# Test conversion
test_convert <- chromosome_names_convert(
    c("chr1", "chr2"),
    "Roman",
    list(validate = TRUE)
)

# Test file validation
test_file <- bigwig_file_validate("test.bw", "exp001")

# Test control filtering
test_controls <- input_controls_filter(
    data.frame(antibody = c("Input", "H3K4me3")),
    list(antibody = "Input")
)

# Test track limits
test_limits <- track_limits_calculate(
    list(track1, track2),
    list(padding = 0.1)
)

# Test sample filtering
test_filter <- comparison_samples_filter(
    data.frame(treatment = c("control", "treated")),
    quote(treatment == "control")
)

# Test track creation
test_track <- sample_track_create(
    data.frame(experiment_number = "001", row.names = 1),
    "data/bigwig",
    list(label = "Sample 1", color = "#fd0036")
)

# Test full analysis
test_analysis <- comparison_analysis_process(
    sample_table,
    list(
        name = "Test Comparison",
        expression = quote(treatment == "control"),
        labels = c("Control 1")
    ),
    list(
        bigwig_directory = "data/bigwig",
        chromosome = "I",
        output_path = "output/test.svg"
    )
)
