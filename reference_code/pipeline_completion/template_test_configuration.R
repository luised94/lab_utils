RUNTIME_CONFIG <- list(
    # Core control flags
    debug_interactive = FALSE,
    debug_verbose = TRUE,
    debug_validate = TRUE,

    # Processing control
    process_single_file = FALSE,
    process_single_comparison = TRUE,
    process_comparison = "comp_1108forNoneAndWT",
    process_chromosome = 14,
    #process_group = 10,
    process_batch = 10,
    process_samples_per_batch = 4,
    #process_samples_per_group = 4,
    process_file_index = 1,
    # Output control
    output_save_plots = FALSE,
    output_dry_run = TRUE,
    output_display_time = 2
)

GENOME_TRACK_CONFIG <- list(
    use_custom_visualization = FALSE,  # Control flag

    # Display dimensions
    display_width = 10,
    display_height = 8,
    
    # Track Creation
    track_points_default = 1000,
    #track_show_title = TRUE,

    # Track defaults
    track_ylim = c(0, 1000),  # Default y-limits, adjust as needed
    track_sampling_rate = 100,  # Points per base pair for empty tracks
    
    # Track colors
    color_placeholder = "#cccccc",
    color_input = "#808080",
    
    # Track naming
    format_sample_track_name = "%s: %s",
    format_control_track_name = "%s: %s - %s",
    format_placeholder_track_name = "%s: %s - %s",
    format_suffix = "(No data)",
    format_genome_axis_track_name = "Chr %s Axis",

    # Labels
    label_always_show = "antibody",
    label_never_show = c("sample_id", "full_name", "short_name", "X__cf_genotype"),
    label_separator = "-",

    # File handling
    file_pattern = "consolidated_.*_sequence\\.fastq$",
    file_sample_id = "consolidated_([0-9]{5,6})_sequence\\.fastq",
    file_sample_id_from_bigwig = "processed_([0-9]{5,6})_sequence_to_S288C_(RPKM|CPM|BPM|RPGC)\\.bw",
    file_genome_pattern = "S288C_refgenome.fna",
    file_genome_directory = file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    file_feature_directory = file.path(Sys.getenv("HOME"), "data", "feature_files"),
    file_feature_pattern = "eaton_peaks",

    # File Names
    filename_format_group_template = "%s_%s_group%02d_chr%s_%s.svg",
    filename_format_comparison_template = "%s_%s_%s_chr%s_%s.svg",
    title_group_template = paste(
        "%s",               # Title
        "Group: %s",   # Comparison ID
        "Chromosome %s (%d samples)", # Chr info
        "%s",               # Additional info
        "Normalization: %s", # Norm method
        sep = "\n"
    ),
    title_comparison_template = paste(
        "%s",               # Title
        "Comparison: %s",   # Comparison ID
        "Chromosome %s (%d samples)", # Chr info
        "%s",               # Additional info
        "Normalization: %s", # Norm method
        sep = "\n"
    ),

    track_defaults_sample = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 1.2,
        background.title = "white",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.6,
        fontface = 1,
        title.width = 1.2
    ),

    track_defaults_placeholder = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 0.8,
        background.title = "white",
        background.panel = "#f5f5f5",    # light gray to indicate "empty"
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9,
        alpha = 0.5,
        grid = FALSE
        #ylim = c(0, 1)                  # fixed range for empty tracks
    ),
    track_defaults_control = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 0.8,
        background.title = "white",
        background.panel = "white",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9
        #alpha = 0.8
    ),
    track_defaults_feature = list(
        showaxis = FALSE,
        showtitle = TRUE,
        size = 0.5,
        background.title = "white",
        background.panel = "#8b7355",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9,
        fill = "#8b4513",
        col = "#8b4513"
    ),
    
    # New Plot Defaults
    plot_defaults = list(
        margin = 15,
        innerMargin = 5,
        spacing = 10,
        extend.left = 0,
        extend.right = 0,
        col.axis = "black",
        cex.axis = 0.8,
        cex.main = 0.8,
        fontface.main = 2,
        background.panel = "transparent"
    )
)
