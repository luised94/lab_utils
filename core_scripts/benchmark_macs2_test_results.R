# Verify required libraries are available.
required_packages <- c("GenomeInfoDb", "IRanges", "GenomicRanges", "rtracklayer", "ggplot2", "Gviz", "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
# add cosmo after manually installing
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
}

# File paths and configuration
INPUT_FILES <- list(
    narrow_peaks = "~/macs2_test_results/macs2_test_241010Bel_peaks.narrowPeak",
    broad_peaks = "~/macs2_test_results/macs2_test_100303Bel_peaks.broadPeak",
    narrow_signal = "~/macs2_test_results/macs2_test_241010Bel_treat_pileup.bdg",
    broad_signal = "~/macs2_test_results/macs2_test_100303Bel_treat_pileup.bdg",
    reference_peaks = list.files(
        file.path(Sys.getenv("HOME"), "data", "feature_files"),
        pattern = "eaton_peaks",
        full.names = TRUE
    )[1]
)
output_dir <- "~/macs2_test_results"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

# Validate file existence
for (file in INPUT_FILES) {
    if (!file.exists(file)) {
        stop(sprintf("File not found: %s", file))
    }
}

# 1. Import data
message("Importing peak files...")
narrow_peaks <- rtracklayer::import(INPUT_FILES$narrow_peaks)
broad_peaks <- rtracklayer::import(INPUT_FILES$broad_peaks)
reference_peaks <- rtracklayer::import(INPUT_FILES$reference_peaks)
# Import signal tracks for visualization
message("Importing signal tracks...")
treat_signal1 <- rtracklayer::import(INPUT_FILES$narrow_signal, format = "BedGraph")
treat_signal2 <- rtracklayer::import(INPUT_FILES$broad_signal, format = "BedGraph")

if (!is.null(reference_peaks)) {
    # Convert to chrRoman format
    GenomeInfoDb::seqlevels(reference_peaks) <- paste0(
        "chr",
        utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(reference_peaks)))
    )
}

# Find a reference origin on chrX for focused view
chromosome_to_plot <- "chrX"
chrX_refs <- reference_peaks[GenomicRanges::seqnames(reference_peaks) == chromosome_to_plot]
example_origin <- chrX_refs[1]
origin_start <- example_origin@ranges@start
origin_end <- origin_start + example_origin@ranges@width

# Expand region around origin for visualization
view_window <- 5000  # 5kb window
region_start <- max(0, origin_start - view_window)
region_end <- origin_end + view_window

# 2. Basic Statistics
message("Calculating basic statistics...")
peak_stats <- data.frame(
    Dataset = c("With_Control", "No_Control"),
    Total_Peaks = c(length(narrow_peaks), length(broad_peaks)),
    Mean_Width = c(mean(narrow_peaks@ranges@width), mean(broad_peaks@ranges@width)),
    Median_Score = c(median(narrow_peaks$score), median(broad_peaks$score))
)

# Chromosome distribution
chr_stats_narrow <- table(GenomicRanges::seqnames(narrow_peaks))
chr_stats_broad <- table(GenomicRanges::seqnames(broad_peaks))

# 4. Create visualization tracks
message("Creating visualization tracks...")
# Create tracks for narrow peaks (with control)
narrow_tracks <- list(
    Gviz::GenomeAxisTrack(),
    Gviz::DataTrack(
        treat_signal1[GenomicRanges::seqnames(treat_signal1) == "chrX"],
        name="With Control Signal",
        type="histogram",
        col.histogram="darkblue",
        fill.histogram="lightblue"
    ),
    Gviz::AnnotationTrack(
        narrow_peaks[GenomicRanges::seqnames(narrow_peaks) == "chrX"],
        name="Narrow Peaks",
        col="darkblue"
    ),
    Gviz::AnnotationTrack(
        reference_peaks[GenomicRanges::seqnames(reference_peaks) == "chrX"],
        name="Reference Origins",
        col="darkred"
    )
)

# Create tracks for broad peaks (no control)
broad_tracks <- list(
    Gviz::GenomeAxisTrack(),
    Gviz::DataTrack(
        treat_signal2[GenomicRanges::seqnames(treat_signal2) == "chrX"],
        name="No Control Signal",
        type="histogram",
        col.histogram="darkgreen",
        fill.histogram="lightgreen"
    ),
    Gviz::AnnotationTrack(
        broad_peaks[GenomicRanges::seqnames(broad_peaks) == "chrX"],
        name="Broad Peaks",
        col="darkgreen"
    ),
    Gviz::AnnotationTrack(
        reference_peaks[GenomicRanges::seqnames(reference_peaks) == "chrX"],
        name="Reference Origins",
        col="darkred"
    )
)

# 5. Generate plots
message("Generating plots...")

# Plot and save focused views
pdf(file.path(output_dir, "narrow_peaks_origin_region.pdf"), 
    width=10, height=6)
Gviz::plotTracks(
    narrow_tracks,
    chromosome="chrX",
    from=region_start,
    to=region_end,
    main="Narrow Peaks - Origin Region"
)
dev.off()

pdf(file.path(output_dir, "broad_peaks_origin_region.pdf"), 
    width=10, height=6)
Gviz::plotTracks(
    broad_tracks,
    chromosome="chrX",
    from=region_start,
    to=region_end,
    main="Broad Peaks - Origin Region"
)
dev.off()

# Plot and save whole chromosome views
pdf(file.path(output_dir, "narrow_peaks_chrX.pdf"), 
    width=15, height=6)
Gviz::plotTracks(
    narrow_tracks,
    chromosome="chrX",
    main="Narrow Peaks - Chromosome X"
)
dev.off()

pdf(file.path(output_dir, "broad_peaks_chrX.pdf"), 
    width=15, height=6)
Gviz::plotTracks(
    broad_tracks,
    chromosome="chrX",
    main="Broad Peaks - Chromosome X"
)
dev.off()

# Print region coordinates for reference
message(sprintf("Origin region plotted - chrX:%d-%d", region_start, region_end))
message(sprintf("Full chromosome plotted - chrX"))

# 6. Overlap Analysis
message("Calculating overlaps...")
narrow_overlaps <- GenomicRanges::findOverlaps(narrow_peaks, reference_peaks)
broad_overlaps <- GenomicRanges::findOverlaps(broad_peaks, reference_peaks)

overlap_stats <- data.frame(
    Dataset = c("With_Control", "No_Control"),
    Total_Peaks = c(length(narrow_peaks), length(broad_peaks)),
    Reference_Peaks = c(length(reference_peaks), length(reference_peaks)),
    Overlapping_Peaks = c(
        length(unique(S4Vectors::queryHits(narrow_overlaps))),
        length(unique(S4Vectors::queryHits(broad_overlaps)))
    ),
    Percent_Overlap = c(
        100 * length(unique(S4Vectors::queryHits(narrow_overlaps))) / length(narrow_peaks),
        100 * length(unique(S4Vectors::queryHits(broad_overlaps))) / length(broad_peaks)
    )
)

# Print results
message("Results summary:")
print(peak_stats)
print(chr_stats_narrow)
print(chr_stats_broad)
print(overlap_stats)
