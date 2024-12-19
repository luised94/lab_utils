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
    narrow_peaks = "macs2_test_results/macs2_test_241010Bel_peaks.narrowPeak",
    broad_peaks = "macs2_test_results/macs2_test_100303Bel_peaks.broadPeak",
    narrow_signal = "macs2_test_results/macs2_test_241010Bel_treat_pileup.bdg",
    broad_signal = "macs2_test_results/macs2_test_100303Bel_treat_pileup.bdg",
    feature_file = list.files(
        file.path(Sys.getenv("HOME"), "data", "feature_files"),
        pattern = "eaton_peaks",
        full.names = TRUE
    )[1]
)

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
reference_peaks <- rtracklayer::import(INPUT_FILES$feature_file)

# 2. Basic Statistics
message("Calculating basic statistics...")
peak_stats <- data.frame(
    Dataset = c("With_Control", "No_Control"),
    Total_Peaks = c(length(narrow_peaks), length(broad_peaks)),
    Mean_Width = c(mean(width(narrow_peaks)), mean(width(broad_peaks))),
    Median_Score = c(median(narrow_peaks$score), median(broad_peaks$score))
)

# Chromosome distribution
chr_stats_narrow <- table(seqnames(narrow_peaks))
chr_stats_broad <- table(seqnames(broad_peaks))

# 3. Import signal tracks for visualization
message("Importing signal tracks...")
treat_signal1 <- rtracklayer::import(INPUT_FILES$narrow_signal)
treat_signal2 <- rtracklayer::import(INPUT_FILES$broad_signal)

# 4. Create visualization tracks
message("Creating visualization tracks...")
gtrack <- Gviz::GenomeAxisTrack()

dtrack1 <- Gviz::DataTrack(
    treat_signal1[seqnames(treat_signal1) == "chrII"],
    name="With Control Signal"
)

dtrack2 <- Gviz::DataTrack(
    treat_signal2[seqnames(treat_signal2) == "chrII"],
    name="No Control Signal"
)

ptrack1 <- Gviz::AnnotationTrack(
    narrow_peaks,
    name="Narrow Peaks"
)

ptrack2 <- Gviz::AnnotationTrack(
    broad_peaks,
    name="Broad Peaks"
)

rtrack <- Gviz::AnnotationTrack(
    reference_peaks,
    name="Reference Peaks"
)

# 5. Generate plots
message("Generating plots...")
Gviz::plotTracks(
    list(gtrack, dtrack1, ptrack1, dtrack2, ptrack2, rtrack),
    chromosome="chrII",
    from=200000,
    to=250000
)

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
