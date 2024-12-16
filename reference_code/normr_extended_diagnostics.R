################################################################################
# REFERENCE CODE: Extended NormR Diagnostics and Statistics
################################################################################
# PURPOSE: Additional diagnostic metrics and detailed statistics for normR peak
#          calling. Currently not in use but potentially useful for debugging
#          or extended analysis.
# AUTHOR: [Your Name]
# DATE: [Current Date]
# DEPENDENCIES: 
# - normr
# - GenomicRanges
# - safe_write_file function from core utilities
# SOURCE: ~/lab_utils/core_scripts/functions_for_file_operations.R
################################################################################

#' Calculate comprehensive statistics from normR results
#' @param enrichment_results Output from normr::enrichR()
#' @param verbose Logical indicating whether to print processing details
#' @return List of calculated statistics
calculate_extended_diagnostics <- function(enrichment_results, verbose = FALSE) {
    # Input validation
    stopifnot(
        "enrichment_results must be a NormRFit object" = inherits(enrichment_results, "NormRFit"),
        "verbose must be logical" = is.logical(verbose) && length(verbose) == 1
    )
    
    stats <- list(
        counts = normr::getCounts(enrichment_results),
        qvals = normr::getQvalues(enrichment_results),
        enrichment_scores = normr::getEnrichment(enrichment_results),
        peak_classes = normr::getClasses(enrichment_results)
    )
    
    # Validate retrieved data
    stopifnot(
        "Failed to retrieve counts" = !is.null(stats$counts),
        "Failed to retrieve q-values" = !is.null(stats$qvals),
        "Failed to retrieve enrichment scores" = !is.null(stats$enrichment_scores),
        "Failed to retrieve peak classes" = !is.null(stats$peak_classes)
    )
    
    # Calculate derived statistics
    stats$total_bins <- length(stats$counts$treatment)
    stats$zero_bins <- sum(stats$counts$treatment == 0 & 
                          stats$counts$control == 0)
    stats$enrichment_by_threshold <- sapply(
        NORMR_CONFIG$fdr_thresholds, 
        function(threshold) {
            sum(!is.na(stats$qvals) & stats$qvals <= threshold)
        }
    )
    stats$sig_regions <- !is.na(stats$qvals) & 
                        stats$qvals <= NORMR_CONFIG$default_fdr
    
    return(stats)
}

#' Export detailed results using safe file writing
#' @param enrichment_results normr enrichment results
#' @param peak_dir Output directory
#' @param chip_id ChIP sample identifier
#' @param input_id Input sample identifier
#' @param verbose Logical indicating whether to print processing details
#' @param interactive Logical indicating whether to prompt for file overwrites
export_detailed_results <- function(enrichment_results, peak_dir, chip_id, input_id,
                                  verbose = FALSE, interactive = TRUE) {
    # Input validation
    stopifnot(
        "enrichment_results must be a NormRFit object" = inherits(enrichment_results, "NormRFit"),
        "peak_dir must exist" = dir.exists(peak_dir),
        "chip_id must be character" = is.character(chip_id) && length(chip_id) == 1,
        "input_id must be character" = is.character(input_id) && length(input_id) == 1,
        "verbose must be logical" = is.logical(verbose) && length(verbose) == 1,
        "interactive must be logical" = is.logical(interactive) && length(interactive) == 1
    )
    
    # Create region information
    ranges <- normr::getRanges(enrichment_results)
    counts <- normr::getCounts(enrichment_results)
    
    region_info <- data.frame(
        chromosome = seqnames(ranges),
        start = start(ranges),
        end = end(ranges),
        treatment_count = counts$treatment,
        control_count = counts$control,
        enrichment = normr::getEnrichment(enrichment_results),
        qvalue = normr::getQvalues(enrichment_results),
        peak_class = normr::getClasses(enrichment_results)
    )
    
    # Export using safe_write_file
    regions_output <- file.path(
        peak_dir,
        sprintf(
            NORMR_CONFIG$region_file_template,
            TIMESTAMPS$full,
            chip_id,
            input_id,
            "peaks"
        )
    )
    
    safe_write_file(
        data = region_info,
        path = regions_output,
        write_fn = write.table,
        verbose = verbose,
        interactive = interactive,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}

################################################################################
# USAGE EXAMPLES
################################################################################
# Example usage with safe_write_file:
#
# results <- export_detailed_results(
#     enrichment_results = enrichment_results,
#     peak_dir = peak_dir,
#     chip_id = chip_id,
#     input_id = input_id,
#     verbose = DEBUG_CONFIG$verbose,
#     interactive = DEBUG_CONFIG$interactive
# )
#
# if (!results) {
#     warning("Failed to export detailed results")
# }
################################################################################
#######

# Initialize data.frame for gathering statistics.
#sample_statistics <- data.frame(
#    timestamp = character(),
#    chip_id = character(),
#    input_id = character(),
#    total_bins = integer(),
#    zero_bins = integer(),
#    zero_bins_percent = numeric(),
#    mean_treatment_counts = numeric(),
#    mean_control_counts = numeric(),
#    strong_peaks = integer(),  # q <= 0.01
#    medium_peaks = integer(),  # q <= 0.05
#    weak_peaks = integer(),    # q <= 0.10
#    mean_enrichment_score = numeric(),
#    enriched_regions = integer(),
#    within_expected_range = logical(),
#    stringsAsFactors = FALSE
#)
#        if (DEBUG_CONFIG$verbose) {
#            # Collect all statistics
#            stats <- list(
#                # Count statistics
#                counts = normr::getCounts(enrichment_results),
#                qvals = normr::getQvalues(enrichment_results),
#                enrichment_scores = normr::getEnrichment(enrichment_results),
#                peak_classes = normr::getClasses(enrichment_results)
#            )
#
#            # Calculate derived statistics
#            stats$total_bins <- length(stats$counts$treatment)
#            stats$zero_bins <- sum(stats$counts$treatment == 0 & stats$counts$control == 0)
#            stats$enrichment_by_threshold <- sapply(NORMR_CONFIG$fdr_thresholds, function(threshold) {
#                sum(!is.na(stats$qvals) & stats$qvals <= threshold)
#            })
#            stats$sig_regions <- !is.na(stats$qvals) & stats$qvals <= NORMR_CONFIG$default_fdr
#
#            # Print comprehensive analysis report
#            message("\n=== normR Peak Calling Analysis Report ===")
#            message("\n1. Sample Information:")
#            message(sprintf("   Treatment: %s", chip_id))
#            message(sprintf("   Control: %s", input_id))
#
#            message("\n2. Count Statistics:")
#            message(sprintf("   Total bins analyzed: %d", stats$total_bins))
#            message(sprintf("   Empty bins (zero in both): %d (%.2f%%)",
#                    stats$zero_bins,
#                    100 * stats$zero_bins/stats$total_bins))
#            message(sprintf("   Mean treatment counts: %.2f", mean(stats$counts$treatment)))
#            message(sprintf("   Mean control counts: %.2f", mean(stats$counts$control)))
#
#            message("\n3. Peak Statistics:")
#            message(sprintf("   Strong peaks (q <= 1%%): %d", stats$enrichment_by_threshold[1]))
#            message(sprintf("   Medium peaks (q <= 5%%): %d", stats$enrichment_by_threshold[2]))
#            message(sprintf("   Weak peaks (q <= 10%%): %d", stats$enrichment_by_threshold[3]))
#
#            if (any(stats$sig_regions)) {
#                message("\n4. Enrichment Statistics:")
#                message(sprintf("   Mean enrichment score (q <= 5%%): %.2f",
#                        mean(stats$enrichment_scores[stats$sig_regions])))
#                enriched_count <- sum(!is.na(stats$peak_classes) & stats$peak_classes == 1)
#                message(sprintf("   Regions classified as enriched: %d", enriched_count))
#            }
#
#            message("\n5. Biological Context:")
#            message("   Expected: 250-400 ORC binding sites in S. cerevisiae")
#            if (stats$enrichment_by_threshold[2] < NORMR_CONFIG$expected_peak_range$min) {
#                message("   WARNING: Fewer peaks than expected for ORC binding")
#            } else if (stats$enrichment_by_threshold[2] > NORMR_CONFIG$expected_peak_range$max) {
#                message("   WARNING: More peaks than expected for ORC binding")
#            } else {
#                message("   Peak count is within expected range for ORC binding")
#            }
#
#            message("\n=========================================\n")
#        }
#        if (!DEBUG_CONFIG$dry_run) {
#            # Get genomic ranges from enrichment results
#            ranges <- normr::getRanges(enrichment_results)
#            counts <- normr::getCounts(enrichment_results)
#
#            # Create detailed region information
#            region_info <- data.frame(
#                chromosome = seqnames(ranges),
#                start = start(ranges),
#                end = end(ranges),
#                treatment_count = counts$treatment,
#                control_count = counts$control,
#                enrichment = normr::getEnrichment(enrichment_results),
#                qvalue = normr::getQvalues(enrichment_results),
#                peak_class = normr::getClasses(enrichment_results)
#            )
#
#            # Export full region information for later analysis
#            regions_output <- file.path(
#                peak_dir,
#                sprintf(
#                    NORMR_CONFIG$region_file_template,
#                    TIMESTAMPS$full,
#                    sample_id_mapping[chip_id],
#                    sample_id_mapping[input_id],
#                    "normr_peaks"
#                )
#            )
#
#            write.table(
#                region_info,
#                file = regions_output,
#                sep = "\t",
#                quote = FALSE,
#                row.names = FALSE
#            )
#
#            # Export significant peaks in bedGraph format for genome browser
#            sig_peaks <- region_info[!is.na(region_info$peak_class) &
#                                   region_info$qvalue <= NORMR_CONFIG$default_fdr, ]
#
#            bedgraph_output <- file.path(
#                peak_dir,
#                sprintf(
#                    NORMR_CONFIG$bedgraph_template,
#                    TIMESTAMPS$full,
#                    sample_id_mapping[chip_id],
#                    sample_id_mapping[input_id],
#                    "normr_peaks"
#                )
#            )
#
#            write.table(
#                sig_peaks[, c("chromosome", "start", "end", "enrichment")],
#                file = bedgraph_output,
#                sep = "\t",
#                quote = FALSE,
#                row.names = FALSE,
#                col.names = FALSE
#            )
#
#            # Export enrichment results in BED format
#            normr::exportR(
#                obj = enrichment_results,
#                filename = output_path,
#                type = "bed"
#            )
#        }
#
#        current_stats <- data.frame(
#            timestamp = TIMESTAMPS$full,
#            chip_id = chip_id,
#            input_id = input_id,
#            total_bins = stats$total_bins,
#            zero_bins = stats$zero_bins,
#            zero_bins_percent = 100 * stats$zero_bins/stats$total_bins,
#            mean_treatment_counts = mean(stats$counts$treatment),
#            mean_control_counts = mean(stats$counts$control),
#            strong_peaks = stats$enrichment_by_threshold[1],
#            medium_peaks = stats$enrichment_by_threshold[2],
#            weak_peaks = stats$enrichment_by_threshold[3],
#            mean_enrichment_score = if(any(stats$sig_regions)) mean(stats$enrichment_scores[stats$sig_regions]) else NA,
#            enriched_regions = sum(!is.na(stats$peak_classes) & stats$peak_classes == 1),
#            within_expected_range = stats$enrichment_by_threshold[2] >= NORMR_CONFIG$expected_peak_range$min &&
#                                  stats$enrichment_by_threshold[2] <= NORMR_CONFIG$expected_peak_range$max
#        )
#
#        sample_statistics <- rbind(sample_statistics, current_stats)
#if (DEBUG_CONFIG$verbose && nrow(sample_statistics) > 0) {
#    message("\n=== Final Analysis Summary ===")
#
#    # Basic statistics
#    message("\n1. Overall Processing Statistics:")
#    message(sprintf("   Total samples processed: %d", nrow(sample_statistics)))
#    message(sprintf("   Samples within expected peak range: %d of %d",
#                   sum(sample_statistics$within_expected_range),
#                   nrow(sample_statistics)))
#
#    # Peak statistics
#    message("\n2. Peak Statistics Summary:")
#    message("   Strong peaks (q <= 1%):")
#    message(sprintf("     Mean: %.1f", mean(sample_statistics$strong_peaks)))
#    message(sprintf("     Range: %d - %d", min(sample_statistics$strong_peaks), max(sample_statistics$strong_peaks)))
#    message("   Medium peaks (q <= 5%):")
#    message(sprintf("     Mean: %.1f", mean(sample_statistics$medium_peaks)))
#    message(sprintf("     Range: %d - %d", min(sample_statistics$medium_peaks), max(sample_statistics$medium_peaks)))
#
#    # Count statistics
#    message("\n3. Count Statistics Summary:")
#    message(sprintf("   Mean zero bins percentage: %.2f%%", mean(sample_statistics$zero_bins_percent)))
#    message(sprintf("   Mean treatment counts: %.2f", mean(sample_statistics$mean_treatment_counts)))
#    message(sprintf("   Mean control counts: %.2f", mean(sample_statistics$mean_control_counts)))
#
#    # Enrichment statistics
#    message("\n4. Enrichment Statistics Summary:")
#    message(sprintf("   Mean enriched regions per sample: %.1f", mean(sample_statistics$enriched_regions)))
#    message(sprintf("   Mean enrichment score: %.2f", mean(sample_statistics$mean_enrichment_score, na.rm=TRUE)))
#
#    # Samples outside expected range
#    if (any(!sample_statistics$within_expected_range)) {
#        message("\n5. Samples Outside Expected Range:")
#        outliers <- sample_statistics[!sample_statistics$within_expected_range, ]
#        for (i in 1:nrow(outliers)) {
#            message(sprintf("   %s vs %s: %d peaks",
#                          outliers$chip_id[i],
#                          outliers$input_id[i],
#                          outliers$medium_peaks[i]))
#        }
#    }
#
#    message("\n===========================\n")
#
#    # Export statistics if not in dry run mode
#    if (!DEBUG_CONFIG$dry_run) {
#        stats_output <- file.path(
#            peak_dir,
#            sprintf(
#                "peak_calling_statistics_%s.tsv",
#                TIMESTAMPS$full
#            )
#        )
#        write.table(
#            sample_statistics,
#            file = stats_output,
#            sep = "\t",
#            quote = FALSE,
#            row.names = FALSE
#        )
#    }
#}
