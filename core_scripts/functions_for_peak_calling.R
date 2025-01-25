#' @title Calculate the overlap between two sets of genomic ranges
#'
#' This function calculates the percentage of overlap between two sets of genomic ranges
#' and returns a list containing the percentage overlap and the unique peaks in the first set.
#'
#' @param peaks A GRanges object representing the first set of peaks.
#' @param reference A GRanges object representing the second set of peaks (reference).
#' @return A list containing:
#'   * percent_overlap: The percentage of overlap between the two sets.
#'   * unique_peaks: A GRanges object containing the peaks in the first set that do not overlap with the reference set.
#' @import GenomicRanges
#' @export
calculate_overlap <- function(peaks, reference) {
  stopifnot("peaks must be a GRanges object" = methods::is(peaks, "GRanges"))
  stopifnot("reference must be a GRanges object" = methods::is(reference, "GRanges"))

  overlap <- GenomicRanges::findOverlaps(peaks, reference)
  percent_overlap <- length(unique(S4Vectors::queryHits(overlap))) / length(peaks) * 100
  unique_peaks <- peaks[-unique(S4Vectors::queryHits(overlap))]
  return(list(percent_overlap = percent_overlap, unique_peaks = unique_peaks))
}
