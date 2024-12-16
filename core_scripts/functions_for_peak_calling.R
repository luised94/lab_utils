#' @title Find the elbow point using geometric approach
#'
#' This function identifies the elbow point in a q-value vs peak count plot
#' using the geometric distance to a line connecting the start and end points.
#'
#' @param qvals A numeric vector of q-values.
#' @return A list containing:
#'   * elbow_qval: The q-value at the elbow point.
#'   * plot_data: A data.frame with q-values, peak counts, and distances
#'     for plotting the elbow curve.
#' @export
find_elbow_point <- function(qvals) {
  stopifnot("qvals must be numeric" = is.numeric(qvals))
  stopifnot("qvals must be between 0 and 1" = all(qvals >= 0 & qvals <= 1, na.rm = TRUE))

  sorted_qvals <- sort(qvals)
  peak_counts <- seq_along(sorted_qvals)

  # Handle cases where all q-values are the same (avoid division by zero)
  if (min(sorted_qvals) == max(sorted_qvals)) {
    return(list(elbow_qval = sorted_qvals[1], plot_data = data.frame(qval = sorted_qvals, peaks = peak_counts, distances = 0)))
  }

  x <- -log10(sorted_qvals)
  x <- (x - min(x)) / (max(x) - min(x))
  y <- peak_counts / max(peak_counts)

  n <- length(x)
  line_vec <- c(x[n] - x[1], y[n] - y[1]) / sqrt(sum(line_vec^2))
  vec_from_first <- cbind(x - x[1], y - y[1])
  distances <- abs(vec_from_first %*% line_vec)

  elbow_idx <- which.max(distances)
  elbow_qval <- sorted_qvals[elbow_idx]

  return(list(
    elbow_qval = elbow_qval,
    plot_data = data.frame(
      qval = sorted_qvals,
      peaks = peak_counts,
      distances = distances
    )
  ))
}

#' @title Check the stability of the elbow point
#'
#' This function checks the stability of the elbow point by recalculating it
#' with slightly smaller and larger datasets (by removing or adding 'window' number of data points).
#'
#' @param qvals A numeric vector of q-values.
#' @param window The number of data points to remove or add for the stability check. Defaults to 10.
#' @return A numeric vector containing the elbow q-values for the lower, original and upper window sizes.
#' @export
stability_check <- function(qvals, window = 10) {
    stopifnot("qvals must be numeric" = is.numeric(qvals))
    stopifnot("Window must be a positive integer" = is.numeric(window) && window > 0 && window == as.integer(window))
    if (length(qvals) <= 2*window){
        stop("Window size is too large relative to the size of qvals")
    }
  original_elbow <- find_elbow_point(qvals)$elbow_qval
  lower_elbow <- find_elbow_point(qvals[1:(length(qvals) - window)])$elbow_qval
  upper_elbow <- find_elbow_point(qvals[1:min(length(qvals) + window, length(qvals)) ])$elbow_qval #Prevents exceeding vector length

  return(c(lower = lower_elbow, original = original_elbow, upper = upper_elbow))
}

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
