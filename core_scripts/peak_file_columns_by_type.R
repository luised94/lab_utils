#' MACS2 Peak File Column Definitions
#'
#' @description
#' This file contains standardized column definitions for various peak file formats
#' produced by MACS2 peak caller, including narrowPeak, broadPeak, gappedPeak,
#' standard BED, and MACS2 XLS output formats. These definitions facilitate consistent
#' parsing and manipulation of MACS2 output files in R analysis pipelines.
#'
#' @author claude3.7, LEMR
#' @date April 9, 2025
#'
#' @details
#' This script defines a list called 'PEAK_FILE_COLUMNS' that contains named vectors
#' for each file type's column structure. The column names follow standard genomics
#' conventions and can be used directly with rtracklayer's import function or
#' when reading files as data frames.
#'
#' File formats included:
#' - BED: Standard 6-column BED format
#' - BED12: Extended 12-column BED format used for features with internal structure
#' - NARROW_PEAK: BED6+4 format for narrow binding sites (MACS2 narrowPeak)
#' - BROAD_PEAK: BED6+3 format for broader regions (MACS2 broadPeak)
#' - GAPPED_PEAK: BED12+3 format for complex regions with gaps (MACS2 gappedPeak)
#' - XLS: MACS2 detailed peak information in tabular format
#'
#' @usage
#' # Source this file in your analysis scripts:
#' source("peak_file_columns.R")
#'
#' # Example usage:
#' # Reading narrowPeak files with correct column names
#' library(rtracklayer)
#' peaks <- import("peaks.narrowPeak", format="narrowPeak")
#'
#' # Or with base R:
#' peaks_df <- read.table("peaks.narrowPeak", sep="\t",
#'                       col.names=PEAK_FILE_COLUMNS$NARROW_PEAK)
#'
#' # For converting to GRanges with proper metadata:
#' library(GenomicRanges)
#' gr <- GRanges(
#'   seqnames = peaks_df$chromosome,
#'   ranges = IRanges(start = peaks_df$start + 1, end = peaks_df$end),  # convert 0-based to 1-based
#'   strand = peaks_df$strand
#' )
#' mcols(gr) <- peaks_df[, 4:ncol(peaks_df)]  # add metadata columns
#' names(mcols(gr)) <- PEAK_FILE_COLUMNS$NARROW_PEAK[4:length(PEAK_FILE_COLUMNS$NARROW_PEAK)]
#'
#' @note
#' Remember that BED format is 0-based for start coordinates, while GRanges objects
#' in Bioconductor are 1-based. Always add 1 to start positions when converting.
#'
#' @version 1.0
#'
#' @seealso
#' MACS2 documentation: https://github.com/macs3-project/MACS
#' ENCODE file formats: https://www.encodeproject.org/data-standards/file-formats/
#'
# Column definitions for different peak file types
PEAK_FILE_COLUMNS <- list(
  # Standard BED file columns (BED6)
  BED = c(
    "chromosome",
    "start",
    "end",
    "name",
    "score"
    #"strand"
  ),

  # BED12 format (used in gappedPeak)
  BED12 = c(
    "chromosome",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts"
  ),

  # narrowPeak columns (BED6+4)
  NARROWPEAK = c(
    "chromosome",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signalValue",
    "pValue",
    "qValue",
    "peak"
  ),

  # broadPeak columns (BED6+3)
  BROADPEAK = c(
    "chromosome",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signalValue",
    "pValue",
    "qValue"
  ),

  # gappedPeak columns (BED12+3)
  GAPPEDPEAK = c(
    "chromosome",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts",
    "signalValue",
    "pValue",
    "qValue"
  ),

  # XLS MACS2 output
  XLS = c(
    "chromosome",
    "start",
    "end",
    "length",
    "summit",
    "tagsCount",
    "pValue",
    "foldEnrichment",
    "qValue",
    "name"
  )
)
