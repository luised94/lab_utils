################################################################################
# MACS2 Peak File Column Definitions
# Author: Claude3.7,Luis | Date: 2025-04-09 | Version: 1.0.0
################################################################################
# PURPOSE: File contains standardized column names for various peak file formats produced by the MACS2 peak caller.
#
# USAGE: source("peak_file_columns.R")
#
# DEPENDENCIES: NONE
#
# OUTPUTS: NONE
################################################################################
PEAK_FILE_COLUMNS <- list(
  # Standard BED file columns (BED6)
  BED = c(
    "chromosome",
    "start",
    "end",
    "name",
    "score"
    "strand"
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
    "pvalue",
    "qvalue",
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
    "pvalue",
    "qvalue"
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
    "pvalue",
    "qvalue"
  ),
  # XLS MACS2 output
  XLS = c(
    "chromosome",
    "start",
    "end",
    "length",
    "pileup",
    "pvalue",
    "fold_enrichment",
    "qvalue",
    "name"
  ),
  # XLS OUTPUT when --call-summits is enabled
  XLS_SUMMITS = c(
    "chromosome",
    "start",
    "end",
    "length",
    "abs_summit",
    "pileup",
    "pvalue",
    "fold_enrichment",
    "qvalue",
    "name"
  )
)
