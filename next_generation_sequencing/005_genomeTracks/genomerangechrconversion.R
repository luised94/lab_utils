#STATUS: REMOVE.

suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(Gviz)
    library(rtracklayer)
    library(ShortRead)
    library(tidyverse)
})
create_chromosome_GRange <- function(refGenome, chromosome_to_plot = 10) {
    cat("Creating chromosome GRange for loading feature, samples, etc\n")
    genomeRange_to_get <- GRanges(seqnames = "chrX",
                                  ranges = IRanges(start = 1, 
                                                   end = 170000,
                                  strand = "*"))
    return(genomeRange_to_get)
}

chr_convert_gr <- function(gr, verbose = FALSE) {
      if (!is(gr, "GRanges")) {
            stop("Input must be a GenomicRanges object, you dimwit!")
          }
      seqnames <- as.character(seqnames(gr))
    
      if (verbose) cat("Original seqnames:", paste(unique(seqnames), collapse = ", "), "\n")
      # Function to convert a single chromosome name
      convert_single_chr <- function(chr) {
            chr <- gsub("^chr", "", chr, ignore.case = TRUE)
            if (chr %in% c("X", "Y", "MT", "M")) {
                  return(paste0("chr", chr))
                }
            if (grepl("^[IVXLCDM]+$", chr)) {
                  return(paste0("chr", chr))  # Already in roman numeral format
                } else if (grepl("^\\d+$", chr)) {
                  return(paste0("chr", as.roman(as.integer(chr))))
                } else {
                  warning(paste("Unable to convert chromosome:", chr))
                  return(paste0("chr", chr))
                }
          }
      # Vectorized conversion
      new_seqnames <- vapply(seqnames, convert_single_chr, character(1))
      if (verbose) {
            cat("Converted seqnames:", paste(unique(new_seqnames), collapse = ", "), "\n")
            cat("Number of seqnames changed:", sum(new_seqnames != seqnames), "\n")
          }
      # Update the GenomicRanges object
      seqlevels(gr) <- unique(new_seqnames)
      seqnames(gr) <- new_seqnames
      return(gr)
}

# Example usage and testing
set.seed(42)  # For reproducibility, you numbskull
test_gr <- GRanges(
      seqnames = c("1", "X", "chrII", "20", "chrY", "M", "chr10", "III", "chrMT"),
      ranges = IRanges(start = sample(1:1000, 9), width = 100)
    )

result <- chr_convert_gr(test_gr, verbose = TRUE)
print(result)

# Test with invalid input
tryCatch(
      chr_convert_gr(data.frame()),
      error = function(e) cat("Error caught:", conditionMessage(e), "\n")
    )

main <- function() {
    genomeRange_to_get <-  create_chromosome_GRange()
    chr_convert_gr(genomeRange_to_get)
    }
main()
