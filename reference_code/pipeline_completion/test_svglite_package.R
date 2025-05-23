#!/usr/bin/env Rscript
###############################################################################
# svglite functionality test script
################################################################################
# PURPOSE: Create some simple plots to ensure svglite package works in the luria cluster environment
# Conclusion: Seems to work in luria. Not errors. Both plots produced succesfully. Opened and had expected content.
# USAGE: 
#   = source("reference_code/pipeline_completion/test_svglite_package.R")
#   = just run directly line by line
# DEPENDENCIES: svglite
# OUTPUT: two svg files in the svglite_test
# AUTHOR: claude 3.7 sonnet
# DATE: 2025-04-08
################################################################################
# 
# Save this as test_svglite.R and run with: Rscript test_svglite.R

cat("Testing svglite functionality...\n")

# Check if svglite is installed
if (!requireNamespace("svglite", quietly = TRUE)) {
  cat("ERROR: svglite package is not installed.\n")
  cat("Install with: install.packages('svglite')\n")
  quit(status = 1)
}

# Load the package
library(svglite)
cat("svglite package loaded successfully.\n")

# Test system dependencies
tryCatch({
  # Create a simple test directory
  test_dir <- "svglite_test"
  dir.create(test_dir, showWarnings = FALSE)
  test_file <- file.path(test_dir, "test_plot.svg")

  cat(sprintf("Attempting to create SVG file: %s\n", test_file))

  # Try creating a simple plot
  svglite(test_file, width = 5, height = 5)
  plot(1:10, main = "Test Plot", xlab = "X Axis", ylab = "Y Axis")
  text(5, 5, "Text rendering test")
  dev.off()

  # Check if file was created with reasonable content
  if (file.exists(test_file)) {
    file_size <- file.info(test_file)$size
    cat(sprintf("SUCCESS: SVG file created, size: %d bytes\n", file_size))

    # Read first few lines to verify it's a valid SVG
    svg_content <- readLines(test_file, n = 5)
    if (any(grepl("<svg", svg_content, fixed = TRUE))) {
      cat("File appears to be a valid SVG.\n")
    } else {
      cat("WARNING: File may not be a valid SVG. Check content manually.\n")
    }

    # Print the file path for manual inspection
    cat(sprintf("File path for manual inspection: %s\n", normalizePath(test_file)))
  } else {
    cat("ERROR: Failed to create SVG file.\n")
  }

  # Test with different parameters that might be used in your actual script
  test_file2 <- file.path(test_dir, "test_plot_large.svg")
  cat(sprintf("Testing with larger dimensions: %s\n", test_file2))

  svglite(test_file2, width = 10, height = 8, pointsize = 12, standalone = TRUE)
  plot(1:100, type = "l", col = "blue", lwd = 2)
  points(1:100, pch = 19, col = "red", cex = 0.5)
  title("Complex Plot Test")
  dev.off()

  if (file.exists(test_file2)) {
    cat("SUCCESS: Complex SVG file created.\n")
  } else {
    cat("ERROR: Failed to create complex SVG file.\n")
  }

  # Overall result
  cat("\nOverall svglite test: PASSED\n")

}, error = function(e) {
  cat(sprintf("ERROR: svglite test failed with error:\n%s\n", e$message))
  quit(status = 1)
})

cat("All tests completed.\n")
