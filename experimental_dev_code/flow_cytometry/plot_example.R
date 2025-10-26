#!/usr/bin/env Rscript
# Configure output
output_file <- "fluorescence_analysis.pdf"
pdf(output_file, width = 14, height = 10) # Professional A4 aspect ratio

# Generate synthetic data
set.seed(2025)
samples <- 3
timepoints <- 5
data <- expand.grid(
  Sample = paste0("S-", LETTERS[1:samples]),
  Timepoint = factor(paste0("T", 0:(timepoints - 1)),
    levels = paste0("T", (timepoints - 1):0)
  ),
  Fluorescence = rnorm(1500, mean = seq(50, 250, length.out = timepoints))
)

# Set up plotting parameters
original_par <- par(no.readonly = TRUE)
on.exit(par(original_par))

par(
  mfcol = c(timepoints, samples),
  oma = c(4, 5, 6, 2),
  mar = c(1.2, 1.2, 0.5, 0.5),
  mgp = c(1.5, 0.4, 0),
  tcl = -0.25
)

# Calculate global parameters
global_breaks <- seq(min(data$Fluorescence), max(data$Fluorescence), length.out = 30)
max_count <- max(hist(data$Fluorescence, breaks = global_breaks, plot = FALSE)$counts)

# Generate plots
for (tp in levels(data$Timepoint)) {
  for (smpl in unique(data$Sample)) {
    subset_data <- data[data$Sample == smpl & data$Timepoint == tp, ]

    h <- hist(subset_data$Fluorescence, breaks = global_breaks, plot = FALSE)

    plot(h,
      col = hcl.colors(length(h$breaks), "TealGrn"),
      main = "", xlab = "", ylab = "",
      ylim = c(0, max_count * 1.1),
      axes = FALSE
    )
    box(col = "gray40")

    # Add mean/variance annotations
    mean_val <- mean(subset_data$Fluorescence)
    sd_val <- sd(subset_data$Fluorescence)
    text(
      x = max(global_breaks) * 0.9, y = max_count * 0.9,
      labels = sprintf("a=%.1f\nN=%.1f", mean_val, sd_val),
      cex = 0.7, adj = 1
    )

    # Conditional axes
    if (smpl == unique(data$Sample)[1]) axis(2, cex.axis = 0.8)
    if (tp == levels(data$Timepoint)[timepoints]) axis(1, cex.axis = 0.8)
  }
}

# Add professional annotations
mtext("Fluorescence Intensity (AU)", side = 1, outer = TRUE, line = 2.5, font = 2)
mtext("Frequency", side = 2, outer = TRUE, line = 3, font = 2)
mtext(unique(data$Sample),
  side = 3, outer = TRUE,
  at = seq(1 / (2 * samples), 1, length.out = samples),
  line = -2, font = 2
)
mtext(levels(data$Timepoint),
  side = 4, outer = TRUE,
  at = seq(1 / (2 * timepoints), 1, length.out = timepoints),
  line = 1, las = 2, font = 2
)

# Add titles
title("Comprehensive Fluorescence Distribution Analysis\nExperimental Time Course (Q1 2025)",
  outer = TRUE, line = 2, cex.main = 1.5
)

# Close device
dev.off()
