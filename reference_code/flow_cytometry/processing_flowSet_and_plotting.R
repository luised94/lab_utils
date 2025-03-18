# 1. Create working copy
filtered_set_density <- subset_flowSet
# Get global FL1-A range across all flowFrames in flowSet
fl1a_ranges <- fsApply(filtered_set_density, function(fr) range(exprs(fr)[,"FL1-A"]))
global_range <- range(fl1a_ranges)

# 2. Trim extreme 0.5% from both tails
trim_level <- 0.005  # Adjust this to control trimming strictness

for(i in seq_along(filtered_set_density)) {
  ff <- filtered_set_density[[i]]
  mat <- exprs(ff)[,"FL1-A"]
  
  # Calculate trim thresholds
  lower <- quantile(mat, trim_level)
  upper <- quantile(mat, 1-trim_level)
  
  # Apply filter
  keep <- mat > lower & mat < upper
  exprs(ff) <- exprs(ff)[keep,]
  filtered_set_density[[i]] <- ff
  
  cat("Sample", sampleNames(filtered_set_density)[i],
      "retained", sum(keep), "/", length(keep), "events\n")
}
ggcyto(filtered_set_density, aes(x = `FL1-A`)) +
  geom_histogram(bins = 50) +
  facet_wrap(~timepoints) +
  geom_vline(
    data = timepoint_medians,
    aes(xintercept = median_FL1A),
    color = "red",
    linetype = "dashed"
  ) +
  coord_cartesian(xlim = global_range) +  # Use computed range
  labs(title = "Density-Trimmed FL1-A Distributions")

# 1. Create working copy to preserve original data
filtered_set_iqr <- subset_flowSet
fl1a_ranges <- fsApply(filtered_set_iqr, function(fr) range(exprs(fr)[,"FL1-A"]))
global_range <- range(fl1a_ranges)

# 2. Process each sample using flowCore's fsApply
for(i in seq_along(filtered_set_iqr)) {
  ff <- filtered_set_iqr[[i]]
  
  # Extract fluorescence data matrix
  mat <- exprs(ff)[,"FL1-A"]
  
  # Calculate quartiles and IQR
  q25 <- quantile(mat, 0.25)
  q75 <- quantile(mat, 0.75)
  iqr <- q75 - q25
  
  # Define valid range
  lower <- q25 - 1.5*iqr
  upper <- q75 + 1.5*iqr
  
  # Filter events and update flowFrame
  keep <- mat > lower & mat < upper
  exprs(ff) <- exprs(ff)[keep,]
  filtered_set_iqr[[i]] <- ff
  
  # Show filtering impact
  cat("Sample", sampleNames(filtered_set_iqr)[i], 
      "retained", sum(keep), "/", length(keep), "events\n")
}

ggcyto(filtered_set_iqr, aes(x = `FL1-A`)) +
    geom_histogram(bins = 50) +
  geom_vline(
    data = timepoint_medians,
    aes(xintercept = median_FL1A),
    color = "red",
    linetype = "dashed"
  ) +
    facet_wrap(~timepoints) +
    coord_cartesian(xlim = global_range) +  # Use computed range
    labs(title = "IQR-Filtered FL1-A Distributions")

# 1. Create WORKING COPY of original flowSet (preserves sample names)
filtered_set <- subset_flowSet

# 2. Process each flowFrame IN PLACE
for(i in seq_along(filtered_set)){
  # Get current flowFrame data
  current_frame <- filtered_set[[i]]
  frame_data <- exprs(current_frame)
  
  # 3. Filter key channels using IQR
  for(channel in c("FSC-A", "SSC-A", "FL1-A")){
    vals <- frame_data[,channel]
    q <- quantile(vals, probs = c(0.25, 0.75), na.rm = TRUE)
    iqr <- q[2] - q[1]
    lower <- q[1] - 1.5 * iqr
    upper <- q[2] + 1.5 * iqr
    
    # Apply filter to THIS channel
    frame_data <- frame_data[vals >= lower & vals <= upper, ]
  }
  
  # 4. UPDATE IN PLACE (preserves sample names)
  exprs(current_frame) <- frame_data
  filtered_set[[i]] <- current_frame
  
  cat("Sample", sampleNames(filtered_set)[i], "now has", nrow(frame_data), "events\n")
}

# 5. Plot (no phenoData reassignment needed)
ggcyto(filtered_set, aes(x = `FL1-A`)) +
  geom_histogram(bins = 50) +
  facet_wrap(~timepoints)

# 1. Create WORKING COPY
filtered_set <- subset_flowSet

# 2. Process each flowFrame IN PLACE
for(i in seq_along(filtered_set)){
  current_frame <- filtered_set[[i]]
  frame_data <- exprs(current_frame)
  
  # 3. Hampel filter for key channels
  for(channel in c("FSC-A", "SSC-A", "FL1-A")){
    vals <- frame_data[,channel]
    med <- median(vals, na.rm = TRUE)
    mad_val <- mad(vals, constant = 1, na.rm = TRUE)
    
    lower <- med - 3 * mad_val
    upper <- med + 3 * mad_val
    
    frame_data <- frame_data[vals >= lower & vals <= upper, ]
  }
  
  # 4. Update IN PLACE
  exprs(current_frame) <- frame_data
  filtered_set[[i]] <- current_frame
  
  cat("Sample", sampleNames(filtered_set)[i], "retained", nrow(frame_data), "events\n")
}

# 5. Plot directly
ggcyto(filtered_set, aes(x = `FSC-A`, y = `SSC-A`)) +
  geom_hex(bins = 100) +
  facet_wrap(~timepoints)
