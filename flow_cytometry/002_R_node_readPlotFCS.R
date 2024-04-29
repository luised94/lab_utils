library(flowCore)
library(ggcyto)

# Single file
fcs_data <- read.FCS("path_to_your_file.fcs")

# Multiple files
fcs_files <- list.files(path = "path_to_fcs_files", pattern = "\\.fcs$", full.names = TRUE)
fcs_set <- read.flowSet(files = fcs_files)

# Plotting a single parameter histogram
p <- ggcyto(fcs_data, aes(x = `FL1-H`)) + geom_histogram(bins = 30)
print(p)

# Plotting a scatter plot of two parameters
p <- ggcyto(fcs_data, aes(x = `FSC-A`, y = `SSC-A`)) + geom_point()
print(p)

# Using facet_wrap with a flowSet
p <- ggcyto(fcs_set, aes(x = `FL1-H`)) + geom_histogram(bins = 30) + facet_wrap(~ name)
print(p)

# Assuming `gs` is a GatingSet with gates applied
p <- ggcyto(gs, aes(x = `FL1-H`)) + geom_histogram(data = "CD3+")  # Subset to CD3+ population
print(p)


