library(flowCore)
library(ggcyto)
library(svglite)
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

path_to_fcs_files <- "/mnt/c/Users/Luis/Dropbox (MIT)/Lab/Experiments/Yeast Genetics/2023_10_04 Cdc6 Overexpression/CL4_Cdc6Oe_CarbonSourceShift_01"
fcs_file_paths <- list.files(path_to_fcs_files, pattern = "\\.fcs$", full.names = TRUE)
fcs_data <- read.FCS(fcs_file_paths[1])
colnames(head(exprs(fcs_data)))
head(exprs(fcs_data)[,11])
head(exprs(fcs_data[,'FITC-Width']))

output_directory <- "/mnt/c/Users/Luis/Dropbox (MIT)/Lab/Experiments/Yeast Genetics/2023_10_04 Cdc6 Overexpression/"

plot_output <- paste0(output_directory, "/", "plot_hist.svg")
svglite(plot_output, width = 6, height = 4)
hist(exprs(fcs_data[,'FITC-Width']),
	  main = "Histogram of FITC-Width",
	  xlab = "FITC-Width",
	  col = "blue",
	  border = "black")
dev.off()

plot_output <- paste0(output_directory, "/", "plot_density.svg")
svglite(plot_output, width = 6, height = 4)
ggplot(fcs_data, aes(x = `FITC-Width`)) +
	geom_histogram(aes(y = after_stat(density)), binwidth = 30, colour = "black", fill = "white") +
	geom_density(alpha = .2, fill = "#FF6666") +
	labs(title = "Fluorescence Intensity Distribution", x = "Fluorescence Intensity", y = "Density")
dev.off()

plot_output <- paste0(output_directory, "/", "plot_count.svg")
svglite(plot_output, width = 6, height = 4)
ggplot(fcs_data, aes(x = `FITC-Width`)) +
	geom_histogram(aes(y = after_stat(count)), binwidth = 30, colour = "black", fill = "white") +
	geom_density(alpha = .2, fill = "#FF6666") +
	labs(title = "Fluorescence Intensity Distribution", x = "Fluorescence Intensity", y = "Density")
dev.off()

summary(fcs_data)
