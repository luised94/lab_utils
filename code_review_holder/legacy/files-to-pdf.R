# Install the webshot package (if not already installed)
# install.packages("webshot")

# Load the required libraries
library(webshot)

# Set the path to the input file
input_file <- "./2021-01-18-reading-tables-from-images-with-magick.html"

# Set the path to the output file
output_file <- "./reading-tables-from-images-with-magick.pdf"

# Use the webshot function to convert the HTML file to PDF
webshot(input_file, output_file, delay = 5, vwidth = 800, vheight = 600)


# Install the epubr and xml2 packages (if not already installed)
# install.packages(c("epubr", "xml2"))

# Load the required libraries
library(epubr)
library(xml2)
library(webshot)

# Set the path to the input file
input_file <- "path/to/input/file.epub"

# Set the path to the output file
output_file <- "path/to/output/file.pdf"

# Use the epub_extract function to extract the contents of the EPUB file
epub_content <- epub_extract(input_file)

# Extract the HTML content from the EPUB file
html_content <- read_html(epub_content$content)

# Use the webshot function to convert the HTML content to PDF
webshot(html_content, output_file, delay = 5, vwidth = 800, vheight = 600)
