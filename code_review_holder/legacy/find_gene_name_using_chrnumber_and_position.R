#from statisticalrecipes blogspot, stephanie hicks https://statisticalrecipes.blogspot.com/2012/08/biomart-find-gene-name-using-chromosome.html
# Load the library
library(biomaRt)

# Define biomart object
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Gives a list of all possible annotations; Currently there are 1668 listed
listAttributes(mart)

# Gives a list of all filters or criteria to search by; Currently there are 333 listed
# I chose to filter by: chromosome_name, start, end
listFilters(mart)

# Read in tab-delimited file with three columns: chromosome number, start position and end position
positions <- read.table("positions.txt")

# Extract HGNC gene symbol
results <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters = c("chromosome_name", "start", "end"), values = list(positions[,1], positions[,2], positions[,3]), mart = mart)