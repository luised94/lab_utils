library(httr)
library(jsonlite)

# Function to retrieve entries from UniProt API
getUniprotEntries <- function(query, format = "json") {
  base_url <- "https://rest.uniprot.org/uniprotkb/search"
  url <- paste0(base_url, "?query=", query, "&format=", format)
  response <- GET(url)
  
  if (status_code(response) == 200) {
    content(response, as = "text")
  } else {
    stop("Error retrieving UniProt entries.")
  }
}

# Function to parse JSON response from UniProt API
parseUniprotJSON <- function(json) {
  entries <- fromJSON(json)
  return(entries$results)
}

# Function to retrieve UniProt IDs of ORC genes in S. cerevisiae
getORCUniProtIDs <- function(taxon_id, format = "json") {
  query <- paste0("(organism_id:", taxon_id, ") AND (gene:ORC)")
  base_url <- "https://rest.uniprot.org/uniprotkb/search"
  url <- paste0(base_url, "?query=", query, "&format=", format)
  response <- GET(url)
  
  if (status_code(response) == 200) {
    content(response, as = "text")
  } else {
    stop("Error retrieving UniProt IDs for ORC genes.")
  }
}

# Function to retrieve protein sequences by UniProt ID
getProteinSequenceByUniProtID <- function(uniprot_id, format = "fasta") {
  base_url <- "https://www.uniprot.org/uniprot/"
  url <- paste0(base_url, uniprot_id, ".", format)
  response <- GET(url)
  
  if (status_code(response) == 200) {
    content(response, as = "text")
  } else {
    stop("Error retrieving protein sequence for UniProt ID.")
  }
}

# Example usage
# query <- "(reviewed:true) AND (organism_id:9606)"
query <- "accession:P54791"
format <- "json"

json_response <- getUniprotEntries(query, format)
entries <- parseUniprotJSON(json_response)

# Print UniProt IDs of retrieved entries
for (entry in entries) {
  print(entry$accession)
}

library(httr)
library(jsonlite)


# Example usage
taxon_id <- "559292"  # S. cerevisiae taxon ID
format <- "json"

# Step 1: Retrieve UniProt IDs of ORC genes in S. cerevisiae
json_response <- getORCUniProtIDs(taxon_id, format)
orc_gene_ids <- fromJSON(json_response)$results$accession

# Step 2: Retrieve protein sequences for each UniProt ID
for (orc_gene_id in orc_gene_ids) {
  fasta_sequence <- getProteinSequenceByUniProtID(orc_gene_id)
  cat(fasta_sequence, "\n")
}
