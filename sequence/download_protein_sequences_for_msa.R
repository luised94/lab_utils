# CONSTANTS BLOCK
# Organism taxonomy IDs
ORGANISM_TAXIDS_int <- c(
  559292,   # S. cerevisiae
  27292,    # S. pastorianus
  1080349,  # S. arboricola
  1080088,  # S. eubayanus
  230603,   # S. uvarum
  27291,    # S. paradoxus
  4954,     # Z. rouxii
  2763761,  # Z. mellis
  37769,    # T. globosa
  5478,     # C. glabrata
  284812,   # S. pombe
  7227,     # D. melanogaster
  6239,     # C. elegans
  10090,    # M. musculus
  9606      # H. sapiens
)

ORGANISM_NAMES_chr <- c(
  "Scer", "Spas", "Sarb", "Seub", "Suva", "Spar",
  "Zrou", "Zmel", "Tglo", "Cgla", "Spom",
  "Dmel", "Cele", "Mmus", "Hsap"
)

# Genes to query
GENE_NAMES_chr <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA1", "TOA2")

# UniProt API
BASE_URL_uniprot_chr <- "https://rest.uniprot.org"
ENDPOINT_search_chr <- "/uniprotkb/search"

# Output
OUTPUT_DIR_path <- "~/data/protein_files"
OUTPUT_PREFIX_chr <- "model_organisms"

# Cache control
FORCE_DOWNLOAD_lgl <- FALSE

# API behavior
RETRY_MAX_int <- 3
RETRY_DELAY_sec <- 2
USER_AGENT_chr <- "ORC_MSA_Analysis/1.0 (R_script)"
