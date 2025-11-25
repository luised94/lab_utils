
# === PATHS ===
INPUT_DIR_path <- "~/data/protein_files"
OUTPUT_DIR_path <- "~/data/protein_files/alignments"

# === GENES TO PROCESS ===
GENE_NAMES_chr <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA1", "TOA2")

# === FILTERING CRITERIA ===
MIN_SEQUENCE_LENGTH_int <- 50  # minimum amino acids
KEEP_ONE_PER_SPECIES_lgl <- TRUE

# === REGIONS OF INTEREST (Manual specification) ===
# Positions refer to S. cerevisiae reference sequence (ungapped)
# TO BE UPDATED after alignment inspection

REGIONS_OF_INTEREST_lst <- list(
  ORC1 = c(100, 250, 500),  # Example - UPDATE AFTER INSPECTION
  ORC2 = c(150, 300),        # Example - UPDATE AFTER INSPECTION
  ORC3 = c(200, 400),        # Example - UPDATE AFTER INSPECTION
  ORC4 = c(180),             # Example - UPDATE AFTER INSPECTION
  ORC5 = c(220, 450),        # Example - UPDATE AFTER INSPECTION
  ORC6 = c(130, 280),        # Example - UPDATE AFTER INSPECTION
  TOA1 = c(100, 200),        # Example - UPDATE AFTER INSPECTION
  TOA2 = c(90, 180)          # Example - UPDATE AFTER INSPECTION
)

# Window around each position for zoomed plots
ZOOM_WINDOW_RESIDUES_int <- 10  # +/- residues around center position

# === ALIGNMENT CONFIGURATION ===
ALIGNMENT_METHOD_chr <- "ClustalOmega"
ALIGNMENT_ORDER_chr <- "input"
FORCE_REALIGNMENT_lgl <- FALSE

# === VISUALIZATION CONFIGURATION ===
PLOT_COLOR_SCHEME_chr <- "Chemistry_AA"
PLOT_FONT_chr <- "TimesNewRoman"
SHOW_SEQUENCE_LOGO_lgl <- TRUE
PLOT_WIDTH_inches <- 12
PLOT_HEIGHT_inches <- 8
