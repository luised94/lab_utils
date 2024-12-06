# Constants
CHROMOSOME_MAPPING <- list(
    numeric_to_roman = c(
        "1" = "I", "2" = "II", "3" = "III", "4" = "IV",
        "5" = "V", "6" = "VI", "7" = "VII", "8" = "VIII",
        "9" = "IX", "10" = "X", "11" = "XI", "12" = "XII",
        "13" = "XIII", "14" = "XIV", "15" = "XV", "16" = "XVI"
    ),
    roman_to_numeric = NULL  # Will be initialized in .onLoad
)

# Initialize reverse mapping
.onLoad <- function(libname, pkgname) {
    CHROMOSOME_MAPPING$roman_to_numeric <- setNames(
        names(CHROMOSOME_MAPPING$numeric_to_roman),
        CHROMOSOME_MAPPING$numeric_to_roman
    )
}

genome_detect_chr_style <- function(chr_names) {
    # Input validation
    if (!is.character(chr_names) || length(chr_names) == 0) {
        stop("Invalid chromosome names input")
    }
    
    # Define style patterns
    style_patterns <- list(
        UCSC = "^chr[0-9]+$",
        Roman = "^chr[IVX]+$",
        Numeric = "^[0-9]+$"
    )
    
    # Check each style
    for (style in names(style_patterns)) {
        if (all(grepl(style_patterns[[style]], chr_names))) {
            return(style)
        }
    }
    
    return("Unknown")
}

genome_convert_chr_names <- function(chr_names, target_style) {
    # Input validation
    if (!is.character(chr_names)) {
        stop("Chromosome names must be character vector")
    }
    if (!target_style %in% c("UCSC", "Roman", "Numeric")) {
        stop("Invalid target style. Must be 'UCSC', 'Roman', or 'Numeric'")
    }
    
    # Remove existing 'chr' prefix
    clean_names <- gsub("^chr", "", chr_names)
    
    # Convert based on target style
    converted_names <- switch(target_style,
        "UCSC" = paste0("chr", clean_names),
        "Roman" = sapply(clean_names, function(x) {
            if (x %in% names(CHROMOSOME_MAPPING$numeric_to_roman)) {
                paste0("chr", CHROMOSOME_MAPPING$numeric_to_roman[x])
            } else {
                paste0("chr", x)
            }
        }),
        "Numeric" = sapply(clean_names, function(x) {
            if (x %in% CHROMOSOME_MAPPING$numeric_to_roman) {
                CHROMOSOME_MAPPING$roman_to_numeric[x]
            } else {
                x
            }
        })
    )
    
    return(unname(converted_names))
}
