#!/usr/bin/env Rscript
## Quick test of stringi parsing of experiment id for argument handling.
#library(stringi)
#experiment_id_string <- "241010Bel,241007Bel,241122Bel"
#print(grepl(",", "241010Bel,241007Bel,241122Bel"))
#split_string <- stri_split_fixed(experiment_id_string, ",")
#
#print(unlist(stri_split_fixed("241010Bel", ",")))
#print(unlist(split_string))
#unlist_strings <- unlist(split_string)
#conform_to_pattern <- sapply(unlist_strings, function(x){
#    grepl("^\\d{6}Bel$", x, perl = TRUE)
#})
#print(all(conform_to_pattern))
#file.path(Sys.getenv("HOME"), "data", unlist_strings)
#
#unlist_strings <- unlist(stri_split_fixed("241010Bel", ","))
#conform_to_pattern <- sapply(unlist_strings, function(x){
#    grepl("^\\d{6}Bel$", x, perl = TRUE)
#})
#
#if (!all(conform_to_pattern) && args$experiment_id != "template") {
#    stop(sprintf(
#        "Invalid experiment-id format.\nExpected: YYMMDD'Bel' or 'template'\nReceived: %s",
#        args$experiment_id
#    ))
#} else {
#    print("Completed conditions.")
#    print(file.path(Sys.getenv("HOME"), "data", unlist_strings))
#
#}
#
#unlist_strings <- unlist(stri_split_fixed(experiment_id_string, ","))
#conform_to_pattern <- sapply(unlist_strings, function(x){
#    grepl("^\\d{6}Bel$", x, perl = TRUE)
#})
#
#if (!all(conform_to_pattern) && args$experiment_id != "template") {
#    stop(sprintf(
#        "Invalid experiment-id format.\nExpected: YYMMDD'Bel' or 'template'\nReceived: %s",
#        args$experiment_id
#    ))
#} else {
#    print("Completed conditions.")
#    print(file.path(Sys.getenv("HOME"), "data", unlist_strings))
#
#}
## Quick test of stringi parsing of experiment id for argument handling.
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#                        DEBUG CONTROL SYSTEM DEMO                        #
#                 Execute with --debug flag to test behavior              #
# Usage:                                                                  #
#   Rscript demo.R --debug=all          # Enable all debug blocks         #
#   Rscript demo.R --debug=1,3          # Enable blocks 1 and 3 by number #
#   Rscript demo.R --debug=metadata     # Enable by tag                   #
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

### [1] Debug Control System Setup #########################################

# Central debug configuration environment
DEBUG_CONFIG <- new.env(parent = emptyenv())

# Active blocks (character vector of block IDs)
DEBUG_CONFIG$active_blocks <- character(0)

# Execution counters (list: block_id -> count)
DEBUG_CONFIG$execution_counts <- list()

# Registry of all blocks (list: block_id -> human-readable tag)
DEBUG_CONFIG$block_registry <- list()

### [2] Core System Functions ##############################################

register_debug_block <- function(block_id, tag = NULL) {
  # Automatically register new blocks when first encountered
  if (!block_id %in% names(DEBUG_CONFIG$block_registry)) {
    DEBUG_CONFIG$block_registry[[block_id]] <- if(!is.null(tag)) tag else block_id
    DEBUG_CONFIG$execution_counts[[block_id]] <- 0
  }
}

execute_debug <- function(block_id, expr, tag = NULL) {
  # Register block (if new)
  register_debug_block(block_id, tag)
  
  # Execute only if block is active
  if (block_id %in% DEBUG_CONFIG$active_blocks) {
    force(expr)
    DEBUG_CONFIG$execution_counts[[block_id]] <- 
      DEBUG_CONFIG$execution_counts[[block_id]] + 1
  }
}

### [3] Command Line Interface Setup #######################################

library(optparse)

option_list <- list(
  make_option(c("--debug"), type = "character", default = "none",
              help = "Enable debug blocks: 'all', 'none', '1,3', or tags")
)

args <- parse_args(OptionParser(option_list = option_list))

# Process --debug argument
if (args$debug == "all") {
  DEBUG_CONFIG$active_blocks <- names(DEBUG_CONFIG$block_registry)
} else if (args$debug != "none") {
  inputs <- unlist(strsplit(args$debug, ","))
  resolved <- sapply(inputs, function(x) {
    # Numeric input (e.g., "1" -> "DB1")
    if (grepl("^\\d+$", x)) return(paste0("DB", x))
    
    # Tag-based lookup
    matches <- names(DEBUG_CONFIG$block_registry)[
      DEBUG_CONFIG$block_registry %in% x
    ]
    if (length(matches) > 0) return(matches)
    
    # Fallback to direct ID
    x
  })
  DEBUG_CONFIG$active_blocks <- unique(unlist(resolved))
}

### [4] Sample Debug Blocks ################################################

# Simple wrapper for cleaner syntax when not capturing output
debug_block <- function(block_id, tag = NULL, ...) {
  execute_debug(block_id, tag = tag, {
    cat(paste0("\n=== DEBUG BLOCK ", block_id, " ===\n"))
    # Properly handle multiple arguments
    args <- list(...)
    if(length(args) > 0) {
      print(args)
    } else {
      cat("(No debug output)\n")
    }
    cat(paste0(rep("~", 40), collapse = ""), "\n")
  })
}
# Demo block 1 - Basic debug output
debug_block("DB1", tag = "metadata",
           "Metadata summary:",
           list(samples = 100, columns = 20))

# Demo block 2 - Conditional processing
debug_block("DB2", tag = "processing",
           "Processing results:",
           list(stats = rnorm(3), status = "complete"))

# Demo block 3 - Performance metrics (simulated)
debug_block("DB3", tag = "performance", {
  Sys.sleep(0.2)  # Simulate work
  list(time_elapsed = Sys.time(), memory = paste0(object.size(x = .GlobalEnv), " bytes"))
})

### [5] System Introspection ##############################################

cat("\n\n=== DEBUG SYSTEM SUMMARY ===")
cat("\nActive blocks:", 
    paste(DEBUG_CONFIG$active_blocks, collapse = ", "))
cat("\nExecution counts:")
print(DEBUG_CONFIG$execution_counts)
cat("\nBlock registry:")
print(DEBUG_CONFIG$block_registry)
