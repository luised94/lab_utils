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
