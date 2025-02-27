tsv_files <- list.files("C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Projects\\automate-the-boring-stuff\\reference_data\\Interpro-tsv\\", pattern = "*.tsv$", full.names = TRUE)
data_list <- lapply(tsv_files, read.delim)
df <- do.call(rbind, data_list)
df <- unique(df)
df_names <- names(df)

yeast_prots <- df[setdiff(grep("cerevisiae", df$Organism),grep("MCM1", df$Entry.Name)), ]
yeast_prots$Entry

writeClipboard(yeast_prots$Entry)
pasteAndCopy(yeast_prots$Entry, separator = " OR ")

lapply(yeast_prots$Entry, function(entry) {
  writeClipboard(entry)
  cat("Copied:", entry, "\n")
  readline("Press Enter to proceed...")
})
writeFasta(yeast_prots, "loading_prots.fasta")
