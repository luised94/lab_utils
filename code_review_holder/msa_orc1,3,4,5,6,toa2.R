library(msa)
library(ggmsa)
library(ggplot2)
library(stringr)
# Function to write protein sequences to a FASTA file
writeFasta <- function(dataframe, output_file) {
  fasta_strings <- character()
  
  for (i in 1:nrow(dataframe)) {
    fasta_strings <- c(
      fasta_strings,
      paste0(">", dataframe$Entry.Name[i], "\n", dataframe$Sequence[i])
    )
  }
  
  writeLines(fasta_strings, output_file)
}
#write the fasta files and read them back to generate an msa. Get the files from uniprot.
writeFastaAndGenerateMSA <- function(dataframe, grep_pattern, output_prefix) {
  # Filter the dataframe using grep
  filtered_df <- dataframe[grep(grep_pattern, dataframe$Entry.Name),]
  
  # Create the output file name
  output_file <- paste0(output_prefix, ".fasta")
  
  # Write protein sequences to FASTA file
  writeFasta(filtered_df, output_file)
  
  # Read protein sequences from the FASTA file
  orc_seqs <- readAAStringSet(filepath = output_file)
  
  # Generate MSA
  msa_res <- msaClustalOmega(orc_seqs)
  
  # Generate MSA LaTeX file
  output_tex_file <- paste0(output_prefix, "_msa.tex")
  msaPrettyPrint(
    x = msa_res,
    output = "tex",
    alFile = paste0(output_prefix, "_msa.fasta"),
    file = output_tex_file,
    showNames = "left",
    askForOverwrite = FALSE,
    shadingMode = "functional",
    shadingModeArg = "structure",
    showNumbering = "right",
    verbose = FALSE
  )
  
  # Compile LaTeX to PDF using TinyTeX
  tinytex::pdflatex(output_tex_file)
}
#Function to generate and write MSA. Use if you have fasta file already.
generateAndWriteMSA <- function(file_path) {
  # Read protein sequences from the FASTA file
  fasta_file <- readAAStringSet(filepath = file_path)
  output_prefix <- strsplit(basename(file_path), "_")[[1]][1]
  output_fasta <- paste0(output_prefix, "_msa.fasta")
  
  if (file.exists(output_fasta)){
    cat(paste("File", output_fasta, "already exists.", sep = " "))
    overwrite_align <- readline(prompt = " Overwrite? (y or n): ")
    if(overwrite_align){
      # Generate MSA
      msa_res <- msaClustalOmega(fasta_file)
      # Write MSA to FASTA file
      writeXStringSet(AAStringSet(AAMultipleAlignment(msa_res)),  output_fasta, format = "fasta") 
      cat(paste("File", output_fasta, "overwritten.", sep = " "))
    }
  } else {
    # Generate MSA
    msa_res <- msaClustalOmega(fasta_file)
    # Write MSA to FASTA file
    writeXStringSet(AAStringSet(AAMultipleAlignment(msa_res)),  output_fasta, format = "fasta") 
    cat(paste("File", output_fasta, "created.", sep = " "))
  }
}
#Function to generate and write MSA. Use if you have a AAStringSet object from biostrings package.
alignAndWriteMSA <- function(AAStringSetfasta, output_prefix) {
  output_fasta <- paste0(output_prefix, "_msa.fasta")
  
  if (file.exists(output_fasta)){
    cat(paste("File", output_fasta, "already exists.\n", sep = " "))
    overwrite_align <- readline(prompt = " Overwrite? (T or F): ")
    if(overwrite_align){
      # Generate MSA
      msa_res <- msaClustalOmega(AAStringSetfasta, order = "input")
      # Write MSA to FASTA file
      writeXStringSet(AAStringSet(AAMultipleAlignment(msa_res)),  output_fasta, format = "fasta") 
      cat(paste("File", output_fasta, "overwritten.\n", sep = " "))
    }
  } else {
    # Generate MSA
    msa_res <- msaClustalOmega(AAStringSetfasta)
    # Write MSA to FASTA file
    writeXStringSet(AAStringSet(AAMultipleAlignment(msa_res)),  output_fasta, format = "fasta") 
    cat(paste("File", output_fasta, "created.\n", sep = " "))
  }
}
#generate and plot msa from fasta file.
generateAndPlotMSA <- function(file_path, center_position, around = 10) {
  # Read protein sequences from the FASTA file
  seqs_to_align <- readAAStringSet(filepath = file_path)
  
  # Generate MSA
  msa_res <- msaClustalOmega(seqs_to_align)
 
  # Plot the MSA
  ggmsa(AAMultipleAlignment(msa_res), center_position-around, center_position+around, 
        color = "Chemistry_AA", font = "TimesNewRoman", 
        char_width = 0.5, seq_name = TRUE) + 
    geom_seqlogo(color = "Chemistry_AA")
}
#Random function to concatenate the elements of list. 
pasteAndCopy <- function(string_list, separator = "") {
  combined_string <- paste(string_list, collapse = separator)
  writeClipboard(combined_string)
}
#Plot msa given center position and residues around it. Use if you have an msa already. 
plotMSA <- function(msa_to_plot, center_position, around = 10){
  ggmsa(AAMultipleAlignment(readAAMultipleAlignment(msa_to_plot)), center_position-around, center_position+around, 
        color = "Chemistry_AA", font = "TimesNewRoman", 
        char_width = 0.5, seq_name = TRUE) + geom_seqlogo(color = "Chemistry_AA") + 
    xlab("Alignment Column Number") + ylab("Protein Sequence")
}
#Rename and reorder fasta 
reorder_fasta <- function(fasta_path) {
  library(Biostrings)
  library(stringr)
  
  # Read the FASTA file
  fasta_to_align <- readAAStringSet(fasta_path)
  
  # Use regular expression to extract the desired substring
  gene_name <- str_extract(names(fasta_to_align), "\\|([^|]+)\\s")
  
  # Remove the leading and trailing "|" characters
  gene_name <- gsub("\\|", "", gene_name)
  
  simple_names <- c()
  
  for (i in 1:length(gene_name)) {
    original_name <- strsplit(gene_name, " ")[[i]][1]
    parts <- strsplit(original_name, "_")[[1]]
    reordered_string <- paste(parts[2], parts[1], sep = "_")
    simple_names <- c(simple_names, reordered_string)
  }
  
  names(fasta_to_align) <- simple_names
  
  fasta_to_align <- fasta_to_align[order(fasta_to_align@ranges@NAMES, decreasing = TRUE)]
  
  return(fasta_to_align)
}

####Uncomment this to generate fastas####
# fasta_to_align <- list.files("C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Inventory\\Protein files\\", pattern = "*.fasta$", full.names = TRUE)
# lapply(fasta_to_align, function(x){
#   output_prefix <- strsplit(basename(x), "_")[[1]][1]
#   alignAndWriteMSA(reorder_fasta(fasta_path = x), output_prefix = output_prefix)
# })

####Generate MSAs####
# Define the msa_to_plot values. Can subset or look elsewhere for fasta files.
msas_to_plot <- list.files(pattern = "_msa.fasta")
# Define the center positions, around values, and the number of times each msa is used
center_positions <- c(806, 896, 946, 822, 767, 586, 704, 754, 712, 253, 364, 16)
around_values <- rep(10, length(center_positions))
msa_repeats <- c(4, 1, 4, 1, 1, 1)  # The number of times each msa_to_plot is repeated


center_position_index <- 1
# Generate plots using lapply
lapply(1:length(msas_to_plot), function(i) {
  msa_to_plot <- msas_to_plot[[i]]
  repeat_count <- msa_repeats[i]
  output_prefix <- strsplit(msa_to_plot, "_")[[1]][1]
  
  for (j in 1:repeat_count) {
    # print(paste("i:", i, "j:", j, "center_position_index:", center_position_index))
    paste0(output_prefix, "_plot_", i, "_", j, ".svg")
    png_file <- paste0("./output/",output_prefix, "_plot_", i, "_", j, ".png")
    svg_file <- paste0("./output/",output_prefix, "_plot_", i, "_", j, ".svg")

    plotMSA(msa_to_plot = msa_to_plot,
            center_position = center_positions[center_position_index],
            around = around_values[center_position_index])
    ggsave(filename = svg_file)
    ggsave(filename = png_file)
    # readline("Press Enter to continue...")
    cat("Plot saved to", svg_file, "\n")
    cat("Plot saved to", png_file, "\n")

    center_position_index <<- center_position_index + 1  # Move to the next position
    }
})



# tsv_files <- list.files("C:\\Users\\Luis\\Dropbox (MIT)\\Lab\\Inventory\\Protein files\\", pattern = "*.tsv$", full.names = TRUE)
# df <- data.frame()
# for (i in 1:length(tsv_files)){
#   data <- read.delim(tsv_files[i])
#   df <- rbind(df, data)
#   df <- unique(df)
# }
# df <- unique(df)
# 
# data_list <- lapply(tsv_files, read.delim)
# df <- do.call(rbind, data_list)
# df <- unique(df)

# msa_to_generate <- c("ORC1", "ORC4")
# 
# for (i in 1:length(msa_to_generate)){
#   writeFastaAndGenerateMSA(df, msa_to_generate[i], tolower(msa_to_generate[i]))
# }
# 
# #Unfortunaly had to use snapgene to figure out what number to use to display AAA+ motifs .
# generateAndPlotMSA("./orc1.fasta", 809)
# generateAndPlotMSA("./orc1.fasta", 896)
# generateAndPlotMSA("./orc1.fasta", 946)
# generateAndPlotMSA("./orc4.fasta", 586)
# generateAndPlotMSA("./orc4.fasta", 704)
# generateAndPlotMSA("./orc4.fasta", 754)

# yeast_prots <- df[setdiff(grep("cerevisiae", df$Organism),grep("ORC6|MCM1|ORC2|ORC5|ORC3", df$Entry.Name)), ]

# rel_json <- drawProteins::get_features(paste(yeast_prots$Entry, collapse = " "))
# df_features<- drawProteins::feature_to_dataframe(rel_json)

# drawProteins::draw_regions(p, df_features)# read.delim("C:\\Users\\Luis\\Downloads\\entry-matching-P09119.tsv")
# df_interpro <- read.delim("C:\\Users\\Luis\\Downloads\\protein-sequences_cerevisiae_interpro.tsv")
# str(df_interpro)
# df_interpro[df_interpro$Accession  %in% yeast_prots$Entry,]
# /api/protein/reviewed/entry/InterPro/taxonomy/uniprot/559292/?page_size=20
# 
# yeast_prots <- yeast_prots[order(yeast_prots$Entry.Name), ]
#
# # Output file path
# output_file <- "yeast_loading_AAA+.fasta"
#
# # Write protein sequences to FASTA file
# writeFasta(yeast_prots, output_file)
# mySeqs <- readAAStringSet(filepath = "./yeast_loading_AAA+.fasta")
# res <- msaClustalW(mySeqs)
# msaPrettyPrint(x=res, output="tex", file="HemoglobinExample.tex",
#                showNames = "left",
#                showLogo = "none",
#                askForOverwrite = FALSE,
#                verbose = FALSE)
# # msaPrettyPrint(x=res, output="tex", file="HemoglobinExample.tex",
# #                paperWidth=15.3, paperHeight=4.5,
# #                shadingMode="functional",
# #                shadingModeArg="structure",
# #                shadingColors="greens",
# #                logoColors="rasmol",
# #                showLogoScale="right",
# #                showLegend=FALSE,
# #                askForOverwrite=FALSE)
# # texshade_path <- system.file("tex", "texshade.sty", package="msa")
# tinytex::pdflatex("HemoglobinExample.tex")


#https://rest.uniprot.org/uniprotkb/accessions?accessions=O16810%2CO43929%2CO74270%2CO88708%2CO93479%2CP09119%2CP11746%2CP24279%2CP29469%2CP29496%2CP30665%2CP32833%2CP38132%2CP38826%2CP50874%2CP53091%2CP54784%2CP54789%2CP54790%2CP54791%2CQ13415%2CQ6EWX1%2CQ6P9Z8%2CQ80Z32%2CQ9SU24%2CQ9Y794%2CQ9Z1N2&compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Csequence&format=tsv
## white space around the resulting plot can be removed using the tool pdfcrop
## system("pdfcrop HemoglobinExample.pdf")

# df[grep("ORC1", df$Entry.Name), ][,"Entry.Name"]
# df[grep("ORC4", df$Entry.Name), ][,"Entry.Name"]
# df[grepl("cerevisiae", df$Organism)&!grepl("ORC6|MCM1", df$Entry.Name), ][,"Entry.Name"]
# 
# grep("cerevisiae", df$Organism)
# grep("ORC6|MCM1", df$Entry.Name)

# for (i in 1:nrow(yeast_prots)) {
#   fasta_strings <- c(
#     fasta_strings,
#     paste0(">", yeast_prots$Entry.Name[i], "\n", yeast_prots$Sequence[i])
#   )
# }
# 
# writeLines(fasta_strings, output_file)

