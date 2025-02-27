# Function to retrieve residue numbers of prolines in a protein sequence
getAminoAcidResidueNumbers <- function(sequence, aminoacid) {
  residue <- gregexpr(aminoacid, sequence)
  return(unlist(residue))
}

# Example usage
protein_sequence <- entries$sequence$value
aminoacid <- "P"
prolines_residue_numbers <- getAminoAcidResidueNumbers(protein_sequence, aminoacid)

# Print residue numbers of prolines
cat("Residue numbers of prolines:", paste(prolines_residue_numbers, collapse = ", "))


residues_for_ddmut <- function(protein_chain, aminoacid, mutation) {
  chain_aminoacid <- paste(protein_chain, aminoacid, sep = " ")
  cat(file = "ddmut_sample.txt", paste(chain_aminoacid, prolines_residue_numbers, mutation, sep ="", "\n"))
}

residues_for_ddmut("D", "P", "L")
