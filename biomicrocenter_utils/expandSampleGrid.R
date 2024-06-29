#Variables
rescue_allele <- c("none", "WT", "4R")
suppressor_allele <- c("none", "1EK", "4PS", "T2")
cell_cycle <- c("G1", "G2")
treatment <- c("No", "Yes")
antibody <- c("ProtG", "INP", "ALFA", "V5", "174", "185")

#Create dataframe and rename columns.
all_combinations <- expand.grid(rescue_allele, suppressor_allele, cell_cycle, treatment, antibody, stringsAsFactors = TRUE)
colnames(all_combinations) <- c("rescue_allele", "suppressor_allele", "cell_cycle", "treatment", "antibody")

#Create the logical vectors to index the dataframe
is_none_sample <- all_combinations[, "rescue_allele"] == "none" & all_combinations[, "suppressor_allele"] == "none" & !(all_combinations[, "antibody"] %in% c("ProtG", "INP"))
is_wildtype_sample <- all_combinations[, "rescue_allele"] == "WT" & all_combinations[, "suppressor_allele"] == "none" & all_combinations[, "treatment"] == "Yes" & !(all_combinations[, "antibody"] %in% c("ProtG", "INP"))
is_fourr_sample <- all_combinations[, "rescue_allele"] == "4R" & all_combinations[, "treatment"] == "Yes" & !(all_combinations[, "antibody"] %in% c("ProtG", "185", "INP"))
is_fourrAnd185_sample <- all_combinations[, "rescue_allele"] == "4R" & all_combinations[, "suppressor_allele"] %in% c("none", "4PS") & all_combinations[, "treatment"] == "Yes" & all_combinations[, "antibody"] == "185"
is_input_sample <- ((all_combinations[, "rescue_allele"] == c("none", "WT") & all_combinations[, "suppressor_allele"] == "none")  | (all_combinations[, "rescue_allele"] == "4R")) & all_combinations[, "antibody"] == "INP" & all_combinations[, "treatment"] == "Yes" & all_combinations[, "cell_cycle"] == "G1"
is_protg_sample <- all_combinations[, "rescue_allele"] == "WT" & all_combinations[, "suppressor_allele"] == "none" & all_combinations[, "treatment"] == "Yes" & all_combinations[, "cell_cycle"] == "G1" & all_combinations[, "antibody"] == "ProtG"
print(all_combinations[is_fourrAnd185_sample,])
print(sum(is_none_sample))
print(sum(is_wildtype_sample))
print(sum(is_fourr_sample))
print(sum(is_fourrAnd185_sample))
print(sum(is_input_sample))
print(sum(is_protg_sample))
sample_set <- rbind(
    all_combinations[is_none_sample,],
    all_combinations[is_wildtype_sample,],
    all_combinations[is_fourr_sample,],
    all_combinations[is_fourrAnd185_sample,],
    all_combinations[is_input_sample,],
    all_combinations[is_protg_sample,]
)
print(sample_set[with(sample_set, order(antibody, rescue_allele, suppressor_allele, cell_cycle, treatment)),])
sample_set <- sample_set[with(sample_set, order(antibody, rescue_allele, suppressor_allele, cell_cycle, treatment)),]
print(nrow(sample_set))
sample_set[ , "full_name"] <- apply(sample_set, 1, function(row) paste(row, collapse = "_"))
print(sample_set[, "full_name"])
##write.table(sample_set, "samples.tsv", sep="\t", row.names=FALSE)
