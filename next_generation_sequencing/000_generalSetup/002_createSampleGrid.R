rm(list = ls())
category_names <- c("rescue_allele", "mcm_tag", "auxin_treatment", "cell_cycle", "antibody")

strain_source <- c("lemr", "oa")
rescue_allele <- c("none", "wt")
mcm_tag <- c("none", "2", "4", "7")
mcm_tag <- c("none", "2", "7")
auxin_treatment <- c("no", "yes")
cell_cycle <- c("G1", "M")
antibody <- c("Input", "ProtG", "ALFA", "HM1108", "74", "CHA", "11HA")

variables_for_sample_grid <- ls()[!grepl("category_names", ls())]
sample_combinations <- expand.grid(lapply(variables_for_sample_grid, get))
names(sample_combinations) <- variables_for_sample_grid
sample_combinations <- sample_combinations[,c("strain_source", "rescue_allele", "mcm_tag", "auxin_treatment", "cell_cycle", "antibody")]

sample_combinations$full_name <- apply(sample_combinations, 1, paste, collapse = "_")
sample_combinations$short_name <- apply(sample_combinations[,!grepl("full_name", colnames(sample_combinations))], 1, function(row) paste0(substr(row, 1, 1), collapse = ""))

is_input <- sample_combinations[, "rescue_allele"] == "none" &
            sample_combinations[, "mcm_tag"] == "none" &
            sample_combinations[, "cell_cycle"] == "M" &
            sample_combinations[, "antibody"] == "Input" &
            ((sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "auxin_treatment"] == "no") | (sample_combinations[, "strain_source"] == "lemr"))

print(sum(is_input))
print(sample_combinations[is_input,])


#print(sample_combinations[is_input,])


is_alfa <- sample_combinations[, "rescue_allele"] == "none" &
            sample_combinations[, "mcm_tag"] == "none" &
            sample_combinations[, "cell_cycle"] == "M" &
            sample_combinations[, "antibody"] == "ALFA" &
            ((sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "auxin_treatment"] == "no") | (sample_combinations[, "strain_source"] == "lemr"))

print(sum(is_alfa))
#print(sample_combinations[is_alfa,])


is_1108 <- sample_combinations[, "rescue_allele"] == "none" &
            sample_combinations[, "mcm_tag"] == "none" &
            sample_combinations[, "cell_cycle"] == "M" &
            sample_combinations[, "antibody"] == "HM1108" &
            ((sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "auxin_treatment"] == "no") | (sample_combinations[, "strain_source"] == "lemr"))


print(sum(is_1108))
#print(sample_combinations[is_1108,])


is_174 <- sample_combinations[, "antibody"] == "74" &
          sample_combinations[, "auxin_treatment"] == "no" &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "rescue_allele"] == "none") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "rescue_allele"] == "wt") &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "mcm_tag"] == "7") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "mcm_tag"] == "2")
print(sum(is_174))
print(sample_combinations[is_174,])

is_cha <- sample_combinations[, "antibody"] == "CHA" &
          sample_combinations[, "auxin_treatment"] == "no" &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "rescue_allele"] == "none") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "rescue_allele"] == "wt") &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "mcm_tag"] == "7") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "mcm_tag"] == "2")
print(sum(is_cha))
print(sample_combinations[is_cha,])


is_11HA <- sample_combinations[, "antibody"] == "11HA" &
          sample_combinations[, "auxin_treatment"] == "no" &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "rescue_allele"] == "none") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "rescue_allele"] == "wt") &
          !(sample_combinations[, "strain_source"] == "lemr" & sample_combinations[, "mcm_tag"] == "7") &
          !(sample_combinations[, "strain_source"] == "oa" & sample_combinations[, "mcm_tag"] == "2")
print(sum(is_11HA))
print(sample_combinations[is_11HA,])

subset_vectors <- ls()[grepl("is_", ls())]
sample_table <- do.call(rbind, lapply(subset_vectors, function(subset_vector){
    sample_combinations[get(subset_vector), ]
}))
print(sample_table)
print(dim(sample_table))

bmc_table <- data.frame(SampleName = sample_table$full_name,
    Vol..uL = rep(10, nrow(sample_table)),
    Conc = rep(NA, nrow(sample_table)),
    Type = rep("ChIP", nrow(sample_table)),
    Genome = rep("Saccharomyces cerevisiae", nrow(sample_table)),
    Notes = rep("none", nrow(sample_table)),
    Pool = rep("A", nrow(sample_table))
)
write.table(sample_table, file = "sample_info.tsv", sep = "\t", row.names = FALSE)
write.table(bmc_table, file = "bmc_info.tsv", sep = "\t", row.names = FALSE)
