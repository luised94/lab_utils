#Description: Configuration file that defines the categories of an experiment, creates the combinations of all the variables and then uses a filter function to grab the combinations.
#USAGE: This is the template for other experiments. Source the sampleGridConfig.R file in the script createSampleGrid.R, not the template file.
# This shows an example setup for BMC CHIP-seq experiment 240808Bel.
# @todo: Consider adding a comprehensive list or an alternative file with all of the variables that is generated programatically.
# To update efficiently use di<character> and yi<character> on the current_experiment, categories and filter_samples variables and function. Update order statement appropriately.
cat("Starting sample grid config.\n")
current_experiment <- "240808Bel"
cat(sprintf("Categories and filter_samples will be configured for %s", current_experiment), "\n")

# Create a list with the different categories and variables in the experiment.
categories <- list(
    strain_source = c("lemr", "oa"),
    rescue_allele = c("none", "wt"),
    mcm_tag = c("none", "2", "7"),
    auxin_treatment = c("no", "yes"),
    cell_cycle = c("G1", "M"),
    antibody = c("Input", "ProtG", "ALFA", "HM1108", "74", "CHA", "11HA")
)

#Define the indexes for filtering all of the combinations of the variables.
# Pick one of the variables and define how it is related to the other variables using conditional expressions. For example, for all of the antibodies, define the other conditions it is used with.
filter_samples <- function(combinations){
    #is_not <- with(combinations,
    #)
    is_input <- with(combinations,
        rescue_allele == "none" &
        cell_cycle == "M" &
        antibody == "Input" &
        ((strain_source == "oa" & auxin_treatment == "no") | (strain_source == "lemr" & auxin_treatment == "no")) &
        !( strain_source == "lemr" & mcm_tag == "7" ) &
        !( strain_source == "lemr" & mcm_tag == "2" ) &
        !( strain_source == "oa" & mcm_tag == "2" ) &
        !( strain_source == "oa" & rescue_allele == "wt" )
    )

    is_protg <- with(combinations,
            rescue_allele == "wt" &
            mcm_tag == "none" &
            cell_cycle == "M" &
            antibody == "ProtG" &
            strain_source == "oa" &
            auxin_treatment == "no"
    )

    is_alfa <- with(combinations,
        rescue_allele == "none" &
        mcm_tag == "none" &
        cell_cycle == "M" &
        antibody == "ALFA" &
        (( strain_source == "oa" & auxin_treatment == "no") | ( strain_source == "lemr"))
    )

    is_1108 <-  with(combinations,
            rescue_allele == "none" &
            mcm_tag == "none" &
            cell_cycle == "M" &
            antibody == "HM1108" &
           (( strain_source == "oa" &  auxin_treatment == "no") | strain_source == "lemr")
    )

    is_174 <- with(combinations,
        antibody == "74" &
        auxin_treatment == "no" &
        !( strain_source == "lemr" &  rescue_allele == "none") &
        !( strain_source == "oa" &  rescue_allele == "wt") &
        !( strain_source == "lemr" &  mcm_tag == "7") &
        !( strain_source == "oa" &  mcm_tag == "2")
    )

    is_cha <- with(combinations,
        antibody == "CHA" &
        auxin_treatment == "no" &
        !(strain_source == "lemr" & rescue_allele == "none") &
        !(strain_source == "oa" & rescue_allele == "wt") &
        !(strain_source == "lemr" & mcm_tag == "7") &
        !(strain_source == "oa" & mcm_tag == "2")
     )

    is_11HA <- with(combinations,
        antibody == "11HA" &
        auxin_treatment == "no" &
        !(strain_source == "lemr" & rescue_allele == "none") &
        !(strain_source == "oa" & rescue_allele == "wt") &
        !(strain_source == "lemr" & mcm_tag == "7") &
        !(strain_source == "oa" & mcm_tag == "2") &
        !(strain_source == "lemr" & mcm_tag == "none" & rescue_allele == "wt" & cell_cycle == "M")
    )
    return(combinations[is_input | is_protg | is_alfa | is_1108 | is_174 | is_cha | is_11HA , ])
}

sample_table <- filter_samples(expand.grid(categories))

sample_table <- sample_table[with(sample_table, order(antibody, strain_source)), ]
cat("Dimensions of sample_table: \n")
dim(sample_table)
cat("Breakdown by antibody")
print(table(sample_table$antibody))
cat("First elements of sample_table:\n")
print(head(sample_table))

sample_table$full_name <- apply(sample_table, 1, paste, collapse = "_")
sample_table$short_name <- apply(sample_table[,!grepl("full_name", colnames(sample_table))], 1, function(row) paste0(substr(row, 1, 1), collapse = ""))


bmc_table <- data.frame(SampleName = sample_table$full_name,
    Vol..uL = 10,
    Conc = NA,
    Type = "ChIP",
    Genome = "Saccharomyces cerevisiae",
    Notes = "none",
    Pool = "A"
)
#print(head(bmc_table))
#print(ls())
print("sampleGridConfig section complete")
