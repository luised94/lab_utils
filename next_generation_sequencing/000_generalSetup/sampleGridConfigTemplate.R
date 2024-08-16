#Description: Configuration file that defines the categories of an experiment, creates the combinations of all the variables and then uses a filter function to grab the combinations.
#USAGE: This is the template for other experiments. Source the sampleGridConfig.R file in the script createSampleGrid.R, not the template file.
# This shows an example setup for BMC CHIP-seq experiment 240808Bel.
# @todo: Consider adding a comprehensive list or an alternative file with all of the variables that is generated programatically.

# @function: Grab the directory of the createSampleGrid script that runs this file. These two files are meant to be in the same directory. Not very flexible. The files accounts for running from interactive repl when testing or from Rscript via cli.
get_script_dir <- function() {
  if (!is.null(sys.frames()[[1]]$ofile)) {
    # 'source'd via R console
    script_dir <- dirname(normalizePath(sys.frames()[[1]]$ofile))
  } else {
    # Rscript
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
      script_dir <- dirname(normalizePath(sub(needle, "", cmdArgs[match])))
    } else {
      stop("Cannot determine script directory")
    }
  }
  return(script_dir)
}
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
        mcm_tag == "none" &
        cell_cycle == "M" &
        antibody == "Input" &
        ((strain_source == "oa" & auxin_treatment == "no") | strain_source == "lemr")
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

sample_grid <- filter_samples(expand.grid(categories))
print(dim(sample_grid))
print(table(sample_grid$antibody))
print(head(sample_grid))
print(get_script_dir())
