INVALID_COMBINATIONS <- list(
    # Group 1: orc1-161 restrictions
    orc1_161_restrictions = quote(
        orc_phenotype == "orc1-161" & 
        (temperature == "23" |                  # no orc1-161 at 23øC
         antibody == "ORC" |                    # no orc1-161 with ORC
         cell_cycle %in% c("async", "alpha"))      # no orc1-161 in async or alpha
    ),
    # Group 2: ORC antibody restrictions
    orc_restrictions = quote(
        antibody == "ORC" &
        (temperature == "23" |                  # no ORC at 23øC
         cell_cycle %in% c("alpha", "async"))      # no ORC in alpha or async
    ),
    # Group 2: ORC antibody restrictions
    orc_restrictions = quote(
        antibody == "ORC" & 
        (temperature == "23" |                  # no ORC at 23øC
         cell_cycle %in% c("alpha", "async"))      # no ORC in alpha or async
    ),
    # Group 3: Nucleosome and temperature restrictions
    nucleosome_temp_restrictions = quote(
        (antibody == "Nucleosomes" & temperature == "23" & cell_cycle %in% c("alpha", "nocodazole")) | # no Nucleosomes at 23øC in alpha/nocodazole
        (temperature == "37" & cell_cycle == "async")    # no async at 37øC
    )
),
