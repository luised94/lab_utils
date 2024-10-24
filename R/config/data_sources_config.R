#' Add to existing configurations
CONFIG <- list(
    SOURCES = list(
        HAWKINS = list(
            URL = "https://ars.els-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx",
            FILE = "hawkins-origins-timing.xlsx",
            TYPE = "timing"
        ),
        
        EATON_PEAKS = list(
            URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM424nnn/GSM424494/suppl/GSM424494_wt_G2_orc_chip_combined.bed.gz",
            FILE = "eaton_peaks.bed.gz",
            TYPE = "peaks"
        ),
        
        SGD = list(
            FEATURES = list(
                URL = "https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab",
                FILE = "SGD_features.tab",
                TYPE = "features"
            ),
            GFF = list(
                URL = "https://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz",
                FILE = "saccharomyces_cerevisiae.gff.gz",
                TYPE = "annotation"
            )
        ),
        
        EATON_ACS = list(
            URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM424nnn/GSM424494/suppl/GSM424494_acs_locations.bed.gz",
            FILE = "eaton_acs.bed.gz",
            TYPE = "acs"
        )
    ),
    
    SGD_COLUMNS = c(
        "Primary_SGDID", "Feature_type", "Feature_qualifier",
        "Feature_name", "Standard_gene_name", "Alias",
        "Parent_feature_name", "Secondary_SGDID", "Chromosome",
        "Start_coordinate", "Stop_coordinate", "Strand",
        "Genetic_position", "Coordinate_version",
        "Sequence_version", "Description"
    )
)
