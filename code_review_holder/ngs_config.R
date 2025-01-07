#' Merge with existing configurations
CONFIG <- list(
    PATHS = list(
        BASE_DIR = Sys.getenv("HOME"),
        SUBDIRS = list(
            ALIGNMENT = "alignment",
            BIGWIG = "bigwig",
            DOCUMENTATION = "documentation"
        )
    ),
    
    PATTERNS = list(
        SAMPLE_TABLE = "sample_table",
        CONTROL_FACTOR_PREFIX = "X__cf_",
        REFERENCE_GENOME = "S288C",
        BAM_SUFFIX = ".bam$"
    ),
    
    CONTROL = list(
        MAX_CONTROLS = 1,
        DEFAULT_INDEX = 1,
        ANTIBODY_COLUMN = "antibody",
        INPUT_VALUE = "Input"
    )
)

#' Add to existing configurations
CONFIG <- list(
    PATHS = list(
        SUBDIRS = list(
            ALIGNMENT = "alignment",
            BIGWIG = "bigwig"
        )
    ),
    
    FILES = list(
        PATTERNS = list(
            BAM = ".bam$",
            REFERENCE = "S288C"
        )
    ),
    
    VALIDATION = list(
        REQUIRED_COLUMNS = c(
            "sample_ID",
            "short_name"
        )
    )
)

#' Genome Configuration Parameters
CONFIG <- list(
    PATHS = list(
        BASE_DIR = file.path(Sys.getenv("HOME"), "data"),
        GENOME_DIR = "REFGENS",
        DOCUMENTATION_DIR = "documentation"
    ),
    
    PATTERNS = list(
        GENOME = "S288C_refgenome.fna",
        SAMPLE_TABLE = "sample_table",
        CONTROL_FACTOR_PREFIX = "X__cf_"
    ),
    
    CHROMOSOME_MAPPING = list(
        STYLES = c("UCSC", "Roman", "Numeric"),
        ROMAN = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                 "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"),
        PREFIX = "chr"
    ),
    
    REQUIRED_PACKAGES = c(
        "QuasR",
        "GenomicAlignments",
        "Gviz",
        "rtracklayer",
        "ShortRead",
        "tidyverse"
    )
)
CONFIG <- list(
    PATHS = list(
        FEATURE_DIR = file.path(Sys.getenv("HOME"), "data", "feature_files")
    ),
    
    FEATURE_TYPES = list(
        PEAKS = "eaton_peaks",
        ORIGINS = "origins",
        GENES = "genes"
    ),
    
    REQUIRED_CATEGORIES = list(
        BASE = "antibody",
        OPTIONAL = c("condition", "treatment", "timepoint")
    ),
    
    LABEL_CONFIG = list(
        SEPARATOR = "_",
        MAX_LENGTH = 50,
        TRUNCATE_SUFFIX = "..."
    )
)
#' Merge with existing configurations
CONFIG <- list(
    GENOME = list(
        DEFAULT_STRAND = "*",
        DEFAULT_CHROMOSOME = "chr1",
        STYLES = list(
            UCSC = "UCSC",
            ENSEMBL = "ENSEMBL"
        )
    ),
    
    VISUALIZATION = list(
        TRACKS = list(
            COLORS = c(
                PRIMARY = "#1f77b4",
                SECONDARY = "#ff7f0e",
                HIGHLIGHT = "#2ca02c"
            ),
            TYPES = list(
                DATA = "l",
                ANNOTATION = "gene"
            )
        ),
        DEFAULTS = list(
            GENOME = "hg19",
            MIN_HEIGHT = 0,
            MAX_HEIGHT = 100
        )
    )
)
#' Add to existing genome configurations
CONFIG <- list(
    CHROMOSOMES = list(
        SPECIAL = c("X", "Y", "MT", "M"),
        PREFIX = "chr",
        DEFAULT = list(
            NAME = "chrX",
            START = 1,
            END = 170000,
            STRAND = "*"
        )
    ),
    
    NAMING = list(
        PATTERNS = list(
            ROMAN = "^[IVXLCDM]+$",
            NUMERIC = "^\\d+$",
            PREFIX = "^chr"
        ),
        STYLES = list(
            UCSC = "UCSC",     # chrI, chrII, etc.
            ENSEMBL = "ENSEMBL", # 1, 2, etc.
            NCBI = "NCBI"      # I, II, etc.
        )
    )
)

#!/usr/bin/env Rscript

#' Project Configuration Definition
#' @export PROJECT_CONFIG
PROJECT_CONFIG <- list(
    SYSTEM = list(
        VERSION = "1.0.0",
        R_MIN_VERSION = "4.2.0",
        PLATFORM = .Platform$OS.type
    ),
    
    PATHS = list(
        ROOT = normalizePath("~/lab_utils"),
        R_BASE = file.path(normalizePath("~/lab_utils"), "R"),
        FUNCTIONS = "functions",
        SCRIPTS = "scripts",
        CONFIG = "config",
        DATA = normalizePath("~/data"),
        LOGS = normalizePath("~/logs")
    ),
    
    LOAD_SEQUENCE = list(
        CRITICAL = c(
            "logging_utils.R",
            "environment_utils.R",
            "validation_utils.R"
        )
    ),
    FILE_TYPES = list(
            NGS = list(
                EXTENSIONS = c(
                    BAM = "\\.bam$",
                    FASTQ = "\\.(fastq|fq)(\\.gz)?$",
                    BIGWIG = "\\.bw$",
                    BED = "\\.bed$",
                    NARROWPEAK = "\\.narrowPeak$",
                    MOTIF = "\\.(meme|pwm|jaspar)$"
                ),
                REQUIRED_INDEX = c(
                    BAM = "\\.bai$",
                    BIGWIG = "\\.bw\\.tbi$"
                )
            )
    ),
    VALIDATION = list(
         LIMITS = list(
                CHROMOSOME = c(min = 1, max = 16),
                READ_LENGTH = c(min = 20, max = 150),
                PEAK_SCORE = c(min = 0, max = 1000)
            ),
            REQUIRED_DIRS = c(
                "peak",
                "alignment",
                "plots",
                "documentation",
                "fastq",
                "processedFastq"
            )
        )
)


#' Package Management Configuration
CONFIG <- list(
    REPOSITORIES = list(
        CRAN = "https://cloud.r-project.org",
        BIOC = BiocManager::repositories()
    ),
    
    PACKAGES = list(
        FASTQ_ANALYSIS = c(
            "ShortRead",
            "Rsubread",
            "Biostrings",
            "dada2"
        ),
        
        TRACK_VISUALIZATION = c(
            "QuasR",
            "GenomicAlignments",
            "Gviz",
            "rtracklayer"
        ),
        
        PEAK_ANALYSIS = c(
            "ChIPseeker",
            "ChIPpeakAnno",
            "DiffBind",
            "normR",
            "mosaics",
            "csaw"
        ),
        
        MOTIF_ANALYSIS = c(
            "motifStack",
            "TFBSTools",
            "JASPAR2020",
            "universalmotif",
            "memes"
        ),
        
        VISUALIZATION = c(
            "ggbio",
            "ComplexHeatmap",
            "EnhancedVolcano",
            "ggplot2",
            "ggcoverage",
            "gggenome"
        ),
        
        STATISTICS = c(
            "DESeq2",
            "edgeR",
            "limma"
        ),
        
        CORE = c(
            "ggplot2",
            "rmarkdown",
            "knitr",
            "tidyverse",
            "furrr"
        )
    ),
    
    RENV = list(
        LOCKFILE = "renv.lock",
        LIBRARY = "renv/library"
    )
)
#' Add to existing package configurations
CONFIG <- list(
    BIOINFORMATICS = list(
        CORE = c(
            "BiocGenerics",
            "MatrixGenerics",
            "GenomeInfoDb",
            "GenomicRanges"
        ),
        
        VISUALIZATION = c(
            "ggplot2",
            "pheatmap",
            "ComplexHeatmap",
            "circlize",
            "ggbio",
            "Gviz"
        ),
        
        SEQUENCING = c(
            "Rqc",
            "ShortRead",
            "QuasR",
            "Rsubread",
            "Rsamtools",
            "Rbowtie",
            "Rbowtie2"
        ),
        
        ANALYSIS = c(
            "DESeq2",
            "RUVSeq",
            "methylKit",
            "ChIPpeakAnno",
            "normr"
        ),
        
        ANNOTATION = c(
            "AnnotationHub",
            "GenomicFeatures",
            "BSgenome",
            "BSgenome.Scerevisiae.UCSC.sacCer3"
        ),
        
        MOTIF = c(
            "MotifDb",
            "TFBSTools",
            "rGADEM",
            "JASPAR2018"
        )
    ),
    
    DEVELOPMENT = list(
        DOCUMENTATION = c(
            "rmarkdown",
            "xaringan",
            "officer",
            "quarto"
        ),
        
        STYLE = c(
            "styler",
            "formatR"
        ),
        
        TOOLS = c(
            "devtools",
            "rcrossref",
            "taskscheduleR",
            "bio3d"
        )
    ),
    
    GITHUB_PACKAGES = list(
        "crsh/citr",
        "paleolimbot/rbbt"
    )
)

#' Add to existing configurations
CONFIG <- list(
    FEATURES = list(
        PATHS = list(
            BASE_DIR = file.path(Sys.getenv("HOME"), "data", "feature_files")
        ),
        
        PATTERNS = list(
            EXCLUDE = c(
                "sample-key.tab",
                "\\.rds$",
                "_converted.bed$"
            ),
            VERIFY = c(
                "\\.rds$",
                "_converted.bed$"
            )
        ),
        
        FILE_TYPES = list(
            NUCLEOSOME = "Nucleosome_calls",
            TIMING = "hawkins",
            TRANSCRIPTION = "Rhee",
            FEATURES = "SGD"
        ),
        
        COLUMNS = list(
            NUCLEOSOME = list(
                EXCLUDE = c("Nucleosome ID", "Nucleosome dyad", "Chromosome"),
                POSITION = "Nucleosome dyad"
            ),
            TIMING = list(
                EXCLUDE = c("Chromosome", "Position"),
                WINDOW = 100  # +/- bp around position
            )
        )
    )
)

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

#' Add to existing visualization configurations
CONFIG <- list(
    MODES = list(
        INTERACTIVE = list(
            DEFAULT_DIR = "240808Bel",
            DEFAULT_CHROMOSOME = 10
        )
    ),
    
    SYNC = list(
        COMMAND = "rsync -nav %s:%s/%s/plots/* %s/%s/plots/",
        USER = "username",
        DOMAIN = "domain",
        REMOTE_BASE = "~/data",
        LOCAL_BASE = "/local/dir"
    )
)

#' Add to existing QC configurations
CONFIG <- list(
    FASTQC = list(
        PATTERNS = list(
            DATA_FILE = "fastqc_data",
            MODULE_START = "^>>",
            MODULE_END = ">>END_MODULE",
            HEADER = "^#"
        ),
        
        OUTPUT = list(
            DATE_FORMAT = "%Y-%m-%d-%H-%M-%S",
            SEPARATOR = "\t",
            EXTENSIONS = list(
                MODULE = ".tab",
                SUMMARY = "_summary.tab"
            )
        ),
        
        COLUMNS = list(
            SUMMARY = c("Stat", "Value")
        )
    ),
    
    PATHS = list(
        BASE_DIR = file.path(Sys.getenv("HOME"), "data"),
        QC_SUBDIR = "qualityControl"
    )
)
#' Add to existing QC configurations
CONFIG <- list(
    BAM_QC = list(
        PATTERNS = list(
            FLAGSTAT = "bamFlagstat",
            SAMPLE_INFO = "sample_info"
        ),
        
        METRICS = list(
            TOTAL_READS = "total.*reads",
            MAPPED_READS = "^mapped$"
        ),
        
        OUTPUT = list(
            DATE_FORMAT = "%Y-%m-%d-%H-%M-%S",
            SEPARATOR = "\t"
        )
    ),
    
    PATHS = list(
        QC_DIR = "qualityControl",
        DOC_DIR = "documentation"
    )
)

#!/usr/bin/env Rscript
# TODO: Do we need to validate the experiment_id format?
# TODO: Should we add validation for category value uniqueness?
# R/config/experiment_config.R

#' Experiment Configuration
#' @description Define experimental setup and validation
#' @export EXPERIMENT_CONFIG


EXPERIMENT_CONFIG <- list(
    METADATA = list(
        EXPERIMENT_ID = "241010Bel",
        EXPECTED_SAMPLES = 65,
        VERSION = "1.0.0"
    ),
    
    CATEGORIES = list(
        rescue_allele = c("NONE", "WT", "4R", "PS"),
        auxin_treatment = c("NO", "YES"),
        time_after_release = c("0", "1", "2"),
        antibody = c("Input", "ProtG", "HM1108", "V5", "ALFA", "UM174")
    ),
    
    INVALID_COMBINATIONS = list(
        rescue_allele_auxin_treatment = quote(rescue_allele %in% c("4R", "PS") & auxin_treatment == "NO"),
        protg_time_after_release = quote(antibody == "ProtG" & time_after_release %in% c("1", "2")),
        input_time_after_release = quote(antibody == "Input" & time_after_release %in% c("1", "2")),
        input_rescue_allele_auxin_treatment = quote(antibody == "Input" & rescue_allele %in% c("NONE", "WT") & auxin_treatment == "YES")
    ),
    
    EXPERIMENTAL_CONDITIONS = list(
        is_input = quote(time_after_release == "0" & antibody == "Input"),
        is_protg = quote(rescue_allele == "WT" & time_after_release == "0" & antibody == "ProtG" & auxin_treatment == "NO"),
        is_v5 = quote(antibody == "V5"),
        is_alfa = quote(antibody == "ALFA"),
        is_1108 = quote(antibody == "HM1108" & time_after_release == "0"),
        is_174 = quote(antibody == "UM174")
    ),
    
    COMPARISONS = list(
        comp_1108forNoneAndWT = quote(antibody == "HM1108" & rescue_allele %in% c("NONE", "WT")),
        comp_1108forNoneAndWT_auxin = quote(antibody == "HM1108" & auxin_treatment == "YES"),
        comp_timeAfterReleaseV5WT = quote(antibody == "V5" & rescue_allele == "WT" & auxin_treatment == "YES"),
        comp_timeAfterReleaseV5NoTag = quote(antibody == "V5" & rescue_allele == "NONE" & auxin_treatment == "YES"),
        comp_V5atTwoHours = quote(antibody == "V5" & time_after_release == "2" & auxin_treatment == "YES"),
        comp_UM174atTwoHours = quote(antibody == "UM174" & time_after_release == "2" & auxin_treatment == "YES"),
        comp_ALFAforNoRescueNoTreat = quote(antibody == "ALFA" & rescue_allele == "NONE" & auxin_treatment == "NO"),
        comp_ALFAforNoRescueWithTreat = quote(antibody == "ALFA" & rescue_allele == "NONE" & auxin_treatment == "YES"),
        comp_ALFAatTwoHoursForAllAlleles = quote(antibody == "ALFA" & time_after_release == "2" & auxin_treatment == "YES"),
        comp_UM174atZeroHoursForAllAlleles = quote(antibody == "UM174" & time_after_release == "0" & auxin_treatment == "YES"),
        comp_AuxinEffectOnUM174 = quote(antibody == "UM174" & time_after_release == "2" & rescue_allele %in% c("NONE", "WT"))
    ),
    
    CONTROL_FACTORS = list(
        genotype = c("rescue_allele")
    ),
    
    COLUMN_ORDER = c("antibody", "rescue_allele", "auxin_treatment", "time_after_release")
)

#' Validate Category Values
validate_experiment_categories <- function(experiment_config) {
    for (category_name in names(experiment_config$CATEGORIES)) {
        values <- experiment_config$CATEGORIES[[category_name]]
        if (!is.character(values)) {
            stop(sprintf(
                "Category '%s' values must be character vectors, got %s",
                category_name, class(values)
            ))
        }
    }
    invisible(TRUE)
}

#' Validate Column References
validate_experiment_column_references <- function(experiment_config) {
    # Get valid column names
    valid_columns <- names(experiment_config$CATEGORIES)
    
    # Check comparison groups
    for (comp_name in names(experiment_config$COMPARISONS)) {
        comp_expr <- experiment_config$COMPARISONS[[comp_name]]
        comp_vars <- all.vars(comp_expr)
        invalid_cols <- setdiff(comp_vars, valid_columns)
        if (length(invalid_cols) > 0) {
            stop(sprintf(
                "Invalid columns in comparison '%s': %s",
                comp_name, paste(invalid_cols, collapse = ", ")
            ))
        }
    }
    
    # Check control factors
    for (factor_name in names(experiment_config$CONTROL_FACTORS)) {
        invalid_cols <- setdiff(
            experiment_config$CONTROL_FACTORS[[factor_name]],
            valid_columns
        )
        if (length(invalid_cols) > 0) {
            stop(sprintf(
                "Invalid columns in control factor '%s': %s",
                factor_name, paste(invalid_cols, collapse = ", ")
            ))
        }
    }
    
    # Check experimental conditions
    for (cond_name in names(experiment_config$EXPERIMENTAL_CONDITIONS)) {
        cond_expr <- experiment_config$EXPERIMENTAL_CONDITIONS[[cond_name]]
        cond_vars <- all.vars(cond_expr)
        invalid_cols <- setdiff(cond_vars, valid_columns)
        if (length(invalid_cols) > 0) {
            stop(sprintf(
                "Invalid columns in condition '%s': %s",
                cond_name, paste(invalid_cols, collapse = ", ")
            ))
        }
    }
    
    invisible(TRUE)
}

#' Validate Column Order
validate_experiment_column_order <- function(experiment_config) {
    if (!identical(
        sort(names(experiment_config$CATEGORIES)),
        sort(experiment_config$COLUMN_ORDER)
    )) {
        stop("Column order must include all category columns")
    }
    invisible(TRUE)
}

#' Generate Experiment Sample Grid
#' @param project_config List Project configuration
#' @param experiment_config List Experiment configuration
#' @param init_logging Logical Whether to initialize logging
#' @return data.frame Filtered experiment grid
generate_experiment_grid <- function(
    experiment_config = EXPERIMENT_CONFIG,
    init_logging = TRUE
) {
    # Initialize logging if requested
    if (init_logging) {
        if (file.exists("~/lab_utils/R/functions/logging_utils.R")) {
            source("~/lab_utils/R/functions/logging_utils.R")
            log_file <- initialize_logging(script_name = "experiment_grid")
            log_info("Starting experiment grid generation", log_file)
        } else {
            warning("logging_utils.R not found. Proceeding without logging.")
        }
    }
    
    tryCatch({
        # Run validations
        validate_experiment_categories(experiment_config)
        validate_experiment_column_references(experiment_config)
        validate_experiment_column_order(experiment_config)

        # Generate combinations
        grid <- do.call(expand.grid, experiment_config$CATEGORIES)
        
        # Filter invalid combinations
        invalid_idx <- Reduce(
            `|`, 
            lapply(experiment_config$INVALID_COMBINATIONS, eval, envir = grid)
        )
        grid <- subset(grid, !invalid_idx)
        
        # Apply experimental conditions
        valid_idx <- Reduce(
            `|`, 
            lapply(experiment_config$EXPERIMENTAL_CONDITIONS, eval, envir = grid)
        )
        grid <- subset(grid, valid_idx)
        
        # Verify sample count
        n_samples <- nrow(grid)
        if (n_samples != experiment_config$METADATA$EXPECTED_SAMPLES) {
            warning(sprintf(
                "Expected %d samples, got %d", 
                experiment_config$METADATA$EXPECTED_SAMPLES, 
                n_samples
            ))
        }
        
        # Sort columns
        grid <- grid[, experiment_config$COLUMN_ORDER]
        
        # Add attributes
        attr(grid, "control_factors") <- experiment_config$CONTROL_FACTORS
        attr(grid, "experiment_id") <- experiment_config$METADATA$EXPERIMENT_ID
        
        if (init_logging) {
            log_info(sprintf("Generated grid with %d samples", nrow(grid)))
        }
        
        return(grid)
        
    }, error = function(e) {
        msg <- sprintf("Failed to generate experiment grid: %s", e$message)
        if (init_logging) {
            log_error(msg)
        }
        stop(msg)
    })
}
#    PATTERNS = list(
#        R_FILES = "\\.R$",
#        EXCLUDE = c("^\\.", "^_", "test_", "example_")
#    ),
#
#    ENVIRONMENT = list(
#        REQUIRED_VARS = c(
#            "HOME",
#            "R_LIBS_USER"
#        ),
#        RENV_SETTINGS = list(
#            AUTO_ACTIVATE = TRUE,
#            SNAPSHOT_INTERVAL = 86400  # 24 hours
