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
