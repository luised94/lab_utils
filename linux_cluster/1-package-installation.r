#STATUS: REMOVE.
#https://rdrr.io/category/biocview/GenomeAnnotation/
#https://uclouvain-cbio.github.io/WSBIM1322/index.html
#https://compgenomr.github.io/book/
#https://rockefelleruniversity.github.io/RU_ChIPseq/index.html
#https://learn.gencore.bio.nyu.edu/

#Loading packages in a clean way. 
# Install packages that will be used most of the time.
# constant_packages <- c("renv","pacman")
# install.packages(constant_packages)
# 
# 
# pacman::p_load(ggplot2, tidyverse,BiocManager,stringr,R.utils)
renv::install("styler")
bioinformatics_packages <- c("BiocGenerics","MatrixGenerics",'qvalue','plot3D','ggplot2','pheatmap','cowplot',
  'cluster', 'NbClust', 'fastICA', 'NMF','matrixStats',
  'Rtsne', 'mosaic', 'knitr', 'genomation',
  'ggbio', 'Gviz', 'DESeq2', 'RUVSeq',
  'gProfileR', 'ggfortify', 'corrplot',
  'gage', 'EDASeq', 'formatR', 'BiocFileCache',
  'svglite', 'Rqc', 'ShortRead', 'QuasR',
  'methylKit','FactoMineR', 'iClusterPlus',
  'enrichR','caret','xgboost','glmnet',
  'DALEX','kernlab','pROC','nnet','RANN',
  'ranger','GenomeInfoDb', 'GenomicRanges',
  'GenomicAlignments', 'ComplexHeatmap', 'circlize', 
  'rtracklayer', 'tidyr', 'dplyr',
  'AnnotationHub', 'GenomicFeatures', 'normr',
  'MotifDb', 'TFBSTools', 'rGADEM', 'JASPAR2018', 
  'BSgenome', 'htmltab', 'usethis',
  'Rsubread', 'Rsamtools', 'Rbowtie', 'Rbowtie2' , "ChIPpeakAnno", 
  "seqinr","GenomeInfoDbData","BSgenome.Scerevisiae.UCSC.sacCer3")

BiocManager::install(bioinformatics_packages)
#Cool way to install Stolen/borrowed from https://www.r-bloggers.com/2020/01/an-efficient-way-to-install-and-load-r-packages-2/
#  
# packages <- c("ggplot2", "readxl", "dplyr", "tidyr", "ggfortify", "DT", "reshape2", "knitr", "lubridate", "pwr", "psy", "car", "doBy", "imputeMissings", "RcmdrMisc", "questionr", "vcd", "multcomp", "KappaGUI", "rcompanion", "FactoMineR", "factoextra", "corrplot", "ltm", "goeveg", "corrplot", "FSA", "MASS", "scales", "nlme", "psych", "ordinal", "lmtest", "ggpubr", "dslabs", "stringr", "assist", "ggstatsplot", "forcats", "styler", "remedy", "snakecaser", "addinslist", "esquisse", "here", "summarytools", "magrittr", "tidyverse", "funModeling", "pander", "cluster", "abind")
# 
# # Install packages not yet installed
# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }

# 
# # Packages loading
# invisible(lapply(packages, library, character.only = TRUE))


#May need to change some variables.
#Run this, after creating virtual environment and directory.
# #Install packages for documentations and parallel processing 
# renv::install(c("tidyverse", "R.utils", "doParallel", "snow","fs","rmarkdown", "gt","formattable","ggplot2", "webshot2",
#               "rmarkdown","xaringan","officer",
#               "quarto","BiocManager"))
# #Install some files for development
# renv::install("devtools","rcrossref","taskscheduleR","bio3d")
# devtools::install_github("crsh/citr")
# remotes::install_github("paleolimbot/rbbt")

#Non Bioconductor package 
#install.packages("formattable")

#Enable SSL on windows: https://answers.microsoft.com/en-us/windows/forum/all/ssl-error-preventing-connection-in-windows-10/192436e5-e37b-4b4c-a4e0-c4ec744b0f5c
#Then install basilisk.utils
# BiocManager::install(c("basilisk.utils", "MACSr"))
# BiocManager::install(c("nanopoRe"))
# update.packages()
# lapply(paste0("package:", names(sessionInfo()$otherPkgs)),   # Unload add-on packages
#        detach,
#        character.only = TRUE, unload = TRUE)


deattachAll <- function(){
  lapply(paste0("package:", names(sessionInfo()$otherPkgs)),   # Unload add-on packages
                detach,
                character.only = TRUE, unload = TRUE)
}



