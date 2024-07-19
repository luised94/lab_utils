
#If running on Bash script, you must add the R and gnu modules 
#Had to run this manually first. May have to include this line next time I try to set it up. 
# r = getOption("repos")
# r["CRAN"] = "https://cloud.r-project.org/"
# options(repos = r)

Sys.time()

install.packages(c("tidyverse", "R.utils", "ggplot2","BiocManager"), 
                 lib = "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2", )

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

library(BiocManager)
BiocManager::install(pkgs=bioinformatics_packages)

package_to_check <- c("Rbowtie2", "BSgenome", "R.utils","ShortRead","QuasR","BiocManager")
package_was_loaded <- unlist(lapply(package_to_check, library, character.only = TRUE, logical.return=TRUE)) #lib.loc = "/home/luised94/R/x86_64-pc-linux-gnu-library/4.2"
if (length(lapply(package_to_check[!package_was_loaded], function(x){
  paste(x, "Package did not install ")
})) == 0) print("All packages loaded.")


Sys.time()
print("package-installation.R complete")

