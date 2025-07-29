#R -e 'install.packages("~/MTD/update_fix/pvr_pkg/BiocManager_1.30.22.tar.gz", repos=NULL, type="source")'
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repos = "http://cran.us.r-project.org")
BiocManager::install(version = "3.14")
#Below fixes the DESeq2 incompatibility issue that causes the pkg fails to install 
BiocManager::install('BiocGenerics', force = TRUE)
BiocManager::install('genefilter', force = TRUE)
#Above fixes the DESeq2 incompatibility issue that causes the pkg fails to install
R -e "install.packages('~/MTD/update_fix/pvr_pkg/Matrix_1.6-5.tar.gz', repos=NULL, type='source')"
BiocManager::install(c("biomaRt","DESeq2","tximeta","limma","phyloseq","glmGamPoi","cmapR","MAST","microbiome","ANCOMBC","Maaslin2","DO.db","clusterProfiler","enrichplot","pathview"))
#R -e 'install.packages("~/MTD/update_fix/pvr_pkg/pacman_0.5.1.tar.gz", repos=NULL, type="source")'
if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(tidyverse,ggrepel,colorspace,RColorBrewer,
                pheatmap,VennDiagram,doParallel,foreach,
                stringi,vegan,ggpubr,reshape2,sctransform,hdf5r,
                ggridges,ggnewscale,ggupset,Seurat)

#included in tidyverse: ggplot2,tidyr,dplyr,stringr,tibble,
