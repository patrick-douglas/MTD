#!/bin/bash -x

R -e 'install.packages("httr",repos = "http://cran.us.r-project.org")'
R -e 'BiocManager::install("SummarizedExperiment",repos = "http://cran.us.r-project.org")'
R -e 'install.packages("~/MTD/update_fix/pvr_pkg/Matrix_1.6-5.tar.gz", repos=NULL, type="source")'
R -e 'BiocManager::install("ggplot2",repos = "http://cran.us.r-project.org")'
R -e 'install.packages("ade4",repos = "http://cran.us.r-project.org")'
R -e 'install.packages("biomformat",repos = "http://cran.us.r-project.org")'
R -e 'BiocManager::install("igraph", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("multtest", force = TRUE, update = FALSE)'
R -e 'install.packages("vegan",repos = "http://cran.us.r-project.org")'
R -e 'install.packages("~/MTD/update_fix/pvr_pkg/RProtoBufLib_2.14.1.tar.gz", repos=NULL, type="source")'
R -e 'remotes::install_github("RGLab/cytolib",repos = "http://cran.us.r-project.org")'
R -e 'remotes::install_github("RGLab/flowcore",repos = "http://cran.us.r-project.org")'
R -e 'BiocManager::install("AnnotationDbi", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("ggfun", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("SingleCellExperiment", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("phyloseq", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("curl", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("annotate", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("BiocFileCache", force = TRUE, update = FALSE)'
R -e 'install.packages("~/MTD/update_fix/pvr_pkg/MASS_7.3-60.tar.gz", repos=NULL, type="source")'
R -e 'BiocManager::install("aplot", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("lme4", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("TMB", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("GenomicFeatures", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("metagenomeSeq", force = TRUE, update = FALSE)'
R -e 'remotes::install_github("YuLab-SMU/ggtree",repos = "http://cran.us.r-project.org")'
R -e 'remotes::install_github("YuLab-SMU/enrichplot",repos = "http://cran.us.r-project.org")'
R -e 'BiocManager::install("ragg", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("sctransform", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("cmapR", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("clusterProfiler", force = TRUE, update = FALSE)'
#fix httpd installations
rm -rf /home/me/miniconda3/envs/R412/lib/R/library/*LOCK-httpuv/
R -e 'install.packages("~/MTD/update_fix/pvr_pkg/httpuv_1.6.0.tar.gz", repos=NULL, type="source")'
R -e 'BiocManager::install("miniUI", force = TRUE, update = FALSE)'
update_fix $ R -e 'BiocManager::install("shiny", force = TRUE, update = FALSE)'
R -e 'BiocManager::install("Seurat", force = TRUE, update = FALSE)'
R -e 'install.packages("Seurat",repos = "http://cran.us.r-project.org")'
R -e 'BiocManager::install("tximeta", force = TRUE, update = FALSE)'