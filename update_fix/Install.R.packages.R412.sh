#!/bin/bash -x

# Número de threads disponíveis
threads=$(nproc)

# Inicializa o Conda no ambiente shell (obrigatório fora de shell interativo)
source ~/miniconda3/etc/profile.d/conda.sh

# Ativa o ambiente Conda
conda activate R412
R -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")'
R -e 'BiocManager::install(version = "3.14", ask = FALSE)'
conda install -y -n R412 -c conda-forge r-httr
R -e 'install.packages("httr", repos = "https://cloud.r-project.org", dependencies = TRUE)'
R -e 'install.packages("lattice",repos = "https://cloud.r-project.org")'

# instalar Matrix compatível com R 4.1.2
R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.tar.gz", repos = NULL, type = "source")'
R -e 'remotes::install_bioc("SummarizedExperiment", version = "3.14")'
R -e "install.packages('~/MTD/update_fix/pvr_pkg/MASS_7.3-60.tar.gz', repos=NULL, type='source')"
R -e "BiocManager::install('ggplot2')"
conda install -n R412 -c conda-forge -y r-ade4
R -e "install.packages('ade4',repos = 'https://cloud.r-project.org')"
R -e "BiocManager::install('biomformat')"
conda install -n R412 -c conda-forge -y r-igraph
R -e "BiocManager::install('multtest', force = TRUE, update = FALSE)"
R -e "install.packages('vegan',repos = 'https://cloud.r-project.org')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/RProtoBufLib_2.14.1.tar.gz', repos=NULL, type='source')"
R -e "remotes::install_github('RGLab/cytolib',repos = 'https://cloud.r-project.org')"
R -e 'BiocManager::install("flowCore", ask = FALSE, update = TRUE)'
R -e "install.packages('~/MTD/update_fix/pvr_pkg/DBI_1.2.3.tar.gz', repos=NULL, type='source')"
R -e "BiocManager::install('AnnotationDbi', force = TRUE, update = FALSE)"
R -e "BiocManager::install('SingleCellExperiment', force = TRUE, update = FALSE)"
R -e "BiocManager::install('phyloseq', force = TRUE, update = FALSE)"
R -e "BiocManager::install('curl', force = TRUE, update = FALSE)"
R -e "BiocManager::install('annotate', force = TRUE, update = FALSE)"
R -e "BiocManager::install('BiocFileCache', force = TRUE, update = FALSE)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/yulab.utils_0.1.5.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/ggfun_0.1.5.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/yulab.utils_0.1.9.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/gtable_0.3.6.tar.gz', repos=NULL, type='source')"
R -e "BiocManager::install('aplot', force = TRUE, update = FALSE)"
R -e "BiocManager::install('lme4', force = TRUE, update = FALSE)"
R -e "BiocManager::install('TMB', force = TRUE, update = FALSE)"
R -e "BiocManager::install('GenomicFeatures', force = TRUE, update = FALSE)"
R -e "BiocManager::install('metagenomeSeq', force = TRUE, update = FALSE)"

conda install -y -n R412 -c conda-forge r-systemfonts
conda install -y -n R412 -c conda-forge r-textshaping
conda install -y -n R412 -c conda-forge r-ragg
R -e 'remotes::install_bioc("biomaRt")'

R -e "install.packages('~/MTD/update_fix/pvr_pkg/yulab.utils_0.1.9.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/ggfun_0.1.7.tar.gz', repos=NULL, type='source')"

conda install -y -n R412 -c conda-forge libcurl freetype libpng jpeg libtiff pkg-config
conda install -y -n R412 -c conda-forge r-curl r-ragg
R -e "remotes::install_github('YuLab-SMU/ggtree',repos = 'https://cloud.r-project.org')"
R -e 'remotes::install_github("ctlab/fgsea")'
R -e "install.packages('~/MTD/update_fix/pvr_pkg/R.methodsS3_1.8.2.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/R.oo_1.27.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/R.utils_2.13.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/AnnotationDbi_1.66.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/GO.db_3.19.1.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/GOSemSim_2.30.2.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/AnnotationDbi_1.66.0.tar.gz', repos=NULL, type='source')"
R -e 'BiocManager::install("DOSE", version="3.14", ask=FALSE, update=FALSE)'
R -e 'install.packages("ggtangle", repos = c("https://guangchuangyu.r-universe.dev", "https://cloud.r-project.org"))'
R -e 'install.packages("scatterpie", repos="https://cloud.r-project.org")'
R -e "install.packages('~/MTD/update_fix/pvr_pkg/GOSemSim_2.30.2.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/ggnewscale_0.5.2.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/GOSemSim_2.34.0.tar.gz', repos=NULL, type='source')"

#conda install -y -n R412 bioconda::bioconductor-enrichplot
conda install -y -n R412 conda-forge::r-tidyverse
#R -e "BiocManager::install('ragg', force = TRUE, update = FALSE)"
R -e "BiocManager::install('sctransform', force = TRUE, update = FALSE)"
R -e "BiocManager::install('cmapR', force = TRUE, update = FALSE)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/GO.db_3.19.1.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/GOSemSim_2.30.2.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/DOSE_3.30.4.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/tidytree_0.4.2.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/tidytree_0.4.5.tar.gz', repos=NULL, type='source')"
R -e "BiocManager::install('clusterProfiler', force = TRUE, update = FALSE)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/gson_0.1.0.tar.gz', repos=NULL, type='source')"
R -e 'install.packages("processx", repos="https://cloud.r-project.org")'
R -e 'BiocManager::install("clusterProfiler", version="3.14", ask=FALSE, update=FALSE)'
R -e "install.packages('~/MTD/update_fix/pvr_pkg/lifecycle_1.0.3.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/tidyselect_1.2.1.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/rlang_1.1.2.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/vctrs_0.6.4.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/cli_3.6.2.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/purrr_1.0.1.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/tibble_3.2.1.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/cmapR_1.16.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/dplyr_1.1.4.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/dbplyr_2.3.4.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/biomaRt_2.60.1.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/tidyr_1.3.0.tar.gz', repos=NULL, type='source')"
conda install -y -n R412 -c conda-forge r-ggpubr
R -e "install.packages('~/MTD/update_fix/pvr_pkg/pbkrtest_0.5.2.tar.gz', repos=NULL, type='source')"
R -e "install.packages('car',repos = 'https://cloud.r-project.org')"
R -e "BiocManager::install('Maaslin2', force = TRUE, update = FALSE)"
#fix httpd installations
rm -rf /home/me/miniconda3/envs/R412/lib/R/library/*LOCK-httpuv/
R -e "install.packages('promises',repos = 'https://cloud.r-project.org')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/httpuv_1.6.0.tar.gz', repos=NULL, type='source')"
R -e "BiocManager::install('miniUI', force = TRUE, update = FALSE)"
R -e "BiocManager::install('shiny', force = TRUE, update = FALSE)"
conda install -y -n R412 conda-forge::r-htmltools
R -e "install.packages('~/MTD/update_fix/pvr_pkg/fastmap_1.2.0.tar.gz', repos=NULL, type='source')"
R -e 'remotes::install_version("Matrix", version="1.5-0", dependencies=TRUE, repos="https://cloud.r-project.org")'
sudo apt-get update
sudo apt-get install libgeos-dev -y
conda install -y -c conda-forge r-rgeos
R -e 'remotes::install_version("SeuratObject", version="4.1.3", repos="https://cloud.r-project.org")'
R -e 'remotes::install_version("Seurat", version="4.3.0", repos="https://cloud.r-project.org")'
R -e "BiocManager::install('tximeta', force = TRUE, update = FALSE)"
R -e "install.packages(c('BiocManager'), repos='https://cloud.r-project.org')"
R -e "install.packages(c('plyr'), repos='https://cloud.r-project.org')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/BiocGenerics_0.50.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/S4Vectors_0.42.1.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/IRanges_2.38.1.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/UCSC.utils_1.0.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/GenomeInfoDbData_1.2.12.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/GenomeInfoDb_1.40.1.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/matrixStats_1.3.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/lambda.r_1.2.4.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/futile.options_1.0.1.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/futile.logger_1.4.3.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/RColorBrewer_1.1-3.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/DO.db_2.9.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/forcats_1.0.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/tidyverse_2.0.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/conflicted_1.1.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/tidytree_0.4.2.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/AnnotationDbi_1.70.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/org.Hs.eg.db_3.21.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/graph_1.86.0.tar.gz', repos=NULL, type='source')"
R -e 'BiocManager::install("Rgraphviz")'
R -e 'BiocManager::install("KEGGgraph")'
R -e "install.packages('~/MTD/update_fix/pvr_pkg/pathview_1.44.0.tar.gz', repos=NULL, type='source')"
R -e "install.packages('prettydoc',repos = 'https://cloud.r-project.org')"
#Fix for DESeq
R -e "BiocManager::install('BiocGenerics', force = TRUE)"
R -e "BiocManager::install('genefilter', force = TRUE)"
R -e "BiocManager::install('DESeq2', force = TRUE, update = FALSE)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/AnnotationForge_1.4.4.tar.gz', repos=NULL, type='source')"

#Extra
R -e 'BiocManager::install("glmGamPoi", version="3.14", ask=FALSE, update=FALSE)'
R -e 'BiocManager::install("MAST", version="3.14", ask=FALSE, update=FALSE)'
R -e 'BiocManager::install("microbiome", version="3.14", ask=FALSE, update=FALSE)'
R -e 'install.packages("pacman", repos="https://cloud.r-project.org")'
R -e 'install.packages("VennDiagram", repos="https://cloud.r-project.org")'
R -e 'install.packages("hdf5r", repos="https://cloud.r-project.org")'
R -e 'install.packages("ggupset", repos="https://cloud.r-project.org")'
R -e 'install.packages("doParallel", repos="https://cloud.r-project.org", ask=FALSE)'
R -e 'BiocManager::install("ANCOMBC", ask = FALSE, update = FALSE)'
