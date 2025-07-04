#!/bin/bash -x
threads=`nproc`

source ~/miniconda3/etc/profile.d/conda.sh

#conda init
#conda activate MTD
conda install -y -n MTD conda-forge::r-textshaping
conda install -y -n MTD conda-forge::r-ragg
conda install -y -n MTD -c conda-forge freetype
conda install -y -n MTD -c conda-forge pkg-config
conda install -y -n MTD -c conda-forge r-tidyverse
conda install -y -n MTD -c conda-forge r-car
conda install -y -n MTD -c conda-forge r-rstatix
conda install -y -n MTD -c conda-forge r-ggpubr
conda install -y -n MTD -c conda-forge r-plyr

R -e "install.packages(c('BiocManager'), repos='http://cran.us.r-project.org', Ncpus=$threads)"
R -e "install.packages(c('plyr'), repos='http://cran.us.r-project.org', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/rlang_1.1.2.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/vctrs_0.6.4.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/BiocGenerics_0.50.0.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/S4Vectors_0.42.1.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/IRanges_2.38.1.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/UCSC.utils_1.0.0.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/GenomeInfoDbData_1.2.12.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/GenomeInfoDb_1.40.1.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/matrixStats_1.3.0.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/formatR_1.14.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/lambda.r_1.2.4.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/futile.options_1.0.1.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/futile.logger_1.4.3.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/RColorBrewer_1.1-3.tar.gz', repos=NULL, type='source', Ncpus=$threads)"

conda deactivate
