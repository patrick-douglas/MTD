#!/bin/bash -x
threads=`nproc`
conda activate MTD
conda install -y -n MTD conda-forge::r-textshaping
conda install -y -n MTD conda-forge::r-ragg
conda install -y -n MTD -c conda-forge freetype
conda install -y -n MTD -c conda-forge pkg-config
conda install -y -n MTD -c conda-forge freetype pkg-config
conda install -y -n MTD -c conda-forge r-tidyverse
conda install -y -n MTD -c conda-forge r-car
conda install -y -n MTD -c conda-forge r-rstatix
conda install -y -n MTD -c conda-forge r-ggpubr
conda install -y -n MTD -c conda-forge r-plyr

R -e "install.packages(c('BiocManager'), repos='http://cran.us.r-project.org', Ncpus=$threads)"
R -e "install.packages(c('plyr'), repos='http://cran.us.r-project.org', Ncpus=$threads)"

conda deactivate
