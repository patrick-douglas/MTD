#!/bin/bash -x
threads=`nproc`
conda activate MTD
conda install -y -n MTD conda-forge::r-textshaping
conda install -y -n MTD conda-forge::r-ragg
conda install -y -n MTD -c conda-forge freetype
conda install -y -n MTD -c conda-forge pkg-config
conda install -y -n MTD -c conda-forge freetype pkg-config
conda install -y -n MTD -c conda-forge r-tidyverse
R -e "install.packages(c('BiocManager'), repos='http://cran.us.r-project.org', Ncpus=$threads)"
conda deactivate
