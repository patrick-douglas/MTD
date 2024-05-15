#!/bin/bash -x

#MTD
conda install -n MTD -y -c conda-forge pkg-config
conda install -n MTD -y -c bioconda metaphlan=3.0.7=pyh7b7c402_0 
#py2
conda install -n py2 -y -c conda-forge pkg-config

#R412
conda install -n R412 -y -c conda-forge pkg-config
conda install -y -n R412 python=3.10
conda run -n R412 ~/MTD/update_fix/Install.R.packages.R412.sh

#halla0820
conda install -n halla0820 -y -c conda-forge pkg-config
#conda run -n R412 ~/MTD/update_fix/Install.R.packages.halla0820.sh
