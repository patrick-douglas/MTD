#!/bin/bash -x
#R412
conda install -n R412 -y -c conda-forge pkg-config
conda install -y -n R412 python=3.10
conda run -n R412 ~/MTD/update_fix/Install.R.packages.R412.sh
conda install -n R412 -y -c conda-forge bioconda::glpk

#py2
conda install -n py2 -y -c conda-forge pkg-config



