#!/bin/bash
conda install -y -c conda-forge pkg-config
conda activate halla0820
./Install.R.packages.sh
conda deactivate
./Install.R.packages.sh
conda activate MTD
./Install.R.packages.sh
conda deactivate
conda activate R412
./Install.R.packages.sh
conda deactivate
