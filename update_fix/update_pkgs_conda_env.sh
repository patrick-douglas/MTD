#!/bin/bash -x

conda activate halla0820
~/MTD/update_fix/Install.R.packages.sh
conda deactivate
exit 1
conda activate MTD
~/MTD/update_fix/Install.R.packages.sh
conda deactivate

conda activate R412
~/MTD/update_fix/Install.R.packages.sh
conda deactivate
