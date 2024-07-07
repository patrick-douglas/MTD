#!/bin/bash
R -e "install.packages('lattice',repos = 'http://cran.us.r-project.org')"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/Matrix_1.6-5.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages('~/MTD/update_fix/pvr_pkg/MASS_7.3-60.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
R -e "install.packages(c('XICOR','mclust','BiocManager'), repos='http://cran.us.r-project.org', Ncpus=$threads)"
R -e "BiocManager::install('preprocessCore', ask = FALSE)"
R -e "install.packages('eva', INSTALL_opts = '--no-lock', repos='http://cran.us.r-project.org', Ncpus=$threads)"

