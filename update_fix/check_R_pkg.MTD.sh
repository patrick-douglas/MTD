#!/bin/bash
w=$(tput sgr0) 
r=$(tput setaf 1)
g=$(tput setaf 2) 
y=$(tput setaf 3) 
p=$(tput setaf 5) 
R_ver=`R --version | grep version | grep R | awk '{print $3}'`
echo "${g}***************************"
echo "${w}R $R_ver packages versions"
echo "${g}Env:${w} MTD"
echo "${g}***************************${g}"

echo -n "${g}tidyverse:            ${w}" ; R --no-restore -e 'packageVersion("tidyverse")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}car:                  ${w}" ; R --no-restore -e 'packageVersion("car")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}rstatix:              ${w}" ; R --no-restore -e 'packageVersion("rstatix")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ggpubr:               ${w}" ; R --no-restore -e 'packageVersion("ggpubr")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}plyr:                 ${w}" ; R --no-restore -e 'packageVersion("plyr")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}textshaping:          ${w}" ; R --no-restore -e 'packageVersion("textshaping")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ragg:                 ${w}" ; R --no-restore -e 'packageVersion("ragg")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}BiocManager:          ${w}" ; R --no-restore -e 'packageVersion("BiocManager")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}rlang:                ${w}" ; R --no-restore -e 'packageVersion("rlang")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}vctrs:                ${w}" ; R --no-restore -e 'packageVersion("vctrs")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}BiocGenerics:         ${w}" ; R --no-restore -e 'packageVersion("BiocGenerics")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}S4Vectors:            ${w}" ; R --no-restore -e 'packageVersion("S4Vectors")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}IRanges:              ${w}" ; R --no-restore -e 'packageVersion("IRanges")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}UCSC.utils:           ${w}" ; R --no-restore -e 'packageVersion("UCSC.utils")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}GenomeInfoDbData:     ${w}" ; R --no-restore -e 'packageVersion("GenomeInfoDbData")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}GenomeInfoDb:         ${w}" ; R --no-restore -e 'packageVersion("GenomeInfoDb")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}matrixStats:          ${w}" ; R --no-restore -e 'packageVersion("matrixStats")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}formatR:              ${w}" ; R --no-restore -e 'packageVersion("formatR")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}lambda.r:             ${w}" ; R --no-restore -e 'packageVersion("lambda.r")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}futile.options:       ${w}" ; R --no-restore -e 'packageVersion("futile.options")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}futile.logger:        ${w}" ; R --no-restore -e 'packageVersion("futile.logger")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}RColorBrewer:         ${w}" ; R --no-restore -e 'packageVersion("RColorBrewer")' | grep packageVersion -A 1 | grep "[1]" | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'

echo "${g}***************************${w}"
