#Check Dependencies installation
#!/bin/bash
w=$(tput sgr0) 
r=$(tput setaf 1)
g=$(tput setaf 2) 
y=$(tput setaf 3) 
p=$(tput setaf 5) 
R_ver=`R --version | grep version | grep R | awk '{print $3}'`
echo "${g}***********************************************"
echo "${w}R $R_ver packages versions"
echo "${g}Env:${w} R412"
echo "${g}***********************************************${g}"
echo -n "${g}tidyverse: 	                       	${w}" ; R --no-restore -e 'packageVersion("tidyverse")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo "${g}***********************************************${w}"
