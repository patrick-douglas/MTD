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
echo "${g}Env:${w} halla0820"
echo "${g}***********************************************${g}"
echo -n "${g}lattice:                      	        ${w}" ; R --no-restore -e 'packageVersion("lattice")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}MASS: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("MASS")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}psychTools: 	                       	${w}" ; R --no-restore -e 'packageVersion("psychTools")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}Matrix:                       	        ${w}" ; R --no-restore -e 'packageVersion("Matrix")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}XICOR: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("XICOR")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}mclust:                       	        ${w}" ; R --no-restore -e 'packageVersion("mclust")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}BiocManager:                   	        ${w}" ; R --no-restore -e 'packageVersion("BiocManager")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}eva: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("eva")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo "${g}***********************************************${w}"
