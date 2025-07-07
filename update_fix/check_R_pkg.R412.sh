#Check Dependencies installation
#!/bin/bash
w=$(tput sgr0) 
r=$(tput setaf 1)
g=$(tput setaf 2) 
y=$(tput setaf 3) 
p=$(tput setaf 5) 

source ~/miniconda3/etc/profile.d/conda.sh

# Ativa o ambiente Conda
conda activate R412

R_ver=`R --version | grep version | grep R | awk '{print $3}'`
echo "${g}***********************************************"
echo "${w}R $R_ver packages versions"
echo "${g}Env:${w} R412"
echo "${g}***********************************************${g}"
echo -n "${g}lattice:                      	        ${w}" ; R --no-restore -e 'packageVersion("lattice")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}MASS: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("MASS")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}Matrix:                       	        ${w}" ; R --no-restore -e 'packageVersion("Matrix")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}BiocManager:                   	        ${w}" ; R --no-restore -e 'packageVersion("BiocManager")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}httr: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("httr")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}SummarizedExperiment:         	        ${w}" ; R --no-restore -e 'packageVersion("SummarizedExperiment")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}Matrix: 	                       	${w}" ; R --no-restore -e 'packageVersion("Matrix")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ggplot2: 	                       	${w}" ; R --no-restore -e 'packageVersion("ggplot2")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ade4: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("ade4")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}biomformat: 	                       	${w}" ; R --no-restore -e 'packageVersion("biomformat")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}igraph: 	                       	${w}" ; R --no-restore -e 'packageVersion("igraph")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}multtest: 	                       	${w}" ; R --no-restore -e 'packageVersion("multtest")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}RProtoBufLib: 	                       	${w}" ; R --no-restore -e 'packageVersion("RProtoBufLib")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}cytolib: 	                       	${w}" ; R --no-restore -e 'packageVersion("cytolib")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}flowCore: 	                       	${w}" ; R --no-restore -e 'packageVersion("flowCore")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}AnnotationDbi: 	                       	${w}" ; R --no-restore -e 'packageVersion("AnnotationDbi")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}igraph: 	                       	${w}" ; R --no-restore -e 'packageVersion("igraph")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ggfun:   	                       	${w}" ; R --no-restore -e 'packageVersion("ggfun")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}SingleCellExperiment:                	${w}" ; R --no-restore -e 'packageVersion("SingleCellExperiment")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}phyloseq: 	                       	${w}" ; R --no-restore -e 'packageVersion("phyloseq")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}curl: 	                          	${w}" ; R --no-restore -e 'packageVersion("curl")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}annotate: 	                       	${w}" ; R --no-restore -e 'packageVersion("annotate")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}BiocFileCache: 	                       	${w}" ; R --no-restore -e 'packageVersion("BiocFileCache")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}MASS: 	                          	${w}" ; R --no-restore -e 'packageVersion("MASS")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}aplot: 	                          	${w}" ; R --no-restore -e 'packageVersion("aplot")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}lme4: 	                          	${w}" ; R --no-restore -e 'packageVersion("lme4")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}TMB: 	                          	${w}" ; R --no-restore -e 'packageVersion("TMB")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}GenomicFeatures: 	                ${w}" ; R --no-restore -e 'packageVersion("GenomicFeatures")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}metagenomeSeq: 	                       	${w}" ; R --no-restore -e 'packageVersion("metagenomeSeq")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ggtree: 	                       	${w}" ; R --no-restore -e 'packageVersion("ggtree")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ragg: 	                          	${w}" ; R --no-restore -e 'packageVersion("ragg")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}sctransform: 	                        ${w}" ; R --no-restore -e 'packageVersion("sctransform")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}cmapR: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("cmapR")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}httpuv: 	                       	${w}" ; R --no-restore -e 'packageVersion("httpuv")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}miniUI: 	                       	${w}" ; R --no-restore -e 'packageVersion("miniUI")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}shiny: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("shiny")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}biomaRt: 	                       	${w}" ; R --no-restore -e 'packageVersion("biomaRt")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}tximeta: 	                       	${w}" ; R --no-restore -e 'packageVersion("tximeta")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}limma: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("limma")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}phyloseq: 	                       	${w}" ; R --no-restore -e 'packageVersion("phyloseq")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}glmGamPoi: 	                       	${w}" ; R --no-restore -e 'packageVersion("glmGamPoi")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}MAST: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("MAST")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}microbiome: 	                       	${w}" ; R --no-restore -e 'packageVersion("microbiome")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ANCOMBC: 	                       	${w}" ; R --no-restore -e 'packageVersion("ANCOMBC")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}Maaslin2: 	                       	${w}" ; R --no-restore -e 'packageVersion("Maaslin2")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}DO.db: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("DO.db")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}clusterProfiler: 	                ${w}" ; R --no-restore -e 'packageVersion("clusterProfiler")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}enrichplot: 	                       	${w}" ; R --no-restore -e 'packageVersion("enrichplot")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}pathview: 	                       	${w}" ; R --no-restore -e 'packageVersion("pathview")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}pacman: 	                       	${w}" ; R --no-restore -e 'packageVersion("pacman")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}tidyverse: 	                       	${w}" ; R --no-restore -e 'packageVersion("tidyverse")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ggrepel: 	                       	${w}" ; R --no-restore -e 'packageVersion("ggrepel")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}colorspace: 	                       	${w}" ; R --no-restore -e 'packageVersion("colorspace")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}RColorBrewer: 	                       	${w}" ; R --no-restore -e 'packageVersion("RColorBrewer")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}pheatmap: 	                       	${w}" ; R --no-restore -e 'packageVersion("pheatmap")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}VennDiagram: 	                       	${w}" ; R --no-restore -e 'packageVersion("VennDiagram")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}doParallel: 	                       	${w}" ; R --no-restore -e 'packageVersion("doParallel")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}foreach: 	                       	${w}" ; R --no-restore -e 'packageVersion("foreach")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}stringi: 	                       	${w}" ; R --no-restore -e 'packageVersion("stringi")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}vegan: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("vegan")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ggpubr: 	                       	${w}" ; R --no-restore -e 'packageVersion("ggpubr")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}reshape2: 	                       	${w}" ; R --no-restore -e 'packageVersion("reshape2")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}sctransform: 	                        ${w}" ; R --no-restore -e 'packageVersion("sctransform")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}hdf5r: 	                       	        ${w}" ; R --no-restore -e 'packageVersion("hdf5r")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ggridges: 	                       	${w}" ; R --no-restore -e 'packageVersion("ggridges")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ggnewscale: 	                        ${w}" ; R --no-restore -e 'packageVersion("ggnewscale")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}ggupset: 	                       	${w}" ; R --no-restore -e 'packageVersion("ggupset")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo -n "${g}Seurat: 	                       	${w}" ; R --no-restore -e 'packageVersion("Seurat")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
echo "${g}***********************************************${w}"
