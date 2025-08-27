#!/bin/bash
# Check R package versions in a nice table for the R412 environment

# Cores
w=$(tput sgr0)
g=$(tput setaf 2)
r=$(tput setaf 1)

# Ativa o ambiente Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate R412

# Versão do R
R_ver=$(R --version | grep version | grep R | awk '{print $3}')

# Lista de pacotes a checar
pkgs=(yulab.utils DBI lattice MASS Matrix BiocManager httr SummarizedExperiment ggplot2 remotes ade4 biomformat igraph multtest 
      RProtoBufLib cytolib flowCore AnnotationDbi ggfun SingleCellExperiment phyloseq curl annotate 
      BiocFileCache aplot lme4 TMB GenomicFeatures metagenomeSeq ggtree ragg sctransform cmapR 
      httpuv miniUI shiny biomaRt tximeta limma glmGamPoi MAST microbiome ANCOMBC Maaslin2 DO.db 
      clusterProfiler enrichplot pathview pacman tidyverse ggrepel colorspace RColorBrewer pheatmap 
      VennDiagram doParallel foreach stringi vegan ggpubr reshape2 hdf5r ggridges ggnewscale ggupset Seurat 
      genefilter DESeq2)

# Cabeçalho da tabela
echo "${g}"
echo "╔══════════════════════╦═══════════════╗"
echo "║ R                    ║ ${w}$R_ver${g}         ║"
echo "║ Conda Environment    ║ ${w}R412${g}          ║"
echo "╠══════════════════════╬═══════════════╣"

# Loop pelos pacotes
for pkg in "${pkgs[@]}"; do
    ver=$(R --no-restore -e "packageVersion(\"${pkg}\")" 2>/dev/null | \
          grep '\[1\]' | awk '{print $2}' | sed -r 's/^.{1}//; s/.$//')

    if [ -z "$ver" ]; then
        ver="${r}not installed${g}"
    fi

    printf "║ %-20s ║ ${w}%-13s${g} ║\n" "$pkg" "$ver"
done

# Rodapé
echo "╚══════════════════════╩═══════════════╝"

