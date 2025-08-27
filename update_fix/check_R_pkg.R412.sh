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
pkgs=(lattice httr BiocManager SummarizedExperiment MASS ggplot2 ade4 biomformat igraph multtest vegan RProtoBufLib remotes cytolib flowCore DBI AnnotationDbi SingleCellExperiment phyloseq curl annotate BiocFileCache yulab.utils ggfun gtable aplot lme4 TMB GenomicFeatures metagenomeSeq biomaRt ggtree fgsea R.methodsS3 R.oo R.utils GO.db GOSemSim DOSE ggtangle scatterpie ggnewscale tidytree processx clusterProfiler gson lifecycle tidyselect rlang vctrs cli purrr tibble cmapR dplyr dbplyr tidyr ggpubr pbkrtest car Maaslin2 promises httpuv miniUI shiny htmltools fastmap rgeos SeuratObject Seurat Matrix tximeta plyr BiocGenerics S4Vectors IRanges UCSC.utils GenomeInfoDbData GenomeInfoDb matrixStats lambda.r futile.options futile.logger RColorBrewer DO.db forcats tidyverse conflicted org.Hs.eg.db graph Rgraphviz KEGGgraph pathview prettydoc AnnotationForge DESeq2 glmGamPoi MAST microbiome ANCOMBC pacman ggrepel colorspace pheatmap VennDiagram doParallel foreach stringi reshape2 hdf5r ggridges ggupset ragg)

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

