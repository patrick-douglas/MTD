R --no-restore -e 'packageVersion("httr")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("SummarizedExperiment")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("Matrix")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("ggplot2")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("httpuv")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("ade4")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("biomformat")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("igraph")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("multtest")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("vegan")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("RProtoBufLib")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("cytolib")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("flowCore")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("AnnotationDbi")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("igraph")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("ggfun")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("SingleCellExperiment")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("phyloseq")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("curl")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("annotate")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("BiocFileCache")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("MASS")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("aplot")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("lme4")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("TMB")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("GenomicFeatures")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("metagenomeSeq")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("ggtree")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("enrichplot")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("ragg")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("sctransform")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("biomaRt")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("DESeq2")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("tximeta")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("limma")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("phyloseq")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("glmGamPoi")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("cmapR")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("MAST")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("microbiome")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("ANCOMBC")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("Maaslin2")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("DO.db")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("clusterProfiler")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("enrichplot")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("pathview")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("pacman")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("tidyverse")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("ggrepel")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("colorspace")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("RColorBrewer")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("pheatmap")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("VennDiagram")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("doParallel")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("foreach")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("stringi")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("vegan")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("ggpubr")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("reshape2")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("sctransform")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("hdf5r")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("ggridges")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("ggnewscale")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("ggupset")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
R --no-restore -e 'packageVersion("Seurat")' | grep packageVersion -A 1 | grep '[1]' | awk {'print $2'} | sed -r 's/^.{1}//' | sed 's/.$//'
