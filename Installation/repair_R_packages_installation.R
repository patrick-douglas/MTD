#!/usr/bin/env Rscript

# ============================================================
# Repair missing R packages for MTD
# R 4.1.2 + Bioconductor 3.14
# ============================================================

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  timeout = 600,
  download.file.method = "libcurl"
)

# Clean possible broken locks
lib <- .libPaths()[1]
locks <- list.files(lib, pattern = "^00LOCK", full.names = TRUE)
if (length(locks) > 0) {
  message("Removing R lock folders:")
  print(locks)
  unlink(locks, recursive = TRUE, force = TRUE)
}

# ------------------------------------------------------------
# BiocManager
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  version = "3.14",
  ask = FALSE,
  update = FALSE
)

options(repos = BiocManager::repositories(version = "3.14"))

# ------------------------------------------------------------
# Critical CRAN/base packages first
# ------------------------------------------------------------

critical_cran <- c(
  "MASS",
  "remotes",
  "pacman",
  "systemfonts",
  "gdtools",
  "htmlwidgets",
  "httpuv",
  "shiny",
  "DT",
  "ggforce",
  "tidytree",
  "treeio",
  "ggiraph",
  "ggraph",
  "ggplotify",
  "shadowtext",
  "aplot",
  "ggfun",
  "yulab.utils"
)

for (pkg in critical_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing critical CRAN package: ", pkg)
    tryCatch(
      {
        install.packages(pkg, dependencies = c("Depends", "Imports", "LinkingTo"))
      },
      error = function(e) {
        message("FAILED: ", pkg, " -> ", e$message)
      }
    )
  } else {
    message(pkg, " already installed.")
  }
}

# ------------------------------------------------------------
# Missing packages from your check
# ------------------------------------------------------------

bioc_missing <- c(
  "tximeta",
  "phyloseq",
  "glmGamPoi",
  "cmapR",
  "microbiome",
  "ANCOMBC",
  "Maaslin2",
  "clusterProfiler",
  "enrichplot"
)

cran_missing <- c(
  "tidyverse",
  "readr",
  "forcats",
  "colorspace",
  "VennDiagram",
  "doParallel",
  "vegan",
  "ggpubr",
  "sctransform",
  "hdf5r",
  "ggridges",
  "ggnewscale",
  "ggupset",
  "Seurat"
)

message("------------------------------------------------------------")
message("Installing missing Bioconductor packages...")
message("------------------------------------------------------------")

tryCatch(
  {
    BiocManager::install(
      bioc_missing,
      ask = FALSE,
      update = FALSE
    )
  },
  error = function(e) {
    message("Bioconductor install failed: ", e$message)
  }
)

message("------------------------------------------------------------")
message("Installing missing CRAN packages...")
message("------------------------------------------------------------")

for (pkg in cran_missing) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing CRAN package: ", pkg)
    tryCatch(
      {
        install.packages(pkg, dependencies = c("Depends", "Imports", "LinkingTo"))
      },
      error = function(e) {
        message("FAILED: ", pkg, " -> ", e$message)
      }
    )
  } else {
    message(pkg, " already installed.")
  }
}

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

all_check <- unique(c(critical_cran, bioc_missing, cran_missing))

status <- sapply(all_check, requireNamespace, quietly = TRUE)

message("------------------------------------------------------------")
message("Repair summary")
message("------------------------------------------------------------")
message("Installed/loadable: ", sum(status))
message("Still missing: ", sum(!status))

if (any(!status)) {
  message("Still missing packages:")
  message(paste(names(status)[!status], collapse = ", "))
} else {
  message("All repair packages are now available.")
}

message("------------------------------------------------------------")
message("Repair script finished.")
