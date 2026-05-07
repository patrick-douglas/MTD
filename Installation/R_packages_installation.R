#!/usr/bin/env Rscript

# ============================================================
# R package installation script for MTD
# Compatible with:
#   R 4.1.2
#   Bioconductor 3.14
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ============================================================
# Local offline packages
# ============================================================

local_biocmanager <- "~/MTD/update_fix/pvr_pkg/BiocManager_1.30.22.tar.gz"
local_pacman      <- "~/MTD/update_fix/pvr_pkg/pacman_0.5.1.tar.gz"
local_remotes     <- "~/MTD/update_fix/pvr_pkg/remotes_2.4.2.tar.gz"

# ============================================================
# Helper functions
# ============================================================

install_from_local_or_cran <- function(pkg, local_file = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(pkg, " is already installed.")
    return(invisible(TRUE))
  }

  if (!is.null(local_file) && file.exists(path.expand(local_file))) {
    message("Installing ", pkg, " from local offline package: ", local_file)
    result <- tryCatch(
      {
        install.packages(
          path.expand(local_file),
          repos = NULL,
          type = "source"
        )
        TRUE
      },
      error = function(e) {
        message("ERROR installing ", pkg, " from local file: ", e$message)
        FALSE
      }
    )
  } else {
    message("Installing ", pkg, " from CRAN...")
    result <- tryCatch(
      {
        install.packages(pkg)
        TRUE
      },
      error = function(e) {
        message("ERROR installing ", pkg, " from CRAN: ", e$message)
        FALSE
      }
    )
  }

  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("WARNING: Failed to install package: ", pkg)
    return(invisible(FALSE))
  }

  invisible(TRUE)
}

safe_bioc_install <- function(pkgs, force = FALSE) {
  message("Installing Bioconductor packages: ", paste(pkgs, collapse = ", "))

  result <- tryCatch(
    {
      BiocManager::install(
        pkgs,
        ask = FALSE,
        update = FALSE,
        force = force
      )
      TRUE
    },
    error = function(e) {
      message("ERROR during Bioconductor installation: ", e$message)
      FALSE
    }
  )

  installed_status <- sapply(pkgs, requireNamespace, quietly = TRUE)
  missing <- names(installed_status)[!installed_status]

  if (length(missing) > 0) {
    message("WARNING: These Bioconductor packages are still missing:")
    message(paste(missing, collapse = ", "))
  }

  invisible(result)
}

install_cran_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing CRAN package: ", pkg)
      tryCatch(
        {
          install.packages(pkg)
        },
        error = function(e) {
          message("ERROR installing CRAN package ", pkg, ": ", e$message)
        }
      )
    } else {
      message(pkg, " is already installed.")
    }
  }
}

# ============================================================
# Check R version
# ============================================================

r_version <- getRversion()
message("Detected R version: ", r_version)

if (r_version < "4.1.0" || r_version >= "4.2.0") {
  warning(
    "This script is configured for R 4.1.x with Bioconductor 3.14. ",
    "Your R version is ", r_version, "."
  )
}

# ============================================================
# Install/load BiocManager
# ============================================================

install_from_local_or_cran(
  pkg = "BiocManager",
  local_file = local_biocmanager
)

# ============================================================
# Set Bioconductor repositories/version
# ============================================================

message("Setting Bioconductor version to 3.14...")

tryCatch(
  {
    BiocManager::install(
      version = "3.14",
      ask = FALSE,
      update = FALSE
    )
  },
  error = function(e) {
    message("WARNING: Could not fully set Bioconductor version: ", e$message)
  }
)

message("Current Bioconductor version: ", BiocManager::version())

# Important: restore Bioconductor repositories plus CRAN
options(repos = BiocManager::repositories(version = "3.14"))

message("Using repositories:")
print(getOption("repos"))

# ============================================================
# Test Bioconductor access
# ============================================================

message("Testing Bioconductor repository access...")

bioc_repos <- BiocManager::repositories(version = "3.14")
bioc_index_url <- paste0(bioc_repos["BioCsoft"], "/src/contrib/PACKAGES")

message("Testing URL: ", bioc_index_url)

bioc_access <- tryCatch(
  {
    con <- url(bioc_index_url)
    open(con)
    close(con)
    TRUE
  },
  error = function(e) {
    message("WARNING: Cannot access Bioconductor repository index.")
    message("Reason: ", e$message)
    FALSE
  }
)

if (!bioc_access) {
  message("------------------------------------------------------------")
  message("Bioconductor repository is not accessible.")
  message("This is a network/DNS/proxy/firewall problem, not an R version problem.")
  message("Bioconductor packages may fail unless you provide local/offline package files.")
  message("------------------------------------------------------------")
}

# ============================================================
# Fixes for DESeq2-related installation issues
# ============================================================

safe_bioc_install(
  c("BiocGenerics", "genefilter"),
  force = TRUE
)

# ============================================================
# Bioconductor packages
# ============================================================

bioc_packages <- c(
  "biomaRt",
  "DESeq2",
  "tximeta",
  "limma",
  "phyloseq",
  "glmGamPoi",
  "cmapR",
  "MAST",
  "microbiome",
  "ANCOMBC",
  "Maaslin2",
  "DO.db",
  "clusterProfiler",
  "enrichplot",
  "pathview"
)

safe_bioc_install(
  bioc_packages,
  force = FALSE
)

# ============================================================
# Install remotes before pacman
# pacman depends on remotes
# ============================================================

install_from_local_or_cran(
  pkg = "remotes",
  local_file = local_remotes
)

# ============================================================
# Install/load pacman
# ============================================================

install_from_local_or_cran(
  pkg = "pacman",
  local_file = local_pacman
)

# ============================================================
# CRAN packages
# ============================================================

cran_packages <- c(
  "tidyverse",
  "ggrepel",
  "colorspace",
  "RColorBrewer",
  "pheatmap",
  "VennDiagram",
  "doParallel",
  "foreach",
  "stringi",
  "vegan",
  "ggpubr",
  "reshape2",
  "sctransform",
  "hdf5r",
  "ggridges",
  "ggnewscale",
  "ggupset",
  "Seurat"
)

if (requireNamespace("pacman", quietly = TRUE)) {
  message("Installing/loading CRAN packages with pacman...")

  pacman::p_load(
    char = cran_packages,
    install = TRUE,
    update = FALSE
  )
} else {
  message("pacman is not available. Installing CRAN packages using install.packages()...")
  install_cran_packages(cran_packages)
}

# ============================================================
# Final check
# ============================================================

all_packages <- c(
  "BiocManager",
  "remotes",
  "pacman",
  bioc_packages,
  cran_packages
)

installed_status <- sapply(
  all_packages,
  requireNamespace,
  quietly = TRUE
)

missing_packages <- names(installed_status)[!installed_status]

message("------------------------------------------------------------")

if (length(missing_packages) == 0) {
  message("All requested packages appear to be installed successfully.")
} else {
  message("Some packages were not installed successfully:")
  message(paste(missing_packages, collapse = ", "))
}

message("------------------------------------------------------------")
message("R version: ", getRversion())
message("Bioconductor version: ", BiocManager::version())
message("Installation script finished.")
