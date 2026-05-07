#!/usr/bin/env Rscript

# ============================================================
# Check R package installation for MTD
# This script does NOT install anything.
# It only verifies if packages are available and loadable.
# ============================================================

# -----------------------------
# Expected Bioconductor version
# -----------------------------

expected_bioc_version <- "3.14"

# -----------------------------
# Packages to check
# -----------------------------

bioc_packages <- c(
  "BiocManager",
  "BiocVersion",
  "BiocGenerics",
  "genefilter",
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

cran_packages <- c(
  "remotes",
  "pacman",
  "tidyverse",
  "ggplot2",
  "tidyr",
  "dplyr",
  "stringr",
  "tibble",
  "readr",
  "purrr",
  "forcats",
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

all_packages <- unique(c(bioc_packages, cran_packages))

# -----------------------------
# Output files
# -----------------------------

out_dir <- "R_package_check_report"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

report_file  <- file.path(out_dir, "R_package_check_report.tsv")
missing_file <- file.path(out_dir, "missing_packages.txt")
loaded_file  <- file.path(out_dir, "successfully_loaded_packages.txt")
summary_file <- file.path(out_dir, "summary.txt")

# -----------------------------
# Helper function
# -----------------------------

check_package <- function(pkg) {
  namespace_available <- requireNamespace(pkg, quietly = TRUE)

  version <- NA_character_
  load_status <- FALSE
  error_message <- NA_character_

  if (namespace_available) {
    version <- tryCatch(
      as.character(packageVersion(pkg)),
      error = function(e) NA_character_
    )

    load_status <- tryCatch(
      {
        suppressPackageStartupMessages(
          library(pkg, character.only = TRUE)
        )
        TRUE
      },
      error = function(e) {
        error_message <<- e$message
        FALSE
      },
      warning = function(w) {
        TRUE
      }
    )
  } else {
    error_message <- "Package namespace not available"
  }

  data.frame(
    package = pkg,
    installed = namespace_available,
    loadable = load_status,
    version = version,
    error = error_message,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Start check
# -----------------------------

message("------------------------------------------------------------")
message("Checking R packages for MTD")
message("------------------------------------------------------------")
message("R version: ", getRversion())
message("R executable: ", R.home())
message("Library paths:")
print(.libPaths())
message("------------------------------------------------------------")

# -----------------------------
# Check Bioconductor version
# -----------------------------

biocmanager_available <- requireNamespace("BiocManager", quietly = TRUE)

current_bioc_version <- NA_character_

if (biocmanager_available) {
  current_bioc_version <- tryCatch(
    as.character(BiocManager::version()),
    error = function(e) NA_character_
  )
}

message("BiocManager available: ", biocmanager_available)
message("Detected Bioconductor version: ", current_bioc_version)
message("Expected Bioconductor version: ", expected_bioc_version)

if (!is.na(current_bioc_version) && current_bioc_version != expected_bioc_version) {
  warning(
    "Bioconductor version mismatch. Expected ",
    expected_bioc_version,
    " but detected ",
    current_bioc_version
  )
}

message("------------------------------------------------------------")

# -----------------------------
# Run package checks
# -----------------------------

results <- do.call(
  rbind,
  lapply(all_packages, check_package)
)

# -----------------------------
# Save reports
# -----------------------------

write.table(
  results,
  file = report_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

missing_packages <- results$package[!results$installed]
not_loadable_packages <- results$package[results$installed & !results$loadable]
success_packages <- results$package[results$installed & results$loadable]

writeLines(missing_packages, missing_file)
writeLines(success_packages, loaded_file)

summary_lines <- c(
  "MTD R package installation check",
  "================================",
  "",
  paste("R version:", getRversion()),
  paste("R home:", R.home()),
  paste("Expected Bioconductor version:", expected_bioc_version),
  paste("Detected Bioconductor version:", current_bioc_version),
  "",
  paste("Total packages checked:", length(all_packages)),
  paste("Installed packages:", sum(results$installed)),
  paste("Missing packages:", length(missing_packages)),
  paste("Installed but not loadable:", length(not_loadable_packages)),
  paste("Successfully loadable:", length(success_packages)),
  "",
  "Missing packages:",
  if (length(missing_packages) > 0) paste("-", missing_packages) else "None",
  "",
  "Installed but not loadable:",
  if (length(not_loadable_packages) > 0) paste("-", not_loadable_packages) else "None",
  "",
  paste("Full report:", report_file),
  paste("Missing package list:", missing_file),
  paste("Successfully loaded package list:", loaded_file)
)

writeLines(summary_lines, summary_file)

# -----------------------------
# Print summary
# -----------------------------

message("Package check finished.")
message("------------------------------------------------------------")
message("Total packages checked: ", length(all_packages))
message("Installed packages: ", sum(results$installed))
message("Missing packages: ", length(missing_packages))
message("Installed but not loadable: ", length(not_loadable_packages))
message("Successfully loadable: ", length(success_packages))
message("------------------------------------------------------------")

if (length(missing_packages) > 0) {
  message("Missing packages:")
  message(paste(missing_packages, collapse = ", "))
  message("------------------------------------------------------------")
}

if (length(not_loadable_packages) > 0) {
  message("Installed but not loadable:")
  message(paste(not_loadable_packages, collapse = ", "))
  message("------------------------------------------------------------")
}

message("Reports saved in: ", normalizePath(out_dir))
message("Main report: ", normalizePath(report_file))
message("Summary: ", normalizePath(summary_file))

# -----------------------------
# Exit status
# -----------------------------
# Exit code 0 = everything OK
# Exit code 1 = some packages missing or not loadable
# -----------------------------

if (length(missing_packages) > 0 || length(not_loadable_packages) > 0) {
  quit(status = 1)
} else {
  quit(status = 0)
}
