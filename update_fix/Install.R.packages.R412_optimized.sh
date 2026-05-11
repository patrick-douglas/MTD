#!/usr/bin/env bash
set -Eeo pipefail

# ============================================================
# MTD - R 4.1.2 / Bioconductor 3.14 package installer
# Version: FINAL v7 - Maaslin2/Seurat dependency completion fixed
# ============================================================
#
# Assumptions:
#   1) The Conda environment was already created with:
#        conda env create -f Installation/R412.yml
#   2) You run this script from inside the MTD folder, or export:
#        dir=/path/to/MTD
#
# This script:
#   - DOES NOT recreate the Conda environment
#   - DOES NOT purge/remove the R library
#   - DOES NOT install libjpeg-turbo
#   - DOES NOT rely on remotes::install_version()
#   - uses direct CRAN Archive URLs for pinned packages
#   - patches ggtree 3.2.1 for ggplot2 compatibility
#   - patches cplm before installing Maaslin2
#   - validates cplm in a fresh R process to avoid stale namespace/cache state
#   - pins sctransform 0.3.5 before installing Seurat
#
# Usage:
#   bash update_fix/Install.R.packages.R412_optimized.sh.sh
#
# Optional variables:
#   ENV_NAME=R412
#   BIOC_VERSION=3.14
#   THREADS=20
#   RUN_GGTREE_PATCH=1
#   INSTALL_SEURAT=1
#   INSTALL_LOCAL_PATCHES=1
#   CONDA_ENSURE=0
#   CRAN_REPO=https://packagemanager.posit.co/cran/2022-05-15
#   CRAN_ARCHIVE_REPO=https://cran.r-project.org
# ============================================================

ENV_NAME="${ENV_NAME:-R412}"
BIOC_VERSION="${BIOC_VERSION:-3.14}"
THREADS="${THREADS:-$(nproc)}"
RUN_GGTREE_PATCH="${RUN_GGTREE_PATCH:-1}"
INSTALL_SEURAT="${INSTALL_SEURAT:-1}"
INSTALL_LOCAL_PATCHES="${INSTALL_LOCAL_PATCHES:-1}"
CONDA_ENSURE="${CONDA_ENSURE:-0}"
CRAN_REPO="${CRAN_REPO:-https://packagemanager.posit.co/cran/2022-05-15}"
CRAN_ARCHIVE_REPO="${CRAN_ARCHIVE_REPO:-https://cran.r-project.org}"
CONDA_SH="${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"

log() {
  printf '\n[%s] %s\n' "$(date '+%F %T')" "$*"
}

warn() {
  printf '\n[%s] WARNING: %s\n' "$(date '+%F %T')" "$*" >&2
}

die() {
  printf '\n[%s] ERROR: %s\n' "$(date '+%F %T')" "$*" >&2
  exit 1
}

# ------------------------------------------------------------
# Locate MTD directory
# ------------------------------------------------------------
SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

if [[ -n "${dir:-}" && -d "${dir:-}" ]]; then
  MTD_DIR="$(readlink -f "$dir")"
elif [[ -f "$PWD/MTD.sh" || -d "$PWD/update_fix" ]]; then
  MTD_DIR="$(readlink -f "$PWD")"
elif [[ "$(basename "$SCRIPT_DIR")" == "update_fix" ]]; then
  MTD_DIR="$(readlink -f "$(dirname "$SCRIPT_DIR")")"
else
  die "Could not detect MTD directory. Run this script from inside the MTD folder or export dir=/path/to/MTD."
fi

PATCH_DIR="${PATCH_DIR:-$MTD_DIR/update_fix/pvr_pkg}"
LOGDIR="${LOGDIR:-$MTD_DIR/update_fix/R412_post_recreate_logs}"
mkdir -p "$LOGDIR"

log "MTD directory: $MTD_DIR"
log "Patch directory: $PATCH_DIR"
log "Log directory: $LOGDIR"
log "Conda env: $ENV_NAME"
log "Bioconductor version: $BIOC_VERSION"
log "CRAN repo/snapshot: $CRAN_REPO"
log "CRAN archive repo: $CRAN_ARCHIVE_REPO"
log "Threads: $THREADS"
log "CONDA_ENSURE: $CONDA_ENSURE"

# ------------------------------------------------------------
# Conda activation
# ------------------------------------------------------------
if [[ ! -f "$CONDA_SH" ]]; then
  die "conda.sh not found: $CONDA_SH"
fi

# Do not use set -u around conda activation.
# shellcheck source=/dev/null
source "$CONDA_SH"
conda activate "$ENV_NAME"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  die "CONDA_PREFIX is empty after activating $ENV_NAME"
fi

log "Conda prefix: $CONDA_PREFIX"
log "R: $(command -v R || true)"
log "Rscript: $(command -v Rscript || true)"

R_VERSION="$(Rscript -e 'cat(as.character(getRversion()))')"
log "Detected R version: $R_VERSION"

if [[ "$R_VERSION" != 4.1.* ]]; then
  die "This installer expects R 4.1.x. Detected: $R_VERSION"
fi

# ------------------------------------------------------------
# Build environment variables
# ------------------------------------------------------------
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:$CONDA_PREFIX/share/pkgconfig:${PKG_CONFIG_PATH:-}"
export CPATH="$CONDA_PREFIX/include:$CONDA_PREFIX/include/freetype2:${CPATH:-}"
export LIBRARY_PATH="$CONDA_PREFIX/lib:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
export MAKEFLAGS="-j$THREADS"

# Help hdf5r detect Conda HDF5 when it is not installed via Conda.
if [[ -x "$CONDA_PREFIX/bin/h5cc" ]]; then
  export HDF5_LIBS="$("$CONDA_PREFIX/bin/h5cc" -show 2>/dev/null | sed 's/^[^ ]* //')"
  export HDF5_CPPFLAGS="-I$CONDA_PREFIX/include"
fi

mkdir -p "$HOME/.R"
cat > "$HOME/.R/Makevars" <<MAKEVARS
CXX11 = x86_64-conda-linux-gnu-c++
CXX14 = x86_64-conda-linux-gnu-c++
CXX17 = x86_64-conda-linux-gnu-c++

CXX11STD = -std=gnu++14
CXX14STD = -std=gnu++14
CXX17STD = -std=gnu++17

CXXFLAGS += -O2
CXX11FLAGS += -O2
CXX14FLAGS += -O2
CXX17FLAGS += -O2
MAKEVARS

export R_MAKEVARS_USER="$HOME/.R/Makevars"

# ------------------------------------------------------------
# Optional Conda ensure step
# ------------------------------------------------------------
# Default is OFF because Installation/R412.yml should already provide
# the correct foundation. This safe list intentionally does NOT include
# libjpeg-turbo.
# ------------------------------------------------------------
if [[ "$CONDA_ENSURE" == "1" ]]; then
  log "Ensuring Conda compiler/system dependencies without libjpeg-turbo"

  conda install -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    c-compiler \
    cxx-compiler \
    fortran-compiler \
    make \
    cmake \
    pkg-config \
    autoconf \
    automake \
    libtool \
    m4 \
    curl \
    libcurl \
    openssl \
    libssh2 \
    libgit2 \
    libxml2 \
    libuv \
    udunits2 \
    hdf5 \
    graphviz \
    fontconfig \
    freetype \
    harfbuzz \
    fribidi \
    cairo \
    pango \
    libpng \
    libtiff \
    jpeg \
    zlib \
    xz \
    bzip2 \
    tar \
    wget \
    rsync \
    sed \
    perl \
    r-hdf5r \
    bioconductor-rgraphviz \
    2>&1 | tee "$LOGDIR/01_conda_ensure.log"
else
  log "Skipping Conda package installation because CONDA_ENSURE=$CONDA_ENSURE"
  log "Assuming Installation/R412.yml already created the correct foundation"
fi

# Remove stale locks only. Do not purge packages.
rm -rf "$CONDA_PREFIX/lib/R/library/00LOCK"*

# ------------------------------------------------------------
# R package installation: phase 1
# ------------------------------------------------------------
R_INSTALL_SCRIPT="$LOGDIR/install_mtd_R412_packages_phase1.R"

cat > "$R_INSTALL_SCRIPT" <<'RSCRIPT'
threads <- as.integer(Sys.getenv("THREADS", "1"))
cran_repo <- Sys.getenv("CRAN_REPO", "https://packagemanager.posit.co/cran/2022-05-15")
cran_archive_repo <- Sys.getenv("CRAN_ARCHIVE_REPO", "https://cran.r-project.org")
bioc_version <- Sys.getenv("BIOC_VERSION", "3.14")
patch_dir <- path.expand(Sys.getenv("PATCH_DIR", "~/MTD/update_fix/pvr_pkg"))
install_seurat <- identical(Sys.getenv("INSTALL_SEURAT", "1"), "1")
install_local_patches <- identical(Sys.getenv("INSTALL_LOCAL_PATCHES", "1"), "1")

options(
  timeout = max(1000, getOption("timeout")),
  Ncpus = threads
)

Sys.setenv(MAKEFLAGS = paste0("-j", threads))

msg <- function(...) cat("\n==== ", sprintf(...), " ====\n", sep = "")
warn <- function(...) warning(sprintf(...), call. = FALSE)

bioc_repos <- c(
  BioCsoft = sprintf("https://bioconductor.org/packages/%s/bioc", bioc_version),
  BioCann = sprintf("https://bioconductor.org/packages/%s/data/annotation", bioc_version),
  BioCexp = sprintf("https://bioconductor.org/packages/%s/data/experiment", bioc_version),
  BioCworkflows = sprintf("https://bioconductor.org/packages/%s/workflows", bioc_version),
  CRAN = cran_repo
)

options(repos = bioc_repos)

pkg_ok <- function(pkg) requireNamespace(pkg, quietly = TRUE)
pkg_ver <- function(pkg) {
  if (pkg_ok(pkg)) as.character(utils::packageVersion(pkg)) else "MISSING"
}
pkg_version_ge <- function(pkg, ver) {
  pkg_ok(pkg) && utils::packageVersion(pkg) >= package_version(ver)
}

fresh_load_ok <- function(pkg) {
  # Check package loading in a brand-new R process.
  # This avoids false negatives caused by the current installer session having
  # stale namespaces after packages like Matrix were upgraded in-place.
  expr <- sprintf(
    "suppressPackageStartupMessages(library(%s)); cat(as.character(packageVersion('%s')))
",
    pkg,
    pkg
  )
  out <- tempfile(sprintf("%s_load_", pkg), fileext = ".log")
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = c("-e", shQuote(expr)),
    stdout = out,
    stderr = out
  )
  txt <- if (file.exists(out)) readLines(out, warn = FALSE) else character()
  if (!identical(status, 0L)) {
    cat(sprintf("\n==== Fresh R load check failed for %s ====\n", pkg))
    cat(paste(txt, collapse = "\n"), "\n")
    return(FALSE)
  }
  msg("Fresh R load check OK for %s: %s", pkg, paste(txt, collapse = " "))
  TRUE
}

install_cran <- function(pkgs, force = FALSE) {
  pkgs <- unique(pkgs)
  todo <- if (force) pkgs else pkgs[!vapply(pkgs, pkg_ok, logical(1))]
  if (!length(todo)) {
    msg("CRAN packages already installed: %s", paste(pkgs, collapse = ", "))
    return(invisible(TRUE))
  }

  options(repos = c(CRAN = cran_repo))
  msg("Installing CRAN packages: %s", paste(todo, collapse = ", "))
  install.packages(todo, dependencies = c("Depends", "Imports", "LinkingTo"))
}

install_bioc <- function(pkgs, force = FALSE) {
  pkgs <- unique(pkgs)

  if (!pkg_ok("BiocManager")) {
    options(repos = c(CRAN = cran_repo))
    install.packages("BiocManager", dependencies = c("Depends", "Imports", "LinkingTo"))
  }

  options(repos = bioc_repos)

  todo <- if (force) pkgs else pkgs[!vapply(pkgs, pkg_ok, logical(1))]
  if (!length(todo)) {
    msg("Bioconductor packages already installed: %s", paste(pkgs, collapse = ", "))
    return(invisible(TRUE))
  }

  msg("Installing Bioconductor packages: %s", paste(todo, collapse = ", "))

  BiocManager::install(
    todo,
    version = bioc_version,
    ask = FALSE,
    update = FALSE,
    force = force,
    Ncpus = threads
  )
}

install_local_tarball <- function(filename) {
  f <- file.path(patch_dir, filename)
  if (!file.exists(f)) {
    warn("Local tarball not found, skipping: %s", f)
    return(invisible(FALSE))
  }
  msg("Installing local tarball: %s", filename)
  install.packages(f, repos = NULL, type = "source", dependencies = FALSE)
  invisible(TRUE)
}

install_url_safe <- function(url, pkg = NULL, version = NULL, force = FALSE) {
  if (!is.null(pkg) && !force && pkg_ok(pkg)) {
    if (is.null(version) || identical(pkg_ver(pkg), version)) {
      msg("%s %s already installed", pkg, pkg_ver(pkg))
      return(invisible(TRUE))
    }
  }

  msg("Installing URL: %s", url)
  install.packages(url, repos = NULL, type = "source", dependencies = FALSE)
  invisible(TRUE)
}

install_archive <- function(pkg, version, force = FALSE) {
  # CRAN Archive filenames use package_version.tar.gz.
  url <- sprintf(
    "%s/src/contrib/Archive/%s/%s_%s.tar.gz",
    sub("/+$", "", cran_archive_repo),
    pkg,
    pkg,
    version
  )
  install_url_safe(url, pkg = pkg, version = version, force = force)
}

install_patched_cplm <- function(version = "0.7-10", force = TRUE) {
  # cplm 0.7-9/0.7-10 can compile against this environment but fail at load time with:
  #   undefined symbol: GET_SLOT
  # We patch src/common.h to include Rdefines.h and define GET_SLOT -> R_do_slot when needed.
  if (!force && pkg_ok("cplm")) {
    ok <- suppressWarnings(suppressPackageStartupMessages(require("cplm", quietly = TRUE, character.only = TRUE)))
    if (isTRUE(ok)) {
      msg("cplm %s already installed and loadable", pkg_ver("cplm"))
      return(invisible(TRUE))
    }
  }

  msg("Installing patched cplm %s", version)

  # Remove any broken cplm installation and stale locks.
  for (lp in .libPaths()) {
    unlink(file.path(lp, "cplm"), recursive = TRUE, force = TRUE)
  }
  unlink(Sys.glob(file.path(.libPaths(), "00LOCK*")), recursive = TRUE, force = TRUE)

  workdir <- tempfile("cplm_patch_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)

  tarball <- file.path(workdir, sprintf("cplm_%s.tar.gz", version))

  # Prefer the CRAN snapshot URL used by this installer. Fall back to CRAN Archive if needed.
  urls <- c(
    sprintf("%s/src/contrib/cplm_%s.tar.gz", sub("/+$", "", cran_repo), version),
    sprintf("%s/src/contrib/Archive/cplm/cplm_%s.tar.gz", sub("/+$", "", cran_archive_repo), version)
  )

  downloaded <- FALSE
  for (u in urls) {
    msg("Trying cplm source: %s", u)
    ok <- tryCatch({
      utils::download.file(u, tarball, mode = "wb", quiet = FALSE)
      TRUE
    }, error = function(e) FALSE)
    if (ok && file.exists(tarball) && file.info(tarball)$size > 0) {
      downloaded <- TRUE
      break
    }
  }

  if (!downloaded) {
    stop("Could not download cplm source tarball")
  }

  utils::untar(tarball, exdir = workdir)
  srcdir <- file.path(workdir, "cplm")
  common_h <- file.path(srcdir, "src", "common.h")

  if (!file.exists(common_h)) {
    stop("Could not find cplm/src/common.h after extraction")
  }

  txt <- paste(readLines(common_h, warn = FALSE), collapse = "\n")

  if (!grepl("#include <Rdefines.h>", txt, fixed = TRUE)) {
    if (grepl("#include <Rinternals.h>", txt, fixed = TRUE)) {
      txt <- sub(
        "#include <Rinternals.h>",
        "#include <Rinternals.h>\n#include <Rdefines.h>",
        txt,
        fixed = TRUE
      )
    } else {
      txt <- paste("#include <Rdefines.h>", txt, sep = "\n")
    }
  }

  fallback <- "#ifndef GET_SLOT\nextern SEXP R_do_slot(SEXP obj, SEXP name);\n#define GET_SLOT R_do_slot\n#endif"

  if (!grepl("#define GET_SLOT R_do_slot", txt, fixed = TRUE)) {
    txt <- sub(
      "#include <Rdefines.h>",
      paste("#include <Rdefines.h>", fallback, sep = "\n"),
      txt,
      fixed = TRUE
    )
  }

  writeLines(txt, common_h)

  cat("\n==== Patched cplm common.h markers ====\n")
  print(grep("Rdefines|GET_SLOT|R_do_slot", readLines(common_h, warn = FALSE), value = TRUE))

  status <- system2(
    file.path(R.home("bin"), "R"),
    args = c("CMD", "INSTALL", srcdir)
  )

  if (!identical(status, 0L)) {
    stop("Patched cplm installation failed")
  }

  # R CMD INSTALL already tests package loading in a child R process.
  # Still, validate one more time in a brand-new Rscript process. Do NOT use
  # requireNamespace()/require() in this current installer session here, because
  # Matrix was upgraded in-place earlier and the parent R session can keep stale
  # namespace/cache state, causing a false failure even after * DONE (cplm).
  if (!fresh_load_ok("cplm")) {
    stop("Patched cplm installed, but a fresh R process could not load it")
  }

  msg("Patched cplm installed and loaded successfully: %s", pkg_ver("cplm"))
  invisible(TRUE)
}

# ------------------------------------------------------------
# Bootstrap: no Suggests, no git2r/mockery/webfakes cascade
# ------------------------------------------------------------
msg("Bootstrap remotes/BiocManager")
install_cran(c("remotes", "BiocManager"))

msg("Setting Bioconductor %s", bioc_version)
BiocManager::install(version = bioc_version, ask = FALSE, update = FALSE)

# ------------------------------------------------------------
# Verify R core/recommended packages from YAML/Conda
# ------------------------------------------------------------
core_pkgs <- c("lattice", "MASS", "Matrix", "mgcv", "nlme", "survival")
msg("Checking R core/recommended packages")
print(sapply(core_pkgs, pkg_ver))
missing_core <- core_pkgs[!vapply(core_pkgs, pkg_ok, logical(1))]
if (length(missing_core)) {
  stop("Core R packages missing before MTD installation: ", paste(missing_core, collapse = ", "))
}

# ------------------------------------------------------------
# Pinned low-level packages
# ------------------------------------------------------------
if (install_local_patches && file.exists(file.path(patch_dir, "Matrix_1.6-5.tar.gz"))) {
  install_local_tarball("Matrix_1.6-5.tar.gz")
} else {
  install_archive("Matrix", "1.6-5", force = !pkg_version_ge("Matrix", "1.6-5"))
}

if (install_local_patches && file.exists(file.path(patch_dir, "MASS_7.3-60.tar.gz"))) {
  install_local_tarball("MASS_7.3-60.tar.gz")
}

# BH is header-only and old versions are often needed by older C++ packages.
install_archive("BH", "1.75.0-0")

# ------------------------------------------------------------
# CRAN foundation, including dependencies that previously broke:
# fs -> yulab.utils
# dplyr/yulab.utils -> ggfun
# lme4 >= 1.1.31 -> pbkrtest 0.5.2
# hdf5r/Rgraphviz -> Maaslin2/pathview stack
# SeuratObject/Seurat extra dependencies
# ------------------------------------------------------------
cran_foundation <- c(
  # general build/runtime
  "fs", "digest", "Rcpp", "RcppEigen", "RcppParallel", "RcppProgress",
  "BH", "cpp11", "pkgconfig", "rlang", "cli", "glue", "lifecycle",
  "vctrs", "pillar", "tibble", "magrittr", "generics", "withr",
  "tidyselect", "purrr", "stringi", "stringr", "tidyr", "dplyr",
  "dbplyr", "forcats", "readr", "readxl", "haven", "modelr",
  "lubridate", "hms", "reprex", "broom", "tidyverse",

  # plotting
  "ggplot2", "gtable", "isoband", "scales", "farver", "labeling",
  "RColorBrewer", "viridisLite", "viridis", "cowplot", "patchwork",
  "ggrepel", "ggplotify", "ggnewscale", "aplot", "scatterpie",
  "ggridges", "ggupset", "VennDiagram", "pheatmap",

  # system/interface
  "curl", "httr", "openssl", "xml2", "DBI", "RSQLite", "systemfonts",
  "textshaping", "ragg", "units", "hdf5r",

  # stats/model dependencies
  # coda/biglm/tweedie are required by patched cplm before Maaslin2.
  # Maaslin2 1.8.0 also requires robustbase, pcaPP, pbapply,
  # chemometrics, hash, pscl, and Bioconductor lpsymphony.
  "numDeriv", "coda", "biglm", "tweedie",
  "robustbase", "pcaPP", "pbapply", "chemometrics", "hash", "pscl",
  "minqa", "nloptr", "statmod", "lme4", "lmerTest",
  "TMB", "glmmTMB", "car", "carData",

  # Seurat stack
  "progressr", "future", "future.apply", "globals", "listenv",
  "parallelly", "fitdistrplus", "ica", "leiden", "lmtest", "plotly",
  "RANN", "RcppAnnoy", "reticulate", "ROCR", "scattermore",
  "spatstat.data", "spatstat.utils", "spatstat.sparse",
  "spatstat.random", "spatstat.geom", "spatstat.explore",
  "deldir", "polyclip", "FNN", "RSpectra", "irlba", "Rtsne",
  "uwot", "spam",

  # MTD/other
  "ade4", "vegan", "biomformat", "TMB", "processx", "promises",
  "miniUI", "shiny", "htmltools", "fastmap", "plyr", "prettydoc",
  "pacman", "doParallel", "foreach", "iterators",
  "R.methodsS3", "R.oo", "R.utils", "gson", "conflicted", "logging",
  "optparse", "data.table"
)

install_cran(cran_foundation)

# Force lme4 new enough for pbkrtest 0.5.2. The 2022 snapshot may provide 1.1-29.
if (!pkg_version_ge("lme4", "1.1.31")) {
  install_archive("lme4", "1.1-35.1", force = TRUE)
}

if (install_local_patches && file.exists(file.path(patch_dir, "pbkrtest_0.5.2.tar.gz"))) {
  install_local_tarball("pbkrtest_0.5.2.tar.gz")
} else {
  install_archive("pbkrtest", "0.5.2", force = !pkg_version_ge("pbkrtest", "0.5.2"))
}

# YuLab helper packages. Prefer local tarballs if present.
if (install_local_patches && file.exists(file.path(patch_dir, "yulab.utils_0.1.9.tar.gz"))) {
  install_local_tarball("yulab.utils_0.1.9.tar.gz")
} else {
  install_archive("yulab.utils", "0.1.9")
}

if (install_local_patches && file.exists(file.path(patch_dir, "ggfun_0.1.6.tar.gz"))) {
  install_local_tarball("ggfun_0.1.6.tar.gz")
} else if (install_local_patches && file.exists(file.path(patch_dir, "ggfun_0.1.7.tar.gz"))) {
  install_local_tarball("ggfun_0.1.7.tar.gz")
} else {
  install_archive("ggfun", "0.1.6")
}

# Optional local CRAN/helper patches, applied after normal CRAN install.
if (install_local_patches) {
  local_patches <- c(
    "gtable_0.3.6.tar.gz",
    "ggnewscale_0.5.2.tar.gz",
    "gson_0.1.0.tar.gz",
    "R.methodsS3_1.8.2.tar.gz",
    "R.oo_1.27.0.tar.gz",
    "R.utils_2.13.0.tar.gz",
    "lifecycle_1.0.3.tar.gz",
    "cli_3.6.2.tar.gz",
    "rlang_1.1.2.tar.gz",
    "vctrs_0.6.4.tar.gz",
    "tidyselect_1.2.1.tar.gz",
    "purrr_1.0.1.tar.gz",
    "tibble_3.2.1.tar.gz",
    "dplyr_1.1.4.tar.gz",
    "dbplyr_2.3.4.tar.gz",
    "tidyr_1.3.0.tar.gz",
    "forcats_1.0.0.tar.gz",
    "tidyverse_2.0.0.tar.gz",
    "conflicted_1.1.0.tar.gz",
    "httpuv_1.6.0.tar.gz",
    "fastmap_1.2.0.tar.gz"
  )
  for (x in local_patches) {
    if (file.exists(file.path(patch_dir, x))) install_local_tarball(x)
  }
}

# ------------------------------------------------------------
# Bioconductor packages - phase 1
# Important: ggtree, enrichplot and clusterProfiler are installed later.
# ------------------------------------------------------------
bioc_pkgs_phase1 <- c(
  "BiocGenerics", "S4Vectors", "IRanges", "XVector",
  "GenomeInfoDbData", "GenomeInfoDb", "GenomicRanges",
  "MatrixGenerics", "DelayedArray", "SummarizedExperiment",
  "Biobase", "BiocParallel", "BiocFileCache", "AnnotationDbi",
  "annotate", "AnnotationForge", "GO.db", "DO.db", "org.Hs.eg.db",
  "GOSemSim", "DOSE", "treeio", "tidytree",
  "fgsea",
  "biomformat", "multtest", "RProtoBufLib", "SingleCellExperiment",
  "phyloseq", "GenomicFeatures", "metagenomeSeq", "biomaRt",
  "graph", "Rgraphviz", "KEGGREST", "KEGGgraph", "pathview",
  "lpsymphony",
  "genefilter", "geneplotter", "limma", "edgeR", "DESeq2",
  "glmGamPoi", "MAST", "microbiome", "mia", "ANCOMBC",
  "tximeta", "cmapR", "flowCore", "cytolib"
)

install_bioc(bioc_pkgs_phase1)

# ------------------------------------------------------------
# cplm/Maaslin2 stack
# ------------------------------------------------------------
# cplm will fail before compiling if these CRAN dependencies are missing.
# Maaslin2 is installed by URL with dependencies=FALSE, so every required
# dependency must be present before installing Maaslin2 itself.
cplm_deps <- c("coda", "biglm", "tweedie")
maaslin2_cran_deps <- c("robustbase", "pcaPP", "pbapply", "chemometrics", "hash", "pscl")
maaslin2_bioc_deps <- c("lpsymphony")

install_cran(c(cplm_deps, maaslin2_cran_deps))
install_bioc(maaslin2_bioc_deps)

missing_cplm_deps <- cplm_deps[!vapply(cplm_deps, pkg_ok, logical(1))]
if (length(missing_cplm_deps)) {
  stop("cplm dependencies still missing: ", paste(missing_cplm_deps, collapse = ", "))
}

missing_maaslin2_cran_deps <- maaslin2_cran_deps[!vapply(maaslin2_cran_deps, pkg_ok, logical(1))]
missing_maaslin2_bioc_deps <- maaslin2_bioc_deps[!vapply(maaslin2_bioc_deps, pkg_ok, logical(1))]
missing_maaslin2_deps <- c(missing_maaslin2_cran_deps, missing_maaslin2_bioc_deps)
if (length(missing_maaslin2_deps)) {
  stop("Maaslin2 dependencies still missing before Maaslin2 install: ", paste(missing_maaslin2_deps, collapse = ", "))
}

install_patched_cplm("0.7-10", force = TRUE)

# Install Maaslin2 directly after patched cplm is present, so BiocManager does not try to reinstall cplm.
install_url_safe(
  sprintf("https://bioconductor.org/packages/%s/bioc/src/contrib/Maaslin2_1.8.0.tar.gz", bioc_version),
  pkg = "Maaslin2",
  force = TRUE
)

if (!fresh_load_ok("Maaslin2")) {
  stop("Maaslin2 installed, but a fresh R process could not load it")
}

# ------------------------------------------------------------
# Seurat stack
# ------------------------------------------------------------
if (install_seurat) {
  # Seurat 4.3.0 requires sctransform >= 0.3.5, and also directly needs pbapply.
  # pbapply is also a Maaslin2 dependency, but keep this guard here so Seurat
  # remains robust even if Maaslin2 installation is disabled later.
  install_cran("pbapply")
  if (!pkg_ok("pbapply")) stop("Seurat dependency still missing: pbapply")

  install_archive("sctransform", "0.3.5", force = !pkg_version_ge("sctransform", "0.3.5"))
  install_archive("SeuratObject", "4.1.3", force = !pkg_version_ge("SeuratObject", "4.1.3"))
  install_archive("Seurat", "4.3.0", force = !pkg_version_ge("Seurat", "4.3.0"))

  if (!fresh_load_ok("Seurat")) {
    stop("Seurat installed, but a fresh R process could not load it")
  }
}

important_versions <- c(
  "Rcpp", "RcppArmadillo", "Matrix", "MASS", "lattice",
  "ggplot2", "ggfun", "yulab.utils", "tidytree", "treeio",
  "DOSE", "fgsea", "DESeq2", "phyloseq", "microbiome", "mia",
  "ANCOMBC", "cplm", "Maaslin2", "Rgraphviz", "hdf5r",
  "cytolib", "flowCore", "cmapR",
  "lme4", "pbkrtest", "sctransform", "SeuratObject", "Seurat"
)

msg("Important package versions before ggtree patch")
print(sapply(important_versions, pkg_ver))
RSCRIPT

log "Running R package installation phase 1"
THREADS="$THREADS" \
CRAN_REPO="$CRAN_REPO" \
CRAN_ARCHIVE_REPO="$CRAN_ARCHIVE_REPO" \
BIOC_VERSION="$BIOC_VERSION" \
PATCH_DIR="$PATCH_DIR" \
INSTALL_SEURAT="$INSTALL_SEURAT" \
INSTALL_LOCAL_PATCHES="$INSTALL_LOCAL_PATCHES" \
Rscript "$R_INSTALL_SCRIPT" 2>&1 | tee "$LOGDIR/03_R_install_phase1.log"

# ------------------------------------------------------------
# ggtree compatibility patch
# ------------------------------------------------------------
patch_ggtree() {
  log "Applying ggtree 3.2.1 compatibility patch"

  local workdir
  workdir="$(mktemp -d /tmp/ggtree_patch_XXXXXX)"
  local tarball="$workdir/ggtree_3.2.1.tar.gz"

  curl -L --retry 3 --fail \
    "https://bioconductor.org/packages/3.14/bioc/src/contrib/ggtree_3.2.1.tar.gz" \
    -o "$tarball"

  tar -xzf "$tarball" -C "$workdir"

  # Avoid importing ggplot2 internal warning_wrap from NAMESPACE.
  sed -i '/warning_wrap/d' "$workdir/ggtree/NAMESPACE" || true

  cat > "$workdir/ggtree/R/ggplot2_compat_warning_wrap.R" <<'EOF2'
warning_wrap <- function(..., call. = FALSE, immediate. = FALSE) {
  msg <- paste0(..., collapse = "")
  warning(msg, call. = call., immediate. = immediate.)
}
EOF2

  # geom_hilight.R used some ggplot2 internals that changed across ggplot2 versions.
  if [[ -f "$workdir/ggtree/R/geom_hilight.R" ]]; then
    sed -i.bak '/warning_wrap <- getFromNamespace("warning_wrap", "ggplot2")/d' "$workdir/ggtree/R/geom_hilight.R" || true
    sed -i.bak '/rect_to_poly <- getFromNamespace("rect_to_poly", "ggplot2")/d' "$workdir/ggtree/R/geom_hilight.R" || true
    sed -i.bak '/new_data_frame <- getFromNamespace("new_data_frame", "ggplot2")/d' "$workdir/ggtree/R/geom_hilight.R" || true
  fi

  cat > "$workdir/ggtree/R/ggplot2_compat_rect_to_poly.R" <<'EOF2'
rect_to_poly <- function(xmin, xmax, ymin, ymax) {
  data.frame(
    x = as.vector(rbind(xmin, xmax, xmax, xmin, NA)),
    y = as.vector(rbind(ymax, ymax, ymin, ymin, NA))
  )
}
EOF2

  cat > "$workdir/ggtree/R/ggplot2_compat_new_data_frame.R" <<'EOF2'
new_data_frame <- function(x = list(), n = NULL) {
  if (!is.list(x)) {
    stop("x must be a list", call. = FALSE)
  }

  if (is.null(n)) {
    n <- if (length(x) == 0) 0L else max(lengths(x))
  }

  for (i in seq_along(x)) {
    if (length(x[[i]]) == 1L && n > 1L) {
      x[[i]] <- rep(x[[i]], n)
    }
  }

  class(x) <- "data.frame"
  attr(x, "row.names") <- .set_row_names(n)
  x
}
EOF2

  find "$workdir/ggtree/R" -type f -name "*.bak" -delete

  rm -rf "$CONDA_PREFIX/lib/R/library/00LOCK"*
  R CMD INSTALL "$workdir/ggtree" 2>&1 | tee "$LOGDIR/04_ggtree_patch.log"
  rm -rf "$CONDA_PREFIX/lib/R/library/00LOCK"*

  rm -rf "$workdir"
}

if [[ "$RUN_GGTREE_PATCH" == "1" ]]; then
  patch_ggtree
else
  log "Skipping ggtree patch because RUN_GGTREE_PATCH=$RUN_GGTREE_PATCH"
fi

# ------------------------------------------------------------
# R package installation: phase 2
# Install enrichplot/clusterProfiler after ggtree is patched.
# ------------------------------------------------------------
R_PHASE2_SCRIPT="$LOGDIR/install_mtd_R412_packages_phase2.R"

cat > "$R_PHASE2_SCRIPT" <<'RSCRIPT'
threads <- as.integer(Sys.getenv("THREADS", "1"))
cran_repo <- Sys.getenv("CRAN_REPO", "https://packagemanager.posit.co/cran/2022-05-15")
bioc_version <- Sys.getenv("BIOC_VERSION", "3.14")

options(
  repos = c(
    BioCsoft = sprintf("https://bioconductor.org/packages/%s/bioc", bioc_version),
    BioCann = sprintf("https://bioconductor.org/packages/%s/data/annotation", bioc_version),
    BioCexp = sprintf("https://bioconductor.org/packages/%s/data/experiment", bioc_version),
    BioCworkflows = sprintf("https://bioconductor.org/packages/%s/workflows", bioc_version),
    CRAN = cran_repo
  ),
  timeout = max(1000, getOption("timeout")),
  Ncpus = threads
)

Sys.setenv(MAKEFLAGS = paste0("-j", threads))

msg <- function(...) cat("\n==== ", sprintf(...), " ====\n", sep = "")
pkg_ok <- function(pkg) requireNamespace(pkg, quietly = TRUE)
pkg_ver <- function(pkg) if (pkg_ok(pkg)) as.character(packageVersion(pkg)) else "MISSING"

if (!pkg_ok("BiocManager")) {
  install.packages("BiocManager", dependencies = c("Depends", "Imports", "LinkingTo"))
}

msg("Checking ggtree before installing enrichplot/clusterProfiler")
print(sapply(c("ggplot2", "treeio", "tidytree", "ggtree", "DOSE"), pkg_ver))

if (!pkg_ok("ggtree")) {
  stop("ggtree is missing before phase 2. The ggtree patch probably failed.")
}

msg("Installing enrichplot and clusterProfiler")
BiocManager::install(
  c("enrichplot", "clusterProfiler"),
  version = bioc_version,
  ask = FALSE,
  update = FALSE,
  force = TRUE,
  Ncpus = threads
)

msg("Phase 2 versions")
print(sapply(c("ggplot2", "ggtree", "DOSE", "enrichplot", "clusterProfiler"), pkg_ver))
RSCRIPT

log "Running R package installation phase 2"
THREADS="$THREADS" \
CRAN_REPO="$CRAN_REPO" \
BIOC_VERSION="$BIOC_VERSION" \
Rscript "$R_PHASE2_SCRIPT" 2>&1 | tee "$LOGDIR/06_R_install_phase2_enrich_clusterProfiler.log"

# ------------------------------------------------------------
# Final validation
# ------------------------------------------------------------
VALIDATION_SCRIPT="$LOGDIR/validate_R412_MTD.R"

cat > "$VALIDATION_SCRIPT" <<'RSCRIPT'
options(
  repos = c(
    BioCsoft = "https://bioconductor.org/packages/3.14/bioc",
    BioCann = "https://bioconductor.org/packages/3.14/data/annotation",
    BioCexp = "https://bioconductor.org/packages/3.14/data/experiment",
    BioCworkflows = "https://bioconductor.org/packages/3.14/workflows",
    CRAN = Sys.getenv("CRAN_REPO", "https://packagemanager.posit.co/cran/2022-05-15")
  )
)

pkgs <- c(
  "lattice", "MASS", "Matrix", "mgcv", "nlme", "survival",
  "ggplot2", "ggfun", "yulab.utils", "tidytree", "treeio", "ggtree",
  "DOSE", "enrichplot", "clusterProfiler", "fgsea",
  "DESeq2", "phyloseq", "microbiome", "mia", "ANCOMBC",
  "cplm", "Maaslin2", "flowCore", "cytolib", "cmapR", "tximeta",
  "hdf5r", "Rgraphviz", "pathview",
  "lme4", "pbkrtest",
  "sctransform", "SeuratObject", "Seurat"
)

versions <- sapply(pkgs, function(p) {
  if (requireNamespace(p, quietly = TRUE)) {
    as.character(packageVersion(p))
  } else {
    "MISSING"
  }
})

cat("\n==== Package versions ====\n")
print(versions)

load_status <- sapply(pkgs, function(p) {
  suppressPackageStartupMessages(
    require(p, character.only = TRUE, quietly = TRUE)
  )
})

cat("\n==== Load status ====\n")
print(load_status)

failed <- names(load_status)[!load_status]
if (length(failed)) {
  stop("Some packages failed to load: ", paste(failed, collapse = ", "))
}

cat("\n==== BiocManager::valid() ====\n")
print(BiocManager::valid())

cat("\nValidation completed successfully.\n")
RSCRIPT

log "Running final validation"
CRAN_REPO="$CRAN_REPO" Rscript "$VALIDATION_SCRIPT" 2>&1 | tee "$LOGDIR/07_validation.log"

log "Done. Logs are in: $LOGDIR"

