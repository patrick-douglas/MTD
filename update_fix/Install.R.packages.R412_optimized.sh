#!/usr/bin/env bash
set -Eeo pipefail

# ============================================================
# MTD - R 4.1.2 / Bioconductor 3.14 package installer
# Version: FINAL v11 - Bioconductor failover and locale-safe TMB ABI validation
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
#   - probes official Bioconductor mirrors before package installation
#   - retries failed Bioconductor installs across healthy mirrors
#   - caches directly downloaded Bioconductor source tarballs
#   - supports BIOC_SMOKE_TEST_ONLY=1 for a focused recovery test
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
#   BIOC_MIRRORS=https://bioconductor.org,https://bioconductor.posit.co,...
#   BIOC_RETRY_ATTEMPTS=3
#   BIOC_RETRY_SLEEP=10
#   BIOC_SMOKE_TEST_ONLY=0
#   BIOC_CACHE_DIR=/path/to/persistent/cache
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
BIOC_MIRRORS="${BIOC_MIRRORS:-https://bioconductor.org,https://bioconductor.posit.co,https://bioconductor.statistik.tu-dortmund.de,https://ftp.gwdg.de/pub/misc/bioconductor}"
BIOC_RETRY_ATTEMPTS="${BIOC_RETRY_ATTEMPTS:-3}"
BIOC_RETRY_SLEEP="${BIOC_RETRY_SLEEP:-10}"
BIOC_SMOKE_TEST_ONLY="${BIOC_SMOKE_TEST_ONLY:-0}"
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
BIOC_CACHE_DIR="${BIOC_CACHE_DIR:-$MTD_DIR/update_fix/bioc_tarball_cache/$BIOC_VERSION}"
BIOC_ACTIVE_MIRROR_FILE="${BIOC_ACTIVE_MIRROR_FILE:-$LOGDIR/bioc_active_mirror.txt}"
mkdir -p "$LOGDIR" "$BIOC_CACHE_DIR"

log "MTD directory: $MTD_DIR"
log "Patch directory: $PATCH_DIR"
log "Log directory: $LOGDIR"
log "Conda env: $ENV_NAME"
log "Bioconductor version: $BIOC_VERSION"
log "CRAN repo/snapshot: $CRAN_REPO"
log "CRAN archive repo: $CRAN_ARCHIVE_REPO"
log "Bioconductor mirrors: $BIOC_MIRRORS"
log "Bioconductor retry attempts: $BIOC_RETRY_ATTEMPTS"
log "Bioconductor retry base sleep: $BIOC_RETRY_SLEEP seconds"
log "Bioconductor tarball cache: $BIOC_CACHE_DIR"
log "Bioconductor smoke-test only: $BIOC_SMOKE_TEST_ONLY"
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
bioc_smoke_test_only <- identical(Sys.getenv("BIOC_SMOKE_TEST_ONLY", "0"), "1")
bioc_cache_dir <- path.expand(Sys.getenv("BIOC_CACHE_DIR", file.path(tempdir(), "mtd_bioc_cache")))
bioc_active_mirror_file <- path.expand(Sys.getenv("BIOC_ACTIVE_MIRROR_FILE", file.path(tempdir(), "mtd_bioc_active_mirror.txt")))
bioc_retry_attempts <- suppressWarnings(as.integer(Sys.getenv("BIOC_RETRY_ATTEMPTS", "3")))
bioc_retry_sleep <- suppressWarnings(as.numeric(Sys.getenv("BIOC_RETRY_SLEEP", "10")))
if (is.na(bioc_retry_attempts) || bioc_retry_attempts < 1L) bioc_retry_attempts <- 3L
if (is.na(bioc_retry_sleep) || bioc_retry_sleep < 0) bioc_retry_sleep <- 10
bioc_mirrors <- trimws(strsplit(
  Sys.getenv(
    "BIOC_MIRRORS",
    paste(c(
      "https://bioconductor.org",
      "https://bioconductor.posit.co",
      "https://bioconductor.statistik.tu-dortmund.de",
      "https://ftp.gwdg.de/pub/misc/bioconductor"
    ), collapse = ",")
  ),
  ",", fixed = TRUE
)[[1]])
bioc_mirrors <- unique(sub("/+$", "", bioc_mirrors[nzchar(bioc_mirrors)]))
dir.create(bioc_cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(bioc_active_mirror_file), recursive = TRUE, showWarnings = FALSE)

options(
  timeout = max(1000, getOption("timeout")),
  Ncpus = threads
)

Sys.setenv(MAKEFLAGS = paste0("-j", threads))

msg <- function(...) cat("\n==== ", sprintf(...), " ====\n", sep = "")
warn <- function(...) warning(sprintf(...), call. = FALSE)

make_bioc_repos <- function(root) {
  root <- sub("/+$", "", root)
  c(
    BioCsoft = sprintf("%s/packages/%s/bioc", root, bioc_version),
    BioCann = sprintf("%s/packages/%s/data/annotation", root, bioc_version),
    BioCexp = sprintf("%s/packages/%s/data/experiment", root, bioc_version),
    BioCworkflows = sprintf("%s/packages/%s/workflows", root, bioc_version),
    CRAN = cran_repo
  )
}

bioc_repos <- make_bioc_repos(bioc_mirrors[[1]])
active_bioc_mirror <- bioc_mirrors[[1]]
healthy_bioc_mirrors <- NULL
options(repos = bioc_repos, BioC_mirror = paste0(active_bioc_mirror, "/"))

probe_bioc_mirror <- function(root) {
  repo <- unname(make_bioc_repos(root)["BioCsoft"])
  urls <- c(
    sprintf("%s/src/contrib/PACKAGES.gz", repo),
    sprintf("%s/src/contrib/PACKAGES", repo)
  )

  for (url in urls) {
    for (attempt in seq_len(bioc_retry_attempts)) {
      tmp <- tempfile("bioc_packages_", fileext = if (grepl("\\.gz$", url)) ".gz" else ".txt")
      curl_bin <- Sys.which("curl")
      if (!nzchar(curl_bin)) stop("curl was not found while probing Bioconductor mirrors")
      status <- tryCatch(
        system2(
          curl_bin,
          args = c(
            "-L", "--fail", "--silent", "--show-error",
            "--connect-timeout", "15", "--max-time", "90",
            "--output", tmp, url
          ),
          stdout = FALSE,
          stderr = FALSE
        ),
        error = function(e) {
          warn("Mirror probe failed for %s: %s", url, conditionMessage(e))
          1L
        }
      )

      ok <- isTRUE(status == 0L) && file.exists(tmp) && isTRUE(file.info(tmp)$size > 20)
      unlink(tmp, force = TRUE)
      if (ok) return(TRUE)

      if (attempt < bioc_retry_attempts && bioc_retry_sleep > 0) {
        Sys.sleep(bioc_retry_sleep * attempt)
      }
    }
  }

  FALSE
}

discover_bioc_mirrors <- function(refresh = FALSE) {
  if (!refresh && !is.null(healthy_bioc_mirrors)) return(healthy_bioc_mirrors)

  msg("Probing Bioconductor %s software repositories", bioc_version)
  healthy <- character()
  for (root in bioc_mirrors) {
    msg("Probing Bioconductor mirror: %s", root)
    if (probe_bioc_mirror(root)) {
      msg("Bioconductor mirror is healthy: %s", root)
      healthy <- c(healthy, root)
    } else {
      warn("Bioconductor mirror is unavailable for release %s: %s", bioc_version, root)
    }
  }

  healthy_bioc_mirrors <<- unique(healthy)
  if (!length(healthy_bioc_mirrors)) {
    stop(
      "No configured Bioconductor mirror exposes the ", bioc_version,
      " software repository. Checked: ", paste(bioc_mirrors, collapse = ", ")
    )
  }

  healthy_bioc_mirrors
}

activate_bioc_mirror <- function(root) {
  active_bioc_mirror <<- sub("/+$", "", root)
  bioc_repos <<- make_bioc_repos(active_bioc_mirror)
  options(repos = bioc_repos, BioC_mirror = paste0(active_bioc_mirror, "/"))
  writeLines(active_bioc_mirror, bioc_active_mirror_file)
  msg("Active Bioconductor mirror: %s", active_bioc_mirror)
  invisible(bioc_repos)
}

download_with_retries <- function(url, dest) {
  part <- paste0(dest, ".part")
  unlink(part, force = TRUE)

  for (attempt in seq_len(bioc_retry_attempts)) {
    msg("Downloading [%d/%d]: %s", attempt, bioc_retry_attempts, url)
    status <- tryCatch(
      suppressWarnings(utils::download.file(
        url,
        part,
        mode = "wb",
        quiet = FALSE,
        method = "libcurl"
      )),
      error = function(e) {
        warn("Download failed for %s: %s", url, conditionMessage(e))
        1L
      }
    )

    if (isTRUE(status == 0L) && file.exists(part) && isTRUE(file.info(part)$size > 0)) {
      if (file.rename(part, dest)) return(TRUE)
      if (file.copy(part, dest, overwrite = TRUE)) {
        unlink(part, force = TRUE)
        return(TRUE)
      }
    }

    unlink(part, force = TRUE)
    if (attempt < bioc_retry_attempts && bioc_retry_sleep > 0) {
      Sys.sleep(bioc_retry_sleep * attempt)
    }
  }

  FALSE
}

install_bioc_tarball <- function(pkg, version, force = FALSE) {
  if (!force && pkg_ok(pkg)) {
    msg("%s %s already installed", pkg, pkg_ver(pkg))
    return(invisible(TRUE))
  }

  filename <- sprintf("%s_%s.tar.gz", pkg, version)
  cached <- file.path(bioc_cache_dir, filename)

  if (!file.exists(cached) || file.info(cached)$size <= 0) {
    mirrors <- discover_bioc_mirrors()
    downloaded <- FALSE
    for (root in mirrors) {
      url <- sprintf(
        "%s/packages/%s/bioc/src/contrib/%s",
        sub("/+$", "", root), bioc_version, filename
      )
      if (download_with_retries(url, cached)) {
        activate_bioc_mirror(root)
        downloaded <- TRUE
        break
      }
    }
    if (!downloaded) stop("Could not download Bioconductor tarball: ", filename)
  } else {
    msg("Using cached Bioconductor tarball: %s", cached)
  }

  msg("Installing cached Bioconductor tarball: %s", cached)
  install.packages(cached, repos = NULL, type = "source", dependencies = FALSE)
  if (!pkg_ok(pkg)) stop("Bioconductor tarball installed but package is unavailable: ", pkg)
  invisible(TRUE)
}

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

# MTD_R412_TMB_MATRIX_ABI_V10
strict_fresh_load_ok <- function(pkg, expected_version) {
  # Matrix/TMB binary incompatibility is reported as a warning during package
  # loading. For this stack a warning is a failed validation, not a soft notice.
  # The exact-version check also prevents R's automatic restoration of an older
  # package from being mistaken for a successful source rebuild.
  expr <- sprintf(
    paste(
      "options(warn = 2)",
      "pkg <- %s",
      "expected <- %s",
      "suppressPackageStartupMessages(library(pkg, character.only = TRUE))",
      "actual <- as.character(packageVersion(pkg))",
      "if (!identical(actual, expected)) {",
      "  stop(sprintf('Unexpected %%s version: expected %%s, found %%s', pkg, expected, actual))",
      "}",
      "cat(sprintf('MTD_STRICT_LOAD_OK:%%s:%%s\\n', pkg, actual))",
      sep = "; "
    ),
    encodeString(pkg, quote = "\""),
    encodeString(expected_version, quote = "\"")
  )

  out <- tempfile(sprintf("%s_strict_load_", pkg), fileext = ".log")
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = c("-e", shQuote(expr)),
    stdout = out,
    stderr = out
  )

  txt <- if (file.exists(out)) readLines(out, warn = FALSE) else character()
  expected_marker <- sprintf(
    "MTD_STRICT_LOAD_OK:%s:%s",
    pkg,
    expected_version
  )
  marker_count <- sum(txt == expected_marker)

  if (!identical(status, 0L) || marker_count != 1L) {
    cat(sprintf("\n==== Strict fresh R load check failed for %s ====\n", pkg))
    cat(paste(txt, collapse = "\n"), "\n")
    return(FALSE)
  }

  msg(
    "Strict fresh R load check OK for %s: %s",
    pkg,
    expected_version
  )
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

  todo <- if (force) pkgs else pkgs[!vapply(pkgs, pkg_ok, logical(1))]
  if (!length(todo)) {
    msg("Bioconductor packages already installed: %s", paste(pkgs, collapse = ", "))
    return(invisible(TRUE))
  }

  mirrors <- discover_bioc_mirrors()
  last_errors <- character()

  for (round in seq_len(bioc_retry_attempts)) {
    for (root in mirrors) {
      activate_bioc_mirror(root)
      msg(
        "Installing Bioconductor packages (round %d/%d): %s",
        round, bioc_retry_attempts, paste(todo, collapse = ", ")
      )

      tryCatch(
        BiocManager::install(
          todo,
          version = bioc_version,
          ask = FALSE,
          update = FALSE,
          force = force,
          Ncpus = threads
        ),
        error = function(e) {
          last_errors <<- c(
            last_errors,
            sprintf("%s: %s", root, conditionMessage(e))
          )
          warn("Bioconductor installation failed through %s: %s", root, conditionMessage(e))
        }
      )

      todo <- pkgs[!vapply(pkgs, pkg_ok, logical(1))]
      if (!length(todo)) {
        msg("Bioconductor package installation completed successfully")
        return(invisible(TRUE))
      }

      warn("Packages still missing after mirror %s: %s", root, paste(todo, collapse = ", "))
    }

    if (round < bioc_retry_attempts && bioc_retry_sleep > 0) {
      delay <- bioc_retry_sleep * round
      msg("Waiting %.0f seconds before the next Bioconductor retry round", delay)
      Sys.sleep(delay)
      mirrors <- discover_bioc_mirrors(refresh = TRUE)
    }
  }

  details <- if (length(last_errors)) paste(unique(last_errors), collapse = " | ") else "no R exception was raised"
  stop(
    "Bioconductor packages still missing after mirror failover: ",
    paste(todo, collapse = ", "),
    ". Errors: ", details
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
activate_bioc_mirror(discover_bioc_mirrors()[[1]])
tryCatch(
  BiocManager::install(version = bioc_version, ask = FALSE, update = FALSE),
  error = function(e) warn("Bioconductor version bootstrap failed: %s", conditionMessage(e))
)
if (!pkg_ok("BiocVersion")) install_bioc("BiocVersion")

if (bioc_smoke_test_only) {
  msg("Running focused Bioconductor recovery smoke test")
  install_bioc(c("BiocVersion", "lpsymphony"))

  smoke_pkgs <- c("BiocVersion", "lpsymphony")
  missing_smoke <- smoke_pkgs[!vapply(smoke_pkgs, pkg_ok, logical(1))]
  if (length(missing_smoke)) {
    stop("Smoke-test packages are still missing: ", paste(missing_smoke, collapse = ", "))
  }
  if (!fresh_load_ok("lpsymphony")) {
    stop("lpsymphony installed, but a fresh R process could not load it")
  }

  msg("Bioconductor recovery smoke test completed successfully")
  quit(save = "no", status = 0L)
}

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
  "car", "carData",

  # Seurat stack
  "progressr", "future", "future.apply", "globals", "listenv",
  "parallelly", "fitdistrplus", "ica", "leiden", "lmtest", "plotly",
  "RANN", "RcppAnnoy", "reticulate", "ROCR", "scattermore",
  "spatstat.data", "spatstat.utils", "spatstat.sparse",
  "spatstat.random", "spatstat.geom", "spatstat.explore",
  "deldir", "polyclip", "FNN", "RSpectra", "irlba", "Rtsne",
  "uwot", "spam",

  # MTD/other
  "ade4", "vegan", "processx", "promises",
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

# MTD_R412_DEPENDENCY_CLEANUP_V9
# MTD_R412_TMB_MATRIX_ABI_V10
# TMB contains compiled code tied to the Matrix binary interface. The
# 2022-05-15 CRAN snapshot provides TMB 1.8.1, which cannot be compiled against
# Matrix 1.6-5 in this R412 stack (undefined GET_SLOT). Pin TMB 1.9.7 from the
# official CRAN Archive, then rebuild the MTD-compatible glmmTMB 1.1.3 release.
tmb_version <- "1.9.7"
glmmtmb_version <- "1.1.3"

msg(
  "Installing TMB %s from CRAN Archive against Matrix %s",
  tmb_version,
  pkg_ver("Matrix")
)
install_archive("TMB", tmb_version, force = TRUE)
if (!strict_fresh_load_ok("TMB", tmb_version)) {
  stop(
    "TMB ", tmb_version,
    " was not installed cleanly against Matrix ", pkg_ver("Matrix"),
    ". Refusing to continue with a restored or ABI-incompatible TMB package."
  )
}

msg(
  "Installing glmmTMB %s from CRAN Archive against TMB %s and Matrix %s",
  glmmtmb_version,
  tmb_version,
  pkg_ver("Matrix")
)
install_archive("glmmTMB", glmmtmb_version, force = TRUE)
if (!strict_fresh_load_ok("glmmTMB", glmmtmb_version)) {
  stop(
    "glmmTMB ", glmmtmb_version,
    " was not installed cleanly against TMB ", tmb_version,
    " and Matrix ", pkg_ver("Matrix")
  )
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
  incompatible_ggnewscale <- file.path(patch_dir, "ggnewscale_0.5.2.tar.gz")
  if (file.exists(incompatible_ggnewscale)) {
    msg(
      "Skipping local ggnewscale 0.5.2 because it requires ggplot2 >= 3.5.0; MTD R412 uses ggplot2 %s",
      pkg_ver("ggplot2")
    )
  }

  local_patches <- c(
    "gtable_0.3.6.tar.gz",
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

if (!fresh_load_ok("ggnewscale")) {
  stop("The CRAN-snapshot ggnewscale installation is missing or incompatible")
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
install_bioc_tarball(
  pkg = "Maaslin2",
  version = "1.8.0",
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
  "ggplot2", "ggnewscale", "ggfun", "yulab.utils", "tidytree", "treeio",
  "DOSE", "fgsea", "DESeq2", "phyloseq", "microbiome", "mia",
  "ANCOMBC", "cplm", "Maaslin2", "Rgraphviz", "hdf5r",
  "cytolib", "flowCore", "cmapR",
  "lme4", "pbkrtest", "TMB", "glmmTMB",
  "sctransform", "SeuratObject", "Seurat"
)

msg("Important package versions before ggtree patch")
print(sapply(important_versions, pkg_ver))
RSCRIPT

log "Running R package installation phase 1"
THREADS="$THREADS" \
CRAN_REPO="$CRAN_REPO" \
CRAN_ARCHIVE_REPO="$CRAN_ARCHIVE_REPO" \
BIOC_VERSION="$BIOC_VERSION" \
BIOC_MIRRORS="$BIOC_MIRRORS" \
BIOC_RETRY_ATTEMPTS="$BIOC_RETRY_ATTEMPTS" \
BIOC_RETRY_SLEEP="$BIOC_RETRY_SLEEP" \
BIOC_SMOKE_TEST_ONLY="$BIOC_SMOKE_TEST_ONLY" \
BIOC_CACHE_DIR="$BIOC_CACHE_DIR" \
BIOC_ACTIVE_MIRROR_FILE="$BIOC_ACTIVE_MIRROR_FILE" \
PATCH_DIR="$PATCH_DIR" \
INSTALL_SEURAT="$INSTALL_SEURAT" \
INSTALL_LOCAL_PATCHES="$INSTALL_LOCAL_PATCHES" \
Rscript "$R_INSTALL_SCRIPT" 2>&1 | tee "$LOGDIR/03_R_install_phase1.log"

if [[ "$BIOC_SMOKE_TEST_ONLY" == "1" ]]; then
  log "Focused Bioconductor recovery smoke test passed"
  log "Active mirror: $(cat "$BIOC_ACTIVE_MIRROR_FILE" 2>/dev/null || printf 'unknown')"
  log "Log: $LOGDIR/03_R_install_phase1.log"
  exit 0
fi

# ------------------------------------------------------------
# Bioconductor source tarball download with mirror failover
# ------------------------------------------------------------
download_bioc_tarball() {
  local package="$1"
  local version="$2"
  local output="$3"
  local filename="${package}_${version}.tar.gz"
  local cached="$BIOC_CACHE_DIR/$filename"
  local mirror url attempt part
  local -a mirrors

  mkdir -p "$BIOC_CACHE_DIR" "$(dirname "$output")"

  if [[ -s "$cached" ]] && tar -tzf "$cached" >/dev/null 2>&1; then
    log "Using cached Bioconductor tarball: $cached"
    cp -f "$cached" "$output"
    return 0
  fi
  rm -f "$cached" "$cached.part"

  IFS=',' read -r -a mirrors <<< "$BIOC_MIRRORS"
  for mirror in "${mirrors[@]}"; do
    mirror="${mirror%/}"
    [[ -n "$mirror" ]] || continue
    url="$mirror/packages/$BIOC_VERSION/bioc/src/contrib/$filename"

    for ((attempt = 1; attempt <= BIOC_RETRY_ATTEMPTS; attempt++)); do
      part="$cached.part"
      rm -f "$part"
      log "Downloading $filename from $mirror (attempt $attempt/$BIOC_RETRY_ATTEMPTS)"

      if curl -L --fail --show-error \
        --connect-timeout 30 --max-time 1800 \
        --retry 2 --retry-delay 3 --retry-all-errors \
        "$url" -o "$part"; then
        if [[ -s "$part" ]] && tar -tzf "$part" >/dev/null 2>&1; then
          mv -f "$part" "$cached"
          cp -f "$cached" "$output"
          printf '%s
' "$mirror" > "$BIOC_ACTIVE_MIRROR_FILE"
          log "Downloaded and cached: $cached"
          return 0
        fi
        warn "Downloaded file is not a valid source tarball: $url"
      fi

      rm -f "$part"
      if (( attempt < BIOC_RETRY_ATTEMPTS && BIOC_RETRY_SLEEP > 0 )); then
        sleep $((BIOC_RETRY_SLEEP * attempt))
      fi
    done
  done

  die "Could not download $filename from any configured Bioconductor mirror"
}

# ------------------------------------------------------------
# ggtree compatibility patch
# ------------------------------------------------------------
patch_ggtree() {
  log "Applying ggtree 3.2.1 compatibility patch"

  local workdir
  workdir="$(mktemp -d /tmp/ggtree_patch_XXXXXX)"
  local tarball="$workdir/ggtree_3.2.1.tar.gz"

  download_bioc_tarball "ggtree" "3.2.1" "$tarball"

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
bioc_retry_attempts <- suppressWarnings(as.integer(Sys.getenv("BIOC_RETRY_ATTEMPTS", "3")))
bioc_retry_sleep <- suppressWarnings(as.numeric(Sys.getenv("BIOC_RETRY_SLEEP", "10")))
if (is.na(bioc_retry_attempts) || bioc_retry_attempts < 1L) bioc_retry_attempts <- 3L
if (is.na(bioc_retry_sleep) || bioc_retry_sleep < 0) bioc_retry_sleep <- 10
bioc_mirrors <- trimws(strsplit(Sys.getenv("BIOC_MIRRORS", "https://bioconductor.org,https://bioconductor.posit.co,https://bioconductor.statistik.tu-dortmund.de,https://ftp.gwdg.de/pub/misc/bioconductor"), ",", fixed = TRUE)[[1]])
bioc_mirrors <- unique(sub("/+$", "", bioc_mirrors[nzchar(bioc_mirrors)]))
active_mirror_file <- path.expand(Sys.getenv("BIOC_ACTIVE_MIRROR_FILE", file.path(tempdir(), "mtd_bioc_active_mirror.txt")))

options(timeout = max(1000, getOption("timeout")), Ncpus = threads)
Sys.setenv(MAKEFLAGS = paste0("-j", threads))

msg <- function(...) cat("\n==== ", sprintf(...), " ====\n", sep = "")
warn <- function(...) warning(sprintf(...), call. = FALSE)
pkg_ok <- function(pkg) requireNamespace(pkg, quietly = TRUE)
pkg_ver <- function(pkg) if (pkg_ok(pkg)) as.character(packageVersion(pkg)) else "MISSING"

make_repos <- function(root) {
  root <- sub("/+$", "", root)
  c(
    BioCsoft = sprintf("%s/packages/%s/bioc", root, bioc_version),
    BioCann = sprintf("%s/packages/%s/data/annotation", root, bioc_version),
    BioCexp = sprintf("%s/packages/%s/data/experiment", root, bioc_version),
    BioCworkflows = sprintf("%s/packages/%s/workflows", root, bioc_version),
    CRAN = cran_repo
  )
}

probe <- function(root) {
  url <- sprintf("%s/src/contrib/PACKAGES.gz", unname(make_repos(root)["BioCsoft"]))
  for (attempt in seq_len(bioc_retry_attempts)) {
    tmp <- tempfile(fileext = ".gz")
    curl_bin <- Sys.which("curl")
    if (!nzchar(curl_bin)) stop("curl was not found while probing Bioconductor mirrors")
    status <- tryCatch(
      system2(
        curl_bin,
        args = c(
          "-L", "--fail", "--silent", "--show-error",
          "--connect-timeout", "15", "--max-time", "90",
          "--output", tmp, url
        ),
        stdout = FALSE,
        stderr = FALSE
      ),
      error = function(e) 1L
    )
    ok <- isTRUE(status == 0L) && file.exists(tmp) && isTRUE(file.info(tmp)$size > 20)
    unlink(tmp, force = TRUE)
    if (ok) return(TRUE)
    if (attempt < bioc_retry_attempts && bioc_retry_sleep > 0) Sys.sleep(bioc_retry_sleep * attempt)
  }
  FALSE
}

if (!pkg_ok("BiocManager")) {
  options(repos = c(CRAN = cran_repo))
  install.packages("BiocManager", dependencies = c("Depends", "Imports", "LinkingTo"))
}

msg("Checking ggtree before installing enrichplot/clusterProfiler")
print(sapply(c("ggplot2", "treeio", "tidytree", "ggtree", "DOSE"), pkg_ver))
if (!pkg_ok("ggtree")) stop("ggtree is missing before phase 2. The ggtree patch probably failed.")

targets <- c("enrichplot", "clusterProfiler")
missing <- targets
healthy <- bioc_mirrors[vapply(bioc_mirrors, probe, logical(1))]
if (!length(healthy)) stop("No healthy Bioconductor mirror found for phase 2")

for (round in seq_len(bioc_retry_attempts)) {
  for (root in healthy) {
    options(repos = make_repos(root), BioC_mirror = paste0(root, "/"))
    writeLines(root, active_mirror_file)
    msg("Installing phase-2 packages through %s", root)
    tryCatch(
      BiocManager::install(
        missing,
        version = bioc_version,
        ask = FALSE,
        update = FALSE,
        force = TRUE,
        Ncpus = threads
      ),
      error = function(e) warn("Phase-2 installation failed through %s: %s", root, conditionMessage(e))
    )
    missing <- targets[!vapply(targets, pkg_ok, logical(1))]
    if (!length(missing)) break
  }
  if (!length(missing)) break
  if (round < bioc_retry_attempts && bioc_retry_sleep > 0) Sys.sleep(bioc_retry_sleep * round)
}

if (length(missing)) stop("Phase-2 Bioconductor packages still missing: ", paste(missing, collapse = ", "))
msg("Phase 2 versions")
print(sapply(c("ggplot2", "ggtree", "DOSE", "enrichplot", "clusterProfiler"), pkg_ver))
RSCRIPT

log "Running R package installation phase 2"
THREADS="$THREADS" \
CRAN_REPO="$CRAN_REPO" \
BIOC_VERSION="$BIOC_VERSION" \
BIOC_MIRRORS="$BIOC_MIRRORS" \
BIOC_RETRY_ATTEMPTS="$BIOC_RETRY_ATTEMPTS" \
BIOC_RETRY_SLEEP="$BIOC_RETRY_SLEEP" \
BIOC_ACTIVE_MIRROR_FILE="$BIOC_ACTIVE_MIRROR_FILE" \
Rscript "$R_PHASE2_SCRIPT" 2>&1 | tee "$LOGDIR/06_R_install_phase2_enrich_clusterProfiler.log"

# ------------------------------------------------------------
# Final validation
# ------------------------------------------------------------
VALIDATION_SCRIPT="$LOGDIR/validate_R412_MTD.R"

cat > "$VALIDATION_SCRIPT" <<'RSCRIPT'
bioc_version <- Sys.getenv("BIOC_VERSION", "3.14")
cran_repo <- Sys.getenv("CRAN_REPO", "https://packagemanager.posit.co/cran/2022-05-15")
active_mirror_file <- path.expand(Sys.getenv("BIOC_ACTIVE_MIRROR_FILE", ""))
bioc_root <- "https://bioconductor.org"
if (nzchar(active_mirror_file) && file.exists(active_mirror_file)) {
  candidate <- trimws(readLines(active_mirror_file, n = 1L, warn = FALSE))
  if (length(candidate) && nzchar(candidate)) bioc_root <- sub("/+$", "", candidate)
}

options(
  repos = c(
    BioCsoft = sprintf("%s/packages/%s/bioc", bioc_root, bioc_version),
    BioCann = sprintf("%s/packages/%s/data/annotation", bioc_root, bioc_version),
    BioCexp = sprintf("%s/packages/%s/data/experiment", bioc_root, bioc_version),
    BioCworkflows = sprintf("%s/packages/%s/workflows", bioc_root, bioc_version),
    CRAN = cran_repo
  ),
  BioC_mirror = paste0(bioc_root, "/")
)

pkgs <- c(
  "lattice", "MASS", "Matrix", "mgcv", "nlme", "survival",
  "ggplot2", "ggnewscale", "ggfun", "yulab.utils", "tidytree", "treeio", "ggtree",
  "DOSE", "enrichplot", "clusterProfiler", "fgsea",
  "DESeq2", "phyloseq", "microbiome", "mia", "ANCOMBC",
  "cplm", "Maaslin2", "flowCore", "cytolib", "cmapR", "tximeta",
  "hdf5r", "Rgraphviz", "pathview",
  "lme4", "pbkrtest", "TMB", "glmmTMB",
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

# MTD_R412_TMB_MATRIX_ABI_V11_FINAL_GUARD
# The TMB/glmmTMB stack is pinned because both packages contain compiled code
# tied to the Matrix ABI. Refuse to validate a silently restored or downgraded
# package even when require() itself returns TRUE.
expected_abi_versions <- c(
  Matrix = "1.6.5",
  TMB = "1.9.7",
  glmmTMB = "1.1.3"
)

abi_mismatch <- names(expected_abi_versions)[
  versions[names(expected_abi_versions)] != expected_abi_versions
]

if (length(abi_mismatch)) {
  details <- paste(
    sprintf(
      "%s expected=%s found=%s",
      abi_mismatch,
      expected_abi_versions[abi_mismatch],
      versions[abi_mismatch]
    ),
    collapse = "; "
  )
  stop("Matrix/TMB ABI version guard failed: ", details)
}

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
valid_result <- tryCatch(
  BiocManager::valid(),
  error = function(e) {
    warning("BiocManager::valid() could not complete its online check: ", conditionMessage(e), call. = FALSE)
    NA
  }
)
print(valid_result)

cat("\nValidation completed successfully.\n")
RSCRIPT

log "Running final validation"
CRAN_REPO="$CRAN_REPO" \
BIOC_VERSION="$BIOC_VERSION" \
BIOC_ACTIVE_MIRROR_FILE="$BIOC_ACTIVE_MIRROR_FILE" \
Rscript "$VALIDATION_SCRIPT" 2>&1 | tee "$LOGDIR/07_validation.log"

log "Done. Logs are in: $LOGDIR"

