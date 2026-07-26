#!/usr/bin/env bash
# ==============================================================================
# MTD Explorer installation checker
# Version: 2026.07.26-r11
#
# Current installation architecture:
#   - repository location is detected from this checker, never assumed as ~/MTD;
#   - installer name can be selected or discovered automatically;
#   - Kraken2 2.17.1 and Bracken are isolated in the MTD_kraken2 environment;
#   - the shared microbiome Kraken2 database is installed by default;
#   - NCBI library caches carry remote-catalog freshness/integrity metadata;
#   - NCBI taxonomy is cached once and copied with a remote checksum manifest;
#   - unmappable plasmid records may only be removed through the audited filter;
#   - the virome is a custom nonredundant RefSeq + Virus-Host DB + SIV collection;
#   - custom-host caches are not mistaken for installed kraken2DB_<TaxID> databases.
#
# The checker is read-only except for its report directory and temporary files.
# ==============================================================================
set -uo pipefail

CHECKER_VERSION="2026.07.26-r11"

SCRIPT_DIR="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &&
    pwd -P
)"

MTD_DIR="$SCRIPT_DIR"
INSTALLER_PATH=""
INSTALL_SH=""
INSTALLER_NAME=""
INSTALLER_RELATIVE=""
CONDA_PATH=""
OFFLINE_DIR=""
READ_LEN=75
MODE="full"
REPORT_DIR=""
STRICT=0
KEEP_TEMP=0
HOST_TAXID=""
NETWORK_CHECK=1

PASS_COUNT=0
WARN_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
TOTAL_COUNT=0

COLOR_RESET=""
COLOR_RED=""
COLOR_GREEN=""
COLOR_YELLOW=""
COLOR_BLUE=""
COLOR_MAGENTA=""

RESULTS_TSV=""
FULL_LOG=""
SUMMARY_TXT=""
TMP_WORK=""

usage() {
    cat <<'USAGE'
MTD Explorer installation checker

Usage:
  bash MTD_check_installation.sh [options]

Options:
  --mtd-dir PATH       MTD Explorer repository/installation directory
                       Default: directory containing this checker

  --installer PATH     Installer script
                       Default: automatically detected inside --mtd-dir

  --conda-path PATH    Miniconda directory
                       Default: value from condaPath, then $HOME/miniconda3

  --offline-dir PATH   Persistent installation cache
                       Default: value from offlineCachePath

  -r, --read-length N  Bracken read length
                       Default: 75

  --hostid TAXID       Check one installed custom-host reference
                       Default: automatically detect numeric host references

  --mode MODE          quick, full, or deep
                       quick: runtime and installed-database essentials
                       full:  quick plus source, packages, cache metadata,
                              HUMAnN databases, and audit contracts
                       deep:  full plus gzip integrity, kraken2-inspect,
                              and safe remote freshness checks

  --no-network         Skip remote freshness checks in deep mode
  --report-dir PATH    Output directory for reports
  --strict             Treat warnings as final failure
  --keep-temp          Keep temporary checker files
  --version            Print checker version
  -h, --help           Show this help

Examples:
  bash MTD_check_installation.sh --mode quick

  bash MTD_check_installation.sh \
      --mode full \
      --strict

  bash MTD_check_installation.sh \
      --mode deep \
      --offline-dir /path/to/installer-cache

Exit status:
  0  no failures; warnings allowed unless --strict is used
  1  one or more failures, or warnings with --strict
  2  invalid checker arguments
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mtd-dir)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --mtd-dir requires a value." >&2
                exit 2
            }
            MTD_DIR="$2"
            shift 2
            ;;
        --mtd-dir=*)
            MTD_DIR="${1#*=}"
            shift
            ;;
        --installer)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --installer requires a value." >&2
                exit 2
            }
            INSTALLER_PATH="$2"
            shift 2
            ;;
        --installer=*)
            INSTALLER_PATH="${1#*=}"
            shift
            ;;
        --conda-path)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --conda-path requires a value." >&2
                exit 2
            }
            CONDA_PATH="$2"
            shift 2
            ;;
        --conda-path=*)
            CONDA_PATH="${1#*=}"
            shift
            ;;
        --offline-dir|-o)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --offline-dir requires a value." >&2
                exit 2
            }
            OFFLINE_DIR="$2"
            shift 2
            ;;
        --offline-dir=*)
            OFFLINE_DIR="${1#*=}"
            shift
            ;;
        -r|--read-length)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --read-length requires a value." >&2
                exit 2
            }
            READ_LEN="$2"
            shift 2
            ;;
        --read-length=*)
            READ_LEN="${1#*=}"
            shift
            ;;
        --hostid)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --hostid requires a value." >&2
                exit 2
            }
            HOST_TAXID="$2"
            shift 2
            ;;
        --hostid=*)
            HOST_TAXID="${1#*=}"
            shift
            ;;
        --mode)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --mode requires a value." >&2
                exit 2
            }
            MODE="$2"
            shift 2
            ;;
        --mode=*)
            MODE="${1#*=}"
            shift
            ;;
        --no-network)
            NETWORK_CHECK=0
            shift
            ;;
        --report-dir)
            [[ $# -ge 2 ]] || {
                echo "ERROR: --report-dir requires a value." >&2
                exit 2
            }
            REPORT_DIR="$2"
            shift 2
            ;;
        --report-dir=*)
            REPORT_DIR="${1#*=}"
            shift
            ;;
        --strict)
            STRICT=1
            shift
            ;;
        --keep-temp)
            KEEP_TEMP=1
            shift
            ;;
        --version)
            printf '%s\n' "$CHECKER_VERSION"
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'ERROR: Unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case "$MODE" in
    quick|full|deep) ;;
    *)
        printf 'ERROR: --mode must be quick, full, or deep. Received: %s\n' \
            "$MODE" >&2
        exit 2
        ;;
esac

if ! [[ "$READ_LEN" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --read-length must be a positive integer." >&2
    exit 2
fi

if [[ -n "$HOST_TAXID" ]] &&
   ! [[ "$HOST_TAXID" =~ ^[1-9][0-9]*$ ]]
then
    echo "ERROR: --hostid must be a positive NCBI Taxonomy ID." >&2
    exit 2
fi

expand_path() {
    local value="$1"
    if [[ "$value" == "~" ]]; then
        value="$HOME"
    elif [[ "$value" == "~/"* ]]; then
        value="$HOME/${value#~/}"
    fi
    readlink -m -- "$value" 2>/dev/null || printf '%s\n' "$value"
}

resolve_installer() {
    local candidate=""
    local -a detected=()

    if [[ -n "$INSTALLER_PATH" ]]; then
        if [[ "$INSTALLER_PATH" != /* ]]; then
            INSTALLER_PATH="$MTD_DIR/$INSTALLER_PATH"
        fi
        INSTALLER_PATH="$(expand_path "$INSTALLER_PATH")"
        if [[ ! -s "$INSTALLER_PATH" ]]; then
            printf 'ERROR: Installer not found or empty: %s\n' \
                "$INSTALLER_PATH" >&2
            return 1
        fi
        printf '%s\n' "$INSTALLER_PATH"
        return 0
    fi

    for candidate in \
        "$MTD_DIR/Install.sh" \
        "$MTD_DIR/MTD_install.sh" \
        "$MTD_DIR/MTD_Explorer_install.sh" \
        "$MTD_DIR/install.sh"
    do
        if [[ -s "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done

    while IFS= read -r -d '' candidate; do
        case "$(basename "$candidate")" in
            MTD_check_installation.sh|MTD_explorer.sh|Create_custom_*.sh)
                continue
                ;;
        esac
        if grep -Fq 'offlineCachePath' "$candidate" 2>/dev/null &&
           grep -Eq '(^|[[:space:]])parse_arguments[[:space:]]*\(\)' \
               "$candidate" 2>/dev/null
        then
            detected+=("$candidate")
        fi
    done < <(
        find "$MTD_DIR" -maxdepth 1 -type f -name '*.sh' -print0 2>/dev/null
    )

    if (( ${#detected[@]} == 1 )); then
        printf '%s\n' "${detected[0]}"
        return 0
    fi

    if (( ${#detected[@]} > 1 )); then
        echo "ERROR: More than one possible installer was detected:" >&2
        printf '  %s\n' "${detected[@]}" >&2
        echo "Use --installer PATH to select one." >&2
    else
        echo "ERROR: No MTD Explorer installer was detected." >&2
    fi
    return 1
}

MTD_DIR="$(expand_path "$MTD_DIR")"
if [[ ! -d "$MTD_DIR" ]]; then
    printf 'ERROR: MTD directory does not exist: %s\n' "$MTD_DIR" >&2
    exit 2
fi

INSTALL_SH="$(resolve_installer)" || exit 2
INSTALLER_NAME="$(basename "$INSTALL_SH")"
if [[ "$INSTALL_SH" == "$MTD_DIR/"* ]]; then
    INSTALLER_RELATIVE="${INSTALL_SH#"$MTD_DIR/"}"
else
    INSTALLER_RELATIVE="$INSTALL_SH"
fi

if [[ -z "$CONDA_PATH" && -s "$MTD_DIR/condaPath" ]]; then
    CONDA_PATH="$(head -n 1 "$MTD_DIR/condaPath" | tr -d '\r\n')"
fi
: "${CONDA_PATH:=$HOME/miniconda3}"
CONDA_PATH="$(expand_path "$CONDA_PATH")"

if [[ -z "$OFFLINE_DIR" && -s "$MTD_DIR/offlineCachePath" ]]; then
    OFFLINE_DIR="$(head -n 1 "$MTD_DIR/offlineCachePath" | tr -d '\r\n')"
fi
if [[ -n "$OFFLINE_DIR" ]]; then
    OFFLINE_DIR="$(expand_path "$OFFLINE_DIR")"
fi

if [[ -z "$REPORT_DIR" ]]; then
    REPORT_DIR="$MTD_DIR/installation_check_$(date +%Y%m%d_%H%M%S)"
fi
REPORT_DIR="$(expand_path "$REPORT_DIR")"

mkdir -p "$REPORT_DIR" || {
    printf 'ERROR: Could not create report directory: %s\n' \
        "$REPORT_DIR" >&2
    exit 2
}

RESULTS_TSV="$REPORT_DIR/MTD_installation_check.tsv"
FULL_LOG="$REPORT_DIR/MTD_installation_check.log"
SUMMARY_TXT="$REPORT_DIR/MTD_installation_summary.txt"
TMP_WORK="$REPORT_DIR/.tmp"
mkdir -p "$TMP_WORK"
printf 'Status\tSection\tCheck\tDetails\n' > "$RESULTS_TSV"
: > "$FULL_LOG"

cleanup() {
    if (( KEEP_TEMP == 0 )); then
        rm -rf -- "$TMP_WORK"
    fi
}

on_signal() {
    cleanup
    exit 130
}

trap cleanup EXIT
trap on_signal INT TERM

init_colors() {
    if [[ -t 1 && -z "${NO_COLOR:-}" ]] &&
       command -v tput >/dev/null 2>&1
    then
        COLOR_RESET="$(tput sgr0 2>/dev/null || true)"
        COLOR_RED="$(tput setaf 1 2>/dev/null || true)"
        COLOR_GREEN="$(tput setaf 2 2>/dev/null || true)"
        COLOR_YELLOW="$(tput setaf 3 2>/dev/null || true)"
        COLOR_BLUE="$(tput setaf 4 2>/dev/null || true)"
        COLOR_MAGENTA="$(tput setaf 5 2>/dev/null || true)"
    fi
}

sanitize_field() {
    printf '%s' "$1" |
        tr '\t\r\n' '   ' |
        sed -E 's/[[:space:]]+/ /g; s/^ //; s/ $//'
}

log_line() {
    printf '%s\n' "$*" | tee -a "$FULL_LOG"
}

record() {
    local status="$1"
    local section="$2"
    local check="$3"
    local details="${4:-}"
    local color=""

    TOTAL_COUNT=$((TOTAL_COUNT + 1))
    case "$status" in
        PASS)
            PASS_COUNT=$((PASS_COUNT + 1))
            color="$COLOR_GREEN"
            ;;
        WARN)
            WARN_COUNT=$((WARN_COUNT + 1))
            color="$COLOR_YELLOW"
            ;;
        FAIL)
            FAIL_COUNT=$((FAIL_COUNT + 1))
            color="$COLOR_RED"
            ;;
        SKIP)
            SKIP_COUNT=$((SKIP_COUNT + 1))
            color="$COLOR_BLUE"
            ;;
    esac

    section="$(sanitize_field "$section")"
    check="$(sanitize_field "$check")"
    details="$(sanitize_field "$details")"

    printf '%s[%s]%s %-23s | %-45s | %s\n' \
        "$color" "$status" "$COLOR_RESET" \
        "$section" "$check" "$details" | tee -a "$FULL_LOG"

    printf '%s\t%s\t%s\t%s\n' \
        "$status" "$section" "$check" "$details" >> "$RESULTS_TSV"
}

file_size() {
    local path="$1"
    if [[ -e "$path" ]]; then
        du -h "$path" 2>/dev/null | awk 'NR == 1 {print $1}'
    fi
}

count_nonempty_lines() {
    local path="$1"
    if [[ -f "$path" ]]; then
        awk 'NF {count++} END {print count + 0}' "$path" 2>/dev/null
    else
        printf '0\n'
    fi
}

capture_command() {
    local output_file="$1"
    shift
    "$@" > "$output_file" 2>&1
}

compact_output() {
    local path="$1"
    local max_lines="${2:-8}"
    [[ -s "$path" ]] || return 0
    sed -E \
        -e 's/\x1B\[[0-9;]*[[:alpha:]]//g' \
        -e 's/\r//g' \
        "$path" |
        tail -n "$max_lines" |
        tr '\n' ' ' |
        sed -E 's/[[:space:]]+/ /g; s/^ //; s/ $//'
}

check_required_file() {
    local section="$1"
    local label="$2"
    local path="$3"
    if [[ -s "$path" ]]; then
        record PASS "$section" "$label" "$path ($(file_size "$path"))"
    elif [[ -e "$path" ]]; then
        record FAIL "$section" "$label" "exists but is empty: $path"
    else
        record FAIL "$section" "$label" "missing: $path"
    fi
}

check_optional_file() {
    local section="$1"
    local label="$2"
    local path="$3"
    if [[ -s "$path" ]]; then
        record PASS "$section" "$label" "$path ($(file_size "$path"))"
    else
        record SKIP "$section" "$label" "optional and absent/empty: $path"
    fi
}

check_required_dir() {
    local section="$1"
    local label="$2"
    local path="$3"
    if [[ -d "$path" ]] &&
       find "$path" -mindepth 1 -print -quit 2>/dev/null | grep -q .
    then
        record PASS "$section" "$label" "$path"
    elif [[ -d "$path" ]]; then
        record FAIL "$section" "$label" "directory is empty: $path"
    else
        record FAIL "$section" "$label" "missing: $path"
    fi
}

check_command() {
    local command_name="$1"
    local resolved=""
    resolved="$(command -v "$command_name" 2>/dev/null || true)"
    if [[ -n "$resolved" ]]; then
        record PASS "System" "command: $command_name" "$resolved"
    else
        record FAIL "System" "command: $command_name" "not found in PATH"
    fi
}

check_shell_syntax() {
    local label="$1"
    local path="$2"
    local tmp="$TMP_WORK/shell_${label//[^A-Za-z0-9_.-]/_}.txt"
    if [[ ! -f "$path" ]]; then
        record FAIL "Source syntax" "$label" "missing: $path"
        return
    fi
    if capture_command "$tmp" bash -n "$path"; then
        record PASS "Source syntax" "$label" "bash -n passed"
    else
        record FAIL "Source syntax" "$label" "$(compact_output "$tmp" 20)"
    fi
}

check_perl_syntax() {
    local label="$1"
    local path="$2"
    local tmp="$TMP_WORK/perl_${label//[^A-Za-z0-9_.-]/_}.txt"
    if [[ ! -f "$path" ]]; then
        record FAIL "Source syntax" "$label" "missing: $path"
        return
    fi
    if capture_command "$tmp" perl -c "$path"; then
        record PASS "Source syntax" "$label" "perl -c passed"
    else
        record FAIL "Source syntax" "$label" "$(compact_output "$tmp" 20)"
    fi
}

check_python_syntax() {
    local label="$1"
    local path="$2"
    local tmp="$TMP_WORK/python_${label//[^A-Za-z0-9_.-]/_}.txt"
    local pycache_root="$TMP_WORK/python_bytecode"
    if [[ ! -f "$path" ]]; then
        record FAIL "Source syntax" "$label" "missing: $path"
        return
    fi
    mkdir -p "$pycache_root"
    if capture_command "$tmp" \
        env PYTHONPYCACHEPREFIX="$pycache_root" \
        python3 -m py_compile "$path"
    then
        record PASS "Source syntax" "$label" "python3 compilation passed"
    else
        record FAIL "Source syntax" "$label" "$(compact_output "$tmp" 20)"
    fi
}

check_text_marker() {
    local section="$1"
    local label="$2"
    local path="$3"
    local marker="$4"
    local severity="${5:-FAIL}"
    if [[ ! -f "$path" ]]; then
        record "$severity" "$section" "$label" "missing: $path"
    elif grep -Fq -- "$marker" "$path"; then
        record PASS "$section" "$label" "marker found: $marker"
    else
        record "$severity" "$section" "$label" "marker missing: $marker"
    fi
}

conda_available() {
    [[ -x "$CONDA_PATH/bin/conda" ]]
}

env_prefix() {
    local env_name="$1"
    if [[ "$env_name" == "base" ]]; then
        printf '%s\n' "$CONDA_PATH"
    else
        printf '%s\n' "$CONDA_PATH/envs/$env_name"
    fi
}

env_exists() {
    local prefix=""
    prefix="$(env_prefix "$1")"
    [[ -d "$prefix/conda-meta" ]]
}

check_conda_env() {
    local env_name="$1"
    local required="${2:-1}"
    local prefix=""
    prefix="$(env_prefix "$env_name")"
    if [[ -d "$prefix/conda-meta" ]]; then
        record PASS "Conda" "environment: $env_name" "$prefix"
    elif (( required == 1 )); then
        record FAIL "Conda" "environment: $env_name" "missing: $prefix"
    else
        record SKIP "Conda" "environment: $env_name" "optional and absent: $prefix"
    fi
}

check_env_command() {
    local env_name="$1"
    local command_name="$2"
    local severity="${3:-FAIL}"
    local tmp="$TMP_WORK/env_${env_name}_${command_name//[^A-Za-z0-9_.-]/_}.txt"
    if ! env_exists "$env_name"; then
        record "$severity" "Tools/$env_name" "$command_name" \
            "environment is unavailable"
        return
    fi
    if capture_command "$tmp" "$CONDA_PATH/bin/conda" run -n "$env_name" \
        bash -c "command -v '$command_name'"
    then
        record PASS "Tools/$env_name" "$command_name" "$(compact_output "$tmp" 2)"
    else
        record "$severity" "Tools/$env_name" "$command_name" "not found"
    fi
}

check_env_version() {
    local env_name="$1"
    local package_name="$2"
    local expected_prefix="${3:-}"
    local severity="${4:-WARN}"
    local tmp="$TMP_WORK/version_${env_name}_${package_name}.json"
    local version=""

    if ! env_exists "$env_name"; then
        record "$severity" "Versions/$env_name" "$package_name" \
            "environment is unavailable"
        return
    fi

    if capture_command "$tmp" "$CONDA_PATH/bin/conda" list -n "$env_name" \
        "$package_name" --json
    then
        version="$(python3 - "$tmp" "$package_name" <<'PY_VERSION'
import json
import sys
path, wanted = sys.argv[1:]
wanted = wanted.lower().replace("_", "-")
try:
    data = json.load(open(path, encoding="utf-8"))
except Exception:
    data = []
for item in data:
    name = str(item.get("name", "")).lower().replace("_", "-")
    if name == wanted:
        print(item.get("version", ""))
        break
PY_VERSION
)"
    fi

    if [[ -z "$version" ]]; then
        record "$severity" "Versions/$env_name" "$package_name" "not installed"
    elif [[ -n "$expected_prefix" && "$version" != "$expected_prefix"* ]]; then
        record "$severity" "Versions/$env_name" "$package_name" \
            "$version (expected ${expected_prefix}*)"
    else
        record PASS "Versions/$env_name" "$package_name" \
            "$version${expected_prefix:+ (expected ${expected_prefix}*)}"
    fi
}

check_repo_r_checker() {
    local env_name="$1"
    local label="$2"
    local script_path="$3"
    local tmp="$TMP_WORK/rcheck_${env_name}_${label//[^A-Za-z0-9_.-]/_}.txt"
    if [[ ! -f "$script_path" ]]; then
        record SKIP "R/$env_name" "$label" "checker script absent"
        return
    fi
    if ! env_exists "$env_name"; then
        record FAIL "R/$env_name" "$label" "environment is unavailable"
        return
    fi
    if capture_command "$tmp" "$CONDA_PATH/bin/conda" run -n "$env_name" \
        bash "$script_path"
    then
        record PASS "R/$env_name" "$label" "$(compact_output "$tmp" 8)"
    else
        record FAIL "R/$env_name" "$label" "$(compact_output "$tmp" 30)"
    fi
}

check_source_tree() {
    local source_file=""
    local script_path=""
    local relative=""
    local manifest_count=0
    local old_path_tmp="$TMP_WORK/old_paths.txt"
    local -a required_source_files=(
        "$INSTALLER_RELATIVE"
        MTD_explorer.sh
        MTD_check_installation.sh
        Create_custom_host.sh
        Create_custom_micro.sh
        HostSpecies.csv
        Installation/MTD.yml
        Installation/MTD_fastp.yml
        Installation/MTD_humann.yml
        Installation/MTD_kraken2.yml
        Installation/MTD_orgdb.yml
        Installation/R412.yml
        Installation/mask_low_complexity.sh
        Installation/filter_kraken2_unmapped_plasmids.py
        Installation/report_kraken2_unmapped.py
        Installation/M33262_SIVMM239.fa
        aux_scripts/Kraken2/build_nonredundant_viral_fasta.py
        aux_scripts/Kraken2/build_refseq_viral_taxid_map.py
        aux_scripts/Kraken2/build_virushost_taxid_map.py
        aux_scripts/manifest_scripts/sync_ncbi_cache.py
    )

    for source_file in "${required_source_files[@]}"; do
        if [[ "$source_file" == /* ]]; then
            check_required_file "Repository source" \
                "$(basename "$source_file")" "$source_file"
        else
            check_required_file "Repository source" "$source_file" \
                "$MTD_DIR/$source_file"
        fi
    done

    check_shell_syntax "$INSTALLER_RELATIVE" "$INSTALL_SH"
    for relative in \
        MTD_explorer.sh \
        MTD_check_installation.sh \
        Create_custom_host.sh \
        Create_custom_micro.sh \
        Installation/mask_low_complexity.sh
    do
        check_shell_syntax "$relative" "$MTD_DIR/$relative"
    done

    while IFS= read -r -d '' script_path; do
        relative="${script_path#"$MTD_DIR/"}"
        manifest_count=$((manifest_count + 1))
        check_shell_syntax "$relative" "$script_path"
    done < <(
        find "$MTD_DIR/aux_scripts/manifest_scripts" \
            -maxdepth 1 -type f -name 'manifest*.sh' -print0 2>/dev/null
    )

    if (( manifest_count >= 5 )); then
        record PASS "Repository source" "manifest script inventory" \
            "$manifest_count manifest script(s) found"
    else
        record FAIL "Repository source" "manifest script inventory" \
            "only $manifest_count manifest script(s) found; expected at least 5"
    fi

    for relative in \
        Installation/filter_kraken2_unmapped_plasmids.py \
        Installation/report_kraken2_unmapped.py \
        aux_scripts/Kraken2/build_nonredundant_viral_fasta.py \
        aux_scripts/Kraken2/build_refseq_viral_taxid_map.py \
        aux_scripts/Kraken2/build_virushost_taxid_map.py \
        aux_scripts/manifest_scripts/sync_ncbi_cache.py
    do
        check_python_syntax "$relative" "$MTD_DIR/$relative"
    done

    for relative in \
        Installation/rsync_from_ncbi.pl \
        Installation/rsync_from_ncbi_archaea.pl \
        Installation/rsync_from_ncbi_bacteria.pl \
        Installation/rsync_from_ncbi_offline.pl
    do
        if [[ -f "$MTD_DIR/$relative" ]]; then
            check_perl_syntax "$relative" "$MTD_DIR/$relative"
        fi
    done

    check_text_marker "Installer contract" "dedicated Kraken environment" \
        "$INSTALL_SH" 'KRAKEN_ENV_NAME="MTD_kraken2"'
    check_text_marker "Installer contract" "Kraken environment file" \
        "$INSTALL_SH" 'Installation/MTD_kraken2.yml'
    check_text_marker "Installer contract" "manifest helper directory" \
        "$INSTALL_SH" 'MANIFEST_SCRIPTS_DIR="$dir/aux_scripts/manifest_scripts"'
    check_text_marker "Installer contract" "Kraken helper directory" \
        "$INSTALL_SH" 'KRAKEN_AUX_DIR="$dir/aux_scripts/Kraken2"'
    check_text_marker "Installer contract" "taxonomy remote manifest" \
        "$INSTALL_SH" 'KRAKEN_TAXONOMY_REMOTE_MANIFEST_NAME'
    check_text_marker "Installer contract" "plasmid audit filter" \
        "$INSTALL_SH" 'MTD_KRAKEN2_PLASMID_UNMAPPED_FILTER_V1'
    check_text_marker "Installer contract" "eukaryote HTTPS cache" \
        "$INSTALL_SH" 'MTD_NCBI_EUKARYOTE_HTTPS_CACHE_V1'
    check_text_marker "Installer contract" "Virus-Host DB pinned mirror" \
        "$INSTALL_SH" 'VIRUSHOST_MIRROR_TAG='
    check_text_marker "Installer contract" "nonredundant virome output" \
        "$INSTALL_SH" 'viral_genomes_combined_nonredundant.fna'
    check_text_marker "Installer contract" "custom virome Kraken insertion" \
        "$INSTALL_SH" '--add-to-library "$viral_library"'
    check_text_marker "Installer contract" "cache path record" \
        "$INSTALL_SH" 'offlineCachePath'
    check_text_marker "Installer contract" "Conda path record" \
        "$INSTALL_SH" 'condaPath'

    : > "$old_path_tmp"
    grep -RInE \
        --include='*.sh' --include='*.py' --include='*.pl' \
        --exclude='*.bak*' --exclude='MTD_check_installation.sh' \
        '(\$HOME|\$\{HOME\})/MTD(/|["'"'"']|$)|/home/[^/]+/MTD(/|["'"'"']|$)|(\$HOME|\$\{HOME\})/MTD-Explorer/' \
        "$INSTALL_SH" \
        "$MTD_DIR/Installation" \
        "$MTD_DIR/aux_scripts/manifest_scripts" \
        > "$old_path_tmp" 2>/dev/null || true

    if [[ -s "$old_path_tmp" ]]; then
        record FAIL "Source path contract" "hard-coded repository paths" \
            "$(compact_output "$old_path_tmp" 8)"
    else
        record PASS "Source path contract" "hard-coded repository paths" \
            "no active ~/MTD or fixed ~/MTD-Explorer path detected"
    fi

    while IFS= read -r -d '' script_path; do
        relative="${script_path#"$MTD_DIR/"}"

        if ! grep -Fq 'sync_ncbi_cache.py' "$script_path" ||
           ! grep -Fq 'MTD_MANIFEST_HELPER' "$script_path"
        then
            record FAIL "Manifest contract" "$relative" \
                "must use the shared sync_ncbi_cache.py helper"
            continue
        fi

        case "$(basename "$script_path")" in
            manifest.plasmid.sh)
                if grep -Fq 'MTD_KRAKEN2_PLASMID_CACHE' "$script_path"; then
                    record PASS "Manifest contract" "$relative" \
                        "uses shared synchronizer and dedicated plasmid cache path"
                else
                    record FAIL "Manifest contract" "$relative" \
                        "must use MTD_KRAKEN2_PLASMID_CACHE"
                fi
                ;;
            *)
                if grep -Fq 'MTD_OFFLINE_FILES_FOLDER' "$script_path"; then
                    record PASS "Manifest contract" "$relative" \
                        "uses shared synchronizer and installation cache root"
                else
                    record FAIL "Manifest contract" "$relative" \
                        "must use MTD_OFFLINE_FILES_FOLDER"
                fi
                ;;
        esac
    done < <(
        find "$MTD_DIR/aux_scripts/manifest_scripts" \
            -maxdepth 1 -type f \
            \( -name 'manifest.bacteria.sh' \
               -o -name 'manifest.archea.sh' \
               -o -name 'manifest.archaea.sh' \
               -o -name 'manifest.virus.sh' \
               -o -name 'manifest.viral.sh' \
               -o -name 'manifest.plasmid.sh' \
               -o -name 'manifest.fungi.sh' \
               -o -name 'manifest.protozoa.sh' \) \
            -print0 2>/dev/null
    )
}

check_system() {
    local command_name=""
    record PASS "Configuration" "checker version" "$CHECKER_VERSION"
    record PASS "Configuration" "mode" "$MODE"
    record PASS "Configuration" "MTD directory" "$MTD_DIR"
    record PASS "Configuration" "installer" "$INSTALL_SH"
    record PASS "Configuration" "Conda path" "$CONDA_PATH"
    if [[ -n "$OFFLINE_DIR" ]]; then
        record PASS "Configuration" "offline cache" "$OFFLINE_DIR"
    else
        record WARN "Configuration" "offline cache" \
            "not supplied and offlineCachePath is unavailable"
    fi

    for command_name in \
        awk bash bzip2 cmp cut find grep gzip perl python3 readlink sed \
        sha256sum sort stat tar tr wc
    do
        check_command "$command_name"
    done

    if [[ "$(uname -s 2>/dev/null || true)" == "Linux" ]]; then
        record PASS "System" "operating system" "Linux"
    else
        record FAIL "System" "operating system" \
            "MTD Explorer installation is supported on Linux"
    fi
}

check_conda_stack() {
    local env_name=""
    local command_name=""
    local expect_orgdb=0

    if ! conda_available; then
        record FAIL "Conda" "executable" "missing: $CONDA_PATH/bin/conda"
        return
    fi

    record PASS "Conda" "executable" "$CONDA_PATH/bin/conda"

    for env_name in \
        base MTD_fastp MTD MTD_kraken2 MTD_humann py2 halla0820 R412
    do
        check_conda_env "$env_name" 1
    done

    if env_exists MTD_orgdb || grep -Fq 'MTD_orgdb' "$INSTALL_SH" 2>/dev/null; then
        expect_orgdb=1
    fi
    check_conda_env MTD_orgdb "$expect_orgdb"

    check_env_version MTD_fastp fastp 1.3.3 WARN
    check_env_version MTD r-base 4.0.3 WARN
    check_env_version MTD hisat2 2.2.1 WARN
    check_env_version MTD_kraken2 kraken2 2.17.1 FAIL
    check_env_version MTD_kraken2 bracken 3.1p1 FAIL
    check_env_version py2 python 2.7 WARN
    check_env_version halla0820 python 3.10 WARN
    check_env_version halla0820 r-base 4.1.2 WARN
    check_env_version R412 r-base 4.1.2 WARN

    for command_name in fastp R Rscript; do
        check_env_command MTD_fastp "$command_name"
    done

    for command_name in \
        python R Rscript hisat2 hisat2-build bowtie2 samtools featureCounts \
        makeblastdb blastdbcmd blastn magicblast diamond emapper.py datasets \
        STAR rsem-calculate-expression nextflow parallel
    do
        check_env_command MTD "$command_name"
    done

    for command_name in \
        kraken2 kraken2-build kraken2-inspect k2 bracken bracken-build
    do
        check_env_command MTD_kraken2 "$command_name"
    done

    for command_name in humann metaphlan bowtie2 diamond python; do
        check_env_command MTD_humann "$command_name"
    done

    for command_name in python hclust2.py; do
        check_env_command py2 "$command_name"
    done

    for command_name in python Rscript halla; do
        check_env_command halla0820 "$command_name"
    done

    for command_name in R Rscript; do
        check_env_command R412 "$command_name"
    done

    if (( expect_orgdb == 1 )); then
        for command_name in R Rscript jq yq; do
            check_env_command MTD_orgdb "$command_name"
        done
    fi

    if [[ "$MODE" != "quick" ]]; then
        check_repo_r_checker MTD "required package checker" \
            "$MTD_DIR/update_fix/check_R_pkg.MTD.sh"
        check_repo_r_checker R412 "required package checker" \
            "$MTD_DIR/update_fix/check_R_pkg.R412.sh"
        check_repo_r_checker halla0820 "required package checker" \
            "$MTD_DIR/update_fix/check_R_pkg.halla0820.sh"
    fi
}

find_kraken_libexec() {
    local prefix="$CONDA_PATH/envs/MTD_kraken2"
    local found=""
    found="$(find "$prefix" -type f -name k2mask -perm -u+x -print \
        2>/dev/null | head -n 1)"
    if [[ -n "$found" ]]; then
        dirname "$found"
    fi
}

check_kraken_helpers() {
    local libexec=""
    local source_path=""
    local active_path=""
    local helper=""

    libexec="$(find_kraken_libexec)"
    if [[ -n "$libexec" ]]; then
        record PASS "Kraken runtime" "dynamic libexec" "$libexec"
        check_required_file "Kraken runtime" "k2mask" "$libexec/k2mask"
        check_required_file "Kraken runtime" "lookup_accession_numbers" \
            "$libexec/lookup_accession_numbers"
    else
        record FAIL "Kraken runtime" "dynamic libexec" \
            "k2mask was not found under $CONDA_PATH/envs/MTD_kraken2"
        return
    fi

    if [[ "$MODE" == "quick" ]]; then
        return
    fi

    for helper in rsync_from_ncbi.pl download_genomic_library.sh; do
        source_path="$MTD_DIR/Installation/$helper"
        active_path="$libexec/$helper"
        if [[ ! -s "$source_path" || ! -s "$active_path" ]]; then
            record FAIL "Kraken helpers" "$helper restored" \
                "source or active helper is missing"
        elif cmp -s -- "$source_path" "$active_path"; then
            record PASS "Kraken helpers" "$helper restored" \
                "active helper matches repository default"
        else
            record FAIL "Kraken helpers" "$helper restored" \
                "active helper differs from repository default"
        fi
    done
}

check_taxonomy_dir() {
    local section="$1"
    local label="$2"
    local taxonomy_dir="$3"
    local file=""

    for file in \
        names.dmp nodes.dmp nucl_gb.accession2taxid nucl_wgs.accession2taxid
    do
        check_required_file "$section" "$label/$file" "$taxonomy_dir/$file"
    done

    if [[ "$MODE" != "quick" ]]; then
        for file in nucl_gb.accession2taxid nucl_wgs.accession2taxid; do
            if [[ -s "$taxonomy_dir/$file" ]] &&
               head -n 1 "$taxonomy_dir/$file" |
                   grep -q $'^accession\taccession.version\ttaxid'
            then
                record PASS "$section" "$label $file header" "valid"
            else
                record FAIL "$section" "$label $file header" \
                    "invalid or unavailable"
            fi
        done
    fi
}

check_kraken_core() {
    local database="$1"
    local label="$2"
    local file=""

    if [[ ! -d "$database" ]]; then
        record FAIL "Kraken DB" "$label directory" "missing: $database"
        return
    fi

    record PASS "Kraken DB" "$label directory" "$database"
    for file in hash.k2d opts.k2d taxo.k2d; do
        check_required_file "Kraken DB" "$label/$file" "$database/$file"
    done
}

check_kraken_library() {
    local database="$1"
    local library="$2"
    local severity="${3:-FAIL}"
    local library_dir="$database/library/$library"
    local file=""

    if [[ ! -d "$library_dir" ]]; then
        record "$severity" "Kraken library" "$library directory" \
            "missing: $library_dir"
        return
    fi

    record PASS "Kraken library" "$library directory" "$library_dir"
    for file in library.fna prelim_map.txt; do
        if [[ -s "$library_dir/$file" ]]; then
            record PASS "Kraken library" "$library/$file" \
                "$library_dir/$file ($(file_size "$library_dir/$file"))"
        else
            record "$severity" "Kraken library" "$library/$file" \
                "missing or empty: $library_dir/$file"
        fi
    done
}

read_virome_summary_metric() {
    local summary="$1"
    local metric="$2"
    awk -F '\t' -v wanted="$metric" '
        NR > 1 && $1 == wanted {
            print $2
            exit
        }
    ' "$summary" 2>/dev/null
}

check_virome_collection() {
    local database="$1"
    local viral_cache=""
    local virushost_cache=""
    local combined=""
    local summary=""
    local details=""
    local metric=""
    local value=""
    local sample_header=""
    local sample_token=""
    local sample_id=""
    local -a required_zero_metrics=(
        records_without_accession
        records_without_taxid
    )
    local -a required_positive_metrics=(
        records_seen
        records_written
        refseq_viral_records_seen
        virushostdb_records_seen
        extra_1_records_seen
    )

    record PASS "Virome architecture" "Kraken library layout" \
        "custom nonredundant collection; database/library/viral is not expected"

    if [[ -z "$OFFLINE_DIR" ]]; then
        record SKIP "Virome source" "persistent collection" \
            "offline cache path unavailable"
        return
    fi

    viral_cache="$OFFLINE_DIR/Kraken2DB_micro/library/viral"
    virushost_cache="$OFFLINE_DIR/Ref_genomes/MTD_virus/official_current"
    combined="$viral_cache/viral_genomes_combined_nonredundant.fna"
    summary="$viral_cache/viral_genomes_combined_nonredundant.summary.tsv"
    details="$viral_cache/viral_genomes_combined_nonredundant.details.tsv"

    check_required_file "Virome source" "combined nonredundant FASTA" "$combined"
    check_required_file "Virome source" "combination summary" "$summary"
    check_required_file "Virome source" "combination details" "$details"
    check_required_file "Virome source" "RefSeq viral aggregate" \
        "$viral_cache/all_viral_genomes.fna"
    check_required_file "Virome source" "RefSeq viral assembly summary" \
        "$viral_cache/assembly_summary_viral.txt"
    check_required_file "Virome source" "RefSeq viral TaxID map" \
        "$viral_cache/refseq_viral_accession2taxid.tsv"
    check_required_file "Virome source" "Virus-Host DB FASTA" \
        "$virushost_cache/virushostdb.genomic.fna.gz"
    check_required_file "Virome source" "Virus-Host DB release" \
        "$virushost_cache/dbrel.txt"
    check_required_file "Virome source" "Virus-Host DB non-segmented list" \
        "$virushost_cache/non-segmented_virus_list.tsv"
    check_required_file "Virome source" "Virus-Host DB segmented list" \
        "$virushost_cache/segmented_virus_list.tsv"
    check_required_file "Virome source" "Virus-Host DB TaxID map" \
        "$virushost_cache/virushostdb_accession2taxid.tsv"
    check_required_file "Virome source" "additional SIV FASTA" \
        "$MTD_DIR/Installation/M33262_SIVMM239.fa"

    if [[ -s "$summary" ]]; then
        if head -n 1 "$summary" | grep -q $'^metric\tvalue$'; then
            record PASS "Virome validation" "summary header" "valid"
        else
            record FAIL "Virome validation" "summary header" \
                "expected metric<TAB>value: $summary"
        fi

        for metric in "${required_zero_metrics[@]}"; do
            value="$(read_virome_summary_metric "$summary" "$metric")"
            if [[ "$value" == "0" ]]; then
                record PASS "Virome validation" "$metric" "$value"
            elif [[ "$value" =~ ^[0-9]+$ ]]; then
                record FAIL "Virome validation" "$metric" \
                    "$value record(s); expected 0"
            else
                record FAIL "Virome validation" "$metric" \
                    "metric missing or invalid in $summary"
            fi
        done

        for metric in "${required_positive_metrics[@]}"; do
            value="$(read_virome_summary_metric "$summary" "$metric")"
            if [[ "$value" =~ ^[1-9][0-9]*$ ]]; then
                record PASS "Virome composition" "$metric" "$value"
            elif [[ "$value" == "0" ]]; then
                record FAIL "Virome composition" "$metric" \
                    "0; expected at least one input/written record"
            else
                record FAIL "Virome composition" "$metric" \
                    "metric missing or invalid in $summary"
            fi
        done
    fi

    if [[ -s "$combined" && -s "$database/seqid2taxid.map" ]]; then
        sample_header="$(grep -m 1 '^>' "$combined" 2>/dev/null || true)"
        sample_token="$(
            printf '%s\n' "$sample_header" |
                sed -E -e 's/^>//' -e 's/[[:space:]].*$//'
        )"
        sample_id="$(
            printf '%s\n' "$sample_token" |
                sed -E -e 's/^kraken:taxid\|[0-9]+\|//'
        )"
        if [[ -z "$sample_token" ]]; then
            record FAIL "Virome validation" "Kraken mapping sample" \
                "could not extract a sequence ID from the combined FASTA"
        elif awk -v token="$sample_token" -v id="$sample_id" '
                $1 == token || $1 == id {found=1; exit}
                END {exit !found}
            ' "$database/seqid2taxid.map"
        then
            record PASS "Virome validation" "Kraken mapping sample" \
                "$sample_id is represented in seqid2taxid.map"
        else
            record WARN "Virome validation" "Kraken mapping sample" \
                "sample ID was not recognized in seqid2taxid.map: $sample_id"
        fi
    elif [[ -s "$combined" ]]; then
        record FAIL "Virome validation" "Kraken mapping sample" \
            "missing or empty: $database/seqid2taxid.map"
    fi
}

check_taxonomy_manifest_alignment() {
    local database="$1"
    local cache_manifest=""
    local database_manifest="$database/taxonomy/.mtd_taxonomy_remote_md5.tsv"

    if [[ -z "$OFFLINE_DIR" ]]; then
        record SKIP "Kraken taxonomy" "remote-manifest alignment" \
            "offline cache path unavailable"
        return
    fi

    cache_manifest="$OFFLINE_DIR/Kraken2_taxonomy_cache/.mtd_taxonomy_remote_md5.tsv"
    check_required_file "Kraken taxonomy" "cache remote manifest" "$cache_manifest"
    check_required_file "Kraken taxonomy" "database remote manifest" \
        "$database_manifest"

    if [[ -s "$cache_manifest" && -s "$database_manifest" ]]; then
        if cmp -s -- "$cache_manifest" "$database_manifest"; then
            record PASS "Kraken taxonomy" "cache/database alignment" \
                "database taxonomy matches the validated cache release"
        else
            record FAIL "Kraken taxonomy" "cache/database alignment" \
                "database taxonomy manifest differs from current cache"
        fi
    fi

    check_required_file "Kraken taxonomy" "cache completion marker" \
        "$OFFLINE_DIR/Kraken2_taxonomy_cache/.mtd_taxonomy_complete"
}

check_plasmid_audit() {
    local database="$1"
    local unmapped="$database/unmapped.txt"
    local prebuild="$database/unmapped.prebuild.txt"
    local origins="$database/unmapped.prebuild.origins.tsv"
    local audit="$database/filtered_plasmid_unmapped.tsv"
    local audit_count=0
    local current_accnum=0
    local denominator=0
    local fraction=""
    local invalid_library=0

    check_required_file "Kraken mapping" "seqid2taxid.map" \
        "$database/seqid2taxid.map"

    if [[ ! -e "$unmapped" || ! -s "$unmapped" ]]; then
        record PASS "Kraken mapping" "post-build unmapped accessions" \
            "absent or empty"
    else
        record FAIL "Kraken mapping" "post-build unmapped accessions" \
            "$(count_nonempty_lines "$unmapped") accession(s): $unmapped"
    fi

    if [[ ! -e "$prebuild" || ! -s "$prebuild" ]]; then
        record PASS "Kraken mapping" "pre-build unmapped accessions" \
            "absent or empty after final preflight"
    else
        record FAIL "Kraken mapping" "pre-build unmapped accessions" \
            "$(count_nonempty_lines "$prebuild") accession(s): $prebuild"
    fi

    if [[ -s "$origins" ]]; then
        record PASS "Kraken mapping" "pre-build origin report" \
            "$origins ($(file_size "$origins"))"
    else
        record SKIP "Kraken mapping" "pre-build origin report" \
            "not retained or empty"
    fi

    if [[ ! -e "$audit" ]]; then
        record PASS "Kraken plasmid audit" "filtered accessions" \
            "no plasmid filtering audit was required"
        return
    fi

    if [[ ! -s "$audit" ]]; then
        record FAIL "Kraken plasmid audit" "audit file" \
            "exists but is empty: $audit"
        return
    fi

    if head -n 1 "$audit" |
       grep -q $'^accession\taccession_version\tsequence_id\tlibrary\treason\tfiltered_at_utc$'
    then
        record PASS "Kraken plasmid audit" "header" "valid"
    else
        record FAIL "Kraken plasmid audit" "header" "invalid: $audit"
    fi

    audit_count="$(awk 'NR > 1 && NF {count++} END {print count + 0}' "$audit")"
    invalid_library="$(awk -F '\t' 'NR > 1 && NF && $4 != "plasmid" {count++} END {print count + 0}' "$audit")"

    if (( invalid_library == 0 )); then
        record PASS "Kraken plasmid audit" "library restriction" \
            "all $audit_count audit row(s) are plasmid-only"
    else
        record FAIL "Kraken plasmid audit" "library restriction" \
            "$invalid_library non-plasmid audit row(s) detected"
    fi

    if (( audit_count <= 1000 )); then
        record PASS "Kraken plasmid audit" "maximum count" \
            "$audit_count filtered accession record(s); limit=1000"
    else
        record FAIL "Kraken plasmid audit" "maximum count" \
            "$audit_count exceeds limit=1000"
    fi

    current_accnum="$(grep -c '^ACCNUM[[:space:]]' \
        "$database/library/plasmid/prelim_map.txt" 2>/dev/null || true)"
    denominator=$((current_accnum + audit_count))
    if (( denominator > 0 )); then
        fraction="$(awk -v n="$audit_count" -v d="$denominator" \
            'BEGIN {printf "%.8f", n/d}')"
        if awk -v f="$fraction" 'BEGIN {exit !(f <= 0.01)}'; then
            record PASS "Kraken plasmid audit" "maximum fraction" \
                "$audit_count/$denominator=$fraction; limit=0.01"
        else
            record FAIL "Kraken plasmid audit" "maximum fraction" \
                "$audit_count/$denominator=$fraction exceeds 0.01"
        fi
    else
        record WARN "Kraken plasmid audit" "maximum fraction" \
            "could not determine plasmid ACCNUM denominator"
    fi
}

check_kraken_database() {
    local database="$MTD_DIR/kraken2DB_micro"
    local library=""
    local inspect_tmp="$TMP_WORK/kraken2_inspect.txt"

    check_kraken_core "$database" "microbiome"
    check_taxonomy_dir "Kraken taxonomy" "microbiome" "$database/taxonomy"

    for library in bacteria archaea protozoa fungi plasmid UniVec_Core; do
        check_kraken_library "$database" "$library" FAIL
    done
    check_virome_collection "$database"

    check_required_file "Bracken DB" "read length $READ_LEN distribution" \
        "$database/database${READ_LEN}mers.kmer_distrib"
    check_required_file "Bracken DB" "database.kraken" \
        "$database/database.kraken"

    if [[ "$MODE" != "quick" ]]; then
        check_taxonomy_manifest_alignment "$database"
        check_plasmid_audit "$database"
    fi

    if [[ "$MODE" == "deep" && -d "$database" ]]; then
        if conda_available && env_exists MTD_kraken2; then
            if capture_command "$inspect_tmp" \
                "$CONDA_PATH/bin/conda" run -n MTD_kraken2 \
                kraken2-inspect --db "$database"
            then
                record PASS "Kraken DB" "kraken2-inspect" \
                    "$(head -n 6 "$inspect_tmp" | tr '\n' ' ' | sed -E 's/[[:space:]]+/ /g')"
            else
                record FAIL "Kraken DB" "kraken2-inspect" \
                    "$(compact_output "$inspect_tmp" 20)"
            fi
        else
            record FAIL "Kraken DB" "kraken2-inspect" \
                "MTD_kraken2 environment unavailable"
        fi
    fi
}

cache_local_dir() {
    local metadata_dir="$1"
    if [[ -d "$metadata_dir/all" ]]; then
        printf '%s\n' "$metadata_dir/all"
    else
        printf '%s\n' "$metadata_dir"
    fi
}

check_catalog_hash() {
    local library="$1"
    local metadata_dir="$2"
    local catalog="$metadata_dir/remote_catalog_${library}.tsv"
    local hash_file="$metadata_dir/remote_catalog_${library}.sha256"
    local expected=""
    local observed=""

    if [[ ! -s "$catalog" || ! -s "$hash_file" ]]; then
        record FAIL "NCBI cache/$library" "remote catalog SHA-256" \
            "catalog or hash file missing"
        return
    fi

    expected="$(awk 'NF {print $1; exit}' "$hash_file")"
    observed="$(sha256sum "$catalog" 2>/dev/null | awk '{print $1}')"
    if [[ -n "$expected" && "$expected" == "$observed" ]]; then
        record PASS "NCBI cache/$library" "remote catalog SHA-256" \
            "$observed"
    else
        record FAIL "NCBI cache/$library" "remote catalog SHA-256" \
            "expected=${expected:-missing}; observed=${observed:-missing}"
    fi
}

check_cache_library() {
    local library="$1"
    local metadata_dir="$OFFLINE_DIR/Kraken2DB_micro/library/$library"
    local local_dir=""
    local names="$metadata_dir/manifest_${library}.names.txt"
    local urls="$metadata_dir/manifest_${library}.list.txt"
    local catalog="$metadata_dir/remote_catalog_${library}.tsv"
    local integrity="$metadata_dir/integrity_${library}.stat.tsv"
    local failed="$metadata_dir/failed_downloads_${library}.txt"
    local obsolete="$metadata_dir/obsolete_local_${library}.txt"
    local marker="$metadata_dir/.mtd_${library}_cache_complete"
    local expected_count=0
    local local_count=0
    local integrity_count=0
    local failed_count=0
    local obsolete_count=0
    local missing_list="$TMP_WORK/cache_${library}_missing.txt"
    local extra_list="$TMP_WORK/cache_${library}_extra.txt"
    local local_names="$TMP_WORK/cache_${library}_local_names.txt"
    local expected_names="$TMP_WORK/cache_${library}_expected_names.txt"

    if [[ ! -d "$metadata_dir" ]]; then
        record FAIL "NCBI cache/$library" "metadata directory" \
            "missing: $metadata_dir"
        return
    fi
    record PASS "NCBI cache/$library" "metadata directory" "$metadata_dir"
    local_dir="$(cache_local_dir "$metadata_dir")"
    check_required_dir "NCBI cache/$library" "local sequence directory" "$local_dir"

    for path_label in \
        "remote catalog:$catalog" \
        "catalog SHA-256:$metadata_dir/remote_catalog_${library}.sha256" \
        "name manifest:$names" \
        "URL manifest:$urls" \
        "integrity state:$integrity" \
        "completion marker:$marker"
    do
        check_required_file "NCBI cache/$library" "${path_label%%:*}" \
            "${path_label#*:}"
    done

    if [[ ! -e "$failed" ]]; then
        record FAIL "NCBI cache/$library" "failed-download record" \
            "missing: $failed"
    elif [[ -s "$failed" ]]; then
        failed_count="$(count_nonempty_lines "$failed")"
        record FAIL "NCBI cache/$library" "failed-download record" \
            "$failed_count failed item(s): $failed"
    else
        record PASS "NCBI cache/$library" "failed-download record" \
            "present and empty"
    fi

    if [[ ! -e "$obsolete" ]]; then
        record FAIL "NCBI cache/$library" "obsolete-local record" \
            "missing: $obsolete"
    elif [[ -s "$obsolete" ]]; then
        obsolete_count="$(count_nonempty_lines "$obsolete")"
        record PASS "NCBI cache/$library" "obsolete-local record" \
            "$obsolete_count obsolete item(s) removed during last synchronization"
    else
        record PASS "NCBI cache/$library" "obsolete-local record" \
            "present and empty"
    fi

    check_catalog_hash "$library" "$metadata_dir"

    if [[ -s "$names" ]]; then
        LC_ALL=C sort -u "$names" > "$expected_names"
        expected_count="$(count_nonempty_lines "$expected_names")"
    fi

    find "$local_dir" -maxdepth 1 -type f -name '*.gz' -printf '%f\n' \
        2>/dev/null | LC_ALL=C sort -u > "$local_names"
    local_count="$(count_nonempty_lines "$local_names")"

    if (( expected_count > 0 && local_count == expected_count )); then
        record PASS "NCBI cache/$library" "catalog/local count" \
            "$local_count local file(s); $expected_count expected"
    elif (( expected_count > 0 )); then
        record FAIL "NCBI cache/$library" "catalog/local count" \
            "$local_count local file(s); $expected_count expected"
    else
        record FAIL "NCBI cache/$library" "catalog/local count" \
            "name manifest is empty or invalid"
    fi

    comm -23 "$expected_names" "$local_names" > "$missing_list" 2>/dev/null || true
    comm -13 "$expected_names" "$local_names" > "$extra_list" 2>/dev/null || true
    if [[ ! -s "$missing_list" ]]; then
        record PASS "NCBI cache/$library" "missing local files" "none"
    else
        record FAIL "NCBI cache/$library" "missing local files" \
            "$(count_nonempty_lines "$missing_list") missing; first: $(head -n 1 "$missing_list")"
    fi
    if [[ ! -s "$extra_list" ]]; then
        record PASS "NCBI cache/$library" "unexpected local files" "none"
    else
        record FAIL "NCBI cache/$library" "unexpected local files" \
            "$(count_nonempty_lines "$extra_list") unexpected; first: $(head -n 1 "$extra_list")"
    fi

    integrity_count="$(count_nonempty_lines "$integrity")"
    if (( expected_count > 0 && integrity_count == expected_count )); then
        record PASS "NCBI cache/$library" "integrity-state coverage" \
            "$integrity_count/$expected_count files recorded"
    else
        record FAIL "NCBI cache/$library" "integrity-state coverage" \
            "$integrity_count/$expected_count files recorded"
    fi

    if [[ -s "$marker" ]] &&
       grep -Fq 'status=complete' "$marker" &&
       grep -Fq "library=$library" "$marker"
    then
        record PASS "NCBI cache/$library" "completion marker contents" \
            "status and library fields are valid"
    else
        record FAIL "NCBI cache/$library" "completion marker contents" \
            "invalid or incomplete marker"
    fi

    if [[ "$library" != "plasmid" ]]; then
        check_required_file "NCBI cache/$library" "assembly summary" \
            "$metadata_dir/assembly_summary_${library}.txt"
    else
        check_optional_file "NCBI cache/$library" "remote release" \
            "$metadata_dir/remote_release_${library}.txt"
    fi

    if [[ "$MODE" == "deep" ]]; then
        check_cache_gzip_integrity "$library" "$local_dir"
        check_remote_cache_freshness "$library" "$metadata_dir"
    fi
}

check_cache_gzip_integrity() {
    local library="$1"
    local local_dir="$2"
    local failures="$TMP_WORK/gzip_failures_${library}.txt"
    local tested=0
    local file=""
    : > "$failures"

    while IFS= read -r -d '' file; do
        tested=$((tested + 1))
        if ! gzip -t "$file" 2>/dev/null; then
            printf '%s\n' "$file" >> "$failures"
        fi
    done < <(
        find "$local_dir" -maxdepth 1 -type f -name '*.gz' -print0 2>/dev/null
    )

    if (( tested == 0 )); then
        record FAIL "NCBI cache/$library" "deep gzip integrity" \
            "no gzip files found"
    elif [[ -s "$failures" ]]; then
        record FAIL "NCBI cache/$library" "deep gzip integrity" \
            "$(count_nonempty_lines "$failures")/$tested failed; first: $(head -n 1 "$failures")"
    else
        record PASS "NCBI cache/$library" "deep gzip integrity" \
            "$tested file(s) passed gzip -t"
    fi
}

manifest_for_library() {
    local library="$1"
    local source_dir="$MTD_DIR/aux_scripts/manifest_scripts"
    local candidate=""
    case "$library" in
        archaea)
            for candidate in manifest.archaea.sh manifest.archea.sh; do
                [[ -f "$source_dir/$candidate" ]] && {
                    printf '%s\n' "$source_dir/$candidate"
                    return 0
                }
            done
            ;;
        viral)
            for candidate in manifest.viral.sh manifest.virus.sh; do
                [[ -f "$source_dir/$candidate" ]] && {
                    printf '%s\n' "$source_dir/$candidate"
                    return 0
                }
            done
            ;;
        *)
            candidate="manifest.${library}.sh"
            [[ -f "$source_dir/$candidate" ]] && {
                printf '%s\n' "$source_dir/$candidate"
                return 0
            }
            ;;
    esac
    return 1
}

check_remote_cache_freshness() {
    local library="$1"
    local metadata_dir="$2"
    local manifest=""
    local tmp="$TMP_WORK/remote_check_${library}.txt"

    if (( NETWORK_CHECK == 0 )); then
        record SKIP "NCBI cache/$library" "remote freshness" \
            "disabled with --no-network"
        return
    fi

    manifest="$(manifest_for_library "$library" 2>/dev/null || true)"
    if [[ -z "$manifest" ]]; then
        record SKIP "NCBI cache/$library" "remote freshness" \
            "source manifest was not found"
        return
    fi

    if ! grep -Eq 'MTD_NCBI_CHECK_ONLY|--check-only' "$manifest"; then
        record SKIP "NCBI cache/$library" "remote freshness" \
            "manifest does not advertise read-only check-only support"
        return
    fi

    if capture_command "$tmp" env \
        MTD_NCBI_CHECK_ONLY=1 \
        MTD_MANIFEST_HELPER="$MTD_DIR/aux_scripts/manifest_scripts/sync_ncbi_cache.py" \
        MTD_OFFLINE_FILES_FOLDER="$OFFLINE_DIR" \
        bash "$manifest"
    then
        if grep -Fq 'check-only found the cache current' "$tmp"; then
            record PASS "NCBI cache/$library" "remote freshness" \
                "live NCBI catalog matches the local cache"
        elif grep -Fq 'check-only found cache changes' "$tmp"; then
            record WARN "NCBI cache/$library" "remote freshness" \
                "remote changes detected; rerun the installer synchronization"
        else
            record WARN "NCBI cache/$library" "remote freshness" \
                "check completed but status text was not recognized: $(compact_output "$tmp" 8)"
        fi
    else
        record WARN "NCBI cache/$library" "remote freshness" \
            "network/check-only command failed: $(compact_output "$tmp" 15)"
    fi
}

check_installation_cache() {
    local library=""
    if [[ -z "$OFFLINE_DIR" ]]; then
        record FAIL "Installation cache" "path" \
            "offline cache path is unavailable"
        return
    fi
    if [[ ! -d "$OFFLINE_DIR" ]]; then
        record FAIL "Installation cache" "path" "missing: $OFFLINE_DIR"
        return
    fi

    record PASS "Installation cache" "path" "$OFFLINE_DIR"
    for library in bacteria archaea protozoa fungi viral plasmid; do
        check_cache_library "$library"
    done

    check_required_dir "Installation cache" "Kraken taxonomy cache" \
        "$OFFLINE_DIR/Kraken2_taxonomy_cache"
    check_required_dir "Installation cache" "HUMAnN cache" \
        "$OFFLINE_DIR/HUMAnN"
    check_required_dir "Installation cache" "eggNOG cache" \
        "$OFFLINE_DIR/eggNOG/emapperdb-5.0.2"
}

check_humann_databases() {
    local root="$MTD_DIR/HUMAnN/ref_database"
    local choco="$root/chocophlan"
    local uniref="$root/uniref"
    local mapping="$root/utility_mapping"
    local metaphlan="$root/metaphlan"
    local count=0

    check_required_dir "HUMAnN" "ChocoPhlAn directory" "$choco"
    check_required_dir "HUMAnN" "UniRef directory" "$uniref"
    check_required_dir "HUMAnN" "utility mapping directory" "$mapping"
    check_required_dir "MetaPhlAn" "database directory" "$metaphlan"
    check_required_file "HUMAnN" "completion marker" \
        "$root/.mtd_humann_databases_complete"
    check_required_file "HUMAnN" "UniRef90 DIAMOND database" \
        "$uniref/uniref90_201901b_full.dmnd"

    count="$(find "$choco" -type f -name '*.ffn.gz' -size +0c \
        2>/dev/null | awk 'END {print NR + 0}')"
    if (( count > 0 )); then
        record PASS "HUMAnN" "ChocoPhlAn files" "$count file(s)"
    else
        record FAIL "HUMAnN" "ChocoPhlAn files" "none found"
    fi

    count="$(find "$mapping" -type f -size +0c 2>/dev/null |
        awk 'END {print NR + 0}')"
    if (( count > 0 )); then
        record PASS "HUMAnN" "utility mapping files" "$count file(s)"
    else
        record FAIL "HUMAnN" "utility mapping files" "none found"
    fi

    count="$(find "$metaphlan" -maxdepth 1 -type f \
        -name 'mpa_vJun23_CHOCOPhlAnSGB_202403*' -size +0c \
        2>/dev/null | awk 'END {print NR + 0}')"
    if (( count > 0 )); then
        record PASS "MetaPhlAn" "pinned database index" \
            "$count matching file(s)"
    else
        record FAIL "MetaPhlAn" "pinned database index" \
            "mpa_vJun23_CHOCOPhlAnSGB_202403 files not found"
    fi
}

check_one_custom_host() {
    local taxid="$1"
    local path=""
    local candidate=""

    for candidate in \
        "$MTD_DIR/kraken2DB_${taxid}" \
        "$MTD_DIR/kraken2DB_host_${taxid}"
    do
        if [[ -d "$candidate" ]]; then
            path="$candidate"
            break
        fi
    done

    if [[ -z "$path" ]]; then
        record WARN "Custom host" "TaxID $taxid" \
            "no installed kraken2DB_${taxid} database was found in $MTD_DIR"
        return
    fi

    record PASS "Custom host" "TaxID $taxid directory" "$path"
    check_kraken_core "$path" "host TaxID $taxid"
    check_required_file "Custom host" "TaxID $taxid genome FASTA" \
        "$path/genome_${taxid}.fa"
    check_required_file "Custom host" "TaxID $taxid inspect summary" \
        "$path/kraken2_inspect_summary.txt"
}

check_installed_custom_hosts() {
    local -A seen=()
    local taxid=""
    local path=""
    local base=""

    if [[ -n "$HOST_TAXID" ]]; then
        check_one_custom_host "$HOST_TAXID"
        return
    fi

    while IFS= read -r -d '' path; do
        base="$(basename "$path")"
        if [[ "$base" =~ ^kraken2DB_([1-9][0-9]*)$ ]]; then
            taxid="${BASH_REMATCH[1]}"
            seen["$taxid"]=1
        elif [[ "$base" =~ ^kraken2DB_host_([1-9][0-9]*)$ ]]; then
            taxid="${BASH_REMATCH[1]}"
            seen["$taxid"]=1
        fi
    done < <(
        find "$MTD_DIR" -mindepth 1 -maxdepth 1 -type d \
            \( -name 'kraken2DB_[1-9]*' \
               -o -name 'kraken2DB_host_[1-9]*' \) \
            -print0 2>/dev/null
    )

    if (( ${#seen[@]} == 0 )); then
        record SKIP "Custom host" "installed references" \
            "none detected; create them separately when needed"
        return
    fi

    while IFS= read -r taxid; do
        check_one_custom_host "$taxid"
    done < <(printf '%s\n' "${!seen[@]}" | sort -n)
}

write_summary() {
    local final_status="PASS"
    local exit_status=0

    if (( FAIL_COUNT > 0 )); then
        final_status="FAIL"
        exit_status=1
    elif (( STRICT == 1 && WARN_COUNT > 0 )); then
        final_status="FAIL (strict warnings)"
        exit_status=1
    elif (( WARN_COUNT > 0 )); then
        final_status="PASS WITH WARNINGS"
    fi

    cat > "$SUMMARY_TXT" <<SUMMARY
MTD Explorer installation checker
Version: $CHECKER_VERSION
Date: $(date --iso-8601=seconds 2>/dev/null || date)
Mode: $MODE
Final status: $final_status

MTD directory: $MTD_DIR
Installer: $INSTALL_SH
Conda path: $CONDA_PATH
Offline cache: ${OFFLINE_DIR:-not available}
Report directory: $REPORT_DIR

Checks: $TOTAL_COUNT
PASS: $PASS_COUNT
WARN: $WARN_COUNT
FAIL: $FAIL_COUNT
SKIP: $SKIP_COUNT

Detailed TSV: $RESULTS_TSV
Full log: $FULL_LOG
SUMMARY

    echo | tee -a "$FULL_LOG"
    log_line "============================================================"
    log_line "MTD Explorer installation check summary"
    log_line "Version: $CHECKER_VERSION"
    log_line "Mode: $MODE"
    log_line "Status: $final_status"
    log_line "PASS=$PASS_COUNT WARN=$WARN_COUNT FAIL=$FAIL_COUNT SKIP=$SKIP_COUNT"
    log_line "Reports: $REPORT_DIR"
    log_line "============================================================"

    return "$exit_status"
}

main() {
    init_colors

    log_line "============================================================"
    log_line "MTD Explorer installation checker $CHECKER_VERSION"
    log_line "============================================================"

    check_system

    if [[ "$MODE" != "quick" ]]; then
        check_source_tree
    fi

    check_conda_stack
    check_kraken_helpers
    check_kraken_database

    if [[ "$MODE" != "quick" ]]; then
        check_installation_cache
        check_humann_databases
    fi

    check_installed_custom_hosts

    if write_summary; then
        return 0
    fi
    return 1
}

main "$@"