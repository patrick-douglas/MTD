#!/usr/bin/env bash
# MTD Create_custom_host benchmark runner
set -Eeuo pipefail

PROGRAM_NAME="$(basename "$0")"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd -P)"

MACHINE=""
TAXID=""
RUN_NUMBER="1"
INTERVAL="5"
OUTPUT_ROOT="${HOME}/MTD_custom_host_benchmarks"
OFFLINE_FOLDER=""
SKIP_ORGDB=0
REBUILD_TAXONOMY=0
NO_CONFIRM=0
EXTRA_ARGS=()

usage() {
    cat <<USAGE
Usage:
  bash benchmark/$PROGRAM_NAME \\
    --machine NAME \\
    --taxid TAXID \\
    [options]

Required:
  --machine NAME              Machine label, e.g. master or s5
  --taxid TAXID               NCBI Taxon ID, e.g. 6526

Options:
  --run-number N              Repetition number [default: 1]
  --interval SECONDS          Resource sampling interval [default: 5]
  --output-root DIR           Benchmark result root
                              [default: \$HOME/MTD_custom_host_benchmarks]
  --offline-folder DIR        Override the cache in offlineCachePath
  --skip-orgdb                Skip only OrgDb package construction
  --rebuild-taxonomy-cache    Also rebuild the large shared Kraken taxonomy cache
  --yes                       Do not request BUILD confirmation
  -h, --help                  Show this help

Additional Create_custom_host.sh arguments may be placed after --.

Examples:
  bash benchmark/$PROGRAM_NAME \\
    --machine master \\
    --taxid 6526 \\
    --run-number 1

  bash benchmark/$PROGRAM_NAME \\
    --machine master \\
    --taxid 6526 \\
    --offline-folder /data/MTD_cache \\
    -- --orgdb-version 1.0.0

Benchmark definition:
  - Existing outputs for this TaxID are rebuilt by Create_custom_host.sh.
  - Species downloads and universal shared caches are preserved unless
    --rebuild-taxonomy-cache is explicitly selected.
  - --force-orgdb is used by default for a complete functional-reference run.
  - Create_custom_host.sh currently uses nproc; the detected logical CPU count
    is recorded in the report.
USAGE
}

die() {
    printf '[STOP] %s\n' "$*" >&2
    exit 1
}

safe_name() {
    printf '%s' "$1" | tr -cs 'A-Za-z0-9._-' '_'
}

human_size() {
    local path="$1"
    if [[ -e "$path" ]]; then
        du -sh "$path" 2>/dev/null | awk '{print $1}'
    else
        printf 'missing\n'
    fi
}

file_count() {
    local path="$1"
    if [[ -d "$path" ]]; then
        find "$path" -type f 2>/dev/null | wc -l
    elif [[ -f "$path" ]]; then
        printf '1\n'
    else
        printf '0\n'
    fi
}

while (($# > 0)); do
    case "$1" in
        --machine)
            (($# >= 2)) || die "--machine requires a value."
            MACHINE="$2"
            shift 2
            ;;
        --taxid)
            (($# >= 2)) || die "--taxid requires a value."
            TAXID="$2"
            shift 2
            ;;
        --run-number)
            (($# >= 2)) || die "--run-number requires a value."
            RUN_NUMBER="$2"
            shift 2
            ;;
        --interval)
            (($# >= 2)) || die "--interval requires a value."
            INTERVAL="$2"
            shift 2
            ;;
        --output-root)
            (($# >= 2)) || die "--output-root requires a value."
            OUTPUT_ROOT="$2"
            shift 2
            ;;
        --offline-folder)
            (($# >= 2)) || die "--offline-folder requires a value."
            OFFLINE_FOLDER="$2"
            shift 2
            ;;
        --skip-orgdb)
            SKIP_ORGDB=1
            shift
            ;;
        --rebuild-taxonomy-cache)
            REBUILD_TAXONOMY=1
            shift
            ;;
        --yes)
            NO_CONFIRM=1
            shift
            ;;
        --)
            shift
            EXTRA_ARGS=("$@")
            break
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

[[ -n "$MACHINE" ]] || die "--machine is required."
[[ "$TAXID" =~ ^[0-9]+$ ]] || die "--taxid must be a positive integer."
[[ "$RUN_NUMBER" =~ ^[0-9]+$ ]] || die "--run-number must be an integer."
[[ "$INTERVAL" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "--interval must be numeric."

BUILDER="$REPO_DIR/Create_custom_host.sh"
MONITOR="$SCRIPT_DIR/MTD_benchmark_install.sh"
REPORTER="$SCRIPT_DIR/MTD_custom_host_benchmark_report.py"
HOST_CSV="$REPO_DIR/HostSpecies.csv"

[[ -s "$BUILDER" ]] || die "Missing Create_custom_host.sh: $BUILDER"
[[ -s "$MONITOR" ]] || die "Missing benchmark monitor: $MONITOR"
[[ -s "$REPORTER" ]] || die "Missing custom-host reporter: $REPORTER"
[[ -s "$HOST_CSV" ]] || die "Missing HostSpecies.csv: $HOST_CSV"

bash -n "$BUILDER"
bash -n "$MONITOR"
bash -n "$0"
python3 -m py_compile "$REPORTER"

if [[ -z "$OFFLINE_FOLDER" ]]; then
    CACHE_POINTER="$REPO_DIR/offlineCachePath"
    [[ -s "$CACHE_POINTER" ]] ||
        die "offlineCachePath is missing. Use --offline-folder DIR."
    IFS= read -r OFFLINE_FOLDER < "$CACHE_POINTER"
    OFFLINE_FOLDER="${OFFLINE_FOLDER%$'\r'}"
fi

[[ -n "$OFFLINE_FOLDER" ]] || die "The persistent cache path is empty."
mkdir -p "$OFFLINE_FOLDER"
OFFLINE_FOLDER="$(readlink -f "$OFFLINE_FOLDER")"
OUTPUT_ROOT="${OUTPUT_ROOT/#\~/$HOME}"
mkdir -p "$OUTPUT_ROOT"
OUTPUT_ROOT="$(readlink -f "$OUTPUT_ROOT")"

LOGICAL_CPUS="$(nproc)"
MACHINE_SAFE="$(safe_name "$MACHINE")"
LABEL="create_custom_host_taxid${TAXID}_${MACHINE_SAFE}_run${RUN_NUMBER}"

SPECIES_CACHE="$OFFLINE_FOLDER/Customized_hosts/$TAXID"
TAXONOMY_CACHE="$OFFLINE_FOLDER/Kraken2_taxonomy_cache"
EGGNOG_CACHE="$OFFLINE_FOLDER/eggNOG/emapperdb-5.0.2"
NCBI_GENE_CACHE="$OFFLINE_FOLDER/NCBI_gene/gene2ensembl.gz"

KRAKEN_OUT="$REPO_DIR/kraken2DB_${TAXID}"
REF_OUT="$REPO_DIR/ref_${TAXID}"
HISAT_OUT="$REPO_DIR/hisat2_index_${TAXID}"
BLAST_OUT="$REPO_DIR/blastdb_${TAXID}"
ORGDB_BUILD="$REPO_DIR/build/orgdb_gold/${TAXID}"
CUSTOM_R_LIB="$REPO_DIR/custom_R_libs"

species_files="$(file_count "$SPECIES_CACHE")"
if (( species_files > 0 )); then
    SPECIES_CACHE_MODE="warm"
else
    SPECIES_CACHE_MODE="cold"
fi

required_taxonomy=(
    names.dmp
    nodes.dmp
    nucl_gb.accession2taxid
    nucl_wgs.accession2taxid
)
taxonomy_complete=1
for required in "${required_taxonomy[@]}"; do
    [[ -s "$TAXONOMY_CACHE/taxonomy/$required" ]] || taxonomy_complete=0
done
if (( taxonomy_complete == 1 )); then
    TAXONOMY_CACHE_MODE="warm"
else
    TAXONOMY_CACHE_MODE="cold_or_incomplete"
fi

if [[ -d "$EGGNOG_CACHE" ]] && [[ "$(file_count "$EGGNOG_CACHE")" -gt 0 ]]; then
    EGGNOG_CACHE_MODE="warm"
else
    EGGNOG_CACHE_MODE="missing"
fi

if [[ -s "$NCBI_GENE_CACHE" ]]; then
    NCBI_GENE_CACHE_MODE="warm"
else
    NCBI_GENE_CACHE_MODE="cold"
fi

CACHE_PROFILE="species-${SPECIES_CACHE_MODE}_taxonomy-${TAXONOMY_CACHE_MODE}_eggnog-${EGGNOG_CACHE_MODE}_ncbi-${NCBI_GENE_CACHE_MODE}"
LABEL="${LABEL}_$(safe_name "$CACHE_PROFILE")"

STATE_DIR="$(mktemp -d "${TMPDIR:-/tmp}/mtd_custom_host_benchmark.XXXXXX")"
trap 'rm -rf "$STATE_DIR"' EXIT
CACHE_STATE_BEFORE="$STATE_DIR/cache_state_before.tsv"

{
    printf 'component\tpath\tstate\tsize\tfile_count\n'
    printf 'species_reference\t%s\t%s\t%s\t%s\n' \
        "$SPECIES_CACHE" "$SPECIES_CACHE_MODE" \
        "$(human_size "$SPECIES_CACHE")" "$(file_count "$SPECIES_CACHE")"
    printf 'kraken_taxonomy\t%s\t%s\t%s\t%s\n' \
        "$TAXONOMY_CACHE" "$TAXONOMY_CACHE_MODE" \
        "$(human_size "$TAXONOMY_CACHE")" "$(file_count "$TAXONOMY_CACHE")"
    printf 'eggnog_database\t%s\t%s\t%s\t%s\n' \
        "$EGGNOG_CACHE" "$EGGNOG_CACHE_MODE" \
        "$(human_size "$EGGNOG_CACHE")" "$(file_count "$EGGNOG_CACHE")"
    printf 'ncbi_gene2ensembl\t%s\t%s\t%s\t%s\n' \
        "$NCBI_GENE_CACHE" "$NCBI_GENE_CACHE_MODE" \
        "$(human_size "$NCBI_GENE_CACHE")" "$(file_count "$NCBI_GENE_CACHE")"
} > "$CACHE_STATE_BEFORE"

HOST_ROW="$(
    python3 - "$HOST_CSV" "$TAXID" <<'PY'
import csv
import sys

csv_path, taxid = sys.argv[1:]
with open(csv_path, newline="", encoding="utf-8-sig") as handle:
    for row in csv.DictReader(handle):
        if str(row.get("Taxon_ID", "")).strip() == taxid:
            scientific = (
                row.get("Scientific_name")
                or row.get("Scientific_Name")
                or row.get("Species")
                or ""
            ).strip()
            reference = str(row.get("Reference_TaxID", "")).strip()
            orgdb = str(row.get("OrgDb", "")).strip()
            print("\t".join([scientific, reference, orgdb]))
            break
PY
)"

if [[ -n "$HOST_ROW" ]]; then
    IFS=$'\t' read -r SCIENTIFIC_NAME REFERENCE_TAXID ORGDB_PACKAGE <<< "$HOST_ROW"
else
    SCIENTIFIC_NAME=""
    REFERENCE_TAXID=""
    ORGDB_PACKAGE=""
fi

printf '\n============================================================\n'
printf 'MTD Create_custom_host benchmark\n'
printf '============================================================\n'
printf 'Repository:             %s\n' "$REPO_DIR"
printf 'Machine:                %s\n' "$MACHINE"
printf 'TaxID:                  %s\n' "$TAXID"
printf 'Scientific name:        %s\n' "${SCIENTIFIC_NAME:-not resolved from CSV}"
printf 'Reference TaxID:        %s\n' "${REFERENCE_TAXID:-not listed}"
printf 'OrgDb package:          %s\n' "${ORGDB_PACKAGE:-not listed}"
printf 'Logical CPUs (nproc):   %s\n' "$LOGICAL_CPUS"
printf 'Sampling interval:      %s s\n' "$INTERVAL"
printf 'Persistent cache:       %s\n' "$OFFLINE_FOLDER"
printf 'Species cache:          %s\n' "$SPECIES_CACHE_MODE"
printf 'Kraken taxonomy cache:  %s\n' "$TAXONOMY_CACHE_MODE"
printf 'eggNOG database:        %s\n' "$EGGNOG_CACHE_MODE"
printf 'NCBI gene cache:        %s\n' "$NCBI_GENE_CACHE_MODE"
printf 'Build OrgDb package:    %s\n' "$([[ "$SKIP_ORGDB" == 1 ]] && echo no || echo yes)"
printf 'Rebuild taxonomy cache: %s\n' "$([[ "$REBUILD_TAXONOMY" == 1 ]] && echo yes || echo no)"
printf 'Benchmark output root:  %s\n' "$OUTPUT_ROOT"
printf '\nExisting TaxID-specific outputs will be removed and rebuilt:\n'
printf '  %s\n' "$KRAKEN_OUT" "$REF_OUT" "$HISAT_OUT" "$BLAST_OUT" "$ORGDB_BUILD"
printf '\nThe persistent species downloads and shared databases are preserved.\n'
if [[ "$REBUILD_TAXONOMY" == 1 ]]; then
    printf '[WARNING] The shared Kraken taxonomy cache will also be removed and rebuilt.\n'
fi
printf '============================================================\n\n'

if [[ "$NO_CONFIRM" != 1 ]]; then
    read -r -p "Type BUILD to start the destructive rebuild benchmark: " answer
    [[ "$answer" == "BUILD" ]] || die "Confirmation was not BUILD; nothing was changed."
fi

COMMAND=(
    bash "$BUILDER"
    --ncbi-taxon-id "$TAXID"
    --offline-folder "$OFFLINE_FOLDER"
    --force-orgdb
)

if [[ "$SKIP_ORGDB" == 1 ]]; then
    COMMAND+=(--skip-orgdb)
fi
if [[ "$REBUILD_TAXONOMY" == 1 ]]; then
    COMMAND+=(--rebuild-kraken-taxonomy-cache)
fi
if ((${#EXTRA_ARGS[@]} > 0)); then
    COMMAND+=("${EXTRA_ARGS[@]}")
fi

before_dirs="$STATE_DIR/before_dirs.txt"
find "$OUTPUT_ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort > "$before_dirs"

set +e
(
    cd "$REPO_DIR"
    bash "$MONITOR" \
        --label "$LABEL" \
        --interval "$INTERVAL" \
        --output-root "$OUTPUT_ROOT" \
        --watch-path "$KRAKEN_OUT" \
        --watch-path "$REF_OUT" \
        --watch-path "$HISAT_OUT" \
        --watch-path "$BLAST_OUT" \
        --watch-path "$ORGDB_BUILD" \
        --watch-path "$SPECIES_CACHE" \
        --watch-path "$TAXONOMY_CACHE" \
        -- "${COMMAND[@]}"
)
benchmark_status=$?
set -e

after_dirs="$STATE_DIR/after_dirs.txt"
new_dirs="$STATE_DIR/new_dirs.txt"
find "$OUTPUT_ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort > "$after_dirs"
comm -13 "$before_dirs" "$after_dirs" > "$new_dirs"

RUN_DIR=""
while IFS= read -r candidate; do
    [[ "$candidate" == "$(safe_name "$LABEL")__"* ]] || continue
    RUN_DIR="$OUTPUT_ROOT/$candidate"
done < "$new_dirs"

if [[ -z "$RUN_DIR" ]]; then
    RUN_DIR="$(
        find "$OUTPUT_ROOT" -mindepth 1 -maxdepth 1 -type d \
            -name "$(safe_name "$LABEL")__*" \
            -printf '%T@\t%p\n' |
        sort -n |
        tail -n 1 |
        cut -f2-
    )"
fi

[[ -n "$RUN_DIR" && -d "$RUN_DIR" ]] ||
    die "Could not locate the benchmark output directory."

cp "$CACHE_STATE_BEFORE" "$RUN_DIR/cache_state_before.tsv"

python3 "$REPORTER" \
    --run-dir "$RUN_DIR" \
    --repo "$REPO_DIR" \
    --taxid "$TAXID" \
    --cache "$OFFLINE_FOLDER" \
    --machine "$MACHINE" \
    --run-number "$RUN_NUMBER" \
    --logical-cpus "$LOGICAL_CPUS" \
    --scientific-name "$SCIENTIFIC_NAME" \
    --reference-taxid "$REFERENCE_TAXID" \
    --orgdb-package "$ORGDB_PACKAGE" \
    --builder-status "$benchmark_status" \
    --skip-orgdb "$SKIP_ORGDB"

BUNDLE_NAME="create_custom_host_benchmark_bundle_$(safe_name "$MACHINE")_taxid${TAXID}_run${RUN_NUMBER}_$(date -u '+%Y%m%dT%H%M%SZ').tar.gz"
BUNDLE_PATH="$OUTPUT_ROOT/$BUNDLE_NAME"

tar -C "$OUTPUT_ROOT" -czf "$BUNDLE_PATH" "$(basename "$RUN_DIR")"

printf '\n============================================================\n'
printf '[RESULT] Benchmark status: %s\n' "$benchmark_status"
printf '[RESULT] Run directory:    %s\n' "$RUN_DIR"
printf '[RESULT] Summary:          %s\n' "$RUN_DIR/custom_host_summary.txt"
printf '[RESULT] Bundle:           %s\n' "$BUNDLE_PATH"
printf '============================================================\n'

exit "$benchmark_status"
