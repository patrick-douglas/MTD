#!/usr/bin/env bash
set -euo pipefail

offline_files_folder=""
offline_files_folder="${MTD_OFFLINE_FILES_FOLDER:-$offline_files_folder}"

if [[ -z "$offline_files_folder" ]]; then
    echo "[FAIL] MTD_OFFLINE_FILES_FOLDER is not configured." >&2
    exit 1
fi

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
helper="${MTD_MANIFEST_HELPER:-$script_dir/sync_ncbi_cache.py}"

[[ -f "$helper" ]] || {
    echo "[FAIL] NCBI cache synchronization helper not found: $helper" >&2
    exit 1
}

metadata_dir="$offline_files_folder/Kraken2DB_micro/library/fungi"

args=(
    assemblies
    --library fungi
    --summary-url "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt"
    --local-dir "$metadata_dir/all"
    --metadata-dir "$metadata_dir"
    --level "Complete Genome"
    --level "Chromosome"
    --min-count 100
    --require-complete
)

[[ "${FULL_GZIP_CHECK:-0}" == "1" ]] && args+=(--full-gzip-check)
[[ "${CHECK_ONLY:-0}" == "1" ]] && args+=(--check-only)

exec python3 "$helper" "${args[@]}"
