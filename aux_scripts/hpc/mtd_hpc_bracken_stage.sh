#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/mtd_hpc_common.sh"

usage() {
    cat <<'USAGE'
Usage:
  mtd_hpc_bracken_stage.sh --hpc-conf FILE --input-manifest TSV \
    --database DIR --read-length INT --threshold INT \
    --work-dir DIR --mtd-root DIR

Input manifest columns (no header):
  sample<TAB>kraken-report<TAB>phylum-output<TAB>genus-output<TAB>species-output
USAGE
}

HPC_CONF=""
INPUT_MANIFEST=""
DATABASE=""
READ_LENGTH=""
THRESHOLD=""
WORK_DIR=""
MTD_ROOT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hpc-conf) HPC_CONF="$2"; shift 2 ;;
        --input-manifest) INPUT_MANIFEST="$2"; shift 2 ;;
        --database) DATABASE="$2"; shift 2 ;;
        --read-length) READ_LENGTH="$2"; shift 2 ;;
        --threshold) THRESHOLD="$2"; shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --mtd-root) MTD_ROOT="$2"; shift 2 ;;
        --help) usage; exit 0 ;;
        *) mtd_hpc_error "Unknown option: $1"; usage; exit 2 ;;
    esac
done

[[ -n "$HPC_CONF" && -n "$INPUT_MANIFEST" && -n "$DATABASE" && \
   -n "$READ_LENGTH" && -n "$THRESHOLD" && -n "$WORK_DIR" && -n "$MTD_ROOT" ]] || {
    usage
    exit 2
}

mtd_hpc_load_config "$HPC_CONF" "$MTD_ROOT"
mtd_hpc_require_file "$INPUT_MANIFEST" "Bracken input manifest"
[[ "$READ_LENGTH" =~ ^[0-9]+$ ]] && (( READ_LENGTH >= 1 )) || \
    mtd_hpc_die "--read-length must be an integer >= 1"
[[ "$THRESHOLD" =~ ^[0-9]+$ ]] || mtd_hpc_die "--threshold must be an integer >= 0"

node_database="$(mtd_hpc_map_path_to_node "$DATABASE" "$MTD_ROOT")"
mkdir -p -- "$WORK_DIR"
manifest="$WORK_DIR/bracken.tasks.tsv"
: > "$manifest"

while IFS=$'\t' read -r sample report phylum genus species extra; do
    [[ -z "$sample" || "$sample" == \#* ]] && continue
    [[ -z "${extra:-}" ]] || mtd_hpc_die "Too many columns in Bracken manifest for sample: $sample"
    mtd_hpc_validate_id "$sample" "sample name"
    mtd_hpc_require_file "$report" "Bracken report for $sample"
    [[ -n "$phylum" && -n "$genus" && -n "$species" ]] || \
        mtd_hpc_die "Missing Bracken output path for sample: $sample"

    command=(
        "$SCRIPT_DIR/mtd_hpc_node_job.sh"
        bracken
        --hpc-conf "$HPC_CONF"
        --sample "$sample"
        --database "$node_database"
        --report "$report"
        --read-length "$READ_LENGTH"
        --threshold "$THRESHOLD"
        --phylum-output "$phylum"
        --genus-output "$genus"
        --species-output "$species"
    )
    printf -v command_string '%q ' "${command[@]}"
    expected="$phylum"$'\n'"$genus"$'\n'"$species"
    mtd_hpc_task_line "bracken_${sample}" "$command_string" "$expected" >> "$manifest"
done < "$INPUT_MANIFEST"

"$SCRIPT_DIR/mtd_hpc_submit_array.sh" \
    --hpc-conf "$HPC_CONF" \
    --manifest "$manifest" \
    --stage bracken \
    --work-dir "$WORK_DIR/slurm" \
    --mtd-root "$MTD_ROOT"
