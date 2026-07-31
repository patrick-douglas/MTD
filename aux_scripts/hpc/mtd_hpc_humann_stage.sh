#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/mtd_hpc_common.sh"

usage() {
    cat <<'EOF'
Usage:
  mtd_hpc_humann_stage.sh --hpc-conf FILE --samples FILE \
    --input-dir DIR --output-dir DIR --work-dir DIR --mtd-root DIR
EOF
}

HPC_CONF=""
SAMPLES_FILE=""
INPUT_DIR=""
OUTPUT_DIR=""
WORK_DIR=""
MTD_ROOT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hpc-conf) HPC_CONF="$2"; shift 2 ;;
        --samples) SAMPLES_FILE="$2"; shift 2 ;;
        --input-dir) INPUT_DIR="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --mtd-root) MTD_ROOT="$2"; shift 2 ;;
        --help) usage; exit 0 ;;
        *) mtd_hpc_error "Unknown option: $1"; usage; exit 2 ;;
    esac
done

[[ -n "$HPC_CONF" && -n "$SAMPLES_FILE" && -n "$INPUT_DIR" && -n "$OUTPUT_DIR" && -n "$WORK_DIR" && -n "$MTD_ROOT" ]] || {
    usage
    exit 2
}

mtd_hpc_load_config "$HPC_CONF" "$MTD_ROOT"
mtd_hpc_require_file "$SAMPLES_FILE" "sample list"
mtd_hpc_require_dir "$INPUT_DIR" "HUMAnN input directory"
mkdir -p -- "$OUTPUT_DIR" "$WORK_DIR"

manifest="$WORK_DIR/humann.tasks.tsv"
: > "$manifest"

while IFS= read -r sample; do
    [[ -z "$sample" || "$sample" == \#* ]] && continue
    mtd_hpc_validate_id "$sample" "sample name"
    input="$INPUT_DIR/${sample}.fq"
    mtd_hpc_require_file "$input" "HUMAnN input for $sample"

    command=(
        "$SCRIPT_DIR/mtd_hpc_node_job.sh"
        humann
        --hpc-conf "$HPC_CONF"
        --sample "$sample"
        --input "$input"
        --output-dir "$OUTPUT_DIR"
    )
    printf -v command_string '%q ' "${command[@]}"
    expected="${OUTPUT_DIR}/${sample}_genefamilies.tsv"$'\n'"${OUTPUT_DIR}/${sample}_pathabundance.tsv"
    mtd_hpc_task_line "humann_${sample}" "$command_string" "$expected" >> "$manifest"
done < "$SAMPLES_FILE"

"$SCRIPT_DIR/mtd_hpc_submit_array.sh" \
    --hpc-conf "$HPC_CONF" \
    --manifest "$manifest" \
    --stage humann \
    --work-dir "$WORK_DIR" \
    --mtd-root "$MTD_ROOT"
