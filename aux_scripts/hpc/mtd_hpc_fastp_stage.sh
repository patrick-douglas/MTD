#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/mtd_hpc_common.sh"

usage() {
    cat <<'USAGE'
Usage:
  mtd_hpc_fastp_stage.sh --hpc-conf FILE --input-manifest TSV \
    --min-length INT --validate-paired-ids 0|1 --pe-max-attempts INT \
    --work-dir DIR --mtd-root DIR

Input manifest columns (no header):
  sample<TAB>layout<TAB>read1<TAB>read2-or--<TAB>output-read1<TAB>output-read2-or--<TAB>html<TAB>json
USAGE
}

HPC_CONF=""
INPUT_MANIFEST=""
MIN_LENGTH=""
VALIDATE_PAIRED_IDS="1"
PE_MAX_ATTEMPTS="2"
WORK_DIR=""
MTD_ROOT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hpc-conf) HPC_CONF="$2"; shift 2 ;;
        --input-manifest) INPUT_MANIFEST="$2"; shift 2 ;;
        --min-length) MIN_LENGTH="$2"; shift 2 ;;
        --validate-paired-ids) VALIDATE_PAIRED_IDS="$2"; shift 2 ;;
        --pe-max-attempts) PE_MAX_ATTEMPTS="$2"; shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --mtd-root) MTD_ROOT="$2"; shift 2 ;;
        --help) usage; exit 0 ;;
        *) mtd_hpc_error "Unknown option: $1"; usage; exit 2 ;;
    esac
done

[[ -n "$HPC_CONF" && -n "$INPUT_MANIFEST" && -n "$MIN_LENGTH" && \
   -n "$WORK_DIR" && -n "$MTD_ROOT" ]] || { usage; exit 2; }

mtd_hpc_load_config "$HPC_CONF" "$MTD_ROOT"
mtd_hpc_require_file "$INPUT_MANIFEST" "fastp input manifest"
[[ "$MIN_LENGTH" =~ ^[0-9]+$ ]] && (( MIN_LENGTH >= 1 )) || \
    mtd_hpc_die "--min-length must be an integer >= 1"
[[ "$VALIDATE_PAIRED_IDS" == "0" || "$VALIDATE_PAIRED_IDS" == "1" ]] || \
    mtd_hpc_die "--validate-paired-ids must be 0 or 1"
[[ "$PE_MAX_ATTEMPTS" =~ ^[0-9]+$ ]] && (( PE_MAX_ATTEMPTS >= 1 )) || \
    mtd_hpc_die "--pe-max-attempts must be an integer >= 1"

mkdir -p -- "$WORK_DIR"
manifest="$WORK_DIR/fastp.tasks.tsv"
: > "$manifest"

while IFS=$'\t' read -r sample layout read1 read2 output_read1 output_read2 html json extra; do
    [[ -z "$sample" || "$sample" == \#* ]] && continue
    [[ -z "${extra:-}" ]] || mtd_hpc_die "Too many columns in fastp manifest for sample: $sample"
    mtd_hpc_validate_id "$sample" "sample name"
    [[ "$layout" == "se" || "$layout" == "pe" ]] || mtd_hpc_die "Invalid layout for $sample: $layout"
    mtd_hpc_require_file "$read1" "fastp R1 for $sample"
    [[ -n "$output_read1" && -n "$html" && -n "$json" ]] || \
        mtd_hpc_die "Missing fastp output path for sample: $sample"

    command=(
        "$SCRIPT_DIR/mtd_hpc_node_job.sh"
        fastp
        --hpc-conf "$HPC_CONF"
        --sample "$sample"
        --layout "$layout"
        --read1 "$read1"
        --output-read1 "$output_read1"
        --html "$html"
        --json "$json"
        --min-length "$MIN_LENGTH"
        --validate-paired-ids "$VALIDATE_PAIRED_IDS"
        --pe-max-attempts "$PE_MAX_ATTEMPTS"
    )
    expected="$output_read1"$'\n'"$html"$'\n'"$json"

    if [[ "$layout" == "pe" ]]; then
        [[ "$read2" != "-" && "$output_read2" != "-" ]] || \
            mtd_hpc_die "Paired fastp sample has no R2 path: $sample"
        mtd_hpc_require_file "$read2" "fastp R2 for $sample"
        command+=(--read2 "$read2" --output-read2 "$output_read2")
        expected+=$'\n'"$output_read2"
    fi

    printf -v command_string '%q ' "${command[@]}"
    mtd_hpc_task_line "fastp_${sample}" "$command_string" "$expected" >> "$manifest"
done < "$INPUT_MANIFEST"

"$SCRIPT_DIR/mtd_hpc_submit_array.sh" \
    --hpc-conf "$HPC_CONF" \
    --manifest "$manifest" \
    --stage fastp \
    --work-dir "$WORK_DIR/slurm" \
    --mtd-root "$MTD_ROOT"
