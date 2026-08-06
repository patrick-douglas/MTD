#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/mtd_hpc_common.sh"

usage() {
    cat <<'USAGE'
Usage:
  mtd_hpc_kraken_stage.sh --hpc-conf FILE --stage NAME --input-manifest TSV \
    --database DIR --confidence FLOAT --minimum-hit-groups INT \
    --input-gzip 0|1 --work-dir DIR --mtd-root DIR

Supported stage names:
  kraken_host
  kraken_micro_raw
  kraken_micro_final

Input manifest columns (no header):
  sample<TAB>layout<TAB>read1<TAB>read2-or--<TAB>report<TAB>kraken-output<TAB>classified-r1<TAB>classified-r2-or--<TAB>unclassified-r1<TAB>unclassified-r2-or--
USAGE
}

HPC_CONF=""
STAGE=""
INPUT_MANIFEST=""
DATABASE=""
CONFIDENCE=""
MINIMUM_HIT_GROUPS=""
INPUT_GZIP="0"
WORK_DIR=""
MTD_ROOT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hpc-conf) HPC_CONF="$2"; shift 2 ;;
        --stage) STAGE="$2"; shift 2 ;;
        --input-manifest) INPUT_MANIFEST="$2"; shift 2 ;;
        --database) DATABASE="$2"; shift 2 ;;
        --confidence) CONFIDENCE="$2"; shift 2 ;;
        --minimum-hit-groups) MINIMUM_HIT_GROUPS="$2"; shift 2 ;;
        --input-gzip) INPUT_GZIP="$2"; shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --mtd-root) MTD_ROOT="$2"; shift 2 ;;
        --help) usage; exit 0 ;;
        *) mtd_hpc_error "Unknown option: $1"; usage; exit 2 ;;
    esac
done

[[ -n "$HPC_CONF" && -n "$STAGE" && -n "$INPUT_MANIFEST" && \
   -n "$DATABASE" && -n "$CONFIDENCE" && -n "$MINIMUM_HIT_GROUPS" && \
   -n "$WORK_DIR" && -n "$MTD_ROOT" ]] || { usage; exit 2; }

case "$STAGE" in
    kraken_host|kraken_micro_raw|kraken_micro_final) ;;
    *) mtd_hpc_die "Unsupported Kraken2 stage: $STAGE" ;;
esac

mtd_hpc_load_config "$HPC_CONF" "$MTD_ROOT"
mtd_hpc_require_file "$INPUT_MANIFEST" "Kraken2 input manifest"
[[ "$MINIMUM_HIT_GROUPS" =~ ^[0-9]+$ ]] && (( MINIMUM_HIT_GROUPS >= 1 )) || \
    mtd_hpc_die "--minimum-hit-groups must be an integer >= 1"
[[ "$INPUT_GZIP" == "0" || "$INPUT_GZIP" == "1" ]] || \
    mtd_hpc_die "--input-gzip must be 0 or 1"

node_database="$(mtd_hpc_map_path_to_node "$DATABASE" "$MTD_ROOT")"
mkdir -p -- "$WORK_DIR"
manifest="$WORK_DIR/${STAGE}.tasks.tsv"
: > "$manifest"

while IFS=$'\t' read -r sample layout read1 read2 report kraken_output classified1 classified2 unclassified1 unclassified2 extra; do
    [[ -z "$sample" || "$sample" == \#* ]] && continue
    [[ -z "${extra:-}" ]] || mtd_hpc_die "Too many columns in Kraken2 manifest for sample: $sample"
    mtd_hpc_validate_id "$sample" "sample name"
    [[ "$layout" == "se" || "$layout" == "pe" ]] || mtd_hpc_die "Invalid layout for $sample: $layout"
    mtd_hpc_require_file "$read1" "Kraken2 R1 for $sample"
    [[ -n "$report" && -n "$kraken_output" && -n "$classified1" && -n "$unclassified1" ]] || \
        mtd_hpc_die "Missing Kraken2 output path for sample: $sample"

    command=(
        "$SCRIPT_DIR/mtd_hpc_node_job.sh"
        kraken2
        --hpc-conf "$HPC_CONF"
        --sample "$sample"
        --layout "$layout"
        --read1 "$read1"
        --database "$node_database"
        --report "$report"
        --kraken-output "$kraken_output"
        --classified-read1 "$classified1"
        --unclassified-read1 "$unclassified1"
        --confidence "$CONFIDENCE"
        --minimum-hit-groups "$MINIMUM_HIT_GROUPS"
        --input-gzip "$INPUT_GZIP"
    )

    expected="$report"$'\n'"$kraken_output"$'\n'"exists:$classified1"$'\n'"exists:$unclassified1"
    if [[ "$layout" == "pe" ]]; then
        [[ "$read2" != "-" && "$classified2" != "-" && "$unclassified2" != "-" ]] || \
            mtd_hpc_die "Paired Kraken2 sample has incomplete R2 paths: $sample"
        mtd_hpc_require_file "$read2" "Kraken2 R2 for $sample"
        command+=(
            --read2 "$read2"
            --classified-read2 "$classified2"
            --unclassified-read2 "$unclassified2"
        )
        expected+=$'\n'"exists:$classified2"$'\n'"exists:$unclassified2"
    fi

    printf -v command_string '%q ' "${command[@]}"
    mtd_hpc_task_line "${STAGE}_${sample}" "$command_string" "$expected" >> "$manifest"
done < "$INPUT_MANIFEST"

"$SCRIPT_DIR/mtd_hpc_submit_array.sh" \
    --hpc-conf "$HPC_CONF" \
    --manifest "$manifest" \
    --stage "$STAGE" \
    --work-dir "$WORK_DIR/slurm" \
    --mtd-root "$MTD_ROOT"
