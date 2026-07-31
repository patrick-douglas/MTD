#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/mtd_hpc_common.sh"

usage() {
    cat <<'EOF'
Usage:
  mtd_hpc_magicblast_stage.sh --hpc-conf FILE --input-manifest TSV \
    --main-database-prefix PREFIX --work-dir DIR --mtd-root DIR

Input manifest columns (no header):
  sample<TAB>layout<TAB>read1<TAB>read2-or--<TAB>final-sam
EOF
}

HPC_CONF=""
INPUT_MANIFEST=""
MAIN_DATABASE=""
WORK_DIR=""
MTD_ROOT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hpc-conf) HPC_CONF="$2"; shift 2 ;;
        --input-manifest) INPUT_MANIFEST="$2"; shift 2 ;;
        --main-database-prefix) MAIN_DATABASE="$2"; shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --mtd-root) MTD_ROOT="$2"; shift 2 ;;
        --help) usage; exit 0 ;;
        *) mtd_hpc_error "Unknown option: $1"; usage; exit 2 ;;
    esac
done

[[ -n "$HPC_CONF" && -n "$INPUT_MANIFEST" && -n "$MAIN_DATABASE" && -n "$WORK_DIR" && -n "$MTD_ROOT" ]] || {
    usage
    exit 2
}

mtd_hpc_load_config "$HPC_CONF" "$MTD_ROOT"
mtd_hpc_require_file "$INPUT_MANIFEST" "Magic-BLAST input manifest"
mkdir -p -- "$WORK_DIR"

node_database="$(mtd_hpc_map_path_to_node "$MAIN_DATABASE" "$MTD_ROOT")"
manifest="$WORK_DIR/magicblast.tasks.tsv"
: > "$manifest"

sample_outputs_dir="$WORK_DIR/sample_outputs"
chunks_root="$WORK_DIR/chunks"
mkdir -p -- "$sample_outputs_dir" "$chunks_root"

while IFS=$'\t' read -r sample layout read1 read2 final_sam; do
    [[ -z "$sample" || "$sample" == \#* ]] && continue
    mtd_hpc_validate_id "$sample" "sample name"
    [[ "$layout" == "se" || "$layout" == "pe" ]] || mtd_hpc_die "Invalid layout for $sample: $layout"
    mtd_hpc_require_file "$read1" "Magic-BLAST R1 for $sample"
    if [[ "$layout" == "pe" ]]; then
        [[ "$read2" != "-" ]] || mtd_hpc_die "Paired sample $sample has no R2"
        mtd_hpc_require_file "$read2" "Magic-BLAST R2 for $sample"
    else
        read2="-"
    fi
    [[ -n "$final_sam" ]] || mtd_hpc_die "Final SAM path is empty for sample: $sample"

    sample_dir="$chunks_root/$sample"
    chunk_manifest="$sample_dir/chunks.tsv"
    split_state="$sample_dir/split.state"
    mkdir -p -- "$sample_dir"

    read1_signature="$(readlink -f -- "$read1")|$(stat -c '%s|%Y' -- "$read1")"
    if [[ "$layout" == "pe" ]]; then
        read2_signature="$(readlink -f -- "$read2")|$(stat -c '%s|%Y' -- "$read2")"
    else
        read2_signature="-"
    fi
    split_signature="layout=$layout;chunk_reads=$MTD_HPC_MAGICBLAST_CHUNK_READS;r1=$read1_signature;r2=$read2_signature"

    reuse_split=0
    if [[ -s "$split_state" && -s "$chunk_manifest" ]] && \
       grep -Fxq "$split_signature" "$split_state"; then
        reuse_split=1
        while IFS=$'\t' read -r _chunk_id existing_r1 existing_r2 _count; do
            [[ -s "$existing_r1" ]] || { reuse_split=0; break; }
            [[ "$existing_r2" == "-" || -s "$existing_r2" ]] || { reuse_split=0; break; }
        done < "$chunk_manifest"
    fi

    if (( reuse_split )); then
        mtd_hpc_info "Reusing existing FASTQ chunks for sample: $sample"
    else
        rm -f -- "$sample_dir"/${sample}.chunk*.R1.fq \
                  "$sample_dir"/${sample}.chunk*.R2.fq \
                  "$chunk_manifest" "$split_state"
        if (( MTD_HPC_MAGICBLAST_CHUNK_READS > 0 )); then
            split_args=(
                "$SCRIPT_DIR/mtd_split_fastq.py"
                --read1 "$read1"
                --output-dir "$sample_dir"
                --sample "$sample"
                --records-per-chunk "$MTD_HPC_MAGICBLAST_CHUNK_READS"
                --manifest "$chunk_manifest"
            )
            [[ "$layout" == "pe" ]] && split_args+=(--read2 "$read2")
            python3 "${split_args[@]}"
        else
            printf '000001\t%s\t%s\t0\n' "$read1" "$read2" > "$chunk_manifest"
        fi
        printf '%s\n' "$split_signature" > "${split_state}.tmp"
        mv -f -- "${split_state}.tmp" "$split_state"
    fi

    sam_list="$sample_outputs_dir/${sample}.sam_chunks.txt"
    : > "$sam_list"

    while IFS=$'\t' read -r chunk_id chunk_r1 chunk_r2 _record_count; do
        [[ -z "$chunk_id" ]] && continue
        chunk_sam="$sample_dir/${sample}.chunk${chunk_id}.sam"
        printf '%s\n' "$chunk_sam" >> "$sam_list"

        command=(
            "$SCRIPT_DIR/mtd_hpc_node_job.sh"
            magicblast
            --hpc-conf "$HPC_CONF"
            --sample "${sample}.chunk${chunk_id}"
            --query "$chunk_r1"
            --database "$node_database"
            --output "$chunk_sam"
            --input-signature "$split_signature"
        )
        [[ "$chunk_r2" == "-" ]] || command+=(--query-mate "$chunk_r2")
        printf -v command_string '%q ' "${command[@]}"
        mtd_hpc_task_line "magicblast_${sample}_${chunk_id}" "$command_string" "$chunk_sam" >> "$manifest"
    done < "$chunk_manifest"

done < "$INPUT_MANIFEST"

"$SCRIPT_DIR/mtd_hpc_submit_array.sh" \
    --hpc-conf "$HPC_CONF" \
    --manifest "$manifest" \
    --stage magicblast \
    --work-dir "$WORK_DIR/slurm" \
    --mtd-root "$MTD_ROOT"

while IFS=$'\t' read -r sample _layout _read1 _read2 final_sam; do
    [[ -z "$sample" || "$sample" == \#* ]] && continue
    sam_list="$sample_outputs_dir/${sample}.sam_chunks.txt"
    python3 "$SCRIPT_DIR/mtd_merge_sam_chunks.py" \
        --chunk-list "$sam_list" \
        --output "$final_sam"
    mtd_hpc_require_file "$final_sam" "merged Magic-BLAST SAM for $sample"
done < "$INPUT_MANIFEST"
