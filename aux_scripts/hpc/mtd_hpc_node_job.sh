#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/mtd_hpc_common.sh"

usage() {
    cat <<'EOF'
Usage:
  mtd_hpc_node_job.sh humann \
    --hpc-conf FILE --sample NAME --input FILE --output-dir DIR

  mtd_hpc_node_job.sh magicblast \
    --hpc-conf FILE --sample NAME --query FILE [--query-mate FILE] \
    --database PREFIX --output FILE
EOF
}

[[ $# -ge 1 ]] || { usage; exit 2; }
MODE="$1"
shift

HPC_CONF=""
SAMPLE=""
INPUT=""
OUTPUT_DIR=""
QUERY=""
QUERY_MATE=""
DATABASE=""
OUTPUT=""
INPUT_SIGNATURE=""

LOCAL_SCRATCH=""
PARTIAL_OUTPUT=""

node_job_exit() {
    local status=$?

    trap - EXIT

    # Shared-storage partial files are disposable. Node-local partial files are
    # preserved together with the scratch directory when failure cleanup is off.
    if (( status != 0 )) &&
       [[ -n "$PARTIAL_OUTPUT" ]] &&
       [[ "$MTD_HPC_STAGE_LOCAL" != "1" ]]
    then
        rm -f -- "$PARTIAL_OUTPUT"
    fi

    mtd_hpc_cleanup_local_scratch "$LOCAL_SCRATCH" "$status" || true
    exit "$status"
}
trap node_job_exit EXIT

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hpc-conf) HPC_CONF="$2"; shift 2 ;;
        --sample) SAMPLE="$2"; shift 2 ;;
        --input) INPUT="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --query) QUERY="$2"; shift 2 ;;
        --query-mate) QUERY_MATE="$2"; shift 2 ;;
        --database) DATABASE="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --input-signature) INPUT_SIGNATURE="$2"; shift 2 ;;
        --help) usage; exit 0 ;;
        *) mtd_hpc_error "Unknown option: $1"; usage; exit 2 ;;
    esac
done

[[ -n "$HPC_CONF" ]] || mtd_hpc_die "--hpc-conf is required"
mtd_hpc_load_config "$HPC_CONF" "${MTD_HPC_MTD_ROOT:-}"
[[ -n "${MTD_NODE_THREADS:-}" ]] || mtd_hpc_export_node_resources

[[ -x "$MTD_HPC_CONDA_BIN" ]] || \
    mtd_hpc_die "Node-local Conda executable not found: $MTD_HPC_CONDA_BIN"
[[ -d "$MTD_HPC_ENV_DIR" ]] || \
    mtd_hpc_die "Node-local environment not found: $MTD_HPC_ENV_DIR"

conda_run=(
    "$MTD_HPC_CONDA_BIN"
    run
    --no-capture-output
    --prefix "$MTD_HPC_ENV_DIR"
)

case "$MODE" in
    humann)
        [[ -n "$SAMPLE" && -n "$INPUT" && -n "$OUTPUT_DIR" ]] || {
            usage
            exit 2
        }

        mtd_hpc_validate_id "$SAMPLE" "sample name"
        mtd_hpc_require_file "$INPUT" "HUMAnN input"
        mtd_hpc_validate_humann_databases \
            "$MTD_HPC_HUMANN_DB_ROOT" \
            "$MTD_HPC_METAPHLAN_INDEX"

        shared_input="$INPUT"
        shared_output_dir="$OUTPUT_DIR"
        humann_input="$shared_input"
        humann_output_dir="$shared_output_dir"

        if [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]]; then
            mtd_hpc_prepare_local_scratch humann "$SAMPLE"
            LOCAL_SCRATCH="$MTD_HPC_TASK_SCRATCH"

            humann_input="$LOCAL_SCRATCH/input/$(basename -- "$shared_input")"
            humann_output_dir="$LOCAL_SCRATCH/output"

            mtd_hpc_stage_in_file \
                "$shared_input" \
                "$humann_input" \
                "HUMAnN input"
        else
            mkdir -p -- "$shared_output_dir"
        fi

        metaphlan_options="--bowtie2db $MTD_HPC_HUMANN_DB_ROOT/metaphlan --index $MTD_HPC_METAPHLAN_INDEX --offline"

        "${conda_run[@]}" env \
            -u PYTHONPATH \
            -u PYTHONHOME \
            PYTHONNOUSERSITE=1 \
            humann \
            --input "$humann_input" \
            --output "$humann_output_dir" \
            --threads "$MTD_NODE_THREADS" \
            --nucleotide-database "$MTD_HPC_HUMANN_DB_ROOT/chocophlan" \
            --protein-database "$MTD_HPC_HUMANN_DB_ROOT/uniref" \
            --metaphlan "$MTD_HPC_ENV_DIR/bin" \
            --metaphlan-options "$metaphlan_options" \
            --verbose

        if [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]]; then
            humann_output_stem="$(basename -- "$humann_input")"
            humann_output_stem="${humann_output_stem%.gz}"
            humann_output_stem="${humann_output_stem%.fastq}"
            humann_output_stem="${humann_output_stem%.fq}"
            humann_output_stem="${humann_output_stem%.fasta}"
            humann_output_stem="${humann_output_stem%.fa}"

            local_genefamilies="$humann_output_dir/${humann_output_stem}_genefamilies.tsv"
            local_pathabundance="$humann_output_dir/${humann_output_stem}_pathabundance.tsv"

            mtd_hpc_require_file \
                "$local_genefamilies" \
                "local HUMAnN gene-family output"
            mtd_hpc_require_file \
                "$local_pathabundance" \
                "local HUMAnN pathway-abundance output"

            mtd_hpc_atomic_stage_out_file \
                "$local_genefamilies" \
                "$shared_output_dir/${SAMPLE}_genefamilies.tsv" \
                "HUMAnN gene-family output"

            mtd_hpc_atomic_stage_out_file \
                "$local_pathabundance" \
                "$shared_output_dir/${SAMPLE}_pathabundance.tsv" \
                "HUMAnN pathway-abundance output"
        fi
        ;;

    magicblast)
        [[ -n "$SAMPLE" && -n "$QUERY" && -n "$DATABASE" && -n "$OUTPUT" ]] || {
            usage
            exit 2
        }

        mtd_hpc_validate_id "$SAMPLE" "Magic-BLAST task sample"
        [[ -n "$INPUT_SIGNATURE" ]] || \
            mtd_hpc_die "--input-signature is required for Magic-BLAST restart safety"

        mtd_hpc_info "Magic-BLAST input signature: $INPUT_SIGNATURE"
        mtd_hpc_require_file "$QUERY" "Magic-BLAST query"
        [[ -z "$QUERY_MATE" ]] || \
            mtd_hpc_require_file "$QUERY_MATE" "Magic-BLAST mate query"

        # Magic-BLAST databases are prefixes, so validate at least one matching file.
        compgen -G "${DATABASE}.*" >/dev/null || \
            mtd_hpc_die "No Magic-BLAST database files match prefix: $DATABASE"

        shared_query="$QUERY"
        shared_query_mate="$QUERY_MATE"
        shared_output="$OUTPUT"

        magicblast_query="$shared_query"
        magicblast_query_mate="$shared_query_mate"
        magicblast_output="$shared_output"

        if [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]]; then
            mtd_hpc_prepare_local_scratch magicblast "$SAMPLE"
            LOCAL_SCRATCH="$MTD_HPC_TASK_SCRATCH"

            magicblast_query="$LOCAL_SCRATCH/input/query.R1.fastq"
            magicblast_output="$LOCAL_SCRATCH/output/result.sam"

            mtd_hpc_stage_in_file \
                "$shared_query" \
                "$magicblast_query" \
                "Magic-BLAST query"

            if [[ -n "$shared_query_mate" ]]; then
                magicblast_query_mate="$LOCAL_SCRATCH/input/query.R2.fastq"

                mtd_hpc_stage_in_file \
                    "$shared_query_mate" \
                    "$magicblast_query_mate" \
                    "Magic-BLAST mate query"
            fi
        else
            mkdir -p -- "$(dirname -- "$shared_output")"
        fi

        PARTIAL_OUTPUT="${magicblast_output}.partial.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID:-0}"
        rm -f -- "$PARTIAL_OUTPUT"

        # Packed Magic-BLAST chunk jobs normally receive one CPU each.
        # Fall back to the detected node CPUs when SLURM_CPUS_PER_TASK is absent.
        magicblast_threads="${SLURM_CPUS_PER_TASK:-$MTD_NODE_THREADS}"
        mtd_hpc_info "Magic-BLAST threads assigned to this task: $magicblast_threads"

        args=(
            -query "$magicblast_query"
            -db "$DATABASE"
            -infmt fastq
            -out "$PARTIAL_OUTPUT"
            -num_threads "$magicblast_threads"
        )
        [[ -z "$magicblast_query_mate" ]] || \
            args+=( -query_mate "$magicblast_query_mate" )

        "${conda_run[@]}" magicblast "${args[@]}"

        mtd_hpc_require_file \
            "$PARTIAL_OUTPUT" \
            "temporary Magic-BLAST SAM"

        mv -f -- "$PARTIAL_OUTPUT" "$magicblast_output"
        PARTIAL_OUTPUT=""

        if [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]]; then
            mtd_hpc_atomic_stage_out_file \
                "$magicblast_output" \
                "$shared_output" \
                "Magic-BLAST SAM"
        fi
        ;;

    *)
        mtd_hpc_error "Unknown node job mode: $MODE"
        usage
        exit 2
        ;;
esac
