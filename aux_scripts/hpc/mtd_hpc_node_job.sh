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

[[ -x "$MTD_HPC_CONDA_BIN" ]] || mtd_hpc_die "Node-local Conda executable not found: $MTD_HPC_CONDA_BIN"
[[ -d "$MTD_HPC_ENV_DIR" ]] || mtd_hpc_die "Node-local environment not found: $MTD_HPC_ENV_DIR"

conda_run=("$MTD_HPC_CONDA_BIN" run --no-capture-output --prefix "$MTD_HPC_ENV_DIR")

case "$MODE" in
    humann)
        [[ -n "$SAMPLE" && -n "$INPUT" && -n "$OUTPUT_DIR" ]] || {
            usage
            exit 2
        }
        mtd_hpc_validate_id "$SAMPLE" "sample name"
        mtd_hpc_require_file "$INPUT" "HUMAnN input"
        mtd_hpc_validate_humann_databases             "$MTD_HPC_HUMANN_DB_ROOT"             "$MTD_HPC_METAPHLAN_INDEX"
        mkdir -p -- "$OUTPUT_DIR"

        metaphlan_options="--bowtie2db $MTD_HPC_HUMANN_DB_ROOT/metaphlan --index $MTD_HPC_METAPHLAN_INDEX --offline"

        "${conda_run[@]}" env             -u PYTHONPATH             -u PYTHONHOME             PYTHONNOUSERSITE=1             humann \
            --input "$INPUT" \
            --output "$OUTPUT_DIR" \
            --threads "$MTD_NODE_THREADS" \
            --nucleotide-database "$MTD_HPC_HUMANN_DB_ROOT/chocophlan" \
            --protein-database "$MTD_HPC_HUMANN_DB_ROOT/uniref" \
            --metaphlan "$MTD_HPC_ENV_DIR/bin" \
            --metaphlan-options "$metaphlan_options" \
            --verbose
        ;;

    magicblast)
        [[ -n "$SAMPLE" && -n "$QUERY" && -n "$DATABASE" && -n "$OUTPUT" ]] || {
            usage
            exit 2
        }
        mtd_hpc_validate_id "$SAMPLE" "Magic-BLAST task sample"
        [[ -n "$INPUT_SIGNATURE" ]] || mtd_hpc_die "--input-signature is required for Magic-BLAST restart safety"
        mtd_hpc_info "Magic-BLAST input signature: $INPUT_SIGNATURE"
        mtd_hpc_require_file "$QUERY" "Magic-BLAST query"
        [[ -z "$QUERY_MATE" ]] || mtd_hpc_require_file "$QUERY_MATE" "Magic-BLAST mate query"

        # Magic-BLAST databases are prefixes, so validate at least one matching file.
        compgen -G "${DATABASE}.*" >/dev/null || \
            mtd_hpc_die "No Magic-BLAST database files match prefix: $DATABASE"

        mkdir -p -- "$(dirname -- "$OUTPUT")"
        temp_output="${OUTPUT}.partial.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID:-0}"
        rm -f -- "$temp_output"
        trap 'rm -f -- "$temp_output"' EXIT

        args=(
            -query "$QUERY"
            -db "$DATABASE"
            -infmt fastq
            -out "$temp_output"
            -num_threads "$MTD_NODE_THREADS"
        )
        [[ -z "$QUERY_MATE" ]] || args+=( -query_mate "$QUERY_MATE" )

        "${conda_run[@]}" magicblast "${args[@]}"
        mtd_hpc_require_file "$temp_output" "temporary Magic-BLAST SAM"
        mv -f -- "$temp_output" "$OUTPUT"
        trap - EXIT
        ;;

    *)
        mtd_hpc_error "Unknown node job mode: $MODE"
        usage
        exit 2
        ;;
esac
