#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/mtd_hpc_common.sh"

usage() {
    cat <<'USAGE'
Usage:
  mtd_hpc_node_job.sh humann \
    --hpc-conf FILE --sample NAME --input FILE --output-dir DIR

  mtd_hpc_node_job.sh magicblast \
    --hpc-conf FILE --sample NAME --query FILE [--query-mate FILE] \
    --database PREFIX --output FILE --input-signature TEXT

  mtd_hpc_node_job.sh fastp \
    --hpc-conf FILE --sample NAME --layout se|pe \
    --read1 FILE --read1-size BYTES [--read2 FILE --read2-size BYTES] \
    --output-read1 FILE [--output-read2 FILE] \
    --html FILE --json FILE --min-length INT \
    --validate-paired-ids 0|1 --pe-max-attempts INT

  mtd_hpc_node_job.sh kraken2 \
    --hpc-conf FILE --sample NAME --layout se|pe --read1 FILE [--read2 FILE] \
    --database DIR --report FILE --kraken-output FILE \
    --classified-read1 FILE [--classified-read2 FILE] \
    --unclassified-read1 FILE [--unclassified-read2 FILE] \
    --confidence FLOAT --minimum-hit-groups INT --input-gzip 0|1

  mtd_hpc_node_job.sh bracken \
    --hpc-conf FILE --sample NAME --database DIR --report FILE \
    --read-length INT --threshold INT --phylum-output FILE \
    --genus-output FILE --species-output FILE
USAGE
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
LAYOUT=""
READ1=""
READ2=""
READ1_SIZE=""
READ2_SIZE=""
OUTPUT_READ1=""
OUTPUT_READ2=""
HTML_OUTPUT=""
JSON_OUTPUT=""
MIN_LENGTH=""
VALIDATE_PAIRED_IDS="1"
PE_MAX_ATTEMPTS="2"
REPORT=""
KRAKEN_OUTPUT=""
CLASSIFIED_READ1=""
CLASSIFIED_READ2=""
UNCLASSIFIED_READ1=""
UNCLASSIFIED_READ2=""
CONFIDENCE=""
MINIMUM_HIT_GROUPS=""
INPUT_GZIP="0"
READ_LENGTH=""
THRESHOLD=""
PHYLUM_OUTPUT=""
GENUS_OUTPUT=""
SPECIES_OUTPUT=""

LOCAL_SCRATCH=""
PARTIAL_OUTPUT=""

node_job_exit() {
    local status=$?

    trap - EXIT

    if (( status != 0 )) &&
       [[ -n "$PARTIAL_OUTPUT" ]] &&
       [[ "${MTD_HPC_STAGE_LOCAL:-0}" != "1" ]]
    then
        rm -f -- "$PARTIAL_OUTPUT"
    fi

    # Release cluster-wide mkdir transfer slots and the node-local stage-out
    # flock explicitly so cancellation/error paths do not leave stale state.
    mtd_hpc_stagein_group_end >/dev/null 2>&1 || true
    mtd_hpc_stageout_lock_release >/dev/null 2>&1 || true

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
        --layout) LAYOUT="$2"; shift 2 ;;
        --read1) READ1="$2"; shift 2 ;;
        --read2) READ2="$2"; shift 2 ;;
        --read1-size) READ1_SIZE="$2"; shift 2 ;;
        --read2-size) READ2_SIZE="$2"; shift 2 ;;
        --output-read1) OUTPUT_READ1="$2"; shift 2 ;;
        --output-read2) OUTPUT_READ2="$2"; shift 2 ;;
        --html) HTML_OUTPUT="$2"; shift 2 ;;
        --json) JSON_OUTPUT="$2"; shift 2 ;;
        --min-length) MIN_LENGTH="$2"; shift 2 ;;
        --validate-paired-ids) VALIDATE_PAIRED_IDS="$2"; shift 2 ;;
        --pe-max-attempts) PE_MAX_ATTEMPTS="$2"; shift 2 ;;
        --report) REPORT="$2"; shift 2 ;;
        --kraken-output) KRAKEN_OUTPUT="$2"; shift 2 ;;
        --classified-read1) CLASSIFIED_READ1="$2"; shift 2 ;;
        --classified-read2) CLASSIFIED_READ2="$2"; shift 2 ;;
        --unclassified-read1) UNCLASSIFIED_READ1="$2"; shift 2 ;;
        --unclassified-read2) UNCLASSIFIED_READ2="$2"; shift 2 ;;
        --confidence) CONFIDENCE="$2"; shift 2 ;;
        --minimum-hit-groups) MINIMUM_HIT_GROUPS="$2"; shift 2 ;;
        --input-gzip) INPUT_GZIP="$2"; shift 2 ;;
        --read-length) READ_LENGTH="$2"; shift 2 ;;
        --threshold) THRESHOLD="$2"; shift 2 ;;
        --phylum-output) PHYLUM_OUTPUT="$2"; shift 2 ;;
        --genus-output) GENUS_OUTPUT="$2"; shift 2 ;;
        --species-output) SPECIES_OUTPUT="$2"; shift 2 ;;
        --help) usage; exit 0 ;;
        *) mtd_hpc_error "Unknown option: $1"; usage; exit 2 ;;
    esac
done

[[ -n "$HPC_CONF" ]] || mtd_hpc_die "--hpc-conf is required"
mtd_hpc_load_config "$HPC_CONF" "${MTD_HPC_MTD_ROOT:-}"
[[ -n "${MTD_NODE_THREADS:-}" ]] || mtd_hpc_export_node_resources

[[ -x "$MTD_HPC_CONDA_BIN" ]] || \
    mtd_hpc_die "Node-local Conda executable not found: $MTD_HPC_CONDA_BIN"

mtd_hpc_conda_run() {
    local environment="$1"
    shift

    [[ -d "$environment/conda-meta" ]] || \
        mtd_hpc_die "Node-local environment not found: $environment" || return 1

    "$MTD_HPC_CONDA_BIN" run --no-capture-output --prefix "$environment" "$@"
}

mtd_hpc_validate_layout() {
    [[ "$LAYOUT" == "se" || "$LAYOUT" == "pe" ]] || \
        mtd_hpc_die "Invalid read layout: $LAYOUT"
}

mtd_hpc_validate_kraken_database() {
    local database="$1"
    mtd_hpc_require_dir "$database" "node-local Kraken2 database" || return 1
    mtd_hpc_require_file "$database/hash.k2d" "Kraken2 hash.k2d" || return 1
    mtd_hpc_require_file "$database/opts.k2d" "Kraken2 opts.k2d" || return 1
    mtd_hpc_require_file "$database/taxo.k2d" "Kraken2 taxo.k2d" || return 1
}

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

            mtd_hpc_stagein_group_begin
            mtd_hpc_stage_in_file \
                "$shared_input" \
                "$humann_input" \
                "HUMAnN input"
            mtd_hpc_stagein_group_end
        else
            mkdir -p -- "$shared_output_dir"
        fi

        metaphlan_options="--bowtie2db $MTD_HPC_HUMANN_DB_ROOT/metaphlan --index $MTD_HPC_METAPHLAN_INDEX --offline"

        mtd_hpc_conda_run "$MTD_HPC_ENV_DIR" env \
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

            mtd_hpc_require_file "$local_genefamilies" "local HUMAnN gene-family output"
            mtd_hpc_require_file "$local_pathabundance" "local HUMAnN pathway-abundance output"

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
        [[ -z "$QUERY_MATE" ]] || mtd_hpc_require_file "$QUERY_MATE" "Magic-BLAST mate query"

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

            mtd_hpc_stagein_group_begin
            mtd_hpc_stage_in_file "$shared_query" "$magicblast_query" "Magic-BLAST query"

            if [[ -n "$shared_query_mate" ]]; then
                magicblast_query_mate="$LOCAL_SCRATCH/input/query.R2.fastq"
                mtd_hpc_stage_in_file \
                    "$shared_query_mate" "$magicblast_query_mate" "Magic-BLAST mate query"
            fi
            mtd_hpc_stagein_group_end
        else
            mkdir -p -- "$(dirname -- "$shared_output")"
        fi

        PARTIAL_OUTPUT="${magicblast_output}.partial.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID:-0}"
        rm -f -- "$PARTIAL_OUTPUT"

        magicblast_threads="${SLURM_CPUS_PER_TASK:-$MTD_NODE_THREADS}"
        mtd_hpc_info "Magic-BLAST threads assigned to this task: $magicblast_threads"

        args=(
            -query "$magicblast_query"
            -db "$DATABASE"
            -infmt fastq
            -out "$PARTIAL_OUTPUT"
            -num_threads "$magicblast_threads"
        )
        [[ -z "$magicblast_query_mate" ]] || args+=( -query_mate "$magicblast_query_mate" )

        mtd_hpc_conda_run "$MTD_HPC_ENV_DIR" magicblast "${args[@]}"
        mtd_hpc_require_file "$PARTIAL_OUTPUT" "temporary Magic-BLAST SAM"

        mv -f -- "$PARTIAL_OUTPUT" "$magicblast_output"
        PARTIAL_OUTPUT=""

        if [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]]; then
            mtd_hpc_atomic_stage_out_file \
                "$magicblast_output" "$shared_output" "Magic-BLAST SAM"
        fi
        ;;

    fastp)
        [[ -n "$SAMPLE" && -n "$READ1" && -n "$OUTPUT_READ1" && \
           -n "$HTML_OUTPUT" && -n "$JSON_OUTPUT" && -n "$MIN_LENGTH" ]] || {
            usage
            exit 2
        }
        mtd_hpc_validate_id "$SAMPLE" "fastp sample"
        mtd_hpc_validate_layout
        [[ "$MIN_LENGTH" =~ ^[0-9]+$ ]] && (( MIN_LENGTH >= 1 )) || \
            mtd_hpc_die "--min-length must be an integer >= 1"
        [[ "$VALIDATE_PAIRED_IDS" == "0" || "$VALIDATE_PAIRED_IDS" == "1" ]] || \
            mtd_hpc_die "--validate-paired-ids must be 0 or 1"
        [[ "$PE_MAX_ATTEMPTS" =~ ^[0-9]+$ ]] && (( PE_MAX_ATTEMPTS >= 1 )) || \
            mtd_hpc_die "--pe-max-attempts must be an integer >= 1"
        [[ "$READ1_SIZE" =~ ^[0-9]+$ ]] && (( READ1_SIZE > 0 )) || \
            mtd_hpc_die "--read1-size must be an integer > 0 for fastp"
        if [[ "$LAYOUT" == "pe" ]]; then
            [[ -n "$READ2" && -n "$OUTPUT_READ2" ]] || mtd_hpc_die "Paired fastp task requires R2 input and output"
            [[ "$READ2_SIZE" =~ ^[0-9]+$ ]] && (( READ2_SIZE > 0 )) || \
                mtd_hpc_die "--read2-size must be an integer > 0 for paired fastp"
        else
            READ2_SIZE="-"
        fi

        if [[ "$MTD_HPC_STAGE_LOCAL" != "1" ]]; then
            mtd_hpc_require_file "$READ1" "fastp R1 input"
            [[ "$LAYOUT" == "se" ]] || mtd_hpc_require_file "$READ2" "fastp R2 input"
        fi

        mtd_hpc_require_local_scratch_capacity_bytes \
            "$MTD_HPC_FASTP_SCRATCH_MULTIPLIER" "fastp" "$READ1_SIZE" "$READ2_SIZE"

        shared_read1="$READ1"
        shared_read2="$READ2"
        shared_output_read1="$OUTPUT_READ1"
        shared_output_read2="$OUTPUT_READ2"
        shared_html="$HTML_OUTPUT"
        shared_json="$JSON_OUTPUT"

        fastp_read1="$shared_read1"
        fastp_read2="$shared_read2"
        fastp_output_read1="$shared_output_read1"
        fastp_output_read2="$shared_output_read2"
        fastp_html="$shared_html"
        fastp_json="$shared_json"

        if [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]]; then
            mtd_hpc_prepare_local_scratch fastp "$SAMPLE"
            LOCAL_SCRATCH="$MTD_HPC_TASK_SCRATCH"
            fastp_read1="$LOCAL_SCRATCH/input/read1.fastq"
            [[ "$shared_read1" == *.gz ]] && fastp_read1+=".gz"
            fastp_output_read1="$LOCAL_SCRATCH/output/read1.fq.gz"
            fastp_html="$LOCAL_SCRATCH/output/report.html"
            fastp_json="$LOCAL_SCRATCH/output/report.json"

            mtd_hpc_stagein_group_begin
            mtd_hpc_stage_in_file \
                "$shared_read1" "$fastp_read1" "fastp R1 input" "$READ1_SIZE"

            if [[ "$LAYOUT" == "pe" ]]; then
                fastp_read2="$LOCAL_SCRATCH/input/read2.fastq"
                [[ "$shared_read2" == *.gz ]] && fastp_read2+=".gz"
                fastp_output_read2="$LOCAL_SCRATCH/output/read2.fq.gz"
                mtd_hpc_stage_in_file \
                    "$shared_read2" "$fastp_read2" "fastp R2 input" "$READ2_SIZE"
            fi
            mtd_hpc_stagein_group_end
        else
            mkdir -p -- \
                "$(dirname -- "$shared_output_read1")" \
                "$(dirname -- "$shared_html")" \
                "$(dirname -- "$shared_json")"
            [[ "$LAYOUT" == "se" ]] || mkdir -p -- "$(dirname -- "$shared_output_read2")"
        fi

        fastp_common_args=(
            --trim_poly_x
            --qualified_quality_phred 15
            --unqualified_percent_limit 40
            --n_base_limit 5
            --cut_front
            --cut_front_window_size 1
            --cut_front_mean_quality 5
            --cut_tail
            --cut_tail_window_size 1
            --cut_tail_mean_quality 5
            --length_required "$MIN_LENGTH"
            --thread "$MTD_NODE_THREADS"
            --html "$fastp_html"
            --json "$fastp_json"
        )

        if [[ "$LAYOUT" == "pe" ]]; then
            fastp_attempt=1
            fastp_ok=0
            while (( fastp_attempt <= PE_MAX_ATTEMPTS )); do
                mtd_hpc_info "fastp paired attempt $fastp_attempt/$PE_MAX_ATTEMPTS for $SAMPLE"
                rm -f -- "$fastp_output_read1" "$fastp_output_read2" "$fastp_html" "$fastp_json"
                if mtd_hpc_conda_run "$MTD_HPC_FASTP_ENV_DIR" fastp \
                    "${fastp_common_args[@]}" \
                    -i "$fastp_read1" -I "$fastp_read2" \
                    -o "$fastp_output_read1" -O "$fastp_output_read2"
                then
                    mtd_hpc_require_file "$fastp_output_read1" "fastp R1 output"
                    mtd_hpc_require_file "$fastp_output_read2" "fastp R2 output"
                    if [[ "$VALIDATE_PAIRED_IDS" == "0" ]] || \
                       mtd_hpc_validate_fastq_pair_ids \
                           "$fastp_output_read1" "$fastp_output_read2" \
                           "fastp output for $SAMPLE"
                    then
                        fastp_ok=1
                        break
                    fi
                fi
                fastp_attempt=$((fastp_attempt + 1))
            done
            (( fastp_ok == 1 )) || mtd_hpc_die "fastp failed to produce synchronized paired output for $SAMPLE"
        else
            mtd_hpc_conda_run "$MTD_HPC_FASTP_ENV_DIR" fastp \
                "${fastp_common_args[@]}" \
                -i "$fastp_read1" -o "$fastp_output_read1"
        fi

        mtd_hpc_require_file "$fastp_output_read1" "fastp read output"
        mtd_hpc_require_file "$fastp_html" "fastp HTML report"
        mtd_hpc_require_file "$fastp_json" "fastp JSON report"
        [[ "$LAYOUT" == "se" ]] || mtd_hpc_require_file "$fastp_output_read2" "fastp R2 output"

        if [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]]; then
            mtd_hpc_atomic_stage_out_file "$fastp_output_read1" "$shared_output_read1" "fastp R1 output"
            if [[ "$LAYOUT" == "pe" ]]; then
                mtd_hpc_atomic_stage_out_file "$fastp_output_read2" "$shared_output_read2" "fastp R2 output"
            fi
            mtd_hpc_atomic_stage_out_file "$fastp_html" "$shared_html" "fastp HTML report"
            mtd_hpc_atomic_stage_out_file "$fastp_json" "$shared_json" "fastp JSON report"
        fi
        ;;

    kraken2)
        [[ -n "$SAMPLE" && -n "$READ1" && -n "$DATABASE" && -n "$REPORT" && \
           -n "$KRAKEN_OUTPUT" && -n "$CLASSIFIED_READ1" && -n "$UNCLASSIFIED_READ1" && \
           -n "$CONFIDENCE" && -n "$MINIMUM_HIT_GROUPS" ]] || {
            usage
            exit 2
        }
        mtd_hpc_validate_id "$SAMPLE" "Kraken2 sample"
        mtd_hpc_validate_layout
        [[ "$INPUT_GZIP" == "0" || "$INPUT_GZIP" == "1" ]] || \
            mtd_hpc_die "--input-gzip must be 0 or 1"
        [[ "$MINIMUM_HIT_GROUPS" =~ ^[0-9]+$ ]] && (( MINIMUM_HIT_GROUPS >= 1 )) || \
            mtd_hpc_die "--minimum-hit-groups must be an integer >= 1"
        mtd_hpc_validate_kraken_database "$DATABASE"
        mtd_hpc_require_file "$READ1" "Kraken2 R1 input"
        if [[ "$LAYOUT" == "pe" ]]; then
            [[ -n "$READ2" && -n "$CLASSIFIED_READ2" && -n "$UNCLASSIFIED_READ2" ]] || \
                mtd_hpc_die "Paired Kraken2 task requires R2 and paired output paths"
            mtd_hpc_require_file "$READ2" "Kraken2 R2 input"
        fi

        mtd_hpc_require_local_scratch_capacity \
            "$MTD_HPC_KRAKEN_SCRATCH_MULTIPLIER" "Kraken2" "$READ1" "$READ2"

        shared_read1="$READ1"
        shared_read2="$READ2"
        shared_report="$REPORT"
        shared_kraken_output="$KRAKEN_OUTPUT"
        shared_classified_read1="$CLASSIFIED_READ1"
        shared_classified_read2="$CLASSIFIED_READ2"
        shared_unclassified_read1="$UNCLASSIFIED_READ1"
        shared_unclassified_read2="$UNCLASSIFIED_READ2"

        kraken_read1="$shared_read1"
        kraken_read2="$shared_read2"
        kraken_report="$shared_report"
        kraken_output="$shared_kraken_output"
        kraken_classified_pattern="$shared_classified_read1"
        kraken_unclassified_pattern="$shared_unclassified_read1"

        if [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]]; then
            mtd_hpc_prepare_local_scratch kraken2 "$SAMPLE"
            LOCAL_SCRATCH="$MTD_HPC_TASK_SCRATCH"
            kraken_read1="$LOCAL_SCRATCH/input/read1.fq"
            [[ "$INPUT_GZIP" == "1" ]] && kraken_read1+=".gz"
            kraken_report="$LOCAL_SCRATCH/output/report.txt"
            kraken_output="$LOCAL_SCRATCH/output/classification.kraken"

            mtd_hpc_stagein_group_begin
            mtd_hpc_stage_in_file "$shared_read1" "$kraken_read1" "Kraken2 R1 input"

            if [[ "$LAYOUT" == "pe" ]]; then
                kraken_read2="$LOCAL_SCRATCH/input/read2.fq"
                [[ "$INPUT_GZIP" == "1" ]] && kraken_read2+=".gz"
                kraken_classified_pattern="$LOCAL_SCRATCH/output/classified#.fq"
                kraken_unclassified_pattern="$LOCAL_SCRATCH/output/unclassified#.fq"
                mtd_hpc_stage_in_file "$shared_read2" "$kraken_read2" "Kraken2 R2 input"
            else
                kraken_classified_pattern="$LOCAL_SCRATCH/output/classified.fq"
                kraken_unclassified_pattern="$LOCAL_SCRATCH/output/unclassified.fq"
            fi
            mtd_hpc_stagein_group_end
        else
            mkdir -p -- \
                "$(dirname -- "$shared_report")" \
                "$(dirname -- "$shared_kraken_output")" \
                "$(dirname -- "$shared_classified_read1")" \
                "$(dirname -- "$shared_unclassified_read1")"
            if [[ "$LAYOUT" == "pe" ]]; then
                kraken_classified_pattern="${shared_classified_read1%_1.fq}#.fq"
                kraken_unclassified_pattern="${shared_unclassified_read1%_1.fq}#.fq"
            fi
        fi

        kraken_args=(
            --db "$DATABASE"
            --use-names
            --confidence "$CONFIDENCE"
            --minimum-hit-groups "$MINIMUM_HIT_GROUPS"
            --report "$kraken_report"
            --output "$kraken_output"
            --threads "$MTD_NODE_THREADS"
        )
        [[ "$INPUT_GZIP" == "0" ]] || kraken_args+=(--gzip-compressed)

        if [[ "$LAYOUT" == "pe" ]]; then
            kraken_args+=(
                --paired
                --classified-out "$kraken_classified_pattern"
                --unclassified-out "$kraken_unclassified_pattern"
                "$kraken_read1" "$kraken_read2"
            )
        else
            kraken_args+=(
                --classified-out "$kraken_classified_pattern"
                --unclassified-out "$kraken_unclassified_pattern"
                "$kraken_read1"
            )
        fi

        mtd_hpc_conda_run "$MTD_HPC_KRAKEN2_ENV_DIR" kraken2 "${kraken_args[@]}"
        mtd_hpc_require_file "$kraken_report" "Kraken2 report"
        mtd_hpc_require_file "$kraken_output" "Kraken2 classification output"

        if [[ "$LAYOUT" == "pe" ]]; then
            local_classified1="${kraken_classified_pattern//#/_1}"
            local_classified2="${kraken_classified_pattern//#/_2}"
            local_unclassified1="${kraken_unclassified_pattern//#/_1}"
            local_unclassified2="${kraken_unclassified_pattern//#/_2}"
            for path in "$local_classified1" "$local_classified2" "$local_unclassified1" "$local_unclassified2"; do
                mtd_hpc_require_path_exists "$path" "Kraken2 category FASTQ"
            done
        else
            local_classified1="$kraken_classified_pattern"
            local_unclassified1="$kraken_unclassified_pattern"
            mtd_hpc_require_path_exists "$local_classified1" "Kraken2 classified FASTQ"
            mtd_hpc_require_path_exists "$local_unclassified1" "Kraken2 unclassified FASTQ"
        fi

        if [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]]; then
            mtd_hpc_atomic_stage_out_file "$kraken_report" "$shared_report" "Kraken2 report"
            mtd_hpc_atomic_stage_out_file "$kraken_output" "$shared_kraken_output" "Kraken2 output"
            mtd_hpc_atomic_stage_out_file "$local_classified1" "$shared_classified_read1" "Kraken2 classified R1" 1
            mtd_hpc_atomic_stage_out_file "$local_unclassified1" "$shared_unclassified_read1" "Kraken2 unclassified R1" 1
            if [[ "$LAYOUT" == "pe" ]]; then
                mtd_hpc_atomic_stage_out_file "$local_classified2" "$shared_classified_read2" "Kraken2 classified R2" 1
                mtd_hpc_atomic_stage_out_file "$local_unclassified2" "$shared_unclassified_read2" "Kraken2 unclassified R2" 1
            fi
        fi
        ;;

    bracken)
        [[ -n "$SAMPLE" && -n "$DATABASE" && -n "$REPORT" && \
           -n "$READ_LENGTH" && -n "$THRESHOLD" && -n "$PHYLUM_OUTPUT" && \
           -n "$GENUS_OUTPUT" && -n "$SPECIES_OUTPUT" ]] || {
            usage
            exit 2
        }
        mtd_hpc_validate_id "$SAMPLE" "Bracken sample"
        [[ "$READ_LENGTH" =~ ^[0-9]+$ ]] && (( READ_LENGTH >= 1 )) || \
            mtd_hpc_die "--read-length must be an integer >= 1"
        [[ "$THRESHOLD" =~ ^[0-9]+$ ]] || mtd_hpc_die "--threshold must be an integer >= 0"
        mtd_hpc_validate_kraken_database "$DATABASE"
        mtd_hpc_require_file "$DATABASE/database${READ_LENGTH}mers.kmer_distrib" \
            "Bracken read-length distribution"
        mtd_hpc_require_file "$REPORT" "Bracken Kraken2 report"

        shared_report="$REPORT"
        shared_phylum="$PHYLUM_OUTPUT"
        shared_genus="$GENUS_OUTPUT"
        shared_species="$SPECIES_OUTPUT"
        bracken_report="$shared_report"
        bracken_phylum="$shared_phylum"
        bracken_genus="$shared_genus"
        bracken_species="$shared_species"

        if [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]]; then
            mtd_hpc_prepare_local_scratch bracken "$SAMPLE"
            LOCAL_SCRATCH="$MTD_HPC_TASK_SCRATCH"
            bracken_report="$LOCAL_SCRATCH/input/report.txt"
            bracken_phylum="$LOCAL_SCRATCH/output/phylum.bracken"
            bracken_genus="$LOCAL_SCRATCH/output/genus.bracken"
            bracken_species="$LOCAL_SCRATCH/output/species.bracken"
            mtd_hpc_stagein_group_begin
            mtd_hpc_stage_in_file "$shared_report" "$bracken_report" "Bracken report"
            mtd_hpc_stagein_group_end
        else
            mkdir -p -- \
                "$(dirname -- "$shared_phylum")" \
                "$(dirname -- "$shared_genus")" \
                "$(dirname -- "$shared_species")"
        fi

        mtd_hpc_conda_run "$MTD_HPC_KRAKEN2_ENV_DIR" bracken \
            -d "$DATABASE" -i "$bracken_report" -o "$bracken_phylum" \
            -r "$READ_LENGTH" -l P -t "$THRESHOLD"
        mtd_hpc_conda_run "$MTD_HPC_KRAKEN2_ENV_DIR" bracken \
            -d "$DATABASE" -i "$bracken_report" -o "$bracken_genus" \
            -r "$READ_LENGTH" -l G -t "$THRESHOLD"
        mtd_hpc_conda_run "$MTD_HPC_KRAKEN2_ENV_DIR" bracken \
            -d "$DATABASE" -i "$bracken_report" -o "$bracken_species" \
            -r "$READ_LENGTH" -l S -t "$THRESHOLD"

        mtd_hpc_require_file "$bracken_phylum" "Bracken phylum output"
        mtd_hpc_require_file "$bracken_genus" "Bracken genus output"
        mtd_hpc_require_file "$bracken_species" "Bracken species output"

        if [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]]; then
            mtd_hpc_atomic_stage_out_file "$bracken_phylum" "$shared_phylum" "Bracken phylum output"
            mtd_hpc_atomic_stage_out_file "$bracken_genus" "$shared_genus" "Bracken genus output"
            mtd_hpc_atomic_stage_out_file "$bracken_species" "$shared_species" "Bracken species output"
        fi
        ;;

    *)
        mtd_hpc_error "Unknown node job mode: $MODE"
        usage
        exit 2
        ;;
esac
