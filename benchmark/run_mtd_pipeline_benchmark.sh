#!/usr/bin/env bash
set -Eeuo pipefail

PROGRAM_NAME="$(basename "$0")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

PIPELINE_SCRIPT="$REPO_DIR/MTD_explorer.sh"
RESOURCE_WRAPPER="$SCRIPT_DIR/MTD_benchmark_install.sh"
REPORTER="$SCRIPT_DIR/MTD_pipeline_benchmark_report.py"

MACHINE_NAME=""
DATASET_LABEL=""
SAMPLESHEET=""
PIPELINE_OUTPUT=""
HOST_ID=""
THREADS="$(nproc)"
# The official benchmark does not extract reads by default. When
# extraction is explicitly enabled, limit it to the top 10 taxa unless
# the user deliberately requests another value.
TOP_N=10
TOP_N_EXPLICIT=0
RUN_NUMBER=1
INTERVAL=5
BENCHMARK_ROOT="$HOME/MTD_pipeline_benchmarks"
READ_LAYOUT="auto"
ANALYSIS_MODE="auto"
USE_BLAST=1
EXTRACT_READS=0
ASSUME_YES=0
EXTRA_MTD_ARGS=()

usage() {
    cat <<USAGE
Usage:
  bash $PROGRAM_NAME [options] [-- extra MTD_explorer.sh arguments]

Required:
  -n, --machine NAME          Machine label, e.g. master or s5
  -i, --input FILE            Path to samplesheet.csv
  -o, --output DIR            Clean MTD Explorer output directory
  -h, --hostid TAXID          Host NCBI TaxID

Recommended:
  -d, --dataset LABEL         Dataset label; default: samplesheet parent directory
      --threads INT           Fixed CPU thread count; default: nproc
      --run-number INT        Repetition number; default: 1

Pipeline options:
      --extract               Enable detected microbiome read extraction
                              for an optional extraction benchmark
      --top-n INT             Taxa extracted from absolute ranking when
                              --extract is enabled; default: 10
                              0 = all detected taxa (potentially very slow)
      --no-extract            Disable read extraction; benchmark default
                              and retained as an explicit compatibility option
      --read-layout MODE      auto, se, or pe; default: auto
      --analysis-mode MODE    auto, comparison, or exploratory; default: auto
      --hisat2                Use HISAT2 instead of Magic-BLAST

Benchmark options:
      --interval SECONDS      Resource sampling interval; default: 5
      --benchmark-root DIR    Benchmark output root
                              default: \$HOME/MTD_pipeline_benchmarks
  -y, --yes                   Delete an existing pipeline output without prompting
      --help                  Show this help

Example for the standard Biomphalaria glabrata benchmark
(read extraction disabled):
  bash benchmark/$PROGRAM_NAME \\
    --machine master \\
    --dataset Bglabrata_PRJNA1306560 \\
    --input /path/to/B.glabrata_fastq/samplesheet.csv \\
    --output "\$HOME/test_MTD_explorer_B.glabrata" \\
    --hostid 6526 \\
    --threads 20 \\
    --run-number 1

Optional functional extraction benchmark using the top 10 taxa:
  bash benchmark/$PROGRAM_NAME \\
    --machine master \\
    --dataset Bglabrata_PRJNA1306560_extract_top10 \\
    --input /path/to/B.glabrata_fastq/samplesheet.csv \\
    --output "\$HOME/test_MTD_explorer_B.glabrata.extract_top10" \\
    --hostid 6526 \\
    --threads 20 \\
    --extract \\
    --top-n 10 \\
    --run-number 1

Stress-test extraction of every detected taxon, which can take many
hours or days and consume substantial storage:
  --extract --top-n 0

Extra MTD Explorer arguments may be added after --, for example:
  -- --metadata /path/to/metadata.csv --pdm spearman
USAGE
}

die() {
    printf '[STOP] %s\n' "$*" >&2
    exit 1
}

warn() {
    printf '[WARN] %s\n' "$*" >&2
}

expand_home() {
    local value="$1"
    case "$value" in
        "~") printf '%s\n' "$HOME" ;;
        "~/"*) printf '%s/%s\n' "$HOME" "${value#~/}" ;;
        *) printf '%s\n' "$value" ;;
    esac
}

safe_label() {
    local value="$1"
    value="$(printf '%s' "$value" | tr -cs 'A-Za-z0-9._-' '_')"
    value="${value##_}"
    value="${value%%_}"
    printf '%s\n' "$value"
}

is_same_or_ancestor() {
    local ancestor="${1%/}"
    local child="${2%/}"
    [[ "$child" == "$ancestor" || "$child" == "$ancestor/"* ]]
}

while (($# > 0)); do
    case "$1" in
        -n|--machine)
            (($# >= 2)) || die "$1 requires a value."
            MACHINE_NAME="$2"
            shift 2
            ;;
        -d|--dataset)
            (($# >= 2)) || die "$1 requires a value."
            DATASET_LABEL="$2"
            shift 2
            ;;
        -i|--input)
            (($# >= 2)) || die "$1 requires a value."
            SAMPLESHEET="$2"
            shift 2
            ;;
        -o|--output)
            (($# >= 2)) || die "$1 requires a value."
            PIPELINE_OUTPUT="$2"
            shift 2
            ;;
        -h|--hostid)
            (($# >= 2)) || die "$1 requires a value."
            HOST_ID="$2"
            shift 2
            ;;
        --threads)
            (($# >= 2)) || die "$1 requires a value."
            THREADS="$2"
            shift 2
            ;;
        --extract)
            EXTRACT_READS=1
            shift
            ;;
        --top-n)
            (($# >= 2)) || die "$1 requires a value."
            TOP_N="$2"
            TOP_N_EXPLICIT=1
            shift 2
            ;;
        --run-number)
            (($# >= 2)) || die "$1 requires a value."
            RUN_NUMBER="$2"
            shift 2
            ;;
        --interval)
            (($# >= 2)) || die "$1 requires a value."
            INTERVAL="$2"
            shift 2
            ;;
        --benchmark-root)
            (($# >= 2)) || die "$1 requires a value."
            BENCHMARK_ROOT="$2"
            shift 2
            ;;
        --read-layout)
            (($# >= 2)) || die "$1 requires a value."
            READ_LAYOUT="${2,,}"
            shift 2
            ;;
        --analysis-mode)
            (($# >= 2)) || die "$1 requires a value."
            ANALYSIS_MODE="${2,,}"
            shift 2
            ;;
        --hisat2)
            USE_BLAST=0
            shift
            ;;
        --no-extract)
            EXTRACT_READS=0
            shift
            ;;
        -y|--yes)
            ASSUME_YES=1
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        --)
            shift
            EXTRA_MTD_ARGS=("$@")
            break
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

# Extraction arguments are managed by this benchmark wrapper so that
# benchmark labels and reports remain consistent with the executed command.
for extra_arg in "${EXTRA_MTD_ARGS[@]}"; do
    case "$extra_arg" in
        --extract-microbiome-reads|--extract-microbiome-reads-top-n|--extract-microbiome-reads-top-n=*)
            die "Do not pass extraction options after --. Use benchmark options --extract and --top-n instead."
            ;;
    esac
done

[[ -n "$MACHINE_NAME" ]] || die "Machine name is required."
[[ -n "$SAMPLESHEET" ]] || die "Samplesheet is required."
[[ -n "$PIPELINE_OUTPUT" ]] || die "Pipeline output directory is required."
[[ -n "$HOST_ID" ]] || die "Host TaxID is required."

[[ "$HOST_ID" =~ ^[0-9]+$ ]] || die "--hostid must be a positive integer."
(( HOST_ID > 0 )) || die "--hostid must be greater than zero."
[[ "$THREADS" =~ ^[0-9]+$ ]] || die "--threads must be a positive integer."
(( THREADS > 0 )) || die "--threads must be greater than zero."
[[ "$TOP_N" =~ ^[0-9]+$ ]] || die "--top-n must be an integer >= 0."

if (( TOP_N_EXPLICIT && ! EXTRACT_READS )); then
    die "--top-n only applies when --extract is enabled. The standard benchmark intentionally disables read extraction."
fi

if (( EXTRACT_READS && TOP_N == 0 )); then
    warn "Read extraction was enabled for all detected taxa (--extract --top-n 0)."
    warn "This stress-test mode may repeatedly scan large Kraken2/FASTQ files, take many hours or days, create hundreds of taxa-specific outputs, and consume substantial disk space."
    warn "For a functional extraction benchmark, prefer --extract --top-n 10."
fi

[[ "$RUN_NUMBER" =~ ^[0-9]+$ ]] || die "--run-number must be a positive integer."
(( RUN_NUMBER > 0 )) || die "--run-number must be greater than zero."
[[ "$INTERVAL" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "--interval must be numeric."
awk -v x="$INTERVAL" 'BEGIN { exit !(x > 0) }' || die "--interval must be greater than zero."

case "$READ_LAYOUT" in
    auto|se|pe) ;;
    *) die "--read-layout must be auto, se, or pe." ;;
esac

case "$ANALYSIS_MODE" in
    auto|comparison|exploratory) ;;
    *) die "--analysis-mode must be auto, comparison, or exploratory." ;;
esac

SAMPLESHEET="$(expand_home "$SAMPLESHEET")"
PIPELINE_OUTPUT="$(expand_home "$PIPELINE_OUTPUT")"
BENCHMARK_ROOT="$(expand_home "$BENCHMARK_ROOT")"

[[ -s "$SAMPLESHEET" ]] || die "Samplesheet not found or empty: $SAMPLESHEET"
SAMPLESHEET="$(readlink -f -- "$SAMPLESHEET")"
PIPELINE_OUTPUT="$(readlink -m -- "$PIPELINE_OUTPUT")"
BENCHMARK_ROOT="$(readlink -m -- "$BENCHMARK_ROOT")"

if [[ -z "$DATASET_LABEL" ]]; then
    DATASET_LABEL="$(basename "$(dirname "$SAMPLESHEET")")"
fi

MACHINE_SAFE="$(safe_label "$MACHINE_NAME")"
DATASET_SAFE="$(safe_label "$DATASET_LABEL")"
[[ -n "$MACHINE_SAFE" ]] || die "Machine label does not contain usable characters."
[[ -n "$DATASET_SAFE" ]] || die "Dataset label does not contain usable characters."

for protected in "/" "$HOME" "$REPO_DIR" "$(dirname "$SAMPLESHEET")" "$BENCHMARK_ROOT"; do
    protected="$(readlink -m -- "$protected")"
    if [[ "$PIPELINE_OUTPUT" == "$protected" ]]; then
        die "Unsafe pipeline output path: $PIPELINE_OUTPUT"
    fi
done

if is_same_or_ancestor "$PIPELINE_OUTPUT" "$SAMPLESHEET"; then
    die "Pipeline output cannot contain the samplesheet."
fi
if is_same_or_ancestor "$PIPELINE_OUTPUT" "$REPO_DIR"; then
    die "Pipeline output cannot contain the repository."
fi
if is_same_or_ancestor "$PIPELINE_OUTPUT" "$BENCHMARK_ROOT" ||
   is_same_or_ancestor "$BENCHMARK_ROOT" "$PIPELINE_OUTPUT"
then
    die "Pipeline output and benchmark root must not contain one another."
fi

command -v python3 >/dev/null 2>&1 || die "python3 is required."
command -v git >/dev/null 2>&1 || die "git is required."
[[ -s "$PIPELINE_SCRIPT" ]] || die "Missing pipeline script: $PIPELINE_SCRIPT"
[[ -s "$RESOURCE_WRAPPER" ]] || die "Missing benchmark wrapper: $RESOURCE_WRAPPER"
[[ -s "$REPORTER" ]] || die "Missing pipeline reporter: $REPORTER"

bash -n "$PIPELINE_SCRIPT" || die "MTD_explorer.sh failed Bash syntax validation."
bash -n "$RESOURCE_WRAPPER" || die "MTD_benchmark_install.sh failed Bash syntax validation."
python3 -m py_compile "$REPORTER" || die "Pipeline reporter failed Python syntax validation."

HELP_OUTPUT="$(bash "$PIPELINE_SCRIPT" --help 2>&1 || true)"
for required_option in \
    "--input" \
    "--output" \
    "--hostid" \
    "--threads"
do
    grep -q -- "$required_option" <<< "$HELP_OUTPUT" || \
        die "Current MTD_explorer.sh does not advertise required option: $required_option"
done

if (( EXTRACT_READS )); then
    for required_option in \
        "--extract-microbiome-reads" \
        "--extract-microbiome-reads-top-n"
    do
        grep -q -- "$required_option" <<< "$HELP_OUTPUT" || \
            die "Current MTD_explorer.sh does not advertise required extraction option: $required_option"
    done
fi

if git -C "$REPO_DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    GIT_COMMIT="$(git -C "$REPO_DIR" rev-parse HEAD)"
    GIT_BRANCH="$(git -C "$REPO_DIR" branch --show-current)"
    if [[ -n "$(git -C "$REPO_DIR" status --short)" ]]; then
        GIT_STATE="dirty"
        warn "Repository contains uncommitted changes; they will be recorded."
    else
        GIT_STATE="clean"
    fi
else
    GIT_COMMIT="not-a-git-repository"
    GIT_BRANCH="NA"
    GIT_STATE="unknown"
fi

if (( USE_BLAST )); then
    ALIGNMENT_LABEL="magicblast"
else
    ALIGNMENT_LABEL="hisat2"
fi

if (( EXTRACT_READS )); then
    if (( TOP_N == 0 )); then
        EXTRACT_LABEL="extract_all"
    else
        EXTRACT_LABEL="extract_top${TOP_N}"
    fi
else
    EXTRACT_LABEL="no_extract"
fi

BENCHMARK_LABEL="${MACHINE_SAFE}_${DATASET_SAFE}_warm_host${HOST_ID}_${ALIGNMENT_LABEL}_${EXTRACT_LABEL}_${READ_LAYOUT}_t${THREADS}_r${RUN_NUMBER}"

CMD=(
    bash "$PIPELINE_SCRIPT"
    --input "$SAMPLESHEET"
    --output "$PIPELINE_OUTPUT"
    --hostid "$HOST_ID"
    --threads "$THREADS"
    --read-layout "$READ_LAYOUT"
    --analysis-mode "$ANALYSIS_MODE"
)

if (( USE_BLAST )); then
    CMD+=(--blast)
fi

if (( EXTRACT_READS )); then
    CMD+=(
        --extract-microbiome-reads
        --extract-microbiome-reads-top-n "$TOP_N"
    )
fi

if ((${#EXTRA_MTD_ARGS[@]} > 0)); then
    CMD+=("${EXTRA_MTD_ARGS[@]}")
fi

mkdir -p -- "$BENCHMARK_ROOT"
BENCHMARK_ROOT="$(readlink -f -- "$BENCHMARK_ROOT")"

printf '%s\n' "============================================================"
printf '%s\n' "MTD EXPLORER PIPELINE BENCHMARK PRECHECK"
printf '%-22s %s\n' "Machine:" "$MACHINE_NAME"
printf '%-22s %s\n' "Dataset:" "$DATASET_LABEL"
printf '%-22s %s\n' "Samplesheet:" "$SAMPLESHEET"
printf '%-22s %s\n' "Pipeline output:" "$PIPELINE_OUTPUT"
printf '%-22s %s\n' "Host TaxID:" "$HOST_ID"
printf '%-22s %s\n' "Threads:" "$THREADS"
printf '%-22s %s\n' "Read layout:" "$READ_LAYOUT"
printf '%-22s %s\n' "Analysis mode:" "$ANALYSIS_MODE"
printf '%-22s %s\n' "Alignment:" "$ALIGNMENT_LABEL"
printf '%-22s %s\n' "Read extraction:" "$EXTRACT_LABEL"
if (( EXTRACT_READS )); then
    if (( TOP_N == 0 )); then
        printf '%-22s %s\n' "Extraction policy:" "all detected taxa; stress-test mode"
    else
        printf '%-22s %s\n' "Extraction policy:" "top $TOP_N taxa from absolute ranking"
    fi
else
    printf '%-22s %s\n' "Extraction policy:" "disabled for the standard benchmark"
fi
printf '%-22s %s\n' "Run number:" "$RUN_NUMBER"
printf '%-22s %s\n' "Sampling interval:" "$INTERVAL s"
printf '%-22s %s\n' "Benchmark label:" "$BENCHMARK_LABEL"
printf '%-22s %s\n' "Benchmark root:" "$BENCHMARK_ROOT"
printf '%-22s %s\n' "Git commit:" "$GIT_COMMIT"
printf '%-22s %s\n' "Git branch/state:" "$GIT_BRANCH / $GIT_STATE"
printf '%s\n' "============================================================"
printf '[COMMAND] '
printf '%q ' "${CMD[@]}"
printf '\n\n'

if [[ -e "$PIPELINE_OUTPUT" || -L "$PIPELINE_OUTPUT" ]]; then
    printf '[REMOVE] Existing pipeline output: %s\n' "$PIPELINE_OUTPUT"
else
    printf '[CREATE] Pipeline output: %s\n' "$PIPELINE_OUTPUT"
fi
printf '[PRESERVE] Input data and installed databases are not deleted.\n'
printf '[PRESERVE] Benchmark results: %s\n' "$BENCHMARK_ROOT"

if (( ! ASSUME_YES )); then
    read -r -p "Type RUN to delete the previous pipeline output and start: " confirmation
    [[ "$confirmation" == "RUN" ]] || die "Cancelled. Nothing was removed."
fi

rm -rf -- "$PIPELINE_OUTPUT"
mkdir -p -- "$(dirname "$PIPELINE_OUTPUT")"

STEP_TEMP="$BENCHMARK_ROOT/.${BENCHMARK_LABEL}.pipeline_steps.$$.tsv"
rm -f -- "$STEP_TEMP"
printf 'timestamp_utc\tepoch_ns\tpercent\tmessage\n' > "$STEP_TEMP"

export MTD_PIPELINE_BENCHMARK_STEPS_TSV="$STEP_TEMP"

BENCHMARKED_CMD=(
    bash -c '
        printf "%s\t%s\t0\tMTD Explorer initialization and input validation\n" \
            "$(date -u "+%Y-%m-%dT%H:%M:%S.%NZ")" \
            "$(date +%s%N)" \
            >> "$MTD_PIPELINE_BENCHMARK_STEPS_TSV"

        set +e
        "$@"
        status=$?
        set -e

        printf "%s\t%s\t100\t__BENCHMARK_WRAPPER_END__\n" \
            "$(date -u "+%Y-%m-%dT%H:%M:%S.%NZ")" \
            "$(date +%s%N)" \
            >> "$MTD_PIPELINE_BENCHMARK_STEPS_TSV"

        exit "$status"
    ' mtd-pipeline-benchmark-command
    "${CMD[@]}"
)

cd "$REPO_DIR"

set +e
bash "$RESOURCE_WRAPPER" \
    --label "$BENCHMARK_LABEL" \
    --interval "$INTERVAL" \
    --output-root "$BENCHMARK_ROOT" \
    --watch-path "$PIPELINE_OUTPUT" \
    -- \
    "${BENCHMARKED_CMD[@]}"
BENCHMARK_STATUS=$?
set -e

LATEST_BENCHMARK="$(
    find "$BENCHMARK_ROOT" \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        -name "${BENCHMARK_LABEL}__*" \
        -printf '%T@ %p\n' 2>/dev/null |
    sort -nr |
    head -n 1 |
    cut -d' ' -f2-
)"

REPORT_STATUS=0
BUNDLE_PATH=""

if [[ -n "$LATEST_BENCHMARK" && -d "$LATEST_BENCHMARK" ]]; then
    cp -- "$STEP_TEMP" "$LATEST_BENCHMARK/pipeline_steps_raw.tsv"

    {
        printf '#!/usr/bin/env bash\n'
        printf '%q ' "${CMD[@]}"
        printf '\n'
    } > "$LATEST_BENCHMARK/pipeline_command.sh"
    chmod +x "$LATEST_BENCHMARK/pipeline_command.sh"

    set +e
    python3 "$REPORTER" \
        --run-dir "$LATEST_BENCHMARK" \
        --pipeline-output "$PIPELINE_OUTPUT" \
        --samplesheet "$SAMPLESHEET" \
        --repo-dir "$REPO_DIR" \
        --dataset-label "$DATASET_LABEL" \
        --machine "$MACHINE_NAME" \
        --hostid "$HOST_ID" \
        --threads "$THREADS" \
        --top-n "$TOP_N" \
        --read-layout "$READ_LAYOUT" \
        --analysis-mode "$ANALYSIS_MODE" \
        --alignment "$ALIGNMENT_LABEL" \
        --extract-reads "$EXTRACT_READS" \
        --exit-status "$BENCHMARK_STATUS"
    REPORT_STATUS=$?
    set -e

    if (( REPORT_STATUS != 0 )); then
        warn "Pipeline-specific report failed with exit status $REPORT_STATUS."
    fi

    BUNDLE_PATH="$BENCHMARK_ROOT/mtd_pipeline_benchmark_bundle_${BENCHMARK_LABEL}_$(date -u '+%Y%m%dT%H%M%SZ').tar.gz"
    tar -czf "$BUNDLE_PATH" \
        -C "$BENCHMARK_ROOT" \
        "$(basename "$LATEST_BENCHMARK")"
    sha256sum "$BUNDLE_PATH" > "${BUNDLE_PATH}.sha256"
else
    warn "Benchmark result directory was not found."
fi

rm -f -- "$STEP_TEMP"
unset MTD_PIPELINE_BENCHMARK_STEPS_TSV

echo
printf '%s\n' "============================================================"
printf '%s\n' "MTD PIPELINE BENCHMARK FINISHED"
printf 'Exit status: %s\n' "$BENCHMARK_STATUS"
if [[ -n "$LATEST_BENCHMARK" && -d "$LATEST_BENCHMARK" ]]; then
    printf 'Results: %s\n' "$LATEST_BENCHMARK"
    [[ -s "$LATEST_BENCHMARK/pipeline_summary.tsv" ]] && \
        printf 'Pipeline summary: %s\n' "$LATEST_BENCHMARK/pipeline_summary.tsv"
    [[ -s "$LATEST_BENCHMARK/pipeline_steps.tsv" ]] && \
        printf 'Stage timings: %s\n' "$LATEST_BENCHMARK/pipeline_steps.tsv"
    [[ -s "$LATEST_BENCHMARK/failure_report.txt" ]] && \
        printf 'Failure report: %s\n' "$LATEST_BENCHMARK/failure_report.txt"
fi
[[ -n "$BUNDLE_PATH" ]] && printf 'Bundle: %s\n' "$BUNDLE_PATH"
printf 'Pipeline output: %s\n' "$PIPELINE_OUTPUT"
printf '%s\n' "============================================================"

if [[ -n "$LATEST_BENCHMARK" && -s "$LATEST_BENCHMARK/pipeline_summary.txt" ]]; then
    echo
    sed -n '1,80p' "$LATEST_BENCHMARK/pipeline_summary.txt"
fi

exit "$BENCHMARK_STATUS"
