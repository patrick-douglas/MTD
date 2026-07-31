#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/mtd_hpc_common.sh"

usage() {
    cat <<'EOF'
Usage:
  mtd_hpc_submit_array.sh \
    --hpc-conf FILE \
    --manifest FILE \
    --stage NAME \
    --work-dir DIR \
    [--mtd-root DIR]
EOF
}

HPC_CONF=""
MANIFEST=""
STAGE=""
WORK_DIR=""
MTD_ROOT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hpc-conf) HPC_CONF="$2"; shift 2 ;;
        --hpc-conf=*) HPC_CONF="${1#*=}"; shift ;;
        --manifest) MANIFEST="$2"; shift 2 ;;
        --manifest=*) MANIFEST="${1#*=}"; shift ;;
        --stage) STAGE="$2"; shift 2 ;;
        --stage=*) STAGE="${1#*=}"; shift ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --work-dir=*) WORK_DIR="${1#*=}"; shift ;;
        --mtd-root) MTD_ROOT="$2"; shift 2 ;;
        --mtd-root=*) MTD_ROOT="${1#*=}"; shift ;;
        --help) usage; exit 0 ;;
        *) mtd_hpc_error "Unknown option: $1"; usage; exit 2 ;;
    esac
done

[[ -n "$HPC_CONF" && -n "$MANIFEST" && -n "$STAGE" && -n "$WORK_DIR" && -n "$MTD_ROOT" ]] || {
    usage
    exit 2
}

mtd_hpc_validate_id "$STAGE" "stage name"
mtd_hpc_load_config "$HPC_CONF" "$MTD_ROOT"
mtd_hpc_require_file "$MANIFEST" "HPC task manifest"

for command in sbatch squeue sacct scancel base64 sha256sum; do
    command -v "$command" >/dev/null 2>&1 || mtd_hpc_die "Required command not found: $command"
done

mkdir -p -- "$WORK_DIR"
WORK_DIR="$(mtd_hpc_realpath "$WORK_DIR")"
MANIFEST="$(mtd_hpc_realpath "$MANIFEST")"
MTD_ROOT="$(mtd_hpc_realpath "$MTD_ROOT")"

log_dir="$WORK_DIR/logs"
success_dir="$WORK_DIR/success"
pending_manifest="$WORK_DIR/${STAGE}.pending.tsv"
mkdir -p -- "$log_dir" "$success_dir"
: > "$pending_manifest"

expected_tasks=0
pending_tasks=0
while IFS=$'\t' read -r task_id task_hash command_b64 expected_b64; do
    [[ -z "$task_id" || "$task_id" == \#* ]] && continue
    mtd_hpc_validate_id "$task_id" "task ID"
    expected_tasks=$((expected_tasks + 1))
    marker="$success_dir/${task_id}.success"
    expected_outputs="$(mtd_hpc_b64_decode "$expected_b64")"

    if [[ "$MTD_HPC_RESUME" == "1" ]] && [[ -s "$marker" ]] && \
       grep -Fxq "hash=$task_hash" "$marker" && \
       mtd_hpc_outputs_exist "$expected_outputs"; then
        mtd_hpc_info "Reusing completed task: $task_id"
        continue
    fi

    rm -f -- "$marker" "$success_dir/${task_id}.failed"
    printf '%s\t%s\t%s\t%s\n' "$task_id" "$task_hash" "$command_b64" "$expected_b64" \
        >> "$pending_manifest"
    pending_tasks=$((pending_tasks + 1))
done < "$MANIFEST"

(( expected_tasks > 0 )) || mtd_hpc_die "Task manifest contains no tasks: $MANIFEST"

if (( pending_tasks == 0 )); then
    mtd_hpc_ok "All $expected_tasks tasks already have matching outputs and success markers."
    exit 0
fi

sbatch_args=(
    --parsable
    --job-name="MTD_${STAGE}"
    --array="1-${pending_tasks}%${MTD_HPC_MAX_PARALLEL}"
    --output="$log_dir/%A_%a.out"
    --error="$log_dir/%A_%a.err"
)

[[ -n "$MTD_HPC_PARTITION" ]] && sbatch_args+=(--partition="$MTD_HPC_PARTITION")
[[ -n "$MTD_HPC_ACCOUNT" ]] && sbatch_args+=(--account="$MTD_HPC_ACCOUNT")
[[ -n "$MTD_HPC_QOS" ]] && sbatch_args+=(--qos="$MTD_HPC_QOS")
[[ -n "$MTD_HPC_CONSTRAINT" ]] && sbatch_args+=(--constraint="$MTD_HPC_CONSTRAINT")
[[ -n "$MTD_HPC_TIME" ]] && sbatch_args+=(--time="$MTD_HPC_TIME")
sbatch_args+=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")

mtd_hpc_info "Submitting $pending_tasks pending tasks ($expected_tasks total)."
job_id="$(sbatch "${sbatch_args[@]}" \
    "$SCRIPT_DIR/mtd_hpc_array_task.sh" \
    "$MTD_HPC_CONF" "$pending_manifest" "$success_dir" "$MTD_ROOT")"
job_id="${job_id%%;*}"
[[ "$job_id" =~ ^[0-9]+$ ]] || mtd_hpc_die "Could not parse Slurm job ID: $job_id"

printf '%s\n' "$job_id" > "$WORK_DIR/${STAGE}.job_id"
mtd_hpc_info "Submitted Slurm array job: $job_id"

cancel_active_job() {
    mtd_hpc_error "Interrupted; cancelling Slurm job $job_id."
    scancel "$job_id" >/dev/null 2>&1 || true
    exit 130
}
trap cancel_active_job INT TERM HUP

while squeue -h -j "$job_id" 2>/dev/null | grep -q .; do
    completed_now=0
    while IFS=$'\t' read -r task_id _ _ _; do
        [[ -z "$task_id" ]] && continue
        [[ -s "$success_dir/${task_id}.success" ]] && completed_now=$((completed_now + 1))
    done < "$pending_manifest"
    mtd_hpc_info "Job $job_id running; new success markers: $completed_now/$pending_tasks"
    sleep "$MTD_HPC_POLL_SECONDS"
done

# sacct can lag briefly after a job leaves squeue.
sacct_output=""
for _ in 1 2 3 4 5; do
    sacct_output="$(sacct -n -P -j "$job_id" --format=JobIDRaw,State,ExitCode 2>/dev/null || true)"
    [[ -n "$sacct_output" ]] && break
    sleep 2
done
printf '%s\n' "$sacct_output" > "$WORK_DIR/${STAGE}.sacct.tsv"
trap - INT TERM HUP

missing=0
while IFS=$'\t' read -r task_id task_hash _ expected_b64; do
    [[ -z "$task_id" || "$task_id" == \#* ]] && continue
    marker="$success_dir/${task_id}.success"
    expected_outputs="$(mtd_hpc_b64_decode "$expected_b64")"
    if [[ ! -s "$marker" ]] || ! grep -Fxq "hash=$task_hash" "$marker" || \
       ! mtd_hpc_outputs_exist "$expected_outputs"; then
        mtd_hpc_error "Task did not produce a matching success marker and outputs: $task_id"
        missing=$((missing + 1))
    fi
done < "$MANIFEST"

if (( missing > 0 )); then
    mtd_hpc_error "$missing task(s) failed or did not finish correctly."
    mtd_hpc_error "Logs: $log_dir"
    [[ -n "$sacct_output" ]] && printf '%s\n' "$sacct_output" >&2
    exit 1
fi

mtd_hpc_ok "All $expected_tasks tasks completed successfully."
