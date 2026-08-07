#!/usr/bin/env bash
set -Eeuo pipefail

# Slurm copies this submitted script into /var/spool/slurm. Therefore,
# BASH_SOURCE[0] no longer points to the shared MTD Explorer checkout.
# Parse MTD_ROOT first and load the common helpers through the NFS path.
[[ $# -ge 4 && $# -le 6 ]] || {
    printf '[MTD-HPC ERROR] Internal array task usage: CONF MANIFEST SUCCESS_DIR MTD_ROOT [STAGE] [ATTEMPT]\n' >&2
    exit 2
}

MTD_HPC_CONF="$1"
MTD_HPC_PENDING_MANIFEST="$2"
MTD_HPC_SUCCESS_DIR="$3"
MTD_HPC_MTD_ROOT="$4"
MTD_HPC_STAGE_NAME="${5:-unknown}"
MTD_HPC_ATTEMPT_LABEL="${6:-unknown}"

MTD_HPC_COMMON_SH="$MTD_HPC_MTD_ROOT/aux_scripts/hpc/mtd_hpc_common.sh"

if [[ ! -s "$MTD_HPC_COMMON_SH" ]]; then
    printf '[MTD-HPC ERROR] Shared HPC helper was not found: %s\n' \
        "$MTD_HPC_COMMON_SH" >&2
    exit 1
fi

# shellcheck source=aux_scripts/hpc/mtd_hpc_common.sh
source "$MTD_HPC_COMMON_SH"

: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

mtd_hpc_load_config "$MTD_HPC_CONF" "$MTD_HPC_MTD_ROOT" || exit 1
mtd_hpc_export_node_resources

line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MTD_HPC_PENDING_MANIFEST")"
[[ -n "$line" ]] || {
    mtd_hpc_error "No task exists for array index $SLURM_ARRAY_TASK_ID"
    exit 2
}

IFS=$'\t' read -r task_id task_hash command_b64 expected_b64 <<< "$line"
[[ -n "$task_id" && -n "$task_hash" && -n "$command_b64" ]] || {
    mtd_hpc_error "Malformed task line at array index $SLURM_ARRAY_TASK_ID"
    exit 2
}
mtd_hpc_validate_id "$task_id" "task ID" || exit 2

command="$(mtd_hpc_b64_decode "$command_b64")"
expected_outputs="$(mtd_hpc_b64_decode "$expected_b64")"

mkdir -p -- "$MTD_HPC_SUCCESS_DIR"
task_node="${SLURMD_NODENAME:-$(hostname -s)}"
success_marker="$MTD_HPC_SUCCESS_DIR/${task_id}.success"
failure_marker="$MTD_HPC_SUCCESS_DIR/${task_id}.failed"
attempt_marker_dir="$MTD_HPC_SUCCESS_DIR/attempts"
mkdir -p -- "$attempt_marker_dir"
attempt_marker_prefix="$attempt_marker_dir/${SLURM_JOB_ID:-unknown}_${SLURM_ARRAY_TASK_ID:-unknown}_${task_id}"
rm -f -- "$success_marker" "$failure_marker"

task_started_at="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
task_start_epoch_ns="$(date +%s%N)"
task_architecture="$(uname -m 2>/dev/null || printf unknown)"
task_kernel="$(uname -r 2>/dev/null || printf unknown)"
task_os_release="$(
    if [[ -r /etc/os-release ]]; then
        sed -n 's/^PRETTY_NAME=//p' /etc/os-release |
            head -n 1 |
            sed 's/^"//; s/"$//'
    fi
)"
: "${task_os_release:=unknown}"
task_memory_total_kb="$(
    awk '/^MemTotal:/ {print $2; exit}' /proc/meminfo 2>/dev/null || true
)"
: "${task_memory_total_kb:=0}"

task_cpu_model="unknown"
task_logical_cpus="$(nproc --all 2>/dev/null || nproc 2>/dev/null || printf 1)"
task_sockets="unknown"
task_cores_per_socket="unknown"
task_threads_per_core="unknown"
if command -v lscpu >/dev/null 2>&1; then
    task_cpu_model="$(
        LC_ALL=C lscpu |
            awk -F: '$1 ~ /^Model name/ {sub(/^[[:space:]]+/, "", $2); print $2; exit}'
    )"
    task_sockets="$(
        LC_ALL=C lscpu |
            awk -F: '$1 ~ /^Socket\(s\)/ {gsub(/[[:space:]]/, "", $2); print $2; exit}'
    )"
    task_cores_per_socket="$(
        LC_ALL=C lscpu |
            awk -F: '$1 ~ /^Core\(s\) per socket/ {gsub(/[[:space:]]/, "", $2); print $2; exit}'
    )"
    task_threads_per_core="$(
        LC_ALL=C lscpu |
            awk -F: '$1 ~ /^Thread\(s\) per core/ {gsub(/[[:space:]]/, "", $2); print $2; exit}'
    )"
fi
: "${task_cpu_model:=unknown}"
: "${task_sockets:=unknown}"
: "${task_cores_per_socket:=unknown}"
: "${task_threads_per_core:=unknown}"

task_finished_at=""
task_finish_epoch_ns=""
task_elapsed_seconds=""

finish_task_timing() {
    [[ -z "$task_finish_epoch_ns" ]] || return 0
    task_finished_at="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    task_finish_epoch_ns="$(date +%s%N)"
    task_elapsed_seconds="$(
        awk -v start="$task_start_epoch_ns" -v finish="$task_finish_epoch_ns" \
            'BEGIN {printf "%.6f", (finish-start)/1000000000}'
    )"
}

write_task_marker() {
    local marker_path="$1"
    local status="$2"
    local reason="$3"
    local temporary_marker="${marker_path}.tmp.$$"

    finish_task_timing
    {
        printf 'hash=%s\n' "$task_hash"
        printf 'task_id=%s\n' "$task_id"
        printf 'stage=%s\n' "$MTD_HPC_STAGE_NAME"
        printf 'attempt=%s\n' "$MTD_HPC_ATTEMPT_LABEL"
        printf 'job_id=%s\n' "${SLURM_JOB_ID:-unknown}"
        printf 'array_task_id=%s\n' "${SLURM_ARRAY_TASK_ID:-unknown}"
        printf 'job_name=%s\n' "${SLURM_JOB_NAME:-unknown}"
        printf 'exit_code=%s\n' "$status"
        printf 'host=%s\n' "$task_node"
        printf 'reason=%s\n' "$reason"
        printf 'threads=%s\n' "$MTD_NODE_THREADS"
        printf 'memory_kb=%s\n' "$MTD_NODE_MEMORY_KB"
        printf 'memory_total_kb=%s\n' "$task_memory_total_kb"
        printf 'slurm_cpus_per_task=%s\n' "${SLURM_CPUS_PER_TASK:-unknown}"
        printf 'slurm_cpus_on_node=%s\n' "${SLURM_CPUS_ON_NODE:-unknown}"
        printf 'architecture=%s\n' "$task_architecture"
        printf 'cpu_model=%s\n' "$task_cpu_model"
        printf 'logical_cpus=%s\n' "$task_logical_cpus"
        printf 'sockets=%s\n' "$task_sockets"
        printf 'cores_per_socket=%s\n' "$task_cores_per_socket"
        printf 'threads_per_core=%s\n' "$task_threads_per_core"
        printf 'kernel=%s\n' "$task_kernel"
        printf 'os_release=%s\n' "$task_os_release"
        printf 'started_at=%s\n' "$task_started_at"
        printf 'started_epoch_ns=%s\n' "$task_start_epoch_ns"
        printf 'finished_at=%s\n' "$task_finished_at"
        printf 'finished_epoch_ns=%s\n' "$task_finish_epoch_ns"
        printf 'elapsed_seconds=%s\n' "$task_elapsed_seconds"
    } > "$temporary_marker"
    mv -f -- "$temporary_marker" "$marker_path"
}

write_task_markers() {
    local canonical_marker="$1"
    local attempt_suffix="$2"
    local status="$3"
    local reason="$4"

    finish_task_timing
    write_task_marker "$canonical_marker" "$status" "$reason"
    write_task_marker "${attempt_marker_prefix}.${attempt_suffix}" "$status" "$reason"
}

record_failure() {
    local status="$1"
    local reason="$2"
    write_task_markers "$failure_marker" failed "$status" "$reason"
}

mtd_hpc_info "Task ID: $task_id"
mtd_hpc_info "Stage: $MTD_HPC_STAGE_NAME"
mtd_hpc_info "Attempt: $MTD_HPC_ATTEMPT_LABEL"
mtd_hpc_info "Command hash: $task_hash"
mtd_hpc_info "Command: $command"

set +e
bash -c "$command"
status=$?
set -e

if (( status != 0 )); then
    record_failure "$status" "command_failed"
    mtd_hpc_error "Task failed with exit code $status: $task_id"
    exit "$status"
fi

if ! mtd_hpc_outputs_exist "$expected_outputs"; then
    record_failure 90 "expected_output_missing_or_empty"
    mtd_hpc_error "At least one expected output is missing or empty for task: $task_id"
    while IFS= read -r expected; do
        [[ -z "$expected" ]] && continue
        if ! mtd_hpc_output_spec_is_valid "$expected"; then
            mtd_hpc_error "Missing output: $(mtd_hpc_output_spec_path "$expected")"
        fi
    done <<< "$expected_outputs"
    exit 90
fi

write_task_markers "$success_marker" success 0 "completed"

mtd_hpc_ok "Task completed: $task_id"
