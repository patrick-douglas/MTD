#!/usr/bin/env bash
set -Eeuo pipefail

# Slurm copies this submitted script into /var/spool/slurm. Therefore,
# BASH_SOURCE[0] no longer points to the shared MTD Explorer checkout.
# Parse MTD_ROOT first and load the common helpers through the NFS path.
[[ $# -eq 4 ]] || {
    printf '[MTD-HPC ERROR] Internal array task usage: CONF MANIFEST SUCCESS_DIR MTD_ROOT\n' >&2
    exit 2
}

MTD_HPC_CONF="$1"
MTD_HPC_PENDING_MANIFEST="$2"
MTD_HPC_SUCCESS_DIR="$3"
MTD_HPC_MTD_ROOT="$4"

MTD_HPC_COMMON_SH="$MTD_HPC_MTD_ROOT/aux_scripts/hpc/mtd_hpc_common.sh"

if [[ ! -s "$MTD_HPC_COMMON_SH" ]]; then
    printf '[MTD-HPC ERROR] Shared HPC helper was not found: %s\n'         "$MTD_HPC_COMMON_SH" >&2
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
rm -f -- "$success_marker" "$failure_marker"

record_failure() {
    local status="$1"
    local reason="$2"
    printf 'hash=%s\nexit_code=%s\nhost=%s\nreason=%s\nfinished_at=%s\n' \
        "$task_hash" "$status" "$task_node" "$reason" "$(date --iso-8601=seconds)" \
        > "$failure_marker"
}

mtd_hpc_info "Task ID: $task_id"
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

printf 'hash=%s\nexit_code=0\nhost=%s\nthreads=%s\nmemory_kb=%s\nfinished_at=%s\n' \
    "$task_hash" "$task_node" "$MTD_NODE_THREADS" "$MTD_NODE_MEMORY_KB" \
    "$(date --iso-8601=seconds)" > "$success_marker"

mtd_hpc_ok "Task completed: $task_id"
