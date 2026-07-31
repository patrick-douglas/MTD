#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=aux_scripts/hpc/mtd_hpc_common.sh
source "$SCRIPT_DIR/mtd_hpc_common.sh"

[[ $# -eq 4 ]] || {
    mtd_hpc_error "Internal array task usage: CONF MANIFEST SUCCESS_DIR MTD_ROOT"
    exit 2
}

MTD_HPC_CONF="$1"
MTD_HPC_PENDING_MANIFEST="$2"
MTD_HPC_SUCCESS_DIR="$3"
MTD_HPC_MTD_ROOT="$4"
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
success_marker="$MTD_HPC_SUCCESS_DIR/${task_id}.success"
failure_marker="$MTD_HPC_SUCCESS_DIR/${task_id}.failed"
rm -f -- "$success_marker" "$failure_marker"

record_failure() {
    local status="$1"
    local reason="$2"
    printf 'hash=%s\nexit_code=%s\nhost=%s\nreason=%s\nfinished_at=%s\n' \
        "$task_hash" "$status" "$(hostname -s)" "$reason" "$(date --iso-8601=seconds)" \
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
        [[ -s "$expected" ]] || mtd_hpc_error "Missing output: $expected"
    done <<< "$expected_outputs"
    exit 90
fi

printf 'hash=%s\nexit_code=0\nhost=%s\nthreads=%s\nmemory_kb=%s\nfinished_at=%s\n' \
    "$task_hash" "$(hostname -s)" "$MTD_NODE_THREADS" "$MTD_NODE_MEMORY_KB" \
    "$(date --iso-8601=seconds)" > "$success_marker"

mtd_hpc_ok "Task completed: $task_id"
