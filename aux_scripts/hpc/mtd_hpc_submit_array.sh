#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/mtd_hpc_common.sh"

usage() {
    cat <<'USAGE'
Usage:
  mtd_hpc_submit_array.sh \
    --hpc-conf FILE \
    --manifest FILE \
    --stage NAME \
    --work-dir DIR \
    [--mtd-root DIR]
USAGE
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

for command in sbatch squeue sacct scancel scontrol base64 sha256sum; do
    command -v "$command" >/dev/null 2>&1 || mtd_hpc_die "Required command not found: $command"
done

mkdir -p -- "$WORK_DIR"
WORK_DIR="$(mtd_hpc_realpath "$WORK_DIR")"
MANIFEST="$(mtd_hpc_realpath "$MANIFEST")"
MTD_ROOT="$(mtd_hpc_realpath "$MTD_ROOT")"

log_root="$WORK_DIR/logs"
success_dir="$WORK_DIR/success"
latest_pending_manifest="$WORK_DIR/${STAGE}.pending.tsv"
job_history="$WORK_DIR/${STAGE}.job_ids.tsv"
sacct_history="$WORK_DIR/${STAGE}.sacct.tsv"
retry_history="$WORK_DIR/${STAGE}.retry.tsv"
mkdir -p -- "$log_root" "$success_dir"
: > "$job_history"
: > "$sacct_history"
: > "$retry_history"

mtd_hpc_task_complete() {
    local task_id="$1"
    local task_hash="$2"
    local expected_b64="$3"
    local marker="$success_dir/${task_id}.success"
    local expected_outputs=""

    [[ -s "$marker" ]] || return 1
    grep -Fxq "hash=$task_hash" "$marker" || return 1
    expected_outputs="$(mtd_hpc_b64_decode "$expected_b64")"
    mtd_hpc_outputs_exist "$expected_outputs"
}

mtd_hpc_count_completed() {
    local manifest="$1"
    local count=0
    local task_id task_hash command_b64 expected_b64

    while IFS=$'\t' read -r task_id task_hash command_b64 expected_b64; do
        [[ -z "$task_id" || "$task_id" == \#* ]] && continue
        if mtd_hpc_task_complete "$task_id" "$task_hash" "$expected_b64"; then
            count=$((count + 1))
        fi
    done < "$manifest"

    printf '%s\n' "$count"
}

mtd_hpc_count_manifest_tasks() {
    local manifest="$1"
    awk -F '\t' 'NF && $1 !~ /^#/ {count++} END {print count + 0}' "$manifest"
}

mtd_hpc_select_resource_policy() {
    case "$STAGE" in
        humann)
            if declare -p MTD_HPC_HUMANN_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
                selected_policy="MTD_HPC_HUMANN_SBATCH_EXTRA_ARGS"
                stage_sbatch_args=("${MTD_HPC_HUMANN_SBATCH_EXTRA_ARGS[@]}")
            else
                selected_policy="MTD_HPC_SBATCH_EXTRA_ARGS (fallback)"
                stage_sbatch_args=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")
            fi
            ;;

        magicblast)
            if declare -p MTD_HPC_MAGICBLAST_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
                selected_policy="MTD_HPC_MAGICBLAST_SBATCH_EXTRA_ARGS"
                stage_sbatch_args=("${MTD_HPC_MAGICBLAST_SBATCH_EXTRA_ARGS[@]}")
            else
                selected_policy="MTD_HPC_SBATCH_EXTRA_ARGS (fallback)"
                stage_sbatch_args=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")
            fi
            ;;

        fastp)
            if declare -p MTD_HPC_FASTP_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
                selected_policy="MTD_HPC_FASTP_SBATCH_EXTRA_ARGS"
                stage_sbatch_args=("${MTD_HPC_FASTP_SBATCH_EXTRA_ARGS[@]}")
            else
                selected_policy="MTD_HPC_SBATCH_EXTRA_ARGS (fallback)"
                stage_sbatch_args=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")
            fi
            ;;

        kraken_host|kraken_micro_raw|kraken_micro_final|bracken)
            if declare -p MTD_HPC_KRAKEN_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
                selected_policy="MTD_HPC_KRAKEN_SBATCH_EXTRA_ARGS"
                stage_sbatch_args=("${MTD_HPC_KRAKEN_SBATCH_EXTRA_ARGS[@]}")
            else
                selected_policy="MTD_HPC_SBATCH_EXTRA_ARGS (fallback)"
                stage_sbatch_args=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")
            fi
            ;;

        *)
            selected_policy="MTD_HPC_SBATCH_EXTRA_ARGS"
            stage_sbatch_args=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")
            ;;
    esac
}

mtd_hpc_failure_host() {
    local task_id="$1"
    local task_hash="$2"
    local marker="$success_dir/${task_id}.failed"
    local host=""

    [[ -s "$marker" ]] || return 1
    grep -Fxq "hash=$task_hash" "$marker" || return 1
    host="$(sed -n 's/^host=//p' "$marker" | head -n 1)"
    [[ "$host" =~ ^[A-Za-z0-9._-]+$ ]] || return 1
    printf '%s\n' "$host"
}

mtd_hpc_resolve_submit_slurm_node() {
    local submit_host_short=""
    local submit_host_fqdn=""
    local submit_host_plain=""
    local line=""
    local field=""
    local node_name=""
    local node_host=""
    local node_addr=""
    local candidate=""
    local candidate_normalized=""
    local target=""
    local target_normalized=""
    local -a candidates=()
    local -A matches=()

    if [[ -n "${SLURMD_NODENAME:-}" ]] && \
       scontrol show node "$SLURMD_NODENAME" -o >/dev/null 2>&1; then
        printf '%s\n' "$SLURMD_NODENAME"
        return 0
    fi

    submit_host_short="$(hostname -s 2>/dev/null || true)"
    submit_host_fqdn="$(hostname -f 2>/dev/null || true)"
    submit_host_plain="$(hostname 2>/dev/null || true)"

    [[ -z "$submit_host_short" ]] || candidates+=("$submit_host_short")
    [[ -z "$submit_host_fqdn" ]] || candidates+=("$submit_host_fqdn")
    [[ -z "$submit_host_plain" ]] || candidates+=("$submit_host_plain")

    while IFS= read -r candidate; do
        [[ -z "$candidate" ]] || candidates+=("$candidate")
    done < <(
        hostname -I 2>/dev/null |
            tr ' ' '\n' |
            awk 'NF && $1 != "127.0.0.1" {print $1}'
    )

    while IFS= read -r line; do
        node_name=""
        node_host=""
        node_addr=""

        for field in $line; do
            case "$field" in
                NodeName=*) node_name="${field#NodeName=}" ;;
                NodeHostName=*) node_host="${field#NodeHostName=}" ;;
                NodeAddr=*) node_addr="${field#NodeAddr=}" ;;
            esac
        done

        [[ -n "$node_name" ]] || continue

        for candidate in "${candidates[@]}"; do
            candidate_normalized="${candidate%.}"
            candidate_normalized="${candidate_normalized,,}"

            for target in "$node_name" "$node_host" "$node_addr"; do
                [[ -n "$target" && "$target" != "(null)" ]] || continue
                target_normalized="${target%.}"
                target_normalized="${target_normalized,,}"

                if [[ "$candidate_normalized" == "$target_normalized" ]]; then
                    matches["$node_name"]=1
                    break 2
                fi
            done
        done
    done < <(scontrol show nodes -o 2>/dev/null)

    if (( ${#matches[@]} == 1 )); then
        printf '%s\n' "${!matches[@]}"
        return 0
    fi

    if (( ${#matches[@]} > 1 )); then
        return 2
    fi

    return 1
}

mtd_hpc_node_supports_partition() {
    local node_name="$1"
    local requested_partitions="$2"
    local node_info=""
    local field=""
    local node_partitions=""
    local requested=""
    local available=""
    local -a requested_list=()
    local -a available_list=()

    [[ -n "$requested_partitions" ]] || return 0

    node_info="$(scontrol show node "$node_name" -o 2>/dev/null)" || return 1
    for field in $node_info; do
        case "$field" in
            Partitions=*)
                node_partitions="${field#Partitions=}"
                break
                ;;
        esac
    done

    [[ -n "$node_partitions" && "$node_partitions" != "(null)" ]] || return 1

    IFS=',' read -r -a requested_list <<< "$requested_partitions"
    IFS=',' read -r -a available_list <<< "$node_partitions"

    for requested in "${requested_list[@]}"; do
        requested="${requested%\*}"
        for available in "${available_list[@]}"; do
            available="${available%\*}"
            [[ "$requested" == "$available" ]] && return 0
        done
    done

    return 1
}

mtd_hpc_build_fallback_stage_args() {
    local argument=""
    local skip_next=0

    fallback_stage_sbatch_args=()

    for argument in "${stage_sbatch_args[@]}"; do
        if (( skip_next == 1 )); then
            skip_next=0
            continue
        fi

        case "$argument" in
            --nodelist|--exclude|-w)
                skip_next=1
                ;;
            --nodelist=*|--exclude=*|-w?*)
                ;;
            *)
                fallback_stage_sbatch_args+=("$argument")
                ;;
        esac
    done
}

expected_tasks=0
initial_pending_manifest="$WORK_DIR/${STAGE}.attempt01.pending.tsv"
: > "$initial_pending_manifest"

while IFS=$'\t' read -r task_id task_hash command_b64 expected_b64; do
    [[ -z "$task_id" || "$task_id" == \#* ]] && continue
    mtd_hpc_validate_id "$task_id" "task ID"
    expected_tasks=$((expected_tasks + 1))

    if [[ "$MTD_HPC_RESUME" == "1" ]] && \
       mtd_hpc_task_complete "$task_id" "$task_hash" "$expected_b64"; then
        mtd_hpc_info "Reusing completed task: $task_id"
        continue
    fi

    rm -f -- "$success_dir/${task_id}.success" "$success_dir/${task_id}.failed"
    printf '%s\t%s\t%s\t%s\n' \
        "$task_id" "$task_hash" "$command_b64" "$expected_b64" \
        >> "$initial_pending_manifest"
done < "$MANIFEST"

(( expected_tasks > 0 )) || mtd_hpc_die "Task manifest contains no tasks: $MANIFEST"

initial_pending_tasks="$(mtd_hpc_count_manifest_tasks "$initial_pending_manifest")"
if (( initial_pending_tasks == 0 )); then
    mtd_hpc_ok "All $expected_tasks tasks already have matching outputs and success markers."
    exit 0
fi

mtd_hpc_select_resource_policy
printf -v stage_sbatch_args_text '%q ' "${stage_sbatch_args[@]}"
mtd_hpc_info "Slurm resource policy: $selected_policy"
mtd_hpc_info "Slurm resource arguments: ${stage_sbatch_args_text% }"
mtd_hpc_info "Automatic attempts per task: $MTD_HPC_MAX_ATTEMPTS"

submit_host_short="$(hostname -s 2>/dev/null || hostname 2>/dev/null || printf unknown)"
submit_host_fqdn="$(hostname -f 2>/dev/null || true)"
submit_slurm_node=""
submit_fallback_unavailable_reason=""

if [[ "$MTD_HPC_FINAL_SUBMIT_NODE_FALLBACK" == "1" ]]; then
    submit_node_resolution_status=0
    submit_slurm_node="$(mtd_hpc_resolve_submit_slurm_node)" || \
        submit_node_resolution_status=$?

    case "$submit_node_resolution_status" in
        0)
            mtd_hpc_validate_id "$submit_slurm_node" "submission Slurm node"
            if ! mtd_hpc_node_supports_partition \
                "$submit_slurm_node" "$MTD_HPC_PARTITION"; then
                submit_fallback_unavailable_reason="submission node $submit_slurm_node is not available in configured partition ${MTD_HPC_PARTITION:-<default>}"
                submit_slurm_node=""
            fi
            ;;
        2)
            submit_fallback_unavailable_reason="submission host matches more than one Slurm NodeName"
            ;;
        *)
            submit_fallback_unavailable_reason="submission host is not registered as a Slurm compute node"
            ;;
    esac

    mtd_hpc_info "Submission host: $submit_host_short"
    [[ -z "$submit_host_fqdn" ]] || mtd_hpc_info "Submission host FQDN: $submit_host_fqdn"

    if [[ -n "$submit_slurm_node" ]]; then
        mtd_hpc_info "Submission Slurm node: $submit_slurm_node"
        mtd_hpc_info \
            "Final submission-node fallback attempts: $MTD_HPC_FINAL_SUBMIT_NODE_ATTEMPTS"
    else
        mtd_hpc_info \
            "Final submission-node fallback is unavailable: $submit_fallback_unavailable_reason"
    fi
fi

mtd_hpc_build_fallback_stage_args

normal_attempt_limit="$MTD_HPC_MAX_ATTEMPTS"
fallback_attempt_limit=0
if [[ "$MTD_HPC_FINAL_SUBMIT_NODE_FALLBACK" == "1" ]] && \
   [[ -n "$submit_slurm_node" ]]; then
    fallback_attempt_limit="$MTD_HPC_FINAL_SUBMIT_NODE_ATTEMPTS"
fi
total_attempt_limit=$((normal_attempt_limit + fallback_attempt_limit))

active_job_id=""
cancel_active_job() {
    if [[ -n "$active_job_id" ]]; then
        mtd_hpc_error "Interrupted; cancelling Slurm job $active_job_id."
        scancel "$active_job_id" >/dev/null 2>&1 || true
    fi
    exit 130
}
trap cancel_active_job INT TERM HUP

declare -A retry_excluded_nodes=()
attempt=1
normal_attempts_run=0
fallback_attempts_run=0
current_manifest="$initial_pending_manifest"

while (( attempt <= total_attempt_limit )); do
    current_tasks="$(mtd_hpc_count_manifest_tasks "$current_manifest")"
    if (( current_tasks == 0 )); then
        break
    fi

    is_fallback=0
    fallback_attempt=0
    round_history_label=""
    round_display=""
    round_tag=""
    round_log_tag=""
    round_job_suffix=""
    round_stage_sbatch_args=()

    if (( attempt <= normal_attempt_limit )); then
        attempt_label="$(printf '%02d' "$attempt")"
        normal_attempts_run="$attempt"
        round_history_label="$attempt"
        round_display="Attempt $attempt/$normal_attempt_limit"
        round_tag="attempt${attempt_label}"
        round_log_tag="attempt_${attempt_label}"
        round_job_suffix="a${attempt_label}"
        round_stage_sbatch_args=("${stage_sbatch_args[@]}")
    else
        is_fallback=1
        fallback_attempt=$((attempt - normal_attempt_limit))
        fallback_label="$(printf '%02d' "$fallback_attempt")"
        fallback_attempts_run="$fallback_attempt"
        round_history_label="fallback_${fallback_attempt}"
        round_display="Submission-node fallback $fallback_attempt/$fallback_attempt_limit"
        round_tag="fallback${fallback_label}"
        round_log_tag="fallback_${fallback_label}"
        round_job_suffix="f${fallback_label}"
        round_stage_sbatch_args=("${fallback_stage_sbatch_args[@]}")
    fi

    cp -- "$current_manifest" "$latest_pending_manifest"
    attempt_log_dir="$log_root/${round_log_tag}"
    mkdir -p -- "$attempt_log_dir"

    sbatch_args=(
        --parsable
        --job-name="MTD_${STAGE}_${round_job_suffix}"
        --array="1-${current_tasks}%${MTD_HPC_MAX_PARALLEL}"
        --output="$attempt_log_dir/%A_%a.out"
        --error="$attempt_log_dir/%A_%a.err"
    )

    [[ -n "$MTD_HPC_PARTITION" ]] && sbatch_args+=(--partition="$MTD_HPC_PARTITION")
    [[ -n "$MTD_HPC_ACCOUNT" ]] && sbatch_args+=(--account="$MTD_HPC_ACCOUNT")
    [[ -n "$MTD_HPC_QOS" ]] && sbatch_args+=(--qos="$MTD_HPC_QOS")
    [[ -n "$MTD_HPC_CONSTRAINT" ]] && sbatch_args+=(--constraint="$MTD_HPC_CONSTRAINT")
    [[ -n "$MTD_HPC_TIME" ]] && sbatch_args+=(--time="$MTD_HPC_TIME")
    sbatch_args+=("${round_stage_sbatch_args[@]}")

    excluded_nodes_csv=""
    placement_record="-"

    if (( is_fallback == 1 )); then
        sbatch_args+=(--nodelist="$submit_slurm_node")
        placement_record="nodelist=$submit_slurm_node"
    elif [[ "$MTD_HPC_RETRY_EXCLUDE_FAILED_NODES" == "1" ]] && \
         (( attempt > 1 )) && (( ${#retry_excluded_nodes[@]} > 0 )); then
        excluded_nodes_csv="$(
            printf '%s\n' "${!retry_excluded_nodes[@]}" |
                LC_ALL=C sort |
                paste -sd, -
        )"
        if [[ -n "$excluded_nodes_csv" ]]; then
            sbatch_args+=(--exclude="$excluded_nodes_csv")
            placement_record="exclude=$excluded_nodes_csv"
        fi
    fi

    mtd_hpc_info \
        "$round_display: submitting $current_tasks task(s)" \
        "($expected_tasks total)."
    if (( is_fallback == 1 )); then
        mtd_hpc_info "Final fallback node: $submit_slurm_node"
    elif [[ -n "$excluded_nodes_csv" ]]; then
        mtd_hpc_info "Retry node exclusion: $excluded_nodes_csv"
    fi

    job_id="$(sbatch "${sbatch_args[@]}" \
        "$SCRIPT_DIR/mtd_hpc_array_task.sh" \
        "$MTD_HPC_CONF" "$current_manifest" "$success_dir" "$MTD_ROOT" \
        "$STAGE" "$round_history_label")"
    job_id="${job_id%%;*}"
    [[ "$job_id" =~ ^[0-9]+$ ]] || mtd_hpc_die "Could not parse Slurm job ID: $job_id"

    active_job_id="$job_id"
    printf '%s\n' "$job_id" > "$WORK_DIR/${STAGE}.job_id"
    printf '%s\n' "$job_id" > "$WORK_DIR/${STAGE}.${round_tag}.job_id"
    printf '%s\t%s\t%s\t%s\n' \
        "$round_history_label" "$job_id" "$current_tasks" "$placement_record" \
        >> "$job_history"
    mtd_hpc_info "Submitted Slurm array job: $job_id"

    while squeue -h -j "$job_id" 2>/dev/null | grep -q .; do
        attempt_completed="$(mtd_hpc_count_completed "$current_manifest")"
        overall_completed="$(mtd_hpc_count_completed "$MANIFEST")"
        mtd_hpc_info \
            "Job $job_id running; overall success: $overall_completed/$expected_tasks; "\
            "$round_display: $attempt_completed/$current_tasks"
        sleep "$MTD_HPC_POLL_SECONDS"
    done

    sacct_output=""
    sacct_fields=""
    for _ in 1 2 3 4 5; do
        for candidate_fields in \
            'JobIDRaw,JobName,State,ExitCode,Submit,Eligible,Start,End,ElapsedRaw,NodeList,NCPUS,CPUTimeRAW,TotalCPU,MaxRSS,MaxVMSize,AveRSS,ReqMem,AllocTRES' \
            'JobIDRaw,JobName,State,ExitCode,Submit,Start,End,ElapsedRaw,NodeList,NCPUS,CPUTimeRAW,TotalCPU,MaxRSS,ReqMem' \
            'JobIDRaw,State,ExitCode'
        do
            sacct_output="$(
                sacct -n -P -j "$job_id" \
                    --format="$candidate_fields" 2>/dev/null || true
            )"
            if [[ -n "$sacct_output" ]]; then
                sacct_fields="$candidate_fields"
                break
            fi
        done
        [[ -n "$sacct_output" ]] && break
        sleep 2
    done

    {
        printf '# fields=%s\n' "${sacct_fields:-unavailable}"
        printf '%s\n' "$sacct_output"
    } > "$WORK_DIR/${STAGE}.${round_tag}.sacct.tsv"
    {
        printf '# round=%s job_id=%s fields=%s\n' \
            "$round_history_label" "$job_id" "${sacct_fields:-unavailable}"
        printf '%s\n' "$sacct_output"
    } >> "$sacct_history"
    active_job_id=""

    next_attempt=$((attempt + 1))
    if (( next_attempt <= normal_attempt_limit )); then
        next_label="$(printf '%02d' "$next_attempt")"
        next_manifest="$WORK_DIR/${STAGE}.attempt${next_label}.pending.tsv"
    else
        next_fallback=$((next_attempt - normal_attempt_limit))
        next_label="$(printf '%02d' "$next_fallback")"
        next_manifest="$WORK_DIR/${STAGE}.fallback${next_label}.pending.tsv"
    fi
    : > "$next_manifest"
    failed_tasks=0

    while IFS=$'\t' read -r task_id task_hash command_b64 expected_b64; do
        [[ -z "$task_id" || "$task_id" == \#* ]] && continue

        if mtd_hpc_task_complete "$task_id" "$task_hash" "$expected_b64"; then
            continue
        fi

        failed_tasks=$((failed_tasks + 1))
        printf '%s\t%s\t%s\t%s\n' \
            "$task_id" "$task_hash" "$command_b64" "$expected_b64" \
            >> "$next_manifest"

        failure_host="$(mtd_hpc_failure_host "$task_id" "$task_hash" || true)"
        failure_reason="unknown"
        failure_code="unknown"
        failure_marker="$success_dir/${task_id}.failed"
        if [[ -s "$failure_marker" ]] && grep -Fxq "hash=$task_hash" "$failure_marker"; then
            failure_reason="$(sed -n 's/^reason=//p' "$failure_marker" | head -n 1)"
            failure_code="$(sed -n 's/^exit_code=//p' "$failure_marker" | head -n 1)"
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$round_history_label" "$job_id" "$task_id" "${failure_host:--}" \
            "${failure_code:-unknown}" "${failure_reason:-unknown}" \
            >> "$retry_history"

        if (( is_fallback == 0 )) && \
           [[ "$MTD_HPC_RETRY_EXCLUDE_FAILED_NODES" == "1" ]] && \
           [[ -n "$failure_host" ]]; then
            retry_excluded_nodes["$failure_host"]=1
        fi
    done < "$current_manifest"

    if (( failed_tasks == 0 )); then
        trap - INT TERM HUP
        if (( fallback_attempts_run > 0 )); then
            mtd_hpc_ok \
                "All $expected_tasks tasks completed after $normal_attempts_run normal attempt(s)" \
                "and $fallback_attempts_run submission-node fallback attempt(s)."
        else
            mtd_hpc_ok \
                "All $expected_tasks tasks completed successfully after" \
                "$normal_attempts_run attempt(s)."
        fi
        exit 0
    fi

    if (( is_fallback == 1 )); then
        mtd_hpc_error \
            "$failed_tasks task(s) failed or did not finish correctly in" \
            "submission-node fallback $fallback_attempt/$fallback_attempt_limit" \
            "on $submit_slurm_node."
    else
        mtd_hpc_error \
            "$failed_tasks task(s) failed or did not finish correctly in attempt" \
            "$attempt/$normal_attempt_limit."
    fi

    if (( attempt >= total_attempt_limit )); then
        if [[ "$MTD_HPC_FINAL_SUBMIT_NODE_FALLBACK" == "1" ]] && \
           (( fallback_attempt_limit == 0 )); then
            mtd_hpc_error \
                "Final submission-node fallback could not run:" \
                "$submit_fallback_unavailable_reason."
        fi
        break
    fi

    if (( attempt == normal_attempt_limit )) && \
       (( fallback_attempt_limit > 0 )); then
        mtd_hpc_info \
            "Normal attempts are exhausted; activating final submission-node fallback."
        mtd_hpc_info \
            "Only the $failed_tasks incomplete task(s) will be submitted to" \
            "$submit_slurm_node."
    fi

    if (( MTD_HPC_RETRY_DELAY_SECONDS > 0 )); then
        if (( attempt == normal_attempt_limit )) && \
           (( fallback_attempt_limit > 0 )); then
            mtd_hpc_info \
                "Starting the final fallback in $MTD_HPC_RETRY_DELAY_SECONDS second(s)."
        else
            mtd_hpc_info \
                "Retrying only the $failed_tasks incomplete task(s) in" \
                "$MTD_HPC_RETRY_DELAY_SECONDS second(s)."
        fi
        sleep "$MTD_HPC_RETRY_DELAY_SECONDS"
    else
        if (( attempt == normal_attempt_limit )) && \
           (( fallback_attempt_limit > 0 )); then
            mtd_hpc_info "Starting the final fallback immediately."
        else
            mtd_hpc_info "Retrying only the $failed_tasks incomplete task(s) immediately."
        fi
    fi

    current_manifest="$next_manifest"
    attempt="$next_attempt"
done

trap - INT TERM HUP

missing=0
while IFS=$'\t' read -r task_id task_hash _ expected_b64; do
    [[ -z "$task_id" || "$task_id" == \#* ]] && continue
    if ! mtd_hpc_task_complete "$task_id" "$task_hash" "$expected_b64"; then
        failure_marker="$success_dir/${task_id}.failed"
        if [[ -s "$failure_marker" ]] && grep -Fxq "hash=$task_hash" "$failure_marker"; then
            failure_host="$(sed -n 's/^host=//p' "$failure_marker" | head -n 1)"
            failure_code="$(sed -n 's/^exit_code=//p' "$failure_marker" | head -n 1)"
            failure_reason="$(sed -n 's/^reason=//p' "$failure_marker" | head -n 1)"
            mtd_hpc_error \
                "Task exhausted all attempts: $task_id "\
                "(host=${failure_host:-unknown}, exit=${failure_code:-unknown}, "\
                "reason=${failure_reason:-unknown})"
        else
            mtd_hpc_error "Task exhausted all attempts without a failure marker: $task_id"
        fi
        missing=$((missing + 1))
    fi
done < "$MANIFEST"

if (( fallback_attempts_run > 0 )); then
    mtd_hpc_error \
        "$missing task(s) still failed after $normal_attempts_run normal attempt(s)" \
        "and $fallback_attempts_run submission-node fallback attempt(s)."
else
    mtd_hpc_error \
        "$missing task(s) still failed after $normal_attempts_run normal attempt(s)."
fi
mtd_hpc_error "Logs: $log_root"
mtd_hpc_error "Retry history: $retry_history"
exit 1
