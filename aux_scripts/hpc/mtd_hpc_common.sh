#!/usr/bin/env bash
# Common functions for the independent MTD Explorer Slurm backend.

mtd_hpc_error() {
    printf '[MTD-HPC ERROR] %s\n' "$*" >&2
}

mtd_hpc_info() {
    printf '[MTD-HPC INFO] %s\n' "$*"
}

mtd_hpc_ok() {
    printf '[MTD-HPC OK] %s\n' "$*"
}

mtd_hpc_die() {
    mtd_hpc_error "$*"
    return 1
}

mtd_hpc_require_file() {
    local path="$1"
    local label="${2:-file}"
    [[ -s "$path" ]] || mtd_hpc_die "$label is missing or empty: $path"
}

mtd_hpc_require_dir() {
    local path="$1"
    local label="${2:-directory}"
    [[ -d "$path" ]] || mtd_hpc_die "$label does not exist: $path"
}

mtd_hpc_realpath() {
    readlink -f -- "$1"
}

mtd_hpc_validate_id() {
    local value="$1"
    local label="${2:-identifier}"
    [[ "$value" =~ ^[A-Za-z0-9._-]+$ ]] || \
        mtd_hpc_die "$label contains unsupported characters: $value"
}

mtd_hpc_load_config() {
    local conf="$1"
    local mtd_root="${2:-}"

    mtd_hpc_require_file "$conf" "HPC configuration" || return 1
    conf="$(mtd_hpc_realpath "$conf")" || return 1

    # This is a trusted Bash configuration file maintained by the user.
    # shellcheck disable=SC1090
    source "$conf"

    : "${MTD_HPC_SCHEDULER:=slurm}"
    : "${MTD_HPC_PREFIX:=/MTD_explorer_HPC}"
    : "${MTD_HPC_ENV_DIR:=${MTD_HPC_PREFIX}/envs/MTD-Explorer-HPC}"
    : "${MTD_HPC_CONDA_BIN:=${MTD_HPC_PREFIX}/miniconda3/bin/conda}"
    : "${MTD_HPC_DATABASE_ROOT:=${MTD_HPC_PREFIX}/databases}"
    : "${MTD_HPC_MTD_DATABASE_ROOT:=${MTD_HPC_DATABASE_ROOT}/MTD-Explorer}"
    : "${MTD_HPC_HUMANN_DB_ROOT:=${MTD_HPC_MTD_DATABASE_ROOT}/HUMAnN/ref_database}"
    : "${MTD_HPC_METAPHLAN_INDEX:=mpa_vJun23_CHOCOPhlAnSGB_202403}"
    : "${MTD_HPC_PARTITION:=}"
    : "${MTD_HPC_ACCOUNT:=}"
    : "${MTD_HPC_QOS:=}"
    : "${MTD_HPC_CONSTRAINT:=}"
    : "${MTD_HPC_TIME:=24:00:00}"
    : "${MTD_HPC_MAX_PARALLEL:=1}"
    : "${MTD_HPC_POLL_SECONDS:=20}"
    : "${MTD_HPC_MAGICBLAST_CHUNK_READS:=1000000}"
    : "${MTD_HPC_RESUME:=1}"
    : "${MTD_HPC_MAX_ATTEMPTS:=3}"
    : "${MTD_HPC_RETRY_DELAY_SECONDS:=20}"
    : "${MTD_HPC_RETRY_EXCLUDE_FAILED_NODES:=1}"
    : "${MTD_HPC_FINAL_SUBMIT_NODE_FALLBACK:=1}"
    : "${MTD_HPC_FINAL_SUBMIT_NODE_ATTEMPTS:=1}"

    # Node-local task staging prevents compute jobs from continuously reading
    # FASTQs or writing large intermediate outputs through the shared NFS.
    : "${MTD_HPC_STAGE_LOCAL:=1}"
    : "${MTD_HPC_LOCAL_SCRATCH_ROOT:=${MTD_HPC_PREFIX}/tmp}"
    : "${MTD_HPC_SERIALIZE_STAGEOUT_PER_NODE:=1}"
    : "${MTD_HPC_CLEAN_LOCAL_ON_SUCCESS:=1}"
    : "${MTD_HPC_CLEAN_LOCAL_ON_FAILURE:=0}"

    if ! declare -p MTD_HPC_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
        MTD_HPC_SBATCH_EXTRA_ARGS=(--exclusive --nodes=1 --mem=0)
    fi
    if ! declare -p MTD_HPC_HUMANN_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
        MTD_HPC_HUMANN_SBATCH_EXTRA_ARGS=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")
    fi
    if ! declare -p MTD_HPC_MAGICBLAST_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
        MTD_HPC_MAGICBLAST_SBATCH_EXTRA_ARGS=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")
    fi
    if ! declare -p MTD_HPC_PATH_MAPS >/dev/null 2>&1; then
        MTD_HPC_PATH_MAPS=()
    fi

    [[ "$MTD_HPC_SCHEDULER" == "slurm" ]] || \
        mtd_hpc_die "Only scheduler=slurm is currently supported. Got: $MTD_HPC_SCHEDULER" || return 1

    [[ "$MTD_HPC_MAX_PARALLEL" =~ ^[0-9]+$ ]] && (( MTD_HPC_MAX_PARALLEL >= 1 )) || \
        mtd_hpc_die "MTD_HPC_MAX_PARALLEL must be an integer >= 1." || return 1

    [[ "$MTD_HPC_POLL_SECONDS" =~ ^[0-9]+$ ]] && (( MTD_HPC_POLL_SECONDS >= 1 )) || \
        mtd_hpc_die "MTD_HPC_POLL_SECONDS must be an integer >= 1." || return 1

    [[ "$MTD_HPC_MAGICBLAST_CHUNK_READS" =~ ^[0-9]+$ ]] || \
        mtd_hpc_die "MTD_HPC_MAGICBLAST_CHUNK_READS must be an integer >= 0." || return 1

    [[ "$MTD_HPC_RESUME" == "0" || "$MTD_HPC_RESUME" == "1" ]] || \
        mtd_hpc_die "MTD_HPC_RESUME must be 0 or 1." || return 1

    [[ "$MTD_HPC_MAX_ATTEMPTS" =~ ^[0-9]+$ ]] && (( MTD_HPC_MAX_ATTEMPTS >= 1 )) || \
        mtd_hpc_die "MTD_HPC_MAX_ATTEMPTS must be an integer >= 1." || return 1

    [[ "$MTD_HPC_RETRY_DELAY_SECONDS" =~ ^[0-9]+$ ]] || \
        mtd_hpc_die "MTD_HPC_RETRY_DELAY_SECONDS must be an integer >= 0." || return 1

    [[ "$MTD_HPC_FINAL_SUBMIT_NODE_ATTEMPTS" =~ ^[0-9]+$ ]] && \
       (( MTD_HPC_FINAL_SUBMIT_NODE_ATTEMPTS >= 1 )) || \
        mtd_hpc_die "MTD_HPC_FINAL_SUBMIT_NODE_ATTEMPTS must be an integer >= 1." || return 1

    local boolean_name boolean_value
    for boolean_name in \
        MTD_HPC_STAGE_LOCAL \
        MTD_HPC_SERIALIZE_STAGEOUT_PER_NODE \
        MTD_HPC_CLEAN_LOCAL_ON_SUCCESS \
        MTD_HPC_CLEAN_LOCAL_ON_FAILURE \
        MTD_HPC_RETRY_EXCLUDE_FAILED_NODES \
        MTD_HPC_FINAL_SUBMIT_NODE_FALLBACK
    do
        boolean_value="${!boolean_name}"
        [[ "$boolean_value" == "0" || "$boolean_value" == "1" ]] || \
            mtd_hpc_die "$boolean_name must be 0 or 1." || return 1
    done

    [[ "$MTD_HPC_LOCAL_SCRATCH_ROOT" == /* ]] || \
        mtd_hpc_die "MTD_HPC_LOCAL_SCRATCH_ROOT must be an absolute path." || return 1

    MTD_HPC_CONF="$conf"
    MTD_HPC_MTD_ROOT="$mtd_root"
    export MTD_HPC_CONF MTD_HPC_MTD_ROOT
    export MTD_HPC_PREFIX MTD_HPC_ENV_DIR MTD_HPC_CONDA_BIN
    export MTD_HPC_DATABASE_ROOT MTD_HPC_MTD_DATABASE_ROOT
    export MTD_HPC_HUMANN_DB_ROOT MTD_HPC_METAPHLAN_INDEX
    export MTD_HPC_STAGE_LOCAL MTD_HPC_LOCAL_SCRATCH_ROOT
    export MTD_HPC_SERIALIZE_STAGEOUT_PER_NODE
    export MTD_HPC_CLEAN_LOCAL_ON_SUCCESS MTD_HPC_CLEAN_LOCAL_ON_FAILURE
    export MTD_HPC_MAX_ATTEMPTS MTD_HPC_RETRY_DELAY_SECONDS
    export MTD_HPC_RETRY_EXCLUDE_FAILED_NODES
    export MTD_HPC_FINAL_SUBMIT_NODE_FALLBACK
    export MTD_HPC_FINAL_SUBMIT_NODE_ATTEMPTS

    return 0
}

mtd_hpc_detect_node_threads() {
    local value=""

    for value in \
        "${SLURM_CPUS_ON_NODE:-}" \
        "${SLURM_CPUS_PER_TASK:-}"
    do
        # Some Slurm installations expose values such as 64(x2) or 64,64.
        value="${value%%(*}"
        value="${value%%,*}"
        if [[ "$value" =~ ^[0-9]+$ ]] && (( value >= 1 )); then
            printf '%s\n' "$value"
            return 0
        fi
    done

    nproc
}

mtd_hpc_detect_node_memory_kb() {
    awk '/^MemAvailable:/ {print $2; found=1; exit} END {if (!found) print 0}' /proc/meminfo
}

mtd_hpc_export_node_resources() {
    MTD_NODE_THREADS="$(mtd_hpc_detect_node_threads)"
    MTD_NODE_MEMORY_KB="$(mtd_hpc_detect_node_memory_kb)"

    [[ "$MTD_NODE_THREADS" =~ ^[0-9]+$ ]] && (( MTD_NODE_THREADS >= 1 )) || MTD_NODE_THREADS=1
    [[ "$MTD_NODE_MEMORY_KB" =~ ^[0-9]+$ ]] || MTD_NODE_MEMORY_KB=0

    MTD_NODE_MEMORY_KB_90=$(( MTD_NODE_MEMORY_KB * 90 / 100 ))
    MTD_NODE_MEMORY_GB_90=$(( MTD_NODE_MEMORY_KB_90 / 1024 / 1024 ))
    (( MTD_NODE_MEMORY_GB_90 >= 1 )) || MTD_NODE_MEMORY_GB_90=1

    export MTD_NODE_THREADS MTD_NODE_MEMORY_KB
    export MTD_NODE_MEMORY_KB_90 MTD_NODE_MEMORY_GB_90
    export OMP_NUM_THREADS="$MTD_NODE_THREADS"
    export OPENBLAS_NUM_THREADS="$MTD_NODE_THREADS"
    export MKL_NUM_THREADS="$MTD_NODE_THREADS"
    export NUMEXPR_NUM_THREADS="$MTD_NODE_THREADS"

    mtd_hpc_info "Node: $(hostname -s)"
    mtd_hpc_info "Detected threads: $MTD_NODE_THREADS"
    mtd_hpc_info "Available memory: ${MTD_NODE_MEMORY_KB} kB"
    mtd_hpc_info "90% memory budget: ${MTD_NODE_MEMORY_GB_90} GiB"
}


mtd_hpc_prepare_local_scratch() {
    local mode="$1"
    local sample="$2"
    local scratch_parent=""
    local scratch_dir=""
    local safe_user="${USER:-unknown}"
    local job_id="${SLURM_JOB_ID:-manual}"
    local task_id="${SLURM_ARRAY_TASK_ID:-0}"

    MTD_HPC_TASK_SCRATCH=""
    export MTD_HPC_TASK_SCRATCH

    [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]] || return 0

    mtd_hpc_validate_id "$mode" "scratch mode" || return 1
    mtd_hpc_validate_id "$sample" "scratch sample" || return 1
    mtd_hpc_validate_id "$safe_user" "scratch user" || safe_user="user"

    # Use the explicitly configured node-local prefix. This avoids relying
    # on SLURM_TMPDIR, which may be unset or backed by shared storage.
    scratch_parent="$MTD_HPC_LOCAL_SCRATCH_ROOT"

    [[ -d "$scratch_parent" ]] || \
        mtd_hpc_die "Node-local scratch directory does not exist: $scratch_parent" || return 1
    [[ -w "$scratch_parent" ]] || \
        mtd_hpc_die "Node-local scratch directory is not writable: $scratch_parent" || return 1

    scratch_dir="$(
        mktemp -d -p "$scratch_parent" \
            "mtd.${safe_user}.${job_id}.${task_id}.${mode}.${sample}.XXXXXX"
    )" || return 1

    chmod 700 "$scratch_dir"
    mkdir -p -- "$scratch_dir/input" "$scratch_dir/output" "$scratch_dir/tmp"

    export TMPDIR="$scratch_dir/tmp"
    export TMP="$TMPDIR"
    export TEMP="$TMPDIR"

    MTD_HPC_TASK_SCRATCH="$scratch_dir"
    export MTD_HPC_TASK_SCRATCH

    mtd_hpc_info "Node-local scratch: $MTD_HPC_TASK_SCRATCH"
}

mtd_hpc_stage_in_file() {
    local source="$1"
    local destination="$2"
    local label="${3:-input}"
    local source_size=0
    local destination_size=0

    mtd_hpc_require_file "$source" "$label" || return 1
    command -v rsync >/dev/null 2>&1 || \
        mtd_hpc_die "rsync is required for node-local input staging." || return 1

    mkdir -p -- "$(dirname -- "$destination")"

    source_size="$(stat -Lc %s -- "$source")"
    mtd_hpc_info "Stage-in $label: $source -> $destination (${source_size} bytes)"

    rsync -aL --whole-file -- "$source" "$destination" || return 1
    mtd_hpc_require_file "$destination" "staged $label" || return 1

    destination_size="$(stat -Lc %s -- "$destination")"
    [[ "$destination_size" == "$source_size" ]] || \
        mtd_hpc_die "Stage-in size mismatch for $label: source=$source_size destination=$destination_size" || return 1
}

mtd_hpc_stageout_lock_acquire() {
    MTD_HPC_STAGEOUT_LOCK_FD=""

    [[ "$MTD_HPC_SERIALIZE_STAGEOUT_PER_NODE" == "1" ]] || return 0
    command -v flock >/dev/null 2>&1 || \
        mtd_hpc_die "flock is required when MTD_HPC_SERIALIZE_STAGEOUT_PER_NODE=1." || return 1

    mkdir -p -- "$MTD_HPC_LOCAL_SCRATCH_ROOT/.locks"
    exec {MTD_HPC_STAGEOUT_LOCK_FD}> \
        "$MTD_HPC_LOCAL_SCRATCH_ROOT/.locks/stageout.lock"
    flock -x "$MTD_HPC_STAGEOUT_LOCK_FD"

    mtd_hpc_info "Acquired node-local stage-out lock."
}

mtd_hpc_stageout_lock_release() {
    [[ -n "${MTD_HPC_STAGEOUT_LOCK_FD:-}" ]] || return 0

    flock -u "$MTD_HPC_STAGEOUT_LOCK_FD" || true
    eval "exec ${MTD_HPC_STAGEOUT_LOCK_FD}>&-"
    MTD_HPC_STAGEOUT_LOCK_FD=""

    mtd_hpc_info "Released node-local stage-out lock."
}

mtd_hpc_atomic_stage_out_file() {
    local source="$1"
    local destination="$2"
    local label="${3:-output}"
    local destination_dir=""
    local temporary_destination=""
    local source_size=0
    local copied_size=0
    local rc=0

    mtd_hpc_require_file "$source" "$label" || return 1
    command -v rsync >/dev/null 2>&1 || \
        mtd_hpc_die "rsync is required for node-local output staging." || return 1

    destination_dir="$(dirname -- "$destination")"
    mkdir -p -- "$destination_dir"

    temporary_destination="$destination.mtd-partial.${SLURM_JOB_ID:-$$}.${SLURM_ARRAY_TASK_ID:-0}.$$"
    rm -f -- "$temporary_destination"

    source_size="$(stat -Lc %s -- "$source")"
    mtd_hpc_info "Stage-out $label: $source -> $destination (${source_size} bytes)"

    mtd_hpc_stageout_lock_acquire || return 1

    if rsync -a --whole-file -- "$source" "$temporary_destination"; then
        if [[ ! -s "$temporary_destination" ]]; then
            mtd_hpc_error "Staged $label is missing or empty: $temporary_destination"
            rc=1
        else
            copied_size="$(stat -Lc %s -- "$temporary_destination")"
            if [[ "$copied_size" != "$source_size" ]]; then
                mtd_hpc_error \
                    "Stage-out size mismatch for $label: source=$source_size destination=$copied_size"
                rc=1
            elif mv -f -- "$temporary_destination" "$destination"; then
                rc=0
            else
                rc=$?
            fi
        fi
    else
        rc=$?
    fi

    [[ "$rc" == "0" ]] || rm -f -- "$temporary_destination"
    mtd_hpc_stageout_lock_release

    (( rc == 0 )) || return "$rc"
    mtd_hpc_require_file "$destination" "staged $label" || return 1
    mtd_hpc_ok "Stage-out completed: $destination"
}

mtd_hpc_cleanup_local_scratch() {
    local scratch_dir="${1:-}"
    local exit_status="${2:-1}"
    local clean=0

    [[ -n "$scratch_dir" && -d "$scratch_dir" ]] || return 0

    if (( exit_status == 0 )); then
        clean="$MTD_HPC_CLEAN_LOCAL_ON_SUCCESS"
    else
        clean="$MTD_HPC_CLEAN_LOCAL_ON_FAILURE"
    fi

    if [[ "$clean" == "1" ]]; then
        rm -rf -- "$scratch_dir"
        mtd_hpc_info "Removed node-local scratch: $scratch_dir"
    else
        mtd_hpc_error "Preserving node-local scratch for diagnosis: $scratch_dir"
    fi
}

mtd_hpc_map_path_to_node() {
    local source_path="$1"
    local mtd_root="${2:-${MTD_HPC_MTD_ROOT:-}}"
    local best_source=""
    local best_target=""
    local mapping src dst

    if [[ -n "$mtd_root" && ( "$source_path" == "$mtd_root" || "$source_path" == "$mtd_root"/* ) ]]; then
        best_source="$mtd_root"
        best_target="$MTD_HPC_MTD_DATABASE_ROOT"
    fi

    for mapping in "${MTD_HPC_PATH_MAPS[@]}"; do
        [[ "$mapping" == *=* ]] || continue
        src="${mapping%%=*}"
        dst="${mapping#*=}"
        [[ -n "$src" && -n "$dst" ]] || continue

        if [[ "$source_path" == "$src" || "$source_path" == "$src"/* ]]; then
            if (( ${#src} > ${#best_source} )); then
                best_source="$src"
                best_target="$dst"
            fi
        fi
    done

    if [[ -z "$best_source" ]]; then
        mtd_hpc_error "No node-local database mapping matches: $source_path"
        return 1
    fi

    printf '%s%s\n' "$best_target" "${source_path#"$best_source"}"
}

mtd_hpc_b64_encode() {
    printf '%s' "$1" | base64 | tr -d '\n'
}

mtd_hpc_b64_decode() {
    printf '%s' "$1" | base64 --decode
}

mtd_hpc_outputs_exist() {
    local outputs="$1"
    local expected

    [[ -n "$outputs" ]] || return 0
    while IFS= read -r expected; do
        [[ -z "$expected" ]] && continue
        [[ -s "$expected" ]] || return 1
    done <<< "$outputs"
    return 0
}

mtd_hpc_validate_humann_databases() {
    local root="${1:-$MTD_HPC_HUMANN_DB_ROOT}"
    local index="${2:-$MTD_HPC_METAPHLAN_INDEX}"
    local utility_count=0
    local suffix

    mtd_hpc_require_dir "$root/chocophlan" "HUMAnN ChocoPhlAn database" || return 1
    mtd_hpc_require_dir "$root/uniref" "HUMAnN UniRef database" || return 1
    mtd_hpc_require_dir "$root/utility_mapping" "HUMAnN utility mapping database" || return 1
    mtd_hpc_require_dir "$root/metaphlan" "MetaPhlAn database" || return 1

    if ! find "$root/chocophlan" -type f -name '*.ffn.gz' -size +0c -print -quit 2>/dev/null | grep -q .; then
        mtd_hpc_die "No usable ChocoPhlAn .ffn.gz files were found under: $root/chocophlan"
        return 1
    fi

    mtd_hpc_require_file         "$root/uniref/uniref90_201901b_full.dmnd"         "HUMAnN UniRef90 DIAMOND database" || return 1

    utility_count="$(find "$root/utility_mapping" -type f -size +0c 2>/dev/null | awk 'END {print NR + 0}')"
    if ! [[ "$utility_count" =~ ^[0-9]+$ ]] || (( utility_count < 10 )); then
        mtd_hpc_die "HUMAnN utility-mapping database appears incomplete: $root/utility_mapping"
        return 1
    fi

    for suffix in pkl 1.bt2l 2.bt2l 3.bt2l 4.bt2l rev.1.bt2l rev.2.bt2l; do
        mtd_hpc_require_file             "$root/metaphlan/${index}.${suffix}"             "MetaPhlAn database component" || return 1
    done
}

mtd_hpc_task_line() {
    local task_id="$1"
    local command="$2"
    local expected_outputs="$3"
    local hash command_b64 expected_b64

    mtd_hpc_validate_id "$task_id" "task ID" || return 1
    hash="$(printf '%s\0%s' "$command" "$expected_outputs" | sha256sum | awk '{print $1}')"
    command_b64="$(mtd_hpc_b64_encode "$command")"
    expected_b64="$(mtd_hpc_b64_encode "$expected_outputs")"

    printf '%s\t%s\t%s\t%s\n' "$task_id" "$hash" "$command_b64" "$expected_b64"
}
