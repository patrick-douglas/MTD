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
    : "${MTD_HPC_FASTP_ENV_DIR:=${MTD_HPC_PREFIX}/envs/MTD_fastp}"
    : "${MTD_HPC_KRAKEN2_ENV_DIR:=${MTD_HPC_PREFIX}/envs/MTD_kraken2}"
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
    : "${MTD_HPC_SCRATCH_RESERVE_GB:=10}"
    : "${MTD_HPC_FASTP_SCRATCH_MULTIPLIER:=3}"
    : "${MTD_HPC_KRAKEN_SCRATCH_MULTIPLIER:=3}"
    : "${MTD_HPC_REMOTE_INPUT_FROM_SUBMIT_HOST:=1}"
    : "${MTD_HPC_SUBMIT_SSH_CONNECT_TIMEOUT:=10}"

    # Cluster-wide transfer throttling. The shared lock directory itself is
    # supplied by mtd_hpc_submit_array.sh for each submitted stage.
    # A value of 0 disables the corresponding global throttle.
    : "${MTD_HPC_STAGEIN_MAX_CONCURRENT:=1}"
    : "${MTD_HPC_STAGEOUT_MAX_CONCURRENT:=1}"
    : "${MTD_HPC_TRANSFER_LOCK_POLL_SECONDS:=1}"
    : "${MTD_HPC_TRANSFER_LOCK_STALE_SECONDS:=300}"

    if ! declare -p MTD_HPC_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
        MTD_HPC_SBATCH_EXTRA_ARGS=(--exclusive --nodes=1 --mem=0)
    fi
    if ! declare -p MTD_HPC_HUMANN_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
        MTD_HPC_HUMANN_SBATCH_EXTRA_ARGS=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")
    fi
    if ! declare -p MTD_HPC_MAGICBLAST_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
        MTD_HPC_MAGICBLAST_SBATCH_EXTRA_ARGS=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")
    fi
    if ! declare -p MTD_HPC_FASTP_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
        MTD_HPC_FASTP_SBATCH_EXTRA_ARGS=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")
    fi
    if ! declare -p MTD_HPC_KRAKEN_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
        MTD_HPC_KRAKEN_SBATCH_EXTRA_ARGS=("${MTD_HPC_SBATCH_EXTRA_ARGS[@]}")
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

    [[ "$MTD_HPC_SCRATCH_RESERVE_GB" =~ ^[0-9]+$ ]] || \
        mtd_hpc_die "MTD_HPC_SCRATCH_RESERVE_GB must be an integer >= 0." || return 1

    local multiplier_name multiplier_value
    for multiplier_name in \
        MTD_HPC_FASTP_SCRATCH_MULTIPLIER \
        MTD_HPC_KRAKEN_SCRATCH_MULTIPLIER
    do
        multiplier_value="${!multiplier_name}"
        [[ "$multiplier_value" =~ ^[0-9]+$ ]] && (( multiplier_value >= 1 )) || \
            mtd_hpc_die "$multiplier_name must be an integer >= 1." || return 1
    done

    local boolean_name boolean_value
    for boolean_name in \
        MTD_HPC_STAGE_LOCAL \
        MTD_HPC_SERIALIZE_STAGEOUT_PER_NODE \
        MTD_HPC_CLEAN_LOCAL_ON_SUCCESS \
        MTD_HPC_CLEAN_LOCAL_ON_FAILURE \
        MTD_HPC_RETRY_EXCLUDE_FAILED_NODES \
        MTD_HPC_FINAL_SUBMIT_NODE_FALLBACK \
        MTD_HPC_REMOTE_INPUT_FROM_SUBMIT_HOST
    do
        boolean_value="${!boolean_name}"
        [[ "$boolean_value" == "0" || "$boolean_value" == "1" ]] || \
            mtd_hpc_die "$boolean_name must be 0 or 1." || return 1
    done

    [[ "$MTD_HPC_SUBMIT_SSH_CONNECT_TIMEOUT" =~ ^[0-9]+$ ]] && \
       (( MTD_HPC_SUBMIT_SSH_CONNECT_TIMEOUT >= 1 )) || \
        mtd_hpc_die "MTD_HPC_SUBMIT_SSH_CONNECT_TIMEOUT must be an integer >= 1." || return 1

    local transfer_limit_name transfer_limit_value
    for transfer_limit_name in \
        MTD_HPC_STAGEIN_MAX_CONCURRENT \
        MTD_HPC_STAGEOUT_MAX_CONCURRENT
    do
        transfer_limit_value="${!transfer_limit_name}"
        [[ "$transfer_limit_value" =~ ^[0-9]+$ ]] || \
            mtd_hpc_die "$transfer_limit_name must be an integer >= 0." || return 1
    done

    [[ "$MTD_HPC_TRANSFER_LOCK_POLL_SECONDS" =~ ^[0-9]+$ ]] && \
       (( MTD_HPC_TRANSFER_LOCK_POLL_SECONDS >= 1 )) || \
        mtd_hpc_die "MTD_HPC_TRANSFER_LOCK_POLL_SECONDS must be an integer >= 1." || return 1

    [[ "$MTD_HPC_TRANSFER_LOCK_STALE_SECONDS" =~ ^[0-9]+$ ]] && \
       (( MTD_HPC_TRANSFER_LOCK_STALE_SECONDS >= 1 )) || \
        mtd_hpc_die "MTD_HPC_TRANSFER_LOCK_STALE_SECONDS must be an integer >= 1." || return 1

    [[ "$MTD_HPC_LOCAL_SCRATCH_ROOT" == /* ]] || \
        mtd_hpc_die "MTD_HPC_LOCAL_SCRATCH_ROOT must be an absolute path." || return 1

    MTD_HPC_CONF="$conf"
    MTD_HPC_MTD_ROOT="$mtd_root"
    export MTD_HPC_CONF MTD_HPC_MTD_ROOT
    export MTD_HPC_PREFIX MTD_HPC_ENV_DIR MTD_HPC_FASTP_ENV_DIR
    export MTD_HPC_KRAKEN2_ENV_DIR MTD_HPC_CONDA_BIN
    export MTD_HPC_DATABASE_ROOT MTD_HPC_MTD_DATABASE_ROOT
    export MTD_HPC_HUMANN_DB_ROOT MTD_HPC_METAPHLAN_INDEX
    export MTD_HPC_STAGE_LOCAL MTD_HPC_LOCAL_SCRATCH_ROOT
    export MTD_HPC_SERIALIZE_STAGEOUT_PER_NODE
    export MTD_HPC_CLEAN_LOCAL_ON_SUCCESS MTD_HPC_CLEAN_LOCAL_ON_FAILURE
    export MTD_HPC_SCRATCH_RESERVE_GB
    export MTD_HPC_FASTP_SCRATCH_MULTIPLIER MTD_HPC_KRAKEN_SCRATCH_MULTIPLIER
    export MTD_HPC_REMOTE_INPUT_FROM_SUBMIT_HOST MTD_HPC_SUBMIT_SSH_CONNECT_TIMEOUT
    export MTD_HPC_STAGEIN_MAX_CONCURRENT MTD_HPC_STAGEOUT_MAX_CONCURRENT
    export MTD_HPC_TRANSFER_LOCK_POLL_SECONDS MTD_HPC_TRANSFER_LOCK_STALE_SECONDS
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


mtd_hpc_require_path_exists() {
    local path="$1"
    local label="${2:-path}"
    [[ -e "$path" ]] || mtd_hpc_die "$label does not exist: $path"
}

mtd_hpc_submit_ssh_candidates() {
    local candidate=""
    local -A seen=()

    for candidate in \
        "${MTD_HPC_SUBMIT_HOST_FQDN:-}" \
        "${MTD_HPC_SUBMIT_HOST_SHORT:-}" \
        "${MTD_HPC_SUBMIT_SLURM_NODE:-}"
    do
        [[ -n "$candidate" && "$candidate" != "unknown" ]] || continue
        [[ "$candidate" =~ ^[A-Za-z0-9._:-]+$ ]] || continue
        [[ -z "${seen[$candidate]:-}" ]] || continue
        seen["$candidate"]=1
        printf '%s\n' "$candidate"
    done
}

mtd_hpc_resolve_submit_ssh_target() {
    local submit_user="${MTD_HPC_SUBMIT_USER:-}"
    local candidate=""
    local target=""

    [[ "$MTD_HPC_REMOTE_INPUT_FROM_SUBMIT_HOST" == "1" ]] || {
        mtd_hpc_error "Remote input staging from the submission host is disabled."
        return 1
    }

    [[ "$submit_user" =~ ^[A-Za-z0-9._-]+$ ]] || {
        mtd_hpc_error \
            "Submission user is unavailable or invalid for remote input staging: ${submit_user:-<unset>}"
        return 1
    }

    command -v ssh >/dev/null 2>&1 || {
        mtd_hpc_error \
            "ssh is required when an input exists only on the host that launched MTD Explorer."
        return 1
    }

    while IFS= read -r candidate; do
        [[ -n "$candidate" ]] || continue
        target="${submit_user}@${candidate}"
        if ssh \
            -o BatchMode=yes \
            -o ConnectTimeout="$MTD_HPC_SUBMIT_SSH_CONNECT_TIMEOUT" \
            -- "$target" true >/dev/null 2>&1
        then
            printf '%s\n' "$target"
            return 0
        fi
    done < <(mtd_hpc_submit_ssh_candidates)

    mtd_hpc_error \
        "Input is not mounted on node ${SLURMD_NODENAME:-$(hostname -s)} and no "\
        "passwordless SSH connection to the host that launched MTD Explorer could be established."
    return 1
}

mtd_hpc_require_local_scratch_capacity_bytes() {
    local multiplier="$1"
    local label="$2"
    shift 2

    local total_input_bytes=0
    local required_bytes=0
    local reserve_bytes=0
    local free_kb=0
    local free_bytes=0
    local size=""

    [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]] || return 0
    [[ "$multiplier" =~ ^[0-9]+$ ]] && (( multiplier >= 1 )) || \
        mtd_hpc_die "Invalid scratch multiplier for $label: $multiplier" || return 1

    for size in "$@"; do
        [[ -z "$size" || "$size" == "-" ]] && continue
        [[ "$size" =~ ^[0-9]+$ ]] && (( size > 0 )) || \
            mtd_hpc_die "Invalid input size for $label scratch precheck: $size" || return 1
        total_input_bytes=$((total_input_bytes + size))
    done

    free_kb="$(df -Pk "$MTD_HPC_LOCAL_SCRATCH_ROOT" | awk 'NR == 2 {print $4}')"
    [[ "$free_kb" =~ ^[0-9]+$ ]] || \
        mtd_hpc_die "Could not determine free scratch space: $MTD_HPC_LOCAL_SCRATCH_ROOT" || return 1

    free_bytes=$((free_kb * 1024))
    reserve_bytes=$((MTD_HPC_SCRATCH_RESERVE_GB * 1024 * 1024 * 1024))
    required_bytes=$((total_input_bytes * multiplier + reserve_bytes))

    mtd_hpc_info \
        "$label scratch precheck: input=${total_input_bytes} bytes, multiplier=${multiplier}," \
        "reserve=${MTD_HPC_SCRATCH_RESERVE_GB} GiB, required=${required_bytes} bytes," \
        "free=${free_bytes} bytes"

    if (( free_bytes < required_bytes )); then
        mtd_hpc_die \
            "Insufficient node-local scratch for $label. Required ${required_bytes} bytes;" \
            "available ${free_bytes} bytes at $MTD_HPC_LOCAL_SCRATCH_ROOT."
        return 1
    fi
}

mtd_hpc_require_local_scratch_capacity() {
    local multiplier="$1"
    local label="$2"
    shift 2

    local total_input_bytes=0
    local required_bytes=0
    local reserve_bytes=0
    local free_kb=0
    local free_bytes=0
    local path=""
    local size=0

    [[ "$MTD_HPC_STAGE_LOCAL" == "1" ]] || return 0
    [[ "$multiplier" =~ ^[0-9]+$ ]] && (( multiplier >= 1 )) || \
        mtd_hpc_die "Invalid scratch multiplier for $label: $multiplier" || return 1

    for path in "$@"; do
        [[ -z "$path" || "$path" == "-" ]] && continue
        mtd_hpc_require_file "$path" "$label input" || return 1
        size="$(stat -Lc %s -- "$path")" || return 1
        total_input_bytes=$((total_input_bytes + size))
    done

    free_kb="$(df -Pk "$MTD_HPC_LOCAL_SCRATCH_ROOT" | awk 'NR == 2 {print $4}')"
    [[ "$free_kb" =~ ^[0-9]+$ ]] || \
        mtd_hpc_die "Could not determine free scratch space: $MTD_HPC_LOCAL_SCRATCH_ROOT" || return 1

    free_bytes=$((free_kb * 1024))
    reserve_bytes=$((MTD_HPC_SCRATCH_RESERVE_GB * 1024 * 1024 * 1024))
    required_bytes=$((total_input_bytes * multiplier + reserve_bytes))

    mtd_hpc_info \
        "$label scratch precheck: input=${total_input_bytes} bytes, multiplier=${multiplier}," \
        "reserve=${MTD_HPC_SCRATCH_RESERVE_GB} GiB, required=${required_bytes} bytes," \
        "free=${free_bytes} bytes"

    if (( free_bytes < required_bytes )); then
        mtd_hpc_die \
            "Insufficient node-local scratch for $label. Required ${required_bytes} bytes;" \
            "available ${free_bytes} bytes at $MTD_HPC_LOCAL_SCRATCH_ROOT."
        return 1
    fi
}

mtd_hpc_validate_fastq_pair_ids() {
    local read1="$1"
    local read2="$2"
    local label="${3:-paired FASTQ}"

    local python_bin="$MTD_HPC_ENV_DIR/bin/python"
    [[ -x "$python_bin" ]] || \
        mtd_hpc_die "Python executable not found in the node-local HPC environment: $python_bin" || return 1

    "$python_bin" - "$read1" "$read2" "$label" <<'PY_FASTQ_ID_CHECK'
import gzip
import sys

read1, read2, label = sys.argv[1:4]

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")

def normalize(header):
    value = header.strip()
    if value.startswith("@"):
        value = value[1:]
    value = value.split()[0]
    if value.endswith("/1") or value.endswith("/2"):
        value = value[:-2]
    return value

def record(handle):
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not quality:
        raise RuntimeError("Incomplete FASTQ record near header: " + header.strip())
    if not header.startswith("@") or not plus.startswith("+"):
        raise RuntimeError("Malformed FASTQ record near header: " + header.strip())
    return header, sequence, plus, quality

checked = 0
mismatches = []
try:
    with open_maybe_gz(read1) as first, open_maybe_gz(read2) as second:
        while True:
            r1 = record(first)
            r2 = record(second)
            if r1 is None and r2 is None:
                break
            if r1 is None or r2 is None:
                raise RuntimeError("Paired FASTQ files have different record counts")
            checked += 1
            id1 = normalize(r1[0])
            id2 = normalize(r2[0])
            if id1 != id2 and len(mismatches) < 20:
                mismatches.append((checked, id1, id2))
except Exception as error:
    print(f"[MTD-HPC ERROR] Could not validate {label}: {error}", file=sys.stderr)
    sys.exit(2)

if mismatches:
    print(f"[MTD-HPC ERROR] Paired FASTQ IDs are desynchronized: {label}", file=sys.stderr)
    for index, id1, id2 in mismatches:
        print(f"{index}\tR1={id1}\tR2={id2}", file=sys.stderr)
    sys.exit(1)

print(f"[MTD-HPC OK] Paired FASTQ IDs validated: {label} ({checked} pairs)")
PY_FASTQ_ID_CHECK
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

mtd_hpc_global_transfer_slot_owner_query_id() {
    if [[ -n "${SLURM_ARRAY_JOB_ID:-}" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        printf '%s_%s\n' "$SLURM_ARRAY_JOB_ID" "$SLURM_ARRAY_TASK_ID"
    else
        printf '%s\n' "${SLURM_JOB_ID:-}"
    fi
}

mtd_hpc_global_transfer_slot_write_owner() {
    local slot_dir="$1"
    local token="$2"
    local owner_tmp="$slot_dir/.owner.$$.${RANDOM}"
    local query_id=""

    query_id="$(mtd_hpc_global_transfer_slot_owner_query_id)"

    {
        printf 'token=%s\n' "$token"
        printf 'host=%s\n' "${SLURMD_NODENAME:-$(hostname -s)}"
        printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID:-}"
        printf 'slurm_array_job_id=%s\n' "${SLURM_ARRAY_JOB_ID:-}"
        printf 'slurm_array_task_id=%s\n' "${SLURM_ARRAY_TASK_ID:-}"
        printf 'slurm_query_id=%s\n' "$query_id"
        printf 'pid=%s\n' "$$"
        printf 'created_epoch=%s\n' "$(date +%s)"
    } > "$owner_tmp" || return 1

    mv -f -- "$owner_tmp" "$slot_dir/owner" || {
        rm -f -- "$owner_tmp"
        return 1
    }
}

mtd_hpc_global_transfer_slot_owner_is_active() {
    local slot_dir="$1"
    local owner_file="$slot_dir/owner"
    local query_id=""
    local state=""

    [[ -r "$owner_file" ]] || return 1

    query_id="$(
        sed -n 's/^slurm_query_id=//p' "$owner_file" 2>/dev/null |
            head -n 1
    )"
    [[ -n "$query_id" ]] || return 1

    if command -v squeue >/dev/null 2>&1; then
        if state="$(squeue -h -j "$query_id" -o '%T' 2>/dev/null)"; then
            [[ -n "$state" ]] && return 0
            return 1
        fi

        # A scheduler query failure is not proof that the owner disappeared.
        # Keep the lock rather than risking two transfers at once.
        return 0
    fi

    if command -v scontrol >/dev/null 2>&1; then
        if scontrol show job "$query_id" >/dev/null 2>&1; then
            return 0
        fi
        return 1
    fi

    # Without a scheduler client on the worker, stale ownership cannot be
    # verified safely. Leave the lock in place rather than stealing it.
    return 0
}

mtd_hpc_global_transfer_slot_reclaim_if_stale() {
    local slot_dir="$1"
    local direction="$2"
    local slot="$3"
    local display_direction="$direction"
    local now=0
    local modified=0
    local age=0
    local quarantine=""

    [[ -d "$slot_dir" ]] || return 1

    case "$direction" in
        stagein) display_direction="stage-in" ;;
        stageout) display_direction="stage-out" ;;
    esac

    modified="$(stat -Lc %Y -- "$slot_dir" 2>/dev/null)" || return 1
    now="$(date +%s)"
    age=$((now - modified))
    (( age >= MTD_HPC_TRANSFER_LOCK_STALE_SECONDS )) || return 1

    if mtd_hpc_global_transfer_slot_owner_is_active "$slot_dir"; then
        # Refresh the directory timestamp after a successful owner check so
        # long transfers do not trigger a scheduler query on every poll once
        # they exceed the stale-age threshold.
        touch -m -- "$slot_dir" >/dev/null 2>&1 || true
        return 1
    fi

    # Rename first instead of deleting the shared path in-place. On NFS this
    # makes stale-lock recovery race-safe: only one waiter can move the old
    # directory, and a newly acquired slot at the original pathname is never
    # removed by the recovery cleanup.
    quarantine="${slot_dir}.reclaimed.${SLURM_JOB_ID:-manual}.${SLURM_ARRAY_TASK_ID:-0}.$$.$RANDOM"
    if mv -- "$slot_dir" "$quarantine" 2>/dev/null; then
        mtd_hpc_info \
            "Recovered stale shared $display_direction transfer slot $slot" \
            "after ${age}s without an active Slurm owner."
        rm -rf -- "$quarantine"
        return 0
    fi

    return 1
}

mtd_hpc_global_transfer_slot_acquire() {
    local direction="$1"
    local max_concurrent="$2"
    local handle_variable="$3"
    local slot_variable="$4"
    local lock_dir="${MTD_HPC_TRANSFER_LOCK_DIR:-}"
    local display_direction="$direction"
    local slot=0
    local slot_dir=""
    local token=""
    local handle=""
    local waited_seconds=0
    local announced_wait=0
    local reclaimed_any=0

    printf -v "$handle_variable" '%s' ""
    printf -v "$slot_variable" '%s' ""

    case "$direction" in
        stagein) display_direction="stage-in" ;;
        stageout) display_direction="stage-out" ;;
    esac

    (( max_concurrent > 0 )) || return 0

    # Direct/manual node-job invocations may not have been launched through the
    # array submitter. Keep those useful for smoke tests instead of failing.
    if [[ -z "$lock_dir" ]]; then
        mtd_hpc_info \
            "Shared $display_direction transfer throttling is unavailable for this standalone task" \
            "(MTD_HPC_TRANSFER_LOCK_DIR is not set)."
        return 0
    fi

    [[ "$lock_dir" == /* ]] || \
        mtd_hpc_die "MTD_HPC_TRANSFER_LOCK_DIR must be an absolute shared path: $lock_dir" || return 1

    mkdir -p -- "$lock_dir/$direction" || return 1
    [[ -w "$lock_dir/$direction" ]] || \
        mtd_hpc_die "Shared $display_direction lock directory is not writable: $lock_dir/$direction" || return 1

    while true; do
        reclaimed_any=0
        for (( slot=1; slot<=max_concurrent; slot++ )); do
            # Use mkdir rather than flock for the cluster-wide semaphore.
            # mkdir is an atomic namespace operation on the shared filesystem
            # and remains mutually exclusive when NFS clients are mounted with
            # nolock/local_lock=all, where flock is only client-local.
            slot_dir="$lock_dir/$direction/slot_${slot}.lockdir"

            if mkdir -- "$slot_dir" 2>/dev/null; then
                token="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}.${SLURM_ARRAY_TASK_ID:-0}.${SLURM_JOB_ID:-0}.$$.$RANDOM"
                if ! mtd_hpc_global_transfer_slot_write_owner "$slot_dir" "$token"; then
                    rm -rf -- "$slot_dir"
                    mtd_hpc_die "Unable to write shared $display_direction transfer lock metadata: $slot_dir" || return 1
                fi

                handle="${slot_dir}"$'\t'"${token}"
                printf -v "$handle_variable" '%s' "$handle"
                printf -v "$slot_variable" '%s' "$slot"

                if (( waited_seconds > 0 )); then
                    mtd_hpc_info \
                        "Acquired shared $display_direction transfer slot $slot/$max_concurrent" \
                        "after waiting ${waited_seconds}s."
                else
                    mtd_hpc_info \
                        "Acquired shared $display_direction transfer slot $slot/$max_concurrent."
                fi
                return 0
            fi

            if mtd_hpc_global_transfer_slot_reclaim_if_stale \
                "$slot_dir" "$direction" "$slot"
            then
                reclaimed_any=1
            fi
        done

        # A stale slot was just removed. Retry immediately instead of adding
        # an unnecessary polling delay before claiming the now-free slot.
        (( reclaimed_any == 0 )) || continue

        if (( announced_wait == 0 )); then
            mtd_hpc_info \
                "Waiting for a shared $display_direction transfer slot" \
                "(limit: $max_concurrent concurrent transfer(s))..."
            announced_wait=1
        fi

        sleep "$MTD_HPC_TRANSFER_LOCK_POLL_SECONDS"
        waited_seconds=$((waited_seconds + MTD_HPC_TRANSFER_LOCK_POLL_SECONDS))
    done
}

mtd_hpc_global_transfer_slot_release() {
    local direction="$1"
    local handle_variable="$2"
    local slot_variable="$3"
    local handle="${!handle_variable:-}"
    local slot="${!slot_variable:-}"
    local display_direction="$direction"
    local slot_dir=""
    local token=""
    local current_token=""
    local released=0

    case "$direction" in
        stagein) display_direction="stage-in" ;;
        stageout) display_direction="stage-out" ;;
    esac

    [[ -n "$handle" ]] || return 0

    slot_dir="${handle%%$'\t'*}"
    token="${handle#*$'\t'}"

    if [[ -r "$slot_dir/owner" ]]; then
        current_token="$(
            sed -n 's/^token=//p' "$slot_dir/owner" 2>/dev/null |
                head -n 1
        )"
    fi

    if [[ -d "$slot_dir" && "$current_token" == "$token" ]]; then
        rm -f -- "$slot_dir/owner"
        if rmdir -- "$slot_dir" 2>/dev/null; then
            released=1
        else
            mtd_hpc_error \
                "Unable to remove shared $display_direction transfer lock directory: $slot_dir"
        fi
    elif [[ -d "$slot_dir" ]]; then
        mtd_hpc_info \
            "Shared $display_direction transfer slot${slot:+ $slot} ownership changed;" \
            "leaving the current lock untouched."
    fi

    printf -v "$handle_variable" '%s' ""
    printf -v "$slot_variable" '%s' ""

    if (( released == 1 )); then
        mtd_hpc_info "Released shared $display_direction transfer slot${slot:+ $slot}."
    fi
}

mtd_hpc_stagein_group_begin() {
    # A group keeps one global slot across all input files belonging to the
    # same task (for example paired R1/R2), so a sample can begin computation
    # as soon as its complete input set has arrived.
    MTD_HPC_STAGEIN_GROUP_ACTIVE=1
    MTD_HPC_STAGEIN_GLOBAL_LOCK_HANDLE=""
    MTD_HPC_STAGEIN_GLOBAL_LOCK_SLOT=""
}

mtd_hpc_stagein_lock_acquire() {
    [[ -n "${MTD_HPC_STAGEIN_GLOBAL_LOCK_HANDLE:-}" ]] && return 0

    mtd_hpc_global_transfer_slot_acquire \
        stagein \
        "$MTD_HPC_STAGEIN_MAX_CONCURRENT" \
        MTD_HPC_STAGEIN_GLOBAL_LOCK_HANDLE \
        MTD_HPC_STAGEIN_GLOBAL_LOCK_SLOT
}

mtd_hpc_stagein_group_end() {
    mtd_hpc_global_transfer_slot_release \
        stagein \
        MTD_HPC_STAGEIN_GLOBAL_LOCK_HANDLE \
        MTD_HPC_STAGEIN_GLOBAL_LOCK_SLOT
    MTD_HPC_STAGEIN_GROUP_ACTIVE=0
}

mtd_hpc_stage_in_file() {
    local source="$1"
    local destination="$2"
    local label="${3:-input}"
    local expected_size="${4:-}"
    local source_size=0
    local destination_size=0
    local remote_target=""
    local rsync_status=0
    local automatic_group=0

    command -v rsync >/dev/null 2>&1 || \
        mtd_hpc_die "rsync is required for node-local input staging." || return 1

    mkdir -p -- "$(dirname -- "$destination")"

    # Most node jobs explicitly group their complete input set. Keep a safe
    # fallback for any older/direct caller that stages one file at a time.
    if [[ "${MTD_HPC_STAGEIN_GROUP_ACTIVE:-0}" != "1" ]]; then
        mtd_hpc_stagein_group_begin
        automatic_group=1
    fi

    if [[ -s "$source" ]]; then
        source_size="$(stat -Lc %s -- "$source")" || return 1
        if [[ -n "$expected_size" && "$source_size" != "$expected_size" ]]; then
            mtd_hpc_die \
                "$label changed after task submission: expected=$expected_size current=$source_size bytes: $source" || return 1
        fi

        mtd_hpc_stagein_lock_acquire || return 1

        mtd_hpc_info "Stage-in $label: $source -> $destination (${source_size} bytes)"
        set +e
        rsync -aL --whole-file -- "$source" "$destination"
        rsync_status=$?
        set -e
    else
        [[ "$expected_size" =~ ^[0-9]+$ ]] && (( expected_size > 0 )) || \
            mtd_hpc_die \
                "$label is not visible on this compute node and no validated source size was supplied: $source" || return 1

        remote_target="$(mtd_hpc_resolve_submit_ssh_target)" || return 1
        source_size="$expected_size"

        mtd_hpc_stagein_lock_acquire || return 1

        mtd_hpc_info \
            "$label is not mounted on node ${SLURMD_NODENAME:-$(hostname -s)};" \
            "pulling it from the pipeline submission host $remote_target."
        mtd_hpc_info \
            "Stage-in $label: ${remote_target}:${source} -> $destination (${source_size} bytes)"

        set +e
        rsync \
            -aL \
            --whole-file \
            --protect-args \
            -e "ssh -o BatchMode=yes -o ConnectTimeout=$MTD_HPC_SUBMIT_SSH_CONNECT_TIMEOUT" \
            -- "${remote_target}:${source}" "$destination"
        rsync_status=$?
        set -e
    fi

    if (( automatic_group == 1 )); then
        mtd_hpc_stagein_group_end
    fi

    (( rsync_status == 0 )) || return "$rsync_status"

    mtd_hpc_require_file "$destination" "staged $label" || return 1

    destination_size="$(stat -Lc %s -- "$destination")"
    [[ "$destination_size" == "$source_size" ]] || \
        mtd_hpc_die "Stage-in size mismatch for $label: source=$source_size destination=$destination_size" || return 1
}

mtd_hpc_stageout_lock_acquire() {
    MTD_HPC_STAGEOUT_LOCK_FD=""
    MTD_HPC_GLOBAL_STAGEOUT_LOCK_HANDLE=""
    MTD_HPC_GLOBAL_STAGEOUT_LOCK_SLOT=""

    mtd_hpc_global_transfer_slot_acquire \
        stageout \
        "$MTD_HPC_STAGEOUT_MAX_CONCURRENT" \
        MTD_HPC_GLOBAL_STAGEOUT_LOCK_HANDLE \
        MTD_HPC_GLOBAL_STAGEOUT_LOCK_SLOT || return 1

    [[ "$MTD_HPC_SERIALIZE_STAGEOUT_PER_NODE" == "1" ]] || return 0
    command -v flock >/dev/null 2>&1 || {
        mtd_hpc_global_transfer_slot_release \
            stageout \
            MTD_HPC_GLOBAL_STAGEOUT_LOCK_HANDLE \
            MTD_HPC_GLOBAL_STAGEOUT_LOCK_SLOT
        mtd_hpc_die "flock is required when MTD_HPC_SERIALIZE_STAGEOUT_PER_NODE=1."
        return 1
    }

    mkdir -p -- "$MTD_HPC_LOCAL_SCRATCH_ROOT/.locks"
    exec {MTD_HPC_STAGEOUT_LOCK_FD}> \
        "$MTD_HPC_LOCAL_SCRATCH_ROOT/.locks/stageout.lock"
    flock -x "$MTD_HPC_STAGEOUT_LOCK_FD"

    mtd_hpc_info "Acquired node-local stage-out lock."
}

mtd_hpc_stageout_lock_release() {
    if [[ -n "${MTD_HPC_STAGEOUT_LOCK_FD:-}" ]]; then
        flock -u "$MTD_HPC_STAGEOUT_LOCK_FD" || true
        eval "exec ${MTD_HPC_STAGEOUT_LOCK_FD}>&-"
        MTD_HPC_STAGEOUT_LOCK_FD=""
        mtd_hpc_info "Released node-local stage-out lock."
    fi

    mtd_hpc_global_transfer_slot_release \
        stageout \
        MTD_HPC_GLOBAL_STAGEOUT_LOCK_HANDLE \
        MTD_HPC_GLOBAL_STAGEOUT_LOCK_SLOT
}

mtd_hpc_atomic_stage_out_file() {
    local source="$1"
    local destination="$2"
    local label="${3:-output}"
    local allow_empty="${4:-0}"
    local destination_dir=""
    local temporary_destination=""
    local source_size=0
    local copied_size=0
    local rc=0

    if [[ "$allow_empty" == "1" ]]; then
        mtd_hpc_require_path_exists "$source" "$label" || return 1
    else
        mtd_hpc_require_file "$source" "$label" || return 1
    fi
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
        if [[ "$allow_empty" == "1" && ! -e "$temporary_destination" ]]; then
            mtd_hpc_error "Staged $label is missing: $temporary_destination"
            rc=1
        elif [[ "$allow_empty" != "1" && ! -s "$temporary_destination" ]]; then
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
    if [[ "$allow_empty" == "1" ]]; then
        mtd_hpc_require_path_exists "$destination" "staged $label" || return 1
    else
        mtd_hpc_require_file "$destination" "staged $label" || return 1
    fi
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
    local expected=""
    local path=""

    [[ -n "$outputs" ]] || return 0
    while IFS= read -r expected; do
        [[ -z "$expected" ]] && continue
        case "$expected" in
            exists:*)
                path="${expected#exists:}"
                [[ -e "$path" ]] || return 1
                ;;
            nonempty:*)
                path="${expected#nonempty:}"
                [[ -s "$path" ]] || return 1
                ;;
            *)
                [[ -s "$expected" ]] || return 1
                ;;
        esac
    done <<< "$outputs"
    return 0
}

mtd_hpc_output_spec_is_valid() {
    local expected="$1"
    local path=""

    case "$expected" in
        exists:*)
            path="${expected#exists:}"
            [[ -e "$path" ]]
            ;;
        nonempty:*)
            path="${expected#nonempty:}"
            [[ -s "$path" ]]
            ;;
        *)
            [[ -s "$expected" ]]
            ;;
    esac
}

mtd_hpc_output_spec_path() {
    local expected="$1"
    case "$expected" in
        exists:*) printf '%s\n' "${expected#exists:}" ;;
        nonempty:*) printf '%s\n' "${expected#nonempty:}" ;;
        *) printf '%s\n' "$expected" ;;
    esac
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
