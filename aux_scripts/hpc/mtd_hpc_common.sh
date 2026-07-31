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

    if ! declare -p MTD_HPC_SBATCH_EXTRA_ARGS >/dev/null 2>&1; then
        MTD_HPC_SBATCH_EXTRA_ARGS=(--exclusive --nodes=1 --mem=0)
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

    MTD_HPC_CONF="$conf"
    MTD_HPC_MTD_ROOT="$mtd_root"
    export MTD_HPC_CONF MTD_HPC_MTD_ROOT
    export MTD_HPC_PREFIX MTD_HPC_ENV_DIR MTD_HPC_CONDA_BIN
    export MTD_HPC_DATABASE_ROOT MTD_HPC_MTD_DATABASE_ROOT
    export MTD_HPC_HUMANN_DB_ROOT MTD_HPC_METAPHLAN_INDEX

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
