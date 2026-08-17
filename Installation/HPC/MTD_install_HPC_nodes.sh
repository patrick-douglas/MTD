#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
REPO_ROOT_DEFAULT="$(readlink -f "$SCRIPT_DIR/../..")"

usage() {
    cat <<'USAGE'
Configure one or more Slurm compute nodes for MTD Explorer HPC execution.

Usage:
  bash Installation/HPC/MTD_install_HPC_nodes.sh \
    --node HOSTNAME [--user USER]

  bash Installation/HPC/MTD_install_HPC_nodes.sh \
    --node-list nodes.txt [--user USER]

Options:
  --node HOSTNAME              Configure one node.
  --node-list FILE             Configure each hostname listed in FILE.
  --user USER                  Non-root owner. Auto-detected unless launched
                               from a root login shell; root is forbidden.
  --prefix DIR                 Node-local installation prefix.
                               Default: /MTD_explorer_HPC
  --repo-root DIR              Main MTD Explorer repository.
  --ssh-user USER              SSH login used for file transfer and ordinary
                               commands. Default: --user. In this release it
                               must equal --user.
  --remote-root-mode MODE      sudo or direct. Default: sudo.
                               sudo: SSH as --user and run sudo -n remotely.
                               direct: SSH to root@HOST for root operations.
  --database SOURCE=RELATIVE   Copy an additional database directory with
                               rsync into PREFIX/databases/RELATIVE. Repeatable.
                               Complete repository-root kraken2DB_* databases
                               are discovered automatically.
  --database-distribution MODE
                               direct or propagate. Default: direct.
                               propagate uses complete nodes from --node-list
                               as one-at-a-time sources for pending nodes.
  --skip-default-databases     Do not auto-copy HUMAnN/ref_database, complete
                               repository-root kraken2DB_* databases, or host
                               Magic-BLAST database directories.
  --force-recreate-env         Remove and recreate all node-local MTD envs.
  --dry-run                    Validate nodes and print mutating commands.
  --help                       Show this help.

Complete HUMAnN/MetaPhlAn, Magic-BLAST, fastp, Kraken2 and Bracken nodes
currently require x86_64 Linux.
USAGE
}

REQUESTED_USER=""
NODE=""
NODE_LIST=""
PREFIX="/MTD_explorer_HPC"
REPO_ROOT="$REPO_ROOT_DEFAULT"
SSH_USER=""
REMOTE_ROOT_MODE="sudo"
DATABASE_DISTRIBUTION="direct"
SKIP_DEFAULT_DATABASES=0
FORCE_RECREATE_ENV=0
DRY_RUN=0
DATABASE_SPECS=()
ORIGINAL_ARGS=("$@")
SSH_OPTIONS=(-o BatchMode=yes -o ConnectTimeout=20 -o ServerAliveInterval=15)

fatal() {
    printf '[ERROR] %s\n' "$*" >&2
    exit 1
}

info() {
    printf '[INFO] %s\n' "$*"
}

ok() {
    printf '[OK] %s\n' "$*"
}

need_value() {
    [[ $# -ge 2 && -n "${2:-}" ]] || fatal "$1 requires a value."
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --node) need_value "$1" "${2:-}"; NODE="$2"; shift 2 ;;
        --node=*) NODE="${1#*=}"; shift ;;
        --node-list) need_value "$1" "${2:-}"; NODE_LIST="$2"; shift 2 ;;
        --node-list=*) NODE_LIST="${1#*=}"; shift ;;
        --user) need_value "$1" "${2:-}"; REQUESTED_USER="$2"; shift 2 ;;
        --user=*) REQUESTED_USER="${1#*=}"; shift ;;
        --prefix) need_value "$1" "${2:-}"; PREFIX="$2"; shift 2 ;;
        --prefix=*) PREFIX="${1#*=}"; shift ;;
        --repo-root) need_value "$1" "${2:-}"; REPO_ROOT="$2"; shift 2 ;;
        --repo-root=*) REPO_ROOT="${1#*=}"; shift ;;
        --ssh-user) need_value "$1" "${2:-}"; SSH_USER="$2"; shift 2 ;;
        --ssh-user=*) SSH_USER="${1#*=}"; shift ;;
        --remote-root-mode) need_value "$1" "${2:-}"; REMOTE_ROOT_MODE="$2"; shift 2 ;;
        --remote-root-mode=*) REMOTE_ROOT_MODE="${1#*=}"; shift ;;
        --database) need_value "$1" "${2:-}"; DATABASE_SPECS+=("$2"); shift 2 ;;
        --database=*) DATABASE_SPECS+=("${1#*=}"); shift ;;
        --database-distribution) need_value "$1" "${2:-}"; DATABASE_DISTRIBUTION="$2"; shift 2 ;;
        --database-distribution=*) DATABASE_DISTRIBUTION="${1#*=}"; shift ;;
        --skip-default-databases) SKIP_DEFAULT_DATABASES=1; shift ;;
        --force-recreate-env) FORCE_RECREATE_ENV=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        --help) usage; exit 0 ;;
        *) fatal "Unknown option: $1" ;;
    esac
done

# Elevate locally, preserving the login owner before sudo changes identity.
if (( EUID != 0 )); then
    detected_user="$(id -un)"
    info "Installation owner detected: $detected_user"
    info "Root privileges are required."
    exec sudo -- /usr/bin/env \
        "MTD_HPC_ORIGINAL_USER=$detected_user" \
        bash "$SCRIPT_PATH" "${ORIGINAL_ARGS[@]}"
fi

if [[ -n "$REQUESTED_USER" ]]; then
    OWNER_USER="$REQUESTED_USER"
elif [[ -n "${MTD_HPC_ORIGINAL_USER:-}" && "$MTD_HPC_ORIGINAL_USER" != "root" ]]; then
    OWNER_USER="$MTD_HPC_ORIGINAL_USER"
elif [[ -n "${SUDO_USER:-}" && "$SUDO_USER" != "root" ]] &&
     [[ "${SUDO_COMMAND:-}" == *"$(basename "$SCRIPT_PATH")"* ]]; then
    # Accept direct sudo invocation, but not an already-open `sudo -i` shell.
    OWNER_USER="$SUDO_USER"
else
    cat >&2 <<'ROOT_ERROR'
[ERROR] The installer was started from a root login shell and the non-root
installation owner cannot be determined safely.

Run it again with an explicit non-root owner, for example:

  bash Installation/HPC/MTD_install_HPC_nodes.sh \
      --node node01 \
      --user me

Do not use --user root.
ROOT_ERROR
    exit 1
fi

[[ "$OWNER_USER" != "root" ]] || fatal "--user must not be root."
getent passwd "$OWNER_USER" >/dev/null || fatal "Local user does not exist: $OWNER_USER"
OWNER_UID="$(id -u "$OWNER_USER")"
OWNER_GID="$(id -g "$OWNER_USER")"
OWNER_GROUP="$(id -gn "$OWNER_USER")"
OWNER_HOME="$(getent passwd "$OWNER_USER" | cut -d: -f6)"
(( OWNER_UID != 0 )) || fatal "The selected account has UID 0."

for command in ssh scp rsync curl sha256sum getent runuser install; do
    command -v "$command" >/dev/null 2>&1 || fatal "Required command not found on the main machine: $command"
done

[[ "$PREFIX" == /* && "$PREFIX" != "/" ]] || fatal "--prefix must be a non-root absolute path."
[[ "$PREFIX" != *$'\n'* && "$PREFIX" != *$'\t'* ]] || fatal "--prefix contains unsupported whitespace."
REPO_ROOT="$(readlink -f "$REPO_ROOT")"
[[ -s "$REPO_ROOT/MTD_explorer.sh" ]] || fatal "MTD Explorer repository not found: $REPO_ROOT"
ENV_YML="$REPO_ROOT/Installation/HPC/MTD-Explorer-HPC.yml"
FASTP_ENV_YML="$REPO_ROOT/Installation/MTD_fastp.yml"
KRAKEN2_ENV_YML="$REPO_ROOT/Installation/MTD_kraken2.yml"
[[ -s "$ENV_YML" ]] || fatal "HPC environment YAML not found: $ENV_YML"
[[ -s "$FASTP_ENV_YML" ]] || fatal "fastp environment YAML not found: $FASTP_ENV_YML"
[[ -s "$KRAKEN2_ENV_YML" ]] || fatal "Kraken2 environment YAML not found: $KRAKEN2_ENV_YML"

case "$REMOTE_ROOT_MODE" in
    sudo|direct) ;;
    *) fatal "--remote-root-mode must be sudo or direct." ;;
esac
case "$DATABASE_DISTRIBUTION" in
    direct|propagate) ;;
    *) fatal "--database-distribution must be direct or propagate." ;;
esac
if [[ "$DATABASE_DISTRIBUTION" == "propagate" && -z "$NODE_LIST" ]]; then
    fatal "--database-distribution propagate requires --node-list."
fi
[[ -n "$SSH_USER" ]] || SSH_USER="$OWNER_USER"
[[ "$SSH_USER" == "$OWNER_USER" ]] || \
    fatal "This release requires --ssh-user to equal --user so copied files retain operational ownership."

nodes=()
if [[ -n "$NODE" ]]; then
    [[ -z "$NODE_LIST" ]] || fatal "Use either --node or --node-list, not both."
    nodes+=("$NODE")
elif [[ -n "$NODE_LIST" ]]; then
    [[ -s "$NODE_LIST" ]] || fatal "Node list missing or empty: $NODE_LIST"
    while IFS= read -r host; do
        host="${host%%#*}"
        host="${host//[[:space:]]/}"
        [[ -n "$host" ]] || continue
        [[ "$host" =~ ^[A-Za-z0-9._-]+$ ]] || fatal "Invalid hostname in node list: $host"
        nodes+=("$host")
    done < "$NODE_LIST"
else
    fatal "One of --node or --node-list is required."
fi
(( ${#nodes[@]} > 0 )) || fatal "No hostnames were found."

# Keep the first occurrence of every hostname.
declare -A seen_nodes=()
unique_nodes=()
for node in "${nodes[@]}"; do
    [[ "$node" =~ ^[A-Za-z0-9._-]+$ ]] || fatal "Invalid hostname: $node"
    [[ -n "${seen_nodes[$node]:-}" ]] && continue
    seen_nodes[$node]=1
    unique_nodes+=("$node")
done
nodes=("${unique_nodes[@]}")

run_as_owner() {
    runuser -u "$OWNER_USER" -- env \
        HOME="$OWNER_HOME" USER="$OWNER_USER" LOGNAME="$OWNER_USER" "$@"
}

remote_query() {
    local node="$1"
    shift
    run_as_owner ssh "${SSH_OPTIONS[@]}" "$SSH_USER@$node" "$@"
}

remote_root() {
    local node="$1"
    shift
    local target
    if [[ "$REMOTE_ROOT_MODE" == "direct" ]]; then
        target="root@$node"
    else
        target="$SSH_USER@$node"
    fi

    if (( DRY_RUN )); then
        printf '[DRY-RUN] ssh %q' "$target"
        [[ "$REMOTE_ROOT_MODE" == "sudo" ]] && printf ' sudo -n'
        printf ' %q' "$@"
        printf '\n'
        return 0
    fi

    if [[ "$REMOTE_ROOT_MODE" == "direct" ]]; then
        run_as_owner ssh "${SSH_OPTIONS[@]}" "$target" "$@"
    else
        run_as_owner ssh "${SSH_OPTIONS[@]}" "$target" sudo -n "$@"
    fi
}

remote_owner() {
    local node="$1"
    shift
    if (( DRY_RUN )); then
        printf '[DRY-RUN] ssh %q' "$SSH_USER@$node"
        printf ' %q' "$@"
        printf '\n'
        return 0
    fi
    remote_query "$node" "$@"
}

copy_as_owner() {
    if (( DRY_RUN )); then
        printf '[DRY-RUN]'
        printf ' %q' "$@"
        printf '\n'
        return 0
    fi
    run_as_owner "$@"
}

# OpenSSH joins remote command arguments with spaces instead of preserving
# their original argument boundaries. Shell-quote every argument when a
# command (notably rsync -e) must retain an argument containing spaces.
remote_owner_quoted() {
    local node="$1"
    shift
    local remote_command
    printf -v remote_command '%q ' "$@"

    if (( DRY_RUN )); then
        printf '[DRY-RUN] ssh %q %s\n' "$SSH_USER@$node" "$remote_command"
        return 0
    fi
    remote_query "$node" "$remote_command"
}

is_complete_kraken2_database() {
    local database_dir="$1"
    [[ -d "$database_dir" ]] &&
        [[ -s "$database_dir/hash.k2d" ]] &&
        [[ -s "$database_dir/opts.k2d" ]] &&
        [[ -s "$database_dir/taxo.k2d" ]]
}

# Build database copy list. Relative targets mirror the main repository.
# Every complete repository-root kraken2DB_* directory is synchronized
# automatically, including custom host databases created after installation.
AUTO_KRAKEN_DB_TARGETS=()
if (( SKIP_DEFAULT_DATABASES == 0 )); then
    [[ -d "$REPO_ROOT/HUMAnN/ref_database" ]] && \
        DATABASE_SPECS+=("$REPO_ROOT/HUMAnN/ref_database=MTD-Explorer/HUMAnN/ref_database")

    shopt -s nullglob
    kraken_candidates=( "$REPO_ROOT"/kraken2DB_* )
    shopt -u nullglob

    for candidate in "${kraken_candidates[@]}"; do
        [[ -d "$candidate" ]] || continue

        if ! is_complete_kraken2_database "$candidate"; then
            fatal "Repository Kraken2 database is incomplete: $candidate (required: hash.k2d, opts.k2d, taxo.k2d)"
        fi

        relative_target="MTD-Explorer/$(basename "$candidate")"
        DATABASE_SPECS+=("$candidate=$relative_target")
        AUTO_KRAKEN_DB_TARGETS+=("$relative_target")
    done

    for candidate in \
        "$REPO_ROOT"/*_blastdb \
        "$REPO_ROOT"/blastdb_*
    do
        [[ -d "$candidate" ]] || continue
        DATABASE_SPECS+=("$candidate=MTD-Explorer/$(basename "$candidate")")
    done
fi

# Validate and deduplicate database mappings before contacting nodes.
declare -A seen_targets=()
validated_database_specs=()
REQUIRED_KRAKEN_DB_TARGETS=()
for spec in "${DATABASE_SPECS[@]}"; do
    [[ "$spec" == *=* ]] || fatal "Invalid --database specification: $spec"
    source_path="${spec%%=*}"
    relative_target="${spec#*=}"
    source_path="$(readlink -f "$source_path")"
    [[ -d "$source_path" ]] || fatal "Database source must be an existing directory: $source_path"
    [[ -n "$relative_target" && "$relative_target" != /* ]] || \
        fatal "Database target must be a non-empty relative path: $relative_target"
    [[ "/$relative_target/" != *"/../"* && "/$relative_target/" != *"/./"* ]] || \
        fatal "Database target is unsafe: $relative_target"
    if [[ -n "${seen_targets[$relative_target]:-}" ]]; then
        info "Skipping duplicate database target: $relative_target"
        continue
    fi

    case "$(basename "$relative_target")" in
        kraken2DB_*)
            if ! is_complete_kraken2_database "$source_path"; then
                fatal "Kraken2 database mapping is incomplete: $source_path -> $relative_target (required: hash.k2d, opts.k2d, taxo.k2d)"
            fi
            REQUIRED_KRAKEN_DB_TARGETS+=("$relative_target")
            ;;
    esac

    seen_targets[$relative_target]=1
    validated_database_specs+=("$source_path=$relative_target")
done
DATABASE_SPECS=("${validated_database_specs[@]}")

info "Installation owner: $OWNER_USER ($OWNER_UID:$OWNER_GID; group $OWNER_GROUP)"
info "Repository: $REPO_ROOT"
info "Node-local prefix: $PREFIX"
info "Nodes: ${nodes[*]}"
info "Database directories: ${#DATABASE_SPECS[@]}"
info "Auto-detected repository Kraken2 databases: ${#AUTO_KRAKEN_DB_TARGETS[@]}"
for relative_target in "${AUTO_KRAKEN_DB_TARGETS[@]}"; do
    info "  Kraken2: $relative_target"
done
info "Required Kraken2 databases for final validation: ${#REQUIRED_KRAKEN_DB_TARGETS[@]}"
info "Database distribution: $DATABASE_DISTRIBUTION"

CACHE_DIR="$REPO_ROOT/Installation/HPC/cache"
mkdir -p -- "$CACHE_DIR"
chown "$OWNER_USER:$OWNER_GROUP" "$CACHE_DIR"

CONDARC_LOCAL="$CACHE_DIR/condarc"
cat > "$CONDARC_LOCAL" <<'CONDARC'
channels:
  - conda-forge
  - bioconda
channel_priority: strict
show_channel_urls: true
CONDARC
chown "$OWNER_USER:$OWNER_GROUP" "$CONDARC_LOCAL"

for node in "${nodes[@]}"; do
    printf '\n============================================================\n'
    printf 'Configuring node: %s\n' "$node"
    printf '============================================================\n'

    remote_query "$node" true || fatal "Cannot connect to $SSH_USER@$node with non-interactive SSH."
    if [[ "$REMOTE_ROOT_MODE" == "sudo" ]]; then
        remote_query "$node" sudo -n true || \
            fatal "Remote passwordless sudo is not available for $SSH_USER on $node. Use --remote-root-mode direct if appropriate."
    else
        run_as_owner ssh "${SSH_OPTIONS[@]}" "root@$node" true || \
            fatal "Cannot connect to root@$node for --remote-root-mode direct."
    fi

    # Node-local stage-in/stage-out requires rsync, while packed chunk
    # stage-out serialization uses flock from util-linux.
    missing_remote_commands=()
    for remote_command in rsync flock; do
        if ! remote_query "$node" command -v "$remote_command" >/dev/null 2>&1; then
            missing_remote_commands+=("$remote_command")
        fi
    done

    if (( ${#missing_remote_commands[@]} > 0 )); then
        info "Installing missing node utilities on $node: ${missing_remote_commands[*]}"

        package_manager="$(remote_query "$node" bash -c '
            if command -v apt-get >/dev/null 2>&1; then echo apt-get;
            elif command -v dnf >/dev/null 2>&1; then echo dnf;
            elif command -v yum >/dev/null 2>&1; then echo yum;
            else echo none; fi
        ' | tr -d '\r\n')"

        case "$package_manager" in
            apt-get)
                remote_root "$node" env DEBIAN_FRONTEND=noninteractive apt-get update
                remote_root "$node" env DEBIAN_FRONTEND=noninteractive \
                    apt-get install -y rsync ca-certificates util-linux
                ;;
            dnf)
                remote_root "$node" dnf install -y rsync ca-certificates util-linux
                ;;
            yum)
                remote_root "$node" yum install -y rsync ca-certificates util-linux
                ;;
            *)
                fatal "Missing node utilities (${missing_remote_commands[*]}) and no supported package manager was found on $node."
                ;;
        esac
    fi

    arch="$(remote_query "$node" uname -m | tr -d '\r\n')"
    case "$arch" in
        x86_64|amd64)
            installer_name="Miniconda3-latest-Linux-x86_64.sh"
            ;;
        aarch64|arm64)
            fatal "Node $node is $arch. The complete runtime currently includes Magic-BLAST 1.7.0, whose Bioconda package is not available for linux-aarch64. Use x86_64 nodes for this first implementation."
            ;;
        *) fatal "Unsupported architecture on $node: $arch" ;;
    esac

    remote_uid="$(remote_query "$node" id -u "$OWNER_USER" | tr -d '\r\n')" || \
        fatal "User $OWNER_USER does not exist on $node."
    remote_gid="$(remote_query "$node" id -g "$OWNER_USER" | tr -d '\r\n')" || \
        fatal "Could not read GID for $OWNER_USER on $node."
    [[ "$remote_uid:$remote_gid" == "$OWNER_UID:$OWNER_GID" ]] || \
        fatal "UID/GID mismatch on $node. Main=$OWNER_UID:$OWNER_GID remote=$remote_uid:$remote_gid"

    installer_url="https://repo.anaconda.com/miniconda/$installer_name"
    installer_local="$CACHE_DIR/$installer_name"
    if [[ ! -s "$installer_local" ]]; then
        info "Downloading $installer_name on the main machine."
        if (( DRY_RUN )); then
            printf '[DRY-RUN] curl -fL %q -o %q\n' "$installer_url" "$installer_local"
        else
            run_as_owner curl -fL --retry 3 --retry-delay 5 \
                "$installer_url" -o "$installer_local"
        fi
    fi

    remote_root "$node" install -d -o "$OWNER_UID" -g "$OWNER_GID" -m 0755 \
        "$PREFIX" "$PREFIX/envs" "$PREFIX/databases" "$PREFIX/config" \
        "$PREFIX/logs" "$PREFIX/cache"

    # Scratch is private to the MTD execution owner. Tasks create unique
    # per-job directories below this node-local path.
    remote_root "$node" install -d -o "$OWNER_UID" -g "$OWNER_GID" -m 0700 \
        "$PREFIX/tmp"

    copy_as_owner scp -q "${SSH_OPTIONS[@]}" \
        "$installer_local" "$SSH_USER@$node:$PREFIX/cache/$installer_name"
    copy_as_owner scp -q "${SSH_OPTIONS[@]}" \
        "$ENV_YML" "$SSH_USER@$node:$PREFIX/config/MTD-Explorer-HPC.yml"
    copy_as_owner scp -q "${SSH_OPTIONS[@]}" \
        "$FASTP_ENV_YML" "$SSH_USER@$node:$PREFIX/config/MTD_fastp.yml"
    copy_as_owner scp -q "${SSH_OPTIONS[@]}" \
        "$KRAKEN2_ENV_YML" "$SSH_USER@$node:$PREFIX/config/MTD_kraken2.yml"
    copy_as_owner scp -q "${SSH_OPTIONS[@]}" \
        "$CONDARC_LOCAL" "$SSH_USER@$node:$PREFIX/config/condarc"

    if (( FORCE_RECREATE_ENV )); then
        remote_root "$node" rm -rf -- \
            "$PREFIX/envs/MTD-Explorer-HPC" \
            "$PREFIX/envs/MTD_fastp" \
            "$PREFIX/envs/MTD_kraken2"
    fi

    if ! remote_query "$node" test -x "$PREFIX/miniconda3/bin/conda"; then
        remote_owner "$node" bash "$PREFIX/cache/$installer_name" -b -p "$PREFIX/miniconda3"
    else
        info "Using existing Miniconda on $node."
    fi

    install_or_update_environment() {
        local environment_name="$1"
        local yaml_name="$2"
        local environment_prefix="$PREFIX/envs/$environment_name"
        local yaml_path="$PREFIX/config/$yaml_name"

        if ! remote_query "$node" test -d "$environment_prefix/conda-meta"; then
            info "Creating $environment_name environment on $node."
            remote_owner "$node" env \
                PYTHONNOUSERSITE=1 \
                CONDARC="$PREFIX/config/condarc" \
                "$PREFIX/miniconda3/bin/conda" env create \
                --prefix "$environment_prefix" \
                --file "$yaml_path"
        else
            info "Updating existing $environment_name environment on $node."
            remote_owner "$node" env \
                PYTHONNOUSERSITE=1 \
                CONDARC="$PREFIX/config/condarc" \
                "$PREFIX/miniconda3/bin/conda" env update \
                --prefix "$environment_prefix" \
                --file "$yaml_path" \
                --prune
        fi
    }

    install_or_update_environment "MTD-Explorer-HPC" "MTD-Explorer-HPC.yml"
    install_or_update_environment "MTD_fastp" "MTD_fastp.yml"
    install_or_update_environment "MTD_kraken2" "MTD_kraken2.yml"

    remote_root "$node" chown -R "$OWNER_UID:$OWNER_GID" "$PREFIX"
done

database_target_is_complete() {
    local node="$1"
    local spec source_path relative_target changes

    (( DRY_RUN == 0 )) || return 1
    for spec in "${DATABASE_SPECS[@]}"; do
        source_path="${spec%%=*}"
        relative_target="${spec#*=}"
        remote_query "$node" test -d "$PREFIX/databases/$relative_target" || return 1
        changes="$(run_as_owner rsync -aHn --delete --itemize-changes \
            -e "ssh ${SSH_OPTIONS[*]}" \
            "$source_path/" "$SSH_USER@$node:$PREFIX/databases/$relative_target/")" || return 1
        [[ -z "$changes" ]] || return 1
    done
}

copy_databases_from_master() {
    local destination="$1"
    local spec source_path relative_target

    for spec in "${DATABASE_SPECS[@]}"; do
        source_path="${spec%%=*}"
        relative_target="${spec#*=}"
        remote_root "$destination" install -d -o "$OWNER_UID" -g "$OWNER_GID" -m 0755 \
            "$PREFIX/databases/$relative_target"
        info "Synchronizing database master -> $destination: $relative_target"
        copy_as_owner rsync -aH --delete --info=progress2 \
            -e "ssh ${SSH_OPTIONS[*]}" \
            "$source_path/" "$SSH_USER@$destination:$PREFIX/databases/$relative_target/"
    done
}

copy_databases_from_node() {
    local source="$1"
    local destination="$2"
    local spec relative_target

    remote_query "$source" ssh "${SSH_OPTIONS[@]}" "$SSH_USER@$destination" true || \
        fatal "Source $source cannot connect non-interactively to $SSH_USER@$destination."
    for spec in "${DATABASE_SPECS[@]}"; do
        relative_target="${spec#*=}"
        remote_root "$destination" install -d -o "$OWNER_UID" -g "$OWNER_GID" -m 0755 \
            "$PREFIX/databases/$relative_target"
        info "Synchronizing database $source -> $destination: $relative_target"
        remote_owner_quoted "$source" rsync -aH --delete --info=progress2 \
            -e "ssh ${SSH_OPTIONS[*]}" \
            "$PREFIX/databases/$relative_target/" \
            "$SSH_USER@$destination:$PREFIX/databases/$relative_target/"
    done
}

if [[ "$DATABASE_DISTRIBUTION" == "direct" ]]; then
    for node in "${nodes[@]}"; do
        copy_databases_from_master "$node"
    done
else
    ready_sources=(master)
    pending_nodes=()
    for node in "${nodes[@]}"; do
        if database_target_is_complete "$node"; then
            info "Using node with complete database targets as a source: $node"
            ready_sources+=("$node")
        else
            pending_nodes+=("$node")
        fi
    done

    while (( ${#pending_nodes[@]} > 0 )); do
        round_destinations=()
        round_sources=()
        round_pids=()
        assignments=${#ready_sources[@]}
        (( assignments > ${#pending_nodes[@]} )) && assignments=${#pending_nodes[@]}

        for (( index=0; index<assignments; index++ )); do
            source="${ready_sources[$index]}"
            destination="${pending_nodes[$index]}"
            round_sources+=("$source")
            round_destinations+=("$destination")
            if [[ "$source" == "master" ]]; then
                copy_databases_from_master "$destination" &
            else
                copy_databases_from_node "$source" "$destination" &
            fi
            round_pids+=("$!")
        done

        round_failed=0
        for (( index=0; index<assignments; index++ )); do
            if ! wait "${round_pids[$index]}"; then
                printf '[ERROR] Database transfer failed: %s -> %s\n' \
                    "${round_sources[$index]}" "${round_destinations[$index]}" >&2
                round_failed=1
            fi
        done
        (( round_failed == 0 )) || fatal "Database propagation stopped after a failed transfer."

        for destination in "${round_destinations[@]}"; do
            if (( DRY_RUN == 0 )); then
                database_target_is_complete "$destination" || \
                    fatal "Database validation failed after transfer to $destination."
            fi
            ready_sources+=("$destination")
            ok "Database targets complete; promoted to source: $destination"
        done
        pending_nodes=("${pending_nodes[@]:assignments}")
    done
fi

for node in "${nodes[@]}"; do
    remote_root "$node" chown -R "$OWNER_UID:$OWNER_GID" "$PREFIX"

    if remote_query "$node" test -d "$PREFIX/databases/MTD-Explorer/HUMAnN/ref_database"; then
        for pair in \
            "nucleotide:$PREFIX/databases/MTD-Explorer/HUMAnN/ref_database/chocophlan" \
            "protein:$PREFIX/databases/MTD-Explorer/HUMAnN/ref_database/uniref" \
            "utility_mapping:$PREFIX/databases/MTD-Explorer/HUMAnN/ref_database/utility_mapping"
        do
            key="${pair%%:*}"
            value="${pair#*:}"
            remote_owner "$node" env PYTHONNOUSERSITE=1 \
                "$PREFIX/miniconda3/bin/conda" run --prefix "$PREFIX/envs/MTD-Explorer-HPC" \
                humann_config --update database_folders "$key" "$value"
        done
    fi

    for version_command in \
        "python --version" \
        "humann --version" \
        "metaphlan --version" \
        "magicblast -version"
    do
        # The command string is fixed by this script, not user input.
        remote_owner "$node" bash -c \
            "PYTHONNOUSERSITE=1 '$PREFIX/miniconda3/bin/conda' run --prefix '$PREFIX/envs/MTD-Explorer-HPC' $version_command"
    done

    remote_owner "$node" env PYTHONNOUSERSITE=1 \
        "$PREFIX/miniconda3/bin/conda" run --prefix "$PREFIX/envs/MTD_fastp" \
        fastp --version
    remote_owner "$node" env PYTHONNOUSERSITE=1 \
        "$PREFIX/miniconda3/bin/conda" run --prefix "$PREFIX/envs/MTD_kraken2" \
        kraken2 --version
    if ! remote_owner "$node" env PYTHONNOUSERSITE=1 \
        "$PREFIX/miniconda3/bin/conda" run --prefix "$PREFIX/envs/MTD_kraken2" \
        bracken -v
    then
        remote_owner "$node" env PYTHONNOUSERSITE=1 \
            "$PREFIX/miniconda3/bin/conda" run --prefix "$PREFIX/envs/MTD_kraken2" \
            bracken -h >/dev/null
    fi

    configured_at="$(date --iso-8601=seconds)"
    manifest_local="$(mktemp)"
    cat > "$manifest_local" <<MANIFEST
node=$node
architecture=$arch
owner=$OWNER_USER
uid=$OWNER_UID
gid=$OWNER_GID
prefix=$PREFIX
environment=$PREFIX/envs/MTD-Explorer-HPC
fastp_environment=$PREFIX/envs/MTD_fastp
kraken2_environment=$PREFIX/envs/MTD_kraken2
database_root=$PREFIX/databases
scratch_root=$PREFIX/tmp
configured_at=$configured_at
required_kraken_db_count=${#REQUIRED_KRAKEN_DB_TARGETS[@]}
MANIFEST
    for relative_target in "${REQUIRED_KRAKEN_DB_TARGETS[@]}"; do
        printf 'required_kraken_db=%s\n' "$relative_target" >> "$manifest_local"
    done

    chown "$OWNER_USER:$OWNER_GROUP" "$manifest_local"
    copy_as_owner scp -q "${SSH_OPTIONS[@]}" \
        "$manifest_local" "$SSH_USER@$node:$PREFIX/config/node_manifest.txt"
    rm -f -- "$manifest_local"

    remote_root "$node" chown -R "$OWNER_UID:$OWNER_GID" "$PREFIX"

    info "Running final node validation: $node"
    checker_args=(
        --node "$node"
        --user "$OWNER_USER"
        --ssh-user "$SSH_USER"
        --prefix "$PREFIX"
    )
    for relative_target in "${REQUIRED_KRAKEN_DB_TARGETS[@]}"; do
        checker_args+=( --require-kraken-db "$relative_target" )
    done

    run_as_owner bash "$SCRIPT_DIR/MTD_check_HPC_nodes.sh" "${checker_args[@]}"

    ok "Node configured: $node"
done

ok "MTD Explorer HPC node configuration completed."
