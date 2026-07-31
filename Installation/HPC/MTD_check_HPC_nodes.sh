#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'USAGE'
Usage:
  MTD_check_HPC_nodes.sh --node HOST [--user USER] [--prefix DIR]
  MTD_check_HPC_nodes.sh --node-list FILE [--user USER] [--prefix DIR]
USAGE
}

NODE=""
NODE_LIST=""
OWNER_USER="${SUDO_USER:-${USER:-}}"
PREFIX="/MTD_explorer_HPC"
SSH_USER=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --node) NODE="$2"; shift 2 ;;
        --node=*) NODE="${1#*=}"; shift ;;
        --node-list) NODE_LIST="$2"; shift 2 ;;
        --node-list=*) NODE_LIST="${1#*=}"; shift ;;
        --user) OWNER_USER="$2"; shift 2 ;;
        --user=*) OWNER_USER="${1#*=}"; shift ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        --prefix=*) PREFIX="${1#*=}"; shift ;;
        --ssh-user) SSH_USER="$2"; shift 2 ;;
        --ssh-user=*) SSH_USER="${1#*=}"; shift ;;
        --help) usage; exit 0 ;;
        *) printf '[ERROR] Unknown option: %s\n' "$1" >&2; usage; exit 2 ;;
    esac
done

[[ -n "$OWNER_USER" && "$OWNER_USER" != "root" ]] || {
    printf '[ERROR] --user must identify a non-root account.\n' >&2
    exit 2
}
[[ -n "$SSH_USER" ]] || SSH_USER="$OWNER_USER"

nodes=()
if [[ -n "$NODE" ]]; then
    [[ -z "$NODE_LIST" ]] || { printf '[ERROR] Use either --node or --node-list.\n' >&2; exit 2; }
    nodes+=("$NODE")
elif [[ -n "$NODE_LIST" ]]; then
    [[ -s "$NODE_LIST" ]] || { printf '[ERROR] Node list missing: %s\n' "$NODE_LIST" >&2; exit 2; }
    while IFS= read -r host; do
        host="${host%%#*}"
        host="${host//[[:space:]]/}"
        [[ -n "$host" ]] && nodes+=("$host")
    done < "$NODE_LIST"
fi
(( ${#nodes[@]} > 0 )) || { usage; exit 2; }

for node in "${nodes[@]}"; do
    printf '\n============================================================\n'
    printf 'Node: %s\n' "$node"
    printf '============================================================\n'
    ssh -o BatchMode=yes -o ConnectTimeout=20 "$SSH_USER@$node" bash -s -- "$PREFIX" "$OWNER_USER" <<'REMOTE'
set -Eeuo pipefail
prefix="$1"
owner="$2"

[[ "$(uname -m)" == "x86_64" || "$(uname -m)" == "amd64" ]]
test -d "$prefix"
test "$(stat -c %U "$prefix")" = "$owner"
test -x "$prefix/miniconda3/bin/conda"
test -d "$prefix/envs/MTD-Explorer-HPC/conda-meta"
test -d "$prefix/databases"
test -s "$prefix/config/node_manifest.txt"

humann_db="$prefix/databases/MTD-Explorer/HUMAnN/ref_database"
metaphlan_index="mpa_vJun23_CHOCOPhlAnSGB_202403"
test -d "$humann_db/chocophlan"
test -d "$humann_db/uniref"
test -d "$humann_db/utility_mapping"
test -d "$humann_db/metaphlan"
find "$humann_db/chocophlan" -type f -name '*.ffn.gz' -size +0c -print -quit | grep -q .
test -s "$humann_db/uniref/uniref90_201901b_full.dmnd"
utility_count="$(find "$humann_db/utility_mapping" -type f -size +0c | awk 'END {print NR + 0}')"
(( utility_count >= 10 ))
for suffix in pkl 1.bt2l 2.bt2l 3.bt2l 4.bt2l rev.1.bt2l rev.2.bt2l; do
    test -s "$humann_db/metaphlan/${metaphlan_index}.${suffix}"
done

conda_bin="$prefix/miniconda3/bin/conda"
env_dir="$prefix/envs/MTD-Explorer-HPC"

"$conda_bin" run --prefix "$env_dir" python --version
"$conda_bin" run --prefix "$env_dir" humann --version
"$conda_bin" run --prefix "$env_dir" metaphlan --version
"$conda_bin" run --prefix "$env_dir" magicblast -version

printf '[INFO] CPUs visible on node: %s\n' "$(nproc)"
awk '/^MemTotal:/ {printf "[INFO] Total memory: %s kB\\n", $2}' /proc/meminfo
printf '[OK] %s: MTD Explorer HPC runtime is ready.\n' "$(hostname -s)"
REMOTE
done
