#!/usr/bin/env bash
set -Eeuo pipefail

REPO="$(pwd)"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --repo) REPO="$2"; shift 2 ;;
        --repo=*) REPO="${1#*=}"; shift ;;
        *) printf '[ERROR] Unknown option: %s\n' "$1" >&2; exit 2 ;;
    esac
done
REPO="$(readlink -f "$REPO")"

required=(
    "$REPO/Installation/HPC/MTD_install_HPC_nodes.sh"
    "$REPO/Installation/HPC/MTD_check_HPC_nodes.sh"
    "$REPO/Installation/HPC/MTD-Explorer-HPC.yml"
    "$REPO/aux_scripts/hpc/mtd_hpc_common.sh"
    "$REPO/aux_scripts/hpc/mtd_hpc_submit_array.sh"
    "$REPO/aux_scripts/hpc/mtd_hpc_array_task.sh"
    "$REPO/aux_scripts/hpc/mtd_hpc_node_job.sh"
    "$REPO/aux_scripts/hpc/mtd_hpc_humann_stage.sh"
    "$REPO/aux_scripts/hpc/mtd_hpc_magicblast_stage.sh"
    "$REPO/aux_scripts/hpc/mtd_split_fastq.py"
    "$REPO/aux_scripts/hpc/mtd_merge_sam_chunks.py"
)

for file in "${required[@]}"; do
    [[ -s "$file" ]] || { printf '[ERROR] Missing: %s\n' "$file" >&2; exit 1; }
done

for script in "$REPO"/Installation/HPC/*.sh "$REPO"/aux_scripts/hpc/*.sh; do
    bash -n "$script"
done
bash -n "$REPO/MTD_explorer.sh"

for script in "$REPO"/aux_scripts/hpc/*.py "$REPO"/update_fix/add_mtd_hpc_backend_20260730.py; do
    python3 -m py_compile "$script"
done

python3 "$REPO/update_fix/add_mtd_hpc_backend_20260730.py" --repo "$REPO" --check >/dev/null

grep -Fq 'MTD_HPC_BACKEND_20260730' "$REPO/MTD_explorer.sh" || {
    printf '[ERROR] MTD_explorer.sh does not contain the HPC integration marker.\n' >&2
    exit 1
}

grep -Fq -- '--hpc-conf' "$REPO/MTD_explorer.sh" || {
    printf '[ERROR] --hpc-conf was not added to MTD_explorer.sh.\n' >&2
    exit 1
}

# Regression checks for the integration boundaries.
grep -Fq 'show_sample_progress "Kraken2"' "$REPO/MTD_explorer.sh" || {
    printf '[ERROR] Kraken2 sample loops were unexpectedly removed by the HPC patch.\n' >&2
    exit 1
}
grep -Fq '# Run HUMAnN once per sample' "$REPO/MTD_explorer.sh" || {
    printf '[ERROR] HUMAnN execution section is missing after integration.\n' >&2
    exit 1
}
grep -Fq 'MTD_HPC_HUMANN_WORK_DIR="$HUMANN_WORK_DIR/hpc_humann"' "$REPO/MTD_explorer.sh" || {
    printf '[ERROR] HPC HUMAnN branch is missing after integration.\n' >&2
    exit 1
}

# FASTQ splitter and SAM merger smoke tests.
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
cat > "$tmp/R1.fq" <<'FASTQ1'
@a/1
ACGT
+
IIII
@b/1
TGCA
+
IIII
@c/1
AAAA
+
IIII
FASTQ1
cat > "$tmp/R2.fq" <<'FASTQ2'
@a/2
ACGT
+
IIII
@b/2
TGCA
+
IIII
@c/2
TTTT
+
IIII
FASTQ2
python3 "$REPO/aux_scripts/hpc/mtd_split_fastq.py" \
    --read1 "$tmp/R1.fq" --read2 "$tmp/R2.fq" \
    --output-dir "$tmp/chunks" --sample smoke \
    --records-per-chunk 2 --manifest "$tmp/chunks.tsv" >/dev/null
[[ "$(wc -l < "$tmp/chunks.tsv")" -eq 2 ]]
[[ "$(grep -c '^@' "$tmp/chunks/smoke.chunk000001.R1.fq")" -eq 2 ]]
[[ "$(grep -c '^@' "$tmp/chunks/smoke.chunk000002.R1.fq")" -eq 1 ]]

cat > "$tmp/a.sam" <<'SAM1'
@HD	VN:1.6
@SQ	SN:ref	LN:10
@PG	ID:magicblast	CL:magicblast -query chunk1.fq
r1	0	ref	1	60	4M	*	0	0	ACGT	IIII
SAM1
cat > "$tmp/b.sam" <<'SAM2'
@HD	VN:1.6
@SQ	SN:ref	LN:10
@PG	ID:magicblast	CL:magicblast -query chunk2.fq
r2	0	ref	2	60	4M	*	0	0	TGCA	IIII
SAM2
printf '%s\n%s\n' "$tmp/a.sam" "$tmp/b.sam" > "$tmp/sams.txt"
python3 "$REPO/aux_scripts/hpc/mtd_merge_sam_chunks.py" \
    --chunk-list "$tmp/sams.txt" --output "$tmp/merged.sam" >/dev/null
[[ "$(grep -c '^@HD' "$tmp/merged.sam")" -eq 1 ]]
[[ "$(grep -c '^@PG' "$tmp/merged.sam")" -eq 1 ]]
[[ "$(grep -vc '^@' "$tmp/merged.sam")" -eq 2 ]]

# A changed reference dictionary must be rejected.
sed 's/LN:10/LN:11/' "$tmp/b.sam" > "$tmp/b_bad.sam"
printf '%s\n%s\n' "$tmp/a.sam" "$tmp/b_bad.sam" > "$tmp/sams_bad.txt"
if python3 "$REPO/aux_scripts/hpc/mtd_merge_sam_chunks.py" \
    --chunk-list "$tmp/sams_bad.txt" --output "$tmp/should_not_exist.sam" >/dev/null 2>&1; then
    printf '[ERROR] SAM merger accepted incompatible reference headers.\n' >&2
    exit 1
fi

# Common resource/config smoke test.
cat > "$tmp/hpc.conf" <<'CONF'
MTD_HPC_MAX_PARALLEL=2
MTD_HPC_MAGICBLAST_CHUNK_READS=10
MTD_HPC_SBATCH_EXTRA_ARGS=(--exclusive --nodes=1 --mem=0)
CONF
# shellcheck source=/dev/null
source "$REPO/aux_scripts/hpc/mtd_hpc_common.sh"
mtd_hpc_load_config "$tmp/hpc.conf" "$REPO"
[[ "$MTD_HPC_MAX_PARALLEL" -eq 2 ]]
[[ "$(mtd_hpc_detect_node_threads)" =~ ^[0-9]+$ ]]

# End-to-end array submission smoke test with a local Slurm command shim.
mkdir -p "$tmp/mock_slurm"
cat > "$tmp/mock_slurm/sbatch" <<'SBATCH'
#!/usr/bin/env bash
set -Eeuo pipefail
array=""
script=""
script_index=0
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
    case "${args[$i]}" in
        --array=*) array="${args[$i]#*=}" ;;
        *.sh)
            script="${args[$i]}"
            script_index=$i
            break
            ;;
    esac
done
[[ -n "$array" && -n "$script" ]]
range="${array%%%*}"
start="${range%-*}"
end="${range#*-}"
printf 'call\n' >> "${MOCK_SBATCH_COUNT:?}"
for ((task_id=start; task_id<=end; task_id++)); do
    SLURM_JOB_ID=9001 SLURM_ARRAY_JOB_ID=9001 SLURM_ARRAY_TASK_ID="$task_id" \
        bash "$script" "${args[@]:$((script_index + 1))}" >/dev/null 2>&1
done
printf '9001\n'
SBATCH
cat > "$tmp/mock_slurm/squeue" <<'SQUEUE'
#!/usr/bin/env bash
exit 0
SQUEUE
cat > "$tmp/mock_slurm/sacct" <<'SACCT'
#!/usr/bin/env bash
printf '9001|COMPLETED|0:0\n'
SACCT
cat > "$tmp/mock_slurm/scancel" <<'SCANCEL'
#!/usr/bin/env bash
exit 0
SCANCEL
chmod +x "$tmp/mock_slurm"/*

cat > "$tmp/submit.conf" <<'CONF'
MTD_HPC_MAX_PARALLEL=2
MTD_HPC_POLL_SECONDS=1
MTD_HPC_MAGICBLAST_CHUNK_READS=10
MTD_HPC_RESUME=1
MTD_HPC_SBATCH_EXTRA_ARGS=()
CONF
: > "$tmp/tasks.tsv"
cmd1="printf '%s\n' alpha > '$tmp/task1.out'"
cmd2="printf '%s\n' beta > '$tmp/task2.out'"
mtd_hpc_task_line task1 "$cmd1" "$tmp/task1.out" >> "$tmp/tasks.tsv"
mtd_hpc_task_line task2 "$cmd2" "$tmp/task2.out" >> "$tmp/tasks.tsv"
: > "$tmp/sbatch.count"
PATH="$tmp/mock_slurm:$PATH" MOCK_SBATCH_COUNT="$tmp/sbatch.count" \
    bash "$REPO/aux_scripts/hpc/mtd_hpc_submit_array.sh" \
        --hpc-conf "$tmp/submit.conf" \
        --manifest "$tmp/tasks.tsv" \
        --stage smoke \
        --work-dir "$tmp/slurm_work" \
        --mtd-root "$REPO" >/dev/null
[[ "$(cat "$tmp/task1.out")" == "alpha" ]]
[[ "$(cat "$tmp/task2.out")" == "beta" ]]
[[ "$(wc -l < "$tmp/sbatch.count")" -eq 1 ]]

# A second identical call must reuse success markers without another submission.
PATH="$tmp/mock_slurm:$PATH" MOCK_SBATCH_COUNT="$tmp/sbatch.count" \
    bash "$REPO/aux_scripts/hpc/mtd_hpc_submit_array.sh" \
        --hpc-conf "$tmp/submit.conf" \
        --manifest "$tmp/tasks.tsv" \
        --stage smoke \
        --work-dir "$tmp/slurm_work" \
        --mtd-root "$REPO" >/dev/null
[[ "$(wc -l < "$tmp/sbatch.count")" -eq 1 ]]

printf '[OK] MTD Explorer HPC addon tests passed.\n'
