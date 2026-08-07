# MTD Explorer pipeline benchmark

This benchmark measures a complete `MTD_explorer.sh` analysis in either local
or Slurm HPC mode. The large scientific output remains separate from the compact
benchmark bundle.

`MTD_benchmark_install.sh` remains the low-level monitor for the controller
machine. The pipeline-specific reporters add MTD stage timings, output
inventory, diagnostics and, in HPC mode, Slurm accounting and per-node hardware
records.

## Files

- `run_mtd_pipeline_benchmark.sh`
  - validates paths and benchmark options;
  - selects `local` or `hpc` execution;
  - removes only the previous analysis output after confirmation;
  - runs the complete MTD Explorer analysis;
  - creates the final benchmark bundle.
- `MTD_benchmark_install.sh`
  - records end-to-end wall time and resource use on the controller machine.
- `MTD_pipeline_benchmark_report.py`
  - creates the common local/HPC pipeline summary;
  - calculates durations between MTD progress stages;
  - inventories output files;
  - creates a one-row comparison table;
  - produces focused diagnostics after a failed run.
- `MTD_slurm_benchmark_report.py`
  - discovers all MTD Slurm stage job histories under the pipeline output;
  - queries `sacct` for jobs, array tasks, queue time, elapsed time, CPU and RSS;
  - combines accounting with MTD task markers;
  - summarizes performance by stage and by compute node.
- `MTD_explorer.sh`
  - writes progress timestamps only when
    `MTD_PIPELINE_BENCHMARK_STEPS_TSV` is set;
  - normal runs are unchanged.

## Execution modes

The default is:

```text
--execution-mode auto
```

Resolution rules:

```text
auto without --hpc-conf  -> local
auto with --hpc-conf     -> hpc
local                     -> rejects --hpc-conf
hpc                       -> requires --hpc-conf
```

`--threads` always describes the controller thread setting passed to
`MTD_explorer.sh`. In local mode this is also the pipeline thread count. In HPC
mode the Slurm jobs determine their allocated/detected resources independently.

## Local benchmark

```bash
cd ~/MTD-Explorer

bash benchmark/run_mtd_pipeline_benchmark.sh \
  --execution-mode local \
  --machine master \
  --dataset Bglabrata_PRJNA1306560 \
  --input /path/to/B.glabrata_fastq/samplesheet.csv \
  --output "$HOME/test_MTD_explorer_B.glabrata.local" \
  --hostid 6526 \
  --threads 20 \
  --run-number 1
```

The old command remains valid because `auto` selects local mode when no HPC
configuration is supplied:

```bash
bash benchmark/run_mtd_pipeline_benchmark.sh \
  --machine master \
  --dataset Bglabrata_PRJNA1306560 \
  --input /path/to/samplesheet.csv \
  --output "$HOME/test_MTD_explorer_B.glabrata.local" \
  --hostid 6526 \
  --threads 20 \
  --run-number 1
```

## HPC benchmark

```bash
cd ~/MTD-Explorer

bash benchmark/run_mtd_pipeline_benchmark.sh \
  --execution-mode hpc \
  --machine LBN_Slurm \
  --dataset Bglabrata_PRJNA1306560 \
  --input /path/to/B.glabrata_fastq/samplesheet.csv \
  --output "$HOME/test_MTD_explorer_B.glabrata.hpc" \
  --hostid 6526 \
  --threads 20 \
  --hpc-conf "$HOME/MTD-Explorer/Installation/HPC/MTD_hpc_slurm.conf" \
  --run-number 1
```

The end-to-end wall time remains the principal user-facing duration. It includes
controller work, Slurm queueing, retries, node-local staging, distributed
processing and final local analyses.

The Slurm compute makespan is reported separately. It spans the earliest start
to the latest end among the collected HPC task attempts and must not be confused
with the complete pipeline wall time.

When Slurm accounting is unavailable or incomplete, the task markers still
provide node, thread, memory and elapsed-time information. Such runs are labeled
`PASS_WITH_PARTIAL_HPC_METRICS` rather than silently presenting incomplete
accounting as a fully collected benchmark.

## Distributed stages

The collector recognizes the stage names written by the current HPC backend:

```text
fastp
kraken_host
kraken_micro_raw
kraken_micro_final
bracken
humann
magicblast
```

Magic-BLAST may produce many task attempts because samples are divided into
FASTQ chunks.

## Read extraction

The standard benchmark disables microbiome-read extraction. To benchmark the
optional extraction function:

```bash
--extract --top-n 10
```

Use `--extract --top-n 0` only as an explicit stress test because every detected
taxon can require repeated scans and substantial storage.

## Main outputs

Every run contains the common files:

```text
benchmark_comparison_row.tsv
console_clean.log
diagnostic_hits.tsv
final_console_tail.txt
input_files.tsv
input_samplesheet.csv
output_extensions.tsv
output_inventory.tsv
pipeline_command.sh
pipeline_steps.tsv
pipeline_steps_raw.tsv
pipeline_summary.tsv
pipeline_summary.txt
failure_report.txt        # only when incomplete or failed
```

An HPC run additionally contains:

```text
hpc_configuration.conf
hpc_summary.tsv
hpc_jobs.tsv
hpc_tasks.tsv
hpc_stage_summary.tsv
hpc_node_summary.tsv
hpc_collection_warnings.txt   # only when accounting is partial
```

The pipeline output also retains per-attempt task markers under each Slurm
stage `success/attempts/` directory. These preserve timings and node hardware
for retries instead of keeping only the final task outcome.

### `hpc_stage_summary.tsv`

Reports, for each distributed stage:

- Slurm jobs and submitted task attempts;
- unique, successful and failed tasks;
- retries and nodes used;
- earliest submit/start and latest end;
- parallel compute makespan;
- sum of task elapsed time;
- allocated and reported consumed CPU time;
- peak task RSS and CPUs per task.

### `hpc_node_summary.tsv`

Reports only the nodes that actually executed task attempts, including CPU
model, architecture, topology, memory, kernel, current Slurm state, stages,
task counts and accumulated processing metrics. The list is discovered from
the per-attempt markers written on the selected compute nodes, with `sacct`
used as a fallback. It is never built from a fixed node list, `sinfo`, or the
set of nodes currently eligible in the partition. Consequently, nodes in
`DRAIN`, `DOWN`, or otherwise unused during that run do not appear; a node that
executed a task remains in the report even if it is drained later. Hardware
values are recorded by the real MTD jobs and supplemented, when available,
with `scontrol show node` metadata.

### `benchmark_comparison_row.tsv`

This one-row file is intended for combining repetitions and comparing local and
HPC scenarios. It includes execution mode, machine/cluster label, controller
threads, nodes used, end-to-end wall time, Slurm compute makespan, status and Git
commit.

## Fair comparisons

For local-versus-HPC documentation:

- use the same Git commit;
- use identical FASTQ files, samplesheet, databases and analysis options;
- use separate clean output directories;
- preserve both successful and failed benchmark bundles;
- report the end-to-end wall time for both scenarios;
- report Slurm compute makespan and node inventory only for HPC;
- do not compare the sum of parallel task durations directly with local wall
  time.

The bundles are created under:

```text
~/MTD_pipeline_benchmarks/
```
