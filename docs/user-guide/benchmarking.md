# Pipeline benchmarking

MTD Explorer includes an end-to-end benchmark workflow for measuring complete
pipeline runs in either ordinary local mode or the optional [Slurm][slurm] HPC
mode.

The benchmark is designed for reproducible performance comparisons. It keeps
large scientific outputs in the normal MTD Explorer output directory and writes
a separate compact benchmark bundle containing commands, timings, resource
summaries, input records, output inventories, and diagnostics.

!!! important "Benchmark the complete pipeline"

    For local-versus-HPC comparisons, use the benchmark wrapper rather than
    timing individual commands manually. The wrapper records the complete
    user-facing wall time and, for HPC runs, collects Slurm-specific metrics
    separately.

## Main benchmark command

The pipeline benchmark is launched with:

```bash
bash benchmark/run_mtd_pipeline_benchmark.sh [options]
```

Display the current options with:

```bash
bash benchmark/run_mtd_pipeline_benchmark.sh --help
```

The benchmark requires:

| Option | Meaning |
|---|---|
| `--machine NAME` | Label for the machine or cluster being measured |
| `--input FILE` | Samplesheet used by the pipeline |
| `--output DIR` | Clean MTD Explorer scientific output directory |
| `--hostid TAXID` | Host NCBI Taxonomy ID |

Frequently used additional options include:

| Option | Default | Meaning |
|---|---:|---|
| `--dataset LABEL` | samplesheet parent directory | Human-readable dataset label |
| `--threads INT` | `nproc` | Controller thread setting passed to `MTD_explorer.sh` |
| `--run-number INT` | `1` | Repetition number |
| `--read-layout MODE` | `auto` | `auto`, `se`, or `pe` |
| `--analysis-mode MODE` | `auto` | `auto`, `comparison`, or `exploratory` |
| `--interval SECONDS` | `5` | Controller resource-sampling interval |
| `--benchmark-root DIR` | `$HOME/MTD_pipeline_benchmarks` | Root for compact benchmark bundles |
| `--execution-mode MODE` | `auto` | `auto`, `local`, or `hpc` |
| `--hpc-conf FILE` | none | Slurm configuration used for an HPC benchmark |
| `--extract` | disabled | Enable optional detected-microbiome read extraction |
| `--top-n INT` | `10` when extraction is enabled | Number of ranked taxa extracted; `0` means all detected taxa |
| `--hisat2` | Magic-BLAST enabled | Use HISAT2 instead of Magic-BLAST |
| `--yes` | disabled | Remove an existing benchmark pipeline output without prompting |

Additional MTD Explorer arguments can be passed after `--`:

```bash
bash benchmark/run_mtd_pipeline_benchmark.sh \
  --machine workstation \
  --input /path/to/samplesheet.csv \
  --output "$HOME/mtd_benchmark_output" \
  --hostid 9606 \
  -- \
  --metadata /path/to/metadata.csv \
  --pdm spearman
```

## Execution modes

The benchmark has three execution modes:

```text
auto
local
hpc
```

The resolution rules are:

```text
auto without --hpc-conf  -> local
auto with --hpc-conf     -> hpc
local                     -> rejects --hpc-conf
hpc                       -> requires --hpc-conf
```

Therefore, an existing benchmark command that does not specify
`--execution-mode` or `--hpc-conf` continues to run locally.

### Meaning of `--threads`

`--threads` always records and passes the controller thread setting to
`MTD_explorer.sh`.

In **local mode**, this is also the main pipeline thread count.

In **HPC mode**, distributed Slurm tasks use the resources allocated or detected
on the compute node according to the HPC configuration. The benchmark therefore
does not assume that `--threads 20`, for example, means that every heterogeneous
compute node receives 20 CPUs.

## Local benchmark

Example:

```bash
cd ~/MTD-Explorer

bash benchmark/run_mtd_pipeline_benchmark.sh \
  --execution-mode local \
  --machine workstation \
  --dataset example_dataset \
  --input /path/to/samplesheet.csv \
  --output "$HOME/MTD_benchmark_local" \
  --hostid 9606 \
  --threads 20 \
  --run-number 1
```

Because `auto` resolves to local mode when no HPC configuration is supplied,
`--execution-mode local` may be omitted when an explicit mode is not needed.

## HPC benchmark

First prepare and validate the HPC runtime as described in
[HPC / Slurm execution](hpc-slurm.md).

Then run the same scientific analysis with an HPC configuration:

```bash
cd ~/MTD-Explorer

bash benchmark/run_mtd_pipeline_benchmark.sh \
  --execution-mode hpc \
  --machine cluster_label \
  --dataset example_dataset \
  --input /path/to/samplesheet.csv \
  --output "$HOME/MTD_benchmark_hpc" \
  --hostid 9606 \
  --threads 20 \
  --hpc-conf "$HOME/MTD-Explorer/Installation/HPC/MTD_hpc_slurm.conf" \
  --run-number 1
```

The HPC collector recognizes the distributed stage names produced by the
current backend:

```text
fastp
kraken_host
kraken_micro_raw
kraken_micro_final
bracken
humann
magicblast
```

Magic-BLAST may contribute many task attempts because synchronized FASTQ input
can be split into multiple chunks before submission.

## Which time should be reported?

For both local and HPC scenarios, the primary comparison metric is the
**end-to-end pipeline wall time**.

For an HPC run, this duration includes:

- controller-side work;
- Slurm queueing;
- retries;
- input stage-in;
- distributed processing;
- output stage-out;
- final pipeline stages that remain local.

The benchmark also reports a **Slurm compute makespan**. This spans the earliest
start to the latest end among collected HPC task attempts. It is useful for
understanding the distributed part of the run, but it is not equivalent to the
complete pipeline wall time.

!!! warning "Do not compare summed parallel task time with local wall time"

    Multiple HPC tasks execute concurrently. The sum of their individual
    elapsed times is therefore not a valid replacement for end-to-end wall
    time when comparing local and distributed execution.

## Read-extraction benchmark behavior

The standard pipeline benchmark **does not extract detected microbiome reads by
default**.

To include the optional extraction feature:

```bash
--extract --top-n 10
```

When extraction is enabled and `--top-n` is not changed, the benchmark uses the
top 10 taxa from the absolute ranking.

Use:

```bash
--extract --top-n 0
```

only as an intentional stress test. Extracting every detected taxon can require
many repeated FASTQ scans, substantially more runtime, and additional storage.

## Benchmark bundle location

By default, compact benchmark bundles are stored under:

```text
~/MTD_pipeline_benchmarks/
```

The scientific pipeline output remains in the directory supplied with
`--output`.

This separation allows the benchmark bundle to remain small enough for
archiving and comparison while preserving the full analysis results elsewhere.

## Common benchmark outputs

Every completed benchmark bundle can contain the following common files:

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
failure_report.txt          # only when incomplete or failed
```

Key files include:

| File | Purpose |
|---|---|
| `pipeline_command.sh` | Exact benchmarked pipeline command |
| `input_samplesheet.csv` | Copy of the benchmark samplesheet |
| `input_files.tsv` | Input-file inventory |
| `pipeline_steps.tsv` | Processed pipeline-stage timing table |
| `pipeline_summary.tsv` | Machine-readable benchmark summary |
| `pipeline_summary.txt` | Human-readable benchmark summary |
| `output_inventory.tsv` | Inventory of produced scientific outputs |
| `benchmark_comparison_row.tsv` | One-row record intended for merging repetitions and comparing scenarios |
| `failure_report.txt` | Focused failure diagnostics when a run is incomplete or fails |

## Additional HPC benchmark outputs

An HPC benchmark additionally records files such as:

```text
hpc_configuration.conf
hpc_summary.tsv
hpc_jobs.tsv
hpc_tasks.tsv
hpc_stage_summary.tsv
hpc_node_summary.tsv
hpc_collection_warnings.txt   # only when collection is partial
```

The pipeline output also retains per-attempt HPC markers below the individual
stage `success/attempts/` directories. These preserve retry history, task timing,
and node hardware information instead of retaining only the final task result.

### `hpc_stage_summary.tsv`

This table summarizes each distributed stage, including:

- Slurm jobs and submitted task attempts;
- unique, successful, and failed tasks;
- retries;
- nodes used;
- earliest submission/start and latest end;
- parallel compute makespan;
- summed task elapsed time;
- allocated and reported CPU time when available;
- peak task RSS;
- CPUs per task.

### `hpc_node_summary.tsv`

This table reports **only compute nodes that actually executed task attempts in
that benchmark run**.

The node list is reconstructed from the per-attempt markers written by real MTD
Explorer jobs, with Slurm accounting used as a fallback. It is not generated
from a fixed node list, `sinfo`, or all nodes currently eligible in a partition.

Consequently:

- a node that remained `DRAIN`, `DOWN`, or otherwise unused does not appear;
- a node that executed a task remains part of the benchmark record even if its
  Slurm state changes later;
- heterogeneous node hardware can be reported according to the nodes that
  actually participated.

Node records can include CPU model, architecture, topology, memory, kernel,
Slurm state, executed stages, task counts, and accumulated processing metrics.

### `benchmark_comparison_row.tsv`

This one-row table is intended for merging benchmark repetitions. It records
fields such as:

- execution mode;
- machine or cluster label;
- controller threads;
- nodes used;
- end-to-end wall time;
- Slurm compute makespan when applicable;
- benchmark status;
- Git commit.

## Partial Slurm accounting

HPC task markers remain useful even when `sacct` is unavailable or returns
incomplete accounting information. The collector can still recover task-level
node, thread, memory, and elapsed-time information from the MTD Explorer HPC
records.

A successful run with incomplete Slurm accounting is labeled:

```text
PASS_WITH_PARTIAL_HPC_METRICS
```

rather than presenting partial metrics as though full accounting had been
collected.

Inspect `hpc_collection_warnings.txt` when that status is reported.

## Fair local-versus-HPC comparisons

Use the following rules when creating performance numbers for documentation,
reports, or manuscripts:

1. use the same MTD Explorer Git commit;
2. use identical FASTQ files;
3. use the same samplesheet and metadata;
4. use identical databases and reference versions;
5. use the same analysis options;
6. use separate clean scientific output directories;
7. preserve benchmark bundles from both successful and failed runs;
8. compare end-to-end wall time for local and HPC scenarios;
9. report Slurm compute makespan and actual node inventory separately for HPC;
10. do not compare summed parallel task durations directly against local wall
    time.

These rules are especially important for heterogeneous clusters, where the set
of compute nodes actually used can differ between repetitions as Slurm node
availability changes.

## Benchmark implementation files

The main pipeline benchmark is implemented by:

```text
benchmark/run_mtd_pipeline_benchmark.sh
benchmark/MTD_benchmark_install.sh
benchmark/MTD_pipeline_benchmark_report.py
benchmark/MTD_slurm_benchmark_report.py
```

`run_mtd_pipeline_benchmark.sh` validates the requested scenario, runs the
complete pipeline, and creates the final bundle.

`MTD_benchmark_install.sh` records controller-side wall time and resource use.

`MTD_pipeline_benchmark_report.py` creates the common local/HPC pipeline
summary, stage timings, output inventory, comparison row, and failure
information.

`MTD_slurm_benchmark_report.py` discovers MTD Explorer Slurm job histories,
combines `sacct` information with the pipeline's own task markers, and produces
stage- and node-level HPC summaries.

Normal pipeline runs are not instrumented with benchmark timestamps unless the
benchmark-specific environment variable used by the wrapper is set.

## Related pages

- [HPC / Slurm execution](hpc-slurm.md)
- [Command-line reference](command-line.md)
- [Methods and reproducibility outputs](methods-reproducibility-outputs.md)
- [Output files](output-files.md)
- [Troubleshooting](../troubleshooting/index.md)

[slurm]: https://slurm.schedmd.com/
