# HPC / Slurm execution

MTD Explorer provides an optional independent [Slurm][slurm] backend for
clusters with multiple Linux compute nodes. The backend distributes selected
sample-level and chunk-level stages while preserving the ordinary local
pipeline when HPC execution is not requested.

Enable the backend with:

```text
--hpc-conf FILE
```

Omitting `--hpc-conf` keeps the standard local execution path.

!!! important "HPC is optional"

    A Slurm cluster is **not** required to use MTD Explorer. Install and run the
    ordinary pipeline normally when `--hpc-conf` is not supplied.

## What the HPC backend distributes

The current backend can submit these stages through Slurm:

| Stage | Current Slurm execution model |
|---|---|
| [fastp][fastp] | One complete SE or PE sample per exclusive node |
| [Kraken2][kraken2] host filtering | One sample per exclusive node |
| Kraken2 raw microbiome classification | One sample per exclusive node |
| Kraken2 final microbiome classification | One sample per exclusive node |
| [Bracken][bracken] | One sample per exclusive node for phylum, genus, and species estimation |
| [HUMAnN][humann] / [MetaPhlAn][metaphlan] | One sample per exclusive node |
| [Magic-BLAST][magic-blast] | Synchronized FASTQ chunks distributed across CPUs and merged after validation |

The remaining pipeline stages continue on the machine orchestrating the MTD
Explorer run. In particular, the optional contaminant/read-extraction step
between raw and final microbiome Kraken2 classification currently remains
local; the final Kraken2 pass can then be distributed again.

The scientific output locations are kept compatible with local execution.
Operational Slurm manifests, logs, success markers, accounting records, and
retry information are stored under the run-specific `hpc/` directory.

## Execution model

For fastp, Kraken2, Bracken, and HUMAnN/MetaPhlAn, a Slurm array task receives
one sample and runs it on one allocated node. Whole-node stages are designed to
use heterogeneous nodes without assuming that every node has the same CPU or
memory capacity.

After allocation, the worker detects the resources exposed by Slurm and uses
node-local software, databases, and scratch space. Inputs are staged to local
scratch, processing is performed locally, and validated final outputs are
copied back atomically to the shared MTD Explorer output directory.

This design reduces sustained FASTQ and temporary-file traffic over shared
network storage.

Magic-BLAST uses a different policy: input FASTQs are split into synchronized
chunks, each chunk is submitted with a smaller resource request, and the SAM
chunks are validated and merged after completion.

## Requirements

In addition to the ordinary MTD Explorer installation, the HPC backend expects:

- [Slurm][slurm] on the cluster;
- x86_64 Linux compute nodes;
- SSH access from the installation/submission machine to the target nodes;
- the same non-root user UID and GID on the submission machine and compute
  nodes;
- a writable node-local prefix, `/MTD_explorer_HPC` by default;
- enough node-local scratch space for staged FASTQ files and temporary outputs;
- the MTD Explorer repository, configuration, work directory, and final output
  directory to be visible where required by the Slurm jobs;
- the required databases synchronized to the same node-local layout on every
  eligible compute node.

Raw FASTQ inputs that are not visible through the shared filesystem can be
pulled from the submission host when remote-input staging is enabled in the HPC
configuration.

!!! note "Heterogeneous clusters"

    The backend does not require every node to have the same CPU count or RAM.
    Whole-node jobs detect the resources made available by Slurm after
    allocation. The important requirement is a consistent software and database
    layout across eligible nodes.

## Repository components

The HPC installation and validation files are stored under:

```text
Installation/HPC/
├── MTD_install_HPC_nodes.sh
├── MTD_check_HPC_nodes.sh
├── MTD-Explorer-HPC.yml
├── MTD_hpc_slurm.conf
└── examples/
    ├── MTD_hpc_slurm.conf
    └── nodes.txt
```

Runtime helpers are stored under:

```text
aux_scripts/hpc/
├── mtd_hpc_common.sh
├── mtd_hpc_submit_array.sh
├── mtd_hpc_array_task.sh
├── mtd_hpc_node_job.sh
├── mtd_hpc_fastp_stage.sh
├── mtd_hpc_kraken_stage.sh
├── mtd_hpc_bracken_stage.sh
├── mtd_hpc_humann_stage.sh
├── mtd_hpc_magicblast_stage.sh
├── mtd_split_fastq.py
└── mtd_merge_sam_chunks.py
```

The HPC backend is implemented directly in MTD Explorer. It does not require
HpcGridRunner or AutoHPC as runtime dependencies.

## Node-local installation layout

The default prefix is:

```text
/MTD_explorer_HPC/
```

A prepared node normally contains:

```text
/MTD_explorer_HPC/
├── miniconda3/
├── envs/
│   ├── MTD-Explorer-HPC/
│   ├── MTD_fastp/
│   └── MTD_kraken2/
├── databases/
│   └── MTD-Explorer/
├── config/
├── cache/
├── logs/
└── tmp/
```

`MTD-Explorer-HPC` contains the HPC runtime used by HUMAnN, MetaPhlAn, and
Magic-BLAST. fastp and Kraken2/Bracken use their dedicated node-local
environments.

## Install the HPC runtime on one node

Complete the ordinary MTD Explorer installation first. Then, from the repository
root, prepare a node with:

```bash
cd "$HOME/MTD-Explorer"

bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node node01 \
  --user your_user
```

The selected account owns the HPC prefix and should be the same account used to
run MTD Explorer and submit Slurm jobs.

When the installer is launched from a `sudo -i` root shell, `--user` is
mandatory. Running the MTD Explorer workload itself as root is not supported.

## Install several nodes

Create a node list with one SSH/Slurm hostname per line:

```text
node01
node02
node03
```

Then run:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user your_user
```

The installer updates the node-local environments and synchronizes databases
incrementally with [rsync][rsync].

## Node installer options

The current installer supports the following main controls:

| Option | Purpose |
|---|---|
| `--node HOST` | Prepare one node |
| `--node-list FILE` | Prepare all nodes listed in a file |
| `--user USER` | Non-root account that owns and runs the HPC installation |
| `--prefix DIR` | Override the default `/MTD_explorer_HPC` prefix |
| `--repo-root DIR` | Explicit MTD Explorer repository root |
| `--ssh-user USER` | Override the SSH account used for remote access |
| `--remote-root-mode MODE` | Select the remote privilege mode |
| `--database SOURCE=RELATIVE` | Synchronize an additional database; may be repeated |
| `--database-distribution direct|propagate` | Select how databases are copied to multiple nodes |
| `--skip-default-databases` | Do not copy the default database set |
| `--force-recreate-env` | Recreate node-local Conda environments |
| `--dry-run` | Validate the proposed installation without changing nodes |

In the default remote mode, the installer connects as the selected user and
uses non-interactive `sudo -n` for privileged operations. Clusters configured
for direct root SSH can use the supported direct-root mode.

## Synchronize databases

When present, the installer can copy the default repository databases used by
the HPC stages, including:

- `HUMAnN/ref_database/`;
- `kraken2DB_micro/`;
- directories matching `*_blastdb` or `blastdb_*`.

Additional host or custom Kraken2 databases can be added with repeated
`--database` arguments. For example:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user your_user \
  --database "$HOME/MTD-Explorer/kraken2DB_6526=MTD-Explorer/kraken2DB_6526"
```

A repository database is mapped automatically to the corresponding relative
path under:

```text
/MTD_explorer_HPC/databases/MTD-Explorer/
```

### Databases outside the repository

For a database stored elsewhere on the submission machine, synchronize it with
an explicit destination:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user your_user \
  --database "/data/custom_db=custom/custom_db"
```

Then map the source path to the node-local path in
`Installation/HPC/MTD_hpc_slurm.conf`:

```bash
MTD_HPC_PATH_MAPS=(
  "/data/custom_db=/MTD_explorer_HPC/databases/custom/custom_db"
)
```

When more than one mapping matches, the longest matching source prefix is used.

For Bracken, every target node must also contain the distribution file matching
the selected read length, for example:

```text
database75mers.kmer_distrib
```

### Direct versus propagated database distribution

`--database-distribution direct` copies the required databases to target nodes
using the configured direct transfer path.

`--database-distribution propagate` can distribute already synchronized data
from prepared nodes to additional nodes, reducing repeated transfer pressure on
a single source. Propagation requires a node list so that the installer can
coordinate the participating nodes.

## Dry-run the node installation

Use `--dry-run` before changing a new cluster:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user your_user \
  --dry-run
```

The dry run checks items such as remote access, privilege handling,
architecture, UID/GID consistency, and the proposed installation state without
performing the full installation.

## Validate HPC nodes

After installation, validate the nodes with:

```bash
bash Installation/HPC/MTD_check_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user your_user
```

The checker validates the node architecture, ownership, node-local scratch,
Miniconda, the three HPC environments, Python, HUMAnN, MetaPhlAn, Magic-BLAST,
fastp, Kraken2, and Bracken. It also inspects node-local Kraken2 databases and
reports available Bracken read-length distributions.

Do not substitute this checker for the ordinary installation checker. The two
validate different parts of the installation:

```text
MTD_check_installation.sh          ordinary/local MTD Explorer installation
Installation/HPC/MTD_check_HPC_nodes.sh   node-local HPC runtime
```

## Configure Slurm

The main configuration file is:

```text
Installation/HPC/MTD_hpc_slurm.conf
```

Start from the provided configuration and review it for the target cluster.
At minimum, the partition normally needs to be changed from the example/default
value to the cluster partition that contains MTD-ready nodes.

Important current defaults include:

| Setting | Default | Purpose |
|---|---:|---|
| `MTD_HPC_SCHEDULER` | `slurm` | Scheduler backend |
| `MTD_HPC_PREFIX` | `/MTD_explorer_HPC` | Node-local installation prefix |
| `MTD_HPC_ENV_DIR` | `${MTD_HPC_PREFIX}/envs/MTD-Explorer-HPC` | HUMAnN/MetaPhlAn/Magic-BLAST node-local environment |
| `MTD_HPC_FASTP_ENV_DIR` | `${MTD_HPC_PREFIX}/envs/MTD_fastp` | fastp node-local environment |
| `MTD_HPC_KRAKEN2_ENV_DIR` | `${MTD_HPC_PREFIX}/envs/MTD_kraken2` | Kraken2/Bracken node-local environment |
| `MTD_HPC_DATABASE_ROOT` | `${MTD_HPC_PREFIX}/databases` | Node-local database root |
| `MTD_HPC_MTD_DATABASE_ROOT` | `${MTD_HPC_DATABASE_ROOT}/MTD-Explorer` | Node-local mirror root for repository databases |
| `MTD_HPC_HUMANN_DB_ROOT` | `${MTD_HPC_MTD_DATABASE_ROOT}/HUMAnN/ref_database` | HUMAnN database root on compute nodes |
| `MTD_HPC_METAPHLAN_INDEX` | `mpa_vJun23_CHOCOPhlAnSGB_202403` | MetaPhlAn index expected by the HPC runtime |
| `MTD_HPC_PARTITION` | `debug` | Slurm partition; review for the target cluster |
| `MTD_HPC_ACCOUNT` | empty | Optional Slurm account |
| `MTD_HPC_QOS` | empty | Optional Slurm QoS |
| `MTD_HPC_CONSTRAINT` | empty | Optional Slurm feature/constraint for MTD-ready nodes |
| `MTD_HPC_TIME` | `48:00:00` | Default job time limit |
| `MTD_HPC_MAX_PARALLEL` | `100` | Maximum concurrent array tasks |
| `MTD_HPC_POLL_SECONDS` | `20` | Internal Slurm/success-marker polling interval |
| `MTD_HPC_PROGRESS_ON_CHANGE` | `1` | Print progress immediately when task-state counts change |
| `MTD_HPC_PROGRESS_HEARTBEAT_SECONDS` | `300` | Print an unchanged-stage heartbeat after this interval; `0` disables it |
| `MTD_HPC_STAGE_LOCAL` | `1` | Stage work through node-local scratch |
| `MTD_HPC_LOCAL_SCRATCH_ROOT` | `${MTD_HPC_PREFIX}/tmp` | Node-local scratch root |
| `MTD_HPC_STAGEIN_MAX_CONCURRENT` | `1` | Cluster-wide concurrent stage-in transfers; `0` disables the throttle |
| `MTD_HPC_STAGEOUT_MAX_CONCURRENT` | `1` | Cluster-wide concurrent stage-out transfers; `0` disables the throttle |
| `MTD_HPC_TRANSFER_LOCK_POLL_SECONDS` | `1` | Retry interval while waiting for a shared transfer slot |
| `MTD_HPC_TRANSFER_LOCK_STALE_SECONDS` | `300` | Minimum age before an abandoned shared transfer lock can be reclaimed after owner validation |
| `MTD_HPC_SERIALIZE_STAGEOUT_PER_NODE` | `1` | Serialize packed same-node stage-out writes in addition to the cluster-wide throttle |
| `MTD_HPC_CLEAN_LOCAL_ON_SUCCESS` | `1` | Remove successful node-local scratch |
| `MTD_HPC_CLEAN_LOCAL_ON_FAILURE` | `0` | Retain failed node-local scratch for debugging |
| `MTD_HPC_SCRATCH_RESERVE_GB` | `10` | Free-space reserve before staging |
| `MTD_HPC_FASTP_SCRATCH_MULTIPLIER` | `3` | fastp scratch-space estimate multiplier |
| `MTD_HPC_KRAKEN_SCRATCH_MULTIPLIER` | `3` | Kraken2 scratch-space estimate multiplier |
| `MTD_HPC_REMOTE_INPUT_FROM_SUBMIT_HOST` | `1` | Fetch raw FASTQs over SSH when they are visible only on the submission host |
| `MTD_HPC_SUBMIT_SSH_CONNECT_TIMEOUT` | `10` | SSH connection timeout for submission-host input fallback |
| `MTD_HPC_MAGICBLAST_CHUNK_READS` | `1000000` | Synchronized SE/PE FASTQ records per Magic-BLAST chunk; `0` disables splitting |
| `MTD_HPC_RESUME` | `1` | Reuse completed validated tasks on rerun when command hashes still match |
| `MTD_HPC_MAX_ATTEMPTS` | `3` | Total normal executions allowed for each incomplete task |
| `MTD_HPC_RETRY_DELAY_SECONDS` | `20` | Delay before resubmitting only incomplete tasks |
| `MTD_HPC_RETRY_EXCLUDE_FAILED_NODES` | `1` | Exclude nodes associated with matching failures during later normal attempts |
| `MTD_HPC_FINAL_SUBMIT_NODE_FALLBACK` | `1` | Enable final launch-host compute-node fallback when the submission host maps to an eligible Slurm node |
| `MTD_HPC_FINAL_SUBMIT_NODE_ATTEMPTS` | `1` | Number of final fallback executions |

Stage-specific `sbatch` argument arrays are also part of the configuration. By
default, fastp, Kraken2/Bracken, and HUMAnN/MetaPhlAn request one exclusive
whole node with `--mem=0`. Magic-BLAST instead requests one task, one CPU, and
8 GiB per chunk so Slurm can pack independent chunks according to node
capacity. `MTD_HPC_SBATCH_EXTRA_ARGS` remains as the generic fallback policy.

`MTD_HPC_PATH_MAPS` is an array of optional source-to-node-local database
mappings. Each entry has the form:

```text
MAIN_MACHINE_PREFIX=NODE_LOCAL_PREFIX
```

The longest matching source prefix wins.

!!! warning "The configuration is executable Bash"

    `MTD_hpc_slurm.conf` is sourced by the runtime and must be treated as trusted
    code. Do not use a configuration file from an untrusted source.

## Slurm node selection

The backend does not hard-code a node list for normal pipeline execution. When
no node list is forced by the Slurm configuration, the scheduler may select any
eligible idle node in the configured partition.

This means that nodes placed in states such as `DRAIN` or otherwise unavailable
to Slurm are naturally excluded from new allocations.

If the partition also contains nodes that do not have the MTD Explorer HPC
runtime, assign a Slurm feature to prepared nodes and use, for example:

```bash
MTD_HPC_CONSTRAINT="mtd_hpc"
```

The documentation and runtime therefore do not assume specific hostnames or a
fixed cluster size.

## Resource policies

The default whole-node stages use policies equivalent to:

```bash
--exclusive
--nodes=1
--mem=0
```

This policy is used for fastp, Kraken2/Bracken, and HUMAnN/MetaPhlAn through
their stage-specific Slurm argument arrays. Kraken2 and Bracken share the
Kraken/Bracken resource policy.

Magic-BLAST is chunk-oriented and uses smaller packed tasks; the default policy
allocates one CPU and 8 GiB of memory per chunk.

Because the worker determines the allocated resources on the compute node,
`--threads` should not be interpreted as a requirement that every HPC node has
the same number of CPUs. `--threads` continues to control local pipeline work
and other stages where the main process uses that setting.

## Scratch-space preflight

Before large fastp or Kraken2 stage-ins, the worker checks available node-local
space. Current defaults are:

```bash
MTD_HPC_SCRATCH_RESERVE_GB=10
MTD_HPC_FASTP_SCRATCH_MULTIPLIER=3
MTD_HPC_KRAKEN_SCRATCH_MULTIPLIER=3
```

For a Kraken2 task with 12 GiB of combined input, the default estimate therefore
requires approximately 46 GiB free before stage-in: three times the input size
plus the 10 GiB reserve.

Adjust these safeguards only after validating the real storage behavior of the
cluster.

## Transfer throttling

Large stage-in and stage-out operations are globally throttled by default:

```bash
MTD_HPC_STAGEIN_MAX_CONCURRENT=1
MTD_HPC_STAGEOUT_MAX_CONCURRENT=1
```

The backend uses an atomic directory-based locking strategy suitable for shared
filesystems to prevent many nodes from simultaneously saturating a constrained
network link. Lock polling and stale-lock recovery are configurable.

## Run MTD Explorer locally

The ordinary command does not change:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output results \
  --hostid 6526 \
  --blast
```

## Run MTD Explorer through Slurm

Add the HPC configuration:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output results \
  --hostid 6526 \
  --blast \
  --hpc-conf Installation/HPC/MTD_hpc_slurm.conf
```

The pipeline validates the configuration before the HPC stages are launched.
The selected configuration path is also shown in the opening parameter summary.

## HPC run outputs

When HPC execution is enabled, operational files are written below:

```text
OUTPUT_DIRECTORY/hpc/
```

Stage-specific subdirectories can include:

```text
hpc/
├── fastp/
├── kraken_host/
├── kraken_micro_raw/
├── kraken_micro_final/
├── bracken/
└── magicblast/
```

HUMAnN HPC controller files are kept with the HUMAnN working structure.
Exact internal filenames may evolve, but the HPC directories are intended to
retain the information needed to diagnose distributed execution, including
input manifests, task logs, job IDs, success/failure markers, accounting data,
and retry history.

These are operational records. The main scientific output files remain in the
same analysis directories used by local execution.

## Resume behavior

With:

```bash
MTD_HPC_RESUME=1
```

rerunning the same stage can reuse tasks whose command hash still matches and
whose expected outputs pass validation. Incomplete, interrupted, invalid, or
missing tasks are resubmitted rather than blindly accepting stale output.

This behavior is specific to the HPC stage controller and should not be confused
with the main pipeline's interactive handling of an already existing top-level
output directory.

## Retry and fallback behavior

The default controller allows three normal attempts for incomplete tasks.
When enabled, nodes associated with matching task failures are excluded from
later normal retry rounds.

After the normal attempts, the backend can request one final Slurm attempt on
the compute node corresponding to the machine that launched MTD Explorer. The
hostname is discovered at runtime rather than hard-coded.

If the launch host is not registered as an eligible compute node in the target
partition, that fallback is skipped safely and the stage reports its final
failure state.

Successful scratch is removed by default. Failed scratch is retained by default
so that the task can be inspected.

## Monitor a run

Use ordinary Slurm commands while a run is active:

```bash
squeue -u "$USER"
```

For completed or failed jobs, inspect accounting information with:

```bash
sacct -u "$USER" --starttime today
```

Also inspect the run-specific `hpc/` directory and the stage log files. The MTD
Explorer controller reports stage progress and periodically checks submitted
jobs until all required outputs have either passed validation or reached a
final failure.

## Common HPC problems

### A job remains pending

Check:

```bash
squeue -u "$USER"
```

Inspect the Slurm pending reason, configured partition, account/QOS, resource
policy, and optional `MTD_HPC_CONSTRAINT`.

### A node is missing an environment or database

Run the HPC node checker again:

```bash
bash Installation/HPC/MTD_check_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user your_user
```

Then synchronize or reinstall only the missing component.

### A custom database works locally but not on a node

Confirm both:

1. the database was copied to every eligible node; and
2. the source path has the correct `MTD_HPC_PATH_MAPS` mapping when it is
   outside the repository.

For Bracken, also verify the read-length-specific `.kmer_distrib` file.

### Stage-in reports insufficient scratch space

Check node-local free space and the configured reserve/multiplier values.
Do not simply disable the safeguard without confirming how much temporary space
the stage actually requires.

### A task repeatedly fails on one node

Inspect the corresponding task log and retained failed scratch. The normal
retry policy can exclude failed nodes from later attempts, but repeated failures
should still be investigated before production runs.

## Recommended deployment sequence

For a new cluster:

1. complete the ordinary MTD Explorer installation;
2. prepare the intended compute nodes with `MTD_install_HPC_nodes.sh`;
3. synchronize all databases required by the planned analyses;
4. run `MTD_check_HPC_nodes.sh` across the prepared nodes;
5. review `MTD_hpc_slurm.conf`, especially partition and any site-specific
   account/QOS/constraint settings;
6. run a small end-to-end smoke test;
7. inspect `sacct`, task logs, stage outputs, and node-local cleanup;
8. only then launch the full dataset.

## Local versus HPC execution

| Behavior | Local mode | HPC mode |
|---|---|---|
| Activation | Default | Add `--hpc-conf FILE` |
| Slurm required | No | Yes |
| fastp | Main machine | Slurm nodes |
| Kraken2 host/microbiome | Main machine | Slurm nodes |
| Bracken | Main machine | Slurm nodes |
| HUMAnN/MetaPhlAn | Main machine | Slurm nodes |
| Magic-BLAST | Main machine when selected | Slurm chunk jobs when selected |
| Other analysis stages | Main machine | Orchestrator/main machine |
| Node-local scratch | Not required by the HPC backend | Used by default |
| Extra operational records | No `hpc/` controller tree | `hpc/` manifests/logs/accounting/retry records |
| Scientific output locations | Standard MTD Explorer paths | Same standard paths |

## Related pages

- [System requirements](../getting-started/requirements.md)
- [Installation](../getting-started/installation.md)
- [Verify installation](../getting-started/verify-installation.md)
- [Quick start](../getting-started/quick-start.md)
- [Command-line reference](command-line.md)
- [Output files](output-files.md)
- [Methods and reproducibility outputs](methods-reproducibility-outputs.md)

[slurm]: https://slurm.schedmd.com/
[fastp]: https://github.com/OpenGene/fastp
[kraken2]: https://github.com/DerrickWood/kraken2
[bracken]: https://github.com/jenniferlu717/Bracken
[humann]: https://github.com/biobakery/humann
[metaphlan]: https://github.com/biobakery/MetaPhlAn
[magic-blast]: https://ncbi.github.io/magicblast/
[rsync]: https://rsync.samba.org/
