# Troubleshooting

This page collects common MTD Explorer installation, input, runtime, database,
and Slurm/HPC problems.

Start with the smallest failing layer. Do not rebuild databases or reinstall all
environments before checking the error message, the installation checker, and
the run-specific logs.

## Before reporting a problem

Record the repository state and basic system information:

```bash
cd /path/to/MTD-Explorer

git status --short
git log -1 --oneline
uname -a
```

For installation problems, also save a checker report:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --report-dir "MTD_check_$(date +%Y%m%d_%H%M%S)"
```

For a pipeline run, preserve:

- the complete terminal output;
- `methods/mtd_methods_run_parameters.csv`;
- `methods/used_samplesheet.csv`;
- `methods/used_fastq_files.tsv`;
- the exact command used to launch MTD Explorer;
- the relevant stage-specific log or Slurm log.

The saved samplesheet and FASTQ manifest are particularly useful because they
identify the exact sample design and input files used by that run.

## Installation checker fails

Run the current checker help first:

```bash
bash MTD_check_installation.sh --help
```

The recommended normal post-installation check is:

```bash
bash MTD_check_installation.sh --mode full
```

Use `quick` for runtime/database essentials and `deep` when you need cache
integrity, `kraken2-inspect`, or remote-freshness checks.

If remote access is unavailable but you still want the deep local checks:

```bash
bash MTD_check_installation.sh --mode deep --no-network
```

See [Verify the installation](../getting-started/verify-installation.md) for the
current options and exit codes.

## Conda command not found

If Conda is installed but the shell does not know the `conda` command, source
the installation explicitly. For the default Miniconda location:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
```

Then confirm:

```bash
conda --version
```

MTD Explorer records the selected Miniconda directory in `condaPath`; the
installation checker can also be pointed to it explicitly:

```bash
bash MTD_check_installation.sh \
  --conda-path /path/to/miniconda3 \
  --mode full
```

## Dedicated environment missing or wrong version

Several production tools are intentionally isolated from the legacy main `MTD`
environment.

The current checker expects, among others:

```text
MTD_fastp          fastp 1.3.6
MTD_featurecounts  Subread/featureCounts 2.1.1
MTD_kraken2        Kraken2 2.17.1 + Bracken package 3.1p1
MTD_humann         HUMAnN / MetaPhlAn runtime
```

If `bracken -v` prints `v3.0.1`, do not use that banner alone to conclude that
the environment is wrong. The validated Bioconda package is `3.1p1`, and the
current executable banner is stale; use the installation checker or Conda
package metadata to verify it.

If the checker reports one of these environments as missing or incompatible,
rerun `Install.sh` with the same persistent cache used for the original
installation, then rerun the checker.

Do not replace the production `MTD_featurecounts` executable with an older
`featureCounts` binary from the main `MTD` environment.

## FASTQ files are not found

Check the permanent input manifest after the run has started:

```text
methods/used_fastq_files.tsv
```

It records the absolute R1/R2 paths that MTD Explorer resolved for each sample.

Before starting a new run, verify that the paths or filenames represented by
the samplesheet correspond to real readable FASTQ files and that compressed
files are complete.

For details about accepted layouts and samplesheet organization, see
[Input files](../user-guide/input-files.md).

## Mixed single-end and paired-end input

One MTD Explorer run must use one sequencing layout consistently across all
samples.

With the default:

```bash
--read-layout auto
```

MTD Explorer detects the effective layout from the resolved FASTQs. A mixture of
SE and PE samples in the same run is rejected.

If you explicitly request:

```bash
--read-layout se
```

or:

```bash
--read-layout pe
```

the resolved files must match that requested layout.

## Paired FASTQs are desynchronized

For paired-end data, MTD Explorer validates R1/R2 record IDs before downstream
processing. A mismatch indicates that the mates do not represent synchronized
records.

Check that R1 and R2 came from the same library and were not independently
filtered, truncated, concatenated, or reordered.

Do not bypass this check by renaming unrelated FASTQ files.

## Output directory already exists

MTD Explorer protects an existing output directory by asking whether it should
be removed and overwritten.

If the directory contains results you need, answer `n`, move or rename the old
output, and launch the new analysis with a different `--output` path.

## Kraken2 host database problem

When `--kraken-host-db` is omitted, the host filtering database is resolved from
`--hostid`.

A custom host Kraken2 database can be selected with:

```bash
--kraken-host-db /path/to/database
```

This override changes only Kraken2 host filtering. It does not replace the host
GTF, HISAT2/Magic-BLAST resources, featureCounts annotation, OrgDb resources, or
other host-specific files associated with `--hostid`.

If a custom database is incomplete, verify the standard Kraken2 database files
and run the installation/reference checks before continuing.

## Kraken2 microbiome or Bracken database problem

A custom microbiome database can be selected with:

```bash
--kraken-micro-db /path/to/database
```

Bracken additionally requires a distribution matching the selected read length.
For example:

```text
database75mers.kmer_distrib
```

must exist for:

```bash
--bracken-read-len 75
```

If the Bracken distribution is absent, rebuild the database for the required
read length or select a read length that was actually prepared.

The default installation may remove the redundant final
`kraken2DB_micro/library/bacteria/all/` download-archive directory after the
completed Kraken2 and Bracken database has been validated. This is expected;
the consolidated bacterial `library.fna` and the persistent installation cache
are retained.

## Custom host reference is incomplete

Run or rebuild the host reference with `Create_custom_host.sh` and verify that
the genome, GTF, and protein FASTA belong to a compatible assembly/release.

The protein FASTA is required for the reference-matched eggNOG/GO functional
resources used by the current ssGSEA `auto` mode.

See [Custom host references](../user-guide/custom-host-references.md).

## ssGSEA fails after using the default options

The current default is:

```bash
--ssgsea-gmt auto
```

This is different from the legacy MSigDB C2 behavior. `auto` expects the
persistent host-specific eggNOG/GO master GMT generated by
`Create_custom_host.sh`, then creates an analysis-specific GMT after
intersecting/filtering it against `host.gct`.

If the selected host reference predates that resource, rebuild or complete the
custom-host functional annotation.

To intentionally request the legacy GMT instead, use:

```bash
--ssgsea-gmt default
```

For custom gene sets, provide a real GMT path:

```bash
--ssgsea-gmt /path/to/custom.gmt
```

See [ssGSEA outputs](../user-guide/ssgsea-outputs.md) for the current generated
files and filtering behavior.

## ssGSEA differential plots are missing

Differential ssGSEA requires at least two valid sample groups.

For multi-group runs, explicit valid `group1`/`group2` pairs in the samplesheet
are used when available; otherwise the ssGSEA plotting step builds all pairwise
contrasts. Inspect:

```text
ssGSEA/ssGSEA_differential_comparisons.tsv
```

If fewer than two groups are available, differential plotting is skipped and a
skip marker is written instead of treating the entire ssGSEA analysis as failed.

## Enhanced volcano plot labels look wrong

Current host and microbiome EnhancedVolcano reprocessing derives the group
labels from the comparison directory, for example:

```text
Liver_vs_Telencephalon
```

becomes:

```text
Liver vs Telencephalon
```

If an older output still shows generic labels, regenerate the volcano outputs
with the current MTD Explorer scripts or rerun the relevant analysis.

## HUMAnN or MetaPhlAn problem

First verify the ordinary runtime and databases:

```bash
bash MTD_check_installation.sh --mode full
```

Then inspect the HUMAnN/MetaPhlAn stage logs and confirm that the expected
HUMAnN databases and MetaPhlAn index are available.

When `--hpc-conf` is not used, this stage runs through the standard local
pipeline path. When HPC is enabled, per-sample HUMAnN/MetaPhlAn work is
submitted through Slurm and node-local resources must also pass the HPC checker.

## HPC mode does not start

HPC mode is enabled only when a valid configuration file is supplied:

```bash
--hpc-conf Installation/HPC/MTD_hpc_slurm.conf
```

Without `--hpc-conf`, local execution is expected.

The configuration file is sourced as trusted Bash and must be readable. Confirm
Slurm itself is available from the submission host:

```bash
sinfo
squeue
sbatch --version
```

Then validate the prepared nodes using the HPC checker described in
[HPC / Slurm execution](../user-guide/hpc-slurm.md).

## HPC job remains pending

Use Slurm to inspect the pending reason:

```bash
squeue -u "$USER" -o '%.18i %.12P %.24j %.8T %.10M %.6D %R'
```

Typical causes include:

- no idle eligible nodes;
- requested partition does not contain usable nodes;
- `MTD_HPC_CONSTRAINT` excludes available nodes;
- a node is `DOWN`, `DRAIN`, or otherwise unavailable;
- requested whole-node resources cannot currently be allocated.

MTD Explorer intentionally does not require a hard-coded `--nodelist` in the
standard HPC configuration; Slurm selects eligible nodes from the configured
partition/constraint.

## HPC node fails validation

Check one node directly with the HPC checker or use the node-list form described
on the HPC page.

Common failures are:

- `/MTD_explorer_HPC` missing or owned by the wrong user;
- UID/GID mismatch across machines;
- required node-local Conda environment missing;
- Kraken2/HUMAnN databases not copied to the expected node-local paths;
- node-local scratch not writable;
- insufficient free scratch space;
- passwordless SSH/sudo requirements not satisfied during installation.

Do not assume that a node is ready only because it is visible to Slurm.

## HPC database works locally but not on compute nodes

The local database path on the submission machine may differ from the
node-local path used by workers.

Databases distributed inside the standard MTD Explorer database tree are
mapped by the HPC installation/runtime layout. For external/custom database
locations, configure the longest matching source-to-node mapping in:

```bash
MTD_HPC_PATH_MAPS=(
  "/main/machine/path=/MTD_explorer_HPC/databases/custom_path"
)
```

Every node eligible for a Kraken2/Bracken task must have the corresponding
usable node-local database.

## HPC scratch space failure

The HPC backend performs a free-space preflight before expensive work when
node-local staging is enabled.

The default configuration reserves additional free space and applies
stage-specific input-size multipliers for fastp and Kraken2. If a worker reports
insufficient scratch:

1. inspect the path printed in the Slurm error log;
2. remove genuinely obsolete scratch data;
3. verify that failed-task scratch was intentionally retained for diagnosis;
4. confirm the node has enough space for the sample and stage;
5. retry after correcting the storage problem.

Do not simply set all scratch safeguards to zero without understanding the
storage requirement.

## HPC transfer appears slow

The default configuration intentionally throttles simultaneous large stage-in
and stage-out transfers. This prevents many workers from splitting a slow shared
or source link at the same time.

Current defaults allow one large stage-in and one large stage-out transfer
concurrently per stage. Computation can overlap with subsequent transfers.

On faster storage/network infrastructure, these limits can be tuned in the HPC
configuration. On a constrained link, increasing them may make total transfer
performance worse rather than better.

## HPC task fails and is retried

The backend records success/failure markers and can retry only incomplete tasks.
The standard configuration allows multiple attempts, can exclude a node that
produced a matching failure marker, and can use a final submission-node fallback
for tasks still incomplete after normal retries.

Inspect the stage-specific `hpc/` logs and retry/accounting files rather than
rerunning the entire MTD Explorer analysis immediately.

See [HPC / Slurm execution](../user-guide/hpc-slurm.md) for the retry and resume
settings.

## Find which nodes actually ran HPC tasks

Do not infer actual participation only from the configured node list or Slurm
partition. Query accounting for the submitted jobs:

```bash
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,NodeList,ExitCode
```

MTD Explorer also retains stage-level job IDs, manifests, success markers, and
logs under the HPC operational output directories.

## A pipeline stage looks idle

Sample-oriented stages print progress in a form such as:

```text
[Kraken2] [3/15] Processing sample_name
```

In HPC mode, the submitter also reports stage/task progress as Slurm state and
success-marker counts change, with periodic heartbeat output when configured.

If no progress changes, inspect the corresponding process/Slurm state rather
than assuming the pipeline has stopped.

## Collect a compact diagnostic package

For a reproducible problem report, collect at minimum:

```text
git log -1 --oneline
git status --short
MTD installation checker summary
exact MTD Explorer command
terminal log
methods/mtd_methods_run_parameters.csv
methods/used_samplesheet.csv
methods/used_fastq_files.tsv
failing stage log
```

For HPC problems, also include:

```text
HPC configuration with credentials/secrets removed
HPC checker output
relevant Slurm .out/.err logs
squeue/sacct status for the failing job IDs
```

Do not remove the original logs before the problem is understood.

## Related pages

- [System requirements](../getting-started/requirements.md)
- [Installation](../getting-started/installation.md)
- [Verify the installation](../getting-started/verify-installation.md)
- [Input files](../user-guide/input-files.md)
- [Command-line reference](../user-guide/command-line.md)
- [HPC / Slurm execution](../user-guide/hpc-slurm.md)
- [Methods and reproducibility outputs](../user-guide/methods-reproducibility-outputs.md)
