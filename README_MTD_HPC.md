# MTD Explorer HPC backend

The independent [Slurm](https://slurm.schedmd.com/) backend distributes the
sample-level and chunk-level stages of
[MTD Explorer](https://github.com/patrick-douglas/MTD-Explorer) across a
heterogeneous Linux cluster.

The current backend supports:

- [fastp](https://github.com/OpenGene/fastp): one complete SE or PE sample per
  exclusive node;
- [Kraken2](https://github.com/DerrickWood/kraken2) host filtering: one sample
  per exclusive node;
- Kraken2 raw and final microbiome classification: one sample per exclusive
  node;
- [Bracken](https://github.com/jenniferlu717/Bracken) phylum, genus and species
  estimation: one sample per exclusive node;
- [HUMAnN](https://github.com/biobakery/humann) and
  [MetaPhlAn](https://github.com/biobakery/MetaPhlAn): one sample per exclusive
  node;
- [Magic-BLAST](https://ncbi.github.io/magicblast/): synchronized FASTQ chunks
  packed across available CPUs and merged after validation.

Omitting `--hpc-conf` preserves the ordinary local pipeline behavior.

The backend does **not** clone, import or require
[HpcGridRunner](https://github.com/patrick-douglas/HpcGridRunner) or
[AutoHPC](https://github.com/patrick-douglas/AutoHPC). Those projects inspired
its centralized orchestration and node-adaptive design, but they are not
runtime dependencies.

## Execution model

For fastp, Kraken2, Bracken and HUMAnN, each Slurm array element processes one
sample on one exclusive node. No fixed CPU count is requested. After allocation,
the task detects the CPUs and available memory exposed by Slurm and uses the
node-local software and databases.

Inputs are copied from the shared filesystem into node-local scratch. Processing
and temporary output creation occur locally. Validated final outputs are copied
atomically back to the shared MTD output directory. This avoids sustained FASTQ
and temporary-file traffic through NFS.

The backend also provides:

- stage-level resume using command hashes, output validation and `.success`
  markers;
- three normal attempts by default, submitting only incomplete tasks;
- exclusion of nodes that produced matching task failures;
- one final Slurm fallback attempt on the compute node corresponding to the host
  that launched MTD Explorer, without hard-coding its hostname;
- task-specific logs and retained local scratch after failures;
- preflight free-space checks before copying large fastp or Kraken2 inputs.

## Repository layout

```text
Installation/HPC/
├── MTD_install_HPC_nodes.sh
├── MTD_check_HPC_nodes.sh
├── MTD-Explorer-HPC.yml
├── MTD_hpc_slurm.conf
└── examples/
    ├── MTD_hpc_slurm.conf
    └── nodes.txt

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

The default node-local structure is:

```text
/MTD_explorer_HPC/
├── miniconda3/
├── envs/
│   ├── MTD-Explorer-HPC/   # HUMAnN, MetaPhlAn and Magic-BLAST
│   ├── MTD_fastp/          # fastp
│   └── MTD_kraken2/        # Kraken2 and Bracken
├── databases/
│   └── MTD-Explorer/
├── config/
├── cache/
├── logs/
└── tmp/
```

The installer runs the privileged filesystem operations as root, but the full
prefix is owned by the non-root account selected with `--user`.

## 1. Configure nodes

After the ordinary MTD Explorer installation, prepare one node with:

```bash
cd "$HOME/MTD-Explorer"

bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node node01 \
  --user me
```

For several nodes, create a file containing one SSH/Slurm hostname per line:

```text
node01
node02
node03
```

Then run:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user me
```

The installer uses [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/)
to create or update all three node-local environments. Re-running it updates
the environments and synchronizes databases incrementally with
[rsync](https://rsync.samba.org/).

When launched from a `sudo -i` root shell, `--user` is mandatory. The selected
account must have the same UID and GID on the submission machine and all nodes.

The default remote mode connects as that user and uses non-interactive
`sudo -n`. Clusters configured for root SSH may use:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node node01 \
  --user me \
  --remote-root-mode direct
```

## 2. Synchronize databases

By default, the node installer copies these repository paths when present:

- `HUMAnN/ref_database/`;
- **every complete repository-root `kraken2DB_*` directory**, including
  `kraken2DB_micro` and all custom host databases such as `kraken2DB_6526`;
- directories matching `*_blastdb` or `blastdb_*`.

A repository-root Kraken2 directory is considered complete only when
`hash.k2d`, `opts.k2d` and `taxo.k2d` are all present and non-empty. If a
`kraken2DB_*` directory exists but is incomplete, the installer stops before
node configuration instead of silently omitting it.

For example, if the main installation contains:

```text
$HOME/MTD-Explorer/
├── kraken2DB_micro/
├── kraken2DB_6526/
└── kraken2DB_59463/
```

the installer automatically synchronizes them as:

```text
/MTD_explorer_HPC/databases/MTD-Explorer/kraken2DB_micro/
/MTD_explorer_HPC/databases/MTD-Explorer/kraken2DB_6526/
/MTD_explorer_HPC/databases/MTD-Explorer/kraken2DB_59463/
```

No repeated `--database` option is required for those repository-root Kraken2
databases.

`--database SOURCE=RELATIVE` remains available for additional databases outside
the automatically discovered repository paths. For example, for a database
outside the repository, synchronize it and declare the matching source-to-node
mapping in `Installation/HPC/MTD_hpc_slurm.conf`:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user me \
  --database "/media/me/custom_db=custom/custom_db"
```

```bash
MTD_HPC_PATH_MAPS=(
  "/media/me/custom_db=/MTD_explorer_HPC/databases/custom/custom_db"
)
```

The longest matching source prefix is used. Every eligible node must contain an
identical database layout. For Bracken, the microbiome database must also
contain the distribution matching the configured read length, for example
`database75mers.kmer_distrib`.

A dry run validates SSH, privileges, architecture, UID/GID and installation
state without changing the nodes:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user me \
  --dry-run
```

## 3. Validate nodes

```bash
bash Installation/HPC/MTD_check_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user me
```

The checker validates ownership, architecture, node-local scratch, Miniconda,
the three environments, Python, HUMAnN, MetaPhlAn, Magic-BLAST, fastp, Kraken2
and Bracken. It also validates every node-local Kraken2 database it finds and
reports available Bracken read-length distributions.

The installer records the full list of Kraken2 databases selected for
synchronization in each node manifest and also passes that list to the final
validation. Final validation therefore fails if any expected database directory,
`hash.k2d`, `opts.k2d` or `taxo.k2d` is missing or empty on a node. Later
standalone checker runs automatically reuse the requirements recorded in the
node manifest, so a deleted or incomplete required database cannot be hidden by
the presence of a different valid Kraken2 database. The checker can also enforce
a database manually with the repeatable option:

```bash
bash Installation/HPC/MTD_check_HPC_nodes.sh \
  --node node01 \
  --user me \
  --require-kraken-db MTD-Explorer/kraken2DB_6526
```

## 4. Configure Slurm

Start from the synchronized example when creating a new configuration:

```bash
cp Installation/HPC/examples/MTD_hpc_slurm.conf \
   Installation/HPC/MTD_hpc_slurm.conf
```

At minimum, review:

```bash
MTD_HPC_PARTITION="your_partition"
MTD_HPC_MAX_PARALLEL=100
```

No `--nodelist` is set by default, so Slurm may use any eligible idle node in the
partition. When that partition contains nodes not prepared for MTD Explorer,
configure a Slurm feature and set:

```bash
MTD_HPC_CONSTRAINT="mtd_hpc"
```

The whole-node stages use independent policies:

```bash
MTD_HPC_FASTP_SBATCH_EXTRA_ARGS=(
  --exclusive
  --nodes=1
  --mem=0
)

MTD_HPC_KRAKEN_SBATCH_EXTRA_ARGS=(
  --exclusive
  --nodes=1
  --mem=0
)

MTD_HPC_HUMANN_SBATCH_EXTRA_ARGS=(
  --exclusive
  --nodes=1
  --mem=0
)
```

Kraken2 and Bracken share the `MTD_HPC_KRAKEN_SBATCH_EXTRA_ARGS` policy. The
Magic-BLAST policy remains packed at one CPU and 8 GiB per chunk.

The free-space preflight defaults are:

```bash
MTD_HPC_SCRATCH_RESERVE_GB=10
MTD_HPC_FASTP_SCRATCH_MULTIPLIER=3
MTD_HPC_KRAKEN_SCRATCH_MULTIPLIER=3
```

For example, a Kraken2 PE task with 12 GiB of combined inputs requires at least
46 GiB free before stage-in: three times the input size plus the 10 GiB reserve.
These values are conservative safeguards and can be adjusted for the node-local
filesystem after observing real runs.

## 5. Run MTD Explorer

Local execution remains unchanged:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output results \
  --hostid 6526 \
  --blast
```

Activate the backend with:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output results \
  --hostid 6526 \
  --blast \
  --hpc-conf Installation/HPC/MTD_hpc_slurm.conf
```

With `--hpc-conf`, fastp, Kraken2 host filtering, both microbiome Kraken2 passes,
Bracken, HUMAnN/MetaPhlAn and optional Magic-BLAST are submitted through Slurm.
The remaining pipeline stages continue on the machine orchestrating the run.
The optional contaminant extraction between the raw and final microbiome
classifications currently remains local; the final Kraken2 classification is
then distributed again.

## Restart and failure behavior

Each stage stores its manifests, logs, success markers, job IDs, accounting
records and retry history below the MTD output directory under `hpc/`.

With `MTD_HPC_RESUME=1`, rerunning the same output directory reuses only tasks
whose command hash matches and whose expected outputs still pass validation.
Tasks that failed, were interrupted or lost outputs are submitted again.

By default, each incomplete task has three normal attempts. A node recorded in a
matching `.failed` marker is excluded from later normal rounds. After those
attempts, one final round is forced through Slurm onto the compute node that
corresponds to the machine that launched MTD Explorer. If that host is not a
registered compute node in the configured partition, the fallback is skipped
safely and the stage reports a definitive failure.

## Limitations

- The complete runtime currently targets **x86_64 Linux nodes**.
- Every node submitted to these stages must have the same local prefix,
  environments and database layout.
- The repository, input files, configuration, work directories and final output
  directories must be visible to every node through the shared filesystem.
- The configuration is sourced as trusted [GNU Bash](https://www.gnu.org/software/bash/)
  and must not come from an untrusted source.
- Before production use, validate the real Conda solves, database copies and a
  small end-to-end Slurm run on the target cluster.

## Design lineage and acknowledgements

The backend was designed from operational patterns previously used in
HpcGridRunner and AutoHPC: centralized command distribution, whole-node use on
heterogeneous hardware, node-local resource detection and orchestration from a
shared filesystem. No source code or runtime component from either repository
is required.
