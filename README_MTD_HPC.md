# MTD Explorer HPC backend

This add-on provides an independent [Slurm](https://slurm.schedmd.com/) backend
for the expensive [HUMAnN](https://github.com/biobakery/humann) / [MetaPhlAn](https://github.com/biobakery/MetaPhlAn)
and [Magic-BLAST](https://ncbi.github.io/magicblast/) stages of
[MTD Explorer](https://github.com/patrick-douglas/MTD-Explorer).

It does **not** clone, import, or require
[HpcGridRunner](https://github.com/patrick-douglas/HpcGridRunner) or
[AutoHPC](https://github.com/patrick-douglas/AutoHPC). Those projects inspired
the whole-node, node-adaptive design and should be acknowledged in the project
documentation, but they are not runtime dependencies.

## What this first implementation does

- `--hpc-conf FILE` activates HPC execution; omitting it preserves local mode.
- One Slurm array task is created per HUMAnN sample.
- Magic-BLAST FASTQs are split into synchronized chunks, distributed as array
  tasks, validated, and merged back into one SAM per sample.
- Each task receives one exclusive node, detects its own CPUs and available RAM,
  and exports `MTD_NODE_THREADS`, `MTD_NODE_MEMORY_KB`, and a 90% RAM budget.
- Matching success markers and non-empty outputs permit stage-level restart
  without repeating completed work while the HPC work directory is preserved.
- Inputs, task manifests, logs, and results remain under the shared NFS-visible
  MTD working paths. Software and databases are read from node-local storage.

## Repository layout

```text
Installation/HPC/
├── MTD_install_HPC_nodes.sh
├── MTD_check_HPC_nodes.sh
├── MTD-Explorer-HPC.yml
└── examples/
    ├── MTD_hpc_slurm.conf
    └── nodes.txt

aux_scripts/hpc/
├── mtd_hpc_common.sh
├── mtd_hpc_submit_array.sh
├── mtd_hpc_array_task.sh
├── mtd_hpc_node_job.sh
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
│   └── MTD-Explorer-HPC/
├── databases/
│   └── MTD-Explorer/
├── config/
├── cache/
├── logs/
└── tmp/
```

The installer runs with root privileges, but the complete prefix is owned by the
non-root account selected with `--user`.

## 1. Add the package to the repository

Extract the downloaded package and run:

```bash
cd /path/to/MTD-Explorer-HPC-addon

bash apply_to_MTD_Explorer.sh \
  --repo "$HOME/MTD-Explorer"
```

The command:

1. checks whether the current `MTD_explorer.sh` matches the expected integration
   points;
2. copies the independent HPC files;
3. creates a timestamped backup of `MTD_explorer.sh`;
4. adds `--hpc-conf` and the two HPC execution branches;
5. runs shell/[Python](https://www.python.org/) syntax and splitter/merger smoke tests.

A preflight without modifying the repository is available:

```bash
bash apply_to_MTD_Explorer.sh \
  --repo "$HOME/MTD-Explorer" \
  --check
```

## 2. Configure one node

Run this from the main MTD Explorer machine after the ordinary `Install.sh`:

```bash
cd "$HOME/MTD-Explorer"

bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node node01 \
  --user me
```

When launched from a normal user terminal, `--user` can be omitted. The script
records the login user before requesting local `sudo`. When launched from a
`sudo -i` root shell, `--user` is mandatory and `root` is refused.

The default remote mode connects as the selected user and requires non-interactive
remote `sudo -n`. A cluster configured for root SSH can instead use:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node node01 \
  --user me \
  --remote-root-mode direct
```

The installer uses [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/)
to create `/MTD_explorer_HPC/envs/MTD-Explorer-HPC`. It installs
[rsync](https://rsync.samba.org/) through a supported node package manager only
when the node does not already provide it.

## 3. Configure several nodes

Create a text file containing one hostname per line:

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

Nodes are configured sequentially. Re-running the command updates the
[Conda](https://docs.conda.io/) environment and synchronizes databases
incrementally with rsync.

The installer verifies that the selected user has the same UID and GID on the
main machine and every node. This is important for the shared NFS home.

### Database synchronization

By default, the installer copies:

- `HUMAnN/ref_database/`;
- every repository directory matching `*_blastdb`;
- every repository directory matching `blastdb_*`.

Additional database directories can be supplied repeatedly:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node node01 \
  --user me \
  --database /main/path/custom_database=custom/custom_database
```

The target is relative to `/MTD_explorer_HPC/databases/`.

To print mutating commands without applying them:

```bash
bash Installation/HPC/MTD_install_HPC_nodes.sh \
  --node node01 \
  --user me \
  --dry-run
```

The dry run still connects to the node to validate SSH, sudo, architecture,
UID/GID, and the existing installation state.

## 4. Validate configured nodes

```bash
bash Installation/HPC/MTD_check_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user me
```

This checks ownership, node architecture, Miniconda, the environment, and the
installed versions of Python, HUMAnN, MetaPhlAn, and Magic-BLAST.

## 5. Configure Slurm

Copy the example:

```bash
cp Installation/HPC/examples/MTD_hpc_slurm.conf \
   Installation/HPC/MTD_hpc_slurm.conf
```

Edit at least:

```bash
MTD_HPC_PARTITION="your_partition"
MTD_HPC_MAX_PARALLEL=6
```

When the Slurm partition also contains unconfigured nodes, set a constraint that
selects only nodes prepared for this backend:

```bash
MTD_HPC_CONSTRAINT="mtd_hpc"
```

The default submission arguments are:

```bash
--exclusive --nodes=1 --mem=0
```

No fixed CPU count is requested. The command running inside the allocated node
uses `SLURM_CPUS_ON_NODE`, falling back to `SLURM_CPUS_PER_TASK` and then
`nproc`.

## 6. Run MTD Explorer

Without `--hpc-conf`, behavior remains local:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output results \
  --hostid 6526 \
  --blast
```

To activate HPC execution:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output results \
  --hostid 6526 \
  --blast \
  --hpc-conf Installation/HPC/MTD_hpc_slurm.conf
```

In this first implementation, HUMAnN/MetaPhlAn and Magic-BLAST use Slurm. The
remaining MTD Explorer stages continue on the machine running the main pipeline.
If `--blast` is omitted, HUMAnN/MetaPhlAn still uses the HPC backend while [HISAT2](https://daehwankimlab.github.io/hisat2/) remains local.

## Important limitations

- The complete environment currently supports **x86_64 Linux nodes**. The
  [Bioconda](https://bioconda.github.io/) Magic-BLAST 1.7.0 package used here is not published for
  `linux-aarch64`; the installer therefore stops clearly on ARM64 nodes rather
  than leaving a partially functional runtime.
- Every submitted node must have the same local prefix and synchronized database
  layout.
- The MTD repository, inputs, HPC configuration, work directory, and outputs must
  be visible from every node through the shared filesystem.
- The configuration file is sourced as trusted [GNU Bash](https://www.gnu.org/software/bash/) and must not come from an
  untrusted source.
- A real Slurm run and full Conda solve must still be validated on the target
  cluster before production use. Start with one configured node and a tiny test
  dataset.


## Design lineage and acknowledgements

The independent backend was designed from operational patterns previously used
in [HpcGridRunner](https://github.com/patrick-douglas/HpcGridRunner) and
[AutoHPC](https://github.com/patrick-douglas/AutoHPC): command distribution,
whole-node execution on heterogeneous hardware, node-local resource detection,
and centralized orchestration from a shared filesystem. No source code or
runtime component from either repository is required by this add-on.
