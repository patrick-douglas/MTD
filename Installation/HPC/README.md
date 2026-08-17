# MTD Explorer node installation

The scripts in this directory prepare node-local software, scratch storage and
databases for the independent [MTD Explorer](https://github.com/patrick-douglas/MTD-Explorer)
[Slurm](https://slurm.schedmd.com/) backend.

The node installer creates separate environments for:

- [fastp](https://github.com/OpenGene/fastp): `MTD_fastp`;
- [Kraken2](https://github.com/DerrickWood/kraken2) and
  [Bracken](https://github.com/jenniferlu717/Bracken): `MTD_kraken2`;
- [HUMAnN](https://github.com/biobakery/humann),
  [MetaPhlAn](https://github.com/biobakery/MetaPhlAn) and
  [Magic-BLAST](https://ncbi.github.io/magicblast/): `MTD-Explorer-HPC`.

It also synchronizes the configured node-local databases and prepares
`/MTD_explorer_HPC/tmp` for stage-in, local processing and atomic stage-out.
Every complete repository-root `kraken2DB_*` directory is discovered
automatically and synchronized to the same relative path on each node.

These scripts do not clone or require
[HpcGridRunner](https://github.com/patrick-douglas/HpcGridRunner) or
[AutoHPC](https://github.com/patrick-douglas/AutoHPC).

See [`README_MTD_HPC.md`](../../README_MTD_HPC.md) for node installation,
automatic Kraken2 database discovery, additional database mappings, Slurm
policies, validation, retry behavior and usage.
