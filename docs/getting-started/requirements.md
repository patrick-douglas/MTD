# System requirements

This page describes the operating system, hardware, storage, and network requirements for installing and running MTD Explorer.

!!! info "Benchmark-based requirements"

    These requirements are based on a complete warm-cache installation benchmark
    performed with the full default database set. The benchmark completed
    successfully on a 20-thread x86-64 computer with 128 GiB of RAM and NVMe
    storage.

    The installation reached approximately 116.6 GiB of process-tree resident
    memory and 20.2 GiB of swap usage. Systems close to the minimum requirements
    must therefore have properly configured swap space.

## Supported operating systems

MTD Explorer is developed and tested on 64-bit GNU/Linux systems.

### Currently tested

* Linux Mint
* Ubuntu-based Linux distributions

### Not currently supported

* Native Microsoft Windows
* Native macOS

Windows users may be able to use a Linux server, virtual machine, or WSL2, but these configurations have not yet been formally validated.

## Hardware requirements

The following values apply to a complete MTD Explorer installation, including
construction of the default Kraken2 and Bracken databases.

| Resource | Minimum | Recommended | Notes |
| --- | ---: | ---: | --- |
| CPU architecture | x86-64 | x86-64 | ARM64 is not currently validated |
| CPU threads | 8 | 20 or more | Fewer threads substantially increase installation and analysis time |
| RAM | 128 GB | 192 GB or more | 256 GB is preferred for additional headroom during large database builds |
| Swap | 32 GB | 64 GB | Swap should preferably be located on SSD or NVMe storage |
| Free storage | 1.5 TB | 2 TB or more | Includes the installation, databases, Conda environments, and reusable cache |
| Installation drive | SSD | NVMe SSD | Database construction performs several terabytes of disk I/O |
| Internet connection | Required | Stable broadband | Used for packages, reference metadata, and database updates |

!!! warning "RAM and swap"

    A full database installation was successfully completed with 128 GiB of RAM,
    but the benchmark reached approximately 116.6 GiB of resident memory and
    used approximately 20.2 GiB of swap.

    A system with 128 GB of RAM should therefore have at least 32 GB of available
    swap. Systems without swap may fail during the most memory-intensive Kraken2
    or Bracken construction steps.

!!! note "Storage layout"

    The benchmark produced approximately:

    * 740 GiB for the MTD Explorer installation and generated databases;
    * 18 GiB for the Conda installation and environments;
    * 373 GiB for the reusable installation cache.

    When the cache and installation are stored on the same filesystem, at least
    1.5 TB of free space should be available before installation. A minimum of
    2 TB is recommended when input data, temporary files, and analysis outputs
    will also be stored on that filesystem.

    The reusable cache may instead be placed on a separate disk. In that
    configuration, plan for approximately 1 TB on the installation disk and
    at least 500 GB on the cache disk.

!!! note "Installation versus analysis"

    These values describe installation of the complete default database set.
    Resource requirements for individual analyses depend on the number of
    samples, FASTQ file sizes, selected databases, and enabled analysis modules.

    A separate official clean pipeline benchmark using the 15-sample
    [PRJNA1306560](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1306560)
    example dataset completed in
    **54 h 58 min** with 20 threads and reached **117.27 GiB** peak process-tree
    RSS. That analysis also produced approximately **1.39 TB** of pipeline output
    with optional top-5 microbiome-read extraction enabled. These are observed
    values for that specific workload, not minimum requirements. See
    [Pipeline benchmarking](../user-guide/benchmarking.md#reference-local-benchmark-prjna1306560)
    for the complete configuration and stage timings.

## Optional Slurm / HPC requirements

The hardware table above describes a complete ordinary MTD Explorer installation. The optional [Slurm][slurm] backend has additional cluster requirements.

For HPC execution, prepared compute nodes currently need x86_64 Linux, Slurm access, a consistent node-local MTD Explorer runtime/database layout, and writable local scratch. Nodes may have different CPU counts and RAM; whole-node workers detect the resources exposed by Slurm after allocation.

The default node-local prefix is `/MTD_explorer_HPC`. The repository, configuration, work directory, and final output directory must be visible where required by the jobs. Raw FASTQs may instead be fetched from the submission host when the remote-input staging option is enabled.

See [HPC / Slurm execution](../user-guide/hpc-slurm.md) for node preparation, database synchronization, Slurm configuration, scratch-space safeguards, and validation.

## Storage planning

Storage should be available for:

1. MTD Explorer source code;
2. Conda environments;
3. Kraken2 databases;
4. host reference genomes and indexes;
5. eggNOG and HUMAnN databases;
6. the optional reusable installation cache;
7. input FASTQ files;
8. analysis outputs and temporary files.

The reusable installation cache may be stored on a separate disk.

## Required permissions

The user should have:

* permission to write inside the MTD Explorer directory;
* permission to write to the selected cache directory;
* permission to create files in the selected output directory;
* `sudo` access when system packages need to be installed.

MTD Explorer itself should not normally be executed as the root user.

## Check your computer

Use the following commands before installation.

### Operating system

```bash
cat /etc/os-release
```

### CPU architecture

```bash
uname -m
```

### Available CPU threads

```bash
nproc
```

### Available memory

```bash
free -h
```

### Available disk space

```bash
df -h
```

### Git availability

```bash
git --version
```

## Next step

Continue to the [installation guide](installation.md).


[slurm]: https://slurm.schedmd.com/
