# Verify the installation

MTD Explorer includes a dedicated verification program that checks whether
the software environments, commands, packages, reference files, indexes,
and databases required by the pipeline are available and usable.

The verification program is:

```text
MTD_check_installation.sh
```

!!! important

    Finishing `Install.sh` does not by itself guarantee that every
    environment, package, index, and reference database is ready for use.

    Run the installation checker before analyzing a real dataset.

## Display the checker help

From the MTD Explorer directory:

```bash
cd ~/MTD-Explorer
bash MTD_check_installation.sh --help
```

The general syntax is:

```text
MTD_check_installation.sh [options]
```

## Basic verification

The checker can usually determine the MTD Explorer directory, Conda
installation, and persistent cache paths automatically.

Run the default full verification with:

```bash
cd ~/MTD-Explorer
bash MTD_check_installation.sh
```

The default verification mode is:

```text
full
```

## Verification modes

MTD Explorer provides three verification levels.

| Mode | Current scope |
|---|---|
| `quick` | Runtime and installed-database essentials |
| `full` | Everything in `quick`, plus source checks, package checks, cache metadata, HUMAnN resources, and installation audit contracts |
| `deep` | Everything in `full`, plus gzip integrity checks, `kraken2-inspect`, and safe remote freshness checks |

### Quick mode

Use quick mode for a fast runtime check:

```bash
bash MTD_check_installation.sh --mode quick
```

This is useful after script updates, permission changes, environment repairs,
or transfer of an existing installation when you primarily want to confirm that
the required runtimes and installed databases are available.

### Full mode

Full mode is the default and is recommended after a normal installation:

```bash
bash MTD_check_installation.sh --mode full
```

The following two commands are therefore equivalent:

```bash
bash MTD_check_installation.sh
```

```bash
bash MTD_check_installation.sh --mode full
```

### Deep mode

Deep mode performs the most comprehensive validation:

```bash
bash MTD_check_installation.sh --mode deep
```

It extends the full check with integrity-oriented and remote-freshness checks.
Because these checks can inspect many cached files and query remote metadata,
deep mode is most useful after a clean installation, interrupted downloads,
cache migration, or suspected cache/database corruption.

To run deep mode without remote freshness checks:

```bash
bash MTD_check_installation.sh --mode deep --no-network
```

## Checker options

| Option | Argument | Description |
|---|---|---|
| `--mtd-dir` | `PATH` | MTD Explorer repository/installation directory; default: directory containing the checker |
| `--installer` | `PATH` | Explicit installer script; default: auto-detected inside `--mtd-dir` |
| `--conda-path` | `PATH` | Miniconda directory; default: `condaPath`, then `$HOME/miniconda3` |
| `-o`, `--offline-dir` | `PATH` | Persistent installation cache; default: `offlineCachePath` |
| `-r`, `--read-length` | `INT` | Bracken read length; default: `75` |
| `--hostid` | `TAXID` | Validate one installed custom-host reference; default: auto-detect numeric host references |
| `--mode` | `quick`, `full`, `deep` | Verification level; default: `full` |
| `--no-network` | — | Skip remote freshness checks in `deep` mode |
| `--report-dir` | `PATH` | Directory where verification reports are written |
| `--strict` | — | Treat warnings as final failure |
| `--keep-temp` | — | Preserve temporary checker files |
| `--version` | — | Display the checker version and exit |
| `-h`, `--help` | — | Display the checker help and exit |

!!! note "No `-m` or `-p` checker aliases"

    In the current checker, `--mtd-dir` and `--conda-path` are long-form
    options. The short flags `-m` and `-p` belong to other MTD Explorer
    commands and should not be used here.

## Exit status

The checker currently uses these process exit codes:

| Exit code | Meaning |
|---|---|
| `0` | No failures; warnings are allowed unless `--strict` is used |
| `1` | One or more failures, or warnings when `--strict` is enabled |
| `2` | Invalid checker arguments |

## Automatic path detection

### MTD Explorer directory

When `--mtd-dir` is not provided, the checker uses the directory containing
`MTD_check_installation.sh`.

For a standard installation, this is normally sufficient:

```bash
cd ~/MTD-Explorer
bash MTD_check_installation.sh
```

A different installation directory can be specified explicitly:

```bash
bash MTD_check_installation.sh \
  --mtd-dir /path/to/MTD-Explorer
```

### Conda installation

When `--conda-path` is not provided, the checker searches for the Conda
installation in the following order:

1. the path recorded in `MTD-Explorer/condaPath`;
2. `$HOME/miniconda3`.

A path can also be supplied explicitly:

```bash
bash MTD_check_installation.sh \
  --conda-path /home/user/miniconda3
```

### Persistent installation cache

When `--offline-dir` is not provided, the checker uses the path recorded in:

```text
MTD-Explorer/offlineCachePath
```

The cache can be specified manually:

```bash
bash MTD_check_installation.sh \
  --offline-dir /path/to/MTD_install_cache
```

An explicit cache path is particularly useful when:

- the cache was moved;
- the installation directory was copied from another computer;
- more than one cache is available;
- `offlineCachePath` is missing or outdated.

## Complete explicit command

All primary paths can be supplied explicitly:

```bash
bash MTD_check_installation.sh \
  --mtd-dir /home/user/MTD-Explorer \
  --conda-path /home/user/miniconda3 \
  --offline-dir /path/to/MTD_install_cache \
  --read-length 75 \
  --mode full
```

Replace the example paths with the paths used by the current installation.

## Bracken read length

The default Bracken read length checked by the program is:

```text
75
```

A different value can be selected with:

```bash
bash MTD_check_installation.sh \
  --read-length 100
```

The value should match the Bracken read length used during installation and
database preparation.

For example, when the installer was run with:

```bash
bash Install.sh \
  -o /path/to/MTD_install_cache \
  -r 100
```

the checker should also use:

```bash
bash MTD_check_installation.sh \
  --read-length 100
```

## Dedicated runtime versions checked

The current checker explicitly validates several isolated production runtimes:

| Environment | Package | Expected version |
|---|---|---|
| `MTD_fastp` | fastp | `1.3.6` |
| `MTD_featurecounts` | Subread / featureCounts | `2.1.1` |
| `MTD_kraken2` | Kraken2 | `2.17.1` |
| `MTD_kraken2` | Bracken package | `3.1p1` |

The Bracken executable included by the validated package can still print an
upstream `v3.0.1` banner. The checker intentionally verifies the Conda package
metadata (`3.1p1`) instead of using that stale banner as the package version.

Production host quantification uses `MTD_featurecounts`; an older
`featureCounts` executable may remain in the main `MTD` environment for
legacy dependency compatibility, but it is not the production quantification
runtime.

If a required isolated runtime is missing or has an incompatible version,
rerun `Install.sh` and repeat the full checker.

## HPC validation is separate

`MTD_check_installation.sh` validates the ordinary MTD Explorer installation. It does not replace validation of node-local Slurm runtimes.

For HPC nodes, use:

```bash
bash Installation/HPC/MTD_check_HPC_nodes.sh \
  --node-list Installation/HPC/examples/nodes.txt \
  --user your_user
```

See [HPC / Slurm execution](../user-guide/hpc-slurm.md) for the complete HPC validation workflow.

## Report directory

A specific output directory can be selected with `--report-dir`:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --report-dir ./MTD_installation_check
```

Using a dedicated report directory is recommended when:

- testing a clean installation;
- comparing different computers;
- preparing benchmark records;
- reporting an installation problem;
- preserving results from multiple checker runs.

A timestamped report directory can be created with:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --report-dir "MTD_check_$(date +%Y%m%d_%H%M%S)"
```

## Preserve temporary test files

Temporary files created during verification are normally removed.

Use `--keep-temp` to preserve them inside the report directory:

```bash
bash MTD_check_installation.sh \
  --mode deep \
  --report-dir ./MTD_deep_check \
  --keep-temp
```

This option is mainly useful for debugging failed tests.

## Strict mode

By default, warnings do not necessarily produce a failing process status.

With `--strict`, warnings are treated as final failure and the checker returns
exit status `1`:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --strict
```

Strict mode is useful for:

- automated validation;
- continuous integration;
- installation benchmarking;
- detecting warnings in shell scripts;
- requiring a completely clean verification report.

The exit status can be inspected immediately after the checker finishes:

```bash
echo $?
```

## Recommended validation workflow

After installation, run the default full verification:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --report-dir ./MTD_check_full
```

For most users, `full` mode is the recommended post-installation check.

Use `quick` mode when you only need a fast structural check after editing
scripts, changing permissions, or moving an existing installation:

```bash
bash MTD_check_installation.sh \
  --mode quick \
  --report-dir ./MTD_check_quick
```

Use `deep` mode after clean installations, interrupted downloads, cache
transfers, or suspected cache corruption:

```bash
bash MTD_check_installation.sh \
  --mode deep \
  --report-dir ./MTD_check_deep
```

This avoids running the most expensive checks unnecessarily while still
providing a clear escalation path for debugging.

## Save the terminal output

The checker output can also be saved with `tee`:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --report-dir ./MTD_check_full \
  2>&1 | tee MTD_check_full.log
```

For a timestamped log:

```bash
bash MTD_check_installation.sh \
  --mode full \
  --report-dir "MTD_check_$(date +%Y%m%d_%H%M%S)" \
  2>&1 | tee "MTD_check_$(date +%Y%m%d_%H%M%S).log"
```

## Understanding the results

The checker reports conditions using statuses such as:

| Status | Meaning |
|---|---|
| `PASS` | The component was found and passed the corresponding check |
| `WARN` | The component may be usable, but the condition requires attention |
| `FAIL` | A required component is missing, invalid, incomplete, or unusable |

A warning should be reviewed before running a real analysis.

A failure should be resolved before using MTD Explorer.

!!! warning "Do not evaluate only the final line"

    Review the complete checker output and generated reports.

    An installation can contain several valid components while still
    having an incomplete environment, database, index, or reference file.

## Reporting installation problems

If the installation checker reports a `FAIL`, or if a `WARN` remains unclear
after reviewing this documentation, please open an issue in the MTD Explorer
GitHub repository:

[Open a GitHub issue](https://github.com/patrick-douglas/MTD-Explorer/issues)

When opening an issue, include the information listed below whenever possible.
This makes it easier to reproduce the problem and identify whether it is
related to the installer, Conda environments, reference databases, file paths,
or local system configuration.

## Information to preserve

When reporting an installation problem, preserve:

- the complete installer log;
- the checker terminal output;
- the checker report directory;
- temporary test files when relevant;
- the MTD Explorer Git commit;
- the current Git working-tree state;
- operating-system information;
- hardware information;
- the Conda installation path;
- the persistent cache path;
- the Bracken read length;
- the verification mode used.

Record the current commit with:

```bash
git log -1 --oneline
```

Record the repository state with:

```bash
git status
```

Record basic system information with:

```bash
cat /etc/os-release
uname -a
nproc
free -h
df -h
```

## Next step

After the required checks pass, continue to the
[Quick start guide](quick-start.md).
