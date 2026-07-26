# MTD Explorer installation benchmark suite

This directory contains the complete benchmark and failure-diagnostic suite used
to measure clean MTD Explorer installations.

## Components

1. `run_mtd_clean_benchmark.sh`
   Clones the current repository, validates every benchmark component, generates
   and validates `Install_profiled.sh`, and only then asks for destructive cleanup
   confirmation. The persistent installation cache is preserved.

2. `MTD_make_instrumented_installer.sh`
   Creates `Install_profiled.sh` without modifying `Install.sh`. It supports both
   the current software-before-databases installer layout and historical snapshots
   containing `prepare_installation_cache()`.

3. `MTD_fix_profiled_locale.sh`
   Repairs older generated profilers that parse a locale-formatted
   `EPOCHREALTIME` value with a decimal comma.

4. `MTD_benchmark_install.sh`
   Runs the installer in the foreground while collecting wall time, CPU, RAM,
   process-tree RSS and swap, disk and network I/O, temperatures, hardware,
   software versions, Git state, logs, and watched-path sizes.

5. `MTD_benchmark_failure_report.py`
   Produces a compact failure diagnosis, contextual log extracts, and stable
   failure metrics that remain mergeable with successful runs.

6. `MTD_benchmark_merge.py`
   Combines benchmark runs from one or more machines into article-ready TSV
   tables and function-level timing summaries.

7. `README_MTD_BENCHMARK.md`
   This guide.

## Safety change after installer reorganization

The current installer separates these stages:

```text
prepare_installation_cache_structure
validate_all_software_before_databases
download_database_caches
```

The benchmark generator now recognizes that layout. More importantly, the clean
runner creates and validates the instrumented installer **before** displaying the
`Type DELETE` prompt. If instrumentation becomes incompatible with a future
`Install.sh`, the runner exits while the existing MTD and Conda installation are
still untouched.

## Recommended clean benchmark command

Keep the runner outside the repository because the clean workflow removes and
replaces `~/MTD-Explorer`:

```bash
cp ~/MTD-Explorer/benchmark/run_mtd_clean_benchmark.sh ~/
chmod +x ~/run_mtd_clean_benchmark.sh

bash ~/run_mtd_clean_benchmark.sh \
    -n mauro \
    -o /media/mago/MTD_Install_cach/MTD_install_cache
```

The runner detects whether the persistent cache is cold or warm and generates a
label ending in `_r2`.

## Manual profiled-installer generation

```bash
cd ~/MTD-Explorer

bash benchmark/MTD_make_instrumented_installer.sh \
    --input Install.sh \
    --output Install_profiled.sh

bash benchmark/MTD_fix_profiled_locale.sh Install_profiled.sh
bash -n Install_profiled.sh
```

For the reorganized installer, the generator should report:

```text
[OK] Installer layout detected: software-before-databases-v2
[OK] Profiler insertion point: main "$@"
```

## Manual benchmark wrapper

```bash
cd ~/MTD-Explorer

bash benchmark/MTD_benchmark_install.sh \
    --label master_warm_native_r2 \
    --interval 5 \
    --output-root "$HOME/MTD_benchmarks" \
    --watch-path "$HOME/miniconda3" \
    --watch-path "/path/to/MTD_install_cache" \
    --watch-path "$HOME/MTD-Explorer" \
    -- \
    bash ./Install_profiled.sh -o "/path/to/MTD_install_cache"
```

Do not pass passwords in command-line arguments. The benchmark records the
executed command in metadata; allow `sudo` to prompt normally.

## Output from each run

```text
console.log
git_state.txt
gnu_time.txt
hardware.txt
metadata.txt
resource_samples.csv
software.txt
steps.tsv
summary.tsv
watch_paths_before.tsv
watch_paths_after.tsv
failure_report.txt
failure_report.tsv
failure_context.log
console_tail.log
diagnostic_matches.tsv
failed_steps.tsv
```

Failure files are generated after both successful and failed clean-runner
executions. For a failure, inspect `failure_report.txt` and
`failure_context.log` first.

Nested function durations overlap. Use `parent_function`, `call_depth`, and the
major installation-stage functions rather than summing every row blindly.

## Merge benchmark runs

```bash
python3 benchmark/MTD_benchmark_merge.py \
    --input "$HOME/MTD_benchmarks_master" \
    --input "$HOME/MTD_benchmarks_secondary" \
    --output "$HOME/MTD_benchmark_merged"
```

The merged directory contains:

```text
runs_wide.tsv
runs_article.tsv
steps_long.tsv
steps_summary.tsv
```

## Path rule

Do not quote a path beginning with a literal tilde:

```bash
# Wrong
-o "~/MTD_install_cache"

# Correct
-o "$HOME/MTD_install_cache"

# Also correct
-o ~/MTD_install_cache
```

## Fair-run checklist

Before each timed run:

- use the same Git commit;
- record whether the persistent cache is cold or warm;
- reboot the machine when practical;
- stop unrelated analyses, backups, and cloud synchronization;
- use the same network type;
- keep the machine connected to AC power;
- preserve successful and failed run directories;
- run `MTD_check_installation.sh --deep` after a successful installation.
