# MTD Explorer pipeline benchmark

This benchmark measures a complete `MTD_explorer.sh` analysis while preserving
the scientific output separately from the compact benchmark bundle.

It reuses `MTD_benchmark_install.sh` as the low-level machine/resource monitor.
Despite its historical filename, that wrapper accepts any command after `--`.
The pipeline runner adds MTD-specific validation, stage timings, output
inventory, diagnostics and a compressed bundle.

## Files

- `run_mtd_pipeline_benchmark.sh`
  - validates the command and paths;
  - removes only the previous analysis output after confirmation;
  - fixes the thread count explicitly;
  - runs the complete MTD Explorer analysis;
  - records hardware and resource use;
  - creates the final benchmark bundle.
- `MTD_pipeline_benchmark_report.py`
  - creates the pipeline summary;
  - calculates durations between MTD progress stages;
  - inventories output files;
  - strips terminal colors from the log;
  - produces focused diagnostics after a failed run.
- `MTD_explorer.sh`
  - receives a tiny opt-in hook inside `show_progress`;
  - writes progress timestamps only when
    `MTD_PIPELINE_BENCHMARK_STEPS_TSV` is set;
  - normal runs are unchanged.

## Biomphalaria glabrata benchmark

The `--extract-microbiome-reads-top-n` option requires an integer.
Use `0` to extract every detected taxon.

```bash
cd ~/MTD-Explorer

bash benchmark/run_mtd_pipeline_benchmark.sh \
  --machine master \
  --dataset Bglabrata_PRJNA1306560 \
  --input /media/me/18TB_BACKUP_LBN/drive.ifpa/LBN_RNA-Seq/Metatranscriptomics/MTD/MTD_Offline_Install_files/personal/BioProject/PRJNA1306560/B.glabrata_fastq/samplesheet.csv \
  --output "$HOME/test_MTD_explorer_B.glabrata" \
  --hostid 6526 \
  --threads 20 \
  --top-n 0 \
  --run-number 1
```

The runner uses Magic-BLAST and microbiome-read extraction by default, matching
the intended command. Use `--hisat2` or `--no-extract` only for a different
benchmark design.

## Repetitions

Use a new repetition number while keeping the dataset, Git commit, options and
thread count identical:

```bash
--run-number 1
--run-number 2
--run-number 3
```

For a fair machine comparison:

- use the same Git commit;
- use identical FASTQ files and samplesheet;
- use the same fixed thread count when comparing per-thread performance;
- alternatively use each machine's full thread count, but label that as native;
- stop unrelated analyses, backups and cloud synchronization;
- preserve failed runs as well as successful runs;
- do not include the large scientific output directory in the bundle.

## Main outputs

Each benchmark run directory contains the generic monitor files plus:

```text
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
failure_report.txt        # only when the run is incomplete or failed
```

The bundle is created under:

```text
~/MTD_pipeline_benchmarks/
```

The large MTD analysis output remains at the path passed with `--output`.
