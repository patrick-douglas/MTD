#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from collections import defaultdict, deque
from typing import Iterable

ANSI_RE = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
FATAL_PATTERNS = [
    ("error_marker", re.compile(r"\[ERROR\]", re.IGNORECASE)),
    ("traceback", re.compile(r"Traceback \(most recent call last\)", re.IGNORECASE)),
    ("no_space", re.compile(r"No space left on device", re.IGNORECASE)),
    ("killed", re.compile(r"(^|\s)Killed(\s|$)", re.IGNORECASE)),
    ("segfault", re.compile(r"Segmentation fault|core dumped", re.IGNORECASE)),
    ("command_not_found", re.compile(r"command not found", re.IGNORECASE)),
    ("memory", re.compile(r"cannot allocate memory|out of memory", re.IGNORECASE)),
    ("explicit_failure", re.compile(r"failed with exit status", re.IGNORECASE)),
]
COMPLETION_MARKERS = (
    "MTD Explorer analysis is finished.",
    "MTD Explorer exploratory run is finished.",
)
HASH_MAX_BYTES = 20 * 1024 * 1024


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate MTD Explorer pipeline-specific benchmark reports."
    )
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--pipeline-output", required=True, type=Path)
    parser.add_argument("--samplesheet", required=True, type=Path)
    parser.add_argument("--repo-dir", required=True, type=Path)
    parser.add_argument("--dataset-label", required=True)
    parser.add_argument("--machine", required=True)
    parser.add_argument("--execution-mode", required=True, choices=("local", "hpc"))
    parser.add_argument("--hpc-conf", required=True)
    parser.add_argument("--hostid", required=True)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--top-n", required=True, type=int)
    parser.add_argument("--read-layout", required=True)
    parser.add_argument("--analysis-mode", required=True)
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--extract-reads", required=True, type=int, choices=(0, 1))
    parser.add_argument("--exit-status", required=True, type=int)
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def command_output(args: list[str], cwd: Path | None = None) -> str:
    try:
        completed = subprocess.run(
            args,
            cwd=cwd,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    except OSError:
        return "NA"
    return completed.stdout.strip() or "NA"


def strip_ansi(text: str) -> str:
    return ANSI_RE.sub("", text)


def iter_compact_console_lines(path: Path, chunk_size: int = 8 * 1024 * 1024):
    """Yield terminal-style logical lines without retaining CR progress history.

    Many command-line tools redraw a progress line using carriage returns.
    Treat each carriage return as a cursor reset and retain only the most
    recent text before the next newline. CRLF is preserved as a normal newline.
    The implementation is streaming and bounded by the current displayed line.
    """
    current = bytearray()
    pending_cr = False

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break

            if pending_cr:
                if chunk.startswith(b"\n"):
                    yield bytes(current)
                    current.clear()
                    chunk = chunk[1:]
                else:
                    current.clear()
                pending_cr = False

            if not chunk:
                continue

            if chunk.endswith(b"\r"):
                chunk = chunk[:-1]
                pending_cr = True

            parts = chunk.split(b"\n")
            last_index = len(parts) - 1

            for index, part in enumerate(parts):
                complete = index < last_index

                if complete and part.endswith(b"\r"):
                    # CRLF: newline, not a progress-line redraw.
                    part = part[:-1]

                if b"\r" in part:
                    # Only the text after the final carriage return remains
                    # visible on a terminal before the newline.
                    current = bytearray(part.rsplit(b"\r", 1)[-1])
                else:
                    current.extend(part)

                if complete:
                    yield bytes(current)
                    current.clear()

    if current:
        yield bytes(current)


def read_summary(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = row.get("metric", "")
            if key:
                values[key] = row.get("value", "")
    return values


def read_metric_rows(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    if not path.is_file():
        return rows
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metric = row.get("metric", "")
            if metric:
                rows.append(
                    (
                        metric,
                        row.get("value", ""),
                        row.get("unit_or_definition", ""),
                    )
                )
    return rows


def read_samplesheet(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel
        reader = csv.reader(handle, dialect)
        rows = [row for row in reader if any(cell.strip() for cell in row)]
    if not rows:
        return [], []
    return rows[0], rows[1:]


def discover_input_files(samplesheet: Path, rows: list[list[str]]) -> list[Path]:
    suffixes = (
        ".fastq",
        ".fq",
        ".fastq.gz",
        ".fq.gz",
        ".fasta",
        ".fa",
        ".fna",
        ".sra",
    )
    found: dict[str, Path] = {}
    base = samplesheet.parent

    for row in rows:
        for raw_cell in row:
            cell = raw_cell.strip().strip('"').strip("'")
            if not cell:
                continue
            lower = cell.lower()
            if not lower.endswith(suffixes):
                continue
            candidate = Path(os.path.expanduser(cell))
            if not candidate.is_absolute():
                candidate = base / candidate
            try:
                resolved = candidate.resolve(strict=False)
            except OSError:
                resolved = candidate.absolute()
            found[str(resolved)] = resolved

    return [found[key] for key in sorted(found)]


def write_input_files(path: Path, files: Iterable[Path]) -> tuple[int, int, int]:
    total = 0
    existing_count = 0
    missing_count = 0
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["path", "exists", "size_bytes", "mtime_utc"])
        for item in files:
            exists = item.is_file()
            size = item.stat().st_size if exists else 0
            if exists:
                existing_count += 1
                total += size
                mtime = dt.datetime.fromtimestamp(
                    item.stat().st_mtime, tz=dt.timezone.utc
                ).isoformat()
            else:
                missing_count += 1
                mtime = "NA"
            writer.writerow([str(item), int(exists), size, mtime])
    return existing_count, missing_count, total


def inventory_output(
    output: Path, destination: Path
) -> tuple[int, int, int, dict[str, list[int]]]:
    file_count = 0
    dir_count = 0
    total_bytes = 0
    ext_summary: dict[str, list[int]] = defaultdict(lambda: [0, 0])

    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "relative_path",
                "type",
                "size_bytes",
                "mtime_utc",
                "sha256_if_le_20MiB",
            ]
        )

        if not output.exists():
            return 0, 0, 0, ext_summary

        for item in sorted(output.rglob("*"), key=lambda p: str(p)):
            try:
                stat_result = item.lstat()
            except OSError:
                continue
            relative = item.relative_to(output)

            if item.is_symlink():
                writer.writerow([str(relative), "symlink", 0, "NA", "NA"])
                continue

            if item.is_dir():
                dir_count += 1
                writer.writerow([str(relative), "directory", 0, "NA", "NA"])
                continue

            if not item.is_file():
                writer.writerow([str(relative), "other", 0, "NA", "NA"])
                continue

            file_count += 1
            size = stat_result.st_size
            total_bytes += size
            mtime = dt.datetime.fromtimestamp(
                stat_result.st_mtime, tz=dt.timezone.utc
            ).isoformat()

            suffixes = "".join(item.suffixes).lower()
            extension = suffixes if suffixes else "[no_extension]"
            ext_summary[extension][0] += 1
            ext_summary[extension][1] += size

            checksum = "SKIPPED_GT_20MiB"
            if size <= HASH_MAX_BYTES:
                try:
                    checksum = sha256_file(item)
                except OSError:
                    checksum = "ERROR"

            writer.writerow([str(relative), "file", size, mtime, checksum])

    return file_count, dir_count, total_bytes, ext_summary


def write_extension_summary(path: Path, summary: dict[str, list[int]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["extension", "file_count", "total_bytes"])
        for extension, (count, size) in sorted(
            summary.items(), key=lambda item: (-item[1][1], item[0])
        ):
            writer.writerow([extension, count, size])


def process_console(run_dir: Path) -> tuple[int, bool]:
    console = run_dir / "console.log"
    clean_console = run_dir / "console_clean.log"
    tail_file = run_dir / "final_console_tail.txt"
    hits_file = run_dir / "diagnostic_hits.tsv"

    hits: list[tuple[int, str, str]] = []
    tail: deque[str] = deque(maxlen=250)
    completion = False

    with clean_console.open("w", encoding="utf-8") as clean_handle:
        if console.is_file():
            for line_number, raw_line in enumerate(
                iter_compact_console_lines(console),
                start=1,
            ):
                line = strip_ansi(
                    raw_line.decode("utf-8", errors="replace")
                )
                clean_handle.write(line + "\n")
                tail.append(line)

                if any(marker in line for marker in COMPLETION_MARKERS):
                    completion = True

                for label, pattern in FATAL_PATTERNS:
                    if pattern.search(line):
                        hits.append((line_number, label, line))
                        break

    tail_file.write_text(
        "\n".join(tail) + ("\n" if tail else ""),
        encoding="utf-8",
    )

    with hits_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["line_number", "pattern", "text"])
        writer.writerows(hits)

    return len(hits), completion


def parse_steps(raw_path: Path, output_path: Path) -> tuple[int, float]:
    events: list[dict[str, str]] = []
    if raw_path.is_file():
        with raw_path.open(newline="", encoding="utf-8", errors="replace") as handle:
            events = list(csv.DictReader(handle, delimiter="\t"))

    rows: list[list[object]] = []
    total = 0.0
    for index in range(max(len(events) - 1, 0)):
        current = events[index]
        following = events[index + 1]
        if current.get("message") == "__BENCHMARK_WRAPPER_END__":
            continue
        try:
            start_ns = int(current["epoch_ns"])
            end_ns = int(following["epoch_ns"])
        except (KeyError, TypeError, ValueError):
            continue
        duration = max((end_ns - start_ns) / 1_000_000_000, 0.0)
        total += duration
        rows.append(
            [
                index + 1,
                current.get("timestamp_utc", "NA"),
                following.get("timestamp_utc", "NA"),
                current.get("percent", "NA"),
                current.get("message", ""),
                f"{duration:.6f}",
            ]
        )

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "stage_index",
                "start_utc",
                "end_utc",
                "start_percent",
                "stage",
                "duration_seconds",
            ]
        )
        writer.writerows(rows)

    return len(rows), total


def write_failure_report(
    path: Path,
    exit_status: int,
    completion: bool,
    output_exists: bool,
    output_file_count: int,
    diagnostic_hits: int,
    clean_console: Path,
) -> None:
    if exit_status == 0 and completion and output_exists and output_file_count > 0:
        if path.exists():
            path.unlink()
        return

    lines: list[str] = [
        "MTD Explorer pipeline benchmark failure report",
        "================================================",
        f"Exit status: {exit_status}",
        f"Completion marker found: {int(completion)}",
        f"Output directory exists: {int(output_exists)}",
        f"Output file count: {output_file_count}",
        f"Fatal-pattern hits in console: {diagnostic_hits}",
        "",
        "Last relevant console context:",
        "------------------------------",
    ]

    if clean_console.is_file():
        console_lines = clean_console.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()
        selected: list[str] = []
        hit_indices: list[int] = []
        for index, line in enumerate(console_lines):
            if any(pattern.search(line) for _, pattern in FATAL_PATTERNS):
                hit_indices.append(index)
        for index in hit_indices[-8:]:
            start = max(index - 4, 0)
            end = min(index + 5, len(console_lines))
            selected.extend(console_lines[start:end])
            selected.append("-----")
        if not selected:
            selected = console_lines[-120:]
        lines.extend(selected)
    else:
        lines.append("console_clean.log was not available.")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    output = args.pipeline_output.resolve(strict=False)
    samplesheet = args.samplesheet.resolve()
    repo = args.repo_dir.resolve()

    run_dir.mkdir(parents=True, exist_ok=True)

    copied_samplesheet = run_dir / "input_samplesheet.csv"
    shutil.copy2(samplesheet, copied_samplesheet)

    header, data_rows = read_samplesheet(samplesheet)
    input_files = discover_input_files(samplesheet, data_rows)
    existing_inputs, missing_inputs, input_total_bytes = write_input_files(
        run_dir / "input_files.tsv", input_files
    )

    output_file_count, output_dir_count, output_total_bytes, ext_summary = (
        inventory_output(output, run_dir / "output_inventory.tsv")
    )
    write_extension_summary(run_dir / "output_extensions.tsv", ext_summary)

    diagnostic_hits, completion = process_console(run_dir)
    stage_count, stage_total_seconds = parse_steps(
        run_dir / "pipeline_steps_raw.tsv",
        run_dir / "pipeline_steps.tsv",
    )

    generic_summary = read_summary(run_dir / "summary.tsv")
    hpc_metric_rows = (
        read_metric_rows(run_dir / "hpc_summary.tsv")
        if args.execution_mode == "hpc"
        else []
    )
    hpc_summary = {metric: value for metric, value, _unit in hpc_metric_rows}
    git_commit = command_output(["git", "rev-parse", "HEAD"], cwd=repo)
    git_branch = command_output(["git", "branch", "--show-current"], cwd=repo)
    git_status = command_output(["git", "status", "--short"], cwd=repo)
    git_dirty = 0 if git_status == "NA" or not git_status.strip() else 1

    mtd_script = repo / "MTD_explorer.sh"
    script_sha = sha256_file(mtd_script) if mtd_script.is_file() else "NA"
    samplesheet_sha = sha256_file(samplesheet)

    output_exists = output.is_dir()
    overall_pass = (
        args.exit_status == 0
        and completion
        and output_exists
        and output_file_count > 0
    )
    hpc_collection_status = (
        hpc_summary.get("hpc_collection_status", "MISSING")
        if args.execution_mode == "hpc"
        else "NA"
    )
    benchmark_data_status = (
        "FAIL"
        if not overall_pass
        else (
            "PASS"
            if args.execution_mode == "local" or hpc_collection_status == "PASS"
            else "PASS_WITH_PARTIAL_HPC_METRICS"
        )
    )

    extraction_mode = (
        "disabled"
        if not args.extract_reads
        else ("all_detected_taxa" if args.top_n == 0 else f"top_{args.top_n}")
    )

    metrics: list[tuple[str, object, str]] = [
        ("pipeline_benchmark_status", "PASS" if overall_pass else "FAIL", "text"),
        ("benchmark_data_status", benchmark_data_status, "text"),
        ("execution_mode", args.execution_mode, "local or hpc"),
        ("hpc_collection_status", hpc_collection_status, "PASS, PARTIAL, NO_JOBS, MISSING, or NA"),
        ("machine", args.machine, "controller machine or cluster label"),
        ("controller_threads", args.threads, "count"),
        ("hpc_configuration", args.hpc_conf if args.execution_mode == "hpc" else "NA", "path"),
        ("dataset_label", args.dataset_label, "text"),
        ("samplesheet", str(samplesheet), "path"),
        ("samplesheet_sha256", samplesheet_sha, "SHA-256"),
        ("samplesheet_columns", ",".join(header), "text"),
        ("samplesheet_data_rows", len(data_rows), "count"),
        ("referenced_input_files", len(input_files), "count"),
        ("existing_referenced_input_files", existing_inputs, "count"),
        ("missing_referenced_input_files", missing_inputs, "count"),
        ("referenced_input_total_bytes", input_total_bytes, "bytes"),
        ("host_taxid", args.hostid, "NCBI TaxID"),
        ("threads", args.threads, "backward-compatible alias for controller_threads"),
        ("requested_read_layout", args.read_layout, "text"),
        ("analysis_mode", args.analysis_mode, "text"),
        ("alignment", args.alignment, "text"),
        ("microbiome_read_extraction", extraction_mode, "text"),
        ("extract_top_n", args.top_n, "0 means all"),
        ("pipeline_exit_status", args.exit_status, "integer"),
        ("completion_marker_found", int(completion), "0/1"),
        ("diagnostic_hit_count", diagnostic_hits, "count"),
        ("pipeline_output", str(output), "path"),
        ("pipeline_output_exists", int(output_exists), "0/1"),
        ("pipeline_output_file_count", output_file_count, "count"),
        ("pipeline_output_directory_count", output_dir_count, "count"),
        ("pipeline_output_total_bytes", output_total_bytes, "bytes"),
        ("pipeline_stage_count", stage_count, "count"),
        ("pipeline_stage_total_seconds", f"{stage_total_seconds:.6f}", "seconds"),
        (
            "generic_wall_time_seconds",
            generic_summary.get("wall_time_seconds", "NA"),
            "seconds",
        ),
        (
            "mean_system_cpu_busy",
            generic_summary.get("mean_system_cpu_busy", "NA"),
            "percent",
        ),
        (
            "peak_process_tree_cpu_one_core",
            generic_summary.get("peak_process_tree_cpu_one_core", "NA"),
            "percent of one core",
        ),
        (
            "peak_process_tree_cpu_total_capacity",
            generic_summary.get("peak_process_tree_cpu_total_capacity", "NA"),
            "percent of total CPU capacity",
        ),
        (
            "peak_process_tree_rss",
            generic_summary.get("peak_process_tree_rss", "NA"),
            "GiB",
        ),
        (
            "peak_system_memory_used",
            generic_summary.get("peak_system_memory_used", "NA"),
            "GiB",
        ),
        (
            "total_physical_disk_read",
            generic_summary.get("total_physical_disk_read", "NA"),
            "GiB",
        ),
        (
            "total_physical_disk_write",
            generic_summary.get("total_physical_disk_write", "NA"),
            "GiB",
        ),
        (
            "total_network_received",
            generic_summary.get("total_network_received", "NA"),
            "GiB",
        ),
        (
            "total_network_transmitted",
            generic_summary.get("total_network_transmitted", "NA"),
            "GiB",
        ),
        (
            "peak_reported_temperature",
            generic_summary.get("peak_reported_temperature", "NA"),
            "degrees Celsius",
        ),
        ("git_commit", git_commit, "text"),
        ("git_branch", git_branch, "text"),
        ("git_dirty", git_dirty, "0/1"),
        ("mtd_explorer_sha256", script_sha, "SHA-256"),
    ]

    existing_metric_names = {metric for metric, _value, _unit in metrics}
    metrics.extend(
        row for row in hpc_metric_rows if row[0] not in existing_metric_names
    )

    with (run_dir / "pipeline_summary.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value", "unit_or_definition"])
        writer.writerows(metrics)

    with (run_dir / "benchmark_comparison_row.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "dataset",
                "execution_mode",
                "machine_or_cluster",
                "host_taxid",
                "controller_threads",
                "slurm_nodes_used",
                "slurm_node_names",
                "end_to_end_wall_seconds",
                "slurm_compute_makespan_seconds",
                "status",
                "git_commit",
            ]
        )
        writer.writerow(
            [
                args.dataset_label,
                args.execution_mode,
                args.machine,
                args.hostid,
                args.threads,
                hpc_summary.get("slurm_nodes_used", "1" if args.execution_mode == "local" else "NA"),
                hpc_summary.get("slurm_node_names", args.machine if args.execution_mode == "local" else "NA"),
                generic_summary.get("wall_time_seconds", "NA"),
                hpc_summary.get("slurm_compute_makespan_seconds", "NA"),
                benchmark_data_status,
                git_commit,
            ]
        )

    summary_lines = [
        "MTD Explorer pipeline benchmark summary",
        "======================================",
        f"Status: {benchmark_data_status}",
        f"Pipeline status: {'PASS' if overall_pass else 'FAIL'}",
        f"Execution mode: {args.execution_mode}",
        f"Machine/cluster: {args.machine}",
        f"Dataset: {args.dataset_label}",
        f"Host TaxID: {args.hostid}",
        f"Controller threads: {args.threads}",
        f"Alignment: {args.alignment}",
        f"Extraction: {extraction_mode}",
        f"Exit status: {args.exit_status}",
        f"Completion marker: {'yes' if completion else 'no'}",
        f"Wall time: {generic_summary.get('wall_time_seconds', 'NA')} seconds",
        (
            "Peak process-tree RSS: "
            f"{generic_summary.get('peak_process_tree_rss', 'NA')} GiB"
        ),
        (
            "Mean system CPU busy: "
            f"{generic_summary.get('mean_system_cpu_busy', 'NA')}%"
        ),
        *(
            [
                f"HPC collection status: {hpc_collection_status}",
                f"Slurm nodes used: {hpc_summary.get('slurm_nodes_used', 'NA')}",
                f"Slurm node names: {hpc_summary.get('slurm_node_names', 'NA')}",
                f"Slurm jobs: {hpc_summary.get('slurm_job_count', 'NA')}",
                f"Slurm task attempts: {hpc_summary.get('slurm_array_task_attempt_count', 'NA')}",
                f"Slurm compute makespan: {hpc_summary.get('slurm_compute_makespan_seconds', 'NA')} seconds",
            ]
            if args.execution_mode == "hpc"
            else []
        ),
        f"Output files: {output_file_count}",
        f"Output size: {output_total_bytes} bytes",
        f"Recorded stages: {stage_count}",
        f"Diagnostic hits: {diagnostic_hits}",
        f"Git commit: {git_commit}",
        "",
        "Important files:",
        f"  {run_dir / 'summary.tsv'}",
        f"  {run_dir / 'pipeline_summary.tsv'}",
        f"  {run_dir / 'benchmark_comparison_row.tsv'}",
        f"  {run_dir / 'pipeline_steps.tsv'}",
        f"  {run_dir / 'resource_samples.csv'}",
        f"  {run_dir / 'output_inventory.tsv'}",
        f"  {run_dir / 'console_clean.log'}",
    ]
    if args.execution_mode == "hpc":
        summary_lines.extend(
            [
                f"  {run_dir / 'hpc_summary.tsv'}",
                f"  {run_dir / 'hpc_stage_summary.tsv'}",
                f"  {run_dir / 'hpc_node_summary.tsv'}",
                f"  {run_dir / 'hpc_jobs.tsv'}",
                f"  {run_dir / 'hpc_tasks.tsv'}",
            ]
        )
    (run_dir / "pipeline_summary.txt").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )

    write_failure_report(
        run_dir / "failure_report.txt",
        args.exit_status,
        completion,
        output_exists,
        output_file_count,
        diagnostic_hits,
        run_dir / "console_clean.log",
    )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(
            f"[ERROR] Could not generate pipeline benchmark report: {exc}",
            file=sys.stderr,
        )
        raise
