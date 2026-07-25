#!/usr/bin/env python3
"""
Generate compact, auditable diagnostics for one MTD installation benchmark run.

Inputs are the files already produced by MTD_benchmark_install.sh. The analyzer
never changes the installation or cache. It writes human-readable and TSV
reports inside the benchmark run directory and appends stable failure metrics
to summary.tsv so failed runs remain mergeable.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

ANSI_RE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
URL_RE = re.compile(r"https?://[^\s<>'\"\])]+")
HTTP_CODE_RE = re.compile(
    r"(?:HTTP(?:/[0-9.]+)?(?:\s+status(?:\s+was)?)?[^0-9]{0,20})([45][0-9]{2})",
    re.IGNORECASE,
)

@dataclass(frozen=True)
class Rule:
    category: str
    severity: str
    priority: int
    rule_id: str
    expression: re.Pattern[str]
    likely_cause: str

@dataclass(frozen=True)
class Match:
    line_number: int
    category: str
    severity: str
    priority: int
    rule_id: str
    text: str
    likely_cause: str

RULES: tuple[Rule, ...] = (
    Rule(
        "network_http", "fatal", 10, "http_5xx",
        re.compile(r"\bHTTP(?:/[0-9.]+)?(?:\s+status(?:\s+was)?)?[^0-9]{0,20}5[0-9]{2}\b|"
                   r"\b5(?:02|03|04)\s+(?:Bad Gateway|Service Unavailable|Gateway Timeout)\b",
                   re.IGNORECASE),
        "A remote HTTP service or mirror returned a server-side 5xx failure.",
    ),
    Rule(
        "network_http", "fatal", 12, "empty_or_short_download",
        re.compile(r"downloaded length\s+0\s*!=\s*reported length|"
                   r"downloaded length .* != reported length|"
                   r"unexpected end of file|empty download",
                   re.IGNORECASE),
        "A download completed with no data or fewer bytes than advertised.",
    ),
    Rule(
        "network_dns", "fatal", 14, "dns_resolution",
        re.compile(r"Could not resolve host|Temporary failure in name resolution|"
                   r"Name or service not known|nodename nor servname provided",
                   re.IGNORECASE),
        "The host name could not be resolved.",
    ),
    Rule(
        "network_http", "fatal", 16, "url_open_failure",
        re.compile(r"cannot open URL|download of package .* failed|"
                   r"curl: \([0-9]+\)|wget: unable|Connection reset by peer|"
                   r"Remote end closed connection|TLS|SSL connect error",
                   re.IGNORECASE),
        "A network transfer or HTTPS session failed.",
    ),
    Rule(
        "network_ftp", "fatal", 18, "ftp_transfer",
        re.compile(r"Net::FTP|PASV:|Transfer aborted|Broken pipe|"
                   r"Unable to close datastream|FTP.*(?:Timeout|Connection closed)|"
                   r"rsync_from_ncbi\.pl: unable to download",
                   re.IGNORECASE),
        "The legacy FTP data channel failed or was interrupted.",
    ),
    Rule(
        "out_of_memory", "fatal", 20, "oom",
        re.compile(r"\bOut of memory\b|\boom-kill\b|\bKilled\b|"
                   r"Cannot allocate memory|std::bad_alloc|MemoryError",
                   re.IGNORECASE),
        "The process or operating system ran out of usable memory.",
    ),
    Rule(
        "disk_space", "fatal", 22, "disk_full",
        re.compile(r"No space left on device|Disk quota exceeded|"
                   r"not enough space|filesystem is full",
                   re.IGNORECASE),
        "The target filesystem ran out of available space or quota.",
    ),
    Rule(
        "permission", "fatal", 24, "permission",
        re.compile(r"Permission denied|Operation not permitted|read-only file system",
                   re.IGNORECASE),
        "The process lacked permission to access or modify a required path.",
    ),
    Rule(
        "integrity", "fatal", 26, "integrity",
        re.compile(r"checksum mismatch|MD5.*mismatch|SHA-?256.*mismatch|"
                   r"invalid gzip|not in gzip format|corrupt(?:ed|ion)?|"
                   r"tar: .*Error|unexpected EOF",
                   re.IGNORECASE),
        "A cached or downloaded archive failed an integrity check.",
    ),
    Rule(
        "conda_solver", "fatal", 30, "conda_solver",
        re.compile(r"UnsatisfiableError|PackagesNotFoundError|LibMambaUnsatisfiableError|"
                   r"failed to solve|Could not solve for environment specs",
                   re.IGNORECASE),
        "Conda could not resolve the requested environment.",
    ),
    Rule(
        "conda_network", "fatal", 32, "conda_network",
        re.compile(r"CondaHTTPError|CondaSSLError|ProxyError|"
                   r"repodata.*(?:failed|error)",
                   re.IGNORECASE),
        "Conda failed while contacting a package channel or downloading metadata.",
    ),
    Rule(
        "kraken_taxonomy", "fatal", 35, "kraken_mapping",
        re.compile(r"accessions remain unmapped|reference mapping validation failed|"
                   r"without taxonomy mapping|unmapped\.prebuild",
                   re.IGNORECASE),
        "Kraken2 reference accessions could not be mapped to TaxIDs.",
    ),
    Rule(
        "missing_dependency", "fatal", 40, "dependency_missing",
        re.compile(r"dependenc(?:y|ies) .* (?:is|are) not available|"
                   r"dependencies still missing|Still missing:|"
                   r"ModuleNotFoundError|No module named|command not found",
                   re.IGNORECASE),
        "One or more required software or package dependencies were unavailable.",
    ),
    Rule(
        "r_package_install", "fatal", 42, "r_install",
        re.compile(r"installation of one or more packages failed|"
                   r"installation of package .* had non-zero exit status|"
                   r"\* removing .*/R/library/|ERROR: dependency|"
                   r"Error: .*package",
                   re.IGNORECASE),
        "An R or Bioconductor package installation failed.",
    ),
    Rule(
        "compiler", "fatal", 45, "build_failure",
        re.compile(r"make(?:\[[0-9]+\])?: \*\*\*|compilation failed|"
                   r"fatal error:|ld: .*error|linker command failed",
                   re.IGNORECASE),
        "Compilation or linking failed.",
    ),
    Rule(
        "syntax", "fatal", 48, "syntax",
        re.compile(r"syntax error|unexpected token|IndentationError|SyntaxError",
                   re.IGNORECASE),
        "A script or generated file contained invalid syntax.",
    ),
    Rule(
        "required_script", "fatal", 55, "required_script",
        re.compile(r"\[ERROR\]\s+Required script failed:|"
                   r"\[ERROR\]\s+Required command failed:|"
                   r"\[ERROR\]\s+.* failed with exit status",
                   re.IGNORECASE),
        "A required installer subcommand or script returned a non-zero status.",
    ),
    Rule(
        "generic_failure", "fatal", 80, "generic_error",
        re.compile(r"^\s*(?:\[ERROR\]|\[FAIL\]|ERROR:|Error:|Execution halted)",
                   re.IGNORECASE),
        "The installer emitted a fatal error that did not match a more specific category.",
    ),
)

SUMMARY_METRICS = {
    "failure_analysis_status",
    "failure_class",
    "failure_primary_message",
    "failure_likely_cause",
    "failure_primary_line",
    "failure_failed_function",
    "failure_required_script",
    "failure_command",
    "failure_http_status_codes",
    "failure_urls",
    "failure_resource_flags",
    "failure_match_count",
    "failure_report_path",
    "failure_context_path",
}

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate diagnostics for one MTD benchmark run."
    )
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--exit-status", type=int)
    parser.add_argument("--context-before", type=int, default=12)
    parser.add_argument("--context-after", type=int, default=24)
    parser.add_argument("--max-context-events", type=int, default=12)
    parser.add_argument("--tail-lines", type=int, default=300)
    return parser.parse_args()

def clean_line(value: str) -> str:
    return ANSI_RE.sub("", value.rstrip("\r\n"))

def read_lines(path: Path) -> list[str]:
    if not path.is_file():
        return []
    return [
        clean_line(line)
        for line in path.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()
    ]

def read_metric_summary(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    if not path.is_file():
        return result
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            name = (row.get("metric") or "").strip()
            if name:
                result[name] = (row.get("value") or "").strip()
    return result

def numeric(value: str | None) -> float | None:
    try:
        return float(value) if value not in (None, "", "NA") else None
    except ValueError:
        return None

def infer_exit_status(
    explicit: int | None,
    summary: dict[str, str],
    run_dir: Path,
) -> int:
    if explicit is not None:
        return explicit
    parsed = numeric(summary.get("exit_status"))
    if parsed is not None:
        return int(parsed)
    metadata = run_dir / "metadata.txt"
    if metadata.is_file():
        match = re.search(
            r"^exit_status\t(-?[0-9]+)\s*$",
            metadata.read_text(encoding="utf-8", errors="replace"),
            flags=re.MULTILINE,
        )
        if match:
            return int(match.group(1))
    return 1

def classify(lines: list[str]) -> list[Match]:
    matches: list[Match] = []
    for number, line in enumerate(lines, 1):
        if not line.strip():
            continue
        for rule in RULES:
            if rule.expression.search(line):
                matches.append(
                    Match(
                        line_number=number,
                        category=rule.category,
                        severity=rule.severity,
                        priority=rule.priority,
                        rule_id=rule.rule_id,
                        text=line.strip(),
                        likely_cause=rule.likely_cause,
                    )
                )
    return matches

def select_primary(matches: list[Match]) -> Match | None:
    if not matches:
        return None
    return min(
        matches,
        key=lambda item: (
            item.priority,
            item.line_number,
            len(item.text),
        ),
    )

def unique_join(values: Iterable[str], limit: int = 10) -> str:
    output: list[str] = []
    seen: set[str] = set()
    for value in values:
        value = value.strip().rstrip(".,;:")
        if not value or value in seen:
            continue
        seen.add(value)
        output.append(value)
        if len(output) >= limit:
            break
    return ", ".join(output)

def extract_last(lines: list[str], expression: str) -> str:
    regex = re.compile(expression, re.IGNORECASE)
    result = ""
    for line in lines:
        match = regex.search(line)
        if match:
            result = match.group(1).strip()
    return result

def analyze_steps(path: Path) -> tuple[str, list[dict[str, str]], list[str]]:
    if not path.is_file():
        return "", [], []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    failed: list[dict[str, str]] = []
    for row in rows:
        value = (row.get("status") or "").strip()
        parsed = numeric(value)
        if parsed is not None and int(parsed) != 0:
            failed.append(row)
    selected = failed[-1] if failed else (rows[-1] if rows else {})
    function = (
        selected.get("function")
        or selected.get("step")
        or selected.get("name")
        or ""
    ).strip()
    return function, failed, fieldnames

def resource_flags(summary: dict[str, str]) -> list[str]:
    flags: list[str] = []
    minimum_available = numeric(summary.get("minimum_system_memory_available"))
    peak_swap = numeric(summary.get("peak_process_tree_swap"))
    peak_rss = numeric(summary.get("peak_process_tree_rss"))
    peak_temp = numeric(summary.get("peak_reported_temperature"))
    if minimum_available is not None:
        if minimum_available < 2:
            flags.append(
                f"critical_memory_pressure:min_available={minimum_available:.3f}GiB"
            )
        elif minimum_available < 8:
            flags.append(
                f"memory_pressure:min_available={minimum_available:.3f}GiB"
            )
    if peak_swap is not None and peak_swap > 0.5:
        flags.append(f"swap_used:peak={peak_swap:.3f}GiB")
    if peak_rss is not None and peak_rss > 100:
        flags.append(f"very_high_process_rss:peak={peak_rss:.3f}GiB")
    if peak_temp is not None and peak_temp >= 95:
        flags.append(f"high_temperature:peak={peak_temp:.1f}C")
    available_values = [
        numeric(value)
        for key, value in summary.items()
        if key.startswith("watch_after_available_bytes::")
    ]
    available_values = [value for value in available_values if value is not None]
    if available_values and min(available_values) < 20 * 1024**3:
        flags.append(
            f"low_filesystem_space:min_available={min(available_values)/1024**3:.3f}GiB"
        )
    return flags

def merge_intervals(
    points: list[int],
    total_lines: int,
    before: int,
    after: int,
) -> list[tuple[int, int]]:
    intervals = sorted(
        (
            max(1, point - before),
            min(total_lines, point + after),
        )
        for point in points
    )
    merged: list[tuple[int, int]] = []
    for start, end in intervals:
        if not merged or start > merged[-1][1] + 2:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged

def write_context(
    path: Path,
    lines: list[str],
    selected_matches: list[Match],
    before: int,
    after: int,
) -> None:
    points = [item.line_number for item in selected_matches]
    intervals = merge_intervals(points, len(lines), before, after)
    matched_lines = {item.line_number for item in selected_matches}
    with path.open("w", encoding="utf-8") as handle:
        if not intervals:
            handle.write("No classified failure context was found.\n")
            return
        for index, (start, end) in enumerate(intervals, 1):
            handle.write(
                f"===== FAILURE CONTEXT {index}: lines {start}-{end} =====\n"
            )
            for line_number in range(start, end + 1):
                marker = ">>" if line_number in matched_lines else "  "
                handle.write(
                    f"{marker} {line_number:9d} | {lines[line_number - 1]}\n"
                )
            handle.write("\n")

def write_matches(path: Path, matches: list[Match]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "line_number",
                "category",
                "severity",
                "priority",
                "rule_id",
                "message",
            ]
        )
        for item in matches:
            writer.writerow(
                [
                    item.line_number,
                    item.category,
                    item.severity,
                    item.priority,
                    item.rule_id,
                    item.text,
                ]
            )

def write_failed_steps(
    path: Path,
    failed_rows: list[dict[str, str]],
    fieldnames: list[str],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        columns = fieldnames or ["status"]
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in failed_rows:
            writer.writerow(row)

def update_summary(path: Path, values: dict[str, tuple[str, str]]) -> None:
    rows: list[dict[str, str]] = []
    if path.is_file():
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if (row.get("metric") or "") not in SUMMARY_METRICS:
                    rows.append(
                        {
                            "metric": row.get("metric", ""),
                            "value": row.get("value", ""),
                            "unit_or_definition": row.get(
                                "unit_or_definition", ""
                            ),
                        }
                    )
    for metric, (value, unit) in values.items():
        rows.append(
            {
                "metric": metric,
                "value": value,
                "unit_or_definition": unit,
            }
        )
    temp = path.with_name(f".{path.name}.failure-analysis.tmp")
    with temp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["metric", "value", "unit_or_definition"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)

def write_report(
    path: Path,
    *,
    run_dir: Path,
    exit_status: int,
    primary: Match | None,
    failure_class: str,
    likely_cause: str,
    failed_function: str,
    required_script: str,
    failed_command: str,
    http_codes: str,
    urls: str,
    flags: list[str],
    match_count: int,
    files: dict[str, Path],
) -> None:
    lines = [
        "MTD BENCHMARK FAILURE DIAGNOSTIC REPORT",
        "=" * 60,
        f"Run directory:        {run_dir}",
        f"Exit status:          {exit_status}",
        f"Classification:       {failure_class}",
        f"Primary log line:     {primary.line_number if primary else 'NA'}",
        f"Primary evidence:     {primary.text if primary else 'No classified fatal line found.'}",
        f"Likely cause:         {likely_cause}",
        f"Failed function:      {failed_function or 'NA'}",
        f"Required script:      {required_script or 'NA'}",
        f"Failed command:       {failed_command or 'NA'}",
        f"HTTP status codes:    {http_codes or 'NA'}",
        f"Related URLs:         {urls or 'NA'}",
        f"Resource flags:       {', '.join(flags) if flags else 'none detected'}",
        f"Diagnostic matches:   {match_count}",
        "",
        "Generated diagnostic files:",
    ]
    for label, file_path in files.items():
        lines.append(f"  {label:20s} {file_path}")
    lines.extend(
        [
            "",
            "Interpretation:",
            (
                "  This report is a heuristic diagnosis. The cited console "
                "lines remain the source of truth."
            ),
            (
                "  Inspect failure_context.log first, then console_tail.log "
                "and diagnostic_matches.tsv."
            ),
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")

def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).expanduser().resolve()
    if not run_dir.is_dir():
        print(f"[ERROR] Benchmark run directory not found: {run_dir}", file=sys.stderr)
        return 2

    console = run_dir / "console.log"
    summary_path = run_dir / "summary.tsv"
    steps_path = run_dir / "steps.tsv"
    lines = read_lines(console)
    summary = read_metric_summary(summary_path)
    exit_status = infer_exit_status(args.exit_status, summary, run_dir)

    matches = classify(lines)
    primary = select_primary(matches)
    failed_function, failed_rows, step_columns = analyze_steps(steps_path)
    flags = resource_flags(summary)

    required_script = extract_last(
        lines, r"\[ERROR\]\s+Required script failed:\s*(.+)$"
    )
    failed_command = extract_last(
        lines, r"\[ERROR\]\s+Command:\s*(.+)$"
    )
    http_codes = unique_join(
        match.group(1)
        for line in lines
        for match in HTTP_CODE_RE.finditer(line)
    )
    urls = unique_join(
        match.group(0)
        for item in matches
        for match in URL_RE.finditer(item.text)
    )

    if exit_status == 0:
        failure_class = "success"
        likely_cause = "The benchmark command completed successfully."
        primary = None
    elif primary:
        failure_class = primary.category
        likely_cause = primary.likely_cause
    else:
        failure_class = "unclassified_failure"
        likely_cause = (
            "The command failed, but no configured failure signature matched."
        )

    ranked = sorted(
        matches,
        key=lambda item: (item.priority, item.line_number),
    )
    selected_matches: list[Match] = []
    selected_lines: set[int] = set()
    for item in ranked:
        if item.line_number in selected_lines:
            continue
        selected_matches.append(item)
        selected_lines.add(item.line_number)
        if len(selected_matches) >= args.max_context_events:
            break

    report_path = run_dir / "failure_report.txt"
    report_tsv = run_dir / "failure_report.tsv"
    context_path = run_dir / "failure_context.log"
    tail_path = run_dir / "console_tail.log"
    matches_path = run_dir / "diagnostic_matches.tsv"
    failed_steps_path = run_dir / "failed_steps.tsv"

    write_context(
        context_path,
        lines,
        selected_matches,
        max(args.context_before, 0),
        max(args.context_after, 0),
    )
    tail_path.write_text(
        "\n".join(lines[-max(args.tail_lines, 1):]) + ("\n" if lines else ""),
        encoding="utf-8",
    )
    write_matches(matches_path, matches)
    write_failed_steps(failed_steps_path, failed_rows, step_columns)

    fields = {
        "analysis_status": "generated",
        "exit_status": str(exit_status),
        "failure_class": failure_class,
        "primary_line": str(primary.line_number if primary else ""),
        "primary_message": primary.text if primary else "",
        "likely_cause": likely_cause,
        "failed_function": failed_function,
        "required_script": required_script,
        "failed_command": failed_command,
        "http_status_codes": http_codes,
        "urls": urls,
        "resource_flags": ",".join(flags),
        "diagnostic_match_count": str(len(matches)),
    }
    with report_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["field", "value"])
        writer.writerows(fields.items())

    generated_files = {
        "human report": report_path,
        "machine summary": report_tsv,
        "error context": context_path,
        "console tail": tail_path,
        "all matches": matches_path,
        "failed steps": failed_steps_path,
    }
    write_report(
        report_path,
        run_dir=run_dir,
        exit_status=exit_status,
        primary=primary,
        failure_class=failure_class,
        likely_cause=likely_cause,
        failed_function=failed_function,
        required_script=required_script,
        failed_command=failed_command,
        http_codes=http_codes,
        urls=urls,
        flags=flags,
        match_count=len(matches),
        files=generated_files,
    )

    update_summary(
        summary_path,
        {
            "failure_analysis_status": ("generated", "text"),
            "failure_class": (failure_class, "text"),
            "failure_primary_message": (
                primary.text if primary else "",
                "text",
            ),
            "failure_likely_cause": (likely_cause, "text"),
            "failure_primary_line": (
                str(primary.line_number if primary else "NA"),
                "console.log line number",
            ),
            "failure_failed_function": (failed_function or "NA", "text"),
            "failure_required_script": (required_script or "NA", "text"),
            "failure_command": (failed_command or "NA", "text"),
            "failure_http_status_codes": (http_codes or "NA", "text"),
            "failure_urls": (urls or "NA", "text"),
            "failure_resource_flags": (
                ",".join(flags) if flags else "none",
                "text",
            ),
            "failure_match_count": (str(len(matches)), "count"),
            "failure_report_path": (str(report_path), "path"),
            "failure_context_path": (str(context_path), "path"),
        },
    )

    print("[DIAGNOSTIC] Failure analysis generated.")
    print(f"[DIAGNOSTIC] Classification: {failure_class}")
    if primary:
        print(
            f"[DIAGNOSTIC] Primary evidence: "
            f"console.log:{primary.line_number}: {primary.text}"
        )
    print(f"[DIAGNOSTIC] Failed function: {failed_function or 'NA'}")
    print(f"[DIAGNOSTIC] Resource flags: {','.join(flags) if flags else 'none'}")
    print(f"[DIAGNOSTIC] Report: {report_path}")
    print(f"[DIAGNOSTIC] Context: {context_path}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
