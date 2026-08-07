#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Iterable

SACCT_CANDIDATES = [
    "JobIDRaw",
    "JobName",
    "State",
    "ExitCode",
    "Submit",
    "Eligible",
    "Start",
    "End",
    "ElapsedRaw",
    "NodeList",
    "NCPUS",
    "CPUTimeRAW",
    "TotalCPU",
    "MaxRSS",
    "MaxVMSize",
    "AveRSS",
    "ReqMem",
    "AllocTRES",
]

JOB_ID_RE = re.compile(r"^[0-9]+$")
ARRAY_TASK_RE = re.compile(r"^(?P<job>[0-9]+)_(?P<task>[0-9]+)$")
ARRAY_STEP_RE = re.compile(r"^(?P<base>[0-9]+_[0-9]+)\.(?P<step>.+)$")
SCONTROL_FIELD_RE = re.compile(r"(?:^|\s)([A-Za-z][A-Za-z0-9]*)=(.*?)(?=\s+[A-Za-z][A-Za-z0-9]*=|$)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collect Slurm accounting and node metadata for an MTD Explorer HPC benchmark."
    )
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--pipeline-output", required=True, type=Path)
    parser.add_argument("--hpc-conf", required=True, type=Path)
    parser.add_argument("--controller-host", required=True)
    parser.add_argument("--controller-threads", required=True, type=int)
    return parser.parse_args()


def run_command(args: list[str]) -> tuple[int, str]:
    try:
        completed = subprocess.run(
            args,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    except OSError as exc:
        return 127, str(exc)
    return completed.returncode, completed.stdout.strip()


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def read_key_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in raw_line:
            continue
        key, value = raw_line.split("=", 1)
        key = key.strip()
        if key:
            values[key] = value.strip()
    return values


def read_manifest(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not path.is_file():
        return rows
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if not fields or not fields[0]:
                continue
            rows.append(
                {
                    "task_id": fields[0],
                    "task_hash": fields[1] if len(fields) > 1 else "",
                }
            )
    return rows


def pending_manifest_for_attempt(work_dir: Path, stage: str, attempt: str) -> Path:
    if attempt.startswith("fallback_"):
        number_text = attempt.split("_", 1)[1]
        try:
            number = int(number_text)
        except ValueError:
            number = 0
        return work_dir / f"{stage}.fallback{number:02d}.pending.tsv"
    try:
        number = int(attempt)
    except ValueError:
        number = 0
    return work_dir / f"{stage}.attempt{number:02d}.pending.tsv"


def parse_datetime(value: str) -> dt.datetime | None:
    text = (value or "").strip()
    if not text or text.upper() in {"NA", "N/A", "NONE", "UNKNOWN"}:
        return None
    text = text.replace("Z", "+00:00")
    try:
        parsed = dt.datetime.fromisoformat(text)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=dt.timezone.utc)
    return parsed.astimezone(dt.timezone.utc)


def iso_or_na(value: dt.datetime | None) -> str:
    if value is None:
        return "NA"
    return value.isoformat().replace("+00:00", "Z")


def seconds_between(start: dt.datetime | None, end: dt.datetime | None) -> float | None:
    if start is None or end is None:
        return None
    return max((end - start).total_seconds(), 0.0)


def to_float(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def to_int(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def parse_slurm_duration(value: str) -> float | None:
    text = (value or "").strip()
    if not text or text.upper() in {"NA", "N/A", "NONE", "UNKNOWN", "UNLIMITED"}:
        return None
    days = 0
    if "-" in text:
        day_text, text = text.split("-", 1)
        try:
            days = int(day_text)
        except ValueError:
            return None
    parts = text.split(":")
    try:
        if len(parts) == 3:
            hours = int(parts[0])
            minutes = int(parts[1])
            seconds = float(parts[2])
        elif len(parts) == 2:
            hours = 0
            minutes = int(parts[0])
            seconds = float(parts[1])
        elif len(parts) == 1:
            hours = 0
            minutes = 0
            seconds = float(parts[0])
        else:
            return None
    except ValueError:
        return None
    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def parse_memory_bytes(value: str) -> int | None:
    text = (value or "").strip()
    if not text or text.upper() in {"NA", "N/A", "NONE", "UNKNOWN"}:
        return None
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)([KMGTP]?)(?:I?B)?", text, re.IGNORECASE)
    if not match:
        return None
    number = float(match.group(1))
    suffix = match.group(2).upper()
    powers = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4, "P": 1024**5}
    return int(number * powers[suffix])


def available_sacct_fields() -> list[str]:
    if shutil.which("sacct") is None:
        return []
    status, output = run_command(["sacct", "--helpformat"])
    if status != 0 or not output:
        return ["JobIDRaw", "State", "ExitCode"]
    available = {token.lower(): token for token in re.findall(r"[A-Za-z][A-Za-z0-9]*", output)}
    selected = [field for field in SACCT_CANDIDATES if field.lower() in available]
    for required in ("JobIDRaw", "State", "ExitCode"):
        if required not in selected:
            selected.append(required)
    return selected


def query_sacct(job_id: str, fields: list[str]) -> tuple[bool, list[dict[str, str]], str]:
    if not fields or shutil.which("sacct") is None:
        return False, [], "sacct command is unavailable"
    status, output = run_command(
        ["sacct", "-n", "-P", "-j", job_id, f"--format={','.join(fields)}"]
    )
    if status != 0:
        return False, [], output or f"sacct exited with status {status}"
    rows: list[dict[str, str]] = []
    for line in output.splitlines():
        if not line.strip():
            continue
        values = line.split("|")
        if values and values[-1] == "":
            values.pop()
        if len(values) < len(fields):
            values.extend([""] * (len(fields) - len(values)))
        rows.append(dict(zip(fields, values)))
    if not rows:
        return False, [], "sacct returned no rows"
    return True, rows, ""


def parse_scontrol_line(line: str) -> dict[str, str]:
    return {match.group(1): match.group(2).strip() for match in SCONTROL_FIELD_RE.finditer(line)}


def query_node(node: str) -> dict[str, str]:
    if not node or shutil.which("scontrol") is None:
        return {}
    status, output = run_command(["scontrol", "show", "node", node, "-o"])
    if status != 0 or not output:
        return {}
    return parse_scontrol_line(output.splitlines()[0])


def normalize_node_name(value: str) -> str:
    """Return one concrete node name, rejecting Slurm placeholders/ranges."""
    node = (value or "").strip()
    if not node or node.upper() in {"NA", "N/A", "NONE", "UNKNOWN", "(NULL)"}:
        return ""
    # Every MTD array task requests exactly one node. A compressed/multi-node
    # expression is not a concrete task placement and must not become a fake
    # entry in hpc_node_summary.tsv.
    if any(token in node for token in (",", "[", "]")):
        return ""
    return node


def marker_hardware(marker: dict[str, str]) -> dict[str, str]:
    return {
        "architecture": marker.get("architecture", ""),
        "cpu_model": marker.get("cpu_model", ""),
        "logical_cpus": marker.get("logical_cpus", ""),
        "sockets": marker.get("sockets", ""),
        "cores_per_socket": marker.get("cores_per_socket", ""),
        "threads_per_core": marker.get("threads_per_core", ""),
        "memory_total_kb": marker.get("memory_total_kb", ""),
        "kernel": marker.get("kernel", ""),
        "os_release": marker.get("os_release", ""),
    }


def main() -> int:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    pipeline_output = args.pipeline_output.resolve(strict=False)
    hpc_conf = args.hpc_conf.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []
    job_rows: list[dict[str, object]] = []
    task_rows: list[dict[str, object]] = []
    stage_context: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "job_count": 0,
            "submitted_task_attempts": 0,
            "unique_tasks": set(),
            "successful_tasks": set(),
            "failed_tasks": set(),
            "retry_events": 0,
            "nodes": set(),
            "submit_times": [],
            "start_times": [],
            "end_times": [],
            "queue_wait": [],
            "elapsed": [],
            "allocated_cpu": [],
            "total_cpu": [],
            "max_rss": [],
            "ncpus": [],
        }
    )
    node_context: dict[str, dict[str, object]] = defaultdict(
        lambda: {
            "stages": set(),
            "task_attempts": 0,
            "successful_task_attempts": 0,
            "failed_task_attempts": 0,
            "elapsed": 0.0,
            "elapsed_observations": 0,
            "allocated_cpu": 0.0,
            "allocated_cpu_observations": 0,
            "total_cpu": 0.0,
            "total_cpu_observations": 0,
            "peak_rss": 0,
            "hardware": {},
        }
    )

    job_history_paths = sorted(pipeline_output.rglob("*.job_ids.tsv")) if pipeline_output.exists() else []
    sacct_fields = available_sacct_fields()
    accounting_jobs = 0

    for history_path in job_history_paths:
        stage = history_path.name.removesuffix(".job_ids.tsv")
        work_dir = history_path.parent
        success_dir = work_dir / "success"
        marker_by_task: dict[str, dict[str, str]] = {}
        marker_status_by_task: dict[str, str] = {}
        attempt_marker_by_job_task: dict[tuple[str, str], dict[str, str]] = {}
        attempt_marker_status: dict[tuple[str, str], str] = {}
        if success_dir.is_dir():
            for marker_path in sorted(success_dir.glob("*.success")):
                task_id = marker_path.name.removesuffix(".success")
                marker_by_task[task_id] = read_key_values(marker_path)
                marker_status_by_task[task_id] = "SUCCESS"
            for marker_path in sorted(success_dir.glob("*.failed")):
                task_id = marker_path.name.removesuffix(".failed")
                if task_id not in marker_by_task:
                    marker_by_task[task_id] = read_key_values(marker_path)
                    marker_status_by_task[task_id] = "FAILED"

            attempt_marker_dir = success_dir / "attempts"
            if attempt_marker_dir.is_dir():
                for marker_status, pattern in (("SUCCESS", "*.success"), ("FAILED", "*.failed")):
                    for marker_path in sorted(attempt_marker_dir.glob(pattern)):
                        marker = read_key_values(marker_path)
                        marker_job_id = marker.get("job_id", "")
                        marker_array_task_id = marker.get("array_task_id", "")
                        if marker_job_id and marker_array_task_id:
                            key = (marker_job_id, marker_array_task_id)
                            attempt_marker_by_job_task[key] = marker
                            attempt_marker_status[key] = marker_status

        original_manifests = [
            work_dir / f"{stage}.tasks.tsv",
            work_dir.parent / f"{stage}.tasks.tsv",
        ]
        for manifest_path in original_manifests:
            for item in read_manifest(manifest_path):
                stage_context[stage]["unique_tasks"].add(item["task_id"])

        retry_path = work_dir / f"{stage}.retry.tsv"
        if retry_path.is_file():
            retry_count = sum(
                1
                for line in retry_path.read_text(encoding="utf-8", errors="replace").splitlines()
                if line.strip() and not line.startswith("#")
            )
            stage_context[stage]["retry_events"] = retry_count

        history_entries: list[tuple[str, str, int, str]] = []
        with history_path.open(encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    warnings.append(f"Malformed job history row: {history_path}: {line}")
                    continue
                attempt = fields[0]
                job_id = fields[1]
                task_count = to_int(fields[2]) or 0
                placement = fields[3] if len(fields) > 3 else ""
                if not JOB_ID_RE.fullmatch(job_id):
                    warnings.append(f"Invalid Slurm job ID in {history_path}: {job_id}")
                    continue
                history_entries.append((attempt, job_id, task_count, placement))

        for attempt, job_id, submitted_tasks, placement in history_entries:
            stage_context[stage]["job_count"] += 1
            stage_context[stage]["submitted_task_attempts"] += submitted_tasks
            pending_manifest = pending_manifest_for_attempt(work_dir, stage, attempt)
            pending_rows = read_manifest(pending_manifest)
            for item in pending_rows:
                stage_context[stage]["unique_tasks"].add(item["task_id"])

            sacct_ok, sacct_rows, sacct_warning = query_sacct(job_id, sacct_fields)
            if sacct_ok:
                accounting_jobs += 1
            else:
                warnings.append(f"Job {job_id} ({stage}, attempt {attempt}): {sacct_warning}")

            by_id = {row.get("JobIDRaw", ""): row for row in sacct_rows}
            parent = by_id.get(job_id, {})
            parent_submit = parse_datetime(parent.get("Submit", ""))
            parent_eligible = parse_datetime(parent.get("Eligible", ""))
            parent_start = parse_datetime(parent.get("Start", ""))
            parent_end = parse_datetime(parent.get("End", ""))
            parent_queue = seconds_between(parent_eligible or parent_submit, parent_start)

            job_rows.append(
                {
                    "stage": stage,
                    "attempt": attempt,
                    "job_id": job_id,
                    "submitted_tasks": submitted_tasks,
                    "placement": placement,
                    "accounting_available": int(sacct_ok),
                    "state": parent.get("State", "NA"),
                    "exit_code": parent.get("ExitCode", "NA"),
                    "submit_utc": iso_or_na(parent_submit),
                    "eligible_utc": iso_or_na(parent_eligible),
                    "start_utc": iso_or_na(parent_start),
                    "end_utc": iso_or_na(parent_end),
                    "elapsed_seconds": parent.get("ElapsedRaw", "NA") or "NA",
                    "queue_wait_seconds": f"{parent_queue:.6f}" if parent_queue is not None else "NA",
                    "node_list": parent.get("NodeList", "NA") or "NA",
                    "ncpus": parent.get("NCPUS", "NA") or "NA",
                    "allocated_cpu_seconds": parent.get("CPUTimeRAW", "NA") or "NA",
                    "total_cpu": parent.get("TotalCPU", "NA") or "NA",
                    "max_rss": parent.get("MaxRSS", "NA") or "NA",
                    "req_mem": parent.get("ReqMem", "NA") or "NA",
                    "alloc_tres": parent.get("AllocTRES", "NA") or "NA",
                    "job_history_file": str(history_path),
                    "pending_manifest": str(pending_manifest),
                }
            )

            base_tasks: dict[str, dict[str, str]] = {}
            step_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
            for row in sacct_rows:
                job_id_raw = row.get("JobIDRaw", "")
                array_match = ARRAY_TASK_RE.fullmatch(job_id_raw)
                if array_match and array_match.group("job") == job_id:
                    base_tasks[job_id_raw] = row
                    continue
                step_match = ARRAY_STEP_RE.fullmatch(job_id_raw)
                if step_match and step_match.group("base").startswith(f"{job_id}_"):
                    step_rows[step_match.group("base")].append(row)

            if not base_tasks and pending_rows:
                for array_index, item in enumerate(pending_rows, start=1):
                    base_tasks[f"{job_id}_{array_index}"] = {}

            for base_id in sorted(base_tasks, key=lambda value: int(value.split("_", 1)[1])):
                base = base_tasks[base_id]
                array_task_id = int(base_id.split("_", 1)[1])
                task_id = pending_rows[array_task_id - 1]["task_id"] if 0 < array_task_id <= len(pending_rows) else ""
                attempt_marker_key = (job_id, str(array_task_id))
                marker = attempt_marker_by_job_task.get(attempt_marker_key, {})
                marker_status = attempt_marker_status.get(attempt_marker_key, "")
                if not marker and task_id:
                    marker = marker_by_task.get(task_id, {})
                    marker_matches_job = not marker.get("job_id") or marker.get("job_id") == job_id
                    if not marker_matches_job:
                        marker = {}
                    else:
                        marker_status = marker_status_by_task.get(task_id, "")

                task_steps = step_rows.get(base_id, [])
                task_submit = parse_datetime(base.get("Submit", "")) or parent_submit
                task_eligible = parse_datetime(base.get("Eligible", "")) or parent_eligible
                task_start = parse_datetime(base.get("Start", "")) or parse_datetime(marker.get("started_at", ""))
                task_end = parse_datetime(base.get("End", "")) or parse_datetime(marker.get("finished_at", ""))
                elapsed = to_float(base.get("ElapsedRaw", ""))
                if elapsed is None:
                    elapsed = to_float(marker.get("elapsed_seconds", ""))
                if elapsed is None:
                    elapsed = seconds_between(task_start, task_end)
                queue_wait = seconds_between(task_eligible or task_submit, task_start)
                # The task marker is written on the compute node itself and is
                # therefore the most exact record of placement. sacct is used
                # only when the task could not write a marker (for example, a
                # hard node failure). No configured/eligible-node inventory is
                # consulted here.
                node = normalize_node_name(marker.get("host", "") or marker.get("node", ""))
                if not node:
                    node = normalize_node_name(base.get("NodeList", ""))
                ncpus = to_int(base.get("NCPUS", ""))
                if ncpus is None:
                    ncpus = to_int(marker.get("slurm_cpus_per_task", "")) or to_int(marker.get("threads", ""))
                allocated_cpu = to_float(base.get("CPUTimeRAW", ""))
                if allocated_cpu is None and elapsed is not None and ncpus is not None:
                    allocated_cpu = elapsed * ncpus

                total_cpu = parse_slurm_duration(base.get("TotalCPU", ""))
                step_cpu_values = [parse_slurm_duration(row.get("TotalCPU", "")) for row in task_steps]
                step_cpu_sum = sum(value for value in step_cpu_values if value is not None)
                if total_cpu is None and step_cpu_sum > 0:
                    total_cpu = step_cpu_sum

                rss_candidates = [parse_memory_bytes(base.get("MaxRSS", ""))]
                rss_candidates.extend(parse_memory_bytes(row.get("MaxRSS", "")) for row in task_steps)
                max_rss_bytes = max((value for value in rss_candidates if value is not None), default=None)

                state = base.get("State", "") or marker_status or "UNKNOWN"
                exit_code = base.get("ExitCode", "") or marker.get("exit_code", "") or "NA"
                success = state.startswith("COMPLETED") or marker_status == "SUCCESS"
                failed = (
                    not success
                    and state not in {"", "UNKNOWN", "PENDING", "RUNNING"}
                ) or marker_status == "FAILED"

                if task_id:
                    stage_context[stage]["unique_tasks"].add(task_id)
                    if success:
                        stage_context[stage]["successful_tasks"].add(task_id)
                    elif failed:
                        stage_context[stage]["failed_tasks"].add(task_id)
                if node:
                    stage_context[stage]["nodes"].add(node)
                    node_context[node]["stages"].add(stage)
                    node_context[node]["task_attempts"] += 1
                    if success:
                        node_context[node]["successful_task_attempts"] += 1
                    elif failed:
                        node_context[node]["failed_task_attempts"] += 1
                    if elapsed is not None:
                        node_context[node]["elapsed"] += elapsed
                        node_context[node]["elapsed_observations"] += 1
                    if allocated_cpu is not None:
                        node_context[node]["allocated_cpu"] += allocated_cpu
                        node_context[node]["allocated_cpu_observations"] += 1
                    if total_cpu is not None:
                        node_context[node]["total_cpu"] += total_cpu
                        node_context[node]["total_cpu_observations"] += 1
                    if max_rss_bytes is not None:
                        node_context[node]["peak_rss"] = max(node_context[node]["peak_rss"], max_rss_bytes)
                    hardware = marker_hardware(marker)
                    for key, value in hardware.items():
                        if value and not node_context[node]["hardware"].get(key):
                            node_context[node]["hardware"][key] = value

                if task_submit:
                    stage_context[stage]["submit_times"].append(task_submit)
                if task_start:
                    stage_context[stage]["start_times"].append(task_start)
                if task_end:
                    stage_context[stage]["end_times"].append(task_end)
                if queue_wait is not None:
                    stage_context[stage]["queue_wait"].append(queue_wait)
                if elapsed is not None:
                    stage_context[stage]["elapsed"].append(elapsed)
                if allocated_cpu is not None:
                    stage_context[stage]["allocated_cpu"].append(allocated_cpu)
                if total_cpu is not None:
                    stage_context[stage]["total_cpu"].append(total_cpu)
                if max_rss_bytes is not None:
                    stage_context[stage]["max_rss"].append(max_rss_bytes)
                if ncpus is not None:
                    stage_context[stage]["ncpus"].append(ncpus)

                task_rows.append(
                    {
                        "stage": stage,
                        "attempt": attempt,
                        "job_id": job_id,
                        "array_task_id": array_task_id,
                        "task_id": task_id or "NA",
                        "state": state or "UNKNOWN",
                        "exit_code": exit_code,
                        "node": node or "NA",
                        "submit_utc": iso_or_na(task_submit),
                        "eligible_utc": iso_or_na(task_eligible),
                        "start_utc": iso_or_na(task_start),
                        "end_utc": iso_or_na(task_end),
                        "elapsed_seconds": f"{elapsed:.6f}" if elapsed is not None else "NA",
                        "queue_wait_seconds": f"{queue_wait:.6f}" if queue_wait is not None else "NA",
                        "ncpus": ncpus if ncpus is not None else "NA",
                        "allocated_cpu_seconds": f"{allocated_cpu:.6f}" if allocated_cpu is not None else "NA",
                        "total_cpu_seconds": f"{total_cpu:.6f}" if total_cpu is not None else "NA",
                        "max_rss_bytes": max_rss_bytes if max_rss_bytes is not None else "NA",
                        "req_mem": base.get("ReqMem", "NA") or "NA",
                        "alloc_tres": base.get("AllocTRES", "NA") or "NA",
                        "marker_status": marker_status or "NA",
                        "detected_threads": marker.get("threads", marker.get("detected_threads", "NA")),
                        "memory_available_kb": marker.get("memory_kb", marker.get("memory_available_kb", "NA")),
                        "success_marker_job_id": marker.get("job_id", "NA"),
                        "pending_manifest": str(pending_manifest),
                    }
                )

        for task_id, status in marker_status_by_task.items():
            stage_context[stage]["unique_tasks"].add(task_id)
            if status == "SUCCESS":
                stage_context[stage]["successful_tasks"].add(task_id)
            elif status == "FAILED":
                stage_context[stage]["failed_tasks"].add(task_id)

    stage_rows: list[dict[str, object]] = []
    all_submit: list[dt.datetime] = []
    all_start: list[dt.datetime] = []
    all_end: list[dt.datetime] = []
    all_queue: list[float] = []
    all_elapsed: list[float] = []
    all_allocated_cpu: list[float] = []
    all_total_cpu: list[float] = []
    all_max_rss: list[int] = []
    all_nodes: set[str] = set()
    all_unique_tasks: set[tuple[str, str]] = set()
    all_success_tasks: set[tuple[str, str]] = set()
    all_failed_tasks: set[tuple[str, str]] = set()
    total_retry_events = 0

    for stage in sorted(stage_context):
        context = stage_context[stage]
        submits = context["submit_times"]
        starts = context["start_times"]
        ends = context["end_times"]
        queue_values = context["queue_wait"]
        elapsed_values = context["elapsed"]
        allocated_cpu_values = context["allocated_cpu"]
        total_cpu_values = context["total_cpu"]
        rss_values = context["max_rss"]
        ncpus_values = context["ncpus"]
        unique_tasks = context["unique_tasks"]
        successful_tasks = context["successful_tasks"]
        failed_tasks = context["failed_tasks"] - successful_tasks
        nodes = context["nodes"]
        stage_makespan = seconds_between(min(starts) if starts else None, max(ends) if ends else None)

        stage_rows.append(
            {
                "stage": stage,
                "job_count": context["job_count"],
                "submitted_task_attempts": context["submitted_task_attempts"],
                "unique_tasks": len(unique_tasks),
                "successful_unique_tasks": len(successful_tasks),
                "failed_unique_tasks": len(failed_tasks),
                "retry_events": context["retry_events"],
                "nodes_used": len(nodes),
                "node_names": ",".join(sorted(nodes)) or "NA",
                "earliest_submit_utc": iso_or_na(min(submits) if submits else None),
                "earliest_start_utc": iso_or_na(min(starts) if starts else None),
                "latest_end_utc": iso_or_na(max(ends) if ends else None),
                "compute_makespan_seconds": f"{stage_makespan:.6f}" if stage_makespan is not None else "NA",
                "sum_queue_wait_seconds": f"{sum(queue_values):.6f}" if queue_values else "NA",
                "max_queue_wait_seconds": f"{max(queue_values):.6f}" if queue_values else "NA",
                "sum_task_elapsed_seconds": f"{sum(elapsed_values):.6f}" if elapsed_values else "NA",
                "allocated_cpu_seconds": f"{sum(allocated_cpu_values):.6f}" if allocated_cpu_values else "NA",
                "total_cpu_seconds": f"{sum(total_cpu_values):.6f}" if total_cpu_values else "NA",
                "peak_task_max_rss_bytes": max(rss_values) if rss_values else "NA",
                "min_ncpus_per_task": min(ncpus_values) if ncpus_values else "NA",
                "max_ncpus_per_task": max(ncpus_values) if ncpus_values else "NA",
            }
        )

        all_submit.extend(submits)
        all_start.extend(starts)
        all_end.extend(ends)
        all_queue.extend(queue_values)
        all_elapsed.extend(elapsed_values)
        all_allocated_cpu.extend(allocated_cpu_values)
        all_total_cpu.extend(total_cpu_values)
        all_max_rss.extend(rss_values)
        all_nodes.update(nodes)
        all_unique_tasks.update((stage, task_id) for task_id in unique_tasks)
        all_success_tasks.update((stage, task_id) for task_id in successful_tasks)
        all_failed_tasks.update((stage, task_id) for task_id in failed_tasks)
        total_retry_events += int(context["retry_events"])

    node_rows: list[dict[str, object]] = []
    for node in sorted(node_context):
        context = node_context[node]
        hardware = context["hardware"]
        scontrol = query_node(node)
        node_rows.append(
            {
                "node": node,
                "node_hostname": scontrol.get("NodeHostName", "NA"),
                "architecture": hardware.get("architecture") or scontrol.get("Arch", "NA"),
                "cpu_model": hardware.get("cpu_model", "NA"),
                "logical_cpus": hardware.get("logical_cpus") or scontrol.get("CPUTot", "NA"),
                "sockets": hardware.get("sockets") or scontrol.get("Sockets", "NA"),
                "cores_per_socket": hardware.get("cores_per_socket") or scontrol.get("CoresPerSocket", "NA"),
                "threads_per_core": hardware.get("threads_per_core") or scontrol.get("ThreadsPerCore", "NA"),
                "memory_total_kb": hardware.get("memory_total_kb") or (
                    str((to_int(scontrol.get("RealMemory", "")) or 0) * 1024)
                    if to_int(scontrol.get("RealMemory", "")) is not None
                    else "NA"
                ),
                "kernel": hardware.get("kernel", "NA"),
                "os_release": hardware.get("os_release") or scontrol.get("OS", "NA"),
                "partitions": scontrol.get("Partitions", "NA"),
                "state_at_collection": scontrol.get("State", "NA"),
                "stages": ",".join(sorted(context["stages"])),
                "task_attempts": context["task_attempts"],
                "successful_task_attempts": context["successful_task_attempts"],
                "failed_task_attempts": context["failed_task_attempts"],
                "sum_task_elapsed_seconds": (
                    f"{context['elapsed']:.6f}" if context["elapsed_observations"] else "NA"
                ),
                "allocated_cpu_seconds": (
                    f"{context['allocated_cpu']:.6f}"
                    if context["allocated_cpu_observations"]
                    else "NA"
                ),
                "total_cpu_seconds": (
                    f"{context['total_cpu']:.6f}"
                    if context["total_cpu_observations"]
                    else "NA"
                ),
                "peak_task_max_rss_bytes": context["peak_rss"] or "NA",
            }
        )

    global_makespan = seconds_between(min(all_start) if all_start else None, max(all_end) if all_end else None)
    failed_final = all_failed_tasks - all_success_tasks
    if not job_history_paths:
        collection_status = "NO_JOBS"
        warnings.append(f"No *.job_ids.tsv files were found under {pipeline_output}")
    elif accounting_jobs == 0:
        collection_status = "PARTIAL"
    elif warnings:
        collection_status = "PARTIAL"
    else:
        collection_status = "PASS"

    summary_rows: list[tuple[str, object, str]] = [
        ("hpc_collection_status", collection_status, "PASS, PARTIAL, or NO_JOBS"),
        ("hpc_controller_host", args.controller_host, "text"),
        ("hpc_controller_threads", args.controller_threads, "count"),
        ("hpc_configuration", str(hpc_conf), "path"),
        ("slurm_accounting_available", int(accounting_jobs > 0), "0/1"),
        ("slurm_accounting_jobs_collected", accounting_jobs, "count"),
        ("slurm_stage_count", len(stage_rows), "count"),
        ("slurm_job_count", len(job_rows), "count"),
        ("slurm_array_task_attempt_count", len(task_rows), "count"),
        ("slurm_unique_task_count", len(all_unique_tasks), "count"),
        ("slurm_successful_unique_task_count", len(all_success_tasks), "count"),
        ("slurm_failed_unique_task_count", len(failed_final), "count"),
        ("slurm_retry_event_count", total_retry_events, "count"),
        ("slurm_nodes_used", len(all_nodes), "count"),
        ("slurm_node_names", ",".join(sorted(all_nodes)) or "NA", "comma-separated"),
        ("slurm_earliest_submit_utc", iso_or_na(min(all_submit) if all_submit else None), "UTC"),
        ("slurm_earliest_start_utc", iso_or_na(min(all_start) if all_start else None), "UTC"),
        ("slurm_latest_end_utc", iso_or_na(max(all_end) if all_end else None), "UTC"),
        ("slurm_compute_makespan_seconds", f"{global_makespan:.6f}" if global_makespan is not None else "NA", "seconds"),
        ("slurm_sum_queue_wait_seconds", f"{sum(all_queue):.6f}" if all_queue else "NA", "seconds across task attempts"),
        ("slurm_max_queue_wait_seconds", f"{max(all_queue):.6f}" if all_queue else "NA", "seconds"),
        ("slurm_sum_task_elapsed_seconds", f"{sum(all_elapsed):.6f}" if all_elapsed else "NA", "seconds across task attempts"),
        ("slurm_allocated_cpu_seconds", f"{sum(all_allocated_cpu):.6f}" if all_allocated_cpu else "NA", "allocated CPU-seconds"),
        ("slurm_total_cpu_seconds", f"{sum(all_total_cpu):.6f}" if all_total_cpu else "NA", "reported consumed CPU-seconds"),
        ("slurm_peak_task_max_rss_bytes", max(all_max_rss) if all_max_rss else "NA", "bytes"),
        ("slurm_collection_warning_count", len(warnings), "count"),
    ]

    write_tsv(
        run_dir / "hpc_jobs.tsv",
        [
            "stage", "attempt", "job_id", "submitted_tasks", "placement",
            "accounting_available", "state", "exit_code", "submit_utc",
            "eligible_utc", "start_utc", "end_utc", "elapsed_seconds",
            "queue_wait_seconds", "node_list", "ncpus", "allocated_cpu_seconds",
            "total_cpu", "max_rss", "req_mem", "alloc_tres",
            "job_history_file", "pending_manifest",
        ],
        job_rows,
    )
    write_tsv(
        run_dir / "hpc_tasks.tsv",
        [
            "stage", "attempt", "job_id", "array_task_id", "task_id", "state",
            "exit_code", "node", "submit_utc", "eligible_utc", "start_utc",
            "end_utc", "elapsed_seconds", "queue_wait_seconds", "ncpus",
            "allocated_cpu_seconds", "total_cpu_seconds", "max_rss_bytes",
            "req_mem", "alloc_tres", "marker_status", "detected_threads",
            "memory_available_kb", "success_marker_job_id", "pending_manifest",
        ],
        task_rows,
    )
    write_tsv(
        run_dir / "hpc_stage_summary.tsv",
        [
            "stage", "job_count", "submitted_task_attempts", "unique_tasks",
            "successful_unique_tasks", "failed_unique_tasks", "retry_events",
            "nodes_used", "node_names", "earliest_submit_utc",
            "earliest_start_utc", "latest_end_utc", "compute_makespan_seconds",
            "sum_queue_wait_seconds", "max_queue_wait_seconds",
            "sum_task_elapsed_seconds", "allocated_cpu_seconds",
            "total_cpu_seconds", "peak_task_max_rss_bytes",
            "min_ncpus_per_task", "max_ncpus_per_task",
        ],
        stage_rows,
    )
    write_tsv(
        run_dir / "hpc_node_summary.tsv",
        [
            "node", "node_hostname", "architecture", "cpu_model", "logical_cpus",
            "sockets", "cores_per_socket", "threads_per_core", "memory_total_kb",
            "kernel", "os_release", "partitions", "state_at_collection",
            "stages", "task_attempts",
            "successful_task_attempts", "failed_task_attempts",
            "sum_task_elapsed_seconds", "allocated_cpu_seconds",
            "total_cpu_seconds", "peak_task_max_rss_bytes",
        ],
        node_rows,
    )
    with (run_dir / "hpc_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value", "unit_or_definition"])
        writer.writerows(summary_rows)

    warning_path = run_dir / "hpc_collection_warnings.txt"
    if warnings:
        warning_path.write_text("\n".join(warnings) + "\n", encoding="utf-8")
    elif warning_path.exists():
        warning_path.unlink()

    print(f"[HPC BENCHMARK] Collection status: {collection_status}")
    print(f"[HPC BENCHMARK] Jobs: {len(job_rows)}")
    print(f"[HPC BENCHMARK] Task attempts: {len(task_rows)}")
    print(f"[HPC BENCHMARK] Nodes: {len(all_nodes)}")
    print(f"[HPC BENCHMARK] Summary: {run_dir / 'hpc_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
