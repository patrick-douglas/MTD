#!/usr/bin/env python3
"""Generate a focused report for a Create_custom_host.sh benchmark run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import re
from typing import Iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--repo", required=True)
    parser.add_argument("--taxid", required=True)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--machine", required=True)
    parser.add_argument("--run-number", required=True)
    parser.add_argument("--logical-cpus", required=True)
    parser.add_argument("--scientific-name", default="")
    parser.add_argument("--reference-taxid", default="")
    parser.add_argument("--orgdb-package", default="")
    parser.add_argument("--builder-status", type=int, required=True)
    parser.add_argument("--skip-orgdb", type=int, choices=(0, 1), required=True)
    return parser.parse_args()


def read_metric_table(path: Path) -> dict[str, str]:
    metrics: dict[str, str] = {}
    if not path.exists():
        return metrics
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = row.get("metric", "")
            if name:
                metrics[name] = row.get("value", "")
    return metrics


def read_cache_state(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def human_seconds(value: str | float | None) -> str:
    try:
        seconds = float(value) if value is not None else 0.0
    except (TypeError, ValueError):
        return "NA"
    whole = int(round(seconds))
    hours, remainder = divmod(whole, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def human_bytes(value: int | None) -> str:
    if value is None:
        return "NA"
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    amount = float(value)
    for unit in units:
        if amount < 1024 or unit == units[-1]:
            return f"{amount:.2f} {unit}"
        amount /= 1024
    return f"{value} B"


def tree_stats(path: Path) -> tuple[bool, int, int]:
    if not path.exists():
        return False, 0, 0
    if path.is_file():
        try:
            return True, 1, path.stat().st_size
        except OSError:
            return True, 1, 0
    count = 0
    total = 0
    for root, _, files in os.walk(path):
        for name in files:
            candidate = Path(root) / name
            count += 1
            try:
                total += candidate.stat().st_size
            except OSError:
                pass
    return True, count, total


def nonempty(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def any_nonempty(directory: Path, patterns: Iterable[str]) -> bool:
    if not directory.is_dir():
        return False
    for pattern in patterns:
        if any(nonempty(path) for path in directory.glob(pattern)):
            return True
    return False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    repo = Path(args.repo).resolve()
    cache = Path(args.cache).resolve()
    taxid = args.taxid

    metrics = read_metric_table(run_dir / "summary.tsv")
    cache_before = read_cache_state(run_dir / "cache_state_before.tsv")

    kraken = repo / f"kraken2DB_{taxid}"
    ref = repo / f"ref_{taxid}"
    hisat = repo / f"hisat2_index_{taxid}"
    blast = repo / f"blastdb_{taxid}"
    orgdb_build = repo / "build" / "orgdb_gold" / taxid
    functional = ref / "functional_annotation"
    species_cache = cache / "Customized_hosts" / taxid
    taxonomy_cache = cache / "Kraken2_taxonomy_cache"
    eggnog_cache = cache / "eggNOG" / "emapperdb-5.0.2"
    ncbi_cache = cache / "NCBI_gene" / "gene2ensembl.gz"

    validations: list[tuple[str, Path, bool, str]] = [
        ("Kraken2 hash", kraken / "hash.k2d", nonempty(kraken / "hash.k2d"), "required"),
        ("Kraken2 options", kraken / "opts.k2d", nonempty(kraken / "opts.k2d"), "required"),
        ("Kraken2 taxonomy", kraken / "taxo.k2d", nonempty(kraken / "taxo.k2d"), "required"),
        (
            "Kraken2 inspect",
            kraken / "kraken2_inspect_summary.txt",
            nonempty(kraken / "kraken2_inspect_summary.txt"),
            "required",
        ),
        (
            "Clean host FASTA",
            kraken / f"genome_{taxid}.fa",
            nonempty(kraken / f"genome_{taxid}.fa"),
            "required",
        ),
        (
            "GTF annotation",
            ref / f"ref_{taxid}.gtf.gz",
            nonempty(ref / f"ref_{taxid}.gtf.gz"),
            "required",
        ),
        (
            "HISAT2 index",
            hisat,
            any_nonempty(hisat, ("genome_tran*.ht2", "genome_tran*.ht2l")),
            "required",
        ),
        (
            "BLAST nucleotide DB",
            blast,
            any_nonempty(blast, (f"blastdb_{taxid}.n*", f"blastdb_{taxid}.*db")),
            "required",
        ),
        (
            "eggNOG master GMT",
            functional / f"custom_taxid_{taxid}_eggNOG_GO_master.gmt",
            nonempty(functional / f"custom_taxid_{taxid}_eggNOG_GO_master.gmt"),
            "required",
        ),
        (
            "Functional manifest",
            functional / "functional_annotation_manifest.tsv",
            nonempty(functional / "functional_annotation_manifest.tsv"),
            "required",
        ),
        (
            "Functional checksums",
            functional / "checksums.sha256",
            nonempty(functional / "checksums.sha256"),
            "required",
        ),
        (
            "OrgDb build directory",
            orgdb_build,
            orgdb_build.is_dir(),
            "required",
        ),
    ]

    if args.skip_orgdb == 0 and args.orgdb_package:
        pkg = repo / "custom_R_libs" / args.orgdb_package
        validations.append(("Installed custom OrgDb", pkg, pkg.is_dir(), "required"))
    elif args.skip_orgdb == 1:
        validations.append(
            ("Installed custom OrgDb", repo / "custom_R_libs", True, "skipped by option")
        )

    validation_path = run_dir / "output_validation.tsv"
    with validation_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["check", "status", "requirement", "path"])
        for name, path, passed, requirement in validations:
            writer.writerow([name, "PASS" if passed else "FAIL", requirement, str(path)])

    inventory_targets = [
        ("kraken2_database", kraken),
        ("reference_annotation", ref),
        ("hisat2_index", hisat),
        ("blast_database", blast),
        ("orgdb_build", orgdb_build),
        ("custom_r_library", repo / "custom_R_libs"),
        ("species_cache", species_cache),
        ("taxonomy_cache", taxonomy_cache),
        ("eggnog_cache", eggnog_cache),
        ("ncbi_gene_cache", ncbi_cache),
    ]
    inventory_path = run_dir / "output_inventory.tsv"
    with inventory_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["component", "exists", "file_count", "size_bytes", "size_human", "path"])
        for name, path in inventory_targets:
            exists, files, size = tree_stats(path)
            writer.writerow([name, int(exists), files, size, human_bytes(size), str(path)])

    cache_after_path = run_dir / "cache_state_after.tsv"
    after_rows = []
    for name, path in (
        ("species_reference", species_cache),
        ("kraken_taxonomy", taxonomy_cache),
        ("eggnog_database", eggnog_cache),
        ("ncbi_gene2ensembl", ncbi_cache),
    ):
        exists, files, size = tree_stats(path)
        after_rows.append((name, path, exists, files, size))
    with cache_after_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["component", "path", "exists", "size_bytes", "size_human", "file_count"])
        for name, path, exists, files, size in after_rows:
            writer.writerow([name, str(path), int(exists), size, human_bytes(size), files])

    console_path = run_dir / "console.log"
    console_text = ""
    if console_path.exists():
        console_text = console_path.read_text(encoding="utf-8", errors="replace")

    diagnostic_patterns = re.compile(
        r"(?i)(\[ERROR\]|\[WARNING\]|traceback|segmentation fault|killed|"
        r"no space left|out of memory|cannot allocate|failed|corrupt|incomplete)"
    )
    diagnostic_lines = [
        f"{number}\t{line}"
        for number, line in enumerate(console_text.splitlines(), start=1)
        if diagnostic_patterns.search(line)
    ]
    (run_dir / "diagnostic_hits.tsv").write_text(
        "line\ttext\n" + "\n".join(diagnostic_lines) + ("\n" if diagnostic_lines else ""),
        encoding="utf-8",
    )

    tail_lines = console_text.splitlines()[-200:]
    (run_dir / "final_console_tail.txt").write_text(
        "\n".join(tail_lines) + ("\n" if tail_lines else ""),
        encoding="utf-8",
    )

    required_failures = [
        name
        for name, _, passed, requirement in validations
        if requirement == "required" and not passed
    ]
    exit_status = metrics.get("exit_status", str(args.builder_status))
    build_ok = str(exit_status) == "0" and args.builder_status == 0 and not required_failures

    wall = metrics.get("wall_time_seconds", "NA")
    summary_rows = [
        ("benchmark_result", "PASS" if build_ok else "FAIL"),
        ("builder_exit_status", str(args.builder_status)),
        ("monitor_exit_status", str(exit_status)),
        ("machine", args.machine),
        ("taxid", taxid),
        ("scientific_name", args.scientific_name or "NA"),
        ("reference_taxid", args.reference_taxid or "NA"),
        ("orgdb_package", args.orgdb_package or "NA"),
        ("run_number", args.run_number),
        ("logical_cpus_nproc", args.logical_cpus),
        ("wall_time_seconds", wall),
        ("wall_time_hh_mm_ss", human_seconds(wall)),
        ("mean_process_tree_cpu_total_capacity", metrics.get("mean_process_tree_cpu_total_capacity", "NA")),
        ("peak_process_tree_cpu_total_capacity", metrics.get("peak_process_tree_cpu_total_capacity", "NA")),
        ("peak_process_tree_cpu_one_core", metrics.get("peak_process_tree_cpu_one_core", "NA")),
        ("peak_process_tree_rss_gib", metrics.get("peak_process_tree_rss", "NA")),
        ("peak_system_memory_used_gib", metrics.get("peak_system_memory_used", "NA")),
        ("minimum_system_memory_available_gib", metrics.get("minimum_system_memory_available", "NA")),
        ("total_physical_disk_read_gib", metrics.get("total_physical_disk_read", "NA")),
        ("total_physical_disk_write_gib", metrics.get("total_physical_disk_write", "NA")),
        ("total_network_received_gib", metrics.get("total_network_received", "NA")),
        ("total_network_transmitted_gib", metrics.get("total_network_transmitted", "NA")),
        ("peak_reported_temperature_c", metrics.get("peak_reported_temperature", "NA")),
        ("required_output_failures", ",".join(required_failures) if required_failures else "none"),
        ("diagnostic_hit_count", str(len(diagnostic_lines))),
        ("create_custom_host_sha256", sha256(repo / "Create_custom_host.sh")),
    ]

    summary_tsv = run_dir / "custom_host_summary.tsv"
    with summary_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerows(summary_rows)

    cache_description = ", ".join(
        f"{row.get('component', '?')}={row.get('state', '?')}" for row in cache_before
    ) or "unavailable"

    summary_text = f"""MTD Create_custom_host benchmark
================================

Result:              {"PASS" if build_ok else "FAIL"}
Machine:             {args.machine}
TaxID:               {taxid}
Scientific name:     {args.scientific_name or "NA"}
Reference TaxID:     {args.reference_taxid or "NA"}
Logical CPUs:        {args.logical_cpus}
Total wall time:     {human_seconds(wall)} ({wall} seconds)
Cache state before:  {cache_description}

Peak/total resources
--------------------
Mean process CPU:    {metrics.get("mean_process_tree_cpu_total_capacity", "NA")} % total capacity
Peak process CPU:    {metrics.get("peak_process_tree_cpu_total_capacity", "NA")} % total capacity
Peak one-core CPU:   {metrics.get("peak_process_tree_cpu_one_core", "NA")} %
Peak process RSS:    {metrics.get("peak_process_tree_rss", "NA")} GiB
Peak system memory:  {metrics.get("peak_system_memory_used", "NA")} GiB
Minimum free memory: {metrics.get("minimum_system_memory_available", "NA")} GiB
Disk read:           {metrics.get("total_physical_disk_read", "NA")} GiB
Disk written:        {metrics.get("total_physical_disk_write", "NA")} GiB
Network received:    {metrics.get("total_network_received", "NA")} GiB
Network transmitted: {metrics.get("total_network_transmitted", "NA")} GiB
Peak temperature:    {metrics.get("peak_reported_temperature", "NA")} C

Output validation
-----------------
Required failures:   {", ".join(required_failures) if required_failures else "none"}
Diagnostic hits:     {len(diagnostic_lines)}

Important files
---------------
Generic metrics:     {run_dir / "summary.tsv"}
Focused metrics:     {summary_tsv}
Output validation:   {validation_path}
Output inventory:    {inventory_path}
Cache before:        {run_dir / "cache_state_before.tsv"}
Cache after:         {cache_after_path}
Console log:         {console_path}
Diagnostics:         {run_dir / "diagnostic_hits.tsv"}
"""
    (run_dir / "custom_host_summary.txt").write_text(summary_text, encoding="utf-8")

    if not build_ok:
        failure = [
            "Create_custom_host benchmark failure report",
            "===========================================",
            f"Builder exit status: {args.builder_status}",
            f"Monitor exit status: {exit_status}",
            f"Missing/failed required outputs: {', '.join(required_failures) if required_failures else 'none'}",
            f"Diagnostic lines found: {len(diagnostic_lines)}",
            "",
            "See:",
            str(run_dir / "diagnostic_hits.tsv"),
            str(run_dir / "final_console_tail.txt"),
            str(run_dir / "console.log"),
        ]
        (run_dir / "failure_report.txt").write_text("\n".join(failure) + "\n", encoding="utf-8")

    print(summary_text)
    return 0 if build_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
