#!/usr/bin/env python3

import argparse
import csv
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


# Microbiome-read extraction is highly I/O intensive because every
# KrakenTools subprocess reads Kraken and FASTA/FASTQ files. Keep the
# internal concurrency conservative to avoid saturating the storage
# device or opening too many large files simultaneously.
MAX_PARALLEL_EXTRACTION_JOBS = 5


KNOWN_RANKING_COLUMNS = {
    "rank_position",
    "taxon",
    "rank",
    "n_samples",
    "prevalence_n",
    "prevalence_percent",
    "mean_abundance",
    "median_abundance",
    "max_abundance",
    "total_abundance",
    "prevalence_class",
    "importance_interpretation",
    "importance_score",
    "prevalence_score",
    "mean_abundance_score",
    "max_abundance_score",
    "total_abundance_score",
    "samples_present",
    "groups_present",
    "n_samples_present",
    "n_groups_present",
    "distribution_scope",
}


def detect_delimiter(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        first = handle.readline()

    counts = {
        "\t": first.count("\t"),
        ",": first.count(","),
        ";": first.count(";"),
    }

    return max(counts, key=counts.get)


def read_table(path: Path):
    delimiter = detect_delimiter(path)

    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = list(reader)
        header = reader.fieldnames or []

    return header, rows


def clean_text(x: str) -> str:
    if x is None:
        return ""
    return str(x).strip().replace("\r", "")


def safe_name(x: str, max_len: int = 90) -> str:
    x = clean_text(x)
    x = x.replace(" ", "_")
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", x)
    x = re.sub(r"_+", "_", x)
    x = x.strip("._-")
    if not x:
        x = "unknown_taxon"
    return x[:max_len]


def to_float(x, default=0.0) -> float:
    try:
        x = clean_text(x)
        if x == "" or x.lower() in {"na", "nan", "none", "null"}:
            return default
        return float(x)
    except Exception:
        return default


def to_int_from_float(x, default=None):
    try:
        return int(float(clean_text(x)))
    except Exception:
        return default


def read_samplesheet_samples(samplesheet: Path):
    header, rows = read_table(samplesheet)
    if not header:
        raise SystemExit(f"[ERROR] Samplesheet has no header: {samplesheet}")

    sample_col = header[0]
    samples = []

    for row in rows:
        sample = clean_text(row.get(sample_col, ""))
        if sample:
            samples.append(sample)

    samples = list(dict.fromkeys(samples))

    if not samples:
        raise SystemExit(f"[ERROR] No samples found in samplesheet: {samplesheet}")

    return samples


def build_taxon_to_taxid(bracken_table: Path):
    header, rows = read_table(bracken_table)

    lower_to_real = {h.lower(): h for h in header}

    name_col = None
    for candidate in ["name", "taxon", "taxa", "scientific_name", "species", "genus"]:
        if candidate in lower_to_real:
            name_col = lower_to_real[candidate]
            break

    taxid_col = None
    for candidate in ["taxonomy_id", "taxid", "tax_id"]:
        if candidate in lower_to_real:
            taxid_col = lower_to_real[candidate]
            break

    if name_col is None or taxid_col is None:
        raise SystemExit(
            "[ERROR] Could not detect name/taxid columns in Bracken table.\n"
            f"File: {bracken_table}\n"
            f"Header: {header}"
        )

    mapping = {}

    for row in rows:
        name = clean_text(row.get(name_col, ""))
        taxid = clean_text(row.get(taxid_col, ""))

        if not name or not taxid:
            continue

        if not re.fullmatch(r"[0-9]+", taxid):
            continue

        if name not in mapping:
            mapping[name] = taxid

    if not mapping:
        raise SystemExit(f"[ERROR] No taxon-to-taxid mapping recovered from: {bracken_table}")

    return mapping


def detect_sequence_format(seq_file: Path):
    with seq_file.open("r", encoding="utf-8", errors="replace") as handle:
        first = handle.readline().strip()

    if first.startswith(">"):
        return "fasta"

    if first.startswith("@"):
        return "fastq"

    return "unknown"


def count_sequences(seq_file: Path, seq_format: str):
    if not seq_file.exists() or seq_file.stat().st_size == 0:
        return 0

    if seq_format == "fasta":
        count = 0
        with seq_file.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(">"):
                    count += 1
        return count

    if seq_format == "fastq":
        line_count = 0
        with seq_file.open("r", encoding="utf-8", errors="replace") as handle:
            for _ in handle:
                line_count += 1
        return line_count // 4

    return 0


def get_sample_columns(ranked_header, samplesheet_samples):
    header_set = set(ranked_header)

    sample_cols = [s for s in samplesheet_samples if s in header_set]

    if sample_cols:
        return sample_cols

    # Fallback: any non-known numeric-looking column after distribution columns.
    sample_cols = []

    for col in ranked_header:
        if col in KNOWN_RANKING_COLUMNS:
            continue
        sample_cols.append(col)

    return sample_cols


def select_ranked_rows(rows, top_n: int):
    def rank_key(row):
        rank = to_int_from_float(row.get("rank_position", ""), default=None)
        if rank is not None:
            return rank

        score = to_float(row.get("importance_score", "0"), default=0.0)
        return -score

    sorted_rows = sorted(rows, key=rank_key)

    if top_n and top_n > 0:
        sorted_rows = sorted_rows[:top_n]

    return sorted_rows


def extract_reads_for_taxon_sample(
    krakentools_script: Path,
    kraken_file: Path,
    seq_file1: Path,
    seq_file2,
    report_file: Path,
    taxid: str,
    out_file1: Path,
    out_file2,
    include_children: bool,
    seq_format: str,
    read_layout: str,
):
    cmd = [
        sys.executable,
        str(krakentools_script),
        "-k",
        str(kraken_file),
        "-s1",
        str(seq_file1),
        "-r",
        str(report_file),
        "-t",
        str(taxid),
        "-o",
        str(out_file1),
    ]

    if read_layout == "pe":
        if seq_file2 is None:
            raise ValueError(
                "Paired-end extraction requires seq_file2."
            )

        if out_file2 is None:
            raise ValueError(
                "Paired-end extraction requires out_file2."
            )

        cmd.extend(
            [
                "-s2",
                str(seq_file2),
                "-o2",
                str(out_file2),
            ]
        )

    if include_children:
        cmd.append("--include-children")

    if seq_format == "fastq":
        cmd.append("--fastq-output")

    env = os.environ.copy()
    env["PYTHONNOUSERSITE"] = "1"
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )

    return result, cmd


def concatenate_files(files, output_file):
    with output_file.open("wb") as out:
        for f in files:
            if not f.exists() or f.stat().st_size == 0:
                continue
            with f.open("rb") as inp:
                out.write(inp.read())


def main():
    parser = argparse.ArgumentParser(
        description="Extract Kraken-classified reads for detected microbiome taxa."
    )

    parser.add_argument("--ranked_table", required=True, help="detected_microbiome_*_ranked_with_samples.tsv")
    parser.add_argument("--bracken_table", required=True, help="Combined Bracken table, e.g. bracken_species_all")
    parser.add_argument("--samplesheet", required=True, help="samplesheet.csv")
    parser.add_argument(
        "--temp_dir",
        required=True,
        help="MTD temp directory containing Kraken and non-host files",
    )
    parser.add_argument(
        "--read-layout",
        "--read_layout",
        dest="read_layout",
        required=True,
        choices=("se", "pe"),
        help="Effective sequencing layout: se or pe",
    )
    parser.add_argument("--krakentools_script", required=True, help="Path to extract_kraken_reads.py")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--top_n", type=int, default=50, help="Top N taxa to extract. Use 0 for all. Default: 50")
    parser.add_argument("--min_abundance", type=float, default=0.0, help="Extract only from samples with abundance > this value. Default: 0")
    parser.add_argument("--include_children", action="store_true", help="Include child taxa from Kraken report")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing extracted read files")

    args = parser.parse_args()

    ranked_table = Path(args.ranked_table)
    bracken_table = Path(args.bracken_table)
    samplesheet = Path(args.samplesheet)
    temp_dir = Path(args.temp_dir)
    krakentools_script = Path(args.krakentools_script)
    output_dir = Path(args.output_dir)

    for label, path in [
        ("ranked_table", ranked_table),
        ("bracken_table", bracken_table),
        ("samplesheet", samplesheet),
        ("temp_dir", temp_dir),
        ("krakentools_script", krakentools_script),
    ]:
        if label == "temp_dir":
            if not path.is_dir():
                raise SystemExit(f"[ERROR] {label} not found or not a directory: {path}")
        else:
            if not path.exists() or path.stat().st_size == 0:
                raise SystemExit(f"[ERROR] {label} not found or empty: {path}")

    output_dir.mkdir(parents=True, exist_ok=True)

    samples = read_samplesheet_samples(samplesheet)
    taxon_to_taxid = build_taxon_to_taxid(bracken_table)

    ranked_header, ranked_rows = read_table(ranked_table)

    if "taxon" not in ranked_header:
        raise SystemExit(
            "[ERROR] ranked_table must contain a 'taxon' column.\n"
            f"Header: {ranked_header}"
        )

    sample_cols = get_sample_columns(ranked_header, samples)

    if not sample_cols:
        raise SystemExit("[ERROR] Could not detect sample columns in ranked_table.")

    selected_rows = select_ranked_rows(ranked_rows, args.top_n)

    available_cpus = os.cpu_count() or 1

    parallel_jobs = max(
        1,
        min(
            MAX_PARALLEL_EXTRACTION_JOBS,
            len(sample_cols),
            available_cpus,
        ),
    )

    planned_operations = sum(
        1
        for selected_row in selected_rows
        for sample in sample_cols
        if to_float(
            selected_row.get(sample, "0"),
            default=0.0,
        ) > args.min_abundance
    )

    summary_file = output_dir / "extraction_summary.tsv"

    summary_fields = [
        "rank_position",
        "taxon",
        "taxid",
        "sample",
        "abundance",
        "seq_format",
        "output_file",
        "extracted_reads",
        "status",
        "message",
        "layout",
        "output_file_r2",
        "extracted_reads_r2",
        "pair_count_match",
    ]

    print("============================================================")
    print("MTD extract reads by detected microbiome")
    print("Ranked table:", ranked_table)
    print("Bracken table:", bracken_table)
    print("Samplesheet:", samplesheet)
    print("Temp dir:", temp_dir)
    print("Output dir:", output_dir)
    print("Read layout:", args.read_layout)
    print("Top N:", args.top_n)
    print("Min abundance:", args.min_abundance)
    print("Include children:", args.include_children)
    print("Samples detected:", ", ".join(samples))
    print("Sample columns used:", ", ".join(sample_cols))
    print("Selected taxa:", len(selected_rows))
    print("Parallel sample jobs:", parallel_jobs)
    print("Planned taxon-sample operations:", planned_operations)
    print("============================================================")

    all_summary_rows = []

    for index, row in enumerate(selected_rows, start=1):
        taxon = clean_text(row.get("taxon", ""))
        if not taxon:
            continue

        taxid = taxon_to_taxid.get(taxon, "")

        rank_position = clean_text(row.get("rank_position", ""))
        if not rank_position:
            rank_position = str(index)

        if not taxid:
            all_summary_rows.append({
                "rank_position": rank_position,
                "taxon": taxon,
                "taxid": "NA",
                "sample": "NA",
                "abundance": "NA",
                "seq_format": "NA",
                "output_file": "NA",
                "extracted_reads": "0",
                "status": "taxid_not_found",
                "message": "Taxon not found in Bracken taxon-to-taxid mapping",
            })
            continue

        taxon_safe = safe_name(taxon)
        taxon_dir = output_dir / f"{int(index):06d}_taxid{taxid}_{taxon_safe}"
        taxon_dir.mkdir(parents=True, exist_ok=True)

        sample_output_files_r1 = []
        sample_output_files_r2 = []
        sample_output_format = None

        sample_work_items = []

        # Prepare and validate all sample-level work first. Immediate
        # results such as missing inputs or existing outputs are stored
        # together with extraction tasks so the final summary retains
        # the original samplesheet order.
        for sample in sample_cols:
            abundance = to_float(
                row.get(sample, "0"),
                default=0.0,
            )

            if abundance <= args.min_abundance:
                continue

            kraken_file = (
                temp_dir
                / f"Report_non-host_{sample}.kraken"
            )

            report_file = (
                temp_dir
                / f"Report_non-host_{sample}.txt"
            )

            if args.read_layout == "pe":
                seq_file1 = (
                    temp_dir
                    / f"{sample}_non-host_1.fq"
                )

                seq_file2 = (
                    temp_dir
                    / f"{sample}_non-host_2.fq"
                )
            else:
                seq_file1 = (
                    temp_dir
                    / f"{sample}_non-host.fq"
                )

                seq_file2 = None

            required_paths = [
                kraken_file,
                report_file,
                seq_file1,
            ]

            if seq_file2 is not None:
                required_paths.append(
                    seq_file2
                )

            missing = [
                str(required_path)
                for required_path in required_paths
                if (
                    not required_path.exists()
                    or required_path.stat().st_size == 0
                )
            ]

            if missing:
                sample_work_items.append(
                    {
                        "summary_row": {
                            "rank_position": rank_position,
                            "taxon": taxon,
                            "taxid": taxid,
                            "sample": sample,
                            "abundance": abundance,
                            "seq_format": "NA",
                            "output_file": "NA",
                            "extracted_reads": "0",
                            "status": "missing_input",
                            "message": (
                                "Missing: "
                                + ";".join(missing)
                            ),
                            "layout": args.read_layout,
                            "output_file_r2": "NA",
                            "extracted_reads_r2": "0",
                            "pair_count_match": "NA",
                        },
                        "usable_output": False,
                    }
                )

                continue

            seq_format = detect_sequence_format(
                seq_file1
            )

            if seq_format == "unknown":
                sample_work_items.append(
                    {
                        "summary_row": {
                            "rank_position": rank_position,
                            "taxon": taxon,
                            "taxid": taxid,
                            "sample": sample,
                            "abundance": abundance,
                            "seq_format": seq_format,
                            "output_file": "NA",
                            "extracted_reads": "0",
                            "status": (
                                "unknown_sequence_format"
                            ),
                            "message": (
                                "Could not detect "
                                "FASTA/FASTQ format: "
                                f"{seq_file1}"
                            ),
                            "layout": args.read_layout,
                            "output_file_r2": "NA",
                            "extracted_reads_r2": "0",
                            "pair_count_match": "NA",
                        },
                        "usable_output": False,
                    }
                )

                continue

            if args.read_layout == "pe":
                seq_format2 = detect_sequence_format(
                    seq_file2
                )

                if seq_format2 != seq_format:
                    sample_work_items.append(
                        {
                            "summary_row": {
                                "rank_position": (
                                    rank_position
                                ),
                                "taxon": taxon,
                                "taxid": taxid,
                                "sample": sample,
                                "abundance": abundance,
                                "seq_format": seq_format,
                                "output_file": "NA",
                                "extracted_reads": "0",
                                "status": (
                                    "mate_format_mismatch"
                                ),
                                "message": (
                                    f"R1 format is {seq_format}; "
                                    f"R2 format is {seq_format2}"
                                ),
                                "layout": args.read_layout,
                                "output_file_r2": "NA",
                                "extracted_reads_r2": "0",
                                "pair_count_match": "NA",
                            },
                            "usable_output": False,
                        }
                    )

                    continue

            ext = (
                "fastq"
                if seq_format == "fastq"
                else "fasta"
            )

            output_prefix = (
                f"{sample}.taxid{taxid}."
                f"{taxon_safe}"
            )

            if args.read_layout == "pe":
                out_file1 = (
                    taxon_dir
                    / f"{output_prefix}.R1.{ext}"
                )

                out_file2 = (
                    taxon_dir
                    / f"{output_prefix}.R2.{ext}"
                )
            else:
                out_file1 = (
                    taxon_dir
                    / f"{output_prefix}.{ext}"
                )

                out_file2 = None

            output_exists = (
                out_file1.exists()
                and out_file1.stat().st_size > 0
            )

            if args.read_layout == "pe":
                output_exists = (
                    output_exists
                    and out_file2.exists()
                    and out_file2.stat().st_size > 0
                )

            if output_exists and not args.overwrite:
                extracted_reads1 = count_sequences(
                    out_file1,
                    seq_format,
                )

                if args.read_layout == "pe":
                    extracted_reads2 = count_sequences(
                        out_file2,
                        seq_format,
                    )

                    pair_count_match = (
                        "yes"
                        if extracted_reads1
                        == extracted_reads2
                        else "no"
                    )
                else:
                    extracted_reads2 = 0
                    pair_count_match = (
                        "not_applicable"
                    )

                if pair_count_match == "no":
                    status = "pair_count_mismatch"

                    message = (
                        "Existing paired outputs have "
                        "different sequence counts."
                    )
                else:
                    status = "already_exists"

                    message = (
                        "Existing file kept. "
                        "Use --overwrite to regenerate."
                    )

                sample_work_items.append(
                    {
                        "summary_row": {
                            "rank_position": rank_position,
                            "taxon": taxon,
                            "taxid": taxid,
                            "sample": sample,
                            "abundance": abundance,
                            "seq_format": seq_format,
                            "output_file": str(
                                out_file1
                            ),
                            "extracted_reads": str(
                                extracted_reads1
                            ),
                            "status": status,
                            "message": message,
                            "layout": args.read_layout,
                            "output_file_r2": (
                                str(out_file2)
                                if out_file2 is not None
                                else "NA"
                            ),
                            "extracted_reads_r2": str(
                                extracted_reads2
                            ),
                            "pair_count_match": (
                                pair_count_match
                            ),
                        },
                        "usable_output": (
                            extracted_reads1 > 0
                            and pair_count_match != "no"
                        ),
                        "out_file1": out_file1,
                        "out_file2": out_file2,
                        "seq_format": seq_format,
                    }
                )

                continue

            sample_work_items.append(
                {
                    "task": {
                        "sample": sample,
                        "abundance": abundance,
                        "seq_format": seq_format,
                        "seq_file1": seq_file1,
                        "seq_file2": seq_file2,
                        "kraken_file": kraken_file,
                        "report_file": report_file,
                        "out_file1": out_file1,
                        "out_file2": out_file2,
                    }
                }
            )

        extraction_item_indices = [
            item_index
            for item_index, item in enumerate(
                sample_work_items
            )
            if "task" in item
        ]

        future_by_item_index = {}

        if extraction_item_indices:
            taxon_jobs = min(
                parallel_jobs,
                len(extraction_item_indices),
            )

            print(
                "[EXTRACT TAXON]",
                f"rank={rank_position}",
                f"taxid={taxid}",
                f"taxon={taxon}",
                (
                    "sample_jobs="
                    f"{len(extraction_item_indices)}"
                ),
                f"parallel_jobs={taxon_jobs}",
                flush=True,
            )

            with ThreadPoolExecutor(
                max_workers=taxon_jobs
            ) as executor:

                # Launch every eligible sample for this taxon
                # before waiting for results.
                for item_index in (
                    extraction_item_indices
                ):
                    task = sample_work_items[
                        item_index
                    ]["task"]

                    future_by_item_index[
                        item_index
                    ] = executor.submit(
                        extract_reads_for_taxon_sample,
                        krakentools_script=(
                            krakentools_script
                        ),
                        kraken_file=(
                            task["kraken_file"]
                        ),
                        seq_file1=task["seq_file1"],
                        seq_file2=task["seq_file2"],
                        report_file=(
                            task["report_file"]
                        ),
                        taxid=taxid,
                        out_file1=task["out_file1"],
                        out_file2=task["out_file2"],
                        include_children=(
                            args.include_children
                        ),
                        seq_format=(
                            task["seq_format"]
                        ),
                        read_layout=args.read_layout,
                    )

                completed_for_taxon = 0

                # Collect results in original sample order.
                # The subprocesses still execute concurrently.
                for item_index, item in enumerate(
                    sample_work_items
                ):
                    if "summary_row" in item:
                        all_summary_rows.append(
                            item["summary_row"]
                        )

                        if item.get(
                            "usable_output",
                            False,
                        ):
                            sample_output_files_r1.append(
                                item["out_file1"]
                            )

                            if args.read_layout == "pe":
                                sample_output_files_r2.append(
                                    item["out_file2"]
                                )

                            sample_output_format = (
                                item["seq_format"]
                            )

                        continue

                    task = item["task"]

                    completed_for_taxon += 1

                    try:
                        result, cmd = (
                            future_by_item_index[
                                item_index
                            ].result()
                        )

                        extraction_exception = None
                    except Exception as exc:
                        result = None
                        cmd = []
                        extraction_exception = exc

                    extracted_reads1 = count_sequences(
                        task["out_file1"],
                        task["seq_format"],
                    )

                    if args.read_layout == "pe":
                        extracted_reads2 = (
                            count_sequences(
                                task["out_file2"],
                                task["seq_format"],
                            )
                        )

                        pair_count_match = (
                            "yes"
                            if extracted_reads1
                            == extracted_reads2
                            else "no"
                        )
                    else:
                        extracted_reads2 = 0
                        pair_count_match = (
                            "not_applicable"
                        )

                    if extraction_exception is not None:
                        status = "failed"

                        message = (
                            "Extraction worker failed: "
                            f"{extraction_exception}"
                        )
                    elif result.returncode == 0:
                        if pair_count_match == "no":
                            status = (
                                "pair_count_mismatch"
                            )

                            message = (
                                "KrakenTools completed, "
                                "but R1 and R2 counts "
                                "differ: "
                                f"R1={extracted_reads1}; "
                                f"R2={extracted_reads2}"
                            )
                        else:
                            status = "ok"

                            if args.read_layout == "pe":
                                message = (
                                    "Extracted paired "
                                    "reads: "
                                    f"{extracted_reads1} "
                                    "pairs"
                                )
                            else:
                                message = (
                                    "Extracted reads: "
                                    f"{extracted_reads1}"
                                )
                    else:
                        status = "failed"

                        message = (
                            result.stderr
                            or result.stdout
                            or (
                                "extract_kraken_reads.py "
                                "failed"
                            )
                        ).replace(
                            "\n",
                            " | ",
                        )

                    if (
                        extracted_reads1 > 0
                        and status == "ok"
                    ):
                        sample_output_files_r1.append(
                            task["out_file1"]
                        )

                        if args.read_layout == "pe":
                            sample_output_files_r2.append(
                                task["out_file2"]
                            )

                        sample_output_format = (
                            task["seq_format"]
                        )

                    all_summary_rows.append(
                        {
                            "rank_position": (
                                rank_position
                            ),
                            "taxon": taxon,
                            "taxid": taxid,
                            "sample": task["sample"],
                            "abundance": (
                                task["abundance"]
                            ),
                            "seq_format": (
                                task["seq_format"]
                            ),
                            "output_file": str(
                                task["out_file1"]
                            ),
                            "extracted_reads": str(
                                extracted_reads1
                            ),
                            "status": status,
                            "message": message,
                            "layout": args.read_layout,
                            "output_file_r2": (
                                str(task["out_file2"])
                                if (
                                    task["out_file2"]
                                    is not None
                                )
                                else "NA"
                            ),
                            "extracted_reads_r2": str(
                                extracted_reads2
                            ),
                            "pair_count_match": (
                                pair_count_match
                            ),
                        }
                    )

                    print(
                        "[EXTRACT DONE]",
                        (
                            f"taxon={index}/"
                            f"{len(selected_rows)}"
                        ),
                        (
                            "sample_job="
                            f"{completed_for_taxon}/"
                            f"{len(extraction_item_indices)}"
                        ),
                        f"sample={task['sample']}",
                        f"taxid={taxid}",
                        f"status={status}",
                        flush=True,
                    )
        else:
            # No new subprocess was required, but existing or
            # skipped results must still be added to the summary
            # and combined-output lists.
            for item in sample_work_items:
                all_summary_rows.append(
                    item["summary_row"]
                )

                if item.get(
                    "usable_output",
                    False,
                ):
                    sample_output_files_r1.append(
                        item["out_file1"]
                    )

                    if args.read_layout == "pe":
                        sample_output_files_r2.append(
                            item["out_file2"]
                        )

                    sample_output_format = (
                        item["seq_format"]
                    )

        # Concatenate sample-level outputs into one file per taxon.
        if sample_output_files_r1:
            ext = (
                "fastq"
                if sample_output_format == "fastq"
                else "fasta"
            )

            if args.read_layout == "pe":
                combined_file1 = (
                    taxon_dir
                    / f"all_samples.taxid{taxid}.{taxon_safe}.R1.{ext}"
                )

                combined_file2 = (
                    taxon_dir
                    / f"all_samples.taxid{taxid}.{taxon_safe}.R2.{ext}"
                )

                concatenate_files(
                    sample_output_files_r1,
                    combined_file1,
                )

                concatenate_files(
                    sample_output_files_r2,
                    combined_file2,
                )

                combined_count1 = count_sequences(
                    combined_file1,
                    sample_output_format,
                )

                combined_count2 = count_sequences(
                    combined_file2,
                    sample_output_format,
                )

                pair_count_match = (
                    "yes"
                    if combined_count1 == combined_count2
                    else "no"
                )

                all_summary_rows.append({
                    "rank_position": rank_position,
                    "taxon": taxon,
                    "taxid": taxid,
                    "sample": "all_samples",
                    "abundance": "NA",
                    "seq_format": sample_output_format,
                    "output_file": str(combined_file1),
                    "extracted_reads": str(combined_count1),
                    "status": (
                        "combined"
                        if pair_count_match == "yes"
                        else "combined_pair_count_mismatch"
                    ),
                    "message": (
                        "Concatenated paired outputs from "
                        f"{len(sample_output_files_r1)} samples"
                    ),
                    "layout": args.read_layout,
                    "output_file_r2": str(combined_file2),
                    "extracted_reads_r2": str(combined_count2),
                    "pair_count_match": pair_count_match,
                })

            else:
                combined_file = (
                    taxon_dir
                    / f"all_samples.taxid{taxid}.{taxon_safe}.{ext}"
                )

                concatenate_files(
                    sample_output_files_r1,
                    combined_file,
                )

                combined_count = count_sequences(
                    combined_file,
                    sample_output_format,
                )

                all_summary_rows.append({
                    "rank_position": rank_position,
                    "taxon": taxon,
                    "taxid": taxid,
                    "sample": "all_samples",
                    "abundance": "NA",
                    "seq_format": sample_output_format,
                    "output_file": str(combined_file),
                    "extracted_reads": str(combined_count),
                    "status": "combined",
                    "message": (
                        f"Concatenated "
                        f"{len(sample_output_files_r1)} sample files"
                    ),
                    "layout": args.read_layout,
                    "output_file_r2": "NA",
                    "extracted_reads_r2": "0",
                    "pair_count_match": "not_applicable",
                })

    with summary_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for row in all_summary_rows:
            writer.writerow(row)

    print("============================================================")
    print("[DONE] Extraction finished")
    print("Summary:", summary_file)
    print("Rows:", len(all_summary_rows))
    print("============================================================")


if __name__ == "__main__":
    main()
