#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ACCESSION_VERSION_RE = re.compile(r"\.\d+$")
TOKEN_RE = re.compile(r"[A-Za-z0-9_]+(?:\.\d+)?")


def normalize(value: str) -> str:
    return ACCESSION_VERSION_RE.sub(
        "",
        value.strip().lstrip(">"),
    )


def read_targets(path: Path) -> set[str]:
    if not path.is_file():
        raise SystemExit(
            f"[ERROR] Unmapped-accession file not found: {path}"
        )

    targets = {
        normalize(line)
        for line in path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()
        if line.strip()
    }

    if not targets:
        raise SystemExit(
            "[ERROR] No unmapped accessions were supplied."
        )

    return targets


def validate_origins(
    path: Path,
    targets: set[str],
) -> None:
    if not path.is_file():
        raise SystemExit(
            f"[ERROR] Origin report not found: {path}"
        )

    origins: dict[str, set[str]] = defaultdict(set)

    with path.open(
        encoding="utf-8",
        errors="replace",
    ) as handle:
        next(handle, None)

        for line in handle:
            fields = line.rstrip("\n").split("\t")

            if not fields or not fields[0].strip():
                continue

            accession = normalize(fields[0])

            if accession not in targets:
                continue

            if len(fields) > 1 and fields[1].strip():
                origins[accession].update(
                    item.strip()
                    for item in fields[1].split(",")
                    if item.strip()
                )

    invalid = {
        accession: origins.get(accession, set())
        for accession in targets
        if origins.get(accession, set()) != {"plasmid"}
    }

    if invalid:
        details = ", ".join(
            (
                f"{accession}="
                f"{','.join(sorted(libraries)) or 'unknown'}"
            )
            for accession, libraries
            in sorted(invalid.items())[:20]
        )
        raise SystemExit(
            "[ERROR] Automatic filtering is permitted only when "
            "every unmapped accession belongs exclusively to the "
            f"plasmid library. Invalid origins: {details}"
        )


def atomic_temp(path: Path) -> tuple[int, Path]:
    descriptor, name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".mtd-filter.tmp",
        dir=str(path.parent),
    )
    return descriptor, Path(name)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Remove a small, audited set of unmappable plasmid "
            "records from a generated Kraken2 library before "
            "database construction."
        )
    )
    parser.add_argument(
        "--database",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--unmapped",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--origins",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--max-count",
        type=int,
        default=1000,
    )
    parser.add_argument(
        "--max-fraction",
        type=float,
        default=0.01,
    )
    args = parser.parse_args()

    database = args.database.expanduser().resolve()
    unmapped_path = args.unmapped.expanduser().resolve()
    origins_path = args.origins.expanduser().resolve()
    output_path = args.output.expanduser().resolve()

    if args.max_count < 1:
        raise SystemExit(
            "[ERROR] --max-count must be at least 1."
        )

    if not 0 < args.max_fraction <= 1:
        raise SystemExit(
            "[ERROR] --max-fraction must be within (0, 1]."
        )

    targets = read_targets(unmapped_path)
    validate_origins(origins_path, targets)

    if len(targets) > args.max_count:
        raise SystemExit(
            "[ERROR] Refusing to filter "
            f"{len(targets):,} plasmid accessions; maximum is "
            f"{args.max_count:,}."
        )

    plasmid_dir = database / "library" / "plasmid"
    prelim_path = plasmid_dir / "prelim_map.txt"
    fasta_path = plasmid_dir / "library.fna"

    for path in (prelim_path, fasta_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(
                f"[ERROR] Required plasmid file is missing: {path}"
            )

    descriptor, prelim_temp = atomic_temp(prelim_path)
    os.close(descriptor)
    descriptor, fasta_temp = atomic_temp(fasta_path)
    os.close(descriptor)

    target_to_seqids: dict[str, set[str]] = defaultdict(set)
    seqid_to_targets: dict[str, set[str]] = defaultdict(set)
    total_accnum = 0
    removed_prelim_lines = 0
    removed_records = 0
    kept_records = 0
    removed_details: dict[
        str,
        set[tuple[str, str]],
    ] = defaultdict(set)

    try:
        with prelim_path.open(
            encoding="utf-8",
            errors="replace",
        ) as source, prelim_temp.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as destination:
            for line in source:
                fields = line.rstrip("\n").split("\t")

                if fields and fields[0] == "ACCNUM":
                    total_accnum += 1
                    matched = {
                        normalize(field)
                        for field in fields[1:]
                        if normalize(field) in targets
                    }

                    if matched:
                        if len(fields) < 2 or not fields[1].strip():
                            raise SystemExit(
                                "[ERROR] An unmapped ACCNUM row "
                                "has no sequence identifier."
                            )

                        seqid = fields[1].strip()

                        for accession in matched:
                            target_to_seqids[accession].add(
                                seqid
                            )
                            seqid_to_targets[seqid].add(
                                accession
                            )

                        removed_prelim_lines += 1
                        continue

                destination.write(line)

        if total_accnum == 0:
            raise SystemExit(
                "[ERROR] No ACCNUM rows were found in the "
                "plasmid prelim map."
            )

        fraction = len(targets) / total_accnum

        if fraction > args.max_fraction:
            raise SystemExit(
                "[ERROR] Refusing to filter "
                f"{len(targets):,}/{total_accnum:,} plasmid "
                f"accessions ({fraction:.6%}); maximum fraction "
                f"is {args.max_fraction:.6%}."
            )

        missing_prelim = targets - set(target_to_seqids)

        if missing_prelim:
            preview = ", ".join(
                sorted(missing_prelim)[:20]
            )
            raise SystemExit(
                "[ERROR] Some unmapped accessions were not found "
                f"in the plasmid prelim map: {preview}"
            )

        current_remove = False

        with fasta_path.open(
            "r",
            encoding="utf-8",
            errors="replace",
            newline="",
        ) as source, fasta_temp.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as destination:
            for line in source:
                if line.startswith(">"):
                    header = line[1:].rstrip("\r\n")
                    seqid = (
                        header.split(None, 1)[0]
                        if header
                        else ""
                    )

                    token_targets: dict[str, str] = {}

                    for token in TOKEN_RE.findall(header):
                        base = normalize(token)

                        if (
                            base in targets
                            and base not in token_targets
                        ):
                            token_targets[base] = token

                    current_targets = set(
                        seqid_to_targets.get(seqid, set())
                    )
                    current_targets.update(token_targets)
                    current_remove = bool(current_targets)

                    if current_remove:
                        removed_records += 1

                        for accession in current_targets:
                            versioned = token_targets.get(
                                accession,
                                accession,
                            )
                            removed_details[accession].add(
                                (seqid, versioned)
                            )
                    else:
                        kept_records += 1
                        destination.write(line)

                    continue

                if not current_remove:
                    destination.write(line)

        if kept_records == 0:
            raise SystemExit(
                "[ERROR] Filtering would remove every plasmid "
                "FASTA record."
            )

        missing_fasta = targets - set(removed_details)

        if missing_fasta:
            preview = ", ".join(
                sorted(missing_fasta)[:20]
            )
            raise SystemExit(
                "[ERROR] Some unmapped accessions were not found "
                f"in the plasmid FASTA: {preview}"
            )

        if (
            removed_records == 0
            or removed_prelim_lines == 0
        ):
            raise SystemExit(
                "[ERROR] No plasmid records were selected "
                "for filtering."
            )

        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        descriptor, report_temp = atomic_temp(
            output_path
        )
        os.close(descriptor)

        try:
            with report_temp.open(
                "w",
                encoding="utf-8",
                newline="",
            ) as report:
                report.write(
                    "accession\taccession_version\tsequence_id\t"
                    "library\treason\tfiltered_at_utc\n"
                )
                timestamp = datetime.now(
                    timezone.utc
                ).strftime(
                    "%Y-%m-%dT%H:%M:%SZ"
                )

                for accession in sorted(targets):
                    for seqid, versioned in sorted(
                        removed_details[accession]
                    ):
                        report.write(
                            f"{accession}\t{versioned}\t"
                            f"{seqid}\tplasmid\t"
                            "absent_from_current_"
                            "accession2taxid_maps\t"
                            f"{timestamp}\n"
                        )

            # Both complete temporary outputs are created before
            # replacing either original generated-library file.
            os.replace(fasta_temp, fasta_path)
            os.replace(prelim_temp, prelim_path)
            os.replace(report_temp, output_path)
        finally:
            report_temp.unlink(missing_ok=True)

    finally:
        prelim_temp.unlink(missing_ok=True)
        fasta_temp.unlink(missing_ok=True)

    print(
        "[MTD-KRAKEN2] Filtered unmappable plasmid "
        f"accessions: {len(targets):,}"
    )
    print(
        "[MTD-KRAKEN2] Removed plasmid FASTA records: "
        f"{removed_records:,}"
    )
    print(
        "[MTD-KRAKEN2] Removed plasmid prelim-map rows: "
        f"{removed_prelim_lines:,}"
    )
    print(
        "[MTD-KRAKEN2] Filtered fraction: "
        f"{len(targets) / total_accnum:.6%}"
    )
    print(
        f"[MTD-KRAKEN2] Audit report: {output_path}"
    )


if __name__ == "__main__":
    main()
