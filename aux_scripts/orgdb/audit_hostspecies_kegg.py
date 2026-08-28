#!/usr/bin/env python3
"""
Audit the KEGG organism-code column in MTD Explorer HostSpecies.csv.

The audit is read-only: it never modifies HostSpecies.csv.

Official KEGG sources
---------------------
- /list/genome
    KEGG genome ID -> organism code/name
- /link/taxonomy/genome
    KEGG genome ID -> exact NCBI Taxonomy ID
- /link/taxonomy/genome/species
    KEGG genome ID -> species-level NCBI Taxonomy ID

The species-level mapping is important for entries where HostSpecies.csv stores a
species TaxID but KEGG annotates a subspecies/strain genome.

Only Python 3 standard-library modules are required.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import socket
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence


SCRIPT_VERSION = "1.0.3"
KEGG_REST = "https://rest.kegg.jp"

REQUIRED_COLUMNS = {
    "Taxon_ID",
    "Scientific_name",
    "Common_name",
    "kegg",
}

REPORT_COLUMNS = (
    "CSV_line",
    "Taxon_ID",
    "Scientific_name",
    "Reference_Taxon_ID",
    "Reference_Scientific_name",
    "Common_name",
    "kegg_csv",
    "kegg_exists",
    "kegg_genome_id",
    "kegg_name",
    "kegg_taxid",
    "kegg_species_taxid",
    "status",
    "severity",
    "candidate_basis",
    "candidate_codes",
    "suggested_kegg",
    "suggested_action",
    "notes",
)

OK_STATUSES = {
    "OK_EXACT_TAXID",
    "OK_REFERENCE_TAXID",
    "OK_WITHIN_SPECIES",
    "OK_REFERENCE_SPECIES",
}

PROBLEM_STATUSES = {
    "WRONG_TAXON",
    "INVALID_OR_RETIRED",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Audit HostSpecies.csv KEGG organism codes against official KEGG "
            "genome-to-NCBI-taxonomy mappings. The input CSV is never modified."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--hostspecies",
        default="HostSpecies.csv",
        help="Input HostSpecies.csv.",
    )
    parser.add_argument(
        "-o",
        "--report",
        default="HostSpecies.kegg_audit.tsv",
        help="Output TSV audit report.",
    )
    parser.add_argument(
        "--cache-dir",
        default=str(Path.home() / ".cache" / "mtd_explorer" / "kegg_audit"),
        help="Directory used to cache KEGG responses.",
    )
    parser.add_argument(
        "--cache-max-age-days",
        type=int,
        default=7,
        help="Reuse cached KEGG responses up to this age; use 0 to refresh.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=90,
        help="Network timeout per request, in seconds.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=5,
        help="Number of attempts for each KEGG request.",
    )
    parser.add_argument(
        "--request-delay",
        type=float,
        default=0.40,
        help="Delay between uncached KEGG requests (KEGG limit: <=3 requests/s).",
    )
    parser.add_argument(
        "--fail-on-problems",
        action="store_true",
        help=(
            "Exit with status 2 when WRONG_TAXON or INVALID_OR_RETIRED entries "
            "are detected. Useful for CI/checks."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=SCRIPT_VERSION,
    )

    args = parser.parse_args()
    if args.cache_max_age_days < 0:
        parser.error("--cache-max-age-days cannot be negative")
    if args.timeout < 1:
        parser.error("--timeout must be at least 1")
    if args.retries < 1:
        parser.error("--retries must be at least 1")
    if args.request_delay < 0:
        parser.error("--request-delay cannot be negative")
    return args


def say(message: str) -> None:
    print(message, flush=True)


def norm_spaces(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip())


def norm_name(value: Any) -> str:
    text = norm_spaces(value).lower().replace("_", " ")
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return norm_spaces(text)


def kegg_base_name(value: str) -> str:
    text = norm_spaces(value)
    return norm_spaces(re.split(r"\s+\(", text, maxsplit=1)[0])


def kegg_common_qualifier(value: str) -> str:
    """Return a final parenthetical KEGG common-name qualifier, if present."""
    text = norm_spaces(value)
    match = re.search(r"\(([^()]*)\)\s*$", text)
    if not match:
        return ""
    return norm_spaces(match.group(1))


def cache_is_fresh(path: Path, max_age_days: int) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    if max_age_days == 0:
        return False
    return (time.time() - path.stat().st_mtime) <= max_age_days * 86400


class Downloader:
    def __init__(
        self,
        cache_dir: Path,
        timeout: int,
        retries: int,
        request_delay: float,
        cache_max_age_days: int,
    ) -> None:
        self.cache_dir = cache_dir
        self.timeout = timeout
        self.retries = retries
        self.request_delay = request_delay
        self.cache_max_age_days = cache_max_age_days
        self.user_agent = f"MTD-Explorer-KEGG-audit/{SCRIPT_VERSION}"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get_text(self, url: str, cache_name: str) -> str:
        cache_path = self.cache_dir / cache_name
        if cache_is_fresh(cache_path, self.cache_max_age_days):
            return cache_path.read_text(encoding="utf-8", errors="replace")

        headers = {
            "User-Agent": self.user_agent,
            "Accept": "text/plain",
        }
        last_error: Optional[BaseException] = None

        for attempt in range(1, self.retries + 1):
            try:
                request = urllib.request.Request(url, headers=headers)
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    payload = response.read()
                if not payload:
                    raise RuntimeError("server returned an empty response")

                temporary = cache_path.with_suffix(cache_path.suffix + ".tmp")
                temporary.write_bytes(payload)
                os.replace(temporary, cache_path)

                if self.request_delay:
                    time.sleep(self.request_delay)

                return payload.decode("utf-8", errors="replace")

            except urllib.error.HTTPError as exc:
                last_error = exc
                if exc.code in {400, 401, 403, 404, 405, 410}:
                    break
            except (
                urllib.error.URLError,
                TimeoutError,
                socket.timeout,
                OSError,
                RuntimeError,
            ) as exc:
                last_error = exc

            if attempt < self.retries:
                time.sleep(min(30, 2 ** (attempt - 1)))

        if cache_path.is_file() and cache_path.stat().st_size > 0:
            say(f"[WARN] KEGG refresh failed for {url}; using stale cache {cache_path}")
            return cache_path.read_text(encoding="utf-8", errors="replace")

        raise RuntimeError(f"Unable to retrieve {url}: {last_error}")


def read_hostspecies(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise RuntimeError(f"HostSpecies.csv not found: {path}")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS.difference(fieldnames)
        if missing:
            raise RuntimeError(
                "HostSpecies.csv is missing required column(s): "
                + ", ".join(sorted(missing))
            )
        rows = [dict(row) for row in reader]

    seen_taxids: set[str] = set()
    duplicates: list[str] = []
    for line_number, row in enumerate(rows, start=2):
        taxid = norm_spaces(row.get("Taxon_ID"))
        scientific_name = norm_spaces(row.get("Scientific_name"))
        if not taxid or not taxid.isdigit():
            raise RuntimeError(
                f"Invalid Taxon_ID {taxid!r} on CSV line {line_number}"
            )
        if not scientific_name:
            raise RuntimeError(f"Missing Scientific_name on CSV line {line_number}")
        if taxid in seen_taxids:
            duplicates.append(taxid)
        seen_taxids.add(taxid)

    if duplicates:
        raise RuntimeError(
            "Duplicate Taxon_ID values in HostSpecies.csv: "
            + ", ".join(sorted(set(duplicates), key=int))
        )

    return fieldnames, rows


def parse_organism_list(text: str) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    """Parse the current /list/genome format and the legacy /list/organism format.

    Current KEGG /list/genome output (June 2026+) is two columns, e.g.::

        T01001\thsa; Homo sapiens (human)

    Older /list/organism output used separate columns for genome ID, organism
    code, name, and lineage. Supporting both formats keeps cached/older data
    usable while avoiding dependence on the currently failing /list/organism
    endpoint.
    """
    by_code: dict[str, dict[str, str]] = {}
    by_genome: dict[str, dict[str, str]] = {}

    for raw_line in text.splitlines():
        if not raw_line.strip():
            continue
        fields = raw_line.rstrip("\n").split("\t")
        genome_id = fields[0].strip() if fields else ""
        code = ""
        name = ""
        lineage = ""

        if len(fields) >= 3:
            # Legacy /list/organism format.
            code = fields[1].strip().lower()
            name = fields[2].strip()
            lineage = fields[3].strip() if len(fields) >= 4 else ""
        elif len(fields) == 2:
            # Current /list/genome format:
            #   T01001<TAB>hsa; Homo sapiens (human)
            definition = fields[1].strip()
            match = re.match(r"^([A-Za-z][A-Za-z0-9]{1,4})\s*;\s*(.+)$", definition)
            if not match:
                continue
            code = match.group(1).lower()
            name = match.group(2).strip()
        else:
            continue

        if not re.fullmatch(r"T\d+", genome_id):
            continue
        if not re.fullmatch(r"[a-z][a-z0-9]{1,4}", code):
            continue
        if not name:
            continue

        entry = {
            "genome_id": genome_id,
            "code": code,
            "name": name,
            "lineage": lineage,
        }
        by_code[code] = entry
        by_genome[genome_id] = entry

    if not by_code:
        raise RuntimeError("KEGG genome catalog returned no usable entries")
    return by_code, by_genome


def parse_genome_taxonomy_links(text: str) -> dict[str, str]:
    """Parse KEGG genome <-> NCBI Taxonomy LINK output.

    Current KEGG output (August 2026) identifies organisms by organism code,
    for example::

        gn:hsa\ttaxid:9606
        gn:ptr\ttaxid:9598

    Older/alternate output may identify the genome with its T-number instead.
    Accept both representations, and also accept either column direction.  The
    returned mapping is therefore keyed by the normalized KEGG organism code
    (``hsa``) or by the T-number (``T01001``), whichever KEGG supplied.
    """
    result: dict[str, str] = {}
    nonempty_lines = [line.strip() for line in text.splitlines() if line.strip()]

    for raw_line in nonempty_lines:
        fields = raw_line.split("\t")
        if len(fields) != 2:
            continue

        genome_field_index: Optional[int] = None
        genome_key = ""
        for index, field in enumerate(fields):
            field = field.strip()

            code_match = re.fullmatch(r"(?:gn:)?([A-Za-z][A-Za-z0-9]{1,4})", field)
            if code_match:
                genome_field_index = index
                genome_key = code_match.group(1).lower()
                break

            tnumber_match = re.fullmatch(r"(?:gn:)?(T\d+)", field, flags=re.IGNORECASE)
            if tnumber_match:
                genome_field_index = index
                genome_key = tnumber_match.group(1).upper()
                break

        if genome_field_index is None:
            continue

        taxonomy_field = fields[1 - genome_field_index].strip()
        taxid_match = re.search(r"(?<!\d)(\d+)(?!\d)", taxonomy_field)
        if not taxid_match:
            continue

        result[genome_key] = taxid_match.group(1)

    if not result:
        preview = " | ".join(repr(line) for line in nonempty_lines[:5])
        if not preview:
            preview = "<empty response>"
        raise RuntimeError(
            "KEGG taxonomy link response returned no usable entries; "
            f"first response line(s): {preview}"
        )
    return result


def enrich_catalog(
    by_code: dict[str, dict[str, str]],
    exact_taxid_by_genome: dict[str, str],
    species_taxid_by_genome: dict[str, str],
) -> tuple[
    dict[str, dict[str, str]],
    dict[str, list[dict[str, str]]],
    dict[str, list[dict[str, str]]],
]:
    exact_index: dict[str, list[dict[str, str]]] = defaultdict(list)
    species_index: dict[str, list[dict[str, str]]] = defaultdict(list)

    for code, entry in by_code.items():
        genome_id = entry["genome_id"]
        # KEGG LINK currently returns gn:<organism_code>, while older output
        # may use the T-number.  Prefer the code and retain the T-number as a
        # compatibility fallback.
        entry["taxid"] = exact_taxid_by_genome.get(
            code,
            exact_taxid_by_genome.get(genome_id, ""),
        )
        entry["species_taxid"] = species_taxid_by_genome.get(
            code,
            species_taxid_by_genome.get(genome_id, ""),
        )

        if entry["taxid"]:
            exact_index[entry["taxid"]].append(entry)
        if entry["species_taxid"]:
            species_index[entry["species_taxid"]].append(entry)

    for values in exact_index.values():
        values.sort(key=lambda item: item["code"])
    for values in species_index.values():
        values.sort(key=lambda item: item["code"])

    return by_code, dict(exact_index), dict(species_index)


def unique_entries(entries: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    seen: set[str] = set()
    result: list[dict[str, str]] = []
    for entry in entries:
        code = entry["code"]
        if code in seen:
            continue
        seen.add(code)
        result.append(entry)
    return result


def candidate_set(
    taxid: str,
    reference_taxid: str,
    exact_index: dict[str, list[dict[str, str]]],
    species_index: dict[str, list[dict[str, str]]],
) -> tuple[list[dict[str, str]], str]:
    """
    Return candidates using a strict priority order.

    Exact TaxID is preferred, followed by the explicitly curated reference
    TaxID, then species-level mappings for the same two IDs.
    """
    levels = (
        ("EXACT_TAXID", exact_index.get(taxid, [])),
        ("REFERENCE_TAXID", exact_index.get(reference_taxid, []) if reference_taxid else []),
        ("SPECIES_TAXID", species_index.get(taxid, [])),
        ("REFERENCE_SPECIES", species_index.get(reference_taxid, []) if reference_taxid else []),
    )
    for basis, entries in levels:
        unique = unique_entries(entries)
        if unique:
            return unique, basis
    return [], ""


def choose_suggestion(
    candidates: Sequence[dict[str, str]],
    scientific_name: str,
    common_name: str,
) -> tuple[str, str]:
    if not candidates:
        return "", ""
    if len(candidates) == 1:
        return candidates[0]["code"], "unique candidate"

    common_norm = norm_name(common_name)
    if common_norm:
        common_matches = [
            entry
            for entry in candidates
            if norm_name(kegg_common_qualifier(entry["name"])) == common_norm
        ]
        if len(common_matches) == 1:
            return common_matches[0]["code"], "unique KEGG common-name match"

    scientific_norm = norm_name(scientific_name)
    scientific_matches = [
        entry
        for entry in candidates
        if norm_name(kegg_base_name(entry["name"])) == scientific_norm
    ]
    if len(scientific_matches) == 1:
        return scientific_matches[0]["code"], "unique scientific-name match"

    return "", "multiple equally plausible candidates"


def classify_existing(
    entry: dict[str, str],
    taxid: str,
    reference_taxid: str,
) -> tuple[str, str, str]:
    exact_taxid = entry.get("taxid", "")
    species_taxid = entry.get("species_taxid", "")

    if exact_taxid and exact_taxid == taxid:
        return "OK_EXACT_TAXID", "OK", "KEGG exact TaxID matches Taxon_ID"
    if reference_taxid and exact_taxid and exact_taxid == reference_taxid:
        return (
            "OK_REFERENCE_TAXID",
            "OK",
            "KEGG exact TaxID matches Reference_Taxon_ID",
        )
    if species_taxid and species_taxid == taxid:
        return (
            "OK_WITHIN_SPECIES",
            "OK",
            "KEGG genome TaxID differs, but its species-level TaxID matches Taxon_ID",
        )
    if reference_taxid and species_taxid and species_taxid == reference_taxid:
        return (
            "OK_REFERENCE_SPECIES",
            "OK",
            "KEGG species-level TaxID matches Reference_Taxon_ID",
        )

    return (
        "WRONG_TAXON",
        "ERROR",
        "KEGG code exists but its exact/species TaxID matches neither Taxon_ID nor Reference_Taxon_ID",
    )


def write_report(path: Path, rows: Sequence[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=REPORT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def audit_rows(
    rows: Sequence[dict[str, str]],
    by_code: dict[str, dict[str, str]],
    exact_index: dict[str, list[dict[str, str]]],
    species_index: dict[str, list[dict[str, str]]],
) -> list[dict[str, str]]:
    reports: list[dict[str, str]] = []

    for csv_line, row in enumerate(rows, start=2):
        taxid = norm_spaces(row.get("Taxon_ID"))
        scientific_name = norm_spaces(row.get("Scientific_name"))
        common_name = norm_spaces(row.get("Common_name"))
        reference_taxid = norm_spaces(row.get("Reference_Taxon_ID"))
        reference_name = norm_spaces(row.get("Reference_Scientific_name"))
        current_code = norm_spaces(row.get("kegg")).lower()

        candidates, candidate_basis = candidate_set(
            taxid,
            reference_taxid,
            exact_index,
            species_index,
        )
        candidate_codes = ",".join(entry["code"] for entry in candidates)
        suggested, suggestion_reason = choose_suggestion(
            candidates,
            scientific_name,
            common_name,
        )

        kegg_entry: Optional[dict[str, str]] = None
        notes: list[str] = []

        if current_code:
            kegg_entry = by_code.get(current_code)
            if kegg_entry is None:
                status = "INVALID_OR_RETIRED"
                severity = "ERROR"
                notes.append("kegg value is absent from current KEGG genome catalog")
            else:
                status, severity, reason = classify_existing(
                    kegg_entry,
                    taxid,
                    reference_taxid,
                )
                notes.append(reason)

            if status in OK_STATUSES:
                # Do not suggest replacing a valid code even if another genome
                # exists for the same TaxID.
                suggested_kegg = current_code
                suggested_action = "KEEP"
            elif suggested:
                suggested_kegg = suggested
                suggested_action = "REPLACE"
                if suggestion_reason:
                    notes.append(suggestion_reason)
            elif status == "WRONG_TAXON":
                suggested_kegg = ""
                suggested_action = "CLEAR_OR_MANUAL_REVIEW"
            else:
                suggested_kegg = ""
                suggested_action = "MANUAL_REVIEW"
        else:
            if not candidates:
                status = "NO_KEGG_ENTRY"
                severity = "INFO"
                suggested_kegg = ""
                suggested_action = "KEEP_EMPTY"
                notes.append("no exact/reference/species-level KEGG candidate found")
            elif suggested:
                if candidate_basis == "EXACT_TAXID":
                    status = "MISSING_CANDIDATE_EXACT_TAXID"
                    severity = "WARN"
                elif candidate_basis == "REFERENCE_TAXID":
                    status = "MISSING_CANDIDATE_REFERENCE_TAXID"
                    severity = "WARN"
                elif candidate_basis == "SPECIES_TAXID":
                    status = "MISSING_CANDIDATE_SPECIES"
                    severity = "WARN"
                else:
                    status = "MISSING_CANDIDATE_REFERENCE_SPECIES"
                    severity = "WARN"
                suggested_kegg = suggested
                suggested_action = "ADD"
                notes.append(suggestion_reason)
            else:
                status = "MULTIPLE_CANDIDATES"
                severity = "WARN"
                suggested_kegg = ""
                suggested_action = "MANUAL_REVIEW"
                notes.append(suggestion_reason)

        reports.append(
            {
                "CSV_line": str(csv_line),
                "Taxon_ID": taxid,
                "Scientific_name": scientific_name,
                "Reference_Taxon_ID": reference_taxid,
                "Reference_Scientific_name": reference_name,
                "Common_name": common_name,
                "kegg_csv": current_code,
                "kegg_exists": "YES" if kegg_entry is not None else ("NO" if current_code else ""),
                "kegg_genome_id": kegg_entry.get("genome_id", "") if kegg_entry else "",
                "kegg_name": kegg_entry.get("name", "") if kegg_entry else "",
                "kegg_taxid": kegg_entry.get("taxid", "") if kegg_entry else "",
                "kegg_species_taxid": kegg_entry.get("species_taxid", "") if kegg_entry else "",
                "status": status,
                "severity": severity,
                "candidate_basis": candidate_basis,
                "candidate_codes": candidate_codes,
                "suggested_kegg": suggested_kegg,
                "suggested_action": suggested_action,
                "notes": "; ".join(notes),
            }
        )

    return reports


def main() -> int:
    args = parse_args()
    input_path = Path(args.hostspecies).expanduser().resolve()
    report_path = Path(args.report).expanduser().resolve()
    cache_dir = Path(args.cache_dir).expanduser().resolve()

    _, rows = read_hostspecies(input_path)
    downloader = Downloader(
        cache_dir=cache_dir,
        timeout=args.timeout,
        retries=args.retries,
        request_delay=args.request_delay,
        cache_max_age_days=args.cache_max_age_days,
    )

    say("============================================================")
    say("MTD Explorer — HostSpecies KEGG auditor")
    say(f"Input:           {input_path}")
    say(f"Rows:            {len(rows)}")
    say(f"Report:          {report_path}")
    say(f"Cache directory: {cache_dir}")
    say("Mode:            READ-ONLY (HostSpecies.csv will not be modified)")
    say("============================================================")

    say("[INFO] Retrieving KEGG genome/organism catalog...")
    organism_text = downloader.get_text(
        f"{KEGG_REST}/list/genome",
        "kegg_list_genome.tsv",
    )
    by_code, _ = parse_organism_list(organism_text)

    say("[INFO] Retrieving KEGG exact genome -> NCBI TaxID mapping...")
    exact_text = downloader.get_text(
        f"{KEGG_REST}/link/taxonomy/genome",
        "kegg_link_taxonomy_genome.tsv",
    )
    exact_taxid_by_genome = parse_genome_taxonomy_links(exact_text)

    say("[INFO] Retrieving KEGG species-level genome -> NCBI TaxID mapping...")
    species_text = downloader.get_text(
        f"{KEGG_REST}/link/taxonomy/genome/species",
        "kegg_link_taxonomy_genome_species.tsv",
    )
    species_taxid_by_genome = parse_genome_taxonomy_links(species_text)

    by_code, exact_index, species_index = enrich_catalog(
        by_code,
        exact_taxid_by_genome,
        species_taxid_by_genome,
    )

    reports = audit_rows(rows, by_code, exact_index, species_index)
    write_report(report_path, reports)

    counts = Counter(row["status"] for row in reports)
    severity_counts = Counter(row["severity"] for row in reports)
    problems = [row for row in reports if row["status"] in PROBLEM_STATUSES]

    say("============================================================")
    say("Audit summary")
    say(f"Rows evaluated:       {len(reports)}")
    say(f"OK:                   {severity_counts.get('OK', 0)}")
    say(f"Warnings:             {severity_counts.get('WARN', 0)}")
    say(f"Errors:               {severity_counts.get('ERROR', 0)}")
    say(f"No KEGG entry (info): {counts.get('NO_KEGG_ENTRY', 0)}")
    say("")
    for status in sorted(counts):
        say(f"  {status:<34} {counts[status]}")

    if problems:
        say("")
        say("High-confidence problems:")
        for row in problems:
            current = row["kegg_csv"] or "-"
            suggestion = row["suggested_kegg"] or "-"
            say(
                f"  line {row['CSV_line']}: {row['Scientific_name']} "
                f"[{row['Taxon_ID']}] kegg={current} -> {row['status']} "
                f"(suggested={suggestion}; action={row['suggested_action']})"
            )

    say("")
    say(f"Report written: {report_path}")
    say("============================================================")

    if args.fail_on_problems and problems:
        return 2
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        say("\n[ERROR] Interrupted by user.")
        raise SystemExit(130)
    except Exception as exc:
        say(f"[ERROR] {exc}")
        raise SystemExit(1)
