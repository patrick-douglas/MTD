#!/usr/bin/env python3

"""Build a persistent host gene-ID to NCBI GeneID mapping.

The primary source is NCBI gene2ensembl for references whose GTF gene IDs are
true Ensembl identifiers.  When the reference uses an external namespace
(e.g. VectorBase BGLAX identifiers), an optional NCBI gene_info file is used as
an exact-match fallback across dbXrefs, locus tags, symbols and synonyms.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
from collections import defaultdict
from pathlib import Path


def normalize_gene_id(value: str) -> str:
    value = str(value).strip()
    return re.sub(r"\.\d+$", "", value)


def read_gtf_gene_ids(path: Path) -> list[str]:
    genes: list[str] = []
    seen: set[str] = set()

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            gene = normalize_gene_id(line)
            if gene and gene not in seen:
                genes.append(gene)
                seen.add(gene)

    return genes


def add_mapping(
    mapping: dict[str, set[str]],
    symbols: dict[str, str],
    gene_id: str,
    entrez: str,
    symbol: str = "",
) -> None:
    gene_id = normalize_gene_id(gene_id)
    entrez = str(entrez).strip()
    symbol = str(symbol).strip()

    if not gene_id or gene_id == "-":
        return
    if not entrez or entrez == "-" or not entrez.isdigit():
        return

    mapping[gene_id].add(entrez)

    if symbol and symbol != "-" and not symbols.get(gene_id):
        symbols[gene_id] = symbol


def read_gene2ensembl(
    path: Path,
    reference_taxid: str,
    wanted: set[str],
) -> tuple[dict[str, set[str]], dict[str, str], int]:
    mapping: dict[str, set[str]] = defaultdict(set)
    symbols: dict[str, str] = {}
    reference_pairs: set[tuple[str, str]] = set()

    with gzip.open(
        path,
        mode="rt",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.reader(handle, delimiter="\t")

        try:
            raw_header = next(reader)
        except StopIteration:
            raise SystemExit("[ERROR] gene2ensembl.gz is empty.")

        header = [column.lstrip("#").strip() for column in raw_header]
        required = {"tax_id", "GeneID", "Ensembl_gene_identifier"}
        missing = required.difference(header)

        if missing:
            raise SystemExit(
                "[ERROR] Missing gene2ensembl columns: "
                + ", ".join(sorted(missing))
            )

        tax_index = header.index("tax_id")
        entrez_index = header.index("GeneID")
        ensembl_index = header.index("Ensembl_gene_identifier")
        maximum_index = max(tax_index, entrez_index, ensembl_index)

        for row in reader:
            if len(row) <= maximum_index:
                continue
            if row[tax_index].strip() != reference_taxid:
                continue

            ensembl = normalize_gene_id(row[ensembl_index])
            entrez = row[entrez_index].strip()

            if not ensembl or ensembl == "-":
                continue
            if not entrez or entrez == "-" or not entrez.isdigit():
                continue

            reference_pairs.add((ensembl, entrez))

            if ensembl in wanted:
                add_mapping(mapping, symbols, ensembl, entrez)

    return mapping, symbols, len(reference_pairs)


def split_pipe_field(value: str) -> list[str]:
    value = str(value).strip()
    if not value or value == "-":
        return []
    return [part.strip() for part in value.split("|") if part.strip()]


def gene_info_candidate_ids(row: dict[str, str]) -> set[str]:
    candidates: set[str] = set()

    for column in ("Symbol", "LocusTag"):
        value = normalize_gene_id(row.get(column, ""))
        if value and value != "-":
            candidates.add(value)

    for synonym in split_pipe_field(row.get("Synonyms", "")):
        synonym = normalize_gene_id(synonym)
        if synonym and synonym != "-":
            candidates.add(synonym)

    for dbxref in split_pipe_field(row.get("dbXrefs", "")):
        # dbXrefs are usually DATABASE:IDENTIFIER.  The GTF commonly uses the
        # identifier component (e.g. VectorBase:BGLAX_050876).
        identifier = dbxref.split(":", 1)[1] if ":" in dbxref else dbxref
        identifier = normalize_gene_id(identifier)
        if identifier and identifier != "-":
            candidates.add(identifier)

    return candidates


def read_gene_info(
    path: Path,
    reference_taxid: str,
    wanted: set[str],
    already_mapped: set[str],
) -> tuple[dict[str, set[str]], dict[str, str], int]:
    mapping: dict[str, set[str]] = defaultdict(set)
    symbols: dict[str, str] = {}
    reference_rows = 0

    with gzip.open(
        path,
        mode="rt",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise SystemExit("[ERROR] gene_info file has no header.")

        reader.fieldnames = [
            field.lstrip("#").strip() if field is not None else ""
            for field in reader.fieldnames
        ]

        required = {"tax_id", "GeneID", "Symbol", "dbXrefs"}
        missing = required.difference(reader.fieldnames)

        if missing:
            raise SystemExit(
                "[ERROR] Missing gene_info columns: "
                + ", ".join(sorted(missing))
            )

        remaining = wanted.difference(already_mapped)

        for row in reader:
            if str(row.get("tax_id", "")).strip() != reference_taxid:
                continue

            reference_rows += 1
            entrez = str(row.get("GeneID", "")).strip()
            symbol = str(row.get("Symbol", "")).strip()

            if not entrez.isdigit():
                continue

            for candidate in gene_info_candidate_ids(row):
                if candidate in remaining:
                    add_mapping(mapping, symbols, candidate, entrez, symbol)

    return mapping, symbols, reference_rows


def unambiguous_rows(
    ordered_gtf_genes: list[str],
    mapping: dict[str, set[str]],
    symbols: dict[str, str],
) -> tuple[list[tuple[str, str, str]], set[str]]:
    rows: list[tuple[str, str, str]] = []
    ambiguous: set[str] = set()

    for gene_id in ordered_gtf_genes:
        geneids = mapping.get(gene_id, set())

        if len(geneids) == 1:
            entrez = next(iter(geneids))
            rows.append((gene_id, symbols.get(gene_id, ""), entrez))
        elif len(geneids) > 1:
            ambiguous.add(gene_id)

    return rows, ambiguous


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene2ensembl", required=True)
    parser.add_argument("--gene-info")
    parser.add_argument(
        "--probe-direct",
        action="store_true",
        help=(
            "Probe gene2ensembl only. Exit with status 3 when the file "
            "contains mappings for the reference TaxID but none match the "
            "selected GTF namespace; this signals the caller to try gene_info."
        ),
    )
    parser.add_argument("--reference-taxid", required=True)
    parser.add_argument("--requested-taxid", required=True)
    parser.add_argument("--gtf-genes", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    gene2ensembl = Path(args.gene2ensembl)
    gene_info = Path(args.gene_info) if args.gene_info else None
    gtf_file = Path(args.gtf_genes)
    output = Path(args.output)
    summary = Path(args.summary)

    if not gene2ensembl.is_file() or gene2ensembl.stat().st_size == 0:
        raise SystemExit(
            f"[ERROR] gene2ensembl cache is missing or empty: {gene2ensembl}"
        )

    if args.probe_direct and gene_info is not None:
        raise SystemExit(
            "[ERROR] --probe-direct cannot be combined with --gene-info."
        )

    if gene_info is not None and (
        not gene_info.is_file() or gene_info.stat().st_size == 0
    ):
        raise SystemExit(f"[ERROR] gene_info cache is missing or empty: {gene_info}")

    if not gtf_file.is_file() or gtf_file.stat().st_size == 0:
        raise SystemExit(f"[ERROR] GTF gene-ID list is missing or empty: {gtf_file}")

    ordered_gtf_genes = read_gtf_gene_ids(gtf_file)
    if not ordered_gtf_genes:
        raise SystemExit("[ERROR] No gene IDs were read from the GTF list.")

    wanted = set(ordered_gtf_genes)

    direct_map, direct_symbols, ncbi_reference_pairs = read_gene2ensembl(
        gene2ensembl,
        args.reference_taxid,
        wanted,
    )

    direct_rows, direct_ambiguous = unambiguous_rows(
        ordered_gtf_genes,
        direct_map,
        direct_symbols,
    )
    direct_mapped = {row[0] for row in direct_rows}

    combined_map: dict[str, set[str]] = defaultdict(set)
    combined_symbols: dict[str, str] = {}

    for gene_id, ids in direct_map.items():
        combined_map[gene_id].update(ids)
    combined_symbols.update(direct_symbols)

    gene_info_reference_rows = 0
    fallback_mapped: set[str] = set()

    if gene_info is not None:
        fallback_map, fallback_symbols, gene_info_reference_rows = read_gene_info(
            gene_info,
            args.reference_taxid,
            wanted,
            already_mapped=direct_mapped.union(direct_ambiguous),
        )

        for gene_id, ids in fallback_map.items():
            combined_map[gene_id].update(ids)
        for gene_id, symbol in fallback_symbols.items():
            if symbol and not combined_symbols.get(gene_id):
                combined_symbols[gene_id] = symbol

        fallback_rows, _ = unambiguous_rows(
            ordered_gtf_genes,
            fallback_map,
            fallback_symbols,
        )
        fallback_mapped = {row[0] for row in fallback_rows}

    rows, ambiguous = unambiguous_rows(
        ordered_gtf_genes,
        combined_map,
        combined_symbols,
    )
    mapped_gtf_genes = {row[0] for row in rows}

    if not mapped_gtf_genes:
        if gene_info is None:
            if args.probe_direct:
                print(
                    "[INFO] No direct gene2ensembl mapping matched the selected "
                    "GTF namespace; gene_info fallback is required."
                )
                raise SystemExit(3)

            raise SystemExit(
                "[WARNING] NCBI mappings were found, but none matched the gene "
                "IDs from the selected GTF. A gene_info fallback may be needed."
            )
        raise SystemExit(
            "[WARNING] Neither gene2ensembl nor gene_info produced a usable "
            "mapping for the selected GTF gene IDs."
        )

    coverage = 100.0 * len(mapped_gtf_genes) / len(wanted)
    unmapped = wanted.difference(mapped_gtf_genes).difference(ambiguous)

    output.parent.mkdir(parents=True, exist_ok=True)
    summary.parent.mkdir(parents=True, exist_ok=True)
    output_tmp = output.with_suffix(output.suffix + ".tmp")
    summary_tmp = summary.with_suffix(summary.suffix + ".tmp")

    with output_tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["ensembl_gene_id", "external_gene_name", "entrezgene_id"])
        writer.writerows(rows)

    with summary_tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["parameter", "value"])
        writer.writerow(["requested_taxid", args.requested_taxid])
        writer.writerow(["reference_taxid", args.reference_taxid])
        writer.writerow(["gtf_unique_genes", len(wanted)])
        writer.writerow(["ncbi_reference_pairs", ncbi_reference_pairs])
        writer.writerow(["gene2ensembl_mapped_genes", len(direct_mapped)])
        writer.writerow(["gene_info_reference_rows", gene_info_reference_rows])
        writer.writerow(["gene_info_fallback_mapped_genes", len(fallback_mapped)])
        writer.writerow(["ambiguous_gtf_genes_excluded", len(ambiguous)])
        writer.writerow(["unmapped_gtf_genes", len(unmapped)])
        writer.writerow(["mapped_gtf_genes", len(mapped_gtf_genes)])
        writer.writerow(["gtf_mapping_coverage_pct", f"{coverage:.4f}"])
        writer.writerow(["gene2ensembl_source", str(gene2ensembl)])
        writer.writerow(["gene_info_source", str(gene_info) if gene_info else ""])

    os.replace(output_tmp, output)
    os.replace(summary_tmp, summary)

    print("[OK] Reference-matched host GID -> Entrez mapping created")
    print(f"  Requested TaxID:          {args.requested_taxid}")
    print(f"  Reference TaxID:          {args.reference_taxid}")
    print(f"  GTF genes:                {len(wanted)}")
    print(f"  gene2ensembl mapped:      {len(direct_mapped)}")
    print(f"  gene_info fallback mapped:{len(fallback_mapped):>7}")
    print(f"  Ambiguous excluded:       {len(ambiguous)}")
    print(f"  Total mapped:             {len(mapped_gtf_genes)}")
    print(f"  Coverage:                 {coverage:.2f}%")
    print(f"  Mapping:                  {output}")
    print(f"  Summary:                  {summary}")


if __name__ == "__main__":
    main()
