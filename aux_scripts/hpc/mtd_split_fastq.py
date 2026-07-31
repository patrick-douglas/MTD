#!/usr/bin/env python3
"""Split single-end or paired FASTQ files into synchronized record chunks."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import IO, Iterator, Tuple

Record = Tuple[str, str, str, str]


def open_text(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    return path.open("rt", encoding="utf-8", errors="strict")


def records(handle: IO[str], label: str) -> Iterator[Record]:
    index = 0
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        index += 1
        if not qual:
            raise ValueError(f"Incomplete FASTQ record {index} in {label}")
        if not header.startswith("@") or not plus.startswith("+"):
            raise ValueError(f"Malformed FASTQ record {index} in {label}")
        yield header, seq, plus, qual


def normalized_id(header: str) -> str:
    token = header[1:].strip().split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def write_record(handle: IO[str], record: Record) -> None:
    handle.writelines(record)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--read1", required=True, type=Path)
    parser.add_argument("--read2", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--records-per-chunk", required=True, type=int)
    parser.add_argument("--manifest", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.records_per_chunk < 1:
        raise SystemExit("--records-per-chunk must be >= 1")
    if not args.read1.is_file():
        raise SystemExit(f"R1 does not exist: {args.read1}")
    if args.read2 is not None and not args.read2.is_file():
        raise SystemExit(f"R2 does not exist: {args.read2}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)

    chunk_index = 0
    records_in_chunk = 0
    total_records = 0
    out1: IO[str] | None = None
    out2: IO[str] | None = None
    chunk_r1: Path | None = None
    chunk_r2: Path | None = None

    def close_outputs() -> None:
        nonlocal out1, out2
        if out1 is not None:
            out1.close()
            out1 = None
        if out2 is not None:
            out2.close()
            out2 = None

    with open_text(args.read1) as h1, args.manifest.open("wt", encoding="utf-8") as manifest:
        r1_iter = records(h1, str(args.read1))
        if args.read2 is None:
            r2_context = None
        else:
            r2_context = open_text(args.read2)

        try:
            r2_iter = records(r2_context, str(args.read2)) if r2_context is not None else None
            while True:
                try:
                    rec1 = next(r1_iter)
                except StopIteration:
                    rec1 = None

                if r2_iter is not None:
                    try:
                        rec2 = next(r2_iter)
                    except StopIteration:
                        rec2 = None
                else:
                    rec2 = None

                if rec1 is None:
                    if rec2 is not None:
                        raise ValueError("R2 contains more records than R1")
                    break
                if r2_iter is not None and rec2 is None:
                    raise ValueError("R1 contains more records than R2")
                if rec2 is not None and normalized_id(rec1[0]) != normalized_id(rec2[0]):
                    raise ValueError(
                        f"Paired FASTQ IDs differ at record {total_records + 1}: "
                        f"{normalized_id(rec1[0])} != {normalized_id(rec2[0])}"
                    )

                if records_in_chunk == 0:
                    chunk_index += 1
                    chunk_r1 = args.output_dir / f"{args.sample}.chunk{chunk_index:06d}.R1.fq"
                    chunk_r2 = (
                        args.output_dir / f"{args.sample}.chunk{chunk_index:06d}.R2.fq"
                        if r2_iter is not None
                        else None
                    )
                    out1 = chunk_r1.open("wt", encoding="utf-8")
                    out2 = chunk_r2.open("wt", encoding="utf-8") if chunk_r2 else None

                assert out1 is not None and chunk_r1 is not None
                write_record(out1, rec1)
                if rec2 is not None:
                    assert out2 is not None
                    write_record(out2, rec2)

                records_in_chunk += 1
                total_records += 1

                if records_in_chunk == args.records_per_chunk:
                    close_outputs()
                    manifest.write(
                        f"{chunk_index:06d}\t{chunk_r1}\t{chunk_r2 or '-'}\t{records_in_chunk}\n"
                    )
                    records_in_chunk = 0

            if records_in_chunk:
                close_outputs()
                assert chunk_r1 is not None
                manifest.write(
                    f"{chunk_index:06d}\t{chunk_r1}\t{chunk_r2 or '-'}\t{records_in_chunk}\n"
                )
        finally:
            close_outputs()
            if r2_context is not None:
                r2_context.close()

    if total_records == 0:
        raise SystemExit(f"FASTQ contains zero records: {args.read1}")

    print(f"sample={args.sample}")
    print(f"records={total_records}")
    print(f"chunks={chunk_index}")
    print(f"manifest={args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
