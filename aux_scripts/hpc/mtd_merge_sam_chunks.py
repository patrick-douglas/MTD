#!/usr/bin/env python3
"""Merge Magic-BLAST SAM chunks while retaining one compatible header.

Magic-BLAST can write a different @PG command line for every chunk. Those
program records are allowed to differ. Reference-defining @HD, @SQ and @RG
records must remain compatible.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chunk-list", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_paths(path: Path) -> list[Path]:
    result: list[Path] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        raw = raw.strip()
        if raw and not raw.startswith("#"):
            result.append(Path(raw))
    return result


def compatibility_header(header: list[str]) -> list[str]:
    return [line for line in header if line.startswith(("@HD", "@SQ", "@RG"))]


def main() -> int:
    args = parse_args()
    chunks = read_paths(args.chunk_list)
    if not chunks:
        raise SystemExit("Chunk list is empty")
    for chunk in chunks:
        if not chunk.is_file() or chunk.stat().st_size == 0:
            raise SystemExit(f"Missing or empty SAM chunk: {chunk}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temp = args.output.with_suffix(args.output.suffix + ".tmp")
    reference_header: list[str] | None = None
    reference_compatibility: list[str] | None = None
    alignments = 0

    try:
        with temp.open("wt", encoding="utf-8") as out:
            for index, chunk in enumerate(chunks):
                header: list[str] = []
                body_started = False
                with chunk.open("rt", encoding="utf-8", errors="strict") as handle:
                    for line in handle:
                        if not body_started and line.startswith("@"):
                            header.append(line)
                            continue
                        body_started = True
                        if index == 0 and reference_header is None:
                            reference_header = header.copy()
                            reference_compatibility = compatibility_header(header)
                            out.writelines(reference_header)
                        out.write(line)
                        alignments += 1

                if index == 0 and reference_header is None:
                    reference_header = header.copy()
                    reference_compatibility = compatibility_header(header)
                    out.writelines(reference_header)

                if index > 0 and compatibility_header(header) != reference_compatibility:
                    raise SystemExit(
                        "SAM reference/read-group headers differ between chunks; "
                        f"first mismatch: {chunk}"
                    )

        if reference_header is None or not reference_header:
            raise SystemExit("No SAM header was found in the first chunk")

        temp.replace(args.output)
    except BaseException:
        temp.unlink(missing_ok=True)
        raise

    print(f"chunks={len(chunks)}")
    print(f"alignments={alignments}")
    print(f"output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
