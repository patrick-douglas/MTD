#!/usr/bin/env python3
"""Safely integrate the independent HPC backend into MTD_explorer.sh."""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile

MARKER = "MTD_HPC_BACKEND_20260730"


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected exactly one match, found {count}")
    return text.replace(old, new, 1)


def regex_replace_once(
    text: str,
    pattern: re.Pattern[str],
    replacement: str,
    label: str,
) -> str:
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(f"{label}: expected exactly one match, found {len(matches)}")
    match = matches[0]
    return text[: match.start()] + replacement + text[match.end() :]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--check", action="store_true")
    return parser.parse_args()


def patch(text: str) -> str:
    if MARKER in text:
        required_tokens = (
            '--hpc-conf FILE',
            'mtd_hpc_humann_stage.sh',
            'mtd_hpc_magicblast_stage.sh',
            'mtd_hpc_load_config "$HPC_CONF" "$MTDIR"',
        )
        missing = [token for token in required_tokens if token not in text]
        if missing:
            raise RuntimeError(
                "HPC marker exists, but the integration is incomplete. Missing: "
                + ", ".join(missing)
            )
        return text

    # Allow one or more blank lines between the layout defaults and show_help().
    default_pattern = re.compile(
        r'READ_LAYOUT_MODE="auto"\n'
        r'READ_LAYOUT=""\n'
        r'(?P<spacing>\n*)'
        r'show_help\(\) \{'
    )
    default_matches = list(default_pattern.finditer(text))
    if len(default_matches) != 1:
        raise RuntimeError(
            f"HPC default variable: expected exactly one match, found {len(default_matches)}"
        )
    default_match = default_matches[0]
    spacing = default_match.group("spacing") or "\n"
    default_replacement = (
        'READ_LAYOUT_MODE="auto"\n'
        'READ_LAYOUT=""\n'
        f'{spacing}'
        f'# {MARKER}: optional independent Slurm backend.\n'
        'HPC_CONF=""\n'
        'show_help() {'
    )
    text = text[: default_match.start()] + default_replacement + text[default_match.end() :]

    text = replace_once(
        text,
        '      --threads INT                        Number of CPU threads\n'
        '                                           Default: nproc = ${threads}\n',
        '      --threads INT                        Number of CPU threads\n'
        '                                           Default: nproc = ${threads}\n'
        '      --hpc-conf FILE                      Enable the independent Slurm backend using FILE.\n'
        '                                           MetaPhlAn/HUMAnN and Magic-BLAST run on configured nodes.\n'
        '                                           Without this option, execution remains local.\n',
        "HPC help option",
    )

    text = replace_once(
        text,
        '        --ssgsea-gmt)\n',
        '        --hpc-conf)\n'
        '            if [[ $# -lt 2 || -z "${2:-}" ]]; then\n'
        '                die "--hpc-conf requires a configuration file"\n'
        '            fi\n'
        '            HPC_CONF="$2"\n'
        '            shift 2\n'
        '            ;;\n'
        '        --hpc-conf=*)\n'
        '            HPC_CONF="${1#*=}"\n'
        '            [[ -n "$HPC_CONF" ]] || die "--hpc-conf requires a configuration file"\n'
        '            shift\n'
        '            ;;\n'
        '        --ssgsea-gmt)\n',
        "HPC argument parsing",
    )

    text = replace_once(
        text,
        'if [[ -z "${hostid:-}" ]]; then\n'
        '    die "Missing required argument: -h or --hostid TAXID"\n'
        'fi\n',
        'if [[ -z "${hostid:-}" ]]; then\n'
        '    die "Missing required argument: -h or --hostid TAXID"\n'
        'fi\n'
        'if [[ -n "$HPC_CONF" ]]; then\n'
        '    require_file "$HPC_CONF" "HPC configuration"\n'
        '    HPC_CONF="$(readlink -f "$HPC_CONF")"\n'
        'fi\n',
        "HPC early validation",
    )

    text = replace_once(
        text,
        'MTDIR=$(dirname "$(readlink -f "$0")")\n',
        'MTDIR=$(dirname "$(readlink -f "$0")")\n'
        'if [[ -n "$HPC_CONF" ]]; then\n'
        '    MTD_HPC_COMMON="$MTDIR/aux_scripts/hpc/mtd_hpc_common.sh"\n'
        '    require_file "$MTD_HPC_COMMON" "MTD HPC common library"\n'
        '    # shellcheck source=aux_scripts/hpc/mtd_hpc_common.sh\n'
        '    source "$MTD_HPC_COMMON"\n'
        '    mtd_hpc_load_config "$HPC_CONF" "$MTDIR" || die "Invalid HPC configuration: $HPC_CONF"\n'
        'fi\n',
        "HPC initialization",
    )

    text = replace_once(
        text,
        'echo "  threads:                        $threads"\n'
        'echo "  host alignment mode:            $blast"\n',
        'echo "  threads:                        $threads"\n'
        'echo "  HPC backend:                    $([[ -n "$HPC_CONF" ]] && echo enabled || echo disabled)"\n'
        'echo "  HPC configuration:              ${HPC_CONF:-none}"\n'
        'echo "  host alignment mode:            $blast"\n',
        "HPC summary",
    )

    # Restrict the HUMAnN replacement to the documented HUMAnN execution section.
    # The previous broad expression could begin at an earlier sample_index loop.
    humann_pattern = re.compile(
        r'(?P<header>'
        r'# ------------------------------------------------------------\n'
        r'# Run HUMAnN once per sample\n'
        r'# ------------------------------------------------------------\n\n'
        r'echo "\$\{g\}Run HUMAnN3\$\{w\}"\n\n'
        r')'
        r'(?P<loop>'
        r'sample_index=0\n'
        r'for i in \$lsn; do\n'
        r'.*?'
        r'echo "\$\{g\}\[OK\] HUMAnN completed for sample:\$\{w\} \$i"\n'
        r'done\n'
        r')'
        r'(?=\ncd "\$HUMANN_WORK_DIR")',
        re.S,
    )
    humann_matches = list(humann_pattern.finditer(text))
    if len(humann_matches) != 1:
        raise RuntimeError(f"HUMAnN loop: expected one match, found {len(humann_matches)}")
    humann_match = humann_matches[0]
    original_humann = humann_match.group("loop")
    indented_humann = "\n".join(
        "    " + line if line else line
        for line in original_humann.rstrip("\n").split("\n")
    )
    new_humann = (
        humann_match.group("header")
        + 'if [[ -n "$HPC_CONF" ]]; then\n'
        + '    MTD_HPC_HUMANN_WORK_DIR="$HUMANN_WORK_DIR/hpc_humann"\n'
        + '    MTD_HPC_SAMPLE_LIST="$MTD_HPC_HUMANN_WORK_DIR/samples.txt"\n'
        + '    mkdir -p -- "$MTD_HPC_HUMANN_WORK_DIR"\n'
        + '    printf \'%s\\n\' $lsn > "$MTD_HPC_SAMPLE_LIST"\n'
        + '    bash "$MTDIR/aux_scripts/hpc/mtd_hpc_humann_stage.sh" \\\n'
        + '        --hpc-conf "$HPC_CONF" \\\n'
        + '        --samples "$MTD_HPC_SAMPLE_LIST" \\\n'
        + '        --input-dir "$HUMANN_INPUT_DIR" \\\n'
        + '        --output-dir "$HUMANN_RESULTS_DIR" \\\n'
        + '        --work-dir "$MTD_HPC_HUMANN_WORK_DIR" \\\n'
        + '        --mtd-root "$MTDIR" || die "HPC HUMAnN stage failed"\n'
        + 'else\n'
        + f'{indented_humann}\n'
        + 'fi\n'
    )
    text = text[: humann_match.start()] + new_humann + text[humann_match.end() :]

    blast_pattern = re.compile(
        r'if \[\[ "\$blast" == "blast" \]\]; then\n'
        r'    echo "\$\{g\}Magic-BLAST\$\{w\}"\n\n'
        r'    for i in \$lsn; do\n.*?'
        r'    done\n\n'
        r'else\n'
        r'    echo "\$\{g\}HISAT2 alignment\$\{w\}"',
        re.S,
    )
    blast_matches = list(blast_pattern.finditer(text))
    if len(blast_matches) != 1:
        raise RuntimeError(f"Magic-BLAST branch: expected one match, found {len(blast_matches)}")
    whole = blast_matches[0].group(0)
    local_body = whole.split('    for i in $lsn; do\n', 1)[1].rsplit('    done\n\nelse\n', 1)[0]
    local_loop = '    for i in $lsn; do\n' + local_body + '    done\n'
    local_indented = "\n".join(
        "    " + line if line else line
        for line in local_loop.rstrip("\n").split("\n")
    )

    new_blast = (
        'if [[ "$blast" == "blast" ]]; then\n'
        '    echo "${g}Magic-BLAST${w}"\n\n'
        '    if [[ -n "$HPC_CONF" ]]; then\n'
        '        MTD_HPC_MAGICBLAST_WORK_DIR="$outputdr/hpc/magicblast"\n'
        '        MTD_HPC_MAGICBLAST_INPUTS="$MTD_HPC_MAGICBLAST_WORK_DIR/inputs.tsv"\n'
        '        mkdir -p -- "$MTD_HPC_MAGICBLAST_WORK_DIR"\n'
        '        : > "$MTD_HPC_MAGICBLAST_INPUTS"\n\n'
        '        for i in $lsn; do\n'
        '            set_host_kraken_fastq_paths "$i"\n'
        '            sam_file="$PIPELINE_TEMP_DIR/${i}.sam"\n'
        '            host_sam_files+=("${i}.sam")\n'
        '            if [[ "$READ_LAYOUT" == "pe" ]]; then\n'
        '                require_file "$HOST_KRAKEN_R1" "Kraken2 host-classified R1 for Magic-BLAST sample $i"\n'
        '                require_file "$HOST_KRAKEN_R2" "Kraken2 host-classified R2 for Magic-BLAST sample $i"\n'
        '                printf \'%s\\t%s\\t%s\\t%s\\t%s\\n\' \\\n'
        '                    "$i" "$READ_LAYOUT" "$HOST_KRAKEN_R1" "$HOST_KRAKEN_R2" "$sam_file" \\\n'
        '                    >> "$MTD_HPC_MAGICBLAST_INPUTS"\n'
        '            else\n'
        '                require_file "$HOST_KRAKEN_R1" "Kraken2 host-classified SE FASTQ for Magic-BLAST sample $i"\n'
        '                printf \'%s\\t%s\\t%s\\t-\\t%s\\n\' \\\n'
        '                    "$i" "$READ_LAYOUT" "$HOST_KRAKEN_R1" "$sam_file" \\\n'
        '                    >> "$MTD_HPC_MAGICBLAST_INPUTS"\n'
        '            fi\n'
        '        done\n\n'
        '        bash "$MTDIR/aux_scripts/hpc/mtd_hpc_magicblast_stage.sh" \\\n'
        '            --hpc-conf "$HPC_CONF" \\\n'
        '            --input-manifest "$MTD_HPC_MAGICBLAST_INPUTS" \\\n'
        '            --main-database-prefix "$DB_blast" \\\n'
        '            --work-dir "$MTD_HPC_MAGICBLAST_WORK_DIR" \\\n'
        '            --mtd-root "$MTDIR" || die "HPC Magic-BLAST stage failed"\n\n'
        '        for i in $lsn; do\n'
        '            require_file "$PIPELINE_TEMP_DIR/${i}.sam" "HPC Magic-BLAST SAM output for sample $i"\n'
        '        done\n'
        '    else\n'
        f'{local_indented}\n'
        '    fi\n\n'
        'else\n'
        '    echo "${g}HISAT2 alignment${w}"'
    )
    text = text[: blast_matches[0].start()] + new_blast + text[blast_matches[0].end() :]

    return text


def validate_bash(text: str, target: Path) -> None:
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=target.parent,
        prefix=f".{target.name}.hpc-check.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        handle.write(text)

    try:
        result = subprocess.run(
            ["bash", "-n", str(temporary)],
            text=True,
            capture_output=True,
        )
        if result.returncode != 0:
            detail = result.stderr.strip() or "unknown Bash syntax error"
            raise RuntimeError(
                "Generated MTD_explorer.sh failed bash -n; the original was not changed. "
                + detail
            )
    finally:
        temporary.unlink(missing_ok=True)


def main() -> int:
    args = parse_args()
    repo = args.repo.resolve()
    target = repo / "MTD_explorer.sh"
    if not target.is_file():
        raise SystemExit(f"MTD_explorer.sh not found: {target}")

    original = target.read_text(encoding="utf-8")

    # Refuse to use a malformed source script as the patch base.
    validate_bash(original, target)
    updated = patch(original)
    validate_bash(updated, target)

    if args.check:
        if updated == original and MARKER in original:
            print("[OK] HPC backend marker is already present.")
        else:
            print("[OK] Current MTD_explorer.sh matches all required patch markers.")
        return 0

    if updated == original:
        print("[OK] HPC backend is already integrated; no changes made.")
        return 0

    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = target.with_name(f"MTD_explorer.sh.backup_before_hpc_{stamp}")
    original_mode = stat.S_IMODE(target.stat().st_mode)
    shutil.copy2(target, backup)

    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=target.parent,
            prefix=f".{target.name}.hpc-install.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            handle.write(updated)

        os.chmod(temporary_path, original_mode | 0o111)
        os.replace(temporary_path, target)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    print(f"[OK] Updated: {target}")
    print(f"[OK] Backup:  {backup}")
    print("[OK] bash -n accepted the generated MTD_explorer.sh")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
