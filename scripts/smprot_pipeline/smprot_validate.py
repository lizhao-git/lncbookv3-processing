#!/usr/bin/env python3
"""Thin wrapper that forwards positional input files as repeated --input-file
arguments to validate_smprot_tsv.py.

CWL cannot natively emit a repeated ``--input-file <f>`` flag per array item, so
this wrapper accepts the files as positional arguments and forwards them
correctly. It keeps the original validate_smprot_tsv.py unchanged.
"""
import argparse
import subprocess
import sys
from pathlib import Path

VALIDATE_SCRIPT = Path(__file__).resolve().parent / "validate_smprot_tsv.py"


def main():
    parser = argparse.ArgumentParser(description="Validate SmProt coordinates (CWL wrapper)")
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()

    cmd = [
        sys.executable,
        str(VALIDATE_SCRIPT),
        "--output-tsv",
        args.output_tsv,
        "--report",
        args.report,
    ]
    for path in args.files:
        cmd += ["--input-file", path]

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
