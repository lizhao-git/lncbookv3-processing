#!/usr/bin/env python3
"""Concatenate pipeline step reports into a single human-readable run log.

Each input file is written under a ``===== <name> =====`` section header, in
the order given on the command line. The result is a plain-text log that
captures every validation report and summary produced by a pipeline run.
"""
import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="Aggregate pipeline step reports into one run log")
    parser.add_argument("--output-log", required=True)
    parser.add_argument("--title", default="lncbookv3-processing pipeline run log")
    parser.add_argument("--input", action="append", required=True,
                        help="Report file to include; repeat for each step (order is preserved)")
    args = parser.parse_args()

    with open(args.output_log, "w", encoding="utf-8") as out:
        out.write(f"# {args.title}\n")
        out.write(f"# steps: {len(args.input)}\n")
        for path in args.input:
            out.write(f"\n{'=' * 20} {os.path.basename(path)} {'=' * 20}\n")
            try:
                with open(path, "r", encoding="utf-8", errors="replace") as fh:
                    content = fh.read()
            except OSError as exc:
                content = f"<unreadable: {exc}>\n"
            out.write(content if content.endswith("\n") else content + "\n")


if __name__ == "__main__":
    main()
