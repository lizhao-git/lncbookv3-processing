#!/usr/bin/env python3
"""Convert a bigBed track back to BED or a bigWig track back to bedGraph.

The input format is detected from the file extension: ``.bigBed``/``.bb``
use ``bigBedToBed`` and ``.bigWig``/``.bw`` use ``bigWigToBedGraph``.
Requires the UCSC kent tools, which are installed in the
lncbookv3-processing Docker image.
"""
import argparse
import os

from genomic_intervals.io import open_text
from genomic_intervals.tooling import require_tool, run_command

BIGBED_EXTENSIONS = (".bigbed", ".bb")
BIGWIG_EXTENSIONS = (".bigwig", ".bw")


def detect_tool(input_path: str) -> str:
    lower = input_path.lower()
    if lower.endswith(BIGBED_EXTENSIONS):
        return "bigBedToBed"
    if lower.endswith(BIGWIG_EXTENSIONS):
        return "bigWigToBedGraph"
    raise SystemExit(
        f"Cannot infer track type from {input_path!r}; "
        "expected a .bigBed/.bb or .bigWig/.bw file"
    )


def bigtrack_to_text(input_path: str, output_bed: str, report_path: str = None):
    tool_name = detect_tool(input_path)
    tool = require_tool(tool_name)
    run_command([tool, input_path, output_bed])

    rows = 0
    with open_text(output_bed) as (fh, _compression):
        for raw in fh:
            if raw.strip():
                rows += 1

    if report_path:
        with open(report_path, "w", encoding="utf-8") as rep:
            rep.write("metric\tvalue\n")
            rep.write(f"tool\t{tool_name}\n")
            rep.write(f"rows\t{rows}\n")
            if os.path.exists(output_bed):
                rep.write(f"output_bytes\t{os.path.getsize(output_bed)}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert bigBed to BED or bigWig to bedGraph")
    parser.add_argument("--input-track", required=True)
    parser.add_argument("--output-bed", required=True,
                        help="output text file: BED for bigBed, bedGraph for bigWig")
    parser.add_argument("--report")
    args = parser.parse_args()
    bigtrack_to_text(args.input_track, args.output_bed, args.report)


if __name__ == "__main__":
    main()
