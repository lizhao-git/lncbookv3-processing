#!/usr/bin/env python3
"""Add a UCSC track line to a BED/bedGraph file.

Existing ``track`` lines are replaced and ``browser`` lines are dropped so
the output carries exactly one track definition, emitted first as UCSC
expects. All other lines pass through unchanged; blank lines are skipped.
"""
import argparse
import re

from genomic_intervals.io import open_text, open_text_write

COLOR_RE = re.compile(r"^\d{1,3},\d{1,3},\d{1,3}$")
VISIBILITIES = ["hide", "dense", "squish", "pack", "full"]


def build_track_line(args) -> str:
    parts = ["track", f'name="{args.track_name}"']
    if args.description:
        parts.append(f'description="{args.description}"')
    if args.track_type:
        parts.append(f"type={args.track_type}")
    if args.color:
        parts.append(f"color={args.color}")
    if args.visibility:
        parts.append(f"visibility={args.visibility}")
    if args.priority is not None:
        parts.append(f"priority={args.priority}")
    return " ".join(parts)


def validate_color(color: str):
    if not COLOR_RE.match(color):
        raise SystemExit(f"Invalid --color {color!r}; expected R,G,B with values 0-255")
    for part in color.split(","):
        if int(part) > 255:
            raise SystemExit(f"Invalid --color {color!r}; each component must be 0-255")


def make_ucsc_track(input_path: str, output_path: str, track_line: str,
                    report_path: str = None):
    input_lines = 0
    track_lines_replaced = 0
    output_lines = 0
    with open_text_write(output_path) as out:
        out.write(track_line + "\n")
        output_lines += 1
        with open_text(input_path) as (fh, _compression):
            for raw in fh:
                line = raw.rstrip("\n")
                if not line:
                    continue
                input_lines += 1
                if line.startswith("track"):
                    track_lines_replaced += 1
                    continue
                if line.startswith("browser"):
                    continue
                out.write(line + "\n")
                output_lines += 1

    if report_path:
        with open(report_path, "w", encoding="utf-8") as rep:
            rep.write("metric\tvalue\n")
            rep.write(f"input_lines\t{input_lines}\n")
            rep.write(f"track_lines_replaced\t{track_lines_replaced}\n")
            rep.write(f"output_lines\t{output_lines}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Add or replace a UCSC track line on a BED/bedGraph file")
    parser.add_argument("--input-bed", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--track-name", required=True)
    parser.add_argument("--description", default="")
    parser.add_argument("--track-type", default="",
                        help="UCSC track type, e.g. bed 6 +, bedGraph")
    parser.add_argument("--color", default="", help="RGB color, e.g. 0,0,178")
    parser.add_argument("--visibility", choices=VISIBILITIES, default="")
    parser.add_argument("--priority", type=int, default=None)
    parser.add_argument("--report")
    args = parser.parse_args()

    if args.color:
        validate_color(args.color)
    track_line = build_track_line(args)
    make_ucsc_track(args.input_bed, args.output_bed, track_line, args.report)


if __name__ == "__main__":
    main()
