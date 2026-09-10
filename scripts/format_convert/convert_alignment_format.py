#!/usr/bin/env python3
"""Convert alignment files between SAM and BAM with samtools.

The target format is taken from ``--to`` (or inferred from the output
extension when set to ``auto``). ``--sort`` runs ``samtools sort`` instead
of a plain conversion, and ``--index`` additionally builds a BAI index
for BAM output (run it on sorted data).
"""
import argparse

from genomic_intervals.tooling import require_tool, run_command


def infer_target(output_path: str) -> str:
    lower = output_path.lower()
    if lower.endswith(".bam"):
        return "bam"
    if lower.endswith(".sam"):
        return "sam"
    raise SystemExit("Cannot infer target format from the output name; pass --to sam|bam.")


def convert_alignment(input_path: str, output_path: str, target: str,
                      sort: bool = False, index: bool = False, threads: int = 0):
    samtools = require_tool("samtools")
    base = [samtools]
    if threads:
        base += ["-@", str(threads)]

    if sort:
        command = base + ["sort", "-O", target, "-o", output_path, input_path]
    else:
        command = base + ["view", "--no-PG", "-o", output_path]
        command.append("-b" if target == "bam" else "-h")
        command.append(input_path)
    run_command(command)

    if index:
        if target != "bam":
            raise SystemExit("--index is only supported for BAM output.")
        run_command(base + ["index", output_path])


def main():
    parser = argparse.ArgumentParser(description="Convert between SAM and BAM (samtools)")
    parser.add_argument("--input-alignment", required=True)
    parser.add_argument("--output-alignment", required=True)
    parser.add_argument("--to", choices=["auto", "sam", "bam"], default="auto")
    parser.add_argument("--sort", action="store_true", help="sort records before writing")
    parser.add_argument("--index", action="store_true", help="build a BAI index (BAM output only)")
    parser.add_argument("--threads", type=int, default=0, help="samtools extra threads (-@)")
    args = parser.parse_args()

    target = infer_target(args.output_alignment) if args.to == "auto" else args.to
    convert_alignment(args.input_alignment, args.output_alignment, target,
                      sort=args.sort, index=args.index, threads=args.threads)


if __name__ == "__main__":
    main()
