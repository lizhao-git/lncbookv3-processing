#!/usr/bin/env python3
"""Stage D: rescale methylation values in a gene x sample matrix.

Typical use: convert percent values (0-100) to beta values (0-1) with
``--scale 0.01``. Comma-separated multi-value cells are handled per token.

Fixes the legacy indentation bug in ``Analyze_and_format/beta_value.py``,
where ``new_content.append(...)`` sat outside the inner ``for`` loop.
"""

import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--scale", type=float, default=0.01,
                    help="multiplier (0.01 for percent -> beta)")
    args = ap.parse_args()

    with open(args.input) as fin, open(args.output, "w") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            out = [parts[0]]
            for cell in parts[1:]:
                new_tokens = []
                for token in cell.split(","):
                    try:
                        new_tokens.append(str(float(token) * args.scale))
                    except ValueError:
                        new_tokens.append(token)
                out.append(",".join(new_tokens))
            fout.write("\t".join(out) + "\n")

    print(f"done: {args.output}")


if __name__ == "__main__":
    main()
