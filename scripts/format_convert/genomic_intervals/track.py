"""Helpers for UCSC track outputs: chromosome sizes and BED/bedGraph sorting.

bigBed and bigWig files require input rows sorted in the chromosome order
declared by a ``chrom.sizes`` file. These helpers parse that file and
produce the sorted, comment-free intermediate both kent tools expect.
"""
from .io import open_text


def read_chrom_sizes(path: str):
    """Parse a UCSC ``chrom.sizes`` file into (ordered names, name -> size)."""
    order = []
    sizes = {}
    with open_text(path) as (fh, _compression):
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise SystemExit(f"{path}: line {line_no}: expected 'chrom<TAB>size'")
            name, size_text = parts[0], parts[1]
            try:
                size = int(size_text)
            except ValueError:
                raise SystemExit(f"{path}: line {line_no}: size is not an integer: {size_text!r}")
            if name not in sizes:
                order.append(name)
            sizes[name] = size
    return order, sizes


def sort_track_file(input_path: str, output_path: str, chrom_order, min_columns: int = 3):
    """Sort tab-separated BED/bedGraph rows by chrom.sizes order, start, end.

    Comment (``#``), ``track`` and ``browser`` lines are dropped. Rows on
    chromosomes absent from the chrom.sizes file abort with a clear message,
    because both kent converters require the reference set to match.
    """
    rank = {name: index for index, name in enumerate(chrom_order)}
    rows = []
    with open_text(input_path) as (fh, _compression):
        for line_no, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track ") or line.startswith("browser "):
                continue
            cols = line.split("\t")
            if len(cols) < min_columns:
                raise SystemExit(
                    f"{input_path}: line {line_no}: expected at least {min_columns} columns, got {len(cols)}"
                )
            chrom, start_text, end_text = cols[0], cols[1], cols[2]
            if chrom not in rank:
                raise SystemExit(
                    f"{input_path}: line {line_no}: chromosome {chrom!r} is not in the chromosome sizes file"
                )
            try:
                start, end = int(start_text), int(end_text)
            except ValueError:
                raise SystemExit(f"{input_path}: line {line_no}: start/end are not integers")
            rows.append((rank[chrom], start, end, line_no, line))
    rows.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
    with open(output_path, "w", encoding="utf-8") as out:
        for _rank, _start, _end, _line_no, line in rows:
            out.write(line + "\n")
