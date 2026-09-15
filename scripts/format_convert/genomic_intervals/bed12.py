"""BED12-style block records and interval arithmetic.

Consumers such as the lncRNA conservation statistics stage read
pslMap-derived group/length/exon tables whose columns follow the BED12
ordering. This module parses those records and exposes the interval algebra
(block expansion, union and overlap lengths, overlap fraction, span and
distance tests) shared across pipelines.
"""

from dataclasses import dataclass


BED12_FIELDS = [
    "chromosome",
    "chromstart",
    "chromend",
    "ids",
    "strand",
    "blockcount",
    "blocksizes",
    "blockstarts",
]


def parse_int(value, default=0):
    try:
        if value in ("", None):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


@dataclass
class BlockRecord:
    chrom: str
    start: int
    end: int
    ident: str
    strand: str
    block_count: int
    block_sizes: list
    block_starts: list

    @classmethod
    def from_row(cls, row):
        return cls(
            chrom=row.get("chromosome") or row.get("chrom") or row.get("chr") or "",
            start=parse_int(row.get("chromstart") or row.get("start")),
            end=parse_int(row.get("chromend") or row.get("end")),
            ident=row.get("ids") or row.get("id") or row.get("name") or "",
            strand=row.get("strand") or ".",
            block_count=parse_int(row.get("blockcount") or row.get("blockCount"), 1),
            block_sizes=parse_list(row.get("blocksizes") or row.get("blockSizes")),
            block_starts=parse_list(row.get("blockstarts") or row.get("blockStarts")),
        )

    def intervals(self):
        if not self.block_sizes:
            return [(self.start, self.end)]
        out = []
        starts = self.block_starts or [0] * len(self.block_sizes)
        for rel_start, size in zip(starts, self.block_sizes):
            left = self.start + rel_start
            out.append((left, left + size))
        return out

    def span_intersects(self, other):
        return not (self.end < other.start or other.end < self.start)

    def distance_within(self, other, max_distance):
        if self.end < other.start:
            return other.start - self.end <= max_distance
        if other.end < self.start:
            return self.start - other.end <= max_distance
        return True


def parse_list(value):
    if value is None:
        return []
    text = str(value).strip().rstrip(",")
    if not text:
        return []
    return [parse_int(item) for item in text.split(",") if item != ""]


def interval_union_length(intervals):
    if not intervals:
        return 0
    ordered = sorted(intervals)
    total = 0
    cur_start, cur_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
    return total + cur_end - cur_start


def overlap_length(a_intervals, b_intervals):
    total = 0
    for a_start, a_end in a_intervals:
        for b_start, b_end in b_intervals:
            left = max(a_start, b_start)
            right = min(a_end, b_end)
            if left < right:
                total += right - left
    return total


def overlap_fraction_of_smaller(a_intervals, b_intervals):
    a_len = interval_union_length(a_intervals)
    b_len = interval_union_length(b_intervals)
    denom = min(a_len, b_len)
    if denom <= 0:
        return 0.0
    return overlap_length(a_intervals, b_intervals) / denom
