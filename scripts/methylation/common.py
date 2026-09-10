#!/usr/bin/env python3
"""Shared helpers for the refactored methylation pipeline.

Indexing conventions (fixed from the legacy scripts):

* BED files are **0-based, half-open** ``[start, end)`` with columns::

      chrom  start  end  value

* Gene info files are **1-based, fully-closed** ``[start, end]`` with columns::

      gene_id  chrom  start  end

  and are converted to 0-based half-open internally before any overlap test.

* The methylation value is always the **4th column** (index 3) of a BED file;
  every format produced by ``preprocess.py`` already guarantees this.
"""

import math


def load_gene_info(path):
    """Return ``{gene_id: (chrom, start0, end0)}`` with ``[start0, end0)``
    0-based half-open."""
    genes = {}
    with open(path) as fh:
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            gid, chrom = parts[0], parts[1]
            try:
                start, end = int(parts[2]), int(parts[3])
            except ValueError:
                continue
            # 1-based closed [start, end] -> 0-based half-open [start-1, end)
            genes[gid] = (chrom, start - 1, end)
    return genes


def mean_over_region(intervals, start0, end0):
    """Mean of the values of every interval overlapping ``[start0, end0)``.

    ``intervals`` must be a list of ``(start, end, value)`` sorted by start.
    Returns ``None`` when no interval overlaps the region.
    """
    total = 0.0
    count = 0
    for (bs, be, value) in intervals:
        if be <= start0:
            continue
        if bs >= end0:
            break
        total += value
        count += 1
    return (total / count) if count else None


def median(values):
    """Median of a sequence of numbers (0.0 for an empty sequence)."""
    vals = sorted(values)
    n = len(vals)
    if n == 0:
        return 0.0
    if n % 2:
        return vals[n // 2]
    return (vals[n // 2 - 1] + vals[n // 2]) / 2.0


def bh_fdr(pvalues):
    """Benjamini-Hochberg FDR q-values, returned in the same order as input."""
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvalues[i])
    q = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        val = pvalues[i] * n / rank
        val = min(1.0, val)
        val = min(val, prev)  # enforce monotonicity
        prev = val
        q[i] = val
    return q


def _normal_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def mann_whitney_u(x, y):
    """Two-sided Mann-Whitney U test.

    Implemented in pure Python using the normal approximation with a tie
    correction (valid for moderate group sizes; the GEO datasets here have
    >= 3 samples per group). Returns ``(U, p_value)``.
    """
    x = [float(v) for v in x]
    y = [float(v) for v in y]
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return 0.0, 1.0

    combined = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    n = len(combined)

    # average ranks with tie handling
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and combined[j + 1][0] == combined[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = avg
        i = j + 1

    r1 = sum(ranks[k] for k in range(n) if combined[k][1] == 0)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    u2 = n1 * n2 - u1
    u = min(u1, u2)

    # tie correction term
    tie_sum = 0.0
    i = 0
    while i < n:
        j = i
        while j + 1 < n and combined[j + 1][0] == combined[i][0]:
            j += 1
        t = j - i + 1
        if t > 1:
            tie_sum += t ** 3 - t
        i = j + 1

    mu = n1 * n2 / 2.0
    var = (n1 * n2 / 12.0) * ((n + 1) - tie_sum / (n * (n - 1)))
    if var <= 0:
        return u, 1.0
    z = (u - mu) / math.sqrt(var)
    p = 2.0 * (1.0 - _normal_cdf(abs(z)))
    return u, p
