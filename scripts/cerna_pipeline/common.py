#!/usr/bin/env python3
"""Shared helpers for the ceRNA analysis pipeline (Python 3, path-parameterized).

Normalized interaction format (produced by ``parse_tools.py`` and consumed by
``predict_ceRNA.py``)::

    miRNA  target  score  energy  start  end

where ``target`` is the lncRNA/transcript ID, ``score``/``energy`` use ``-`` when
the tool does not report them, and ``start``/``end`` are the 1-based binding-site
coordinates on the target.
"""


def interaction_key(miRNA, target):
    return f"{miRNA}!{target}"


def load_interactions(path):
    """Load a normalized interaction file into
    ``{mirna!target: [(score, energy, start, end), ...]}``."""
    result = {}
    if not path:
        return result
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            miRNA, target = parts[0], parts[1]
            score = parts[2] if parts[2] != "-" else "-"
            energy = parts[3] if parts[3] != "-" else "-"
            result.setdefault(interaction_key(miRNA, target),
                              []).append((score, energy, parts[4], parts[5]))
    return result
