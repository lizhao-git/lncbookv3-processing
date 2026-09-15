import csv
import json
import os
from collections import OrderedDict


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def resolve_path(path, data_root=None, manifest_dir=None):
    if not path:
        return None
    if os.path.isabs(path):
        return path
    for base in (data_root, manifest_dir):
        if base:
            candidate = os.path.join(base, path)
            if os.path.exists(candidate):
                return candidate
    return os.path.join(data_root or manifest_dir or ".", path)


def load_manifest(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def read_tsv(path, header=True, names=None):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return []
    if header:
        fields = rows[0]
        data_rows = rows[1:]
    else:
        fields = names or [f"V{i + 1}" for i in range(len(rows[0]))]
        data_rows = rows
    out = []
    for row in data_rows:
        padded = row + [""] * max(0, len(fields) - len(row))
        out.append(dict(zip(fields, padded[:len(fields)])))
    return out


def write_tsv(path, rows, fields=None):
    ensure_dir(os.path.dirname(path) or ".")
    if fields is None:
        fields = []
        seen = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    fields.append(key)
                    seen.add(key)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def read_fasta(path):
    records = OrderedDict()
    name = None
    seq = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        records[name] = "".join(seq)
    return records


def prefixed(row, prefix, skip=()):
    skip = set(skip)
    return {f"{prefix}{key}": value for key, value in row.items() if key not in skip}


def parse_number(value, default=0.0):
    try:
        if value in ("", None):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value, default=0):
    try:
        if value in ("", None):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default

