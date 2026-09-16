"""Lift BED/bedGraph intervals across assemblies with the UCSC kent liftOver binary."""
import subprocess

from .io import open_text
from .tooling import require_tool


def count_records(path: str) -> int:
    """Count non-comment, non-track lines in a BED-family text file."""
    records = 0
    with open_text(path) as (fh, _compression):
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("track ") or line.startswith("browser "):
                continue
            records += 1
    return records


def run_liftover(input_track: str, chain: str, output_bed: str, unmapped_bed: str, min_match: float = 0.95):
    """Run kent liftOver and return the metrics used in the report."""
    lift = require_tool("liftOver")
    command = [lift]
    if min_match is not None:
        command.append(f"-minMatch={min_match}")
    command += [input_track, chain, output_bed, unmapped_bed]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()
        raise SystemExit(f"liftOver failed ({result.returncode}): {detail}")

    input_records = count_records(input_track)
    mapped = count_records(output_bed)
    unmapped = count_records(unmapped_bed)
    rate = (mapped / input_records) if input_records else 0.0
    return {
        "input_records": input_records,
        "mapped_records": mapped,
        "unmapped_records": unmapped,
        "mapping_rate": f"{rate:.4f}",
    }


def write_report(report_path: str, input_track: str, chain: str, min_match: float, metrics: dict):
    with open(report_path, "w", encoding="utf-8") as report:
        report.write("metric\tvalue\n")
        report.write(f"input_track\t{input_track}\n")
        report.write(f"chain\t{chain}\n")
        report.write(f"min_match\t{min_match}\n")
        for key, value in metrics.items():
            report.write(f"{key}\t{value}\n")
