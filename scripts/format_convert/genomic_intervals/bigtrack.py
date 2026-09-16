"""Validate bigWig/bigBed tracks through the kent info tools (stdlib only)."""
import os
import shutil
import subprocess

from .tooling import require_tool

BIGBED_EXTENSIONS = (".bigbed", ".bb")
BIGWIG_EXTENSIONS = (".bigwig", ".bw")

# Metrics worth surfacing from bigWigInfo / bigBedInfo output.
BIGWIG_METRICS = ("version", "chromCount", "basesCovered", "mean", "min", "max")
BIGBED_METRICS = ("version", "chromCount", "itemCount", "validItemCount", "basesCovered")


def detect_info_tool(path: str):
    """Return (info tool, track type) based on the file extension."""
    lower = path.lower()
    if lower.endswith(BIGBED_EXTENSIONS):
        return "bigBedInfo", "bigBed"
    if lower.endswith(BIGWIG_EXTENSIONS):
        return "bigWigInfo", "bigWig"
    raise SystemExit(
        f"Cannot infer track type from {path!r}; "
        "expected a .bigBed/.bb or .bigWig/.bw file"
    )


def parse_info(output: str) -> dict:
    fields = {}
    for line in output.splitlines():
        key, sep, value = line.partition(":")
        if sep:
            fields[key.strip()] = value.strip()
    return fields


def validate_bigtrack(input_path: str, output_path: str, report_path: str):
    if os.path.realpath(input_path) == os.path.realpath(output_path):
        raise SystemExit("Input and output track paths must differ")

    tool_name, track_type = detect_info_tool(input_path)
    tool = require_tool(tool_name)
    result = subprocess.run([tool, input_path], capture_output=True, text=True)
    fields = parse_info(result.stdout)

    errors = []
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip().splitlines()
        errors.append(detail[-1] if detail else f"{tool_name} exited with {result.returncode}")

    wanted = BIGBED_METRICS if track_type == "bigBed" else BIGWIG_METRICS
    with open(report_path, "w", encoding="utf-8") as report:
        report.write("metric\tvalue\n")
        report.write(f"track_type\t{track_type}\n")
        report.write(f"info_tool\t{tool_name}\n")
        for key in wanted:
            if key in fields:
                report.write(f"{key}\t{fields[key]}\n")
        report.write(f"error_count\t{len(errors)}\n")
        if errors:
            report.write("errors\t" + " | ".join(errors[:100]) + "\n")

    if errors:
        raise SystemExit("Big track validation failed. See report for details.")

    shutil.copyfile(input_path, output_path)
