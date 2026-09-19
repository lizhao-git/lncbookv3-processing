#!/usr/bin/env python3
"""Generate Markdown and HTML debug reports from pipeline step reports.

Input files are the metric/value TSV reports written by the validation and
conversion modules (plus any plain header+rows TSV such as summaries or
annotation tables). The script parses each file, extracts anything that looks
like an error or a warning, and renders:

* a Markdown report (``--output-md``) suitable for issue trackers and repos;
* a self-contained HTML report (``--output-html``) with status badges, inline
  CSS and no external dependencies.

Both outputs are deterministic (no timestamps), so re-running the same
pipeline does not invalidate Nextflow task caching.
"""
import argparse
import html
import os

MAX_TABLE_ROWS = 50  # per-section cap for wide tables in the report body

STATUS_ORDER = {"FAIL": 0, "WARN": 1, "OK": 2, "INFO": 3}


def parse_tsv(path: str) -> dict:
    """Parse a report file into a normalized dict.

    metric/value reports (header ``metric\\tvalue``) become ``kind=metric``;
    everything else is treated as a generic table.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first = None
        for line in fh:
            if line.strip():
                first = line.rstrip("\n").split("\t")
                break
        if first is None:
            return {"kind": "table", "header": [], "rows": [], "row_count": 0}

        if first[:2] == ["metric", "value"] and len(first) == 2:
            metrics = []
            for line in fh:
                if not line.strip():
                    continue
                row = line.rstrip("\n").split("\t")
                if len(row) >= 2:
                    metrics.append((row[0], "\t".join(row[1:])))
                else:
                    metrics.append((row[0], ""))
            return {"kind": "metric", "metrics": metrics}

        # Generic tables (e.g. annotation tables) can hold millions of rows:
        # keep only the first MAX_TABLE_ROWS for display and count the rest.
        rows = []
        row_count = 0
        for line in fh:
            if not line.strip():
                continue
            row_count += 1
            if len(rows) < MAX_TABLE_ROWS:
                rows.append(line.rstrip("\n").split("\t"))
        return {"kind": "table", "header": first, "rows": rows, "row_count": row_count}


def issue_rows(data: dict) -> list[tuple[str, str]]:
    """Return (key, value) pairs that carry error/warning content."""
    if data["kind"] != "metric":
        return []
    issues = []
    for key, value in data["metrics"]:
        key_l = key.lower()
        if (key_l.startswith("error") or key_l.startswith("warning")) and value not in ("", "0"):
            issues.append((key, value))
    return issues


def step_status(data: dict) -> str:
    issues = issue_rows(data)
    if any(key.lower().startswith("error") for key, _ in issues):
        return "FAIL"
    if issues:
        return "WARN"
    return "OK"


def step_records(data: dict) -> str:
    if data["kind"] == "metric":
        for key in ("records", "parsed_records"):
            for mkey, value in data["metrics"]:
                if mkey == key:
                    return value
        return "-"
    return str(data["row_count"])


def md_escape(cell: str) -> str:
    return cell.replace("|", "\\|").replace("\n", "<br>")


def render_markdown(title: str, steps: list[tuple[str, dict]]) -> str:
    lines = [f"# {title}", ""]

    # ---- overview ----
    lines += ["## Overview", "", "| # | step | records | errors | warnings | status |", "| --- | --- | --- | --- | --- | --- |"]
    for idx, (name, data) in enumerate(steps, start=1):
        data_kind = data["kind"]
        if data_kind == "metric":
            get = lambda key: next((v for k, v in data["metrics"] if k == key), "0")
            errors = get("error_count")
            warnings = get("warning_count")
        else:
            errors = warnings = "-"
        lines.append(f"| {idx} | `{name}` | {step_records(data)} | {errors} | {warnings} | {step_status(data)} |")

    # ---- errors & warnings ----
    failing = [(name, issue_rows(data)) for name, data in steps if issue_rows(data)]
    lines += ["", "## Errors & Warnings", ""]
    if not failing:
        lines.append("_No errors or warnings recorded in any step._")
        lines.append("")
    else:
        for name, issues in failing:
            lines.append(f"### `{name}`")
            lines.append("")
            for key, value in issues:
                lines.append(f"- **{key}**: {md_escape(value)}")
            lines.append("")

    # ---- step details ----
    lines += ["## Step details", ""]
    for name, data in steps:
        lines.append(f"### `{name}`")
        lines.append("")
        if data["kind"] == "metric":
            lines += ["| metric | value |", "| --- | --- |"]
            for key, value in data["metrics"]:
                lines.append(f"| {md_escape(key)} | {md_escape(value)} |")
        else:
            header = data["header"]
            rows = data["rows"]
            truncated = data["row_count"] > MAX_TABLE_ROWS
            shown = rows[:MAX_TABLE_ROWS]
            if header:
                lines += ["| " + " | ".join(md_escape(c) for c in header) + " |",
                          "|" + " --- |" * len(header)]
                for row in shown:
                    padded = row + [""] * (len(header) - len(row))
                    lines.append("| " + " | ".join(md_escape(c) for c in padded[: len(header)]) + " |")
            if truncated:
                lines.append("")
                lines.append(f"_... {data['row_count'] - MAX_TABLE_ROWS} more rows omitted._")
        lines.append("")

    return "\n".join(lines) + "\n"


def render_html(title: str, steps: list[tuple[str, dict]]) -> str:
    esc = html.escape

    badge = {
        "OK": '<span class="badge ok">OK</span>',
        "WARN": '<span class="badge warn">WARN</span>',
        "FAIL": '<span class="badge fail">FAIL</span>',
        "INFO": '<span class="badge info">INFO</span>',
    }

    def metric_table(data: dict) -> str:
        rows_html = []
        for key, value in data["metrics"]:
            cls = ""
            key_l = key.lower()
            if key_l.startswith("error") and value not in ("", "0"):
                cls = ' class="err"'
            elif key_l.startswith("warning") and value not in ("", "0"):
                cls = ' class="warn"'
            rows_html.append(f"<tr><td>{esc(key)}</td><td{cls}>{esc(value)}</td></tr>")
        return '<table class="kv"><tr><th>metric</th><th>value</th></tr>' + "".join(rows_html) + "</table>"

    def generic_table(data: dict) -> str:
        header = data["header"]
        rows = data["rows"]
        if not header:
            return "<p><em>(empty file)</em></p>"
        thead = "<tr>" + "".join(f"<th>{esc(c)}</th>" for c in header) + "</tr>"
        shown = rows[:MAX_TABLE_ROWS]
        body = "".join(
            "<tr>" + "".join(f"<td>{esc(c)}</td>" for c in (row + [""] * (len(header) - len(row)))[: len(header)]) + "</tr>"
            for row in shown
        )
        table = f'<table class="data">{thead}{body}</table>'
        if data["row_count"] > MAX_TABLE_ROWS:
            table += f"<p><em>... {data['row_count'] - MAX_TABLE_ROWS} more rows omitted.</em></p>"
        return table

    # ---- overview ----
    overview_rows = []
    for idx, (name, data) in enumerate(steps, start=1):
        status = step_status(data)
        if data["kind"] == "metric":
            get = lambda key: next((v for k, v in data["metrics"] if k == key), "0")
            errors = esc(get("error_count"))
            warnings = esc(get("warning_count"))
        else:
            errors = warnings = "-"
        overview_rows.append(
            f"<tr><td>{idx}</td><td><code>{esc(name)}</code></td>"
            f"<td>{esc(step_records(data))}</td><td>{errors}</td><td>{warnings}</td>"
            f"<td>{badge[status]}</td></tr>"
        )

    # ---- errors & warnings ----
    failing = [(name, issue_rows(data)) for name, data in steps if issue_rows(data)]
    if failing:
        issue_html = "".join(
            f"<h3><code>{esc(name)}</code></h3><ul>"
            + "".join(f"<li><strong>{esc(key)}</strong>: <code>{esc(value)}</code></li>" for key, value in issues)
            + "</ul>"
            for name, issues in failing
        )
    else:
        issue_html = "<p><em>No errors or warnings recorded in any step.</em></p>"

    # ---- step details ----
    details_html = "".join(
        f"<h2><code>{esc(name)}</code> {badge[step_status(data)]}</h2>"
        + (metric_table(data) if data["kind"] == "metric" else generic_table(data))
        for name, data in steps
    )

    total_fail = sum(1 for _, d in steps if step_status(d) == "FAIL")
    total_warn = sum(1 for _, d in steps if step_status(d) == "WARN")
    if total_fail:
        summary = f"{len(steps)} steps, {total_fail} failing, {total_warn} with warnings"
    elif total_warn:
        summary = f"{len(steps)} steps, no failures, {total_warn} with warnings"
    else:
        summary = f"{len(steps)} steps, all clean"

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{esc(title)}</title>
<style>
  :root {{ color-scheme: light dark; }}
  body {{ font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
         margin: 2rem auto; max-width: 60rem; padding: 0 1rem; line-height: 1.5; }}
  h1 {{ border-bottom: 2px solid #8884; padding-bottom: .3rem; }}
  h2 {{ border-bottom: 1px solid #8884; padding-bottom: .2rem; margin-top: 2rem; }}
  h3 {{ margin-bottom: .3rem; }}
  code {{ background: #88888822; padding: .05rem .3rem; border-radius: 3px; }}
  table {{ border-collapse: collapse; margin: .6rem 0 1.2rem; width: 100%; font-size: .92rem; }}
  th, td {{ border: 1px solid #8885; padding: .3rem .55rem; text-align: left; vertical-align: top;
           word-break: break-word; }}
  th {{ background: #88888822; }}
  tr:hover td {{ background: #88888811; }}
  td.err {{ color: #b30000; font-weight: 600; }}
  td.warn {{ color: #9a6700; }}
  .badge {{ display: inline-block; padding: .05rem .5rem; border-radius: 999px;
           font-size: .78rem; font-weight: 700; letter-spacing: .03em; }}
  .badge.ok {{ background: #1a7f37; color: #fff; }}
  .badge.warn {{ background: #9a6700; color: #fff; }}
  .badge.fail {{ background: #b30000; color: #fff; }}
  .badge.info {{ background: #57606a; color: #fff; }}
  .summary {{ color: #57606a; }}
</style>
</head>
<body>
<h1>{esc(title)}</h1>
<p class="summary">{summary}</p>
<h2>Overview</h2>
<table>
<tr><th>#</th><th>step</th><th>records</th><th>errors</th><th>warnings</th><th>status</th></tr>
{''.join(overview_rows)}
</table>
<h2>Errors &amp; Warnings</h2>
{issue_html}
<h2>Step details</h2>
{details_html}
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description="Generate Markdown and HTML debug reports from pipeline step reports")
    parser.add_argument("--input", action="append", required=True,
                        help="Report file to include; repeat for each step (order is preserved)")
    parser.add_argument("--output-md", required=True)
    parser.add_argument("--output-html", required=True)
    parser.add_argument("--title", default="Pipeline debug report")
    args = parser.parse_args()

    steps = []
    for path in args.input:
        name = os.path.basename(path)
        try:
            data = parse_tsv(path)
        except OSError as exc:
            data = {"kind": "metric", "metrics": [("unreadable", str(exc))]}
        steps.append((name, data))

    # FAIL first, then WARN, then OK — makes problems visible at the top
    steps.sort(key=lambda item: (STATUS_ORDER.get(step_status(item[1]), 9), item[0]))

    with open(args.output_md, "w", encoding="utf-8") as out:
        out.write(render_markdown(args.title, steps))
    with open(args.output_html, "w", encoding="utf-8") as out:
        out.write(render_html(args.title, steps))


if __name__ == "__main__":
    main()
