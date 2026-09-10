"""Helpers for invoking external command-line tools (stdlib only)."""
import shutil
import subprocess


def require_tool(name: str) -> str:
    """Return the absolute path of an external tool or abort with guidance."""
    path = shutil.which(name)
    if not path:
        raise SystemExit(
            f"Required tool '{name}' was not found on PATH. "
            "Run inside the lncbookv3-processing Docker image or install it locally."
        )
    return path


def run_command(command, **kwargs):
    """Run a command, aborting on non-zero exit status."""
    subprocess.run(command, check=True, **kwargs)
