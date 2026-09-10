import gzip
import io
import zipfile
from contextlib import contextmanager


def detect_compression(path: str) -> str:
    lower = path.lower()
    if lower.endswith(".gz"):
        return "gzip"
    if lower.endswith(".zip"):
        return "zip"
    try:
        with open(path, "rb") as fh:
            signature = fh.read(4)
    except OSError:
        return "plain"
    if signature[:2] == b"\x1f\x8b":
        return "gzip"
    if signature == b"PK\x03\x04":
        return "zip"
    return "plain"


@contextmanager
def open_text(path: str, preferred_exts=()):
    compression = detect_compression(path)
    if compression == "gzip":
        fh = gzip.open(path, "rt", encoding="utf-8", errors="replace")
        try:
            yield fh, compression
        finally:
            fh.close()
        return

    if compression == "zip":
        zf = zipfile.ZipFile(path, "r")
        members = [name for name in zf.namelist() if not name.endswith("/")]
        if not members:
            zf.close()
            raise SystemExit(f"ZIP archive is empty: {path}")
        preferred = [
            name for name in members
            if any(name.lower().endswith(ext.lower()) for ext in preferred_exts)
        ]
        member = preferred[0] if preferred else members[0]
        raw = zf.open(member, "r")
        fh = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
        try:
            yield fh, compression
        finally:
            fh.close()
            raw.close()
            zf.close()
        return

    fh = open(path, "r", encoding="utf-8", errors="replace")
    try:
        yield fh, compression
    finally:
        fh.close()


@contextmanager
def open_text_write(path: str):
    """Open a path for text writing, gzip-compressing when it ends with .gz."""
    if path.lower().endswith(".gz"):
        fh = gzip.open(path, "wt", encoding="utf-8")
    else:
        fh = open(path, "w", encoding="utf-8")
    try:
        yield fh
    finally:
        fh.close()

