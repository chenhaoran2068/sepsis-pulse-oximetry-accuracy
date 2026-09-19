"""Create a candidate-only file/hash inventory and basic disclosure screen.

This is a screening aid, not proof that source code is licensed for release or
that tabular example records are nonclinical. The owner must review its output.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXCLUDE_DIRS = {"runs", ".git", "__pycache__"}
ALLOWED_SUFFIXES = {"", ".R", ".py", ".md", ".tsv", ".cff", ".gitignore"}
SENSITIVE_PATTERNS = {
    "absolute_windows_path": re.compile(r"(?i)\b[a-z]:[\\/]"),
    "email_address": re.compile(r"[\w.+-]+@[\w.-]+\.[a-zA-Z]{2,}"),
    "github_token": re.compile(r"(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,})"),
    "api_key": re.compile(r"\bsk-[A-Za-z0-9_-]{20,}\b"),
    "private_key": re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path, help="new TSV outside the candidate tree")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"inventory exists; refusing overwrite: {output}")
    if ROOT in output.parents or output == ROOT:
        parser.error("inventory must be outside the candidate tree")
    files = sorted(
        (path for path in ROOT.rglob("*") if path.is_file()
         and not any(part in EXCLUDE_DIRS for part in path.relative_to(ROOT).parts)),
        key=lambda path: path.relative_to(ROOT).as_posix().lower(),
    )
    if not files:
        parser.error("no candidate files")
    records = []
    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        if path.suffix not in ALLOWED_SUFFIXES and path.name != ".gitignore":
            parser.error(f"unexpected public-file extension: {relative}")
        content = path.read_bytes()
        if b"\0" in content:
            parser.error(f"binary public file: {relative}")
        text = content.decode("utf-8")
        flags = [name for name, pattern in SENSITIVE_PATTERNS.items() if pattern.search(text)]
        records.append({
            "relative_path": relative,
            "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest(),
            "disclosure_screen_flags": ",".join(flags) if flags else "none",
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print("INVENTORY_CREATED", output)
    print("CANDIDATE_FILES", len(records))
    print("DISCLOSURE_FLAGGED_FILES", sum(row["disclosure_screen_flags"] != "none" for row in records))


if __name__ == "__main__":
    main()
