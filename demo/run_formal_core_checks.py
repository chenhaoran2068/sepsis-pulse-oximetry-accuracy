"""Exercise the copied formal scientific cores on invented data only.

This is a portability check, not the paper-analysis entrypoint.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import uuid
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def assert_all_pass(path: Path, expected: str) -> None:
    if not path.is_file():
        raise RuntimeError(f"FORMAL_CORE_QA_MISSING: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise RuntimeError(f"FORMAL_CORE_QA_EMPTY: {path}")
    for row in rows:
        result = row.get("pass", row.get("status", ""))
        if result not in {"TRUE", "PASS"}:
            raise RuntimeError(f"FORMAL_CORE_QA_FAIL: {path}: {row}")
    print(f"{expected}_QA_PASS checks={len(rows)}", flush=True)


def run_one(rscript: str, script: Path, output: Path, log: Path, env: dict[str, str]) -> None:
    with log.open("w", encoding="utf-8") as stream:
        result = subprocess.run([rscript, str(script), str(output)], cwd=ROOT,
                                env=env, stdout=stream, stderr=subprocess.STDOUT,
                                check=False)
    if result.returncode != 0:
        raise RuntimeError(f"FORMAL_CORE_TEST_FAILED: {script.name}; log={log}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--run-id")
    parser.add_argument("--extra-r-library", help="optional existing R package library")
    args = parser.parse_args()
    rscript = shutil.which(args.rscript) or (args.rscript if Path(args.rscript).is_file() else None)
    if rscript is None:
        raise SystemExit("FORMAL_CORE_RSCRIPT_NOT_FOUND")
    run_id = args.run_id or uuid.uuid4().hex[:12]
    if re.fullmatch(r"[a-z0-9-]{1,32}", run_id) is None:
        raise SystemExit("FORMAL_CORE_RUN_ID_INVALID")
    run_root = ROOT / "runs" / f"formal-core-{run_id}"
    if run_root.exists():
        raise SystemExit(f"FORMAL_CORE_RUN_ID_EXISTS: {run_id}")
    run_root.mkdir(parents=True)
    env = os.environ.copy()
    if args.extra_r_library:
        library = Path(args.extra_r_library).resolve(strict=True)
        if not library.is_dir():
            raise SystemExit("FORMAL_CORE_LIBRARY_NOT_DIRECTORY")
        env["R9_EXTRA_R_LIB"] = str(library)
    targets = (
        ("mi", ROOT / "tests" / "formal_mi_30x20_synthetic.R",
         Path("aggregate/synthetic_formal_checks.tsv")),
        ("lmm", ROOT / "tests" / "formal_lmm_synthetic.R",
         Path("aggregate/synthetic_checks.tsv")),
        ("rf", ROOT / "tests" / "formal_rf_synthetic.R",
         Path("end_to_end_synthetic_results.tsv")),
    )
    try:
        for name, script, checks in targets:
            output = run_root / name
            run_one(rscript, script, output, run_root / f"{name}.log", env)
            assert_all_pass(output / checks, name.upper())
    except Exception as exc:
        raise SystemExit(f"{exc}; run={run_root}") from exc
    print(f"FORMAL_CORE_SYNTHETIC_CHECKS_PASS run={run_root}")
    print("NOT_A_STANDARDIZED_INPUT_OR_PAPER_RESULT_REPRODUCTION_RUN")


if __name__ == "__main__":
    main()
