"""Exercise refusal paths using only invented-data demo commands."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit("usage: run_negative_guards.py RSCRIPT OUTPUT_TSV UPSTREAM_RUN FOUR_COHORT_RUN")
    rscript = Path(sys.argv[1]).resolve(strict=True)
    output = Path(sys.argv[2]).resolve(strict=False)
    if output.parent != (ROOT / "qa").resolve(strict=True) or output.exists():
        raise ValueError("NEW_QA_OUTPUT_REQUIRED")
    upstream_existing = Path(sys.argv[3]).resolve(strict=True)
    four = Path(sys.argv[4]).resolve(strict=True)
    run_parent = (ROOT / "runs").resolve(strict=True)
    if upstream_existing.parent != run_parent or four.parent != run_parent:
        raise ValueError("INPUT_RUN_OUTSIDE_CANDIDATE")
    upstream = ROOT / "src" / "run_upstream_synthetic.R"
    mi_lmm = ROOT / "src" / "run_mi_lmm_synthetic.R"
    rf = ROOT / "src" / "run_rf_synthetic.R"
    cases = [
        ("invalid_patient_count", [str(rscript), str(upstream), str(ROOT / "runs" / "negative-invalid-count"), "10"], ROOT / "runs" / "negative-invalid-count"),
        ("overwrite_refused", [str(rscript), str(upstream), str(upstream_existing), "200"], None),
        ("outside_run_parent", [str(rscript), str(upstream), str(ROOT / "negative-outside"), "200"], ROOT / "negative-outside"),
        ("rf_upstream_gate", [str(rscript), str(rf), str(four), str(ROOT / "runs" / "negative-rf-gate")], ROOT / "runs" / "negative-rf-gate"),
        ("mi_lmm_upstream_gate", [str(rscript), str(mi_lmm), str(four), str(ROOT / "runs" / "negative-mi-lmm-gate")], ROOT / "runs" / "negative-mi-lmm-gate"),
    ]
    results = []
    for name, command, forbidden_output in cases:
        completed = subprocess.run(command, text=True, capture_output=True, timeout=60, check=False)
        pass_case = completed.returncode != 0 and (forbidden_output is None or not forbidden_output.exists())
        results.append((name, completed.returncode, pass_case))
        if not pass_case:
            raise AssertionError("NEGATIVE_GUARD_FAILED_" + name)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("case", "exit_code", "passed"))
        writer.writerows((name, code, str(passed).upper()) for name, code, passed in results)
    print("SYNTHETIC_NEGATIVE_GUARDS_PASS", len(results))


if __name__ == "__main__":
    main()
