"""Independently check the invented clinical-input adapter demonstration."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValueError("QA_FILE_ABSENT")
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify_clinical_input_demo.py RUN_DIRECTORY")
    run = Path(sys.argv[1]).resolve(strict=True)
    if (run / "independent_qa_status.txt").exists():
        raise ValueError("QA_REFUSE_OVERWRITE")
    if (run / "status.txt").read_text(encoding="utf-8").strip() != "CLINICAL_INPUT_DEMO_PASS":
        raise ValueError("RUNNER_NOT_PASS")
    aggregate = rows(run / "aggregate.tsv")
    expected = {"60": 80, "5": 60}
    if len(aggregate) != 2 or {r["window_minutes"] for r in aggregate} != set(expected):
        raise ValueError("AGGREGATE_WINDOWS_INVALID")
    for row in aggregate:
        if int(row["stay_n"]) != 20 or int(row["patient_n"]) != 19:
            raise ValueError("AGGREGATE_ROSTER_INVALID")
        if int(row["analytic_pair_n"]) != expected[row["window_minutes"]]:
            raise ValueError("AGGREGATE_PAIR_COUNT_INVALID")
    checks = rows(run / "checks.tsv")
    required = {
        "stay_row_cardinality", "independent_pair_cardinality", "outcome_and_background",
        "stay_level_value_consistency", "wrong_stay_mapping_rejected", "missing_pair_rejected",
        "duplicate_pair_uid_rejected", "unit_mismatch_rejected",
        "state_value_mismatch_rejected", "missing_input_rejected", "window_mismatch_rejected",
    }
    if {r["check"] for r in checks} != required or any(r["pass"] != "TRUE" for r in checks):
        raise ValueError("RUNNER_CHECKS_INVALID")
    input_dir = run / "input"
    if len(rows(input_dir / "stay_clinical.tsv")) != 20:
        raise ValueError("STAY_INPUT_COUNT_INVALID")
    for window, count in expected.items():
        clinical = rows(input_dir / f"pair_clinical_{window}.tsv")
        if len(clinical) != count:
            raise ValueError("PAIR_INPUT_COUNT_INVALID")
        if len({r["pair_key_internal"] for r in clinical}) != count:
            raise ValueError("PAIR_UID_DUPLICATE")
    (run / "independent_qa_status.txt").write_text(
        "CLINICAL_INPUT_INDEPENDENT_QA_PASS\n", encoding="utf-8"
    )
    print("CLINICAL_INPUT_INDEPENDENT_QA_PASS")


if __name__ == "__main__":
    main()
