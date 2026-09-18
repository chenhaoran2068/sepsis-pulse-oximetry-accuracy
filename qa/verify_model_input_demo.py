"""Independent aggregate-only checks for the invented model-input interface."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


EXPECTED_CHECKS = {
    "lmm_anchor_60", "lmm_structural_60", "lmm_stay_age_60",
    "coverage_ledger_60", "rf_pair_identity", "rf_no_time_or_source_predictors",
    "rf_formula_identity", "rf_5_minute_rejected", "lmm_anchor_5",
    "lmm_structural_5", "lmm_stay_age_5", "coverage_ledger_5",
    "rf_sensitivity_window_rejected",
}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify_model_input_demo.py RUN_DIR")
    root = Path(__file__).resolve().parents[1]
    runs = root / "runs"
    target = Path(sys.argv[1]).resolve(strict=True)
    if target.parent != runs.resolve(strict=True):
        raise ValueError("RUN_PATH_OUTSIDE_CANDIDATE")
    if (target / "status.txt").read_text(encoding="utf-8").strip() != "MODEL_INPUT_DEMO_PASS":
        raise ValueError("MODEL_INPUT_STATUS_NOT_PASS")
    aggregate = rows(target / "aggregate.tsv")
    if len(aggregate) != 2 or {int(r["window_minutes"]) for r in aggregate} != {5, 60}:
        raise ValueError("WINDOWS_INCORRECT")
    if {int(r["window_minutes"]): int(r["pair_n"]) for r in aggregate} != {5: 60, 60: 80}:
        raise ValueError("PAIR_COUNTS_INCORRECT")
    if any(int(r["patient_n"]) != 19 or int(r["lmm_field_n"]) < 25 for r in aggregate):
        raise ValueError("MODEL_INPUT_SHAPE_INCORRECT")
    checks = rows(target / "checks.tsv")
    if {r["check"] for r in checks} != EXPECTED_CHECKS or any(r["pass"] != "TRUE" for r in checks):
        raise ValueError("CHECK_SET_INCORRECT")
    (target / "independent_qa_status.txt").write_text("MODEL_INPUT_INDEPENDENT_QA_PASS\n", encoding="utf-8")
    print("MODEL_INPUT_INDEPENDENT_QA_PASS")


if __name__ == "__main__":
    main()
