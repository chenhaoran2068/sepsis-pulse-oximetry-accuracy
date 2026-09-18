"""Independent, aggregate-only check of the invented standardized-event demo."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def require(condition: bool, code: str) -> None:
    if not condition:
        raise AssertionError(code)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify_standardized_event_demo.py OUTPUT_DIRECTORY")
    root = Path(sys.argv[1]).resolve(strict=True)
    require((root / "status.txt").read_text(encoding="utf-8").strip() == "STANDARDIZED_EVENT_DEMO_PASS", "RUNNER_STATUS")
    stays = read_tsv(root / "input" / "stays.tsv")
    spo2 = read_tsv(root / "input" / "spo2_events.tsv")
    sao2 = read_tsv(root / "input" / "sao2_events.tsv")
    require(len(stays) == 20 and len({r["stay_key_internal"] for r in stays}) == 20, "STAY_COUNT")
    require(len({r["patient_key_internal"] for r in stays}) == 19, "MULTI_STAY_PATIENT")
    require(len(spo2) == 81 and len(sao2) == 81, "EVENT_COUNT")
    for events in (spo2, sao2):
        require(sum(r["analysis_eligible_range"] == "TRUE" for r in events) == 80, "ELIGIBLE_EVENT_COUNT")
        require(all(r["stay_key_internal"] in {s["stay_key_internal"] for s in stays} for r in events), "ROSTER_LINK")
    checks = read_tsv(root / "checks.tsv")
    require(len(checks) >= 15 and all(r["pass"] == "TRUE" for r in checks), "NEGATIVE_AND_POSITIVE_CHECKS")
    aggregate = read_tsv(root / "aggregate.tsv")
    require({int(r["window_minutes"]) for r in aggregate} == {60, 5}, "WINDOWS")
    by_window = {int(r["window_minutes"]): r for r in aggregate}
    require(int(by_window[60]["final_pair_n"]) == 80 and int(by_window[5]["final_pair_n"]) == 60, "PAIR_COUNTS")
    require(all(int(r["stay_n"]) == 20 and int(r["patient_n"]) == 19 for r in aggregate), "DENOMINATORS")
    (root / "independent_qa_status.txt").write_text("STANDARDIZED_EVENT_INDEPENDENT_QA_PASS\n", encoding="utf-8")
    print("STANDARDIZED_EVENT_INDEPENDENT_QA_PASS")


if __name__ == "__main__":
    main()
