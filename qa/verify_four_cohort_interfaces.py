"""Independent aggregate-only checks of the four synthetic cohort interfaces."""

import csv
import hashlib
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
EXPECTED_FILES = {"aggregate.tsv", "checks.tsv", "source_receipt.tsv", "status.txt"}
EXPECTED_COHORTS = {"Amsterdam", "eICU", "SICDB", "Lianyungang"}
EXPECTED_CHECKS = {
    "bridge_unique_and_cohort_scoped", "no_raw_fields", "duplicate_and_range",
    "artifact_mode", "invalid_staged_schema_rejected", "independent_pair_windows",
    "range_excluded", "one_to_one_60", "one_to_one_5",
    "pairing_deterministic", "cross_modality_patient_conflict_rejected",
}


def rows(run, name):
    with (run / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify_four_cohort_interfaces.py RUN_DIRECTORY")
    run = Path(sys.argv[1]).resolve(strict=True)
    if run.parent != (ROOT / "runs").resolve(strict=True):
        raise AssertionError("RUN_OUTSIDE_CANDIDATE")
    if {p.name for p in run.iterdir()} != EXPECTED_FILES:
        raise AssertionError("UNEXPECTED_OR_MISSING_OUTPUT_FILE")
    if (run / "status.txt").read_text(encoding="utf-8").strip() != "FOUR_COHORT_SYNTHETIC_PASS_PENDING_INDEPENDENT_QA":
        raise AssertionError("RUN_STATUS_NOT_PASS")

    aggregate = rows(run, "aggregate.tsv")
    if len(aggregate) != 4 or {x["cohort"] for x in aggregate} != EXPECTED_COHORTS:
        raise AssertionError("AGGREGATE_COHORT_SET")
    for row in aggregate:
        for field, expected in {"synthetic_stays": 2, "eligible_spo2_events": 6,
                                "eligible_sao2_events": 6, "pairs_60_min": 6,
                                "pairs_5_min": 4}.items():
            if int(row[field]) != expected:
                raise AssertionError("AGGREGATE_COUNT_" + row["cohort"] + "_" + field)
        expected_mode = "NOT_ASSESSABLE_TEMPORAL_RESOLUTION" if row["cohort"] == "SICDB" else "P99_15_MIN_THREE_POINT"
        if row["artifact_qc_mode"] != expected_mode:
            raise AssertionError("ARTIFACT_MODE_" + row["cohort"])

    checks = rows(run, "checks.tsv")
    if not checks or any(row["pass"] != "TRUE" for row in checks):
        raise AssertionError("FAILED_SYNTHETIC_CHECK")
    by_cohort = {}
    for row in checks:
        by_cohort.setdefault(row["cohort"], set()).add(row["check"])
    if set(by_cohort) != EXPECTED_COHORTS | {"ALL"}:
        raise AssertionError("CHECK_COHORT_SET")
    for cohort in EXPECTED_COHORTS:
        if not EXPECTED_CHECKS <= by_cohort[cohort]:
            raise AssertionError("MISSING_CHECK_" + cohort)
    if {"fraction_normalized_only_for_item12311", "explicit_venous_comments_detected"} - by_cohort["Amsterdam"]:
        raise AssertionError("AMSTERDAM_SOURCE_CHECKS")
    if "icu_relative_minute_interface" not in by_cohort["eICU"]:
        raise AssertionError("EICU_SOURCE_CHECK")
    if "seconds_to_icu_minutes" not in by_cohort["SICDB"]:
        raise AssertionError("SICDB_SOURCE_CHECK")
    if {"label_rules", "source_time_and_value_parse"} - by_cohort["Lianyungang"]:
        raise AssertionError("LIANYUNGANG_SOURCE_CHECKS")
    if by_cohort["ALL"] != {"distinct_internal_stay_keys_between_cohorts"}:
        raise AssertionError("CROSS_COHORT_KEY_CHECK")

    module_root = ROOT / "code" / "cores"
    receipts = rows(run, "source_receipt.tsv")
    if {r["name"] for r in receipts} != {"R04", "R05"}:
        raise AssertionError("SOURCE_RECEIPT_MODULES")
    for row in receipts:
        path = module_root / row["relative_path"]
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest().upper() != row["expected_sha256"] or row["observed_sha256"] != row["expected_sha256"]:
            raise AssertionError("SOURCE_HASH_MISMATCH_" + row["name"])

    source = ROOT / "src" / "run_four_cohort_interfaces_synthetic.R"
    if not source.is_file():
        raise AssertionError("TEST_SOURCE_MISSING")
    (run / "independent_qa_status.txt").write_text(
        "FOUR_COHORT_INTERFACES_INDEPENDENT_QA_PASS\n", encoding="utf-8")
    with (run / "independent_qa_receipt.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("component", "sha256"))
        for name, path in (("test_script", source), ("aggregate", run / "aggregate.tsv"),
                           ("checks", run / "checks.tsv"), ("source_receipt", run / "source_receipt.tsv")):
            writer.writerow((name, hashlib.sha256(path.read_bytes()).hexdigest().upper()))
    print("FOUR_COHORT_INTERFACES_INDEPENDENT_QA_PASS", len(checks))


if __name__ == "__main__":
    main()
