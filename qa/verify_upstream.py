"""Independent aggregate-only checks for the private synthetic R04–R05–403–404 run."""

import csv
import pathlib
import sys


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: verify_upstream.py OUTPUT_DIRECTORY")
    root = pathlib.Path(sys.argv[1]).resolve(strict=True)
    check((root / "status.txt").read_text(encoding="utf-8").strip() == "SYNTHETIC_UPSTREAM_PASS", "runner status")
    receipts = rows(root / "source_receipt.tsv")
    check(len(receipts) == 6, "six source receipts")
    check(all(r["expected_sha256"] == r["observed_sha256"] for r in receipts), "source hashes")
    checks = rows(root / "checks.tsv")
    check(len(checks) >= 20 and all(r["pass"] == "TRUE" for r in checks), "runner checks")
    spec = rows(root / "fixture_spec.tsv")
    check(len(spec) == 1, "fixture specification")
    patient_n = int(spec[0]["synthetic_patient_n"])
    check(20 <= patient_n <= 500, "synthetic patient range")
    check((root / "synthetic_pairs_INTERNAL_ONLY.rds").is_file(), "synthetic downstream input")
    aggs = rows(root / "aggregate.tsv")
    check({r["window"] for r in aggs} == {"M60", "M5"}, "both independent windows")
    for row in aggs:
        n = int(row["final_pair_n"])
        expected_pairs = int(spec[0]["synthetic_pairs_per_patient_" + row["window"].lower()]) * patient_n
        check(n == expected_pairs, "expected pair count")
        check(int(row["patient_n"]) == patient_n, "patient count")
        check(int(row["low_sao2_n"]) == patient_n, "low SaO2 count")
        check(0 <= int(row["severe_occult_n"]) <= int(row["occult_n"]) <= patient_n, "forward subset")
        check(0 <= int(row["reverse88_numerator_n"]) <= int(row["reverse88_denominator_n"]) <= n, "reverse subset")
        check(0 <= int(row["affected_patient_n"]) <= patient_n, "patient subset")
        lower, mean, upper = (float(row[k]) for k in ("lower_loa", "mean_bias", "upper_loa"))
        check(lower < mean < upper and float(row["arms"]) >= 0, "agreement bounds")
    (root / "independent_qa_status.txt").write_text("SYNTHETIC_UPSTREAM_INDEPENDENT_QA_PASS\n", encoding="utf-8")
    print("SYNTHETIC_UPSTREAM_INDEPENDENT_QA_PASS")


if __name__ == "__main__":
    main()
