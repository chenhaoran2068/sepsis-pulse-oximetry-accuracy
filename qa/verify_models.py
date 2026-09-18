"""Independent aggregate-only interface checks for private synthetic model runs."""

import csv
import math
import pathlib
import sys


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def require(condition, reason):
    if not condition:
        raise AssertionError(reason)


def verify_receipts(root, expected_count):
    receipts = rows(root / "source_receipt.tsv")
    require(len(receipts) == expected_count, "module receipt count")
    require(all(r["expected_sha256"] == r["observed_sha256"] for r in receipts), "source hash mismatch")


def verify_common(root, kind, receipt_count, check_count):
    require((root / "status.txt").read_text(encoding="utf-8").strip() == f"SYNTHETIC_{kind}_PASS", "run status")
    verify_receipts(root, receipt_count)
    checks = rows(root / "checks.tsv")
    require(len(checks) == check_count and all(r["pass"] == "TRUE" for r in checks), "runner checks")
    aggregate = rows(root / "aggregate.tsv")
    require(len(aggregate) == 1, "one aggregate row")
    allowed = {"aggregate.tsv", "checks.tsv", "source_receipt.tsv", "status.txt", "independent_qa_status.txt"}
    if kind == "MI_LMM":
        allowed.add("imputation_logged_event_summary.tsv")
    require({p.name for p in root.iterdir()} <= allowed, "unexpected model-run file")
    return aggregate[0]


def verify_mi_lmm(root):
    row = verify_common(root, "MI_LMM", 3, 10)
    require((int(row["synthetic_patient_n"]), int(row["synthetic_pair_n"])) == (60, 240), "MI fixture cardinality")
    require((int(row["pilot_imputations"]), int(row["pilot_iterations"])) == (3, 5), "pilot design")
    require(int(row["imputed_cell_n"]) > 0 and int(row["imputation_domain_fail_n"]) == 0, "imputed domains")
    require(int(row["imputation_logged_event_n"]) == 0, "pilot MI logged events")
    require(int(row["imputation_internal_singular_message_n"]) == 0, "intermediate MI singular message")
    require(int(row["lmm_singular_fit_n"]) == 0, "final LMM singular fit")
    require(rows(root / "imputation_logged_event_summary.tsv") == [], "MI logged event summary")
    require(int(row["lmm_fit_n"]) == 3 and int(row["lmm_term_n"]) > 0, "LMM fit count")
    require(row["primary_all_finite"] == row["balanced_all_finite"] == "TRUE", "pooled coefficients")


def verify_rf(root):
    row = verify_common(root, "RF", 5, 13)
    require((int(row["synthetic_patient_n"]), int(row["synthetic_pair_n"])) == (200, 800), "RF fixture cardinality")
    require((int(row["development_patient_n"]), int(row["development_pair_n"])) == (160, 640), "development split")
    require((int(row["test_patient_n"]), int(row["test_pair_n"])) == (40, 160), "test split")
    require(int(row["inner_task_n"]) == 60, "inner CV task count")
    require(int(row["outer_retained_clinical_n"]) > 0 and int(row["predictor_n"]) > 1, "training recipe")
    require(1 <= int(row["selected_grid_index"]) <= 12, "selected tuning rule")
    require(math.isfinite(float(row["test_rmse"])) and math.isfinite(float(row["test_mae"])), "test metrics")
    require(int(row["bootstrap_success_n"]) == 50, "bootstrap")
    require((int(row["shap_patient_n"]), int(row["shap_group_n"])) == (15, 10), "SHAP summary")


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in {"mi_lmm", "rf"}:
        raise SystemExit("usage: verify_models.py mi_lmm|rf OUTPUT_DIRECTORY")
    root = pathlib.Path(sys.argv[2]).resolve(strict=True)
    if sys.argv[1] == "mi_lmm":
        verify_mi_lmm(root)
        label = "SYNTHETIC_MI_LMM_INDEPENDENT_QA_PASS"
    else:
        verify_rf(root)
        label = "SYNTHETIC_RF_INDEPENDENT_QA_PASS"
    (root / "independent_qa_status.txt").write_text(label + "\n", encoding="utf-8")
    print(label)


if __name__ == "__main__":
    main()
