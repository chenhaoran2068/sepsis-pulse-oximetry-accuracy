"""Independent QA of a 30-display synthetic technical bundle."""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DISPLAY_IDS = (
    *(f"main-table-{n}" for n in (1, 2)),
    *(f"main-figure-{n}" for n in (1, 2, 3)),
    *(f"supplementary-table-{n}" for n in range(10, 25)),
    *(f"supplementary-figure-{n}" for n in range(2, 12)),
)


def run_r(rscript: str, script: str, *arguments: Path | str) -> None:
    command = [rscript, str(ROOT / "tests" / script), *(str(x) for x in arguments)]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"{script} failed:\n{result.stdout}\n{result.stderr}")
    print(result.stdout.strip().splitlines()[-1])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--bundle", required=True, type=Path)
    p.add_argument("--run-map", required=True, type=Path)
    p.add_argument("--rscript", default="Rscript")
    args = p.parse_args()
    bundle = args.bundle.resolve()
    rscript = shutil.which(args.rscript) or (
        args.rscript if Path(args.rscript).is_file() else None
    )
    if rscript is None:
        p.error("Rscript unavailable")
    if (bundle / "status.txt").read_text(encoding="utf-8").strip() != "TECHNICAL_DISPLAY_CHAIN_PASS_30_OF_40":
        raise RuntimeError("BUNDLE_STATUS_NOT_PASS")
    with (bundle / "display_inventory.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 30 or {row["display"] for row in rows} != set(DISPLAY_IDS):
        raise RuntimeError("BUNDLE_DISPLAY_SET_INVALID")
    for row in rows:
        label = row["display"]
        pdf = ROOT / row["pdf_path"]
        if pdf.resolve() != (bundle / f"{label}.pdf").resolve():
            raise RuntimeError(f"BUNDLE_PATH_INVALID {label}")
        if not pdf.is_file() or pdf.stat().st_size < 1000 or not pdf.read_bytes().startswith(b"%PDF"):
            raise RuntimeError(f"BUNDLE_PDF_INVALID {label}")
        if hashlib.sha256(pdf.read_bytes()).hexdigest() != row["sha256"]:
            raise RuntimeError(f"BUNDLE_PDF_HASH_MISMATCH {label}")
        if not (bundle / f"{label}-receipt.tsv").is_file():
            raise RuntimeError(f"BUNDLE_RECEIPT_MISSING {label}")
    run_r(rscript,"verify_display_run_map.R",args.run_map.resolve())
    def f(name: str) -> Path:
        return bundle / name
    for label,w in (("main-table-1",60),("supplementary-table-11",5)):
        run_r(rscript,"verify_agreement_table_candidate.R",
              f(f"agreement_manifest_{w}.tsv"),str(w),
              f(f"{label}.pdf"),f(f"{label}-rows.tsv"),f(f"{label}-receipt.tsv"))
    run_r(rscript,"verify_threshold_table_candidate.R",f("threshold_manifest.tsv"),
          f("main-table-2.pdf"),f("main-table-2-rows.tsv"),f("main-table-2-receipt.tsv"))
    run_r(rscript,"verify_pair_characteristics_table_candidate.R",
          f("pair_manifest.tsv"),f("supplementary-table-10-rows.tsv"),
          f("supplementary-table-10-receipt.tsv"))
    run_r(rscript,"verify_stratified_agreement_table_candidate.R",
          f("agreement_manifest_60.tsv"),f("agreement_manifest_5.tsv"),
          f("supplementary-table-13-rows.tsv"),f("supplementary-table-13-receipt.tsv"))
    run_r(rscript,"verify_agreement_diagnostic_bundle.R",
          f("diagnostic_manifest.tsv"),f("supplementary-table-14-rows.tsv"),
          f("supplementary-table-15-rows.tsv"))
    for number,kind in ((16,"reverse"),(17,"patient")):
        label=f"supplementary-table-{number}"
        run_r(rscript,"verify_threshold_aux_tables_candidate.R",
              f("threshold_manifest.tsv"),kind,f(f"{label}.pdf"),
              f(f"{label}-rows.tsv"),f(f"{label}-receipt.tsv"))
    for number,w in ((18,60),(19,5)):
        label=f"supplementary-table-{number}"
        run_r(rscript,"verify_icu_day_table_candidate.R",f(f"day_manifest_{w}.tsv"),
              str(w),f(f"{label}.pdf"),f(f"{label}-rows.tsv"),
              f(f"{label}-receipt.tsv"))
    for number,script,with_receipt in (
        (20,"verify_lmm_characteristics_table_candidate.R",False),
        (21,"verify_lmm_diagnostics_table_candidate.R",False),
        (22,"verify_lmm_balance_table_candidate.R",True),
    ):
        label=f"supplementary-table-{number}"
        arguments=[f("lmm_manifest_60.tsv"),f("lmm_manifest_5.tsv"),
                   f(f"{label}-rows.tsv")]
        if with_receipt:
            arguments.append(f(f"{label}-receipt.tsv"))
        run_r(rscript,script,*arguments)
    for number,script in ((23,"verify_rf_performance_table_candidate.R"),
                          (24,"verify_rf_importance_table_candidate.R")):
        label=f"supplementary-table-{number}"
        run_r(rscript,script,f("rf_manifest.tsv"),f(f"{label}-rows.tsv"),
              f(f"{label}-receipt.tsv"))
    run_r(rscript,"verify_rf_shap_figures_candidate.R",
          f("rf_manifest.tsv"),bundle,"bundle")
    print("THIRTY_DISPLAY_BUNDLE_INDEPENDENT_QA_PASS displays=30")


if __name__=="__main__":
    main()
