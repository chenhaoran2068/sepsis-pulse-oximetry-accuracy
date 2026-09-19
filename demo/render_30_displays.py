"""Render the 30 supported technical display identities from five completed cohort runs.

The ten source/cohort-flow/baseline displays are intentionally not claimed here.
No patient-level input or generated output belongs in a public commit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
import subprocess
import uuid
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COHORTS = ("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
COLUMNS = (
    "cohort", "analysis_dir", "agreement_60_dir", "agreement_5_dir",
    "threshold_60_dir", "threshold_5_dir", "lmm_60_dir", "lmm_5_dir", "rf_60_dir",
)


def read_map(path: Path) -> dict[str, dict[str, Path]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != COLUMNS:
            raise ValueError("RUN_MAP_COLUMNS_INVALID: expected " + ", ".join(COLUMNS))
        rows = list(reader)
    if len(rows) != 5 or {row["cohort"] for row in rows} != set(COHORTS):
        raise ValueError("RUN_MAP_COHORT_SET_INVALID")
    result: dict[str, dict[str, Path]] = {}
    for row in rows:
        co = row["cohort"]
        if co in result:
            raise ValueError("RUN_MAP_DUPLICATE_COHORT")
        result[co] = {}
        for field in COLUMNS[1:]:
            raw = row[field]
            if not raw:
                raise ValueError(f"RUN_MAP_EMPTY_PATH {co} {field}")
            p = Path(raw)
            if not p.is_absolute():
                p = path.parent / p
            p = p.resolve()
            if not p.is_dir():
                raise ValueError(f"RUN_MAP_DIRECTORY_MISSING {co} {field}: {p}")
            result[co][field] = p
    required = {
        "analysis_dir": ("matched_pairs_60.rds", "matched_pairs_5.rds",
                         "analytic_pairs_60.rds", "analytic_pairs_5.rds"),
        "agreement_60_dir": ("agreement_overall.tsv", "agreement_by_spo2_stratum.tsv",
                             "agreement_proportional_bias.tsv", "agreement_heteroscedasticity.tsv",
                             "agreement_patient_balanced.tsv"),
        "agreement_5_dir": ("agreement_overall.tsv", "agreement_by_spo2_stratum.tsv"),
        "threshold_60_dir": ("threshold_from_low_sao2.tsv", "threshold_day1_7_strata.tsv"),
        "threshold_5_dir": ("threshold_from_low_sao2.tsv", "threshold_day1_7_strata.tsv"),
        "lmm_60_dir": ("aggregate/lmm_primary_pooled.tsv",
                       "aggregate/lmm_patient_balanced_pooled.tsv"),
        "lmm_5_dir": ("aggregate/lmm_primary_pooled.tsv",
                      "aggregate/lmm_patient_balanced_pooled.tsv"),
        "rf_60_dir": ("rf_shap_summary.tsv", "rf_shap_sample.tsv",
                      "rf_shap_feature_values.tsv", "rf_shap_grouped.rds"),
    }
    for co in COHORTS:
        for field, names in required.items():
            for name in names:
                p = result[co][field] / name
                if not p.is_file():
                    raise ValueError(f"RUN_MAP_SOURCE_MISSING {co} {field}: {p}")
    return result


def write_manifest(path: Path, header: tuple[str, ...], rows: list[tuple[object, ...]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-map", type=Path, required=True,
                        help="five-row TSV of already completed, matching cohort runs")
    parser.add_argument("--run-id", help="new lowercase letters/digits/dashes ID")
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    rscript = shutil.which(args.rscript) or (
        args.rscript if Path(args.rscript).is_file() else None
    )
    if not rscript:
        parser.error("Rscript not found; use --rscript PATH")
    if not args.run_map.is_file():
        parser.error("run map not found")
    try:
        sources = read_map(args.run_map.resolve())
    except ValueError as exc:
        parser.error(str(exc))
    provenance = subprocess.run(
        [str(rscript), str(ROOT / "tests/verify_display_run_map.R"),
         str(args.run_map.resolve())],
        cwd=ROOT, check=False,
    )
    if provenance.returncode:
        parser.error("cohort-run provenance check failed; no display output was created")
    run_id = args.run_id or uuid.uuid4().hex[:12]
    if not re.fullmatch(r"[a-z0-9-]{1,32}", run_id):
        parser.error("run ID must contain only lowercase letters, digits, or dashes")
    output = ROOT / "runs" / f"displays-{run_id}"
    if output.exists():
        parser.error(f"run ID exists; refusing overwrite: {output}")
    output.mkdir(parents=True)
    status = output / "status.txt"
    status.write_text("TECHNICAL_DISPLAY_CHAIN_RUNNING\n", encoding="utf-8")
    try:
        manifests = {}
        for w in (60, 5):
            ba = output / f"agreement_manifest_{w}.tsv"
            day = output / f"day_manifest_{w}.tsv"
            lmm = output / f"lmm_manifest_{w}.tsv"
            write_manifest(ba, ("cohort", "analysis_rds", "agreement_tsv"), [
                (co, sources[co]["analysis_dir"] / f"analytic_pairs_{w}.rds",
                 sources[co][f"agreement_{w}_dir"] / "agreement_overall.tsv")
                for co in COHORTS])
            write_manifest(day, ("cohort", "day_tsv"), [
                (co, sources[co][f"threshold_{w}_dir"] / "threshold_day1_7_strata.tsv")
                for co in COHORTS])
            write_manifest(lmm, ("cohort", "lmm_tsv"), [
                (co, sources[co][f"lmm_{w}_dir"] / "aggregate/lmm_primary_pooled.tsv")
                for co in COHORTS])
            manifests[f"ba{w}"], manifests[f"day{w}"], manifests[f"lmm{w}"] = ba, day, lmm
        threshold = output / "threshold_manifest.tsv"
        write_manifest(threshold, ("cohort", "window_minutes", "low_sao2_tsv"), [
            (co, w, sources[co][f"threshold_{w}_dir"] / "threshold_from_low_sao2.tsv")
            for co in COHORTS for w in (60, 5)])
        diagnostic = output / "diagnostic_manifest.tsv"
        write_manifest(diagnostic, ("cohort", "window_minutes", "diagnostic_dir"), [
            (co, w, sources[co][f"agreement_{w}_dir"])
            for co in COHORTS for w in (60, 5)])
        rf = output / "rf_manifest.tsv"
        write_manifest(rf, ("cohort", "rf_dir"), [
            (co, sources[co]["rf_60_dir"]) for co in COHORTS])
        pair = output / "pair_manifest.tsv"
        write_manifest(pair, ("cohort", "analysis_dir"), [
            (co, sources[co]["analysis_dir"]) for co in COHORTS])
        jobs: list[tuple[str, str, list[str], bool]] = []
        def add(label: str, script: str, inputs: list[object], rows: bool = False) -> None:
            jobs.append((label, script, [str(x) for x in inputs], rows))
        add("main-table-1","render_agreement_table_candidate.R",[manifests["ba60"],60],True)
        add("main-table-2","render_threshold_table_candidate.R",[threshold],True)
        add("main-figure-1","render_bland_altman_candidate.R",[manifests["ba60"],60])
        add("main-figure-2","render_icu_day_candidate.R",[manifests["day60"],60])
        add("main-figure-3","render_lmm_forest_candidate.R",[manifests["lmm60"],60],True)
        add("supplementary-table-10","render_pair_characteristics_table_candidate.R",[pair],True)
        add("supplementary-table-11","render_agreement_table_candidate.R",[manifests["ba5"],5],True)
        add("supplementary-table-12","render_patient_balance_agreement_table_candidate.R",[manifests["ba60"]],True)
        add("supplementary-table-13","render_stratified_agreement_table_candidate.R",[manifests["ba60"],manifests["ba5"]],True)
        add("supplementary-table-14","render_agreement_diagnostic_tables_candidate.R",[diagnostic,"proportional"],True)
        add("supplementary-table-15","render_agreement_diagnostic_tables_candidate.R",[diagnostic,"heteroscedastic"],True)
        add("supplementary-table-16","render_threshold_aux_tables_candidate.R",[threshold,"reverse"],True)
        add("supplementary-table-17","render_threshold_aux_tables_candidate.R",[threshold,"patient"],True)
        add("supplementary-table-18","render_icu_day_table_candidate.R",[manifests["day60"],60],True)
        add("supplementary-table-19","render_icu_day_table_candidate.R",[manifests["day5"],5],True)
        add("supplementary-table-20","render_lmm_characteristics_table_candidate.R",[manifests["lmm60"],manifests["lmm5"]],True)
        add("supplementary-table-21","render_lmm_diagnostics_table_candidate.R",[manifests["lmm60"],manifests["lmm5"]],True)
        add("supplementary-table-22","render_lmm_balance_table_candidate.R",[manifests["lmm60"],manifests["lmm5"]],True)
        add("supplementary-table-23","render_rf_performance_table_candidate.R",[rf],True)
        add("supplementary-table-24","render_rf_importance_table_candidate.R",[rf],True)
        add("supplementary-figure-2","render_bland_altman_candidate.R",[manifests["ba5"],5])
        add("supplementary-figure-3","render_stratified_bland_altman_candidate.R",[manifests["ba60"]])
        add("supplementary-figure-4","render_spo2_stratum_candidate.R",[manifests["ba60"]])
        add("supplementary-figure-5","render_icu_day_candidate.R",[manifests["day5"],5])
        add("supplementary-figure-6","render_lmm_forest_candidate.R",[manifests["lmm5"],5],True)
        for number, co in enumerate(COHORTS, start=7):
            add(f"supplementary-figure-{number}","render_rf_shap_figure_candidate.R",
                [sources[co]["rf_60_dir"],co,number])
        if len(jobs)!=30:
            raise RuntimeError("DISPLAY_JOB_COUNT_INVALID")
        inventory = []
        for label, script, inputs, has_rows in jobs:
            pdf = output / f"{label}.pdf"
            receipt = output / f"{label}-receipt.tsv"
            command = [str(rscript),str(ROOT/"src"/script),*inputs,str(pdf)]
            if has_rows:
                command.append(str(output/f"{label}-rows.tsv"))
            command.append(str(receipt))
            print("RUN",label,flush=True)
            subprocess.run(command,cwd=ROOT,check=True)
            if pdf.stat().st_size < 1000 or not receipt.is_file():
                raise RuntimeError(f"DISPLAY_OUTPUT_INCOMPLETE {label}")
            inventory.append((label,str(pdf.relative_to(ROOT)),hashlib.sha256(pdf.read_bytes()).hexdigest(),
                              pdf.stat().st_size))
        write_manifest(output/"display_inventory.tsv",
                       ("display","pdf_path","sha256","bytes"),inventory)
        status.write_text("TECHNICAL_DISPLAY_CHAIN_PASS_30_OF_40\n",encoding="utf-8")
        print(status.read_text(encoding="utf-8"),end="")
        print("OUTPUT_DIRECTORY",output)
    except (subprocess.CalledProcessError,RuntimeError) as exc:
        status.write_text(f"TECHNICAL_DISPLAY_CHAIN_FAIL: {exc}\n",encoding="utf-8")
        raise SystemExit(status.read_text(encoding="utf-8")) from exc


if __name__ == "__main__":
    main()
