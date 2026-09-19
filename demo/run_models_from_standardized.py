"""One-command standardized-input analysis with isolated outputs and independent QA.

No database extraction, cleaning, baseline display, or final figure renderer is implied.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import uuid
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COHORTS = ("MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang")
REQUIRED_INPUTS = (
    "stays.tsv", "spo2_events.tsv", "sao2_events.tsv", "stay_clinical.tsv",
    "pair_clinical_60.tsv", "pair_clinical_5.tsv",
)


def run(command: list[str], env: dict[str, str]) -> None:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, env=env, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    scope = parser.add_mutually_exclusive_group(required=True)
    scope.add_argument("--synthetic", action="store_true", help="generate a new MIMIC-shaped invented fixture")
    scope.add_argument("--input", type=Path, help="six-file standardized-input directory")
    parser.add_argument("--cohort", choices=COHORTS, default="MIMIC")
    parser.add_argument("--mode", choices=("paper", "smoke"), default="paper")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--extra-r-library", type=Path, help="optional existing library; not needed if CRAN installation works")
    parser.add_argument("--run-id", help="fresh lowercase letters, digits, and dashes")
    args = parser.parse_args()
    if args.synthetic and args.cohort != "MIMIC":
        parser.error("--synthetic currently generates only the MIMIC-shaped interface")
    rscript = shutil.which(args.rscript) or (args.rscript if Path(args.rscript).is_file() else None)
    if rscript is None:
        parser.error("Rscript not found; use --rscript PATH")
    run_id = args.run_id or uuid.uuid4().hex[:12]
    if re.fullmatch(r"[a-z0-9-]{1,32}", run_id) is None:
        parser.error("run ID must contain only lowercase letters, digits, or dashes")
    parent = ROOT / "runs"
    prefix = f"models-{run_id}"
    stage = lambda name: parent / f"{prefix}-{name}"
    used = list(parent.glob(f"{prefix}-*")) if parent.exists() else []
    if used:
        parser.error(f"run ID exists; refusing overwrite: {prefix}")
    if args.input and not args.input.is_dir():
        parser.error(f"standardized input directory not found: {args.input}")
    if args.input:
        missing = [name for name in REQUIRED_INPUTS if not (args.input / name).is_file()]
        if missing:
            parser.error("standardized input missing required files: " + ", ".join(missing))
    if args.extra_r_library and not args.extra_r_library.is_dir():
        parser.error(f"extra R library not found: {args.extra_r_library}")
    parent.mkdir(exist_ok=True)
    env = os.environ.copy()
    library = stage("r_library")
    env["R_LIBS"] = os.pathsep.join(str(path) for path in
                                     (library, args.extra_r_library) if path) + os.pathsep + env.get("R_LIBS", "")
    env["R9_EXTRA_R_LIB"] = str(library)
    status = stage("status.txt")
    status.write_text("FORMAL_MODEL_CHAIN_RUNNING\n", encoding="utf-8")
    try:
        run([rscript, "demo/check_r_packages.R", str(library),
             str(stage("r_package_versions.tsv"))], env)
        if args.synthetic:
            run([rscript, "src/run_upstream_synthetic.R", str(stage("invented_pairs")), "300"], env)
            run([rscript, "tests/make_rich_synthetic_standardized_input.R",
                 str(stage("invented_pairs") / "synthetic_pairs_INTERNAL_ONLY.rds"),
                 str(stage("standardized"))], env)
            standardized = stage("standardized")
        else:
            standardized = args.input.resolve()
        analysis = stage("analysis_inputs")
        run([rscript, "src/build_analysis_inputs.R", args.cohort, str(standardized), str(analysis)], env)
        run([rscript, "tests/verify_analysis_inputs.R", str(analysis),
             str(stage("analysis_input_qa"))], env)
        for window in (60, 5):
            agreement = stage(f"agreement_{window}")
            run([rscript, "src/run_agreement_from_input.R", args.cohort,
                 str(window), str(analysis), str(agreement), args.mode], env)
            run([rscript, "tests/verify_descriptive.R", "agreement", str(analysis),
                 str(agreement), str(stage(f"agreement_{window}_qa"))], env)
            threshold = stage(f"threshold_{window}")
            run([rscript, "src/run_threshold_from_input.R", args.cohort,
                 str(window), str(analysis), str(threshold), args.mode], env)
            run([rscript, "tests/verify_descriptive.R", "threshold", str(analysis),
                 str(threshold), str(stage(f"threshold_{window}_qa"))], env)
            out = stage(f"mi_lmm_{window}")
            run([rscript, "src/run_formal_mi_lmm_from_input.R", args.cohort,
                 str(window), str(analysis), str(out)], env)
            run([rscript, "tests/verify_formal_mi_lmm.R", str(analysis), str(out),
                 str(stage(f"mi_lmm_{window}_qa"))], env)
        rf = stage("rf_60")
        run([rscript, "src/run_formal_rf_from_input.R", args.cohort, str(analysis),
             str(rf), args.mode], env)
        run([rscript, "tests/verify_formal_rf.R", str(analysis), str(rf),
             str(stage("rf_60_qa"))], env)
    except subprocess.CalledProcessError as exc:
        status.write_text(f"FORMAL_MODEL_CHAIN_FAIL exit={exc.returncode}\n", encoding="utf-8")
        raise SystemExit(f"FORMAL_MODEL_CHAIN_FAIL: {prefix}; step exit {exc.returncode}") from exc
    status.write_text("ANALYSIS_CHAIN_PASS_SYNTHETIC" if args.synthetic else
                      "ANALYSIS_CHAIN_PASS_INPUT_QA_ONLY", encoding="utf-8")
    print(status.read_text(encoding="utf-8"))
    print("OUTPUT_PREFIX", parent / prefix)


if __name__ == "__main__":
    main()
