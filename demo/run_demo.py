"""Run the synthetic demonstration from one entrypoint with isolated run IDs."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import uuid
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str], env: dict[str, str]) -> None:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, env=env, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable")
    parser.add_argument("--run-id", help="optional fresh run ID, lowercase letters, digits, and dashes")
    args = parser.parse_args()
    rscript = shutil.which(args.rscript) or (args.rscript if Path(args.rscript).is_file() else None)
    if rscript is None:
        raise SystemExit("DEMO_RSCRIPT_NOT_FOUND: install R or pass --rscript PATH")
    run_id = args.run_id or uuid.uuid4().hex[:12]
    if re.fullmatch(r"[a-z0-9-]{1,32}", run_id) is None:
        raise SystemExit("DEMO_RUN_ID_INVALID: use lowercase letters, digits, or dashes")
    run_parent = ROOT / "runs"
    run_parent.mkdir(exist_ok=True)
    paths = {name: run_parent / f"demo-{run_id}-{name}" for name in
             ("upstream", "mi", "rf", "cohorts", "event", "clinical", "model")}
    run_library = run_parent / f"demo-{run_id}-r_library"
    receipt = run_parent / f"demo-{run_id}-r_package_versions.tsv"
    if any(path.exists() for path in (*paths.values(), run_library, receipt)):
        raise SystemExit(f"DEMO_RUN_ID_EXISTS: {run_id}; choose a fresh run ID")
    r = {key: str(value) for key, value in paths.items()}
    env = os.environ.copy()
    env["R_LIBS"] = str(run_library) + os.pathsep + env.get("R_LIBS", "")
    try:
        run([rscript, "demo/check_r_packages.R", str(run_library), str(receipt)], env)
        run([rscript, "src/run_upstream_synthetic.R", r["upstream"], "200"], env)
        run([sys.executable, "qa/verify_upstream.py", r["upstream"]], env)
        run([rscript, "src/run_mi_lmm_synthetic.R", r["upstream"], r["mi"]], env)
        run([sys.executable, "qa/verify_models.py", "mi_lmm", r["mi"]], env)
        run([rscript, "src/run_rf_synthetic.R", r["upstream"], r["rf"]], env)
        run([sys.executable, "qa/verify_models.py", "rf", r["rf"]], env)
        run([rscript, "src/run_four_cohort_interfaces_synthetic.R", r["cohorts"]], env)
        run([sys.executable, "qa/verify_four_cohort_interfaces.py", r["cohorts"]], env)
        run([rscript, "src/run_standardized_event_demo.R", r["event"]], env)
        run([sys.executable, "qa/verify_standardized_event_demo.py", r["event"]], env)
        run([rscript, "src/run_clinical_input_demo.R", r["event"], r["clinical"]], env)
        run([sys.executable, "qa/verify_clinical_input_demo.py", r["clinical"]], env)
        run([rscript, "src/run_model_input_demo.R", r["event"], r["clinical"], r["model"]], env)
        run([sys.executable, "qa/verify_model_input_demo.py", r["model"]], env)
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"DEMO_FAILED step_exit={exc.returncode} run_id={run_id} runs={run_parent}") from exc
    print("SYNTHETIC_DEMO_PASS", run_id)
    print("OUTPUT_PARENT", run_parent)


if __name__ == "__main__":
    main()
