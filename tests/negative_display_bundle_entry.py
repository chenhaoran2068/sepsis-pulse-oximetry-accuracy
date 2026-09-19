"""Reject a wrong-cohort/window run map and refuse reuse of an existing run ID."""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--run-map",type=Path,required=True)
    p.add_argument("--existing-bundle",type=Path,required=True)
    p.add_argument("--rscript",default="Rscript")
    args=p.parse_args()
    rscript=shutil.which(args.rscript) or (
        args.rscript if Path(args.rscript).is_file() else None)
    if not rscript:
        p.error("Rscript unavailable")
    with args.run_map.open(encoding="utf-8",newline="") as handle:
        reader=csv.DictReader(handle,delimiter="\t")
        fields=reader.fieldnames
        rows=list(reader)
    if not fields or len(rows)!=5:
        raise RuntimeError("TEST_RUN_MAP_INVALID")
    inventory=args.existing_bundle/"display_inventory.tsv"
    before=hashlib.sha256(inventory.read_bytes()).hexdigest()
    trial=ROOT/"runs"/f"negative-display-{uuid.uuid4().hex[:12]}"
    trial.mkdir(parents=True)
    for row in rows:
        for field in fields[1:]:
            source=Path(row[field])
            row[field]=str((source if source.is_absolute() else
                            args.run_map.parent/source).resolve())
        if row["cohort"]=="SICDB":
            row["agreement_60_dir"]=row["agreement_5_dir"]
    wrong=trial/"wrong-window-map.tsv"
    with wrong.open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",
                              lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    result=subprocess.run([str(rscript),str(ROOT/"tests/verify_display_run_map.R"),
                           str(wrong)],cwd=ROOT,capture_output=True,text=True)
    if result.returncode==0 or "DISPLAY_RUN_MAP_PROVENANCE_MISMATCH" not in (
        result.stdout+result.stderr):
        raise RuntimeError("WRONG_WINDOW_NOT_REJECTED")
    reused=args.existing_bundle.name.removeprefix("displays-")
    result=subprocess.run([sys.executable,str(ROOT/"demo/render_30_displays.py"),
                           "--run-map",str(args.run_map.resolve()),
                           "--run-id",reused,"--rscript",str(rscript)],
                          cwd=ROOT,capture_output=True,text=True)
    if result.returncode==0 or "run ID exists; refusing overwrite" not in (
        result.stdout+result.stderr):
        raise RuntimeError("RUN_ID_REUSE_NOT_REJECTED")
    if hashlib.sha256(inventory.read_bytes()).hexdigest()!=before:
        raise RuntimeError("PRIOR_BUNDLE_MUTATED")
    print("DISPLAY_BUNDLE_NEGATIVE_PASS wrong_window_rejected=1 overwrite_rejected=1")


if __name__=="__main__":
    main()
