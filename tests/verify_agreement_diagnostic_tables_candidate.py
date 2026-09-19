"""Independent invented-data and negative checks for Tables S14/S15 renderers."""
import csv
import pathlib
import subprocess
import sys


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    if len(sys.argv) != 4:
        raise SystemExit("Usage: Rscript.exe renderer.R NEW_TEST_DIR")
    rscript, renderer, root = sys.argv[1], pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
    if root.exists():
        raise SystemExit("TEST_OUTPUT_EXISTS_REFUSE_OVERWRITE")
    root.mkdir(parents=True)
    cohorts = ["MIMIC", "Amsterdam", "eICU", "SICDB", "Lianyungang"]
    manifest = []
    expected = {}
    for i, cohort in enumerate(cohorts):
        for window in (60, 5):
            folder = root / "input" / f"{cohort}_{window}"
            folder.mkdir(parents=True)
            manifest.append({"cohort": cohort, "window_minutes": window,
                             "diagnostic_dir": str(folder)})
            index = i + (0 if window == 60 else 5)
            slope = (index - 4) / 100
            delta = -(index + 1) * 0.013
            intercept = 0.6 + index / 10
            write_tsv(folder / "agreement_proportional_bias.tsv",
                      ["cohort", "window_minutes", "slope", "ci_lower", "ci_upper",
                       "model_estimated_bias_at_90", "singular_fit"],
                      [{"cohort": cohort, "window_minutes": window, "slope": slope,
                        "ci_lower": slope - 0.02, "ci_upper": slope + 0.02,
                        "model_estimated_bias_at_90": intercept, "singular_fit": "FALSE"}])
            write_tsv(folder / "agreement_heteroscedasticity.tsv",
                      ["cohort", "window_minutes", "variance_function_parameter",
                       "ci_lower", "ci_upper", "residual_sd_at_90", "likelihood_ratio", "p_value"],
                      [{"cohort": cohort, "window_minutes": window,
                        "variance_function_parameter": delta,
                        "ci_lower": delta - 0.004, "ci_upper": delta + 0.004,
                        "residual_sd_at_90": 2 + index / 10,
                        "likelihood_ratio": 5 + index, "p_value": 0.0005 if i < 3 else 0.02}])
            expected[(cohort, window)] = (slope, intercept, delta, 2 + index / 10)
    manifest_path = root / "manifest.tsv"
    write_tsv(manifest_path, ["cohort", "window_minutes", "diagnostic_dir"], manifest)
    checks = 0
    for kind in ("proportional", "heteroscedastic"):
        out = root / kind
        out.mkdir()
        cmd = [rscript, str(renderer), str(manifest_path), kind,
               str(out / "table.pdf"), str(out / "rows.tsv"), str(out / "receipt.tsv")]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode:
            raise AssertionError(result.stdout + result.stderr)
        if not (out / "table.pdf").read_bytes().startswith(b"%PDF"):
            raise AssertionError("PDF_MISSING_OR_INVALID")
        rows = read_tsv(out / "rows.tsv")
        if len(rows) != 10:
            raise AssertionError("TEN_ROWS_EXPECTED")
        checks += 2
        pretty = {"MIMIC": "MIMIC", "Amsterdam": "AmsterdamUMCdb", "eICU": "eICU",
                  "SICDB": "SICdb", "Lianyungang": "Lianyungang"}
        for j, (window, cohort) in enumerate((w, c) for w in (60, 5) for c in cohorts):
            row = rows[j]
            if row["cohort"] != pretty[cohort] or row["window"] != f"{window} min":
                raise AssertionError("ROW_ORDER_OR_IDENTITY_INVALID")
            slope, intercept, delta, sigma = expected[(cohort, window)]
            if kind == "proportional":
                if row["slope_with_95ci"] != f"{slope:.2f} ({slope-0.02:.2f} to {slope+0.02:.2f})":
                    raise AssertionError("SLOPE_OR_CI_WRONG")
                if row["model_estimated_bias_at_90"] != f"{intercept:.2f}":
                    raise AssertionError("INTERCEPT_WRONG")
            else:
                if row["variance_parameter_with_95ci"] != f"{delta:.3f} ({delta-0.004:.3f} to {delta+0.004:.3f})":
                    raise AssertionError("VARIANCE_PARAMETER_WRONG")
                if row["residual_sd_at_90"] != f"{sigma:.2f}":
                    raise AssertionError("RESIDUAL_SD_WRONG")
            checks += 4
        duplicate = subprocess.run(cmd, capture_output=True, text=True)
        if duplicate.returncode == 0 or "OUTPUT_EXISTS_REFUSE_OVERWRITE" not in duplicate.stderr:
            raise AssertionError("OVERWRITE_GUARD_FAILED")
        checks += 1
    print(f"AGREEMENT_DIAGNOSTIC_TABLES_INDEPENDENT_QA_PASS checks={checks}")


if __name__ == "__main__":
    main()
