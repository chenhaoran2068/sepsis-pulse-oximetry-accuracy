# Pulse Oximetry Accuracy in Adult Sepsis — code and synthetic demonstration

**Release:** `v1.0.1`, published at <https://github.com/chenhaoran2068/sepsis-pulse-oximetry-accuracy/releases/tag/v1.0.1>.
**Profile:** `code_with_synthetic_demo`.

This candidate is intended to accompany *Pulse Oximetry Accuracy in Adult Sepsis: A Five-Cohort Retrospective Study Across North America, Europe, and East Asia*. It supplies a documented standardized-input interface and executable analysis modules tested with entirely invented records. The interface begins **after** database-specific extraction, cohort selection, source review, unit harmonization, and event-level quality control. It contains no clinical records or manuscript estimates.

## What this candidate can and cannot do

The repository boundary is fixed as follows:

| Included | Not included |
| --- | --- |
| Specifications for six standardized TSV inputs per cohort | Source-database credentials or patient data |
| An entirely invented six-file example | Database-specific extraction, sepsis-cohort construction, and raw-data cleaning |
| Agreement, threshold, MI/LMM, patient-balanced, random-forest, and grouped-SHAP modules | Real study inputs, imputed datasets, model objects, aggregates, or manuscript estimates |
| One-command synthetic demonstrations and independent QA | A claim that invented data reproduce the paper's numerical results |
| Technical renderers for 30 display identities | Supplementary Tables 1-9, Supplementary Figure 1, or submission-ready typesetting |

There are two distinct one-command analysis entrypoints. `demo/run_demo.py` is a small technical demonstration of event validation, independent 60- and 5-minute one-to-one matching, selected descriptive calculations, reduced model interfaces, and five-cohort-shaped event checks. `demo/run_models_from_standardized.py` runs the integrated analysis modules from the six standardized input files: independent pairing and clinical alignment, agreement and threshold aggregates in both windows, 30-imputation × 20-iteration mixed-effects and patient-balanced models in both windows, and 60-minute random forest with held-out performance and grouped SHAP outputs. `--mode smoke` reduces only the descriptive bootstrap and forest computations. The imputation/LMM portion still runs 30 × 20. All synthetic estimates have no clinical meaning.

The copied formal calculation cores and portable runners passed a clean-directory MIMIC-shaped synthetic chain with independent checks. Five cohort-shaped invented inputs also passed formal 30 × 20 MI/LMM quality gates in both pairing windows and 60-minute RF smoke checks. A separate one-command technical rendering chain generates **30 of the 40 current display identities** from five completed cohort runs; its 30 PDFs passed independent synthetic-data QA. The ten excluded displays require upstream fields that are deliberately outside this standardized-input boundary. See [display coverage](expected/paper-output-map.md) and the [R2 display status](DISPLAY_INTEGRATION_STATUS.md).

## Repository layout and code roles

| Path | Role |
| --- | --- |
| `STANDARDIZED_INPUT_V1.md` | Oxygen-event and stay-roster input schema |
| `STANDARDIZED_CLINICAL_INPUT_V1.md` | Stay- and pair-level clinical-variable input schema |
| `examples/synthetic_mimic/` | Static, entirely invented six-file input example |
| `code/cores/` | Calculation cores retained from the reviewed formal implementations |
| `src/` | Portable validators, adapters, analysis runners, and technical display renderers |
| `demo/` | User-facing one-command orchestration entrypoints |
| `qa/` | Independent checks used by the analysis entrypoints and release-inventory screen |
| `tests/` | Synthetic success, failure, provenance, arithmetic, and display checks |
| `expected/` | Output contracts and the 30-of-40 display-coverage map |
| `runs/` | Generated local outputs; Git-ignored and never part of the release |

`code/cores/` and `src/` contain analysis implementation. `demo/` defines the supported execution order. Users should invoke the documented `demo/` entrypoints instead of guessing a sequence of individual scripts. The QA and test code validates interfaces and calculations; it is not a second source of manuscript results.

## Requirements and one-command synthetic run

Install R (tested with R 4.5.1) and Python 3. Network access to CRAN is needed on the first run if R packages are missing. The command checks all required R packages, installs missing packages into an **isolated library under `runs/`**, records actual package versions, then executes the stages and independent aggregate checks in order. It does not install into a global R library. Package installation may take several minutes. Other operating systems have not yet been tested.

From the repository root, use Python 3 (on Windows, the `py -3` launcher is supported):

```text
py -3 demo/run_demo.py
```

On systems without the Windows `py` launcher, use `python3 demo/run_demo.py`. If `Rscript` is not on `PATH`, add `--rscript "PATH_TO_RSCRIPT"`. The entrypoint assigns a fresh run ID by default. The optional `--run-id` must be unused and contain only lowercase letters, digits, and dashes. Each stage writes to its own `runs/demo-RUN_ID-STAGE` path. The command refuses to overwrite any prior stage or package receipt and prints `SYNTHETIC_DEMO_PASS` only after all steps and QA pass. Failures state the run ID and failing step. `runs/` is ignored by Git and should not be published. These are demonstrations, not formal model checkpoint-resume runs.

The narrower commands and individual file contracts are in [standardized event input](STANDARDIZED_INPUT_V1.md), [standardized clinical input](STANDARDIZED_CLINICAL_INPUT_V1.md), and [expected outputs](expected/expected-output-contract.md). The code uses repository-relative paths or caller-provided paths.

The separate `demo/run_formal_core_checks.py` command tests copied scientific cores on invented frames. It is not the integrated standardized-input chain. The following command generated 300 invented MIMIC-shaped ICU stays and passed the integrated clean-directory chain and independent QA on Windows with R 4.5.1. It writes only under a fresh, Git-ignored `runs/models-RUN_ID-*` prefix:

```text
py -3 demo/run_models_from_standardized.py --synthetic --mode smoke --run-id NEW_RUN_ID --rscript "PATH_TO_RSCRIPT"
```

If packages are already installed in another isolated library, `--extra-r-library PATH_TO_LIBRARY` may be supplied. For a locally authorized standardized six-file directory, replace `--synthetic` with `--input PATH_TO_STANDARDIZED_INPUT --cohort COHORT`, where `COHORT` is `MIMIC`, `Amsterdam`, `eICU`, `SICDB`, or `Lianyungang`. `--mode paper` uses the original descriptive bootstrap and forest settings; it does **not** validate a clinical result merely by finishing. The script refuses an existing run ID, reports the failing stage, and does not overwrite prior runs. The forest stage has guarded checkpoint/resume within its own runner; the integrated entrypoint and formal MI/LMM stage do not yet resume an interrupted run. Use a new run ID after a failure. The 20-stay static example in `examples/synthetic_mimic/` demonstrates the file contract but is too small and homogeneous for the 30 × 20 model fit; use `--synthetic` for that fit.

For a locally authorized standardized input directory containing all six TSV files described in the two input contracts, the isolated intermediate-input command is:

```text
Rscript src/build_analysis_inputs.R MIMIC PATH_TO_STANDARDIZED_INPUT NEW_OUTPUT_DIRECTORY
```

Replace `MIMIC` with the applicable cohort key. The command refuses to overwrite its output. Its RDS output may contain patient-level data when run by an authorized researcher and must never be committed to the repository.

The six files in `examples/synthetic_mimic/` are a static, entirely invented standardized-input example. They were generated by the repository's synthetic event and clinical-input scripts, not extracted or sampled from a study database. A smoke test can use `examples/synthetic_mimic` as `PATH_TO_STANDARDIZED_INPUT` and a fresh directory under ignored `runs/` as `NEW_OUTPUT_DIRECTORY`.

## One-command technical rendering from completed cohort runs

After running the analysis entrypoint for **each of the five cohorts** with a fresh run ID, prepare a five-row UTF-8 TSV run map with these columns in this exact order: `cohort`, `analysis_dir`, `agreement_60_dir`, `agreement_5_dir`, `threshold_60_dir`, `threshold_5_dir`, `lmm_60_dir`, `lmm_5_dir`, `rf_60_dir`. Use the stage directories printed by each completed analysis run; no clinical data or run map should be committed. The display entrypoint verifies each source file and the SHA-256 identity of each stage against its cohort's analysis input **before** creating output. It refuses an existing run ID.

```text
py -3 demo/render_30_displays.py --run-map PATH_TO_FIVE_COHORT_RUN_MAP.tsv --run-id NEW_RUN_ID --rscript "PATH_TO_RSCRIPT"
```

This creates 30 technical PDFs, their receipts, an output inventory, and a status file under Git-ignored `runs/displays-NEW_RUN_ID/`. Its scope is exactly the 30 identities marked in [display coverage](expected/paper-output-map.md); it does not render Supplementary Tables 1–9 or Supplementary Figure 1. The separate second-pass command checks source provenance, output hashes, and independent values for the implemented tables and SHAP figures:

```text
py -3 tests/verify_30_display_bundle.py --bundle runs/displays-NEW_RUN_ID --run-map PATH_TO_FIVE_COHORT_RUN_MAP.tsv --rscript "PATH_TO_RSCRIPT"
```

A wrong cohort/window source, missing output, or reused run ID fails clearly. The rendering test used entirely invented inputs and RF smoke settings; it did not regenerate the approved clinical estimates or submission-ready figure styling. The historical [formal-code integration status](FORMAL_CODE_INTEGRATION_STATUS.md) predates this display extension.

## Data and reuse

No real study records or real-result row-level files are distributed. See [data access](DATA_ACCESS.md) for the five source datasets and the distinction between generated demonstration data and clinical results. The original code in this candidate is proposed under the [MIT license](LICENSE). No third-party data, article text, or manuscript figure is licensed by that code license.

## Citation and maintenance

The [citation file](CITATION.cff) records the author order, paper title, repository, and fixed release version. Cite the tagged release rather than a moving branch. Corrections after release should use a new commit and tag rather than changing an already cited version. This project is a research-code supplement, not a clinical decision-support product.
