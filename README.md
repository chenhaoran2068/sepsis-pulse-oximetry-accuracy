# Pulse Oximetry Accuracy in Adult Sepsis — code and synthetic demonstration

**Status:** release candidate under review; no tagged release exists yet.  
**Profile:** `code_with_synthetic_demo`.

This repository is intended to accompany *Pulse Oximetry Accuracy in Adult Sepsis: A Five-Cohort Retrospective Study Across North America, Europe, and East Asia*. It supplies selected analysis functions, a documented standardized-input interface, and executable tests using newly generated, entirely invented records. The repository begins **after** database-specific extraction, cohort selection, and source-data cleaning. It contains no clinical records.

## What this candidate can and cannot do

The single demonstration command exercises source-event validation, one-to-one 60- and 5-minute matching, agreement and threshold calculations, a small 3-imputation × 5-iteration mixed-model pilot, a reduced random-forest/SHAP interface, five-cohort-shaped synthetic event checks, and standardized clinical-variable alignment to model input schemas. The imputation and forest demonstrations are small technical tests whose estimates have no clinical meaning.

This candidate **cannot reproduce the paper's numerical results** from its scripts and invented records alone. It does not include database-specific extraction or cleaning, the formal 30-imputation × 20-iteration run, the final five-cohort mixed-effects and random-forest fits, or final rendering code for all 40 manuscript tables and figures. The [output-coverage map](expected/paper-output-map.md) records this boundary item by item. Do not cite a synthetic estimate as a study result.

## Requirements and one-command synthetic run

Install R (tested with R 4.5.1) and Python 3. Network access to CRAN is needed on the first run if R packages are missing. The command checks all required R packages, installs missing packages into an **isolated library under `runs/`**, records actual package versions, then executes the stages and independent aggregate checks in order. It does not install into a global R library. Package installation may take several minutes. Other operating systems have not yet been tested.

From the repository root, use Python 3 (on Windows, the `py -3` launcher is supported):

```text
py -3 demo/run_demo.py
```

On systems without the Windows `py` launcher, use `python3 demo/run_demo.py`. If `Rscript` is not on `PATH`, add `--rscript "PATH_TO_RSCRIPT"`. The entrypoint assigns a fresh run ID by default. The optional `--run-id` must be unused and contain only lowercase letters, digits, and dashes. Each stage writes to its own `runs/demo-RUN_ID-STAGE` path. The command refuses to overwrite any prior stage or package receipt and prints `SYNTHETIC_DEMO_PASS` only after all steps and QA pass. Failures state the run ID and failing step. `runs/` is ignored by Git and should not be published. These are demonstrations, not formal model checkpoint-resume runs.

The narrower commands and individual file contracts are in [standardized event input](STANDARDIZED_INPUT_V1.md), [standardized clinical input](STANDARDIZED_CLINICAL_INPUT_V1.md), and [expected outputs](expected/expected-output-contract.md). The code uses repository-relative paths or caller-provided paths.

## Data and reuse

No real study records or real-result row-level files are distributed. See [data access](DATA_ACCESS.md) for the five source datasets and the distinction between generated demonstration data and clinical results. The original code in this candidate is proposed under the [MIT license](LICENSE). No third-party data, article text, or manuscript figure is licensed by that code license.

## Citation and maintenance

The [citation file](CITATION.cff) records the author order and paper title but deliberately omits a release version and publication identifier until a reviewed tag exists. Use the exact tagged release, not this candidate name, when citing code. The planned repository owner is Haoran Chen. Corrections after release should use a new commit and tag rather than changing an already cited version. This project is a research-code supplement, not a clinical decision-support product.
