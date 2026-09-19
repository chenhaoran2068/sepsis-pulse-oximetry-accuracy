# Formal-code integration status — private candidate, 2026-09-19

This preserved R1 integration record predates the R2 technical display work. Its statements that all display renderers are absent describe that earlier stage. See `DISPLAY_INTEGRATION_STATUS.md` and `expected/paper-output-map.md` for the current private candidate boundary.

This is a technical candidate, not a released paper-result reproduction. Every reported test used entirely invented records in a new `runs/` path. No accepted clinical result, model, figure, manifest, or QA record was overwritten. The current manuscript sources were checked against `09_manuscript/CURRENT.md` before this integration.

## Controlled scientific cores and one corrected source-version mismatch

The copied scientific cores below match the named controlled source versions by SHA-256. The portable standardized-input runners are new wrappers and are **not** claimed to be byte-identical to the controlled real-data wrappers.

| Candidate core | SHA-256 |
| --- | --- |
| `stage405_multilevel_core.R` | `2A08798D3327E2718A63A3CDCC5A7D71EF43975674B73E32ACD4A8E0C9DCA93F` |
| `stage405_formal_convergence_adapter.R` | `01ECCA9A244EAEF86D3BD2C1137DB914B9135C7A20ABBC7E1DAC876C8A3F4A7B` |
| `stage405_formal_imputation.R` | `0C2EB110D15E5F7DF07FF920C36BB612850399DF4934599F4D71030C4614AB20` |
| `stage405_formal_lmm_stream.R` | `6EF950E7FFCAF3765248D138AE23928F43B70790D8716F955E4AD205F17B1EE0` |
| `stage405_rf_amended_formal_core.R` | `94AAF8F9A20E70DA2D6BAEEECCBA0BA62761035B112FC9B106EF247ED47DC09A` |
| `stage405_rf_amended_fitting_core.R` | `53F5AB508FC806DF6D8F0770255FEA91DAE05839B7ADD1F8CF81D69183A23163` |
| `r211e_core.R` | `32E2C4DC0B71E6B2267017A7D0C3F0A0DED4DE9492DFFD79B64896F46A334013` |
| `stage4_variable_validity_core.R` | `B8FF66DB1B90D4B954370C135C943772EE0E39B7C33AE8BB13BEB288E144C777` |

The initial private copy of `stage4_variable_validity_core.R` came from the superseded pre-R1 source (`ECB58AE...`). The current formal multilevel core and approved baseline source reference the R1 validity core (`B8FF66...`). The candidate was corrected to the R1 source before the fresh post-correction clean run. The only line-level difference is removal of an obsolete severe-hypothermia review lower threshold. This did not change the invented fixture values but would have made the public code scientifically out of version if left unfixed. The controlled originals remain untouched.

## Passed technical checks

- Copied-core synthetic checks: 14 for formal imputation, 25 for LMM, and 9 for RF.
- A single clean-directory command, `demo/run_models_from_standardized.py --synthetic --mode smoke`, passed from six newly generated invented standardized TSVs through both independent pairing windows, agreement and threshold aggregates, 30 × 20 MI/LMM in each window, and the 60-minute RF. At that point independent QA passed 36 input checks, 12 agreement checks and 19 threshold checks per window, 27 MI/LMM checks per window, and 58 RF smoke checks. The complete chain passed again under the fresh `sep19-postvalidity-r1` run ID after correcting the copied validity core. The agreement runner was subsequently extended to write patient-balanced aggregates and its independent QA increased to 16 checks; this later version passed separate fresh MIMIC-shaped runs in both windows.
- The separate `--input` route also completed from a clean directory using the same six already standardized invented files and a new `sep19-external-input-r1` run ID. The route reported `ANALYSIS_CHAIN_PASS_INPUT_QA_ONLY`. That wording deliberately does not assert that clinical results match the manuscript.
- After the patient-balanced agreement output and input preflight were added, a final new-directory **paper-settings** command (`sep19-final-paper-r1`) passed the entire synthetic chain. It used 2,000 descriptive bootstrap replicates, 30 × 20 MI in each window, the full 60-minute RF tuning/final-fit/bootstrap/SHAP settings, and independent QA: input 36, agreement 16 per window, threshold 19 per window, MI/LMM 27 per window, and RF 59 checks. Its success status is `ANALYSIS_CHAIN_PASS_SYNTHETIC`, not clinical-result reproduction.
- The checkpointed RF runner also passed full paper computational settings on invented MIMIC-shaped input, followed by 59 independent checks. An intentional interruption after five inner-fold tasks was resumable with an identity-matched checkpoint. A completed run was refused for resume.
- All five cohort-shaped invented inputs passed the **current** builder with 36 independent checks each. Agreement and threshold calculations passed 16 and 19 independent checks, respectively, in **both** pairing windows for each cohort-shaped input. The checkpointed 60-minute RF smoke run passed 58 checks for each of the five cohort representations. The eICU fixture used a source-like two-stay-per-patient structure. These checks test structural interfaces, not observed clinical effects.
- All five cohort-shaped fixtures completed the 30 × 20 formal MI/LMM runner **in each independently constructed 60- and 5-minute window**, passed the final-fit quality gate, and passed 27 independent checks per cohort-window. The earlier eICU, SICdb, and Lianyungang synthetic-gate failures were retained rather than overwritten.
- Reusing an integrated run ID was refused before any output was replaced. The default static 20-stay invented example remains an input-contract example, not a suitable 30 × 20 model-fitting fixture.
- An input directory missing the six required TSVs was rejected before package installation or output creation, with the missing filenames stated explicitly.

## Failure evidence and unresolved boundaries

- The first portable LMM independent QA omitted the final-fit convergence, Hessian, singularity, variance-floor, design, VIF, and Monte Carlo error gates. It could report PASS despite two nonconverged eICU-shaped invented fits. This was a QA defect in the private wrapper, **not evidence that the accepted clinical-data LMM failed**. The wrapper now stops on the original formal LMM quality gate; the verifier independently checks the written diagnostics. The previously passed 21-check QA is superseded by 27-check QA.
- A revised eICU repeated-stay invented fixture exposed a synthetic-data error: all patients were male after naively assigning two stays per patient. The variant generator was corrected without altering the scientific core. All 30 resulting eICU LMM fits then converged, passed the Hessian/singularity/variance/design checks, but one PaCO2 pooled Monte Carlo error/SE ratio was 0.1044 against the existing 0.10 gate. This invented-data run was correctly rejected. A new fixture retained missing PaCO2 values at a lower rate and passed the same unmodified gate and independent QA. Both test outcomes are preserved.
- The integrated entrypoint has no restart of an interrupted 30 × 20 MI/LMM. Its RF sub-runner has guarded checkpoint/resume, but the top-level command refuses an existing run ID and requires a new isolated run after a failure. The README states this plainly.
- Five-cohort, both-window formal-model **synthetic** QA is closed for the current interface. This does **not** establish numeric equivalence to the accepted clinical results; full real-data reruns were deliberately out of scope. Nor does it validate the final manuscript display renderers, which are absent.
- The code produces analysis aggregates, not the final 2 main tables, 3 main figures, 24 supplementary tables, and 11 supplementary figures. None of those 40 display files has a current renderer in this candidate. Database-specific extraction, cleaning, and source-cohort construction are also outside the standardized-input boundary. See `expected/paper-output-map.md`.

## Publication boundary

The existing public GitHub repository is a limited earlier demonstration, not this private candidate. No reviewed tag or release exists for this candidate. Before any push or release, provide the owner with an exact file-and-hash inventory, disclosure scan, test receipts, explicit unresolved-limit list, and proposed public scope. The owner must confirm that precise list. Only after a fixed public release is verified may the main manuscript gain a bilingual Code availability statement and a newly built LaTeX/PDF candidate. Do not claim that a public synthetic run reproduces the paper's numerical estimates.
