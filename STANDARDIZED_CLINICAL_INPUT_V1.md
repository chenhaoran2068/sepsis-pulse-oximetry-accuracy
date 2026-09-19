# Standardized clinical-input contract (draft v1)

Status: candidate technical interface to the oxygen-event contract. It does not extract or clean source-database records. Its formal model runners have passed a MIMIC-shaped paper-settings synthetic integration test and five cohort-shaped 30 × 20 MI/LMM tests; the five cohort-shaped random-forest tests used smoke settings. These tests do not reproduce manuscript results.

For **each cohort**, prepare `stay_clinical.tsv` once for the eligible ICU-stay roster. Prepare `pair_clinical_60.tsv` and `pair_clinical_5.tsv` separately after independent 60-minute and 5-minute one-to-one matching. All three files are UTF-8 tab-separated, with one header row and no extra columns. The package contains only newly invented demonstration rows. Authorized source data and local configurations are never distributed.

## Stay-level clinical input

`stay_clinical.tsv` has exactly one row per eligible ICU stay and these fields in order: `patient_key_internal`, `stay_key_internal`, `age_years_numeric`, `age_interval`, `age_model_representation`, `sex_class`, `cci_value`. The patient–stay mapping must match `stays.tsv` exactly. Age, sex, and CCI are not repeated with conflicting values across pairs. An empty value is an unavailable source value, not zero. Age and CCI must not be negative, and an observed CCI must be an integer. Source age intervals are retained as intervals; no arbitrary midpoint is created. `age_model_representation` is `numeric` or `source_interval` when known. An observed sex is `female` or `male`; `unknown_or_missing` is treated as unavailable by the imputation core.

## Pair-level clinical input

Each pair file has one row per pair formed **in that window**. The first fields are `patient_key_internal`, `stay_key_internal`, `spo2_stable_event_ordinal`, `sao2_stable_event_ordinal`, and a unique `pair_key_internal`. The four event keys must cover the matching output exactly. A 5-minute table is not made by filtering the 60-minute table.

Next are `sofa_total` and `sofa_status`, followed by four fields for each of `pao2`, `paco2`, `ph`, `lactate`, and `hb`: `_candidate_value`, `_canonical_unit`, `_linkage_class`, and `_value_state`. The canonical units are mmHg, mmHg, pH, mmol/L, and g/dL respectively. An observed candidate value has state `eligible_pending_candidate_build` and the specified unit and a nonblank linkage class. A held/unavailable value has no numeric value. Source-specific validity, specimen linkage, and SaO₂-time anchoring are resolved **before** this interface. The reader checks their resulting states, but cannot reconstruct source evidence from a standardized table.

Finally, each pair has `_state` and `_scope` for `concurrent_vasoactive_medication_use` and `current_invasive_mechanical_ventilation`. A permitted scope has state `yes` or `unknown`; an excluded scope has no state. Unknown must never be recoded as no treatment.

The adapter aligns these tables to paired events without dropping pairs. It calculates signed `bias_spo2_minus_sao2` and centered paired mean saturation from the observed pair. These are never imputed. Both original saturation values remain available for audit but must not enter the fourth-section mixed-effects or random-forest predictors separately. Downstream runners implement the 70% LMM availability decision, 30 × 20 imputation, mixed-effects fits, patient-balanced sensitivity, and random-forest fitting. The separate technical display driver covers 30 display identities. Supplementary Tables 1-9 and Supplementary Figure 1 remain outside this contract because they require upstream source-audit, cohort-flow, and baseline-display fields. Nothing in this contract authorizes export of completed or row-level study data.
