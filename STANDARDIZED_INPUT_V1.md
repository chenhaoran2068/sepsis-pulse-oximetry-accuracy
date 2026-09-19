# Standardized input contract: oxygen-event entrypoint (draft v1)

Status: candidate technical contract for analyses beginning at standardized inputs. It begins **after** each source database has been accessed, the eligible sepsis ICU stays identified, oxygen records selected, units and times harmonized, and event-level source and artifact rules applied. No raw database extraction or source-specific cleaning code is in this module.

The three UTF-8 tab-separated files below belong to **one cohort at a time**. They are supplied locally by a researcher with authorized data access. The package contains only files generated from invented records. Never commit real input files or local run outputs.

## `stays.tsv`

Exactly two columns, in order: `stay_key_internal`, `patient_key_internal`. One row per eligible ICU stay. Stay keys must be unique and nonblank. A patient may contribute more than one stay where the cohort rules allow it. Both keys are local analysis identifiers, not names or public identifiers.

## `spo2_events.tsv` and `sao2_events.tsv`

Exactly these columns, in order: `stay_key_internal`, `patient_key_internal`, `event_time_min`, `saturation_percent`, `source_family`, `source_priority`, `range_status`, `analysis_eligible_range`, `artifact_qc_status`, `possible_transient_artifact`, `stable_event_ordinal`.

- Each row is one standardized event within an eligible ICU stay. Event time is finite, nonnegative minutes from ICU admission. Saturation is a finite percentage in the source's verified unit. Out-of-range events may be retained with `below_70` or `above_100` status for audit but cannot enter pairing.
- `range_status` is `below_70`, `eligible_70_100`, or `above_100`, consistent with the numeric value. `analysis_eligible_range` is exactly `TRUE` for `eligible_70_100` and `FALSE` otherwise.
- `source_family` is a nonblank locally defined source label. `source_priority` is a positive integer following the applicable source rule. These are **already resolved upstream**; the generic adapter does not infer an arterial specimen or choose a database field.
- `stable_event_ordinal` is a positive integer unique within each stay and modality, assigned by a deterministic upstream ordering rule. The same ordinal may occur in different stays or modalities.
- For SaO₂, `artifact_qc_status` is `not_applicable_non_spo2` and `possible_transient_artifact` is `FALSE`. For SpO₂, the flag is `TRUE` exactly when status is `possible_transient_artifact`; other status strings reflect the upstream temporal-resolution assessment.
- Every event's patient–stay mapping must match `stays.tsv`. An empty file, missing column, invalid numeric or logical value, unknown stay, or inconsistent mapping is an explicit error rather than an empty result.

The generic reader validates this interface, then the existing exact one-to-one pairing core independently forms 60-minute and 5-minute pairs. It does not apply the source-specific ICU selection, sepsis definition, specimen checks, unit conversions, or original event-cleaning algorithms. Stay- and pair-level clinical inputs for the mixed-effects and random-forest analyses are defined separately in `STANDARDIZED_CLINICAL_INPUT_V1.md`. Invented demonstration output is not a reproduction of the paper's observed estimates or approved displays.
