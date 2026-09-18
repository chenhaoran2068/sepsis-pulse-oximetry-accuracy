# Data access and scope

This code candidate distributes **no real patient, ICU-stay, event, pair, or model-output records**. Running the demonstration generates invented identifiers, saturation measurements, and clinical values from programmatic examples. They were not sampled from, perturbed from, or calibrated to the five clinical cohorts. Demonstration statistics are interface checks, not estimates in the manuscript.

The study used MIMIC-IV 3.1, AmsterdamUMCdb 1.0.2, eICU Collaborative Research Database 2.0, SICdb 1.0.8, and a local Lianyungang cohort. Users who are independently eligible to access a source database must obtain it from its own custodian under the applicable conditions. This code does not grant access. The Lianyungang source data are not public; requests can be directed to the paper's corresponding authors. The repository includes neither credentials nor a copy of any source dataset.

The documented code boundary starts with already harmonized standardized inputs. Database-specific extraction, sepsis-cohort construction, quality review of raw source fields, and actual clinical-data processing remain outside this candidate. Generated `runs/` files are local test output and are not distributed as study data.
