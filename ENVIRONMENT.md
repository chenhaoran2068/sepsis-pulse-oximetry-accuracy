# Tested software environment

The synthetic demonstration was checked on Windows with R 4.5.1 (ucrt) and Python 3. The observed R packages and versions were: `data.table` 1.17.8, `igraph` 2.2.1, `digest` 0.6.39, `lubridate` 1.9.4, `stringr` 1.6.0, `lme4` 1.1-37, `Matrix` 1.7-3, `nlme` 3.1-168, `ranger` 0.17.0, `fastshap` 0.1.1, `mice` 3.18.0, `miceadds` 3.20.10, `car` 3.1-3, and `sandwich` 3.1-1. These are observed versions, not a tested compatibility range.

Install the needed packages in a user-controlled R library. The repository does not include a lockfile or automatically restore dependencies, and no non-Windows system has yet passed a clean run. `python demo/run_demo.py --rscript PATH` accepts an explicit Rscript executable path when it is not on `PATH`. No local library, Python environment, or clinical data are bundled.
