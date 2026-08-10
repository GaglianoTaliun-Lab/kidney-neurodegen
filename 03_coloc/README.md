# Stage 03 — Bayesian colocalization (coloc.abf)

This stage tests whether PD and a kidney-related trait share a single causal variant at each
non-MHC conjFDR-identified locus, using Bayesian colocalization
(`coloc.abf`). `coloc.abf` derives colocalization posteriors directly from each trait's
marginal association statistics under a single-causal-variant assumption. 
Prior probabilities are set conservatively
(p1 = p2 = 1×10⁻⁴, p12 = 5×10⁻⁶); the sensitivity of the calls to
this choice is assessed across a grid of priors.

## Software used

- **coloc** (R package, v5.2.3) — `coloc.abf`. https://github.com/chr1swallace/coloc (CRAN package `coloc`)
- **R 4.3.1**

## Scripts (run in order)

- `coloc_abf_common.R` — shared configuration and per-locus data loading (sourced by the runners).
- `07_extract_locus_regions.R` — extracts the ±500 kb summary-statistic window per locus and trait pair.
- `run_coloc_abf_single_locus.R` — runs `coloc.abf` at one locus for a trait pair (SLURM-array element).
- `02_submit_coloc_abf_ALL.sh` — submits `coloc.abf` across all non-MHC conjFDR loci.
- `03_aggregate_coloc_abf.R` — collates per-locus posteriors into the master results table.
- `04_dedup_loci.R` — collapses loci within 500 kb of each other into distinct genomic loci for reporting.
- `run_coloc_abf_sensitivity.R` + `05_submit_sensitivity.sh` — repeat `coloc.abf` across a grid of priors (prior-sensitivity; Supplementary Table).
