# Stage 03 — Bayesian colocalization (coloc.susie / coloc.abf)
 
This stage tests whether PD and a kidney-related trait share a single causal variant at
each conjFDR-identified locus, using Bayesian colocalization. 
 
## Software used
 
- **coloc** (R) — Bayesian colocalization (`coloc.susie` and `coloc.abf`).
  Repository: https://github.com/chr1swallace/coloc (CRAN package `coloc`).
- **susieR** (R) — Sum of Single Effects fine-mapping, used by `coloc.susie`.
  Repository: https://github.com/stephenslab/susieR (CRAN package `susieR`).
- **R 4.3.1**; **PLINK 1.9** (1.9b_6.21) for computing per-locus LD matrices.
