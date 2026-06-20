# Genome-wide cross-trait analysis of Parkinson's disease and kidney-related traits

Analysis code and pipeline for the manuscript:

> **Genome-wide cross-trait analysis characterizes shared genetic architecture between Parkinson's disease and kidney-related traits**
> Le Chang\*, Sadaf Gawhary\*, Lyza Maameri, Wiame Belbellaj, Frida Lona-Durazo, Sarah A. Gagliano Taliun
> \*Equal contribution. Corresponding author: Sarah A. Gagliano Taliun (sarah.gagliano-taliun@umontreal.ca)
>

This repository contains the custom code used to characterize shared genetic
architecture between Parkinson's disease (PD) and five kidney-related traits —
estimated glomerular filtration rate (eGFR), urinary albumin-to-creatinine ratio
(uACR), hematuria, urinary potassium-to-creatinine ratio (uK/Cr), and urinary
sodium-to-creatinine ratio (uNa/Cr) — in individuals of European genetic ancestry.
All coordinates are on genome build GRCh37 (hg19).

## Repository structure

```
.
├── 01_ldsc/              # Cross-trait LD score regression 
├── 02_lava/              # Local genetic correlation across LD blocks 
├── 03_conjfdr/           # conjFDR / pleioFDR across 15 trait pairs
├── 04_coloc_susie/       # Bayesian colocalization 
├── 05_eqtl_annotation/   # eQTL lookup in GTEx v11 and NephQTL2
├── 06_enrichment/        # FUMA GENE2FUNC + NetworkAnalyst 3.0 
├── 07_figures/           # Figures
├── data/                 # Data pointers only (see data/README.md) -- no raw data committed
├── env/                  # Environment capture (renv.lock, sessionInfo, module list)
├── CITATION.cff
├── LICENSE
└── README.md
```

## Pipeline overview

1. **LDSC** (`01_ldsc/`) — genome-wide heritability and pairwise genetic correlation
   between PD and each kidney trait.
2. **LAVA** (`02_lava/`) — local genetic correlation across approximately independent
   LD blocks. 
3. **conjFDR / pleioFDR** (`03_conjfdr/`) — identification of pleiotropic loci jointly
   associated with both traits at conjFDR < 0.05. 
4. **coloc.susie** (`04_coloc_susie/`) — colocalization at conjFDR loci using
   conservative priors (p1 = p2 = 1e-4, p12 = 5e-6), with coloc.abf as a fallback
   where SuSiE did not converge. 
5. **eQTL annotation** (`05_eqtl_annotation/`) — candidate gene assignment for
   colocalized loci using GTEx v11 (brain tissues) and NephQTL2 (kidney).
6. **Enrichment** (`06_enrichment/`) — pathway enrichment (FUMA GENE2FUNC) and PPI
   network enrichment (NetworkAnalyst 3.0). 
7. **Figures** (`07_figures/`) — figures for the manuscrpt.

## Software environment

- **R 4.3.1**. Key packages: `data.table`, `ggplot2`, `ggrepel`, `patchwork`,
  `cowplot`, `locuszoomr` (0.3.8), `coloc`, `susieR`, `openxlsx`,
  `EnsDb.Hsapiens.v75`. Exact versions are pinned in `env/` (`renv.lock` /
  `sessionInfo.txt`).
- **PLINK 1.9** (1.9b_6.21) for LD operations.
- **MATLAB** for the pleioFDR/conjFDR step (see `03_conjfdr/`).
- Analyses were run on the Digital Research Alliance of Canada HPC under a SLURM
  scheduler. Module loads are recorded in `env/modules.txt`.

## License

Code is released under the MIT License (see `LICENSE`). The published article is
licensed separately under CC BY 4.0.

## Contact

Questions about the code: open an issue or contact the corresponding author,
Sarah A. Gagliano Taliun (sarah.gagliano-taliun@umontreal.ca).