# Stage 06 — SLALOM (DENTIST-S) LD-consistency diagnostic (Supplementary)

Because both discovery GWAS are themselves meta-analyses, and fine-mapping of meta-analysis
summary statistics can be miscalibrated, this stage screens every reported colocalized locus
for the summary-statistic/LD inconsistency that flags this problem, using **SLALOM**
(Kanai et al., 2022) with its simplified DENTIST test (DENTIST-S). Following Kanai et al.,
the LD reference is used **only** to diagnose summary-statistic/LD inconsistency, not for
fine-mapping. Because a gnomAD reference was not available, an in-sample **UK Biobank** LD
reference is built from the genotyping array (unrelated individuals of white-British ancestry).

## Software used

- **SLALOM** / **DENTIST-S** (Kanai et al., 2022) [confirm version/commit]
- **Python 3.11**, **R 4.3.1**, **PLINK 1.9** (1.9b_6.21)

## Scripts

- `build_ukb_custom_ld.py` — builds the in-sample UK Biobank LD reference for each reported locus.
- `build_slalom_snp.R` — prepares the per-locus `.snp` input file for SLALOM.
- `run_slalom_loci.sh` / `run_slalom_custom.sh` — run SLALOM/DENTIST-S at each reported colocalized locus and at the *SCARB2/FAM47E-STBD1* distinct-variant control, for both traits.
- `collect_slalom_summaries.R` — collates the per-locus SLALOM output into the Supplementary Table.
