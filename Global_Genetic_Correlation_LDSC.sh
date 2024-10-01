#!/bin/bash

# Define variable and set paths

# Path to summary statisctics files
TRAIT1="/path_to_sumstats/PD_males.tsv.gz"          # Summary statisctic of neurodegenrative disease (AD, PD or Noproxy PD)
TRAIT2="/path_to_sumstats/hematuria_males.tsv.gz"   # Summary statistic of Kidney traits

# Path to output files     
TRAIT1_MUNGE="/path_to_munge_output/sumstats_reformatted/PD_males"             # Reformatted neurodegenrative disease
TRAIT2_MUNGE="/path_to_munge_output/sumstats_reformatted/hematuria_males"      # Reformatted kidney traits
OUTPUT_RG="/path_to_ldsc_output/ldscore_regression/PD_hematuria_male"          # Genetic correlation output from ldscy.py

# Path to LD Score files and HapMap3 SNP list
LDSCORE_DIR="/path_to_ldsc/ldsc"                        # Path to LDSC directory
REF_LD_SCORE="/path_to_ldscores/eur_w_ld_chr"           # Reference LD score path
HM3_SNPLIST="/path_to_ldscores/w_hm3.snplist"           # HapMap3 SNPs file used for filtering


# Reformatting the summary statistics files
python /LDSCORE_DIR/munge_sumstats.py \
    --sumstats /TRAIT1 \
    --out /TRAIT1_MUNGE \
    --merge_alleles /HM3_SNPLIST \
    --signed-sumstats beta,0 \
    --a1 effect_allele \
    --a2 other_allele

# Estimating Genetic Correlation
python /LDSCORE_DIR/ldsc.py \
    --rg /TRAIT1_MUNGE,/TRAIT2_MUNGE \
    --ref-ld-chr /REF_LD_SCORE \
    --w-ld-chr /REF_LD_SCORE \
    --out OUTPUT_RG
