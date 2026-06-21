#!/bin/bash
###############################################################################
# extract_fig3_windows_rorqual.sh
#
# Extracts +/-500 kb windows around each Figure 3 index SNP from GWAS sumstats.
# RORQUAL paths.
#
# Run:
#   bash extract_fig3_windows_rorqual.sh
###############################################################################

set -euo pipefail

SUMSTATS_DIR="${HOME}/links/projects/def-gsarah/lchang24/github/nf-genetic-correlations-260504/nf-genetic-correlations/data/sumstats"
OUT_DIR="${HOME}/links/projects/def-gsarah/lchang24/github/coloc-susie/figure3/windows"
mkdir -p "${OUT_DIR}"

WINDOW=500000

echo "=========================================="
echo "Extracting Figure 3 windows (+/- ${WINDOW} bp)"
echo "Sumstats: ${SUMSTATS_DIR}"
echo "Output:   ${OUT_DIR}"
echo "=========================================="

extract_locus() {
  local name="$1" chr="$2" center="$3" pd_file="$4" egfr_file="$5"
  local from=$((center - WINDOW)) to=$((center + WINDOW))
  [ "${from}" -lt 0 ] && from=0

  echo ""
  echo "--- ${name}: chr${chr}:${from}-${to} ---"

  awk -F'\t' -v c="${chr}" -v lo="${from}" -v hi="${to}" \
    'NR==1 || ($2==c && $3>=lo && $3<=hi)' \
    "${SUMSTATS_DIR}/${pd_file}" > "${OUT_DIR}/${name}_PD.tsv"
  echo "  PD   (${pd_file}): $(($(wc -l < "${OUT_DIR}/${name}_PD.tsv") - 1)) SNPs"

  awk -F'\t' -v c="${chr}" -v lo="${from}" -v hi="${to}" \
    'NR==1 || ($2==c && $3>=lo && $3<=hi)' \
    "${SUMSTATS_DIR}/${egfr_file}" > "${OUT_DIR}/${name}_eGFR.tsv"
  echo "  eGFR (${egfr_file}): $(($(wc -l < "${OUT_DIR}/${name}_eGFR.tsv") - 1)) SNPs"
}

extract_locus "RAB19" 7 140125361 "PD_metaSS_noProxy.tsv" "eGFR_meta.tsv"
extract_locus "MSRA"  8 10193772  "PD_metaSS_noProxy.tsv" "eGFR_meta.tsv"
extract_locus "BIN3"  8 22519815  "PD_metaSS_noProxy.tsv" "eGFR_meta.tsv"
extract_locus "TNK2"  3 195634993 "PD_female_noProxy.tsv" "eGFR_female.tsv"

echo ""
echo "Done. Files:"
ls -la "${OUT_DIR}"
