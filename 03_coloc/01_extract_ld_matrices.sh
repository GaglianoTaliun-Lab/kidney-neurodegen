#!/bin/bash
###############################################################################
# 01_extract_ld_matrices.sh
#
# Extracts per-locus LD matrices AND allele frequencies from 1000G Phase 3 EUR.
# Generates .ld (LD matrix), .snplist (SNP order), .bim_positions (positions),
# and .frq (MAF) for each locus.
#
# Input file columns (tab-separated):
#   1:pair 2:pd_stratum 3:kidney_trait 4:locus_num 5:chr 6:position
#   7:conjfdr 8:pval_PD 9:nearest_gene
#
# Usage:
#   bash 01_extract_ld_matrices.sh <loci_file> <ld_ref_prefix> <output_dir>
###############################################################################

set -euo pipefail

LOCI_FILE="${1:?Usage: $0 <loci_file> <ld_ref_prefix> <output_dir>}"
LD_REF="${2:?Provide 1000G EUR PLINK prefix}"
OUTPUT_DIR="${3:?Provide output directory}"

WINDOW=500000

mkdir -p "${OUTPUT_DIR}"

module load StdEnv/2020 plink/1.9b_6.21 2>/dev/null || true

# Deduplicate: extract unique chr:position pairs (exclude MHC)
tail -n +2 "${LOCI_FILE}" | \
  awk -F'\t' '!($5 == 6 && $6 >= 25000000 && $6 <= 34000000) {print $5"\t"$6}' | \
  sort -u > "${OUTPUT_DIR}/_unique_positions.txt"

N_UNIQUE=$(wc -l < "${OUTPUT_DIR}/_unique_positions.txt")
echo "============================================"
echo "Extracting LD matrices + MAF for ${N_UNIQUE} unique positions"
echo "Reference: ${LD_REF}"
echo "Window: +/- ${WINDOW} bp"
echo "Output: ${OUTPUT_DIR}"
echo "============================================"
echo ""

COUNT=0
while IFS=$'\t' read -r CHR POS; do
  COUNT=$((COUNT + 1))

  FROM=$((POS - WINDOW))
  TO=$((POS + WINDOW))
  [ "${FROM}" -lt 0 ] && FROM=0

  PREFIX="${OUTPUT_DIR}/ld_chr${CHR}_${POS}"

  # Skip if already done (check both .ld and .frq exist)
  if [ -f "${PREFIX}.ld" ] && [ -f "${PREFIX}.frq" ]; then
    echo "[${COUNT}/${N_UNIQUE}] chr${CHR}:${POS} — already exists, skipping"
    continue
  fi

  echo "[${COUNT}/${N_UNIQUE}] chr${CHR}:${FROM}-${TO}"

  # Extract LD matrix
  plink --bfile "${LD_REF}" \
    --chr "${CHR}" \
    --from-bp "${FROM}" \
    --to-bp "${TO}" \
    --r square \
    --write-snplist \
    --out "${PREFIX}" \
    --memory 4000 \
    --threads 1 \
    2>/dev/null

  # Extract allele frequencies for the SAME window
  plink --bfile "${LD_REF}" \
    --chr "${CHR}" \
    --from-bp "${FROM}" \
    --to-bp "${TO}" \
    --freq \
    --out "${PREFIX}" \
    --memory 4000 \
    --threads 1 \
    2>/dev/null

  # Save BIM positions for SNP-to-position mapping
  awk -v chr="${CHR}" -v from="${FROM}" -v to="${TO}" \
    '$1 == chr && $4 >= from && $4 <= to' \
    "${LD_REF}.bim" > "${PREFIX}.bim_positions"

  # Clean up PLINK temp files
  rm -f "${PREFIX}.nosex" "${PREFIX}.log"

  if [ -f "${PREFIX}.ld" ] && [ -f "${PREFIX}.frq" ]; then
    N_SNPS=$(wc -l < "${PREFIX}.snplist")
    echo "  -> ${N_SNPS} SNPs (LD + MAF extracted)"
  else
    echo "  -> WARNING: extraction failed for chr${CHR}:${POS}"
  fi
done < "${OUTPUT_DIR}/_unique_positions.txt"

rm -f "${OUTPUT_DIR}/_unique_positions.txt"

echo ""
echo "Done. LD files: $(ls ${OUTPUT_DIR}/*.ld 2>/dev/null | wc -l)"
echo "      FRQ files: $(ls ${OUTPUT_DIR}/*.frq 2>/dev/null | wc -l)"
