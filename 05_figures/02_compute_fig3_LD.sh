#!/bin/bash
###############################################################################
# 02_compute_fig3_LD.sh
#
# Computes LD r^2 between each locus's index SNP and all SNPs in the window,
# using the 1000G Phase 3 EUR reference. Output: one .ld file per locus that
# locuszoomr reads for r^2 coloring.
#
# Run on Rorqual:
#   bash 02_compute_fig3_LD.sh
###############################################################################

set -euo pipefail

LD_REF="${HOME}/links/projects/def-gsarah/lchang24/data/EUR"
OUT_DIR="${HOME}/links/projects/def-gsarah/lchang24/github/coloc-susie/figure3/ld"
mkdir -p "${OUT_DIR}"

WINDOW=500000

module load StdEnv/2020 plink/1.9b_6.21 2>/dev/null || \
  module load plink/1.9b_6.21 2>/dev/null || true

echo "=========================================="
echo "Computing LD for Figure 3 loci"
echo "Reference: ${LD_REF}"
echo "Output:    ${OUT_DIR}"
echo "=========================================="

# index SNP rsIDs (these are the top shared causal variants from coloc)
# RAB19: rs12216582  (7:140125361)
# MSRA:  rs6982308   (8:10193772)
# BIN3:  rs7005025   (8:22519815)
# TNK2:  rs13059257  (3:195634993)

compute_ld() {
  local name="$1" chr="$2" center="$3" index_rsid="$4"
  local from=$((center - WINDOW)) to=$((center + WINDOW))
  [ "${from}" -lt 0 ] && from=0

  echo ""
  echo "--- ${name}: index ${index_rsid} (chr${chr}:${center}) ---"

  # Compute r^2 between index SNP and every SNP in the window.
  # --r2 with --ld-snp gives pairwise r^2 to the index variant.
  plink --bfile "${LD_REF}" \
    --chr "${chr}" \
    --from-bp "${from}" \
    --to-bp "${to}" \
    --ld-snp "${index_rsid}" \
    --r2 \
    --ld-window-kb 1000 \
    --ld-window 999999 \
    --ld-window-r2 0 \
    --out "${OUT_DIR}/${name}" \
    --memory 4000 --threads 1 \
    2>/dev/null

  if [ -f "${OUT_DIR}/${name}.ld" ]; then
    local n=$(($(wc -l < "${OUT_DIR}/${name}.ld") - 1))
    echo "  -> ${n} SNP pairs with r^2"
  else
    echo "  -> WARNING: LD computation failed (index SNP may not be in reference)"
    echo "     Check that ${index_rsid} exists in ${LD_REF}.bim"
  fi

  # clean up
  rm -f "${OUT_DIR}/${name}.nosex" "${OUT_DIR}/${name}.log"
}

compute_ld "RAB19" 7 140125361 "rs12216582"
compute_ld "MSRA"  8 10193772  "rs6982308"
compute_ld "BIN3"  8 22519815  "rs7005025"
compute_ld "TNK2"  3 195634993 "rs13059257"

echo ""
echo "=========================================="
echo "Done. LD files:"
ls -la "${OUT_DIR}"/*.ld 2>/dev/null || echo "No .ld files produced — check warnings above"
echo "=========================================="
echo ""
echo "NOTE: If any index SNP was not found in the reference, find a proxy:"
echo "  grep <rsid> ${LD_REF}.bim"
echo "  Or use the chr:pos to find the rsID in the reference at that position."
