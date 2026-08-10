#!/bin/bash
###############################################################################
# fig3_02_compute_LD.sh   (mirrors the manuscript figure3/02_compute_fig3_LD.sh)
#
# Computes LD r^2 between each locus's index SNP and every SNP in the +/-500 kb
# window, from the 1000G Phase 3 EUR reference, with PLINK. One .ld file per
# locus for the plot step (fig3_03_plot.R) to colour by r^2. No internet needed.
#
# Index SNP = coloc.abf candidate shared variant (or, for SCARB2, the PD
# sentinel, so its two panels visibly peak on different variants).
#
# Run on narval (login or compute):
#   bash fig3_02_compute_LD.sh
###############################################################################
set -euo pipefail

LD_REF="/home/lchang24/projects/def-gsarah/lchang24/data/1kg/EUR"     # 1000G Ph3 EUR .bed/.bim/.fam
BASE="/home/lchang24/projects/def-gsarah/lchang24/github/coloc-abf"
OUT_DIR="${BASE}/ld"
mkdir -p "${OUT_DIR}"
WINDOW=500000

module load StdEnv/2020 plink/1.9b_6.21 2>/dev/null || module load plink/1.9b_6.21 2>/dev/null || true

echo "Reference: ${LD_REF}"
echo "Output:    ${OUT_DIR}"

compute_ld() {
  local name="$1" chr="$2" center="$3" index_rsid="$4"
  local from=$((center - WINDOW)) to=$((center + WINDOW)); [ "${from}" -lt 0 ] && from=0
  echo ""
  echo "--- ${name}: index ${index_rsid} (chr${chr}:${center}) ---"
  plink --bfile "${LD_REF}" \
    --chr "${chr}" --from-bp "${from}" --to-bp "${to}" \
    --ld-snp "${index_rsid}" --r2 \
    --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 0 \
    --out "${OUT_DIR}/${name}" --memory 4000 --threads 1 2>/dev/null || true
  if [ -f "${OUT_DIR}/${name}.ld" ]; then
    echo "  -> $(($(wc -l < "${OUT_DIR}/${name}.ld") - 1)) SNP pairs with r^2"
  else
    echo "  -> WARNING: LD failed. Is ${index_rsid} in the reference?"
    echo "     check:  grep -w ${index_rsid} ${LD_REF}.bim"
  fi
  rm -f "${OUT_DIR}/${name}.nosex" "${OUT_DIR}/${name}.log"
}

# name    chr  window-center (matches 07's extraction centers)  index rsID
compute_ld "RAB19"  7 140117809 "rs12216582"
compute_ld "MSRA"   8 10264000  "rs6982308"
compute_ld "BIN3"   8 22501830  "rs13264187"    # used by both BIN3 panels (eGFR + eGFRcys)
compute_ld "TNK2"   3 195624393 "rs13059257"
compute_ld "SCARB2" 4 77108306  "rs4859429"     # PD sentinel (distinct control)

echo ""
echo "Done. LD files:"
ls -la "${OUT_DIR}"/*.ld 2>/dev/null || echo "No .ld files — check warnings above"