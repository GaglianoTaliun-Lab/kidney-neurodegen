#!/usr/bin/env bash
# run_all_loci.sh — run 04_define_loci_genes.py over every pair.
# Each pair runs in seconds, so a plain serial loop is fine (no SLURM array needed).
#
# Usage:
#   module load StdEnv/2020 plink/1.9b_6.21
#   cd ~/links/projects/def-gsarah/lchang24/github/pd-kidney-nc
#   bash scripts/run_all_loci.sh
set -euo pipefail

if ! command -v plink >/dev/null 2>&1; then
  echo "ERROR: plink not on PATH. Run:  module load StdEnv/2020 plink/1.9b_6.21" >&2
  exit 1
fi

BASE=~/links/projects/def-gsarah/lchang24
CONJ_BASE=$BASE/github/pleiofdr/conjfdr_clean/results
EUR=$BASE/data/EUR
GENES=$BASE/github/nf-genetic-correlations-260504/nf-genetic-correlations/protein_coding_genes_grch37.tsv
SCRIPT=scripts/04_define_loci_genes.py
OUTDIR=results/fuma_loci
mkdir -p "$OUTDIR"


PAIRS=(
  PD_metaSS_noProxy__eGFR_meta
  PD_metaSS_noProxy__uACR_sexComb
  PD_metaSS_noProxy__hematuria_sexComb
  PD_metaSS_noProxy__uKCr_sexComb
  PD_metaSS_noProxy__uNaCr_sexComb
  PD_female_noProxy__eGFR_female
  PD_female_noProxy__hematuria_female
  PD_male_noProxy__eGFR_male
  PD_male_noProxy__uNaCr_male
)

SUMMARY="$OUTDIR/_run_summary.tsv"
echo -e "pair\tsig_snps\tmhc_dropped\tmatched\tmerged_loci\tgenes_candidate\tgenes_ext500kb\tstatus" > "$SUMMARY"

for PAIR in "${PAIRS[@]}"; do
  FILEPREFIX="${PAIR/__/_}"   
  CSV="$CONJ_BASE/$PAIR/${FILEPREFIX}_zscore_conjfdr_0.05_all.csv"

  echo "===== $PAIR ====="
  if [ ! -f "$CSV" ]; then
    echo "  MISSING input: $CSV"
    echo -e "${PAIR}\tNA\tNA\tNA\tNA\tNA\tNA\tMISSING_INPUT" >> "$SUMMARY"
    continue
  fi

  LOG="$OUTDIR/${PAIR}.log"
  if python3 "$SCRIPT" \
        --conjfdr-csv "$CSV" \
        --bfile "$EUR" \
        --gene-file "$GENES" \
        --pair "$PAIR" \
        --outdir "$OUTDIR" \
        --exclude-mhc 2> "$LOG"; then
    # scrape the numbers from the stderr log
    sig=$(grep -oP 'Significant SNPs.*?: \K[0-9]+' "$LOG" || echo NA)
    mhc=$(grep -oP 'MHC dropped: \K[0-9]+' "$LOG" || echo NA)
    matched=$(grep -oP 'Matched to .*by position: \K[0-9]+' "$LOG" || echo NA)
    loci=$(grep -oP 'merged loci: \K[0-9]+' "$LOG" || echo NA)
    gcand=$(grep -oP 'genes \(candidate span\): \K[0-9]+' "$LOG" || echo NA)
    gext=$(grep -oP 'genes \(\+/-500kb\): \K[0-9]+' "$LOG" || echo NA)
    echo -e "${PAIR}\t${sig}\t${mhc}\t${matched}\t${loci}\t${gcand}\t${gext}\tOK" >> "$SUMMARY"
    echo "  OK: $loci loci, $gcand genes (candidate), $gext genes (ext500kb)"
  else
    echo "  FAILED — see $LOG"
    echo -e "${PAIR}\tNA\tNA\tNA\tNA\tNA\tNA\tFAILED" >> "$SUMMARY"
  fi
done

echo
echo "Done. Summary table:"
column -t -s $'\t' "$SUMMARY"
