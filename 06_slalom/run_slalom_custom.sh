#!/bin/bash
###############################################################################
# run_slalom_custom.sh -- run SLALOM with a CUSTOM UK Biobank LD panel
# (no gnomAD, no GCP). Consumes manifest_custom.tsv from build_ukb_custom_ld.py.
#
#   Usage: ./run_slalom_custom.sh <manifest_custom.tsv> <results_dir>
#
# Each locus is run with EXACT allele matching (the panel's variant index stores
# ref=allele1/alt=allele2 exactly as the .snp, and the BM is aligned to allele2),
# so --align-alleles is NOT used and no gnomAD sites table is needed.
###############################################################################
set -uo pipefail
MANIFEST="${1:?usage: run_slalom_custom.sh <manifest_custom.tsv> <results_dir>}"
RESULTS="${2:?usage: run_slalom_custom.sh <manifest_custom.tsv> <results_dir>}"
mkdir -p "$RESULTS"
R2="${R2:-0.6}"; NLOG10P="${NLOG10P:-4}"
command -v slalom >/dev/null 2>&1 || { echo "ERROR: 'slalom' not on PATH. pip install slalom-qc"; exit 1; }

n_ok=0; n_fail=0; first=1
# cols: snp_file bm_path vi_path locus_idx chr pos gene trait trait_type n_panel
while IFS=$'\t' read -r snp bm vi locus_idx chr pos gene trait ttype npanel; do
  if [[ "$first" == "1" ]]; then first=0; continue; fi
  [[ -z "$snp" ]] && continue
  base="$(basename "${snp%.snp}")"
  out="$RESULTS/${base}.slalom.txt"; summ="$RESULTS/${base}.summary.txt"
  cc=(); [[ "$ttype" == "cc" ]] && cc=(--case-control)
  echo ">> locus $locus_idx  $gene  [$trait]  (UKB panel: $npanel variants)"
  if slalom --snp "$snp" --out "$out" --out-summary "$summ" \
        --reference-genome GRCh37 \
        --ld-reference custom --custom-ld-path "$bm" \
        --custom-ld-variant-index-path "$vi" --custom-ld-label ukb \
        --dentist-s --abf --summary \
        --r2-threshold "$R2" --nlog10p-dentist-s-threshold "$NLOG10P" \
        "${cc[@]}" > "$RESULTS/${base}.log" 2>&1; then
    printf 'locus_idx\tgene\ttrait\ttrait_type\tchr\tpos\n' > "$RESULTS/${base}.meta"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$locus_idx" "$gene" "$trait" "$ttype" "$chr" "$pos" >> "$RESULTS/${base}.meta"
    n_ok=$((n_ok+1)); echo "   OK -> $summ"
  else
    n_fail=$((n_fail+1)); echo "   FAILED (see $RESULTS/${base}.log):"; tail -3 "$RESULTS/${base}.log" | sed 's/^/      /'
  fi
done < "$MANIFEST"

echo; echo "Done: $n_ok ok, $n_fail failed."
echo "Next: Rscript collect_slalom_summaries.R $RESULTS slalom_ukb_supplementary_table.tsv"
