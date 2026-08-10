#!/bin/bash
###############################################################################
# run_slalom_loci.sh  --  run SLALOM (slalom-qc) on every per-locus .snp file
# in a manifest, for Reviewer-2 Comment 1 (meta-analysis fine-mapping
# miscalibration; Kanai et al. 2022).
#
#   Usage: ./run_slalom_loci.sh <manifest.tsv> <results_dir>
#
# Reference data (gnomAD): the LD BlockMatrix is read from the PUBLIC bucket
# gs://gcp-public-data--gnomad (free, but needs GCS network). The LD variant
# index and (for --align-alleles) the gnomAD sites table are requester-pays.
# Pick ONE of the two modes below by setting env vars:
#
#   Mode A (requester-pays, simplest if you have a GCP project with billing):
#       export GCP_PROJECT=my-gcp-project
#     -> per-locus reads are billed to your project (cents total).
#
#   Mode B (local bundles, no per-read billing; BM still from public bucket):
#       export LD_INDEX_DIR=/path/to/ld_index          # extracted ldcov bundle
#       export SITES_PARQUET=/path/to/gnomad.genomes.r2.1.1.sites.most_severe.b37.parquet
#       export CUP_PARQUET=/path/to/FASTA_BED.ALL_GRCh37.cups.parquet   # only if ANNOTATE_CUPS=1
#     (download the bundles once with gcloud + a billing project; see the guide.)
#
# Toggles (defaults in []):
#   ALIGN=1        [1]  --align-alleles (reconcile ref/alt vs gnomAD; needs SITES).
#                        Set ALIGN=0 for the leanest run (LD index only, no sites),
#                        but then allele1/allele2 must already be correct b37 ref/alt.
#   ANNOTATE=0     [0]  add --annotate-consequence --annotate-gnomad-freq (needs SITES).
#   ANNOTATE_CUPS=0[0]  add --annotate-cups (needs CUP parquet).
#   R2=0.6         [0.6]  --r2-threshold  (Kanai default)
#   NLOG10P=4      [4]    --nlog10p-dentist-s-threshold  (P<1e-4, Kanai default)
###############################################################################
set -uo pipefail

MANIFEST="${1:?usage: run_slalom_loci.sh <manifest.tsv> <results_dir>}"
RESULTS="${2:?usage: run_slalom_loci.sh <manifest.tsv> <results_dir>}"
mkdir -p "$RESULTS"

ALIGN="${ALIGN:-1}"; ANNOTATE="${ANNOTATE:-0}"; ANNOTATE_CUPS="${ANNOTATE_CUPS:-0}"
R2="${R2:-0.6}"; NLOG10P="${NLOG10P:-4}"

command -v slalom >/dev/null 2>&1 || { echo "ERROR: 'slalom' not on PATH. pip install slalom-qc"; exit 1; }

# ---- assemble the reference-data flags for the chosen mode -------------------
REF_FLAGS=()
if [[ -n "${LD_INDEX_DIR:-}" ]]; then
  echo "== Reference mode B: local bundles =="
  REF_FLAGS+=(--ld-variant-index-dir "$LD_INDEX_DIR")
  [[ -n "${SITES_PARQUET:-}" ]] && REF_FLAGS+=(--gnomad-sites-parquet "$SITES_PARQUET")
  [[ -n "${CUP_PARQUET:-}"   ]] && REF_FLAGS+=(--cup-parquet "$CUP_PARQUET")
elif [[ -n "${GCP_PROJECT:-}" ]]; then
  echo "== Reference mode A: requester-pays (project=$GCP_PROJECT) =="
  REF_FLAGS+=(--storage-options "{\"requester_pays\": true, \"project\": \"$GCP_PROJECT\"}")
else
  echo "ERROR: set GCP_PROJECT (mode A) or LD_INDEX_DIR (mode B). See header."; exit 1
fi

# In local-bundle mode, --align-alleles needs the sites parquet; without it SLALOM would
# fall back to the hosted requester-pays sites path with no billing project and fail.
if [[ -n "${LD_INDEX_DIR:-}" && "$ALIGN" == "1" && -z "${SITES_PARQUET:-}" ]]; then
  echo "ERROR: ALIGN=1 in bundle mode needs SITES_PARQUET (the gnomAD sites parquet from the"
  echo "       slalom.b37 bundle). Set SITES_PARQUET, or set ALIGN=0 to skip allele alignment."
  exit 1
fi

# common statistic flags (the DENTIST-S verdict + ABF max-PIP + per-locus summary)
STAT_FLAGS=(--reference-genome GRCh37 --dentist-s --abf --summary
            --r2-threshold "$R2" --nlog10p-dentist-s-threshold "$NLOG10P")
[[ "$ALIGN" == "1" ]] && STAT_FLAGS+=(--align-alleles)
if [[ "$ANNOTATE" == "1" ]]; then STAT_FLAGS+=(--annotate-consequence --annotate-gnomad-freq); fi
[[ "$ANNOTATE_CUPS" == "1" ]] && STAT_FLAGS+=(--annotate-cups)

# ---- loop the manifest (tab-separated; trait names contain spaces) -----------
# cols: snp_file locus_idx chr pos gene trait trait_type n_variants
n_ok=0; n_fail=0; first=1
while IFS=$'\t' read -r snp locus_idx chr pos gene trait ttype nvar; do
  if [[ "$first" == "1" ]]; then first=0; continue; fi        # skip header
  [[ -z "$snp" ]] && continue
  base="$(basename "${snp%.snp}")"
  out="$RESULTS/${base}.slalom.txt"
  summ="$RESULTS/${base}.summary.txt"
  cc=(); [[ "$ttype" == "cc" ]] && cc=(--case-control)

  echo ">> locus $locus_idx  $gene  [$trait]  ($ttype, $nvar variants)"
  if slalom --snp "$snp" --out "$out" --out-summary "$summ" \
        "${STAT_FLAGS[@]}" "${cc[@]}" "${REF_FLAGS[@]}" > "$RESULTS/${base}.log" 2>&1; then
    # tag the summary with locus/trait provenance for the collector
    echo -e "locus_idx\tgene\ttrait\ttrait_type\tchr\tpos" > "$RESULTS/${base}.meta"
    echo -e "${locus_idx}\t${gene}\t${trait}\t${ttype}\t${chr}\t${pos}" >> "$RESULTS/${base}.meta"
    n_ok=$((n_ok+1)); echo "   OK -> $summ"
  else
    n_fail=$((n_fail+1)); echo "   FAILED (see $RESULTS/${base}.log):"; tail -3 "$RESULTS/${base}.log" | sed 's/^/      /'
  fi
done < "$MANIFEST"

echo
echo "Done: $n_ok ok, $n_fail failed. Summaries in $RESULTS/*.summary.txt"
echo "Next: Rscript collect_slalom_summaries.R $RESULTS slalom_supplementary_table.tsv"
