#!/bin/bash

# Define directories and common paths
SUMSTATS_DIR="/home/sadafgy/projects/def-gsarah/sadafgy/ldsc_results/sumstats_reformatted"
OUTPUT_DIR="/home/sadafgy/projects/def-gsarah/sadafgy/ldsc_results/heritability"
LDSCORE_DIR="/home/sadafgy/projects/def-gsarah/group_writable/ldsc_ldscores/eur_w_ld_chr"

# Path to the ldsc.py script
LDSC_SCRIPT="/home/sadafgy/projects/def-gsarah/sadafgy/ldsc/ldsc.py"

# Loop through all `.sumstats.gz` files in the directory
for SUMSTATS_FILE in "$SUMSTATS_DIR"/*.sumstats.gz; do
  # Extract the base name of the file (e.g., "AD_sexcombined" from "AD_sexcombined.sumstats.gz")
  BASE_NAME=$(basename "$SUMSTATS_FILE" .sumstats.gz)

  # Define output file path
  OUTPUT_FILE="$OUTPUT_DIR/${BASE_NAME}"

  # Run ldsc.py
  python "$LDSC_SCRIPT" \
    --h2 "$SUMSTATS_FILE" \
    --ref-ld-chr "$LDSCORE_DIR/" \
    --w-ld-chr "$LDSCORE_DIR/" \
    --out "$OUTPUT_FILE"

  # Print a message for progress tracking
  echo "Processed: $SUMSTATS_FILE -> $OUTPUT_FILE"
done

echo "All files processed!"
