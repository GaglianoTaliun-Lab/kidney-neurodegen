#!/bin/bash
###############################################################################
# 02_submit_coloc_susie.sh
#
# SLURM array job: runs coloc.susie for all conjFDR loci in parallel.
#
# Before running:
#   1. Run 01_extract_ld_matrices.sh
#   2. Set correct array range:
#      tail -n +2 $LOCI_FILE | awk -F'\t' '!($5==6 && $6>=25e6 && $6<=34e6)' | wc -l
#
# Usage: sbatch 02_submit_coloc_susie.sh
###############################################################################

#SBATCH --account=def-gsarah
#SBATCH --job-name=coloc_susie
#SBATCH --array=17,47,48,49,50,51,71,72,75,76,122
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/coloc_susie_%A_%a.out
#SBATCH --error=logs/coloc_susie_%A_%a.err

# -------------- PATHS ------------------------------------------------------
BASE_DIR="${HOME}/links/projects/def-gsarah/lchang24"

LOCI_FILE="${BASE_DIR}/github/pleiofdr/conjfdr_clean/supp_tables/SuppTable_conjfdr_all_loci.tsv"
SUMSTATS_DIR="${BASE_DIR}/github/nf-genetic-correlations-260504/nf-genetic-correlations/data/sumstats"
LD_DIR="${BASE_DIR}/github/coloc-susie/ld_matrices"
OUTPUT_DIR="${BASE_DIR}/github/coloc-susie/results"
RSCRIPT="${BASE_DIR}/github/coloc-susie/scripts/run_coloc_susie_single_locus.R"
# ---------------------------------------------------------------------------

mkdir -p "${OUTPUT_DIR}" logs
module load StdEnv/2023 r/4.3.1

echo "=== Task ${SLURM_ARRAY_TASK_ID} | $(date) ==="
Rscript "${RSCRIPT}" "${SLURM_ARRAY_TASK_ID}" "${LOCI_FILE}" "${SUMSTATS_DIR}" "${LD_DIR}" "${OUTPUT_DIR}"
echo "Done: $(date)"
