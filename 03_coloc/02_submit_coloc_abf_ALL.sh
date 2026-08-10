#!/bin/bash
###############################################################################
# 02_submit_coloc_abf_ALL.sh  —  coloc.abf (LD-free), FULL RUN
# narval; account def-gsarah.
#
# ALL pairs: every kidney trait x PD in the combined loci file
#   (eGFR, uACR, uK/Cr, uNa/Cr, hematuria, eGFRcys  x  PD, incl. sex strata).
#
# The number of loci is read FROM THE FILE at submit time, so you never
# hard-code an array size. Submit with the run command in the guide:
#   N=$(( $(wc -l < "$LOCI") - 1 ));  sbatch --array=1-${N}%20  <this script>
###############################################################################

#SBATCH --account=def-gsarah
#SBATCH --job-name=coloc_abf_all
#SBATCH --array=1-175%20
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/lchang24/projects/def-gsarah/lchang24/github/coloc-abf/logs/coloc_abf_all_%A_%a.out
#SBATCH --error=/home/lchang24/projects/def-gsarah/lchang24/github/coloc-abf/logs/coloc_abf_all_%A_%a.err

# ========================= active settings ==================================
REPO="/home/lchang24/projects/def-gsarah/lchang24/github/coloc-abf"
SUMSTATS_DIR="/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen"

LOCI_FILE="${REPO}/scripts/SuppTable_conjfdr_all_loci_with_eGFRcys.tsv"  # combined loci
RSCRIPT="${REPO}/scripts/run_coloc_abf_single_locus.R"                   # main analysis
OUTPUT_DIR="${REPO}/results_abf_all"
# ============================================================================

mkdir -p "${OUTPUT_DIR}" "${REPO}/logs"
module load StdEnv/2023 r/4.4.0            # must match the R your packages were built under

echo "=== Task ${SLURM_ARRAY_TASK_ID} | $(hostname) | $(date) ==="
echo "loci : ${LOCI_FILE}"
echo "out  : ${OUTPUT_DIR}"
Rscript "${RSCRIPT}" "${SLURM_ARRAY_TASK_ID}" "${LOCI_FILE}" "${SUMSTATS_DIR}" "${OUTPUT_DIR}"
echo "Done: $(date)"
