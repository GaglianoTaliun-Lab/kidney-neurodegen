#!/bin/bash
###############################################################################
# 05_submit_sensitivity.sh  —  coloc.abf PRIOR-SENSITIVITY sweep
# narval; account def-gsarah.
#
# Sweeps p12 = {1e-6, 5e-6(primary), 1e-5(coloc default), 5e-5, Frida 1e-4}
# for the loci that actually carry a call, so R2's "is it prior-dependent?"
# is answered on exactly the loci we claim:
#
#   colocalized / Yang rows :  1  6  20  24  65  71  72  79  148
#     1  ZZZ3  (eGFR female)      65  RAB19 (eGFR meta)
#     6  TNK2  (eGFR female)      71  MSRA  (eGFR meta)
#     20 RAB19 (eGFR male)        72  PRSS51(eGFR meta)
#     24 BIN3  (eGFR male)        79  BIN3  (eGFR meta)
#                                 148 BIN3  (eGFRcys meta)
#   distinct control        :  53  SCARB2 (eGFR meta)  -- should stay distinct
###############################################################################

#SBATCH --account=def-gsarah
#SBATCH --job-name=coloc_abf_sens
#SBATCH --array=1,6,20,24,53,65,71,72,79,148
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/lchang24/projects/def-gsarah/lchang24/github/coloc-abf/logs/coloc_abf_sens_%A_%a.out
#SBATCH --error=/home/lchang24/projects/def-gsarah/lchang24/github/coloc-abf/logs/coloc_abf_sens_%A_%a.err

# ========================= active settings ==================================
REPO="/home/lchang24/projects/def-gsarah/lchang24/github/coloc-abf"
SUMSTATS_DIR="/home/lchang24/projects/def-gsarah/sumstats/kidney_neurodegen"

LOCI_FILE="${REPO}/scripts/SuppTable_conjfdr_all_loci_with_eGFRcys.tsv"   # same file/indices as the main run
RSCRIPT="${REPO}/scripts/run_coloc_abf_sensitivity.R"                    # prior sweep
OUTPUT_DIR="${REPO}/results_abf_sensitivity"
# ============================================================================

mkdir -p "${OUTPUT_DIR}" "${REPO}/logs"
module load StdEnv/2023 r/4.4.0

echo "=== Sensitivity task ${SLURM_ARRAY_TASK_ID} | $(hostname) | $(date) ==="
Rscript "${RSCRIPT}" "${SLURM_ARRAY_TASK_ID}" "${LOCI_FILE}" "${SUMSTATS_DIR}" "${OUTPUT_DIR}"
echo "Done: $(date)"
