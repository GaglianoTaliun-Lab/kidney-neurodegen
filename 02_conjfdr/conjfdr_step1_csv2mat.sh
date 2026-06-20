#!/bin/bash
#SBATCH --job-name=conjfdr_step1_csv2mat
#SBATCH --time=06:00:00
#SBATCH --account=def-gsarah
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=conjfdr_clean/logs/step1_csv2mat_%j.out
#SBATCH --error=conjfdr_clean/logs/step1_csv2mat_%j.err

###############################################################################
# STEP 1: Convert all CSV sumstats to .mat format for pleioFDR
#
# Input:  csv_convert/*.csv (18 files)
# Output: conjfdr_clean/mat_files/*.mat (18 files)
###############################################################################

module load StdEnv/2023
module load python/3.11

BASEDIR=/home/lchang24/links/projects/def-gsarah/lchang24/github/pleiofdr
CSVDIR=${BASEDIR}/csv_convert
MATDIR=${BASEDIR}/conjfdr_clean/mat_files
CONVERT=${BASEDIR}/python_convert/sumstats.py

mkdir -p ${MATDIR}

echo "========================================"
echo "STEP 1: CSV -> MAT conversion"
echo "Started: $(date)"
echo "========================================"

# Verify python_convert exists
if [ ! -f "${CONVERT}" ]; then
    echo "ERROR: ${CONVERT} not found"
    exit 1
fi

TOTAL=0
FAIL=0

for csv in ${CSVDIR}/*.csv; do
    base=$(basename "$csv" .csv)
    matfile="${MATDIR}/${base}.mat"

    echo ""
    echo "--- Converting: ${base} ---"
    python3 "${CONVERT}" mat \
        --sumstats "${csv}" \
        --ref "${BASEDIR}/9545380.ref" \
        --out "${matfile}" \
        --force \
        --without-n

    if [ -f "${matfile}" ]; then
        SIZE=$(ls -lh "${matfile}" | awk '{print $5}')
        echo "  OK: ${SIZE}"
        TOTAL=$((TOTAL + 1))
    else
        echo "  FAILED: mat file not created"
        FAIL=$((FAIL + 1))
    fi
done

echo ""
echo "========================================"
echo "STEP 1 COMPLETE: ${TOTAL} converted, ${FAIL} failed"
echo "Finished: $(date)"
echo "========================================"
echo ""
echo "Mat files:"
ls -lh ${MATDIR}/*.mat

if [ ${FAIL} -eq 0 ] && [ ${TOTAL} -eq 18 ]; then
    echo ""
    echo "ALL 18 FILES CONVERTED. Ready for step 2:"
    echo "  sbatch conjfdr_step2_array.sh"
else
    echo ""
    echo "WARNING: Expected 18 files, got ${TOTAL}. Fix errors before step 2."
fi