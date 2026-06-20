#!/bin/bash
#SBATCH --job-name=conjfdr_step2
#SBATCH --time=04:00:00
#SBATCH --account=def-gsarah
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --array=1-15%5
#SBATCH --output=conjfdr_clean/logs/conjfdr_%A_%a.out
#SBATCH --error=conjfdr_clean/logs/conjfdr_%A_%a.err

###############################################################################
# STEP 2: Run conjFDR for all 15 sex-concordant PD × kidney pairs
#
# Pairs (3 PD × 5 kidney = 15):
#   PD_metaSS_noProxy  × eGFR_meta, uACR_sexComb, hematuria_sexComb, uKCr_sexComb, uNaCr_sexComb
#   PD_female_noProxy  × eGFR_female, uACR_female, hematuria_female, uKCr_female, uNaCr_female
#   PD_male_noProxy    × eGFR_male, uACR_male, hematuria_male, uKCr_male, uNaCr_male
#
# Each array task:
#   1. Reads its pair from PAIRS array
#   2. Generates config on the fly
#   3. Copies all pleioFDR files to $SLURM_TMPDIR (avoids symlink issues)
#   4. Runs MATLAB from $SLURM_TMPDIR
#   5. Results write to absolute output path
###############################################################################

module load StdEnv/2023
module load matlab/2023b.2

BASEDIR=/home/lchang24/links/projects/def-gsarah/lchang24/github/pleiofdr
MATDIR=${BASEDIR}/conjfdr_clean/mat_files
RESULTSDIR=${BASEDIR}/conjfdr_clean/results
CONFIGDIR=${BASEDIR}/conjfdr_clean/configs

mkdir -p ${RESULTSDIR} ${CONFIGDIR}

# ---- Define all 15 pairs ----
# Format: PD_TRAIT:KIDNEY_TRAIT
PAIRS=(
    "PD_metaSS_noProxy:eGFR_meta"
    "PD_metaSS_noProxy:uACR_sexComb"
    "PD_metaSS_noProxy:hematuria_sexComb"
    "PD_metaSS_noProxy:uKCr_sexComb"
    "PD_metaSS_noProxy:uNaCr_sexComb"
    "PD_female_noProxy:eGFR_female"
    "PD_female_noProxy:uACR_female"
    "PD_female_noProxy:hematuria_female"
    "PD_female_noProxy:uKCr_female"
    "PD_female_noProxy:uNaCr_female"
    "PD_male_noProxy:eGFR_male"
    "PD_male_noProxy:uACR_male"
    "PD_male_noProxy:hematuria_male"
    "PD_male_noProxy:uKCr_male"
    "PD_male_noProxy:uNaCr_male"
)

# ---- Get this task's pair ----
IDX=$((SLURM_ARRAY_TASK_ID - 1))
PAIR_STR=${PAIRS[$IDX]}
PD_TRAIT=$(echo $PAIR_STR | cut -d: -f1)
KIDNEY_TRAIT=$(echo $PAIR_STR | cut -d: -f2)
PAIR_NAME="${PD_TRAIT}__${KIDNEY_TRAIT}"

echo "========================================"
echo "Task ${SLURM_ARRAY_TASK_ID}: ${PAIR_NAME}"
echo "Started: $(date)"
echo "========================================"

# ---- Verify mat files exist ----
PD_MAT="${MATDIR}/${PD_TRAIT}.mat"
KIDNEY_MAT="${MATDIR}/${KIDNEY_TRAIT}.mat"

if [ ! -f "${PD_MAT}" ]; then
    echo "ERROR: PD mat file not found: ${PD_MAT}"
    exit 1
fi
if [ ! -f "${KIDNEY_MAT}" ]; then
    echo "ERROR: Kidney mat file not found: ${KIDNEY_MAT}"
    exit 1
fi

# ---- Create output directory ----
OUTDIR="${RESULTSDIR}/${PAIR_NAME}"
mkdir -p "${OUTDIR}"

# ---- Set up isolated workspace in SLURM_TMPDIR ----
WORKDIR=${SLURM_TMPDIR}/pleiofdr_${PAIR_NAME}
mkdir -p ${WORKDIR}/mat_files

# Copy pleioFDR MATLAB scripts
cp ${BASEDIR}/*.m ${WORKDIR}/

# Copy reference file
cp ${BASEDIR}/ref9545380_1kgPhase3eur_LDr2p1.mat ${WORKDIR}/

# Copy the two mat files needed for this pair
cp "${PD_MAT}" ${WORKDIR}/mat_files/
cp "${KIDNEY_MAT}" ${WORKDIR}/mat_files/

# ---- Generate config ----
cat > ${WORKDIR}/config.txt << EOF
reffile=ref9545380_1kgPhase3eur_LDr2p1.mat
traitfolder=.
traitfile1=mat_files/${PD_TRAIT}.mat
traitname1=${PD_TRAIT}
traitfiles={'mat_files/${KIDNEY_TRAIT}.mat'}
traitnames={'${KIDNEY_TRAIT}'}
outputdir=${OUTDIR}
randprune_n=500
stattype=conjfdr
fdrthresh=0.05
EOF

# Also save config for reference
cp ${WORKDIR}/config.txt ${CONFIGDIR}/config_${PAIR_NAME}.txt

echo "Config:"
cat ${WORKDIR}/config.txt
echo ""

# ---- Run MATLAB ----
cd ${WORKDIR}
echo "Running MATLAB from: ${WORKDIR}"
echo "Output to: ${OUTDIR}"
echo ""

matlab -nodisplay -nosplash < runme.m

echo ""
echo "========================================"
echo "Completed: ${PAIR_NAME}"
echo "Finished: $(date)"
echo "========================================"
echo ""
echo "Results:"
ls -lh ${OUTDIR}/
