#!/bin/bash
#SBATCH --job-name=coloc_conjfdr
#SBATCH --time=24:00:00
#SBATCH --account=def-gsarah
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=conjfdr_clean/logs/coloc_conjfdr_%j.out
#SBATCH --error=conjfdr_clean/logs/coloc_conjfdr_%j.err

module load StdEnv/2023
module load r/4.3.1

export R_LIBS=~/.local/R/4.3.1/

cd /home/lchang24/links/projects/def-gsarah/lchang24/github/pleiofdr

# Check coloc is installed
Rscript -e 'if (!require("coloc")) quit(status=1)' 2>/dev/null
if [ $? -ne 0 ]; then
    echo "ERROR: coloc not installed. Install with:"
    echo "  Rscript -e 'install.packages(\"coloc\", repos=\"https://cloud.r-project.org\", lib=\"~/.local/R/4.3.1/\")'"
    exit 1
fi

SUPP_DIR="conjfdr_clean/supp_tables"
SUMSTATS_DIR="/home/lchang24/links/projects/def-gsarah/lchang24/github/nf-genetic-correlations-260504/nf-genetic-correlations/data/sumstats"
OUTDIR="conjfdr_clean/coloc_results_noMHC"

mkdir -p ${OUTDIR}

echo "============================================"
echo "Running coloc on conjFDR loci"
echo "Supp tables: ${SUPP_DIR}"
echo "Sumstats: ${SUMSTATS_DIR}"
echo "Output: ${OUTDIR}"
echo "============================================"
echo ""

Rscript run_coloc_conjfdr_loci_rmv_MHC.R "${SUPP_DIR}" "${SUMSTATS_DIR}" "${OUTDIR}"
