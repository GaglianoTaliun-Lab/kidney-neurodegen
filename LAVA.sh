#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --time=70:00:00
#SBATCH --job-name=PD_sod_male
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=sadaf.gawhary@umontreal.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G

# Load necessary modules
module load r/4.4.0
module load gcc/12.3

# Set R library path
export R_LIBS_USER=/scratch/sadafgy/R/x86_64-pc-linux-gnu-library/4.4

# Define project and script directories
PROJECT_DIR="/home/sadafgy/projects/def-gsarah/sadafgy"
OUT_DIR="$PROJECT_DIR/LAVA/lava_results"
SCRIPT_PATH="/scratch/sadafgy/PD_sod_male.R"

# Create output directory if it doesn't exist
mkdir -p $OUT_DIR

# Run the R script with the current array index as the block number
Rscript $SCRIPT_PATH 
