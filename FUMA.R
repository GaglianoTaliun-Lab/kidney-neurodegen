## Step 1: Combine all the kidney traits combinations for each female, male, sex-combined and sex-combined meta-analyzed 

library(dplyr)
library(readr)

# Define paths
input_folder <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/lava_results/bivariate_results/NP_PD"  # Folder with LAVA results
output_file <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/FUMA/joined_lava_file/NP_PD_female.txt"  # Output file for significant loci
pattern <- c("NP_PD_female")  # Pattern to search for in filenames

# Get all files containing NP_PD_male
files <- list.files(input_folder, pattern = pattern, full.names = TRUE)

# Initialize empty data frame for results
significant_loci <- data.frame()

# Loop through files
for (file in files) {
  # Read the LAVA results file
  lava_data <- read_tsv(file, col_types = cols())
  
  # Ensure required columns exist
  required_cols <- c("locus", "chr", "start", "stop", "phen1", "phen2", "p")
  if (!all(required_cols %in% colnames(lava_data))) {
    cat("Skipping file due to missing columns:", file, "\n")
    next
  }
  
  # Calculate Bonferroni threshold
  n_loci <- nrow(lava_data)
  bonferroni_threshold <- 0.05/n_loci
  
  # Filter loci that pass the threshold
  significant <- lava_data %>%
    filter(p < bonferroni_threshold) %>%
    select(locus, chr, start, stop, phen1, phen2, p)  # Keep necessary columns
  
  # Append results to the combined data frame
  significant_loci <- bind_rows(significant_loci, significant)
}

# Save significant loci to the output file
if (nrow(significant_loci) > 0) {
  write_tsv(significant_loci, output_file)
  cat("Significant loci saved to:", output_file, "\n")
} else {
  cat("No significant loci found across all files.\n")
}


## Step 2: Find FUMA compatible genes from loci obtained with LAVA

library(GenomicRanges)
library(dplyr)
library(rtracklayer)
library(readr)
library(stringr)

# Define input and output
lava_file <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/FUMA/joined_lava_file/NP_PD_female.txt"  # Replace with the path to your LAVA results file
gtf_file <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/FUMA/gencode.v37lift37.annotation.gtf.gz"  # Replace with the path to your GTF file
output_file <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/FUMA/genes_for_FUMA/NP_PD_female.txt"  # Path to save FUMA-compatible gene list

# Load LAVA results
lava_data <- read_tsv(lava_file, col_types = cols())

# Ensure LAVA results have necessary columns
required_cols <- c("locus", "chr", "start", "stop")
if (!all(required_cols %in% colnames(lava_data))) {
  stop("The LAVA results file must contain the columns: locus, chr, start, stop.")
}

# Convert LAVA chromosomes to match GTF format if necessary
lava_data <- lava_data %>%
  mutate(chr = ifelse(str_detect(chr, "^chr"), chr, str_c("chr", chr)))

# Convert LAVA loci to GRanges
locus_gr <- GRanges(
  seqnames = lava_data$chr,
  ranges = IRanges(start = lava_data$start, end = lava_data$stop),
  locus = lava_data$locus
)

# Load reference gene annotation (GTF format)
ref <- rtracklayer::import(gtf_file)
ref <- ref[ref$type == "gene"]


# Extract overlapping genes for all loci
overlap <- GenomicRanges::findOverlaps(locus_gr, ref)
genes_by_locus <- data.frame(
  locus = lava_data$locus[queryHits(overlap)],
  gene_id = ref[subjectHits(overlap)]$gene_id,
  gene_name = ref[subjectHits(overlap)]$gene_name
) %>%
  distinct(gene_name)  # Ensure unique gene names

# Save genes in FUMA-compatible format
writeLines(genes_by_locus$gene_name, con = output_file)

cat("Gene list saved to:", output_file, "\n")
