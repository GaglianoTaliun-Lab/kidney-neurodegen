# Load necessary libraries
library(tidyverse)
library(data.table)

# Read the VCF file
vcf_data <- fread(file_path, skip = "#CHROM", header = TRUE)
head(vcf_data)

# Extract the required fields from the FORMAT column, and split ES, SE, LP, AF, ID
vcf_processed <- vcf_data %>%
  separate('UKB-b-323', into = c("ES", "SE", "LP", "AF", "ID"), sep = ":") %>%
  mutate(ES = as.numeric(ES),
         SE = as.numeric(SE),
         LP = as.numeric(LP),
         AF = as.numeric(AF),
         Z = ES / SE,                    # Calculate Z-score from effect size and standard error
         P = 10^(-LP),                   # Convert LP (log p-value) back to p-value
         MAF = pmin(AF, 1 - AF)) %>%     # Calculate MAF from AF
  select('#CHROM', POS, ID, REF, ALT, SE, Z, P, MAF)

# Rename the columns to match standard GWAS summary stats format
vcf_final <- vcf_processed %>%
  rename(SNP = ID,
         CHR = '#CHROM',
         BP = POS,
         a1 = ALT,  # Reference allele
         a2 = REF,  # Alternate allele
         pval = P)  

# Filter out rows where values are missing or invalid
vcf_final <- vcf_final %>%
  filter(!is.na(Z), !is.na(pval), !is.na(SNP))  # Ensure Z, pval, and SNP are not missing

# Write the processed summary statistics to a file
write.table(vcf_final, "/path_to_output/processed_sumstats/processed_sumstats.txt", sep = "\t", row.names = FALSE, quote = FALSE)





