# libraries needed to use the sample script
# https://github.com/RHReynolds/colochelpR
library(BiocManager)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37) # use this library if sumstats are in GRCh37
library(colochelpR)
library(tidyverse)
library(data.table)

# rsids positions in GRCh37
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37

# path to summary statistics in your directory
file_path <- "/path_to_GWASsummstats/" 

# Load only the necessary columns (variant, n_complete_samples, beta, se, pval, minor_AF, low_confidence_variant)
sumstats <- fread(file_path, select = c("variant", "n_complete_samples", "beta", "se", "pval", "minor_AF", "low_confidence_variant"))
head(sumstats)

# Filter out low-confidence variants and keep relevant columns
sumstats_filtered <- sumstats %>%
  filter(low_confidence_variant == FALSE) %>%  # Remove low-confidence variants
  select(variant, n_complete_samples, beta, pval, se, minor_AF)  # Keep only essential columns

# Separate the 'variant' column into 'CHR', 'BP', 'a1', and 'a2'
sumstats_filtered <- sumstats_filtered %>%
  tidyr::separate(variant, into = c("CHR", "BP", "a1", "a2"), sep = ":") %>%
  dplyr::mutate(CHR = as.integer(CHR), BP = as.integer(BP))

# Filter out rows where CHR or BP is missing or invalid
sumstats_filtered <- sumstats_filtered %>%
  dplyr::filter(!is.na(CHR), !is.na(BP))

# need to have the columns 'CHR' and 'BP' to be able to annotate rsids.
# Convert CHR and BP to rsID using colochelpR's convert_loc_to_rs
sumstats_with_rs <- colochelpR::convert_loc_to_rs(df = sumstats_filtered, dbsnp_144)

# then need to remove positions with more than one rsid.
sumstats <- sumstats_with_rs %>%
  dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
  dplyr::group_by(CHR_BP) %>%
  dplyr::filter(!any(row_number() > 1))  %>%
  dplyr::ungroup() %>%
  dplyr::select(-CHR_BP)

# Write the processed summary statistics to a file
write.table(sumstats, "/path_to_output/processed_sumstats.txt", sep = "\t", row.names = FALSE, quote = FALSE)


