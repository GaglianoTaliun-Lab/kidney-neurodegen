# load packages
library(tidyverse)
library(stringr)
library(data.table)
library(readr)

# Replace this with your actual project directory
project_dir <- "/home/sadafgy/projects/def-gsarah/sadafgy"
gwas_dir <- file.path(project_dir, "ldscore_regression")
out_dir <- file.path(project_dir, "LAVA", "ldsc_corr")

# number of lines to skip in the ldsc output:
n_skip = 60
# number of lines to read from the output:
n_read = 1

###### Extracting LDSC results ######

# Extract all outputs into a single file
file_paths <- list.files(gwas_dir, pattern = ".log", full.names = TRUE)

files <- vector(mode = "list", length = length(file_paths))

for(i in seq_along(file_paths)) {
  files[[i]] <- read_table(file = file_paths[i], skip = n_skip, n_max = n_read)
}

# Clean the p1 and p2 columns to remove file extensions and unwanted characters
all_rg <- rbindlist(files, fill = TRUE) %>%
  mutate(
    p1 = basename(p1) %>% str_remove(".sumstats.gz"),
    p2 = basename(p2) %>% str_remove(".sumstats.gz")
  ) %>%
  select(
    p1, p2, rg, se, z, p, h2_obs, h2_obs_se, h2_int, h2_int_se, gcov_int, gcov_int_se
  )

# Save data 
write.table(all_rg, file = file.path(out_dir, "ldsc_correlations.txt"), sep = "\t", row.names = FALSE, quote = FALSE)



###### Creating sample overlap matrix by extracting the intercept from LDSC results ######

# Load required libraries
library(tidyverse)
library(data.table)

# Load the LDSC results
ldsc_results <- read.table("/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/ldsc_corr/ldsc_correlations.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract unique phenotypes for proxies and non-proxies
phenotypes_prox <- c("AD_female", "AD_male6", "AD_modified-v2", "PD_females", "PD_males", "PD_sexcomb")
phenotypes_noprox <- c("NP_PD_female", "NP_PD_male", "NP_PD_sexcomb", "NP_AD_sexcomb")

# Directories for saving outputs
output_dir <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/sample_overlap"

# Create a function to generate correlation matrix for given phenotype set and save in specific directory
generate_covar_matrix <- function(phenotypes_set, output_name, output_dir) {
  
  # Create an empty correlation matrix
  n <- length(phenotypes_set)
  covar_matrix <- matrix(NA, n, n, dimnames = list(phenotypes_set, phenotypes_set))
  
  # Fill the matrix with gcov_int values from ldsc_results
  for (i in phenotypes_set) {
    for (j in phenotypes_set) {
      value <- ldsc_results %>%
        filter((p1 == i & p2 == j) | (p1 == j & p2 == i)) %>%
        select(gcov_int) %>%
        pull()
      if (length(value) == 1 && !is.na(value)) {
        covar_matrix[i, j] <- value
      }
    }
  }
  
  # Fill missing diagonals with 1 (self-correlations)
  diag(covar_matrix) <- 1
  
  # Write the matrix to a CSV file in the specific output directory
  write.csv(covar_matrix, file.path(output_dir, paste0(output_name, "_sample_overlap.csv")), row.names = TRUE)
}

# Generate sample overlap matrices for each phenotype set
generate_covar_matrix(phenotypes_prox, "proxies", file.path(output_dir, "AD_output"))
generate_covar_matrix(phenotypes_noprox, "non_proxies", file.path(output_dir, "no_proxies_output"))

# Generate individual matrices for specific phenotype combinations and save in their respective directories
phenotype_combinations <- list(
  # AD combinations
  list(combo = c("AD_female", "creatinine_females"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_male6", "creatinine_male"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_modified-v2", "creatinine_sexcombined"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_female", "microalbumine_females"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_male6", "microalbumine_males"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_modified-v2", "microalbumin_sexcomb"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_female", "GFR_females"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_male6", "GFR_males"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_modified-v2", "GFR_sexcombined"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_female", "potassium_females"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_male6", "potassium_males"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_modified-v2", "potassium_sexcombined"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_female", "sodium_female"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_male6", "sodium_male"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_modified-v2", "sodium_sexcombined"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_female", "hematuria_females"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_male6", "hematuria_males"), dir = file.path(output_dir, "AD_output")),
  list(combo = c("AD_modified-v2", "hematuria_sexcombined"), dir = file.path(output_dir, "AD_output")),
  
  # PD combinations
  list(combo = c("PD_females", "creatinine_females"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_males", "creatinine_male"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_sexcomb", "creatinine_sexcombined"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_females", "microalbumine_females"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_males", "microalbumine_males"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_sexcomb", "microalbumin_sexcomb"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_females", "GFR_females"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_males", "GFR_males"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_sexcomb", "GFR_sexcombined"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_females", "potassium_females"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_males", "potassium_males"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_sexcomb", "potassium_sexcombined"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_females", "sodium_female"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_males", "sodium_male"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_sexcomb", "sodium_sexcombined"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_females", "hematuria_females"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_males", "hematuria_males"), dir = file.path(output_dir, "PD_output")),
  list(combo = c("PD_sexcomb", "hematuria_sexcombined"), dir = file.path(output_dir, "PD_output")),
  
  # NP_PD (no proxies) combinations
  list(combo = c("NP_PD_female", "creatinine_females"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_male", "creatinine_male"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_female", "microalbumine_females"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_male", "microalbumine_males"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_female", "GFR_females"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_male", "GFR_males"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_female", "potassium_females"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_male", "potassium_males"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_female", "sodium_female"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_male", "sodium_male"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_female", "hematuria_females"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_male", "hematuria_males"), dir = file.path(output_dir, "no_proxies_output")),
  
  # NP_PD_sexcomb and NP_AD_sexcomb
  list(combo = c("NP_PD_sexcomb", "creatinine_sexcombined"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_sexcomb", "microalbumine_sexcombined"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_sexcomb", "potassium_sexcombined"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_sexcomb", "sodium_sexcombined"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_sexcomb", "GFR_sexcombined"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_PD_sexcomb", "hematuria_sexcombined"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_AD_sexcomb", "creatinine_sexcombined"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_AD_sexcomb", "microalbumine_sexcombined"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_AD_sexcomb", "potassium_sexcombined"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_AD_sexcomb", "sodium_sexcombined"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_AD_sexcomb", "GFR_sexcombined"), dir = file.path(output_dir, "no_proxies_output")),
  list(combo = c("NP_AD_sexcomb", "hematuria_sexcombined"), dir = file.path(output_dir, "no_proxies_output"))
)

# Loop through the phenotype combinations and generate matrices
for (combo_info in phenotype_combinations) {
  generate_covar_matrix(combo_info$combo, paste0(combo_info$combo[1], "_vs_", combo_info$combo[2]), combo_info$dir)
}
