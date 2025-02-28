## Step 1 : Find the number of loci in the bivariate results that respects the p-value thrsehold of 0.05 and bonferroni thrsehold. 

# Load necessary library
library(dplyr)

# Directories for input and output
bivar_files <- c(
  "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/lava_results/bivariate_results/PD_filtered"
)
output_less_0_05 <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/phenoGram_plot/passing_bivariate_0.05"
output_bonferroni <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/phenoGram_plot/bivariate_bonferroni_threshold"

# Thresholds
threshold_0_05 <- 0.05

# Function to process individual files
process_file <- function(file_path, bonferroni_threshold, output_dir_005, output_dir_bonferroni) {
  # Read the bivariate results
  bivariate_results <- read.delim(file_path, header = TRUE, sep = "\t")
  
  # Filter loci with p-values less than 0.05
  loci_less_0_05 <- bivariate_results %>%
    filter(p < threshold_0_05)
  
  # Filter for significant results for Bonferroni threshold
  significant_results <- bivariate_results %>%
    filter(p < bonferroni_threshold) %>%
    arrange(p) # Sort by p-value
  
  # Get filename
  file_name <- basename(file_path)
  
  # Write results to output directories
  write.table(loci_less_0_05, 
              file = file.path(output_dir_005, file_name), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  write.table(significant_results, 
              file = file.path(output_dir_bonferroni, file_name), 
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Process all files
for (bivar_dir in bivar_files) {
  file_list <- list.files(bivar_dir, pattern = "*.lava.tsv", full.names = TRUE)
  
  for (file_path in file_list) {
    # Read the file to calculate Bonferroni threshold
    bivariate_results <- read.delim(file_path, header = TRUE, sep = "\t")
    num_tests <- nrow(bivariate_results)
    bonferroni_threshold <- threshold_0_05 / num_tests
    
    # Call the process_file function
    process_file(file_path, bonferroni_threshold, output_less_0_05, output_bonferroni)
  }
}

cat("Processing completed. Results are saved in the specified directories.\n")



## Step 2 : Combine male, female, sex-combined and sex-combined meta-analyzed loci (with and without proxies) into one file (input for phenogram plot)

# Load necessary library
library(dplyr)

# Read the input file
input_file <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/phenoGram_plot/passing_bivariate_0.05"  #Change directory to passing bonferroni or passing univariate test 

proxy_file_female <- file.path(input_file, "PD", "PD_females_GFR_females_bivar_combined.lava.tsv")  # Replace with the actual file path
no_proxy_file_female <- file.path(input_file, "NP_PD", "NP_PD_female_GFR_females_bivar_combined.lava.tsv")
proxy_file_male <- file.path(input_file, "PD", "PD_males_GFR_males_bivar_combined.lava.tsv") 
no_proxy_file_male <- file.path(input_file, "NP_PD", "NP_PD_male_GFR_males_bivar_combined.lava.tsv")
proxy_file_comb <- file.path(input_file, "PD_filtered", "PD_sexcomb_GFR_sexcombined_bivar_combined.lava.tsv") 
no_proxy_file_comb <- file.path(input_file, "NP_PD", "NP_PD_sexcomb_GFR_sexcombined_bivar_combined.lava.tsv")
proxy_file_meta <- file.path(input_file, "PD_meta", "PD_sexcombined_meta_GFR_sexcombined_bivar_combined.lava.tsv") 
no_proxy_file_meta <- file.path(input_file, "NP_PD", "NP_PD_sexcombined_meta_GFR_sexcombined_bivar_combined.lava.tsv")

proxy_female <- read.delim(proxy_file_female, header = TRUE)
no_proxy_female <- read.delim(no_proxy_file_female, header = TRUE)
proxy_male <- read.delim(proxy_file_male, header = TRUE)
no_proxy_male <- read.delim(no_proxy_file_male, header = TRUE)
proxy_comb <- read.delim(proxy_file_comb, header = TRUE)
no_proxy_comb <- read.delim(no_proxy_file_comb, header = TRUE)
proxy_meta <- read.delim(proxy_file_meta, header = TRUE)
no_proxy_meta <- read.delim(no_proxy_file_meta, header = TRUE)

create_phenotype_data <- function(data, phenotype_name) {
  data %>%
    select(locus, chr, start) %>%       # Select required columns
    rename(pos = start) %>%            # Rename 'start' to 'pos'
    mutate(phenotype = phenotype_name) # Add the phenotype column
}

proxy_female <- create_phenotype_data(proxy_female, "Female Proxy")
no_proxy_female <- create_phenotype_data(no_proxy_female, "Female No Proxy")
proxy_male <- create_phenotype_data(proxy_male, "Male Proxy")
no_proxy_male <- create_phenotype_data(no_proxy_male, "Male No proxy")
proxy_comb <- create_phenotype_data(proxy_comb, "Sex-combined Proxy")
no_proxy_comb <- create_phenotype_data(no_proxy_comb, "Sex-combined No proxy")
proxy_meta <- create_phenotype_data(proxy_meta, "Sex-combined meta-analyzed with Proxy")
no_proxy_meta <- create_phenotype_data(no_proxy_meta, "Sex-combined meta-analyzed No proxy")

combined_data <- bind_rows(proxy_female, no_proxy_female, proxy_male,no_proxy_male, proxy_comb, no_proxy_comb, proxy_meta, no_proxy_meta)

# Ensure only the desired columns are kept
combined_data <- combined_data %>%
  select(locus, chr, pos, phenotype)  # Retain only the 4 specified columns

# Specify the output directory
output_file <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/phenoGram_plot/ phenogram_data_005_thrsehold" 

# Concatenate the directory path with the desired file name
output_path <- file.path(output_file, "PD_GFR.tsv")

# Write the data to the specified file
write.table(combined_data, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)

