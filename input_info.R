# Load necessary libraries
library(dplyr)
library(here)

# Set the project and summary statistics directories
project_dir <- "~/projects/def-gsarah/sadafgy/LAVA/info_files"
sumstats_dir <- "/scratch/sadafgy/processed_sumstats"

# Read the combined case/control information
N_samples <- read.table(
  file = here(project_dir, "case_control_input.csv"), 
  sep = ",", 
  header = TRUE
) %>%
  arrange(trait) %>%
  dplyr::rename(phenotype = trait, cases = n_cases, controls = n_controls)

# Get the list of summary statistic files in the sumstats directory
file_paths <- list.files(
  path = here(sumstats_dir),
  full.names = TRUE,
  recursive = TRUE,
  pattern = "\\.tsv$"  # Match only .tsv files
)

# Create a data frame for file paths and extract phenotype names from filenames
file_path_df <- data.frame(
  filename = file_paths,
  phenotype = basename(file_paths) %>%
    gsub("\\.tsv$", "", .)  # Remove the file extension to match phenotype names
)

# Merge the case/control information with the file paths by phenotype
input_info <- left_join(N_samples, file_path_df, by = "phenotype")
output_file <- file.path(project_dir, "input_info.txt")

# Save the combined input info file
write.table(
  input_info,
  file = output_file,
  sep = "\t", 
  row.names = FALSE, 
  quote = FALSE
)

# Output for confirmation
cat("Combined input info file created successfully at:", output_file)
