# Limit OpenMP threads to 1 
Sys.setenv(OMP_NUM_THREADS = 1)

# Limit the number of threads in data.table, if it’s being used
data.table::setDTthreads(1)

# Set number of threads for RcppParallel (used by packages like tidyverse or LAVA)
RcppParallel::setThreadOptions(numThreads = 1)

# Load necessary packages
library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Define base project directory
project_dir <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA"

# Define trait pairs
all_traits <- list(
  c("AD_male6", "microalbumine_males"),
  c("AD_female", "microalbumine_females"),
  c("AD_modified-v2", "microalbumine_sexcombined"),
  c("AD_male6", "creatinine_male"),
  c("AD_female", "creatinine_female"),
  c("AD_modified-v2", "creatinine_sexcombined"),
  c("AD_male6", "potassium_males"),
  c("AD_female", "potassium_females"),
  c("AD_modified-v2", "potassium_sexcombined"),
  c("AD_male6", "sodium_male"),
  c("AD_female", "sodium_female"),
  c("AD_modified-v2", "sodium_sexcombined"),
  c("AD_male6", "hematuria_males"),
  c("AD_female", "hematuria_females"),
  c("AD_modified-v2", "hematuria_sexcombined"),
  c("AD_male6", "GFR_males"),
  c("AD_female", "GFR_females"),
  c("AD_modified-v2", "GFR_sexcombined"),
  c("PD_males", "microalbumine_males"),
  c("PD_females", "microalbumine_females"),
  c("PD_sexcomb", "microalbumine_sexcombined"),
  c("PD_males", "creatinine_male"),
  c("PD_females", "creatinine_females"),
  c("PD_sexcomb", "creatinine_sexcombined"),
  c("PD_males", "potassium_males"),
  c("PD_females", "potassium_females"),
  c("PD_sexcomb", "potassium_sexcombined"),
  c("PD_males", "sodium_male"),
  c("PD_females", "sodium_female"),
  c("PD_sexcomb", "sodium_sexcombined"),
  c("PD_males", "hematuria_males"),
  c("PD_females", "hematuria_females"),
  c("PD_sexcomb", "hematuria_sexcombined"),
  c("PD_males", "GFR_males"),
  c("PD_females", "GFR_females"),
  c("PD_sexcomb", "GFR_sexcombined"),
  c("NP_PD_male", "microalbumine_males"),
  c("NP_PD_female", "microalbumine_females"),
  c("NP_PD_male", "creatinine_male"),
  c("NP_PD_female", "creatinine_females"),
  c("NP_PD_male", "potassium_males"),
  c("NP_PD_female", "potassium_females"),
  c("NP_PD_male", "sodium_male"),
  c("NP_PD_female", "sodium_female"),
  c("NP_PD_male", "hematuria_males"),
  c("NP_PD_female", "hematuria_females"),
  c("NP_PD_male", "GFR_males"),
  c("NP_PD_female", "GFR_females"),
  c("NP_PD_sexcomb", "microalbumine_sexcombined"),
  c("NP_PD_sexcomb", "creatinine_sexcombined"),
  c("NP_PD_sexcomb", "potassium_sexcombined"),
  c("NP_PD_sexcomb", "sodium_sexcombined"),
  c("NP_PD_sexcomb", "hematuria_sexcombined"),
  c("NP_PD_sexcomb", "GFR_sexcombined"),
  c("NP_AD_sexcomb", "microalbumine_sexcombined"),
  c("NP_AD_sexcomb", "creatinine_sexcombined"),
  c("NP_AD_sexcomb", "potassium_sexcombined"),
  c("NP_AD_sexcomb", "sodium_sexcombined"),
  c("NP_AD_sexcomb", "hematuria_sexcombined"),
  c("NP_AD_sexcomb", "GFR_sexcombined"))

# Select which trait pair to process
selected_index <- 3  # Change this to process a different trait pair
traits <- all_traits[[selected_index]]

# Determine the directory for sample overlap files based on trait names
get_overlap_directory <- function(trait) {
  if (str_detect(trait, "^AD")) {
    return("AD_output")
  } else if (str_detect(trait, "^PD")) {
    return("no_proxies_output")
  } else {
    return("no_proxies_output")  # Default case
  }
}

# Construct the sample overlap file path
overlap_dir <- get_overlap_directory(traits[1])
sample_overlap_file <- file.path(project_dir, "sample_overlap", overlap_dir, paste0(traits[1], "_vs_", traits[2], "_sample_overlap.txt"))

# Set up arguments
args <- list(
  ref_prefix = file.path(project_dir, "g1000_eur/g1000_eur"),
  loc_file = file.path(project_dir, "blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"),
  info_file = file.path(project_dir, "info_files/input_info.txt"), 
  sample_overlap_file = sample_overlap_file,
  phenotypes = traits,
  output_filename = paste(traits, collapse = "_")
)

# Load the complete input_info file and filter it for the selected traits
input_info_data <- read.table(args$info_file, sep = "\t", header = TRUE)
filtered_info <- input_info_data %>% 
  filter(phenotype %in% args$phenotypes) %>%
  rename(cases = cases, controls = controls)

# Check if both traits are present in the filtered data
if (nrow(filtered_info) < 2) {
      stop("Both selected traits were not found in the input_info file.")
}

# Prepare temporary filtered input file for LAVA
temp_info_file <- tempfile()
write.table(filtered_info, temp_info_file, sep = "\t", row.names = FALSE, quote = FALSE)

# Load loci and input data for selected traits
loci <- LAVA::read.loci(args$loc_file)
n_loci <- nrow(loci)
input <- LAVA::process.input(
  input.info.file = temp_info_file,
  sample.overlap.file = args$sample_overlap_file,
  ref.prefix = args$ref_prefix,
  phenos = args$phenotypes
)

# Parameters for block processing
bloc_size <- 100   # Number of loci per block
n_blocks <- ceiling(n_loci / bloc_size)
univar_threshold <- 0.05/n_loci

# Define output directories for each trait pair and block results
output_block_dir <- file.path(project_dir, "lava_results", "partial", args$output_filename)
combined_univar_dir <- file.path(project_dir, "lava_results", "univar_results", "PD")
combined_bivar_dir <- file.path(project_dir, "lava_results", "bivar_results", "PD")
dir.create(output_block_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(combined_bivar_dir, showWarnings = FALSE, recursive = TRUE)


# Initialize lists to collect all results across blocks
all_univar_results <- list()
all_bivar_results <- list()
successful_loci <- 0  # Counter for successful loci

# Process loci in blocks
for (bloc_num in 1:n_blocks) {
  bloc_start <- (bloc_num - 1) * bloc_size + 1
  bloc_end <- min(bloc_num * bloc_size, n_loci)
  loci_bloc <- loci[bloc_start:bloc_end, ]
  
  print(paste("Processing block", bloc_num, "of", n_blocks, "(", bloc_start, "-", bloc_end, ")"))
    
  # Initialize lists to collect block-specific results
  univar_bloc <- list()
  bivar_bloc <- list()

  for (i in seq_len(nrow(loci_bloc))) {
    locus <- LAVA::process.locus(loci_bloc[i, ], input)
    
    if (!is.null(locus)) {
      successful_loci <- successful_loci + 1
      loc_info <- data.frame(
        locus = locus$id, chr = locus$chr,
        start = locus$start, stop = locus$stop,
        n_snps = locus$n.snps, n_pcs = locus$K
      )
      
      # Run univariate test
      univ_out <- LAVA::run.univ(locus)
      
      # Store univariate result
      univar_bloc[[i]] <- dplyr::bind_cols(loc_info, univ_out)
      
      # Run bivariate test if univariate results pass threshold
      if (all(univ_out$p < univar_threshold, na.rm = TRUE)) {
        bivar_out <- LAVA::run.bivar(locus)
        
        # Store bivariate results if they exist
        if (!is.null(bivar_out)) {
          bivar_bloc[[i]] <- dplyr::bind_cols(loc_info, bivar_out)
        }
      }
    }
  }
  
  # Combine results for the current block
  univar_bloc_df <- bind_rows(univar_bloc)
  bivar_bloc_df <- bind_rows(bivar_bloc)
    
  # Save the block-specific results
  write.table(
    univar_bloc_df,
    file = file.path(output_block_dir, paste0(args$output_filename, "_univ_block_", bloc_num, ".lava.tsv")),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  
  if (nrow(bivar_bloc_df) > 0) {
    write.table(
      bivar_bloc_df,
      file = file.path(output_block_dir, paste0(args$output_filename, "_bivar_block_", bloc_num, ".lava.tsv")),
      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
    )
  }

  # Append block results to cumulative lists
  all_univar_results[[bloc_num]] <- univar_bloc_df
  all_bivar_results[[bloc_num]] <- bivar_bloc_df
}

# Combine all blocks' results into final data frames
final_univar_df <- bind_rows(all_univar_results)
final_bivar_df <- bind_rows(all_bivar_results)

# Save the final combined results to output files
write.table(
  final_univar_df,
  file = file.path(combined_univar_dir, paste0(args$output_filename, "_univ_combined.lava.tsv")),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
)

if (nrow(final_bivar_df) > 0) {
  write.table(
    final_bivar_df,
    file = file.path(combined_bivar_dir, paste0(args$output_filename, "_bivar_combined.lava.tsv")),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
}

# Final summary
print(paste("Analysis complete! Final results saved in", output_block_dir))
print(paste("Total number of loci successfully tested:", successful_loci))

# Remove temporary file
unlink(temp_info_file)
