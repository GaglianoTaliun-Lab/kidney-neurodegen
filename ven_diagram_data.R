## 1. Number of loci passing the univariate test with a Bonferroni thrsehold
# Load necessary libraries
library(dplyr)

# Define file paths and output directory
file1 <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/lava_results/bivariate_results/AD/AD_sexcomb_hematuria_sexcombined_bivar_combined.lava.tsv"  # Replace with the path to the first file
file2 <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/lava_results/bivariate_results/NP_AD/NP_AD_sexcomb_hematuria_sexcombined_bivar_combined.lava.tsv"  # Replace with the path to the second file
output_dir <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/ven_diagram"  # Replace with the path to the output directory

# Define function to compare loci between two files
compare_loci <- function(file1, file2, output_dir) {
  
  # Load the two datasets
  data1 <- read.table(file1, header = TRUE)
  data2 <- read.table(file2, header = TRUE)
  
  # Initialize data frames to store results
  similar <- data.frame()
  unique_file1 <- data.frame()
  unique_file2 <- data.frame()
  
  # Compare loci from file1 against file2
  for (i in seq_len(nrow(data1))) {
    locus1 <- data1[i, ]
    match <- FALSE
    
    for (j in seq_len(nrow(data2))) {
      locus2 <- data2[j, ]
      
      # Check if loci are the same (based on chr, start, and stop)
      if (locus1$chr == locus2$chr &&
          locus1$start == locus2$start &&
          locus1$stop == locus2$stop) {
        # Add to the similar list
        similar <- bind_rows(similar, locus1)
        match <- TRUE
        break
      }
    }
    
    # If no match found, add to unique_file1
    if (!match) {
      unique_file1 <- bind_rows(unique_file1, locus1)
    }
  }
  
  # Identify loci in file2 that are not in the similar list
  for (k in seq_len(nrow(data2))) {
    locus2 <- data2[k, ]
    match <- FALSE
    
    for (m in seq_len(nrow(similar))) {
      similar_locus <- similar[m, ]
      
      if (locus2$chr == similar_locus$chr &&
          locus2$start == similar_locus$start &&
          locus2$stop == similar_locus$stop) {
        match <- TRUE
        break
      }
    }
    
    # If no match found, add to unique_file2
    if (!match) {
      unique_file2 <- bind_rows(unique_file2, locus2)
    }
  }
  
  # Save results to files
  write.table(similar, file = file.path(output_dir, "similar_loci/AD_hematuria_sexcomb.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(unique_file1, file = file.path(output_dir, "with_proxy/unique_AD_hematuria_sexcomb.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(unique_file2, file = file.path(output_dir, "no_proxy/unique_NP_AD_hematuria_sexcomb.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Print summary
  cat("Number of similar loci:", nrow(similar), "\n")
  cat("Number of unique loci in file1:", nrow(unique_file1), "\n")
  cat("Number of unique loci in file2:", nrow(unique_file2), "\n")
}

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Run the function
compare_loci(file1, file2, output_dir)


## 2. Number of loci passing the bivariate test after applying a p-value thrsehold of p<0.05

# Load necessary libraries 
library(dplyr)

# Define file paths and output directory
file1 <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/lava_results/bivariate_results/AD/AD_sexcombined_microalbumine_sexcombined_bivar_combined.lava.tsv"  # Replace with the path to the first file
file2 <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/lava_results/bivariate_results/NP_AD/NP_AD_sexcomb_microalbumine_sexcombined_bivar_combined.lava.tsv"  # Replace with the path to the second file
output_dir <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/ven_diagram/p_value_0.05"  # Replace with the path to the output directory

# Define function to compare loci between two files
compare_loci <- function(file1, file2, output_dir, p_threshold=0.05) {
  
  # Load the two datasets
  data1 <- read.table(file1, header = TRUE, sep="\t")
  data2 <- read.table(file2, header = TRUE, sep="\t")
  
    # Filter datasets based on p-value threshold
  data1 <- data1 %>% filter(p <= p_threshold)
  data2 <- data2 %>% filter(p <= p_threshold)

  # Initialize data frames to store results
  similar <- data.frame()
  unique_file1 <- data.frame()
  unique_file2 <- data.frame()
  
  # Compare loci from file1 against file2
  for (i in seq_len(nrow(data1))) {
    locus1 <- data1[i, ]
    match <- FALSE
    
    for (j in seq_len(nrow(data2))) {
      locus2 <- data2[j, ]
      
      # Check if loci are the same (based on chr, start, and stop)
      if (locus1$chr == locus2$chr &&
          locus1$start == locus2$start &&
          locus1$stop == locus2$stop) {
        # Add to the similar list
        similar <- bind_rows(similar, locus1)
        match <- TRUE
        break
      }
    }
    
    # If no match found, add to unique_file1
    if (!match) {
      unique_file1 <- bind_rows(unique_file1, locus1)
    }
  }
  
  # Identify loci in file2 that are not in the similar list
  for (k in seq_len(nrow(data2))) {
    locus2 <- data2[k, ]
    match <- FALSE
    
    for (m in seq_len(nrow(similar))) {
      similar_locus <- similar[m, ]
      
      if (locus2$chr == similar_locus$chr &&
          locus2$start == similar_locus$start &&
          locus2$stop == similar_locus$stop) {
        match <- TRUE
        break
      }
    }
    
    # If no match found, add to unique_file2
    if (!match) {
      unique_file2 <- bind_rows(unique_file2, locus2)
    }
  }
  
  # Save results to files
  write.table(similar, file = file.path(output_dir, "similar_loci/AD_hematuria_sexcomb.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(unique_file1, file = file.path(output_dir, "with_proxy/unique_AD_hematuria_sexcomb.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(unique_file2, file = file.path(output_dir, "no_proxy/unique_NP_AD_hematuria_sexcomb.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Print summary
  cat("Number of similar loci:", nrow(similar), "\n")
  cat("Number of unique loci in file1:", nrow(unique_file1), "\n")
  cat("Number of unique loci in file2:", nrow(unique_file2), "\n")
}

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Run the function
compare_loci(file1, file2, output_dir, p_threshold=0.05)



## 3. Number of loci passing the bivariate test after applying a Bonferroni threshold (0.05/N_loci tested)

# Load necessary libraries
library(dplyr)

# Define file paths and output directory
file1 <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/lava_results/bivariate_results/AD/AD_sexcombined_potassium_sexcombined_bivar_combined.lava.tsv"
file2 <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/lava_results/bivariate_results/NP_AD/NP_AD_sexcomb_potassium_sexcombined_bivar_combined.lava.tsv"
output_dir <- "/home/sadafgy/projects/def-gsarah/sadafgy/LAVA/ven_diagram/p_value_nloci"

# Define function to compare loci between two files
compare_loci <- function(file1, file2, output_dir) {
  
  # Load the two datasets
  data1 <- read.table(file1, header = TRUE, sep = "\t")
  data2 <- read.table(file2, header = TRUE, sep = "\t")
  
  n_loci_file1 <- nrow(data1)
  n_loci_file2 <- nrow(data2)

  # Calculate the bivariate threshold: 0.05 / n_loci
  bivar_threshold1 <- 0.05 / n_loci_file1
  bivar_threshold2 <- 0.05 / n_loci_file2
  
  cat("Number of loci passing univariate threshold in file1 (PD):", n_loci_file1, "\n")
  cat("Number of loci passing univariate threshold in file2 (NP_PD):", n_loci_file2, "\n")
  cat("Bivariate threshold for file1:", bivar_threshold1, "\n")
  cat("Bivariate threshold for file2:", bivar_threshold2, "\n")
  
  # Filter loci based on the bivariate thresholds
  data1_filtered <- data1 %>% filter(p <= bivar_threshold1)
  data2_filtered <- data2 %>% filter(p <= bivar_threshold2)
  
  cat("Number of loci passing bivariate threshold in file1:", nrow(data1_filtered), "\n")
  cat("Number of loci passing bivariate threshold in file2:", nrow(data2_filtered), "\n")
  
  # Initialize data frames to store results
  similar <- data.frame()
  unique_file1 <- data.frame()
  unique_file2 <- data.frame()
  
  # Compare loci from file1 against file2
  for (i in seq_len(nrow(data1_filtered))) {
    locus1 <- data1_filtered[i, ]
    match <- FALSE
    
    for (j in seq_len(nrow(data2_filtered))) {
      locus2 <- data2_filtered[j, ]
      
      # Check if loci are the same (based on chr, start, and stop)
      if (locus1$chr == locus2$chr &&
          locus1$start == locus2$start &&
          locus1$stop == locus2$stop) {
        # Add to the similar list
        similar <- bind_rows(similar, locus1)
        match <- TRUE
        break
      }
    }
    
    # If no match found, add to unique_file1
    if (!match) {
      unique_file1 <- bind_rows(unique_file1, locus1)
    }
  }
  
  # Identify loci in file2 that are not in the similar list
  for (k in seq_len(nrow(data2_filtered))) {
    locus2 <- data2_filtered[k, ]
    match <- FALSE
    
    for (m in seq_len(nrow(similar))) {
      similar_locus <- similar[m, ]
      
      if (locus2$chr == similar_locus$chr &&
          locus2$start == similar_locus$start &&
          locus2$stop == similar_locus$stop) {
        match <- TRUE
        break
      }
    }
    
    # If no match found, add to unique_file2
    if (!match) {
      unique_file2 <- bind_rows(unique_file2, locus2)
    }
  }
 
  
  # Save results to files
  write.table(similar, file = file.path(output_dir, "similar_loci/similar_AD_hematuria_sexcomb.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(unique_file1, file = file.path(output_dir, "with_proxy/unique_AD_hematuria_sexcomb.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(unique_file2, file = file.path(output_dir, "no_proxy/unique_NP_AD_hematuria_sexcomb.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Print summary
  cat("Number of similar loci:", nrow(similar), "\n")
  cat("Number of unique loci in file1:", nrow(unique_file1), "\n")
  cat("Number of unique loci in file2:", nrow(unique_file2), "\n")
}

# Run the function
compare_loci(file1, file2, output_dir)
