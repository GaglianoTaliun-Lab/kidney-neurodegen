#!/usr/bin/env Rscript
###############################################################################
# run_coloc_conjfdr_loci.R
#
# Runs Bayesian colocalization (coloc.abf) on conjFDR loci (conjFDR < 0.05)
# for all PD-kidney pairs with signal.
#
# Input:
#   - conjFDR supp tables directory (with annotated loci TSVs)
#   - Summary statistics directory
#   - Output directory
#
# Usage:
#   Rscript run_coloc_conjfdr_loci.R <supp_tables_dir> <sumstats_dir> <output_dir>
###############################################################################

suppressPackageStartupMessages({
  library(coloc)
  library(dplyr)
})

args <- commandArgs(TRUE)
if (length(args) < 3) stop("Usage: Rscript run_coloc_conjfdr_loci.R <supp_tables_dir> <sumstats_dir> <output_dir>")

supp_dir    <- args[1]
sumstats_dir <- args[2]
output_dir  <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("Colocalization analysis for conjFDR loci\n")
cat("============================================\n\n")

# ============================================================================
# Trait file mapping
# ============================================================================
trait_file_map <- c(
  "PD (meta)"            = "PD_metaSS_noProxy.tsv",
  "PD (female)"          = "PD_female_noProxy.tsv",
  "PD (male)"            = "PD_male_noProxy.tsv",
  "eGFR (meta)"          = "eGFR_meta.tsv",
  "eGFR (female)"        = "eGFR_female.tsv",
  "eGFR (male)"          = "eGFR_male.tsv",
  "uACR"                 = "uACR_sexComb.tsv",
  "uACR (female)"        = "uACR_female.tsv",
  "uACR (male)"          = "uACR_male.tsv",
  "Hematuria"            = "hematuria_sexComb.tsv",
  "Hematuria (female)"   = "hematuria_female.tsv",
  "Hematuria (male)"     = "hematuria_male.tsv",
  "uK/Cr"                = "uKCr_sexComb.tsv",
  "uK/Cr (female)"       = "uKCr_female.tsv",
  "uK/Cr (male)"         = "uKCr_male.tsv",
  "uNa/Cr"               = "uNaCr_sexComb.tsv",
  "uNa/Cr (female)"      = "uNaCr_female.tsv",
  "uNa/Cr (male)"        = "uNaCr_male.tsv"
)

trait_type_map <- c(
  "PD (meta)" = "cc", "PD (female)" = "cc", "PD (male)" = "cc",
  "eGFR (meta)" = "quant", "eGFR (female)" = "quant", "eGFR (male)" = "quant",
  "uACR" = "quant", "uACR (female)" = "quant", "uACR (male)" = "quant",
  "Hematuria" = "cc", "Hematuria (female)" = "cc", "Hematuria (male)" = "cc",
  "uK/Cr" = "quant", "uK/Cr (female)" = "quant", "uK/Cr (male)" = "quant",
  "uNa/Cr" = "quant", "uNa/Cr (female)" = "quant", "uNa/Cr (male)" = "quant"
)

# Case proportions for case-control traits
case_props <- c(
  "PD (meta)" = 20967 / 60633,
  "PD (female)" = 7947 / 28277,
  "PD (male)" = 13020 / 32356,
  "Hematuria" = 16409 / 396345,
  "Hematuria (female)" = 7061 / 220094,
  "Hematuria (male)" = 9774 / 187440
)

# ============================================================================
# Read master loci table
# ============================================================================
cat("Reading conjFDR loci...\n")
all_loci <- read.delim(file.path(supp_dir, "SuppTable_conjfdr_all_loci.tsv"),
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# After reading all_loci, add:
all_loci <- all_loci %>%
  filter(!(chr == 6 & position >= 25000000 & position <= 34000000))
                          
cat("  Total loci to colocalize:", nrow(all_loci), "\n\n")

# ============================================================================
# Define coloc window: +/- 500kb around each lead SNP
# ============================================================================
WINDOW <- 500000  # 500 kb on each side

# ============================================================================
# Helper: read and subset sumstats for a window
# ============================================================================
read_window_sumstats <- function(trait_label, chr_num, center_pos) {
  fname <- trait_file_map[trait_label]
  if (is.na(fname)) stop(paste("No file mapping for trait:", trait_label))
  
  fpath <- file.path(sumstats_dir, fname)
  if (!file.exists(fpath)) stop(paste("File not found:", fpath))
  
  ss <- read.delim(fpath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  ss_window <- ss %>%
    filter(chromosome == chr_num,
           base_pair_location >= (center_pos - WINDOW),
           base_pair_location <= (center_pos + WINDOW)) %>%
    filter(!is.na(beta), !is.na(standard_error), !is.na(p_value),
           standard_error > 0) %>%
    distinct(base_pair_location, .keep_all = TRUE)
  
  return(ss_window)
}

# ============================================================================
# Run coloc for each locus
# ============================================================================
results <- list()
total <- nrow(all_loci)
completed <- 0
failed <- 0

for (i in seq_len(total)) {
  row <- all_loci[i, ]
  
  cat(sprintf("[%d/%d] %s x %s | chr%s:%s | gene: %s\n",
              i, total, row$pd_stratum, row$kidney_trait,
              row$chr, row$position, row$nearest_gene))
  
  tryCatch({
    ss1 <- read_window_sumstats(row$pd_stratum, row$chr, row$position)
    ss2 <- read_window_sumstats(row$kidney_trait, row$chr, row$position)
    
    # Find shared SNPs by position
    shared_pos <- intersect(ss1$base_pair_location, ss2$base_pair_location)
    
    if (length(shared_pos) < 50) {
      cat(sprintf("  SKIP: too few shared SNPs (%d)\n", length(shared_pos)))
      failed <- failed + 1
      next
    }
    
    ss1 <- ss1 %>% filter(base_pair_location %in% shared_pos) %>% arrange(base_pair_location)
    ss2 <- ss2 %>% filter(base_pair_location %in% shared_pos) %>% arrange(base_pair_location)
    
    # Fix missing SNP IDs: use chr:pos as identifier
    ss1$snp_id <- ifelse(is.na(ss1$variant_id) | ss1$variant_id == "" | ss1$variant_id == ".",
                         paste0(row$chr, ":", ss1$base_pair_location),
                         ss1$variant_id)
    ss2$snp_id <- ifelse(is.na(ss2$variant_id) | ss2$variant_id == "" | ss2$variant_id == ".",
                         paste0(row$chr, ":", ss2$base_pair_location),
                         ss2$variant_id)
    
    # Use position-based IDs to ensure no NAs
    snp_ids <- paste0(row$chr, ":", ss1$base_pair_location)
    
    n1 <- median(ss1$N, na.rm = TRUE)
    n2 <- median(ss2$N, na.rm = TRUE)
    
    type1 <- trait_type_map[row$pd_stratum]
    type2 <- trait_type_map[row$kidney_trait]
    
    dataset1 <- list(
      beta = ss1$beta,
      varbeta = ss1$standard_error^2,
      snp = snp_ids,
      position = ss1$base_pair_location,
      type = type1,
      N = n1
    )
    
    dataset2 <- list(
      beta = ss2$beta,
      varbeta = ss2$standard_error^2,
      snp = snp_ids,
      position = ss2$base_pair_location,
      type = type2,
      N = n2
    )
    
    # For quantitative traits, coloc needs sdY or MAF
    # Estimate MAF from beta and SE: MAF ≈ 1/(2*N*SE^2) for standardized betas
    # Use a simple approximation: set sdY=1 (assumes phenotype is standardized)
    if (type1 == "quant") {
      dataset1$sdY <- 1
    }
    if (type2 == "quant") {
      dataset2$sdY <- 1
    }
    
    # Add case proportion for case-control traits
    if (type1 == "cc" && row$pd_stratum %in% names(case_props)) {
      dataset1$s <- case_props[row$pd_stratum]
    }
    if (type2 == "cc" && row$kidney_trait %in% names(case_props)) {
      dataset2$s <- case_props[row$kidney_trait]
    }
    
    # Run coloc.abf
    res <- suppressMessages(coloc.abf(dataset1, dataset2))
    
    cat(sprintf("  SNPs: %d | PP.H3=%.3f | PP.H4=%.3f\n",
                length(shared_pos),
                res$summary["PP.H3.abf"],
                res$summary["PP.H4.abf"]))
    
    results[[length(results) + 1]] <- data.frame(
      pair = paste0(row$pd_stratum, " x ", row$kidney_trait),
      pd_stratum = row$pd_stratum,
      kidney_trait = row$kidney_trait,
      locus_num = row$locus_num,
      chr = row$chr,
      position = row$position,
      conjfdr = row$conjfdr,
      nearest_gene = row$nearest_gene,
      n_shared_snps = length(shared_pos),
      PP.H0 = round(res$summary["PP.H0.abf"], 4),
      PP.H1 = round(res$summary["PP.H1.abf"], 4),
      PP.H2 = round(res$summary["PP.H2.abf"], 4),
      PP.H3 = round(res$summary["PP.H3.abf"], 4),
      PP.H4 = round(res$summary["PP.H4.abf"], 4),
      stringsAsFactors = FALSE
    )
    
    completed <- completed + 1
    
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    failed <<- failed + 1
  })
}

# ============================================================================
# Save results
# ============================================================================
if (length(results) > 0) {
  results_df <- do.call(rbind, results)
  rownames(results_df) <- NULL
  
  # Save full results
  write.table(results_df,
              file.path(output_dir, "coloc_conjfdr_all_results.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save summary for loci with PP.H4 > 0.5
  h4_loci <- results_df %>% filter(PP.H4 > 0.5)
  write.table(h4_loci,
              file.path(output_dir, "coloc_conjfdr_H4_significant.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Summary by pair
  pair_coloc_summary <- results_df %>%
    group_by(pair) %>%
    summarize(
      n_loci_tested = n(),
      n_H4_gt_0.5 = sum(PP.H4 > 0.5),
      n_H4_gt_0.8 = sum(PP.H4 > 0.8),
      n_H3_gt_0.5 = sum(PP.H3 > 0.5),
      mean_PP.H4 = round(mean(PP.H4), 3),
      .groups = "drop"
    )
  write.table(pair_coloc_summary,
              file.path(output_dir, "coloc_conjfdr_pair_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("\n============================================\n")
  cat("COLOC RESULTS SUMMARY\n")
  cat("============================================\n")
  cat(sprintf("Total loci tested: %d\n", completed))
  cat(sprintf("Failed: %d\n", failed))
  cat(sprintf("PP.H4 > 0.5 (shared causal variant): %d\n", sum(results_df$PP.H4 > 0.5)))
  cat(sprintf("PP.H4 > 0.8 (strong evidence): %d\n", sum(results_df$PP.H4 > 0.8)))
  cat(sprintf("PP.H3 > 0.5 (distinct variants): %d\n", sum(results_df$PP.H3 > 0.5)))
  
  cat("\n--- Per-pair summary ---\n")
  print(as.data.frame(pair_coloc_summary), row.names = FALSE)
  
  if (nrow(h4_loci) > 0) {
    cat("\n--- Loci with PP.H4 > 0.5 ---\n")
    print(h4_loci[, c("pair", "chr", "position", "nearest_gene", "conjfdr", "PP.H3", "PP.H4")],
          row.names = FALSE)
  }
  
} else {
  cat("\nNo coloc results generated.\n")
}

cat("\nDone! Results saved to:", output_dir, "\n")