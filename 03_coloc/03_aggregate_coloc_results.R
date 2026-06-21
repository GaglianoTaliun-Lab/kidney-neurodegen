#!/usr/bin/env Rscript
###############################################################################
# 03_aggregate_coloc_results.R
#
# Aggregates per-locus coloc.susie results into final summary tables.
# Run after all SLURM array tasks complete.
#
# Usage:
#   Rscript 03_aggregate_coloc_results.R <results_dir> <final_output_dir>
###############################################################################

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(TRUE)
if (length(args) < 2) stop("Usage: Rscript 03_aggregate_coloc_results.R <results_dir> <final_output_dir>")

results_dir <- args[1]
final_dir   <- args[2]
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("Aggregating coloc.susie results\n")
cat("============================================\n\n")

# ============================================================================
# Read all per-locus results
# ============================================================================
result_files <- list.files(results_dir, pattern = "locus_.*_result\\.tsv$", full.names = TRUE)
skip_files   <- list.files(results_dir, pattern = "locus_.*_SKIP\\.tsv$", full.names = TRUE)
error_files  <- list.files(results_dir, pattern = "locus_.*_ERROR\\.tsv$", full.names = TRUE)

cat(sprintf("Result files: %d\n", length(result_files)))
cat(sprintf("Skipped: %d\n", length(skip_files)))
cat(sprintf("Errors: %d\n\n", length(error_files)))

if (length(result_files) == 0) {
  stop("No result files found. Check that SLURM jobs completed.")
}

# Combine all results
all_results <- do.call(rbind, lapply(result_files, function(f) {
  read.delim(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}))
all_results <- all_results %>% arrange(locus_idx)

cat(sprintf("Total loci with results: %d\n", nrow(all_results)))

# ============================================================================
# Save full results
# ============================================================================
write.table(all_results,
            file.path(final_dir, "coloc_susie_all_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================================
# Apply Yang et al. selection criteria
# ============================================================================
yang_pass <- all_results %>% filter(passes_yang == TRUE)
yang_h4   <- all_results %>% filter(PP.H4 > 0.8)
yang_h3   <- all_results %>% filter(PP.H3 > 0.8)

write.table(yang_pass,
            file.path(final_dir, "coloc_susie_yang_criteria_PASS.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(yang_h4,
            file.path(final_dir, "coloc_susie_PPH4_gt_0.8.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================================
# Aggregate credible sets
# ============================================================================
cs_files <- list.files(results_dir, pattern = "locus_.*_credset\\.tsv$", full.names = TRUE)
if (length(cs_files) > 0) {
  all_cs <- do.call(rbind, lapply(cs_files, function(f) {
    cs <- read.delim(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    # Extract locus index from filename
    idx <- as.integer(gsub(".*locus_(\\d+)_credset.*", "\\1", basename(f)))
    # Merge with locus metadata
    locus_meta <- all_results %>% filter(locus_idx == idx)
    if (nrow(locus_meta) > 0) {
      cs$locus_idx    <- idx
      cs$pair         <- locus_meta$pair[1]
      cs$chr          <- locus_meta$chr[1]
      cs$locus_center <- locus_meta$position[1]
      cs$nearest_gene <- locus_meta$nearest_gene[1]
      cs$PP.H4_locus  <- locus_meta$PP.H4[1]
    }
    cs
  }))

  write.table(all_cs,
              file.path(final_dir, "coloc_susie_all_credible_sets.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Top SNP per locus → for GTEx eQTL lookup
  top_snps <- all_cs %>%
    group_by(locus_idx) %>%
    slice_max(SNP.PP.H4, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(PP.H4_locus))

  write.table(top_snps,
              file.path(final_dir, "coloc_susie_top_SNPs_for_eQTL.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  cat(sprintf("\nCredible sets aggregated: %d loci\n", length(cs_files)))
  cat(sprintf("Top shared SNPs saved for eQTL lookup: %d\n", nrow(top_snps)))
}

# ============================================================================
# Per-pair summary
# ============================================================================
pair_summary <- all_results %>%
  group_by(pair) %>%
  summarize(
    n_tested        = n(),
    n_susie         = sum(method == "coloc.susie"),
    n_abf_fallback  = sum(method == "coloc.abf"),
    n_H4_gt_0.5     = sum(PP.H4 > 0.5),
    n_H4_gt_0.8     = sum(PP.H4 > 0.8),
    n_H3_gt_0.8     = sum(PP.H3 > 0.8),
    n_yang_pass     = sum(passes_yang == TRUE),
    mean_PP.H4      = round(mean(PP.H4), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(n_yang_pass))

write.table(pair_summary,
            file.path(final_dir, "coloc_susie_pair_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================================
# Print summary
# ============================================================================
cat("\n============================================\n")
cat("FINAL SUMMARY\n")
cat("============================================\n")
cat(sprintf("Total loci tested:              %d\n", nrow(all_results)))
cat(sprintf("Method: coloc.susie:            %d\n", sum(all_results$method == "coloc.susie")))
cat(sprintf("Method: coloc.abf (fallback):   %d\n", sum(all_results$method == "coloc.abf")))
cat(sprintf("PP.H4 > 0.5:                    %d\n", sum(all_results$PP.H4 > 0.5)))
cat(sprintf("PP.H4 > 0.8:                    %d\n", sum(all_results$PP.H4 > 0.8)))
cat(sprintf("PP.H3 > 0.8:                    %d\n", sum(all_results$PP.H3 > 0.8)))
cat(sprintf("Yang criteria (H3+H4>=0.8 & ratio>=5): %d\n", sum(all_results$passes_yang == TRUE)))
cat(sprintf("Skipped loci:                   %d\n", length(skip_files)))
cat(sprintf("Error loci:                     %d\n", length(error_files)))

cat("\n--- Per-pair summary ---\n")
print(as.data.frame(pair_summary), row.names = FALSE)

if (nrow(yang_pass) > 0) {
  cat("\n--- Loci passing Yang et al. criteria ---\n")
  print(yang_pass[, c("pair", "chr", "position", "nearest_gene",
                       "method", "PP.H3", "PP.H4", "ratio_H4H3",
                       "top_snp", "top_snp_pp")],
        row.names = FALSE)
}

cat("\nAll results saved to:", final_dir, "\n")
cat("\nNEXT STEP: Look up SNPs in coloc_susie_top_SNPs_for_eQTL.tsv\n")
cat("  on GTEx Portal (gtexportal.org) or Open Targets (genetics.opentargets.org)\n")
cat("  to identify which gene each shared causal variant regulates.\n")
