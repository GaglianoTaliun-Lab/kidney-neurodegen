#!/usr/bin/env Rscript
###############################################################################
# 03_aggregate_coloc_abf.R
#
# Aggregates per-locus coloc.abf results AND the prior-sensitivity sweep.
#
# Usage:
#   Rscript 03_aggregate_coloc_abf.R <results_dir> <final_output_dir>
#
# <results_dir> may contain primary results (locus_*_result.tsv) and/or
# sensitivity results (locus_*_sensitivity.tsv). Point it at whichever you ran
# (or a directory holding both).
###############################################################################

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(TRUE)
if (length(args) < 2) stop("Usage: Rscript 03_aggregate_coloc_abf.R <results_dir> <final_output_dir>")
results_dir <- args[1]
final_dir   <- args[2]
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv_safe <- function(f) read.delim(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# ============================================================================
# PART A: primary coloc.abf results
# ============================================================================
result_files <- list.files(results_dir, pattern = "locus_.*_result\\.tsv$", full.names = TRUE)
skip_files   <- list.files(results_dir, pattern = "locus_.*_SKIP\\.tsv$",   full.names = TRUE)
error_files  <- list.files(results_dir, pattern = "locus_.*_ERROR\\.tsv$",  full.names = TRUE)

if (length(result_files) > 0) {
  all_results <- do.call(rbind, lapply(result_files, read_tsv_safe)) %>% arrange(locus_idx)
  write.table(all_results, file.path(final_dir, "coloc_abf_all_results.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  write.table(all_results %>% filter(passes_yang == TRUE),
              file.path(final_dir, "coloc_abf_yang_criteria_PASS.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(all_results %>% filter(PP.H4 > 0.8),
              file.path(final_dir, "coloc_abf_PPH4_gt_0.8.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(all_results %>% filter(PP.H4 >= 0.70),
              file.path(final_dir, "coloc_abf_PPH4_ge_0.70.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # credible sets -> top SNP per locus for eQTL lookup
  cs_files <- list.files(results_dir, pattern = "locus_.*_credset\\.tsv$", full.names = TRUE)
  if (length(cs_files) > 0) {
    all_cs <- do.call(rbind, lapply(cs_files, function(f) {
      cs  <- read_tsv_safe(f)
      idx <- as.integer(gsub(".*locus_(\\d+)_credset.*", "\\1", basename(f)))
      meta <- all_results %>% filter(locus_idx == idx)
      if (nrow(meta) > 0) {
        cs$locus_idx <- idx; cs$pair <- meta$pair[1]; cs$chr <- meta$chr[1]
        cs$locus_center <- meta$position[1]; cs$nearest_gene <- meta$nearest_gene[1]
        cs$PP.H4_locus <- meta$PP.H4[1]
      }
      cs
    }))
    write.table(all_cs, file.path(final_dir, "coloc_abf_all_credible_sets.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    top_snps <- all_cs %>% group_by(locus_idx) %>%
      slice_max(SNP.PP.H4, n = 1, with_ties = FALSE) %>% ungroup() %>%
      arrange(desc(PP.H4_locus))
    write.table(top_snps, file.path(final_dir, "coloc_abf_top_SNPs_for_eQTL.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }

  pair_summary <- all_results %>% group_by(pair) %>%
    summarize(n_tested = n(),
              n_H4_gt_0.5  = sum(PP.H4 > 0.5),
              n_H4_ge_0.70 = sum(PP.H4 >= 0.70),
              n_H4_gt_0.80 = sum(PP.H4 > 0.80),
              n_H3_gt_0.80 = sum(PP.H3 > 0.80),
              n_yang_pass  = sum(passes_yang == TRUE),
              mean_PP.H4   = round(mean(PP.H4), 4), .groups = "drop") %>%
    arrange(desc(n_yang_pass))
  write.table(pair_summary, file.path(final_dir, "coloc_abf_pair_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  cat("=== coloc.abf primary ===\n")
  cat(sprintf("Loci with results: %d | SKIP: %d | ERROR: %d\n",
              nrow(all_results), length(skip_files), length(error_files)))
  cat(sprintf("PP.H4 > 0.80: %d | PP.H4 >= 0.70: %d | PP.H3 > 0.80: %d | Yang pass: %d\n",
              sum(all_results$PP.H4 > 0.8), sum(all_results$PP.H4 >= 0.70),
              sum(all_results$PP.H3 > 0.8), sum(all_results$passes_yang == TRUE)))
} else {
  cat("No primary result files (locus_*_result.tsv) found in results_dir.\n")
}

# ============================================================================
# PART B: prior-sensitivity sweep -> stability summary
# ============================================================================
sens_files <- list.files(results_dir, pattern = "locus_.*_sensitivity\\.tsv$", full.names = TRUE)
if (length(sens_files) > 0) {
  all_sens <- do.call(rbind, lapply(sens_files, read_tsv_safe))
  write.table(all_sens, file.path(final_dir, "coloc_abf_sensitivity_all.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  sens_summary <- all_sens %>% group_by(locus_idx, pair, chr, position, nearest_gene) %>%
    summarize(n_priors        = n(),
              PP.H4_min        = round(min(PP.H4), 4),
              PP.H4_max        = round(max(PP.H4), 4),
              PP.H4_median     = round(median(PP.H4), 4),
              n_yang_pass      = sum(passes_yang == TRUE),
              n_H4_ge_0.70     = sum(PPH4_ge_0.70 == TRUE),
              robust_yang_all  = all(passes_yang == TRUE),   # passes under EVERY prior
              robust_yang_none = all(passes_yang == FALSE),  # fails under EVERY prior
              .groups = "drop") %>%
    arrange(desc(PP.H4_median))
  write.table(sens_summary, file.path(final_dir, "coloc_abf_sensitivity_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  cat("\n=== coloc.abf prior sensitivity ===\n")
  cat(sprintf("Loci swept: %d | robust (Yang pass under ALL priors): %d | Yang under NO prior: %d\n",
              nrow(sens_summary), sum(sens_summary$robust_yang_all),
              sum(sens_summary$robust_yang_none)))
}

cat("\nSaved to:", final_dir, "\n")
