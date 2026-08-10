#!/usr/bin/env Rscript
###############################################################################
# run_coloc_abf_single_locus.R
#
# coloc.abf (LD-free, single-causal-variant) for ONE locus. Drop-in for the
# SLURM array. Config + data loading live in coloc_abf_common.R (sourced below).
#
# Usage:
#   Rscript run_coloc_abf_single_locus.R <locus_index> <loci_file> <sumstats_dir> <output_dir>
###############################################################################

suppressPackageStartupMessages(library(dplyr))

# --- locate + source the shared helper (must sit next to this script) --------
get_script_dir <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f[1])) else getwd()
}
common <- file.path(get_script_dir(), "coloc_abf_common.R")
if (!file.exists(common)) stop("coloc_abf_common.R not found next to this script: ", common)
source(common)

args <- commandArgs(TRUE)
if (length(args) < 4) stop("Usage: Rscript run_coloc_abf_single_locus.R <locus_index> <loci_file> <sumstats_dir> <output_dir>")
locus_idx    <- as.integer(args[1])
loci_file    <- args[2]
sumstats_dir <- args[3]
output_dir   <- args[4]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

all_loci <- read.delim(loci_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  filter(!(chr == 6 & position >= 25000000 & position <= 34000000))
if (locus_idx < 1 || locus_idx > nrow(all_loci)) stop("Locus index out of range. Total: ", nrow(all_loci))
row <- all_loci[locus_idx, ]
cat(sprintf("=== Locus %d/%d: %s x %s (chr%s:%s | %s) ===\n",
            locus_idx, nrow(all_loci), row$pd_stratum, row$kidney_trait, row$chr, row$position, row$nearest_gene))

write_status <- function(status, extra = list()) {
  df <- data.frame(locus_idx = locus_idx, pair = paste0(row$pd_stratum, " x ", row$kidney_trait),
                   chr = row$chr, position = row$position, nearest_gene = row$nearest_gene,
                   status = status, stringsAsFactors = FALSE)
  for (nm in names(extra)) df[[nm]] <- extra[[nm]]
  suffix <- if (grepl("^ERROR", status)) "ERROR" else "SKIP"
  write.table(df, file.path(output_dir, sprintf("locus_%04d_%s.tsv", locus_idx, suffix)),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

tryCatch({
  type1 <- unname(trait_type_map[row$pd_stratum]); type2 <- unname(trait_type_map[row$kidney_trait])
  if (is.na(type1)) stop("No type for trait: ", row$pd_stratum)
  if (is.na(type2)) stop("No type for trait: ", row$kidney_trait)

  ss1 <- read_window_sumstats(row$pd_stratum,   row$chr, row$position, sumstats_dir)
  ss2 <- read_window_sumstats(row$kidney_trait, row$chr, row$position, sumstats_dir)
  ss1$snp_id <- paste0(row$chr, ":", ss1$base_pair_location)
  ss2$snp_id <- paste0(row$chr, ":", ss2$base_pair_location)
  cat(sprintf("  SNPs in window: trait1=%d trait2=%d\n", nrow(ss1), nrow(ss2)))

  shared <- sort(unique(intersect(ss1$snp_id, ss2$snp_id)))
  if (length(shared) < MIN_SNPS) { write_status("SKIPPED_FEW_SNPS", list(n_shared = length(shared))); quit(save = "no", status = 0) }
  ss1 <- ss1[match(shared, ss1$snp_id), ]; ss2 <- ss2[match(shared, ss2$snp_id), ]

  # MAF (only for quant traits not in trait_sdY_map) + joint validity filter
  valid <- rep(TRUE, length(shared)); maf1 <- maf2 <- NULL
  if (trait_uses_maf(type1, row$pd_stratum)) {
    maf1 <- get_maf(ss1); if (is.null(maf1)) stop(sprintf("Quant trait '%s' has no frequency column.", row$pd_stratum))
    valid <- valid & is.finite(maf1) & maf1 > 0 & maf1 <= 0.5
  }
  if (trait_uses_maf(type2, row$kidney_trait)) {
    maf2 <- get_maf(ss2); if (is.null(maf2)) stop(sprintf("Quant trait '%s' has no frequency column.", row$kidney_trait))
    valid <- valid & is.finite(maf2) & maf2 > 0 & maf2 <= 0.5
  }
  ss1 <- ss1[valid, ]; ss2 <- ss2[valid, ]
  if (!is.null(maf1)) maf1 <- maf1[valid]
  if (!is.null(maf2)) maf2 <- maf2[valid]
  shared <- shared[valid]
  cat(sprintf("  Shared SNPs after QC/MAF filter: %d\n", length(shared)))
  if (length(shared) < MIN_SNPS) { write_status("SKIPPED_FEW_SNPS_POST_FILTER", list(n_shared = length(shared))); quit(save = "no", status = 0) }

  dataset1 <- build_dataset(ss1, type1, row$pd_stratum,   maf1)
  dataset2 <- build_dataset(ss2, type2, row$kidney_trait, maf2)
  invisible(check_dataset(dataset1, warn.minp = 1e-8))
  invisible(check_dataset(dataset2, warn.minp = 1e-8))

  res <- suppressMessages(coloc.abf(dataset1, dataset2, p1 = P1, p2 = P2, p12 = P12))
  s <- res$summary
  pp_h0 <- unname(s["PP.H0.abf"]); pp_h1 <- unname(s["PP.H1.abf"]); pp_h2 <- unname(s["PP.H2.abf"])
  pp_h3 <- unname(s["PP.H3.abf"]); pp_h4 <- unname(s["PP.H4.abf"])
  sum_h3h4 <- pp_h3 + pp_h4; ratio_h4h3 <- ifelse(pp_h3 > 0, pp_h4 / pp_h3, Inf)
  passes_yang <- (sum_h3h4 >= 0.8) & (ratio_h4h3 >= 5)
  cat(sprintf("  PP.H3=%.4f PP.H4=%.4f | Yang=%s | H4>=0.70=%s | H4>0.80=%s\n",
              pp_h3, pp_h4, passes_yang, pp_h4 >= 0.70, pp_h4 > 0.80))

  top_snp_id <- NA; top_snp_pos <- NA; top_snp_pp <- NA; cs_size <- NA
  if (pp_h4 > PP_H4_CS && "results" %in% names(res)) {
    cs_df <- extract_credible_set(res)
    top_snp_id <- cs_df$snp[1]
    if ("position" %in% names(cs_df)) top_snp_pos <- cs_df$position[1]
    top_snp_pp <- cs_df$SNP.PP.H4[1]; cs_size <- nrow(cs_df)
    write.table(cs_df, file.path(output_dir, sprintf("locus_%04d_credset.tsv", locus_idx)),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }

  result_df <- data.frame(
    locus_idx = locus_idx, pair = paste0(row$pd_stratum, " x ", row$kidney_trait),
    pd_stratum = row$pd_stratum, kidney_trait = row$kidney_trait,
    locus_num = row$locus_num, chr = row$chr, position = row$position,
    conjfdr = row$conjfdr, nearest_gene = row$nearest_gene,
    method = "coloc.abf", n_shared_snps = length(shared), n_cs_pairs = 1,
    PP.H0 = round(pp_h0, 6), PP.H1 = round(pp_h1, 6), PP.H2 = round(pp_h2, 6),
    PP.H3 = round(pp_h3, 6), PP.H4 = round(pp_h4, 6),
    sum_H3H4 = round(sum_h3h4, 6), ratio_H4H3 = round(min(ratio_h4h3, 9999), 2),
    passes_yang = passes_yang, PPH4_ge_0.70 = pp_h4 >= 0.70, PPH4_gt_0.80 = pp_h4 > 0.80,
    top_snp = top_snp_id, top_snp_pos = top_snp_pos,
    top_snp_pp = round(ifelse(is.na(top_snp_pp), 0, top_snp_pp), 6), cs_size_95 = cs_size,
    p1 = P1, p2 = P2, p12 = P12, stringsAsFactors = FALSE
  )
  write.table(result_df, file.path(output_dir, sprintf("locus_%04d_result.tsv", locus_idx)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved.\n")

}, error = function(e) { cat(sprintf("FATAL: %s\n", e$message)); write_status(paste("ERROR:", e$message)) })

cat("Done.\n")