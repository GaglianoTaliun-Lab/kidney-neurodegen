#!/usr/bin/env Rscript
###############################################################################
# run_coloc_abf_sensitivity.R
#
# Prior-sensitivity analysis for coloc.abf at ONE locus. Config + data loading
# come from coloc_abf_common.R (sourced below), so results are directly
# comparable to run_coloc_abf_single_locus.R.
#
# Outputs per locus:
#   locus_XXXX_sensitivity.tsv        -- coloc.abf across a grid of priors
#   plots/locus_XXXX_sensitivity.pdf  -- coloc's sensitivity() vs p12, using the
#                                        Yang rule, for candidate shared loci only
#
# Usage:
#   Rscript run_coloc_abf_sensitivity.R <locus_index> <loci_file> <sumstats_dir> <output_dir>
###############################################################################

suppressPackageStartupMessages(library(dplyr))

get_script_dir <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f[1])) else getwd()
}
common <- file.path(get_script_dir(), "coloc_abf_common.R")
if (!file.exists(common)) stop("coloc_abf_common.R not found next to this script: ", common)
source(common)

args <- commandArgs(TRUE)
if (length(args) < 4) stop("Usage: Rscript run_coloc_abf_sensitivity.R <locus_index> <loci_file> <sumstats_dir> <output_dir>")
locus_idx    <- as.integer(args[1]); loci_file <- args[2]
sumstats_dir <- args[3]; output_dir <- args[4]
plots_dir <- file.path(output_dir, "plots"); dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Selection rule in coloc sensitivity() syntax = the Yang criterion
SENS_RULE <- "H4 > 5*H3 & H3+H4 > 0.8"
# Prior grid: sweep p12 at p1=p2=1e-4 (all < p1 to avoid the boundary warning),
# plus coloc's default and the Frida point (p1=p2=1e-3, p12=1e-4).
PRIOR_GRID <- data.frame(
  label = c("p12_1e-6", "p12_5e-6(primary)", "p12_1e-5(coloc default)",
            "p12_5e-5", "frida_p1p2_1e-3_p12_1e-4"),
  p1  = c(1e-4, 1e-4, 1e-4, 1e-4, 1e-3),
  p2  = c(1e-4, 1e-4, 1e-4, 1e-4, 1e-3),
  p12 = c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4),
  stringsAsFactors = FALSE
)

all_loci <- read.delim(loci_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  filter(!(chr == 6 & position >= 25000000 & position <= 34000000))
if (locus_idx < 1 || locus_idx > nrow(all_loci)) stop("Locus index out of range.")
row <- all_loci[locus_idx, ]
cat(sprintf("=== Sensitivity | Locus %d/%d: %s x %s (chr%s:%s) ===\n",
            locus_idx, nrow(all_loci), row$pd_stratum, row$kidney_trait, row$chr, row$position))

result <- tryCatch({
  type1 <- unname(trait_type_map[row$pd_stratum]); type2 <- unname(trait_type_map[row$kidney_trait])
  ss1 <- read_window_sumstats(row$pd_stratum,   row$chr, row$position, sumstats_dir)
  ss2 <- read_window_sumstats(row$kidney_trait, row$chr, row$position, sumstats_dir)
  ss1$snp_id <- paste0(row$chr, ":", ss1$base_pair_location)
  ss2$snp_id <- paste0(row$chr, ":", ss2$base_pair_location)
  shared <- sort(unique(intersect(ss1$snp_id, ss2$snp_id)))
  if (length(shared) < MIN_SNPS) stop(sprintf("Too few shared SNPs (%d)", length(shared)))
  ss1 <- ss1[match(shared, ss1$snp_id), ]; ss2 <- ss2[match(shared, ss2$snp_id), ]

  valid <- rep(TRUE, length(shared)); maf1 <- maf2 <- NULL
  if (trait_uses_maf(type1, row$pd_stratum)) { maf1 <- get_maf(ss1); if (is.null(maf1)) stop("trait1 needs frequency"); valid <- valid & is.finite(maf1) & maf1 > 0 & maf1 <= 0.5 }
  if (trait_uses_maf(type2, row$kidney_trait)) { maf2 <- get_maf(ss2); if (is.null(maf2)) stop("trait2 needs frequency"); valid <- valid & is.finite(maf2) & maf2 > 0 & maf2 <= 0.5 }
  ss1 <- ss1[valid, ]; ss2 <- ss2[valid, ]
  if (!is.null(maf1)) maf1 <- maf1[valid]
  if (!is.null(maf2)) maf2 <- maf2[valid]
  shared <- shared[valid]
  if (length(shared) < MIN_SNPS) stop(sprintf("Too few SNPs after filter (%d)", length(shared)))

  dataset1 <- build_dataset(ss1, type1, row$pd_stratum,   maf1)
  dataset2 <- build_dataset(ss2, type2, row$kidney_trait, maf2)

  grid_rows <- lapply(seq_len(nrow(PRIOR_GRID)), function(i) {
    g <- PRIOR_GRID[i, ]
    s <- suppressMessages(coloc.abf(dataset1, dataset2, p1 = g$p1, p2 = g$p2, p12 = g$p12))$summary
    h3 <- unname(s["PP.H3.abf"]); h4 <- unname(s["PP.H4.abf"])
    data.frame(locus_idx = locus_idx, pair = paste0(row$pd_stratum, " x ", row$kidney_trait),
               chr = row$chr, position = row$position, nearest_gene = row$nearest_gene,
               prior_label = g$label, p1 = g$p1, p2 = g$p2, p12 = g$p12,
               PP.H0 = round(unname(s["PP.H0.abf"]), 6), PP.H1 = round(unname(s["PP.H1.abf"]), 6),
               PP.H2 = round(unname(s["PP.H2.abf"]), 6), PP.H3 = round(h3, 6), PP.H4 = round(h4, 6),
               sum_H3H4 = round(h3 + h4, 6), ratio_H4H3 = round(ifelse(h3 > 0, min(h4 / h3, 9999), 9999), 2),
               passes_yang = (h3 + h4 >= 0.8) & (ifelse(h3 > 0, h4 / h3, Inf) >= 5),
               PPH4_ge_0.70 = h4 >= 0.70, PPH4_gt_0.80 = h4 > 0.80, stringsAsFactors = FALSE)
  })
  sens_df <- do.call(rbind, grid_rows)
  write.table(sens_df, file.path(output_dir, sprintf("locus_%04d_sensitivity.tsv", locus_idx)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  PP.H4 %.3f-%.3f across %d priors; Yang-pass %d/%d\n",
              min(sens_df$PP.H4), max(sens_df$PP.H4), sum(sens_df$passes_yang), nrow(sens_df)))

  primary <- suppressMessages(coloc.abf(dataset1, dataset2, p1 = P1, p2 = P2, p12 = P12))
  if (unname(primary$summary["PP.H4.abf"]) > PP_H4_CS) {
    pdf(file.path(plots_dir, sprintf("locus_%04d_sensitivity.pdf", locus_idx)), width = 8, height = 6)
    tryCatch(sensitivity(primary, rule = SENS_RULE),
             error = function(e) cat(sprintf("  sensitivity() plot failed: %s\n", e$message)))
    dev.off()
  }
  "OK"
}, error = function(e) {
  cat(sprintf("FATAL: %s\n", e$message))
  write.table(data.frame(locus_idx = locus_idx, pair = paste0(row$pd_stratum, " x ", row$kidney_trait),
                         status = paste("ERROR:", e$message)),
              file.path(output_dir, sprintf("locus_%04d_sensitivity_ERROR.tsv", locus_idx)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  "ERROR"
})
cat(sprintf("Done (%s).\n", result))
