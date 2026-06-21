#!/usr/bin/env Rscript
###############################################################################
# run_coloc_susie_single_locus.R
#
# Runs coloc.susie for a SINGLE locus (called by SLURM array).
#
# Implements Yang et al. (Brain Communications 2026) approach:
#   - SuSiE fine-mapping with LD matrix from 1000G Phase 3 EUR
#   - MAF from 1000G EUR for sdY estimation (quantitative traits)
#   - Conservative priors: p1=p2=1e-4, p12=5e-6
#   - Selection: PP.H3+PP.H4 >= 0.8 AND PP.H4/PP.H3 >= 5
#   - PIP extraction for top shared causal variant
#
# Usage:
#   Rscript run_coloc_susie_single_locus.R \
#     <locus_index> <loci_file> <sumstats_dir> <ld_dir> <output_dir>
###############################################################################

suppressPackageStartupMessages({
  library(coloc)
  library(susieR)
  library(dplyr)
})

args <- commandArgs(TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript run_coloc_susie_single_locus.R <locus_index> <loci_file> <sumstats_dir> <ld_dir> <output_dir>")
}

locus_idx    <- as.integer(args[1])
loci_file    <- args[2]
sumstats_dir <- args[3]
ld_dir       <- args[4]
output_dir   <- args[5]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Parameters
# ============================================================================
P1     <- 1e-4
P2     <- 1e-4
P12    <- 5e-6
WINDOW <- 500000

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

case_props <- c(
  "PD (meta)"          = 20967 / 60633,
  "PD (female)"        = 7947  / 28277,
  "PD (male)"          = 13020 / 32356,
  "Hematuria"          = 16409 / 396345,
  "Hematuria (female)" = 7061  / 220094,
  "Hematuria (male)"   = 9774  / 187440
)

# ============================================================================
# Read locus info
# ============================================================================
all_loci <- read.delim(loci_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
all_loci <- all_loci %>%
  filter(!(chr == 6 & position >= 25000000 & position <= 34000000))

if (locus_idx < 1 || locus_idx > nrow(all_loci)) {
  stop(paste("Locus index", locus_idx, "out of range. Total loci:", nrow(all_loci)))
}

row <- all_loci[locus_idx, ]

cat(sprintf("=== Locus %d/%d ===\n", locus_idx, nrow(all_loci)))
cat(sprintf("Pair: %s x %s\n", row$pd_stratum, row$kidney_trait))
cat(sprintf("Position: chr%s:%s\n", row$chr, row$position))
cat(sprintf("Gene: %s\n\n", row$nearest_gene))

# ============================================================================
# Helper: read sumstats for a window
# ============================================================================
read_window_sumstats <- function(trait_label, chr_num, center_pos) {
  fname <- trait_file_map[trait_label]
  fpath <- file.path(sumstats_dir, fname)
  if (!file.exists(fpath)) stop(paste("File not found:", fpath))

  ss <- read.delim(fpath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  ss %>%
    filter(chromosome == chr_num,
           base_pair_location >= (center_pos - WINDOW),
           base_pair_location <= (center_pos + WINDOW)) %>%
    filter(!is.na(beta), !is.na(standard_error), !is.na(p_value),
           standard_error > 0) %>%
    distinct(base_pair_location, .keep_all = TRUE) %>%
    arrange(base_pair_location)
}

# ============================================================================
# Helper: read LD matrix + MAF for a locus
# ============================================================================
read_ld_and_maf <- function(chr_num, center_pos) {
  prefix <- file.path(ld_dir, sprintf("ld_chr%s_%s", chr_num, center_pos))

  # --- LD matrix ---
  ld_file  <- paste0(prefix, ".ld")
  if (!file.exists(ld_file)) stop(paste("LD file not found:", ld_file))
  ld_mat <- as.matrix(read.table(ld_file, header = FALSE))

  # --- BIM positions (CHR, SNP, CM, BP, A1, A2) ---
  bim_file <- paste0(prefix, ".bim_positions")
  if (!file.exists(bim_file)) stop(paste("BIM file not found:", bim_file))
  bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("chr", "snp_plink", "cm", "bp", "a1", "a2"))

  # Create chr:pos IDs for LD matrix
  ld_pos_ids <- paste0(bim$chr, ":", bim$bp)

  # --- FIX: Remove duplicate positions (multi-allelic sites in 1000G) ---
  dup_mask <- !duplicated(ld_pos_ids)
  n_dups   <- sum(!dup_mask)
  if (n_dups > 0) {
    cat(sprintf("  Removing %d duplicate positions from LD matrix\n", n_dups))
    ld_mat     <- ld_mat[dup_mask, dup_mask]
    bim        <- bim[dup_mask, ]
    ld_pos_ids <- ld_pos_ids[dup_mask]
  }

  colnames(ld_mat) <- rownames(ld_mat) <- ld_pos_ids

  # --- MAF from .frq file (PLINK --freq output) ---
  frq_file <- paste0(prefix, ".frq")
  maf_lookup <- NULL

  if (file.exists(frq_file)) {
    # PLINK .frq format: CHR SNP A1 A2 MAF NCHROBS (space-delimited, leading whitespace)
    frq <- read.table(frq_file, header = TRUE, stringsAsFactors = FALSE)
    # frq$SNP matches bim$snp_plink in the same order
    # Create position-based MAF lookup
    if (nrow(frq) == length(dup_mask)) {
      # Apply same deduplication as LD matrix
      frq <- frq[dup_mask, ]
    }
    if (nrow(frq) == nrow(bim)) {
      maf_lookup <- data.frame(
        pos_id = ld_pos_ids,
        bp     = bim$bp,
        MAF    = frq$MAF,
        stringsAsFactors = FALSE
      )
      cat(sprintf("  MAF loaded: %d SNPs (median MAF=%.3f)\n",
                  nrow(maf_lookup), median(maf_lookup$MAF, na.rm = TRUE)))
    } else {
      warning("FRQ and BIM row counts don't match; MAF not used")
    }
  } else {
    cat("  Note: no .frq file found; sdY=1 will be used for quantitative traits\n")
  }

  cat(sprintf("  LD matrix: %d x %d SNPs\n", nrow(ld_mat), ncol(ld_mat)))

  return(list(ld = ld_mat, maf = maf_lookup))
}

# ============================================================================
# Helper: extract 95% credible set
# ============================================================================
extract_credible_set <- function(coloc_res) {
  snp_res <- coloc_res$results
  o  <- order(snp_res$SNP.PP.H4, decreasing = TRUE)
  cs <- cumsum(snp_res$SNP.PP.H4[o])
  w  <- which(cs > 0.95)[1]
  snp_res[o[1:w], c("snp", "position", "SNP.PP.H4")]
}

# ============================================================================
# Main analysis
# ============================================================================
tryCatch({

  # Read sumstats
  ss1 <- read_window_sumstats(row$pd_stratum,  row$chr, row$position)
  ss2 <- read_window_sumstats(row$kidney_trait, row$chr, row$position)
  cat(sprintf("  Trait 1 SNPs in window: %d\n", nrow(ss1)))
  cat(sprintf("  Trait 2 SNPs in window: %d\n", nrow(ss2)))

  # Read LD + MAF
  ld_maf <- read_ld_and_maf(row$chr, row$position)
  ld_mat     <- ld_maf$ld
  maf_lookup <- ld_maf$maf

  # Create position-based IDs for sumstats
  ss1$snp_id <- paste0(row$chr, ":", ss1$base_pair_location)
  ss2$snp_id <- paste0(row$chr, ":", ss2$base_pair_location)

  # Find SNPs shared across sumstats + LD (deduplicated)
  shared_snps <- Reduce(intersect, list(ss1$snp_id, ss2$snp_id, colnames(ld_mat)))
  shared_snps <- unique(sort(shared_snps))
  cat(sprintf("  SNPs shared (sumstats + LD): %d\n", length(shared_snps)))

  if (length(shared_snps) < 100) {
    cat("  SKIP: too few shared SNPs\n")
    write.table(
      data.frame(locus_idx = locus_idx, pair = paste0(row$pd_stratum, " x ", row$kidney_trait),
                 chr = row$chr, position = row$position, nearest_gene = row$nearest_gene,
                 status = "SKIPPED_FEW_SNPS", n_shared = length(shared_snps)),
      file.path(output_dir, sprintf("locus_%04d_SKIP.tsv", locus_idx)),
      sep = "\t", row.names = FALSE, quote = FALSE)
    quit(save = "no", status = 0)
  }

  # Subset and deduplicate sumstats
  ss1 <- ss1 %>% filter(snp_id %in% shared_snps)
  ss2 <- ss2 %>% filter(snp_id %in% shared_snps)
  ss1 <- ss1[!duplicated(ss1$snp_id), ]
  ss2 <- ss2[!duplicated(ss2$snp_id), ]

  # Re-intersect after dedup
  shared_snps <- Reduce(intersect, list(ss1$snp_id, ss2$snp_id, colnames(ld_mat)))
  shared_snps <- unique(sort(shared_snps))

  # Align everything using match() for correct ordering
  ss1    <- ss1[match(shared_snps, ss1$snp_id), ]
  ss2    <- ss2[match(shared_snps, ss2$snp_id), ]
  ld_idx <- match(shared_snps, colnames(ld_mat))
  ld_mat <- ld_mat[ld_idx, ld_idx]

  # Verify alignment
  if (!all(ss1$snp_id == colnames(ld_mat)) || !all(ss2$snp_id == colnames(ld_mat))) {
    stop("SNP alignment failed after deduplication — check for remaining duplicates")
  }
  cat(sprintf("  Aligned: %d SNPs\n", length(shared_snps)))

  # Merge MAF by position using match()
  maf_vec <- NULL
  if (!is.null(maf_lookup)) {
    maf_idx <- match(shared_snps, maf_lookup$pos_id)
    if (all(!is.na(maf_idx))) {
      maf_vec <- maf_lookup$MAF[maf_idx]
      cat(sprintf("  MAF merged: %d SNPs matched\n", length(maf_vec)))
    } else {
      n_miss <- sum(is.na(maf_idx))
      cat(sprintf("  Warning: %d SNPs missing MAF; sdY=1 fallback\n", n_miss))
    }
  }

  n1    <- median(ss1$N, na.rm = TRUE)
  n2    <- median(ss2$N, na.rm = TRUE)
  type1 <- trait_type_map[row$pd_stratum]
  type2 <- trait_type_map[row$kidney_trait]

  # Build dataset 1
  dataset1 <- list(
    beta     = ss1$beta,
    varbeta  = ss1$standard_error^2,
    snp      = ss1$snp_id,
    position = ss1$base_pair_location,
    type     = type1,
    N        = n1,
    LD       = ld_mat
  )

  # Build dataset 2
  dataset2 <- list(
    beta     = ss2$beta,
    varbeta  = ss2$standard_error^2,
    snp      = ss2$snp_id,
    position = ss2$base_pair_location,
    type     = type2,
    N        = n2,
    LD       = ld_mat
  )

  # --- MAF for quantitative traits ---
  if (type1 == "quant") {
    if (!is.null(maf_vec)) {
      dataset1$MAF <- maf_vec
      cat(sprintf("  Trait 1 (%s): MAF provided for sdY estimation\n", row$pd_stratum))
    } else {
      dataset1$sdY <- 1
      cat(sprintf("  Trait 1 (%s): no MAF, using sdY=1\n", row$pd_stratum))
    }
  }
  if (type2 == "quant") {
    if (!is.null(maf_vec)) {
      dataset2$MAF <- maf_vec
      cat(sprintf("  Trait 2 (%s): MAF provided for sdY estimation\n", row$kidney_trait))
    } else {
      dataset2$sdY <- 1
      cat(sprintf("  Trait 2 (%s): no MAF, using sdY=1\n", row$kidney_trait))
    }
  }

  # Case proportions for case-control traits
  if (type1 == "cc" && row$pd_stratum %in% names(case_props)) {
    dataset1$s <- case_props[row$pd_stratum]
  }
  if (type2 == "cc" && row$kidney_trait %in% names(case_props)) {
    dataset2$s <- case_props[row$kidney_trait]
  }

  # ------------------------------------------------------------------
  # Run SuSiE fine-mapping
  # ------------------------------------------------------------------
  cat("  Running SuSiE for trait 1...\n")
  S1 <- tryCatch(runsusie(dataset1, max_iter = 200),
                 error = function(e) { cat(sprintf("  SuSiE failed (trait 1): %s\n", e$message)); NULL })

  cat("  Running SuSiE for trait 2...\n")
  S2 <- tryCatch(runsusie(dataset2, max_iter = 200),
                 error = function(e) { cat(sprintf("  SuSiE failed (trait 2): %s\n", e$message)); NULL })

  # ------------------------------------------------------------------
  # Run colocalization
  # ------------------------------------------------------------------
  use_susie <- !is.null(S1) && !is.null(S2)

  if (use_susie) {
    cat("  Running coloc.susie...\n")
    res <- tryCatch(coloc.susie(S1, S2, p1 = P1, p2 = P2, p12 = P12),
                    error = function(e) { cat(sprintf("  coloc.susie failed: %s\n", e$message)); NULL })

    if (!is.null(res) && !is.null(res$summary) && nrow(res$summary) > 0) {
      method <- "coloc.susie"
      best   <- res$summary[which.max(res$summary$PP.H4.abf), ]
      pp_h0 <- best$PP.H0.abf; pp_h1 <- best$PP.H1.abf
      pp_h2 <- best$PP.H2.abf; pp_h3 <- best$PP.H3.abf
      pp_h4 <- best$PP.H4.abf
      hit1  <- best$hit1; hit2 <- best$hit2
      n_cs_pairs <- nrow(res$summary)
      cat(sprintf("  coloc.susie: %d credible set pairs | best PP.H4=%.4f\n", n_cs_pairs, pp_h4))
    } else {
      use_susie <- FALSE
      cat("  coloc.susie: no credible sets found, falling back\n")
    }
  }

  # Fallback to coloc.abf
  if (!use_susie) {
    cat("  Running coloc.abf (fallback)...\n")
    method <- "coloc.abf"
    res <- suppressMessages(coloc.abf(dataset1, dataset2, p1 = P1, p2 = P2, p12 = P12))
    pp_h0 <- res$summary["PP.H0.abf"]; pp_h1 <- res$summary["PP.H1.abf"]
    pp_h2 <- res$summary["PP.H2.abf"]; pp_h3 <- res$summary["PP.H3.abf"]
    pp_h4 <- res$summary["PP.H4.abf"]
    hit1 <- NA; hit2 <- NA; n_cs_pairs <- NA
  }

  cat(sprintf("  PP.H3=%.4f | PP.H4=%.4f | method=%s\n", pp_h3, pp_h4, method))

  # Yang et al. criteria
  sum_h3h4    <- pp_h3 + pp_h4
  ratio_h4h3  <- ifelse(pp_h3 > 0, pp_h4 / pp_h3, Inf)
  passes_yang <- (sum_h3h4 >= 0.8) & (ratio_h4h3 >= 5)
  cat(sprintf("  Yang criteria: H3+H4=%.4f, ratio=%.1f -> %s\n",
              sum_h3h4, min(ratio_h4h3, 9999), ifelse(passes_yang, "PASS", "FAIL")))

  # Extract top shared causal SNP
  top_snp_id <- NA; top_snp_pos <- NA; top_snp_pp <- NA; cs_size <- NA

  if (pp_h4 > 0.5 && method == "coloc.abf" && "results" %in% names(res)) {
    cs_df <- extract_credible_set(res)
    top_snp_id  <- cs_df$snp[1]
    top_snp_pos <- cs_df$position[1]
    top_snp_pp  <- cs_df$SNP.PP.H4[1]
    cs_size     <- nrow(cs_df)
    write.table(cs_df, file.path(output_dir, sprintf("locus_%04d_credset.tsv", locus_idx)),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  Top SNP: %s (PP=%.3f, CS size=%d)\n", top_snp_id, top_snp_pp, cs_size))
  } else if (pp_h4 > 0.5 && method == "coloc.susie") {
    top_snp_id <- hit2; top_snp_pp <- pp_h4
    write.table(res$summary, file.path(output_dir, sprintf("locus_%04d_susie_summary.tsv", locus_idx)),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  Top hit: %s (PP.H4=%.3f)\n", hit2, pp_h4))
  }

  # Save result
  result_df <- data.frame(
    locus_idx = locus_idx,
    pair = paste0(row$pd_stratum, " x ", row$kidney_trait),
    pd_stratum = row$pd_stratum, kidney_trait = row$kidney_trait,
    locus_num = row$locus_num, chr = row$chr, position = row$position,
    conjfdr = row$conjfdr, nearest_gene = row$nearest_gene,
    method = method, n_shared_snps = length(shared_snps),
    n_cs_pairs = ifelse(is.na(n_cs_pairs), 1, n_cs_pairs),
    PP.H0 = round(pp_h0, 6), PP.H1 = round(pp_h1, 6),
    PP.H2 = round(pp_h2, 6), PP.H3 = round(pp_h3, 6), PP.H4 = round(pp_h4, 6),
    sum_H3H4 = round(sum_h3h4, 6),
    ratio_H4H3 = round(min(ratio_h4h3, 9999), 2),
    passes_yang = passes_yang,
    top_snp = top_snp_id, top_snp_pos = top_snp_pos,
    top_snp_pp = round(ifelse(is.na(top_snp_pp), 0, top_snp_pp), 6),
    cs_size_95 = cs_size,
    stringsAsFactors = FALSE
  )
  write.table(result_df, file.path(output_dir, sprintf("locus_%04d_result.tsv", locus_idx)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved.\n")

}, error = function(e) {
  cat(sprintf("FATAL: %s\n", e$message))
  write.table(
    data.frame(locus_idx = locus_idx, pair = paste0(row$pd_stratum, " x ", row$kidney_trait),
               chr = row$chr, position = row$position, nearest_gene = row$nearest_gene,
               status = paste("ERROR:", e$message)),
    file.path(output_dir, sprintf("locus_%04d_ERROR.tsv", locus_idx)),
    sep = "\t", row.names = FALSE, quote = FALSE)
})

cat("Done.\n")
