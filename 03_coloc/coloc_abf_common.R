###############################################################################
# coloc_abf_common.R
#
# Shared config + data loading for the coloc.abf pipeline.
# Sourced by run_coloc_abf_single_locus.R and run_coloc_abf_sensitivity.R.
#
# Handles heterogeneous sumstats via column-name normalisation + chr:pos parsing:
#   - CKDGen METAL (eGFR/eGFRcys): MarkerName Allele1 Allele2 n Freq1 Effect StdErr P.value ... chr pos ...
#   - regenie (eGFRcrea strat, uK/Cr, uNa/Cr): CHROM GENPOS ... A1FREQ INFO N ... BETA SE ... LOG10P [P]
#   - harmonised tsv (uACR): CHR BP SNP A2 A1 N BETA P MAF SE
#   - PD sex-combined (.tbl): SNP CHR BP A1 A2 BETA SE A1_AF P N
#   - PD sex-stratified (METAL): MarkerName ... Effect StdErr P-value ... (chr:pos in MarkerName; NO N column)
###############################################################################

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)   # fast reader; handles .gz + whitespace/tab
})

# ---- Primary priors ---------------------------------------------------------
P1 <- 1e-4; P2 <- 1e-4; P12 <- 5e-6

# ---- Window / QC ------------------------------------------------------------
WINDOW   <- 500000
MIN_SNPS <- 100
INFO_MIN <- 0.8        # drop imputed variants with INFO < this (if INFO present)
PP_H4_CS <- 0.5

# ---- sdY / N / frequency handling -------------------------------------------
# RINT'd quantitative traits -> sdY = 1 (uACR used --apply-rint).
trait_sdY_map <- c("uACR" = 1, "uACR (female)" = 1, "uACR (male)" = 1)
SDY_QUANT_DEFAULT <- NULL   # forbid silent assumptions for other quant traits
# Total N for traits whose files carry NO per-SNP N column (PD sex-stratified METAL).
trait_N_map <- c("PD (female)" = 28277, "PD (male)" = 32356)
FREQ_COL_OVERRIDE <- NULL

# ---- Column-name candidates (case-insensitive) ------------------------------
CHR_CANDS    <- c("chromosome","CHROM","CHR","chr","#CHROM","Chr")
POS_CANDS    <- c("base_pair_location","GENPOS","BP","POS","pos","position","Pos_b37","POS_b37")
MARKER_CANDS <- c("MarkerName","markername","marker","MarkerID","variant_id","variant","SNPID")
BETA_CANDS   <- c("beta","BETA","B","Effect","EFFECT","effect")
SE_CANDS     <- c("standard_error","SE","se","StdErr","stderr")
P_CANDS      <- c("p_value","P","p","PVAL","pval","P.value","P-value","Pvalue")
LOG10P_CANDS <- c("LOG10P","log10p","LOG10_P","neglog10p")
N_CANDS      <- c("N","n","n_total_sum","N_total","n_total","Neff","TotalSampleSize")
INFO_CANDS   <- c("INFO","info","Rsq","R2","impinfo")
FREQ_CANDS   <- c("A1FREQ","a1freq","effect_allele_frequency","eaf","EAF",
                  "Freq1","FREQ1","freq1","AF1","af1","MAF","maf","freq","FRQ",
                  "A1_FREQ","ALT_FREQS","AF","minor_AF")

# ---- Trait -> file map ------------------------------------------------------
# Bare filenames resolve against SUMSTATS_DIR (the coloc submit arg = kidney_neurodegen).
# Absolute paths (starting with / or ~) are used as-is  -> PD lives in a different dir.
PD_DIR <- "/home/lchang24/projects/def-gsarah/sumstats/PD-Blauwendraat-2021"
trait_file_map <- c(
  "PD (meta)"       = file.path(PD_DIR, "SEXCOMBINED_PD_meta-analysis_NO_PROXY.tbl"),
  "PD (female)"     = file.path(PD_DIR, "FEMALE_PD_filtered_sumstats_NO_PROXY_no_multi_allelics_RSID.txt.gz"),
  "PD (male)"       = file.path(PD_DIR, "MALE_PD_filtered_sumstats_NO_PROXY_MAF_no_multi_allelics_RSID.txt.gz"),
  "eGFR (meta)"     = "metal_eGFR_meta_ea1.TBL.map.annot.gc.gz",
  "eGFRcys (meta)"  = "metal_eGFRcys_meta1.TBL.map.annot.gc.gz",
  "eGFR (female)"   = "eGFRcrea_female_HRC_imputed_regenie_allchrs.regenie",
  "eGFR (male)"     = "eGFRcrea_male_HRC_imputed_regenie_allchrs.regenie",
  "uACR"            = "uacr_sex_combined.tsv",
  "uACR (female)"   = "uacr_female.tsv",
  "uACR (male)"     = "uacr_male.tsv",
  "uK/Cr"           = "uk_cr_combined_uKCr_all_chr.txt.gz",
  "uK/Cr (female)"  = "uk_cr_female_uKCr_all_chr.txt.gz",
  "uK/Cr (male)"    = "uk_cr_male_uKCr_all_chr.txt.gz",
  "uNa/Cr"          = "una_cr_combined_uNaCr_all_chr_with_P.txt.gz",
  "uNa/Cr (female)" = "una_cr_female_uNaCr_all_chr_with_P.txt.gz",
  "uNa/Cr (male)"   = "una_cr_male_uNaCr_all_chr_with_P.txt.gz",
  "Hematuria"       = "hematuria_sexcombined.tsv",                              # <-- CONFIRM columns before the full run
  "Hematuria (female)" = "X593.FEMALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz",    # <-- CONFIRM
  "Hematuria (male)"   = "X593.MALES.ukb_v3.SAIGE.MAC_20.INFO_0.4.txt.gz"       # <-- CONFIRM
)

trait_type_map <- c(
  "PD (meta)" = "cc", "PD (female)" = "cc", "PD (male)" = "cc",
  "eGFR (meta)" = "quant", "eGFR (female)" = "quant", "eGFR (male)" = "quant",
  "eGFRcys (meta)" = "quant",
  "uACR" = "quant", "uACR (female)" = "quant", "uACR (male)" = "quant",
  "Hematuria" = "cc", "Hematuria (female)" = "cc", "Hematuria (male)" = "cc",
  "uK/Cr" = "quant", "uK/Cr (female)" = "quant", "uK/Cr (male)" = "quant",
  "uNa/Cr" = "quant", "uNa/Cr (female)" = "quant", "uNa/Cr (male)" = "quant"
)

case_props <- c(
  "PD (meta)" = 20967/60633, "PD (female)" = 7947/28277, "PD (male)" = 13020/32356,
  "Hematuria" = 16409/396345, "Hematuria (female)" = 7061/220094, "Hematuria (male)" = 9774/187440
)

# ============================================================================
# Column helpers
# ============================================================================
.find_col <- function(nm, cands) {
  for (c in cands) { i <- which(tolower(nm) == tolower(c)); if (length(i)) return(nm[i[1]]) }
  NA_character_
}

# chromosome + position, from explicit columns OR parsed from a chr:pos marker
.get_chr_pos <- function(ss) {
  nm <- names(ss)
  chr_c <- .find_col(nm, CHR_CANDS); pos_c <- .find_col(nm, POS_CANDS)
  if (!is.na(chr_c) && !is.na(pos_c)) {
    chr <- suppressWarnings(as.integer(gsub("^chr", "", as.character(ss[[chr_c]]), ignore.case = TRUE)))
    pos <- suppressWarnings(as.numeric(ss[[pos_c]]))
  } else {
    mk_c <- .find_col(nm, MARKER_CANDS); if (is.na(mk_c)) mk_c <- nm[1]
    v <- as.character(ss[[mk_c]])
    chk <- v[seq_len(min(200L, length(v)))]
    if (mean(grepl("^[0-9XYMxym]+:[0-9]+", chk), na.rm = TRUE) < 0.9)
      stop(sprintf("No chr/pos columns and '%s' is not chr:pos formatted.", mk_c))
    chr <- suppressWarnings(as.integer(gsub("^chr", "", sub(":.*$", "", v), ignore.case = TRUE)))
    pos <- suppressWarnings(as.numeric(sub("^[^:]*:([0-9]+).*$", "\\1", v)))
  }
  list(chr = chr, pos = pos)
}

normalize_cols <- function(ss, trait_label) {
  nm <- names(ss)
  cp <- .get_chr_pos(ss)
  beta_c <- .find_col(nm, BETA_CANDS); se_c <- .find_col(nm, SE_CANDS)
  if (is.na(beta_c) || is.na(se_c))
    stop(sprintf("Trait '%s': missing beta/se column. Header: %s", trait_label, paste(nm, collapse = ", ")))
  out <- data.frame(
    chromosome         = cp$chr,
    base_pair_location = cp$pos,
    beta               = suppressWarnings(as.numeric(ss[[beta_c]])),
    standard_error     = suppressWarnings(as.numeric(ss[[se_c]])),
    stringsAsFactors   = FALSE
  )
  n_c <- .find_col(nm, N_CANDS)
  out$N <- if (!is.na(n_c)) suppressWarnings(as.numeric(ss[[n_c]])) else NA_real_
  p_c <- .find_col(nm, P_CANDS); l_c <- .find_col(nm, LOG10P_CANDS)
  if (!is.na(p_c))      out$p_value <- suppressWarnings(as.numeric(ss[[p_c]]))
  else if (!is.na(l_c)) out$p_value <- 10^(-suppressWarnings(as.numeric(ss[[l_c]])))
  else                  out$p_value <- NA_real_
  f_c <- if (!is.null(FREQ_COL_OVERRIDE) && FREQ_COL_OVERRIDE %in% nm) FREQ_COL_OVERRIDE else .find_col(nm, FREQ_CANDS)
  if (!is.na(f_c)) out$freq <- suppressWarnings(as.numeric(ss[[f_c]]))
  i_c <- .find_col(nm, INFO_CANDS)
  if (!is.na(i_c)) out$INFO <- suppressWarnings(as.numeric(ss[[i_c]]))
  out
}

# ============================================================================
# Read a +/- WINDOW region for one trait (window-filter first, then normalise)
# ============================================================================
read_window_sumstats <- function(trait_label, chr_num, center_pos, sumstats_dir) {
  fname <- trait_file_map[trait_label]
  if (is.na(fname)) stop(paste("No file mapping for trait label:", trait_label))
  fpath <- if (grepl("^(/|~)", fname)) path.expand(fname) else file.path(sumstats_dir, fname)
  if (!file.exists(fpath)) stop(paste("File not found:", fpath))

  is_gz <- grepl("\\.gz$", fpath, ignore.case = TRUE)
  base_cmd <- if (is_gz) paste("zcat", shQuote(fpath)) else paste("cat", shQuote(fpath))
  # Peek at the header. Some regenie files (the uNa/Cr *_with_P) are space-delimited
  # but have a stray TAB before the appended P column; fread auto-detects that lone
  # tab as THE separator and collapses the row into 1-2 columns. When the header
  # holds BOTH a tab and a space, squeeze all blank runs to a single space, read sep=" ".
  con <- if (is_gz) gzfile(fpath) else file(fpath)
  hdr <- tryCatch(readLines(con, n = 1L), error = function(e) character(0), finally = close(con))
  mixed_delim <- length(hdr) == 1L && grepl("\t", hdr) && grepl(" ", hdr)
  dt <- if (mixed_delim)
          data.table::fread(cmd = paste(base_cmd, "| tr -s '[:blank:]' ' '"),
                            sep = " ", showProgress = FALSE)
        else if (is_gz)
          data.table::fread(cmd = base_cmd, showProgress = FALSE)
        else
          data.table::fread(fpath, showProgress = FALSE)

  cp <- .get_chr_pos(as.data.frame(dt))
  keep <- !is.na(cp$chr) & !is.na(cp$pos) & cp$chr == chr_num &
          cp$pos >= (center_pos - WINDOW) & cp$pos <= (center_pos + WINDOW)
  dt <- dt[which(keep)]                       # small window subset
  if (nrow(dt) == 0) return(dt[0])

  ss <- normalize_cols(as.data.frame(dt), trait_label)
  ss <- ss[!is.na(ss$beta) & !is.na(ss$standard_error) & ss$standard_error > 0, , drop = FALSE]
  if ("INFO" %in% names(ss)) ss <- ss[!is.na(ss$INFO) & ss$INFO >= INFO_MIN, , drop = FALSE]
  ss <- ss[!duplicated(ss$base_pair_location), , drop = FALSE]
  ss[order(ss$base_pair_location), , drop = FALSE]
}

# ============================================================================
# Dataset construction (NO LD)
# ============================================================================
trait_uses_maf <- function(type, trait_label) type == "quant" && !(trait_label %in% names(trait_sdY_map))
get_maf <- function(ss) if ("freq" %in% names(ss)) pmin(ss$freq, 1 - ss$freq) else NULL

build_dataset <- function(ss, type, trait_label, maf_vec = NULL) {
  n <- suppressWarnings(median(ss$N, na.rm = TRUE))
  if (is.na(n) && trait_label %in% names(trait_N_map)) n <- trait_N_map[[trait_label]]
  ds <- list(beta = ss$beta, varbeta = ss$standard_error^2, snp = ss$snp_id,
             position = ss$base_pair_location, type = type)
  if (!is.na(n)) ds$N <- n
  if (type == "quant") {
    if (trait_label %in% names(trait_sdY_map))      ds$sdY <- as.numeric(trait_sdY_map[[trait_label]])
    else if (!is.null(maf_vec))                     ds$MAF <- maf_vec
    else if (!is.null(SDY_QUANT_DEFAULT))           ds$sdY <- SDY_QUANT_DEFAULT
    else stop(sprintf("Quant trait '%s': no sdY entry and no frequency column.", trait_label))
    if (is.na(n)) stop(sprintf("Quant trait '%s': N is required (for sdY estimation) but was not found.", trait_label))
  }
  if (type == "cc" && trait_label %in% names(case_props)) ds$s <- unname(case_props[trait_label])
  ds
}

extract_credible_set <- function(coloc_res) {
  r <- coloc_res$results
  o <- order(r$SNP.PP.H4, decreasing = TRUE)
  w <- which(cumsum(r$SNP.PP.H4[o]) > 0.95)[1]
  cols <- intersect(c("snp", "position", "SNP.PP.H4"), names(r))
  r[o[1:w], cols, drop = FALSE]
}