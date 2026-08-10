#!/usr/bin/env Rscript
###############################################################################
# build_slalom_snp.R
#
# Build SLALOM (Kanai et al. 2022) per-locus .snp inputs from the REPORTED
# coloc.abf loci + the GWAS summary statistics, so we can run the reviewer's own
# tool as QC for Comment 1 (meta-analysis fine-mapping miscalibration).
#
# One .snp file per locus x trait. SLALOM QC is per-trait: it picks that trait's
# lead (min-p) in the window and asks whether every other variant's Z is
# consistent with the lead given local LD (DENTIST-S). So each GWAS is its own
# .snp file / its own SLALOM run.
#
# SLALOM .snp format (whitespace-delimited, header row). Required columns:
#   chromosome position allele1 allele2 beta se
# We also write p (SLALOM's default lead = min-p) and, best-effort, rsid maf
# n_samples n_cases (enables SLALOM's min_neff_r2/max_neff_r2 power columns).
#
# CRITICAL allele convention (README):
#   * allele2 is ALWAYS the EFFECT allele (the one beta/z is expressed for).
#   * allele1 is the OTHER allele.
#   * --align-alleles will reconcile ref/alt against gnomAD and sign-flip beta as
#     needed, but the effect stays on allele2 -- so we only need to put the effect
#     allele in allele2, not know the true b37 ref.
#   * chromosome is a BARE INTEGER (1..22) for GRCh37 (gnomAD v2 b37 contig
#     convention). Run SLALOM with --reference-genome GRCh37.
#
# Reuses coloc_abf_common.R (trait_file_map, WINDOW, trait_N_map, case_props,
# .find_col, .get_chr_pos, and the *_CANDS column-name vectors) so trait handling
# is identical to the coloc run.
#
# Usage:
#   Rscript build_slalom_snp.R <loci_file> <sumstats_dir> <out_dir> [locus_indices]
#     <loci_file>    SuppTable_conjfdr_all_loci_with_eGFRcys.tsv  (same file/indices
#                    as run_coloc_abf_sensitivity.R -- MHC filtered, then indexed)
#     <sumstats_dir> kidney_neurodegen sumstats dir (PD files use absolute paths)
#     <out_dir>      where .snp files + manifest.tsv are written
#     [locus_indices] optional comma-separated override (default = the 10 reported)
###############################################################################

suppressPackageStartupMessages(library(data.table))

get_script_dir <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f[1])) else getwd()
}
source(file.path(get_script_dir(), "coloc_abf_common.R"))  # trait maps, WINDOW, helpers, *_CANDS

# ------------------------------------------------------------------ config ---
# Reported loci to QC = the loci we make claims about + the SCARB2 distinct
# control. SAME indices as 05_submit_sensitivity.sh / run_coloc_abf_sensitivity.R
# (resolved AFTER the MHC filter, exactly as that script does).
TARGET_LOCI_DEFAULT <- c(1, 6, 20, 24, 53, 65, 71, 72, 79, 148)

# EFFECT-allele candidates (allele beta/z is for) and the OTHER allele, per the
# formats in this pipeline: METAL Allele1=effect, regenie ALLELE1=effect,
# PD .tbl A1=effect, harmonised uACR A1=effect.
A1E_CANDS  <- c("Allele1","ALLELE1","A1","effect_allele","EA","ALT","a1")   # EFFECT
A2E_CANDS  <- c("Allele2","ALLELE0","A2","other_allele","NEA","REF","a2")    # OTHER
RSID_CANDS <- c("rsid","RSID","SNP","MarkerName","ID","variant_id","rs_id","SNPID","snp")

args <- commandArgs(TRUE)
loci_file    <- if (length(args) >= 1) args[1] else stop("give <loci_file>")
sumstats_dir <- if (length(args) >= 2) args[2] else stop("give <sumstats_dir>")
out_dir      <- if (length(args) >= 3) args[3] else "slalom_snp"
TARGET_LOCI  <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else TARGET_LOCI_DEFAULT
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------- window reader ---
# Read the +/- WINDOW window for one trait, keeping alleles + p + (best-effort)
# rsid/freq/N. Mirrors read_window_sumstats() but keeps the columns SLALOM needs.
read_window_slalom <- function(trait_label, chr_num, center_pos, sumstats_dir) {
  fname <- trait_file_map[trait_label]
  if (is.na(fname)) stop("no file mapping for trait: ", trait_label)
  fpath <- if (grepl("^(/|~)", fname)) path.expand(fname) else file.path(sumstats_dir, fname)
  if (!file.exists(fpath)) stop("file not found: ", fpath)
  is_gz <- grepl("\\.gz$", fpath, ignore.case = TRUE)
  base_cmd <- if (is_gz) paste("zcat", shQuote(fpath)) else paste("cat", shQuote(fpath))
  con <- if (is_gz) gzfile(fpath) else file(fpath)
  hdr <- tryCatch(readLines(con, n = 1L), error = function(e) character(0), finally = close(con))
  mixed <- length(hdr) == 1L && grepl("\t", hdr) && grepl(" ", hdr)
  dt <- if (mixed) fread(cmd = paste(base_cmd, "| tr -s '[:blank:]' ' '"), sep = " ", showProgress = FALSE)
        else if (is_gz) fread(cmd = base_cmd, showProgress = FALSE)
        else fread(fpath, showProgress = FALSE)
  df <- as.data.frame(dt); nm <- names(df)

  cp  <- .get_chr_pos(df)                                   # chr + pos (coloc_abf_common)
  a1c <- .find_col(nm, A1E_CANDS); a2c <- .find_col(nm, A2E_CANDS)
  bc  <- .find_col(nm, BETA_CANDS); sc <- .find_col(nm, SE_CANDS)
  if (any(is.na(c(a1c, a2c, bc, sc))))
    stop(sprintf("%s: missing effect/other-allele or beta/se column. header: %s",
                 trait_label, paste(nm, collapse = ",")))

  keep <- !is.na(cp$chr) & !is.na(cp$pos) & cp$chr == chr_num &
          cp$pos >= center_pos - WINDOW & cp$pos <= center_pos + WINDOW
  ki <- which(keep); df <- df[ki, , drop = FALSE]; cp <- lapply(cp, function(v) v[ki])

  beta <- suppressWarnings(as.numeric(df[[bc]]))
  se   <- suppressWarnings(as.numeric(df[[sc]]))
  eff  <- toupper(as.character(df[[a1c]]))   # EFFECT allele -> SLALOM allele2
  oth  <- toupper(as.character(df[[a2c]]))   # OTHER allele  -> SLALOM allele1

  # p: prefer file P, else 10^-log10p, else two-sided from z.
  pc <- .find_col(nm, P_CANDS); lc <- .find_col(nm, LOG10P_CANDS)
  if (!is.na(pc))      p <- suppressWarnings(as.numeric(df[[pc]]))
  else if (!is.na(lc)) p <- 10^(-suppressWarnings(as.numeric(df[[lc]])))
  else                 p <- 2 * pnorm(-abs(beta / se))

  out <- data.frame(
    chromosome = cp$chr,                       # bare integer (b37)
    position   = cp$pos,
    allele1    = oth,                           # OTHER allele
    allele2    = eff,                           # EFFECT allele (beta's allele)
    beta       = beta,
    se         = se,
    p          = p,
    stringsAsFactors = FALSE
  )

  # optional: rsid
  rc <- .find_col(nm, RSID_CANDS)
  if (!is.na(rc)) {
    rs <- as.character(df[[rc]])
    # only keep if it actually looks like an rsID (MarkerName can be chr:pos)
    if (mean(grepl("^rs[0-9]+$", rs[seq_len(min(200L, length(rs)))]), na.rm = TRUE) > 0.5)
      out$rsid <- rs
  }
  # optional: maf (from freq column)
  fc <- .find_col(nm, FREQ_CANDS)
  if (!is.na(fc)) {
    fr <- suppressWarnings(as.numeric(df[[fc]]))
    out$maf <- pmin(fr, 1 - fr)
  }
  # optional: n_samples (per-SNP N if present, else fixed total from trait_N_map)
  nc <- .find_col(nm, N_CANDS)
  if (!is.na(nc)) out$n_samples <- suppressWarnings(as.numeric(df[[nc]]))
  else if (trait_label %in% names(trait_N_map)) out$n_samples <- as.numeric(trait_N_map[[trait_label]])
  # optional: n_cases for case-control traits (constant case fraction x n_samples)
  if (trait_label %in% names(case_props) && "n_samples" %in% names(out))
    out$n_cases <- round(unname(case_props[trait_label]) * out$n_samples)

  # clean: finite, positive se, biallelic SNVs, one row per position
  ok <- is.finite(out$beta) & is.finite(out$se) & out$se > 0 &
        nchar(out$allele1) == 1 & nchar(out$allele2) == 1 &
        out$allele1 %in% c("A","C","G","T") & out$allele2 %in% c("A","C","G","T")
  out <- out[ok, , drop = FALSE]
  out <- out[!duplicated(out$position), , drop = FALSE]
  out[order(out$position), , drop = FALSE]
}

# ------------------------------------------------------------------- main ----
all_loci <- read.delim(loci_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# SAME MHC filter + positional indexing as run_coloc_abf_sensitivity.R
all_loci <- all_loci[!(all_loci$chr == 6 & all_loci$position >= 25000000 & all_loci$position <= 34000000), ]

manifest <- list()
for (li in TARGET_LOCI) {
  if (li < 1 || li > nrow(all_loci)) { cat(sprintf("[skip] index %d out of range\n", li)); next }
  row  <- all_loci[li, ]
  gene <- gsub("[^A-Za-z0-9]", "", as.character(row$nearest_gene)); if (is.na(gene) || gene == "" || gene == "NA") gene <- "NA"
  for (tr in c(row$pd_stratum, row$kidney_trait)) {
    trec <- gsub("[^A-Za-z0-9]", "", tr)
    snp  <- file.path(out_dir, sprintf("locus%03d_chr%s_%s_%s.snp", li, row$chr, gene, trec))
    ss <- tryCatch(read_window_slalom(tr, row$chr, row$position, sumstats_dir),
                   error = function(e) { cat(sprintf("  [%s @ locus %d] %s\n", tr, li, conditionMessage(e))); NULL })
    if (is.null(ss) || nrow(ss) < 100) {
      cat(sprintf("  [%s @ locus %d] too few variants (%s) -- skipped\n", tr, li,
                  if (is.null(ss)) "read failed" else nrow(ss))); next
    }
    write.table(ss, snp, sep = " ", row.names = FALSE, quote = FALSE)
    ttype <- if (tr %in% names(case_props)) "cc" else "quant"
    manifest[[length(manifest) + 1]] <- data.frame(
      snp_file = normalizePath(snp), locus_idx = li, chr = row$chr, pos = row$position,
      gene = gene, trait = tr, trait_type = ttype, n_variants = nrow(ss), stringsAsFactors = FALSE)
    cat(sprintf("  wrote %-45s  n=%d  (%s)\n", basename(snp), nrow(ss), ttype))
  }
}
man <- do.call(rbind, manifest)
write.table(man, file.path(out_dir, "manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\n%d .snp files + manifest.tsv written to %s\n", nrow(man), out_dir))
cat("Next: run_slalom_loci.sh <out_dir>/manifest.tsv <results_dir>\n")
