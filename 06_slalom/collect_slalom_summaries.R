#!/usr/bin/env Rscript
###############################################################################
# collect_slalom_summaries.R  --  aggregate SLALOM per-locus summaries into ONE
# Supplementary Table for Reviewer-2 Comment 1.
#
#   Usage: Rscript collect_slalom_summaries.R <results_dir> [out.tsv]
#
# Reads every <base>.summary.txt written by run_slalom_loci.sh, joins the
# locus/trait provenance from the paired <base>.meta, and emits one row per
# locus x trait with the SLALOM verdict.
#
# SLALOM verdict (Kanai et al. 2022): a locus is "suspicious" for meta-analysis
# fine-mapping if it has >=1 DENTIST-S outlier variant -- a variant in LD with
# the lead (r^2 > 0.6) whose Z is inconsistent with the lead given LD
# (DENTIST-S P < 1e-4). n_dentist_s_outlier == 0  =>  NOT suspicious (clean).
# Base R only (no extra packages).
###############################################################################

args <- commandArgs(TRUE)
results_dir <- if (length(args) >= 1) args[1] else "results_slalom"
out_tsv     <- if (length(args) >= 2) args[2] else "slalom_supplementary_table.tsv"

summ_files <- list.files(results_dir, pattern = "\\.summary\\.txt$", full.names = TRUE)
if (!length(summ_files)) stop("no *.summary.txt in ", results_dir)

rd1 <- function(f) read.delim(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                              check.names = FALSE, colClasses = "character")

rows <- list()
for (sf in summ_files) {
  base <- sub("\\.summary\\.txt$", "", basename(sf))
  s <- tryCatch(rd1(sf), error = function(e) NULL)
  if (is.null(s) || !nrow(s)) { cat("[skip] unreadable/empty:", basename(sf), "\n"); next }
  mf <- file.path(results_dir, paste0(base, ".meta"))
  m  <- if (file.exists(mf)) rd1(mf) else data.frame(locus_idx = NA, gene = NA, trait = NA,
                                                     trait_type = NA, chr = NA, pos = NA)
  rows[[length(rows) + 1]] <- cbind(m, s, file = base, stringsAsFactors = FALSE)
}
if (!length(rows)) stop("no summaries could be parsed")

# union columns across rows, then rbind
allcols <- unique(unlist(lapply(rows, names)))
rows <- lapply(rows, function(d) { d[setdiff(allcols, names(d))] <- NA; d[allcols] })
tab <- do.call(rbind, rows)

num <- function(x) suppressWarnings(as.numeric(x))
nout <- num(tab[["n_dentist_s_outlier"]])
tab$verdict <- ifelse(is.na(nout), "no LD / not evaluable",
                ifelse(nout == 0, "consistent (no suspicious variant)",
                       paste0("FLAGGED (", nout, " DENTIST-S outlier", ifelse(nout == 1, "", "s"), ")")))

# order + select the reporting columns (keep neff / nonsyn if present)
want <- c("locus_idx","gene","trait","trait_type","chr","pos",
          "lead_pip_variant","n_total","n_r2","n_na","n_dentist_s_outlier","fraction","max_pip",
          "min_neff_r2","max_neff_r2","n_nonsyn","n_nonsyn_outlier","cs_nonsyn","verdict")
want <- want[want %in% names(tab)]
tab <- tab[, c(want, setdiff(names(tab), c(want, "file"))), drop = FALSE]
ord <- order(num(tab$locus_idx), tab$trait); tab <- tab[ord, , drop = FALSE]

write.table(tab, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# ---- console report ---------------------------------------------------------
nflag <- sum(!is.na(nout) & nout > 0); nclean <- sum(!is.na(nout) & nout == 0); nna <- sum(is.na(nout))
cat(sprintf("\n=== SLALOM summary: %d locus x trait runs ===\n", nrow(tab)))
cat(sprintf("  consistent (0 outliers): %d | FLAGGED (>=1): %d | not evaluable: %d\n", nclean, nflag, nna))
cat(sprintf("  -> %s\n", out_tsv))
show <- c("locus_idx","gene","trait","n_r2","n_dentist_s_outlier","max_pip","verdict")
show <- show[show %in% names(tab)]
cat("\n"); print(tab[, show], row.names = FALSE)
if (nflag > 0) cat("\nNOTE: FLAGGED rows are the loci to inspect/caveat in the response. ",
                   "An eGFR/eGFRcys-side flag is the meaningful one; a PD-side flag on a ",
                   "UKB-proxy LD panel is more ambiguous.\n", sep = "")
