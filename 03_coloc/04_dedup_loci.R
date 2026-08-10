#!/usr/bin/env Rscript
###############################################################################
# 04_dedup_loci.R  —  collapse row-level coloc.abf results to DISTINCT LOCI,
#                     each identified by SNP (not by nearest gene).
#
# Each locus is labelled by:
#   lead_snp  / lead_rsid  = the conjFDR sentinel variant (defines the locus)
#   cand_snp  / cand_rsid  = coloc credible-set top SNP (candidate SHARED variant)
# nearest_gene is retained only as a secondary annotation.
#
# rsIDs are filled if the input already carries lead_rsid/cand_rsid columns
# (i.e. you ran 06_annotate_rsid.R first); otherwise those columns are NA and
# only chr:pos is populated.
#
# Merge rule: same chromosome, +/-500 kb window overlap (gap <= MERGE_BP).
# Representative row per locus = the one with the highest PP.H4.
#
# Usage:
#   Rscript 04_dedup_loci.R final_abf_all/coloc_abf_all_results.tsv       final_abf_all
#   Rscript 04_dedup_loci.R final_abf_all/coloc_abf_all_results_rsid.tsv  final_abf_all   # after 06
###############################################################################

suppressPackageStartupMessages(library(dplyr))

args    <- commandArgs(TRUE)
infile  <- if (length(args) >= 1) args[1] else "final_abf_all/coloc_abf_all_results.tsv"
out_dir <- if (length(args) >= 2) args[2] else dirname(infile)
MERGE_BP <- 500000

x <- read.delim(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
for (c in c("chr","position","PP.H3","PP.H4","top_snp_pos"))
  if (c %in% names(x)) x[[c]] <- suppressWarnings(as.numeric(x[[c]]))
x$passes_yang <- as.logical(x$passes_yang)

# ensure identifier columns exist so summarise can always reference them
if (!"top_snp"   %in% names(x)) x$top_snp   <- ifelse(is.na(x$top_snp_pos), NA, paste0(x$chr, ":", x$top_snp_pos))
if (!"lead_rsid" %in% names(x)) x$lead_rsid <- NA_character_
if (!"cand_rsid" %in% names(x)) x$cand_rsid <- NA_character_
x$lead_snp <- paste0(x$chr, ":", x$position)          # conjFDR sentinel, chr:pos

cluster_ids <- function(chr, pos, gap) {
  o <- order(chr, pos); cid <- integer(length(pos)); k <- 0L
  pchr <- NA_real_; ppos <- NA_real_
  for (i in o) {
    if (is.na(pchr) || chr[i] != pchr || (pos[i] - ppos) > gap) k <- k + 1L
    cid[i] <- k; pchr <- chr[i]; ppos <- pos[i]
  }
  cid
}
classify <- function(maxH4, maxH3, any_yang)
  ifelse(maxH4 > 0.80, "colocalized",
    ifelse(any_yang,   "colocalized_Yang_0.70-0.80",
      ifelse(maxH3 > 0.80, "distinct", "inconclusive")))
gene_label <- function(g) {
  g <- g[!is.na(g) & trimws(g) != "" & toupper(g) != "NA"]
  if (!length(g)) "NA" else paste(unique(g), collapse = "/")
}

# ============================ Level A: within pair ==========================
xa <- x %>% group_by(pair) %>% mutate(loc = cluster_ids(chr, position, MERGE_BP)) %>% ungroup()

perpair_loci <- xa %>% group_by(pair, pd_stratum, kidney_trait, loc) %>%
  summarize(chr = chr[1],
            lead_snp  = lead_snp[which.max(PP.H4)],  lead_rsid = lead_rsid[which.max(PP.H4)],
            cand_snp  = top_snp[which.max(PP.H4)],   cand_rsid = cand_rsid[which.max(PP.H4)],
            maxPP.H3 = max(PP.H3), maxPP.H4 = max(PP.H4), any_yang = any(passes_yang),
            nearest_gene = gene_label(nearest_gene), n_rows = n(), .groups = "drop") %>%
  mutate(class = classify(maxPP.H4, maxPP.H3, any_yang))

perpair_summary <- perpair_loci %>% group_by(pair) %>%
  summarize(n_distinct_loci = n(), n_colocalized = sum(maxPP.H4 > 0.80),
            n_yang = sum(any_yang), n_distinct_sig = sum(class == "distinct"), .groups = "drop") %>%
  arrange(desc(n_colocalized))

# ============================ Level B: across pairs =========================
xb <- x %>% mutate(loc = cluster_ids(chr, position, MERGE_BP))

global_loci <- xb %>% group_by(loc) %>%
  summarize(chr = chr[1],
            lead_snp  = lead_snp[which.max(PP.H4)],  lead_rsid = lead_rsid[which.max(PP.H4)],
            cand_snp  = top_snp[which.max(PP.H4)],   cand_rsid = cand_rsid[which.max(PP.H4)],
            maxPP.H4 = max(PP.H4), maxPP.H3 = max(PP.H3), any_yang = any(passes_yang),
            best_pair_H4 = paste0(pd_stratum," x ",kidney_trait)[which.max(PP.H4)],
            coloc_traits = paste(sort(unique(kidney_trait[PP.H4 > 0.80])), collapse = "; "),
            n_pair_tests = n_distinct(pair),
            pairs_tested = paste(sort(unique(paste0(pd_stratum," x ",kidney_trait))), collapse = "; "),
            nearest_gene = gene_label(nearest_gene), .groups = "drop") %>%
  mutate(class = classify(maxPP.H4, maxPP.H3, any_yang)) %>%
  arrange(desc(maxPP.H4))

write.table(perpair_loci,    file.path(out_dir,"coloc_abf_perpair_distinct_loci.tsv"),  sep="\t", row.names=FALSE, quote=FALSE)
write.table(perpair_summary, file.path(out_dir,"coloc_abf_perpair_distinct_summary.tsv"),sep="\t", row.names=FALSE, quote=FALSE)
write.table(global_loci,     file.path(out_dir,"coloc_abf_global_distinct_loci.tsv"),    sep="\t", row.names=FALSE, quote=FALSE)

cat(sprintf("=== DISTINCT-LOCUS de-dup (merge <= %g kb), SNP-identified ===\n", MERGE_BP/1e3))
cat(sprintf("Row-level results: %d | within-pair loci: %d | global loci: %d\n",
            nrow(x), nrow(perpair_loci), nrow(global_loci)))
cat(sprintf("colocalized (PP.H4>0.80): %d | Yang: %d | distinct(PP.H3>0.80): %d\n",
            sum(global_loci$maxPP.H4>0.80), sum(global_loci$any_yang), sum(global_loci$class=="distinct")))
have_rsid <- any(!is.na(global_loci$lead_rsid))
cat(sprintf("rsIDs present: %s (run 06_annotate_rsid.R first to fill them)\n", ifelse(have_rsid,"yes","no - chr:pos only")))
cat("\n--- Global COLOCALIZED loci (PP.H4>0.80) ---\n")
print(as.data.frame(global_loci %>% filter(maxPP.H4>0.80) %>%
      select(lead_snp, lead_rsid, cand_snp, cand_rsid, maxPP.H4, best_pair_H4, coloc_traits, nearest_gene)), row.names=FALSE)