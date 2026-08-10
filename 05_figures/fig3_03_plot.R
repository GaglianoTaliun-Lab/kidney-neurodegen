#!/usr/bin/env Rscript
###############################################################################
# fig3_03_plot.R   (mirrors the manuscript figure3/03_plot_figure3.R)
#
# Paired regional association plots (PD on top, kidney below) for the coloc.abf
# results, using locuszoomr gg_scatter + gg_genetracks, SNPs coloured by PLINK
# r^2 to the index (candidate shared) variant, with a gene track.
#
# Reads windows from plot_data/ (written by 07_extract_locus_regions.R; columns
# SNP CHR BP P) and LD from ld/ (written by fig3_02_compute_LD.sh).
#
# Panels:
#   MAIN (2x2):  RAB19 | BIN3xeGFR | BIN3xeGFRcys | SCARB2 (distinct control)
#   SUPP (1x2):  MSRA/PRSS51 | TNK2
#
# Run on narval:
#   module load StdEnv/2023 r/4.4.0
#   Rscript fig3_03_plot.R
###############################################################################
suppressPackageStartupMessages({
  library(locuszoomr); library(EnsDb.Hsapiens.v75); library(ggplot2); library(cowplot)
})

BASE    <- "/home/lchang24/projects/def-gsarah/lchang24/github/coloc-abf"
WIN_DIR <- file.path(BASE, "plot_data")
LD_DIR  <- file.path(BASE, "ld")
OUT_DIR <- file.path(BASE, "figure3_output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# cairo devices so the ♂/♀ glyphs render on the cluster (default device fails: mbcsToSbcs)
save_pdf <- function(file, plot, w, h) { cairo_pdf(file, width=w, height=h); print(plot); dev.off() }
save_png <- function(file, plot, w, h) {
  ok <- tryCatch({ png(file, width=w, height=h, units="in", res=300, type="cairo"); TRUE },
                 error=function(e) FALSE)
  if (!ok) png(file, width=w, height=h, units="in", res=300)
  print(plot); dev.off()
}

# ---- loci (name, files, LD key, region, index rsID, PP.H4, sex, kidney, kind) ----
loci <- list(
  RAB19 = list(title="RAB19", pd="RAB19__PD", kid="RAB19__KID", ld="RAB19",
               chr=7, center=140117809, rsid="rs12216582", pph4=0.96, sex="both", kidney="eGFR", kind="shared"),
  BIN3_eGFR = list(title="BIN3 × eGFR", pd="BIN3_eGFR__PD", kid="BIN3_eGFR__KID", ld="BIN3",
               chr=8, center=22501830, rsid="rs13264187", pph4=0.87, sex="both", kidney="eGFR", kind="shared"),
  BIN3_eGFRcys = list(title="BIN3 × eGFRcys", pd="BIN3_eGFRcys__PD", kid="BIN3_eGFRcys__KID", ld="BIN3",
               chr=8, center=22501830, rsid="rs13264187", pph4=0.87, sex="both", kidney="eGFRcys", kind="shared"),
  SCARB2 = list(title="SCARB2", pd="SCARB2_distinct__PD", kid="SCARB2_distinct__KID", ld="SCARB2",
               chr=4, center=77108306, rsid="rs4859429", pph4=NA, sex="both", kidney="eGFR", kind="distinct"),
  MSRA_PRSS51 = list(title="MSRA / PRSS51", pd="MSRA_PRSS51__PD", kid="MSRA_PRSS51__KID", ld="MSRA",
               chr=8, center=10264000, rsid="rs6982308", pph4=0.93, sex="both", kidney="eGFR", kind="shared"),
  TNK2 = list(title="TNK2", pd="TNK2__PD", kid="TNK2__KID", ld="TNK2",
               chr=3, center=195624393, rsid="rs13059257", pph4=0.83, sex="female", kidney="eGFR", kind="shared")
)

sex_symbol <- function(s) switch(s, both="♂♀", female="♀", male="♂", "")

read_ss <- function(path) {
  df <- read.delim(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  df[!is.na(df$P) & df$P > 0, ]
}
read_ld <- function(path) {          # PLINK --r2 --ld-snp: cols CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
  if (!file.exists(path)) { warning("LD file not found: ", path, " - plotting grey"); return(NULL) }
  ld <- read.table(path, header=TRUE, stringsAsFactors=FALSE)
  data.frame(rsid = ld$SNP_B, r2 = ld$R2, stringsAsFactors=FALSE)
}
make_locus <- function(ss, ld_df, chr, center, window=5e5) {
  if (!is.null(ld_df)) { ss$ld <- ld_df$r2[match(ss$SNP, ld_df$rsid)]; ss$ld[is.na(ss$ld)] <- 0 }
  locus(data=ss, seqname=chr, xrange=c(center-window, center+window),
        ens_db="EnsDb.Hsapiens.v75", chrom="CHR", pos="BP", p="P", labs="SNP",
        LD = if (!is.null(ld_df)) "ld" else NULL)
}

build_stack <- function(L) {
  pd  <- read_ss(file.path(WIN_DIR, paste0(L$pd,  ".tsv")))
  kid <- read_ss(file.path(WIN_DIR, paste0(L$kid, ".tsv")))
  ld  <- read_ld(file.path(LD_DIR, paste0(L$ld,  ".ld")))
  loc_pd  <- make_locus(pd,  ld, L$chr, L$center)
  loc_kid <- make_locus(kid, ld, L$chr, L$center)
  sym <- sex_symbol(L$sex)

  # ---- print the index SNP (to add as a label BY HAND) + each panel's peak ----
  ipos <- pd$BP[match(L$rsid, pd$SNP)]; if (is.na(ipos)) ipos <- kid$BP[match(L$rsid, kid$SNP)]
  tpd <- pd[which.min(pd$P), ]; tkd <- kid[which.min(kid$P), ]
  cat(sprintf("\n=== %-13s  INDEX SNP: %s  (chr%d:%s)  %s ===\n", L$title, L$rsid, L$chr,
              ifelse(is.na(ipos), "?", format(ipos, scientific=FALSE)),
              if (L$kind=="shared") sprintf("PP.H4 = %.2f", L$pph4) else "distinct (control)"))
  cat(sprintf("    >> LABEL BY HAND:  %s\n", L$rsid))
  cat(sprintf("    PD    peak: %-12s  -log10P = %5.2f  (chr%d:%d)\n", tpd$SNP, -log10(tpd$P), L$chr, tpd$BP))
  cat(sprintf("    %-5s peak: %-12s  -log10P = %5.2f  (chr%d:%d)\n", L$kidney, tkd$SNP, -log10(tkd$P), L$chr, tkd$BP))

  # ---- panels: trait name TOP-LEFT; NO title banner; index SNP drawn as the
  #      purple diamond but NOT auto-labelled (you add the rsID text yourself) ----
  p_pd <- gg_scatter(loc_pd, index_snp=L$rsid, legend_pos=NULL) +
    annotate("text", x=-Inf, y=Inf, label=sprintf("PD %s", sym), hjust=-0.10, vjust=1.4, size=4, fontface="bold") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  p_kid <- gg_scatter(loc_kid, index_snp=L$rsid, legend_pos=NULL) +
    annotate("text", x=-Inf, y=Inf, label=sprintf("%s %s", L$kidney, sym), hjust=-0.10, vjust=1.4, size=4, fontface="bold") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  p_gene <- gg_genetracks(loc_pd, maxrows=4, gene_col="grey40", cex.text=0.5)

  plot_grid(p_pd, p_kid, p_gene, ncol=1, align="v", axis="lr", rel_heights=c(1, 1, 0.6))
}

# ---- build every panel + save individually ----
stacks <- list()
for (nm in names(loci)) {
  stacks[[nm]] <- build_stack(loci[[nm]])
  save_pdf(file.path(OUT_DIR, sprintf("Fig_%s.pdf", nm)), stacks[[nm]], 6, 7.5)
}

# ---- hand-built horizontal r^2 legend (same style as the manuscript) ----
legend_items <- data.frame(
  label = c("Missing","0.0 - 0.2","0.2 - 0.4","0.4 - 0.6","0.6 - 0.8","0.8 - 1.0","Index SNP"),
  color = c("#BEBEBE","#0000FF","#00FFFF","#00CC00","#FFA500","#FF0000","#9933FF"), stringsAsFactors=FALSE)
legend_items$x <- seq_len(nrow(legend_items)); legend_items$label <- factor(legend_items$label, levels=legend_items$label)
legend_plot <- ggplot(legend_items, aes(x=x, y=1)) +
  geom_point(aes(fill=label), shape=21, size=3, color="black", stroke=0.3, show.legend=FALSE) +
  scale_fill_manual(values=setNames(legend_items$color, legend_items$label)) +
  geom_text(aes(label=label), hjust=0, nudge_x=0.1, size=2.3) +
  annotate("text", x=0.45, y=1, label=expression(r^2), size=2.8) +
  scale_x_continuous(limits=c(0.3, nrow(legend_items)+1.0), expand=c(0,0)) +
  ylim(0.9,1.1) + coord_fixed(ratio=0.25) + theme_void() + theme(plot.margin=margin(2,2,2,2))

# ---- MAIN figure: 2x2 (shared, cross-biomarker, distinct control) ----
main_grid <- plot_grid(stacks[["RAB19"]], stacks[["BIN3_eGFR"]],
                       stacks[["MSRA_PRSS51"]], stacks[["TNK2"]],
                       ncol=2, labels=c("a","b","c","d"), label_size=16)
main_fig <- plot_grid(main_grid, legend_plot, ncol=1, rel_heights=c(1, 0.03))
save_pdf(file.path(OUT_DIR, "Figure_coloc_main.pdf"), main_fig, 12, 15.5)
save_png(file.path(OUT_DIR, "Figure_coloc_main.png"), main_fig, 12, 15.5)

# ---- SUPP figure: 1x2 — the two rebuttal-support panels ----
# a = BIN3 x eGFRcys (cross-biomarker replication); b = SCARB2 (distinct control)
supp_grid <- plot_grid(stacks[["BIN3_eGFRcys"]], stacks[["SCARB2"]], ncol=2, labels=c("a","b"), label_size=16)
supp_fig <- plot_grid(supp_grid, legend_plot, ncol=1, rel_heights=c(1, 0.05))
save_pdf(file.path(OUT_DIR, "Figure_coloc_supp.pdf"), supp_fig, 12, 8)
save_png(file.path(OUT_DIR, "Figure_coloc_supp.png"), supp_fig, 12, 8)

cat(sprintf("\nDone -> %s\n  Figure_coloc_main.{pdf,png}  (a RAB19, b BIN3, c MSRA/PRSS51, d TNK2 - all PD x eGFR)\n  Figure_coloc_supp.{pdf,png}  (a BIN3 x eGFRcys, b SCARB2 distinct control)\n", OUT_DIR))