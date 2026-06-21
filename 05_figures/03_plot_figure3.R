#!/usr/bin/env Rscript
###############################################################################
# 03_plot_figure3.R
#
# Generates Figure 3: paired regional association plots (PD on top, eGFR on
# bottom) for four colocalized loci, using locuszoomr v0.3.8.
#
# Each locus: two stacked panels sharing the x-axis, SNPs colored by LD r^2
# to the index (top shared causal) variant, with a gene track.
#
# Run on Rorqual:
#   module load StdEnv/2023 r/4.3.1
#   Rscript 03_plot_figure3.R
###############################################################################

suppressPackageStartupMessages({
  library(locuszoomr)
  library(EnsDb.Hsapiens.v75)
  library(ggplot2)
  library(cowplot)
})

# ============================================================================
# Paths
# ============================================================================
BASE     <- file.path(Sys.getenv("HOME"),
                      "links/projects/def-gsarah/lchang24/github/coloc-susie/figure3")
WIN_DIR  <- file.path(BASE, "windows")
LD_DIR   <- file.path(BASE, "ld")
OUT_DIR  <- file.path(BASE, "output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Locus definitions
#   name, chr, index position (hg19), index rsID, PP.H4, pair label, gene
# ============================================================================
loci <- list(
  list(name = "RAB19", chr = 7, pos = 140125361, rsid = "rs12216582",
       pph4 = 0.97, gene = "RAB19", sex = "both"),
  list(name = "MSRA",  chr = 8, pos = 10193772,  rsid = "rs6982308",
       pph4 = 0.90, gene = "MSRA", sex = "both"),
  list(name = "BIN3",  chr = 8, pos = 22519815,  rsid = "rs7005025",
       pph4 = 0.87, gene = "BIN3", sex = "both"),
  list(name = "TNK2",  chr = 3, pos = 195634993, rsid = "rs13059257",
       pph4 = 0.83, gene = "TNK2", sex = "female")
)

# Sex symbols: both = male+female combined; female = female-only
sex_symbol <- function(sex) {
  if (sex == "both")   return("\u2642\u2640")   # ♂♀
  if (sex == "female") return("\u2640")          # ♀
  if (sex == "male")   return("\u2642")          # ♂
  ""
}

# ============================================================================
# Helper: read a windowed sumstats file
# ============================================================================
read_sumstats <- function(path) {
  df <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  # Standard columns: variant_id chromosome base_pair_location effect_allele
  #                   other_allele beta standard_error p_value N
  df <- df[!is.na(df$p_value) & df$p_value > 0, ]
  df
}

# ============================================================================
# Helper: read PLINK .ld file and return data.frame(rsid, r2)
#   PLINK --r2 with --ld-snp output columns:
#   CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
# ============================================================================
read_ld <- function(path) {
  if (!file.exists(path)) {
    warning(sprintf("LD file not found: %s — plotting without r^2 coloring", path))
    return(NULL)
  }
  ld <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  # SNP_B is the partner SNP, R2 is r^2 to the index (SNP_A)
  data.frame(rsid = ld$SNP_B, r2 = ld$R2, stringsAsFactors = FALSE)
}

# ============================================================================
# Helper: build a locus object for one trait at one region
# ============================================================================
make_locus <- function(sumstats, ld_df, chr, index_pos, window = 5e5) {
  # Merge LD r^2 into sumstats by rsID
  if (!is.null(ld_df)) {
    sumstats$ld <- ld_df$r2[match(sumstats$variant_id, ld_df$rsid)]
    sumstats$ld[is.na(sumstats$ld)] <- 0
  }

  # locuszoomr::locus() builds the locus object.
  # We pass explicit column names matching our sumstats.
  loc <- locus(
    data       = sumstats,
    seqname    = chr,
    xrange     = c(index_pos - window, index_pos + window),
    ens_db     = "EnsDb.Hsapiens.v75",
    chrom      = "chromosome",
    pos        = "base_pair_location",
    p          = "p_value",
    labs       = "variant_id",
    LD         = if (!is.null(ld_df)) "ld" else NULL
  )
  loc
}

# ============================================================================
# Generate each locus panel-pair and save individually
# ============================================================================
panel_list <- list()

for (L in loci) {
  cat(sprintf("\n=== %s (chr%d:%d) ===\n", L$name, L$chr, L$pos))

  # File paths
  pd_file   <- file.path(WIN_DIR, sprintf("%s_PD.tsv", L$name))
  egfr_file <- file.path(WIN_DIR, sprintf("%s_eGFR.tsv", L$name))
  ld_file   <- file.path(LD_DIR, sprintf("%s.ld", L$name))

  # Read data
  pd_ss   <- read_sumstats(pd_file)
  egfr_ss <- read_sumstats(egfr_file)
  ld_df   <- read_ld(ld_file)

  cat(sprintf("  PD SNPs: %d | eGFR SNPs: %d | LD SNPs: %s\n",
              nrow(pd_ss), nrow(egfr_ss),
              if (is.null(ld_df)) "none" else nrow(ld_df)))

  # Build locus objects
  loc_pd   <- make_locus(pd_ss,   ld_df, L$chr, L$pos)
  loc_egfr <- make_locus(egfr_ss, ld_df, L$chr, L$pos)

  # Sex symbol for this locus
  sym <- sex_symbol(L$sex)
  pd_label   <- sprintf("PD %s", sym)
  egfr_label <- sprintf("eGFR %s", sym)

  # --- Plot PD panel (top): scatter, label inside, no r2 legend ---
  p_pd <- gg_scatter(
    loc_pd,
    index_snp = L$rsid,
    legend_pos = NULL
  ) +
    annotate("text", x = -Inf, y = Inf, label = pd_label,
             hjust = -0.15, vjust = 1.5, size = 4, fontface = "bold") +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank())

  # --- Plot eGFR panel (middle): scatter, label inside, no r2 legend ---
  p_egfr <- gg_scatter(
    loc_egfr,
    index_snp = L$rsid,
    legend_pos = NULL
  ) +
    annotate("text", x = -Inf, y = Inf, label = egfr_label,
             hjust = -0.15, vjust = 1.5, size = 4, fontface = "bold") +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank())

  # --- Gene track (bottom): fewer rows + smaller labels to avoid overlap ---
  p_genes <- gg_genetracks(loc_pd, maxrows = 6, gene_col = "grey40",
                           cex.text = 0.5)

  # Stack the three panels for this locus (PD top, eGFR middle, genes bottom)
  locus_stack <- plot_grid(
    p_pd, p_egfr, p_genes,
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1, 1, 0.5)
  )

  panel_list[[L$name]] <- locus_stack

  # Save individual locus plot
  out_single <- file.path(OUT_DIR, sprintf("Figure3_%s.pdf", L$name))
  ggsave(out_single, locus_stack, width = 6, height = 7, units = "in")
  cat(sprintf("  Saved: %s\n", out_single))
}

# ============================================================================
# Assemble final Figure 3 — 2 x 2 (each panel wide enough for gene labels)
# ============================================================================
cat("\n=== Assembling final Figure 3 (2 x 2) ===\n")

figure3_grid <- plot_grid(
  panel_list[["RAB19"]], panel_list[["BIN3"]],
  panel_list[["MSRA"]],  panel_list[["TNK2"]],
  ncol = 2,
  labels = c("a", "b", "c", "d"),
  label_size = 16
)

# ----------------------------------------------------------------------------
# Build a standalone horizontal r^2 legend with circular points,
# matching locuszoomr's style: Missing | 0.0-0.2 | ... | 0.8-1.0 | Index SNP
# ----------------------------------------------------------------------------
# locuszoomr default LD colors (blue -> cyan -> green -> orange -> red)
legend_items <- data.frame(
  label = c("Missing", "0.0 - 0.2", "0.2 - 0.4", "0.4 - 0.6",
            "0.6 - 0.8", "0.8 - 1.0", "Index SNP"),
  color = c("#BEBEBE", "#0000FF", "#00FFFF", "#00CC00",
            "#FFA500", "#FF0000", "#9933FF"),
  stringsAsFactors = FALSE
)
legend_items$x <- seq_len(nrow(legend_items))
legend_items$label <- factor(legend_items$label, levels = legend_items$label)

legend_plot <- ggplot(legend_items, aes(x = x, y = 1)) +
  geom_point(aes(fill = label), shape = 21, size = 3, color = "black",
             stroke = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = setNames(legend_items$color, legend_items$label)) +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.1, size = 2.3) +
  annotate("text", x = 0.45, y = 1, label = expression(r^2), size = 2.8) +
  scale_x_continuous(limits = c(0.3, nrow(legend_items) + 1.0),
                     expand = c(0, 0)) +
  ylim(0.9, 1.1) +
  coord_fixed(ratio = 0.25) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 2))

# Stack figure above the horizontal legend
figure3 <- plot_grid(
  figure3_grid,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.03)
)

# Save as PDF (vector, for publication) and PNG (raster, for preview)
out_pdf <- file.path(OUT_DIR, "Figure3_combined.pdf")
out_png <- file.path(OUT_DIR, "Figure3_combined.png")

# 2x2 canvas + legend strip
ggsave(out_pdf, figure3, width = 12, height = 14.8, units = "in")
ggsave(out_png, figure3, width = 12, height = 14.8, units = "in", dpi = 300)

cat(sprintf("\nFinal figure saved:\n  %s\n  %s\n", out_pdf, out_png))
cat("\nDone.\n")