#!/usr/bin/env Rscript
###############################################################################
# plot_ldsc_forest_final.R
#
# Publication-ready LDSC forest plot for PD × kidney sex-concordant pairs
#
# Features:
#   - Three significance levels: not sig, p<0.05, Bonferroni (p<0.05/15)
#   - ♂/♀ symbols in y-axis labels
#   - Sex-combined labeled as such (not "meta")
#   - No title, clean minimal design
#   - 300 DPI, landscape
#
# Usage:
#   Rscript plot_ldsc_forest_final.R <rg_file> <output_dir>
#
# rg_file: tab-separated with columns p1, p2, rg, se, z, p
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(TRUE)
if (length(args) < 2) stop("Usage: Rscript plot_ldsc_forest_final.R <rg_file> <output_dir>")

rg_file    <- args[1]
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Read data
# ============================================================================
cat("Reading LDSC rg results...\n")
rg <- read.delim(rg_file, sep = "\t", header = TRUE, colClasses = "character")
rg$rg <- as.numeric(rg$rg)
rg$se <- as.numeric(rg$se)
rg$z  <- as.numeric(rg$z)
rg$p  <- as.numeric(rg$p)
rg$p[is.na(rg$p)] <- 0

# ============================================================================
# Define traits and sex-concordant pairs
# ============================================================================
pd_traits <- c("PD_metaSS_noProxy", "PD_female_noProxy", "PD_male_noProxy")

# Filter to PD × kidney pairs
all_traits <- unique(c(rg$p1, rg$p2))
kidney_traits <- setdiff(all_traits, pd_traits)

rg_pd <- rg %>%
  filter(
    (p1 %in% pd_traits & p2 %in% kidney_traits) |
    (p2 %in% pd_traits & p1 %in% kidney_traits)
  ) %>%
  mutate(
    pd_trait = ifelse(p1 %in% pd_traits, p1, p2),
    kidney_trait = ifelse(p1 %in% pd_traits, p2, p1)
  )

# Sex-concordant mapping
concordant <- data.frame(
  pd_trait = c(rep("PD_metaSS_noProxy", 5), rep("PD_female_noProxy", 5), rep("PD_male_noProxy", 5)),
  kidney_trait = c("eGFR_meta", "uACR_sexComb", "hematuria_sexComb", "uKCr_sexComb", "uNaCr_sexComb",
                   "eGFR_female", "uACR_female", "hematuria_female", "uKCr_female", "uNaCr_female",
                   "eGFR_male", "uACR_male", "hematuria_male", "uKCr_male", "uNaCr_male"),
  stringsAsFactors = FALSE
)

rg_concordant <- rg_pd %>%
  inner_join(concordant, by = c("pd_trait", "kidney_trait"))

cat("  Sex-concordant pairs:", nrow(rg_concordant), "\n")

# ============================================================================
# Bonferroni threshold: 0.05 / 15 sex-concordant pairs = 0.00333
# ============================================================================
n_concordant_tests <- 15
bonf_threshold <- 0.05 / n_concordant_tests
cat("  Bonferroni threshold: p <", bonf_threshold, "(0.05 /", n_concordant_tests, ")\n")

# ============================================================================
# Create labels with ♂/♀ symbols
# ============================================================================
# Kidney trait group (without sex suffix)
rg_concordant$kidney_group <- case_when(
  grepl("eGFR", rg_concordant$kidney_trait) ~ "eGFR",
  grepl("uACR", rg_concordant$kidney_trait) ~ "uACR",
  grepl("hematuria", rg_concordant$kidney_trait) ~ "Hematuria",
  grepl("uKCr", rg_concordant$kidney_trait) ~ "uK/Cr",
  grepl("uNaCr", rg_concordant$kidney_trait) ~ "uNa/Cr"
)

# Sex stratum
rg_concordant$sex_stratum <- case_when(
  rg_concordant$pd_trait == "PD_metaSS_noProxy" ~ "sex-combined",
  rg_concordant$pd_trait == "PD_female_noProxy" ~ "female",
  rg_concordant$pd_trait == "PD_male_noProxy" ~ "male"
)

# Y-axis label with symbols
rg_concordant$pair_label <- case_when(
  rg_concordant$sex_stratum == "sex-combined" ~ paste0("PD \u00d7 ", rg_concordant$kidney_group, " (\u2640\u2642)"),
  rg_concordant$sex_stratum == "female" ~ paste0("PD \u00d7 ", rg_concordant$kidney_group, " (\u2640)"),
  rg_concordant$sex_stratum == "male" ~ paste0("PD \u00d7 ", rg_concordant$kidney_group, " (\u2642)")
)

# Significance categories
rg_concordant$sig_level <- case_when(
  rg_concordant$p < bonf_threshold ~ "Bonferroni",
  rg_concordant$p < 0.05 ~ "Nominal",
  TRUE ~ "Not significant"
)

# CI
rg_concordant$ci_lo <- rg_concordant$rg - 1.96 * rg_concordant$se
rg_concordant$ci_hi <- rg_concordant$rg + 1.96 * rg_concordant$se

# ============================================================================
# Order: group by kidney trait, within each group order: sex-combined, ♀, ♂
# ============================================================================
kidney_order <- c("eGFR", "uACR", "Hematuria", "uK/Cr", "uNa/Cr")
sex_order <- c("sex-combined", "female", "male")

rg_concordant$kidney_group <- factor(rg_concordant$kidney_group, levels = kidney_order)
rg_concordant$sex_stratum <- factor(rg_concordant$sex_stratum, levels = sex_order)

rg_concordant <- rg_concordant %>% arrange(kidney_group, sex_stratum)
rg_concordant$pair_label <- factor(rg_concordant$pair_label,
                                    levels = rev(rg_concordant$pair_label))

# ============================================================================
# Color mapping by sex stratum
# ============================================================================
strata_colors <- c(
  "sex-combined" = "#B8860B",  # muted gold
  "female" = "#6B8E23",         # olive green
  "male" = "#8B5E3C"            # warm brown
)

# Shape mapping by significance
sig_shapes <- c(
  "Bonferroni" = 16,    # Filled circle
  "Nominal" = 17,       # Filled triangle
  "Not significant" = 1  # Open circle
)

# ============================================================================
# Plot
# ============================================================================
cat("Generating forest plot...\n")

p <- ggplot(rg_concordant, aes(x = rg, y = pair_label)) +
  # Zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  # Horizontal grouping lines (subtle)
  geom_hline(yintercept = seq(3.5, 12.5, 3), color = "grey90", linewidth = 0.3) +
  # CI bars
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, color = sex_stratum),
                 height = 0, linewidth = 0.6) +
  # Points
  geom_point(aes(color = sex_stratum, shape = sig_level), size = 3.5) +
  # Scales
  scale_color_manual(
    values = strata_colors,
    labels = c(
      "sex-combined" = "Sex-combined",
      "female" = "Female",
      "male" = "Male"
    ),
    name = "Sex stratum"
  ) +
  scale_shape_manual(
    values = sig_shapes,
    labels = c(
      "Bonferroni" = "Bonferroni",
      "Nominal" = "Nominal",
      "Not significant" = "Not sig."
    ),
    name = "Significance"
  ) +
#  scale_shape_manual(
#    values = sig_shapes,
#    labels = c(
#      "Bonferroni" = paste0("p < ", format(bonf_threshold, digits = 2, scientific = FALSE)),
#      "Nominal" = "p < 0.05",
#      "Not significant" = "Not significant"
#    ),
#    name = "Significance"
#  ) +
  # Axis
  labs(
    x = expression("Genetic correlation (" * r[g] * ")"),
    y = NULL
  ) +
  # Theme
  theme_minimal(base_size = 13) +
  theme(
    # No title
    plot.title = element_blank(),
    # Typography
    text = element_text(family = "sans"),
    axis.title.x = element_text(size = 13, margin = margin(t = 10)),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11),
    # Grid
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    # Legend
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.margin = margin(t = 5),
    # Margins
    plot.margin = margin(15, 20, 10, 10),
    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 4)),
    shape = guide_legend(order = 2, override.aes = list(size = 4))
  )

# ============================================================================
# Save at high resolution
# ============================================================================
# 4800 x 3000 pixels at 300 DPI = 16 x 10 inches
ggsave(file.path(output_dir, "ldsc_forest_concordant_final.pdf"), p,
       width = 10, height = 7, dpi = 300, device = cairo_pdf)
ggsave(file.path(output_dir, "ldsc_forest_concordant_final.png"), p,
       width = 10, height = 7, dpi = 300, type = "cairo")

cat("  Saved: ldsc_forest_concordant_final.pdf/png\n")

# ============================================================================
# Print significance summary
# ============================================================================
cat("\n=== Significance summary ===\n")
cat("Bonferroni threshold: p <", bonf_threshold, "\n")
cat("Bonferroni significant:\n")
bonf_pairs <- rg_concordant %>% filter(sig_level == "Bonferroni")
if (nrow(bonf_pairs) > 0) {
  print(bonf_pairs[, c("pair_label", "rg", "se", "p", "sig_level")], row.names = FALSE)
} else {
  cat("  None\n")
}
cat("\nNominally significant (p < 0.05):\n")
nom_pairs <- rg_concordant %>% filter(sig_level == "Nominal")
if (nrow(nom_pairs) > 0) {
  print(nom_pairs[, c("pair_label", "rg", "se", "p", "sig_level")], row.names = FALSE)
} else {
  cat("  None\n")
}

cat("\nDone!\n")
