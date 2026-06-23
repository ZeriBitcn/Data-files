###############################################################################
# 03_anova_posthoc.R
#
# ANOVA and post-hoc tests for:
#   Asrat et al. "Distribution and Regeneration Potential of the Soil Seed Bank
#   of Ericaceous Vegetation in the Bale Mountains"
#   Journal of Vegetation Science
#
# This script performs:
#   1. One-way ANOVA comparing means across PES for all response variables
#   2. Tukey's HSD post-hoc pairwise comparisons
#   3. Welch's t-tests (standing vegetation vs seed bank diversity)
#   4. Shapiro-Wilk normality tests
#   5. Regression analyses for Figure 4 (Erica cover, seedlings, germination
#      vs. moss cover and disturbance intensity)
#
# Required R packages: dplyr, readr
# R version: 4.3.0+
###############################################################################

## ── Load packages ──────────────────────────────────────────────────────────
library(dplyr)
library(readr)

set.seed(42)

## ── Load data ──────────────────────────────────────────────────────────────
data_dir <- "data"

env <- read_csv(file.path(data_dir, "Environmental_data.csv"))
germ <- read_csv(file.path(data_dir, "SOIL_seed-germ_dat.csv"))

# Clean column names
env <- env %>%
  rename(plot = `Unnamed:_0`,
         elevation = Altitude, slope = Slop, soil = Soil, rock = Rock,
         deadwood = DW, litter = Litter,
         erica_tree = Ericatree, erica_seedling = `Erika_kindergarten`,
         shrub_cover = Shrubcov, moss = Bryo, lichen = Lich,
         disturbance = Disturbance)

germ <- germ %>%
  rename(plot = Plot, moss_sum = `Moss_sum`, lichen_cover = `Lycon_cover`,
         germ_erica = `Germ_Erica`, germ_total = `Ger_cout`)

# Merge
df <- env %>%
  left_join(germ, by = "plot") %>%
  mutate(erica_cover = erica_tree + erica_seedling,
         PES = case_when(
           elevation >= 3800 ~ "PES-P",
           elevation < 3700 ~ "PES-S",
           TRUE ~ "PES-N"
         ),
         PES = factor(PES, levels = c("PES-N", "PES-S", "PES-P")))

cat("Data loaded:", nrow(df), "plots\n")

## ── 1. Shapiro-Wilk normality tests ────────────────────────────────────────
cat("\n=== SHAPIRO-WILK NORMALITY TESTS ===\n")

variables <- c("erica_cover", "erica_seedling", "germ_erica", "germ_total",
                "moss", "litter", "elevation", "rock", "disturbance")

for (v in variables) {
  x <- df[[v]]
  if (length(x[!is.na(x)]) >= 3) {
    sw <- shapiro.test(x)
    normal <- if (sw$p.value >= 0.05) "normal" else "NON-normal"
    cat(sprintf("  %-18s  W = %.4f  p = %8.4g  → %s\n", v, sw$statistic, sw$p.value, normal))
  }
}

## ── 2. One-way ANOVA across PES ────────────────────────────────────────────
cat("\n=== ONE-WAY ANOVA ACROSS PES ===\n")

anova_vars <- c("erica_cover", "erica_seedling", "germ_erica", "germ_total",
                "moss", "litter", "elevation")

for (v in anova_vars) {
  formula_str <- paste0(v, " ~ PES")
  m <- aov(as.formula(formula_str), data = df)
  s <- summary(m)[[1]]
  cat(sprintf("  %-18s  F = %6.2f  p = %8.4g  df = %d, %d\n",
              v, s$`F value`[1], s$`Pr(>F)`[1], s$Df[1], s$Df[2]))
}

## ── 3. Tukey HSD post-hoc tests ────────────────────────────────────────────
cat("\n=== TUKEY HSD POST-HOC TESTS ===\n")

for (v in anova_vars) {
  formula_str <- paste0(v, " ~ PES")
  m <- aov(as.formula(formula_str), data = df)
  tukey <- TukeyHSD(m)
  cat(sprintf("\n  %s:\n", v))
  for (i in 1:nrow(tukey[[1]])) {
    row <- tukey[[1]][i, ]
    sig <- if (row["p adj"] < 0.001) "***" else if (row["p adj"] < 0.01) "**" else if (row["p adj"] < 0.05) "*" else "n.s."
    cat(sprintf("    %-12s  diff = %8.3f  CI = [%.3f, %.3f]  p-adj = %.4f  %s\n",
                rownames(tukey[[1]])[i], row["diff"], row["lwr"], row["upr"], row["p adj"], sig))
  }
}

## ── 4. Welch's t-tests (standing veg vs seed bank) ─────────────────────────
cat("\n=== WELCH'S T-TESTS (standing vegetation vs seed bank) ===\n")
cat("(Note: Run after 01_diversity_ordination.R produces alpha_diversity.csv)\n")

alpha_file <- "output/alpha_diversity.csv"
if (file.exists(alpha_file)) {
  alpha <- read_csv(alpha_file)
  t_rich <- t.test(alpha$veg_richness, alpha$seed_richness)
  t_shan <- t.test(alpha$veg_shannon, alpha$seed_shannon)
  t_even <- t.test(alpha$veg_evenness, alpha$seed_evenness)
  cat(sprintf("  Richness:  t = %6.2f  p = %8.4g\n", t_rich$statistic, t_rich$p.value))
  cat(sprintf("  Shannon:   t = %6.2f  p = %8.4g\n", t_shan$statistic, t_shan$p.value))
  cat(sprintf("  Evenness:  t = %6.2f  p = %8.4g\n", t_even$statistic, t_even$p.value))
} else {
  cat("  (alpha_diversity.csv not found — run 01_diversity_ordination.R first)\n")
}

## ── 5. Regression analyses for Figure 4 ────────────────────────────────────
cat("\n=== REGRESSION ANALYSES (FIGURE 4) ===\n")

# Panel (a): Erica cover vs Moss cover
m_a <- lm(erica_cover ~ moss, data = df)
cat(sprintf("\n  (a) Erica cover vs Moss cover:\n"))
cat(sprintf("      R² = %.3f  p = %.4f  slope = %.4f  n = %d\n",
            summary(m_a)$r.squared, summary(m_a)$coef[2, 4],
            summary(m_a)$coef[2, 1], nrow(df)))

# Panel (b): Erica cover vs Disturbance
m_b <- lm(erica_cover ~ disturbance, data = df)
cat(sprintf("\n  (b) Erica cover vs Disturbance:\n"))
cat(sprintf("      R² = %.3f  p = %.4f  slope = %.4f  n = %d\n",
            summary(m_b)$r.squared, summary(m_b)$coef[2, 4],
            summary(m_b)$coef[2, 1], nrow(df)))

# Panel (c): Erica seedling vs Disturbance
m_c <- lm(erica_seedling ~ disturbance, data = df)
cat(sprintf("\n  (c) Erica seedling vs Disturbance:\n"))
cat(sprintf("      R² = %.3f  p = %.4f  slope = %.4f  n = %d\n",
            summary(m_c)$r.squared, summary(m_c)$coef[2, 4],
            summary(m_c)$coef[2, 1], nrow(df)))

# Panel (d): Erica germination vs Disturbance
m_d <- lm(germ_erica ~ disturbance, data = df)
cat(sprintf("\n  (d) Erica germination vs Disturbance:\n"))
cat(sprintf("      R² = %.3f  p = %.4f  slope = %.4f  n = %d\n",
            summary(m_d)$r.squared, summary(m_d)$coef[2, 4],
            summary(m_d)$coef[2, 1], nrow(df)))

## ── 6. Additional regression analyses (elevation patterns, Figure 5) ───────
cat("\n=== ELEVATION PATTERNS (FIGURE 5) ===\n")

# Cubic regression for Erica cover vs elevation
m_cubic <- lm(erica_cover ~ poly(elevation, 3), data = df)
cat(sprintf("\n  Erica cover ~ elevation (cubic):\n"))
cat(sprintf("      R² = %.3f  p = %.4f\n", summary(m_cubic)$r.squared,
            summary(m_cubic)$fstatistic[1] %>% pf(3, summary(m_cubic)$df[2], lower.tail = FALSE)))

# Cubic regression for germination vs elevation
m_cubic_g <- lm(germ_erica ~ poly(elevation, 3), data = df)
cat(sprintf("\n  Erica germination ~ elevation (cubic):\n"))
cat(sprintf("      R² = %.3f  p = %.4f\n", summary(m_cubic_g)$r.squared,
            summary(m_cubic_g)$fstatistic[1] %>% pf(3, summary(m_cubic_g)$df[2], lower.tail = FALSE)))

## ── 7. Seed density vs germination relationship (Figure 6) ─────────────────
cat("\n=== SEED DENSITY vs GERMINATION (FIGURE 6) ===\n")

# Note: In this dataset, erica_tree is a cover proxy, not a direct seed count
# The R² = 0.04 relationship reported in the manuscript refers to the
# physical seed extraction counts vs greenhouse germination
m_seed_germ <- lm(germ_erica ~ erica_cover, data = df)
cat(sprintf("\n  Erica germination vs Erica cover (proxy for seed density):\n"))
cat(sprintf("      R² = %.3f  p = %.4f  n = %d\n",
            summary(m_seed_germ)$r.squared, summary(m_seed_germ)$coef[2, 4], nrow(df)))

## ── 8. Disturbance effects on vegetation attributes (Table 4) ──────────────
cat("\n=== DISTURBANCE EFFECTS ON VEGETATION (TABLE 4) ===\n")

dist_vars <- c("erica_cover", "moss", "litter", "erica_seedling", "shrub_cover")
for (v in dist_vars) {
  formula_str <- paste0(v, " ~ as.factor(disturbance)")
  m <- aov(as.formula(formula_str), data = df)
  s <- summary(m)[[1]]
  cat(sprintf("  %-18s  F = %6.2f  p = %8.4g\n", v, s$`F value`[1], s$`Pr(>F)`[1]))
}

## ── Save results ───────────────────────────────────────────────────────────
fig4_results <- data.frame(
  panel = c("(a)", "(b)", "(c)", "(d)"),
  x_variable = c("Moss cover", "Disturbance", "Disturbance", "Disturbance"),
  y_variable = c("Erica cover", "Erica cover", "Erica seedlings", "Erica germination"),
  R_squared = c(summary(m_a)$r.squared, summary(m_b)$r.squared,
                summary(m_c)$r.squared, summary(m_d)$r.squared),
  p_value = c(summary(m_a)$coef[2, 4], summary(m_b)$coef[2, 4],
              summary(m_c)$coef[2, 4], summary(m_d)$coef[2, 4]),
  slope = c(summary(m_a)$coef[2, 1], summary(m_b)$coef[2, 1],
            summary(m_c)$coef[2, 1], summary(m_d)$coef[2, 1]),
  n = nrow(df)
)
write_csv(fig4_results, "output/figure4_regression_results.csv")

cat("\n=== ANOVA AND POST-HOC ANALYSIS COMPLETE ===\n")
cat("Results saved to output/ directory\n")
