###############################################################################
# 02_glmm_models.R
#
# Generalized Linear Mixed Models (GLMMs) for:
#   Asrat et al. "Distribution and Regeneration Potential of the Soil Seed Bank
#   of Ericaceous Vegetation in the Bale Mountains"
#   Journal of Vegetation Science
#
# This script performs:
#   1. GLMM for Erica seed density (negative binomial)
#   2. GLMM for Erica germination potential (negative binomial)
#   3. Model selection via AICc
#   4. Post-hoc power analysis for fire-history / disturbance effects
#
# Required R packages: glmmTMB, lme4, MuMIn, DHARMa, dplyr, readr
# R version: 4.3.0+
###############################################################################

## ── Load packages ──────────────────────────────────────────────────────────
library(glmmTMB)
library(lme4)
library(MuMIn)      # for AICc
library(DHARMa)     # for residual diagnostics
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
         elevation = Altitude,
         slope = Slop,
         soil = Soil,
         rock = Rock,
         deadwood = DW,
         litter = Litter,
         erica_tree = Ericatree,
         erica_seedling = `Erika_kindergarten`,
         shrub_cover = Shrubcov,
         moss = Bryo,
         lichen = Lich,
         disturbance = Disturbance)

germ <- germ %>%
  rename(plot = Plot,
         moss_sum = `Moss_sum`,
         lichen_cover = `Lycon_cover`,
         germ_erica = `Germ_Erica`,
         germ_total = `Ger_cout`)

# Merge
df <- env %>%
  left_join(germ, by = "plot") %>%
  mutate(
    erica_cover = erica_tree + erica_seedling,
    PES = case_when(
      grepl("South", as.character(plot)) ~ "PES-S",  # placeholder; replace with actual PES mapping
      TRUE ~ "PES-other"
    ),
    # Scale elevation for model convergence
    elevation_scaled = scale(elevation)[, 1]
  )

# Assign PES based on elevation ranges from the manuscript:
# PES-N (northeast): ~3400-3850 m
# PES-S (southwest): ~3000-3800 m
# PES-P (plateau):   ~3800-4300 m
# Since we don't have explicit PES in the env data, use elevation as proxy
# (In the actual study, PES was assigned in the field)
df <- df %>%
  mutate(PES = case_when(
    elevation >= 3800 ~ "PES-P",
    elevation < 3700 ~ "PES-S",
    TRUE ~ "PES-N"
  ))

cat("Data loaded:", nrow(df), "plots\n")
cat("PES distribution:\n")
print(table(df$PES))

## ── 1. GLMM for Erica seed density ─────────────────────────────────────────
cat("\n=== GLMM: ERICA SEED DENSITY ===\n")

# Note: In this dataset, erica_tree and erica_seedling are cover values,
# not seed counts. The seed density data comes from the physical extraction
# (not available as a separate column here). We model erica_cover as a proxy
# for abundance, and germ_erica as germination potential.

# Model: Erica cover ~ elevation + disturbance + stoniness + moss + (1|PES)
# Using negative binomial for overdispersed count-like data
m_seed <- glmmTMB(
  erica_cover ~ elevation_scaled + disturbance + rock + moss +
    (1 | PES),
  family = nbinom2(),
  data = df
)

cat("GLMM for Erica cover (proxy for seed density/abundance):\n")
print(summary(m_seed))

# AICc
cat("\nAICc:", round(AICc(m_seed), 2), "\n")

# Residual diagnostics
sim_res <- simulateResiduals(m_seed)
# plot(sim_res)  # uncomment to view diagnostic plots

## ── 2. GLMM for Erica germination ──────────────────────────────────────────
cat("\n=== GLMM: ERICA GERMINATION POTENTIAL ===\n")

# Model: Germ_Erica ~ elevation + disturbance + stoniness + moss + (1|PES)
m_germ <- glmmTMB(
  germ_erica ~ elevation_scaled + disturbance + rock + moss +
    (1 | PES),
  family = nbinom2(),
  data = df
)

cat("GLMM for Erica germination potential:\n")
print(summary(m_germ))

cat("\nAICc:", round(AICc(m_germ), 2), "\n")

## ── 3. Model selection via AICc ────────────────────────────────────────────
cat("\n=== MODEL SELECTION ===\n")

# Compare candidate models for Erica cover
m_full <- glmmTMB(erica_cover ~ elevation_scaled + disturbance + rock + moss +
                    (1 | PES), family = nbinom2(), data = df)
m_no_elev <- glmmTMB(erica_cover ~ disturbance + rock + moss +
                       (1 | PES), family = nbinom2(), data = df)
m_no_dist <- glmmTMB(erica_cover ~ elevation_scaled + rock + moss +
                       (1 | PES), family = nbinom2(), data = df)
m_no_moss <- glmmTMB(erica_cover ~ elevation_scaled + disturbance + rock +
                       (1 | PES), family = nbinom2(), data = df)
m_null <- glmmTMB(erica_cover ~ 1 + (1 | PES), family = nbinom2(), data = df)

model_comparison <- AICc(m_full, m_no_elev, m_no_dist, m_no_moss, m_null)
model_comparison$delta <- model_comparison$AICc - min(model_comparison$AICc)
model_comparison <- model_comparison[order(model_comparison$delta), ]
cat("\nModel comparison (Erica cover):\n")
print(model_comparison)

# Compare candidate models for Erica germination
m_germ_full <- glmmTMB(germ_erica ~ elevation_scaled + disturbance + rock + moss +
                         (1 | PES), family = nbinom2(), data = df)
m_germ_no_elev <- glmmTMB(germ_erica ~ disturbance + rock + moss +
                            (1 | PES), family = nbinom2(), data = df)
m_germ_no_dist <- glmmTMB(germ_erica ~ elevation_scaled + rock + moss +
                            (1 | PES), family = nbinom2(), data = df)
m_germ_null <- glmmTMB(germ_erica ~ 1 + (1 | PES), family = nbinom2(), data = df)

model_comparison_germ <- AICc(m_germ_full, m_germ_no_elev, m_germ_no_dist, m_germ_null)
model_comparison_germ$delta <- model_comparison_germ$AICc - min(model_comparison_germ$AICc)
model_comparison_germ <- model_comparison_germ[order(model_comparison_germ$delta), ]
cat("\nModel comparison (Erica germination):\n")
print(model_comparison_germ)

## ── 4. Polynomial regression comparison (linear vs quadratic) ──────────────
cat("\n=== POLYNOMIAL REGRESSION (elevation) ===\n")

# Linear vs quadratic for Erica cover vs elevation
m_linear <- lm(erica_cover ~ elevation, data = df)
m_quad <- lm(erica_cover ~ elevation + I(elevation^2), data = df)
cat("Erica cover ~ elevation:\n")
cat("  Linear AICc:", round(AICc(m_linear), 2), "\n")
cat("  Quadratic AICc:", round(AICc(m_quad), 2), "\n")
cat("  Delta AICc:", round(AICc(m_quad) - AICc(m_linear), 2), "\n")
if (AICc(m_quad) - AICc(m_linear) < -4) {
  cat("  → Quadratic term retained (ΔAICc < -4)\n")
} else {
  cat("  → Linear model retained (ΔAICc ≥ -4)\n")
}

## ── 5. Post-hoc power analysis ─────────────────────────────────────────────
cat("\n=== POST-HOC POWER ANALYSIS ===\n")

# ANOVA for disturbance effect on Erica cover
aov_erica <- aov(erica_cover ~ as.factor(disturbance), data = df)
aov_summary <- summary(aov_erica)[[1]]
ss_effect <- aov_summary$`Sum Sq`[1]
ss_total <- sum(aov_summary$`Sum Sq`)
eta_sq <- ss_effect / ss_total
cat("ANOVA — Erica cover vs disturbance:\n")
cat("  F =", round(aov_summary$`F value`[1], 2), "\n")
cat("  p =", round(aov_summary$`Pr(>F)`[1], 4), "\n")
cat("  η² =", round(eta_sq, 4), "\n")

# Power calculation using Cohen's f
cohen_f <- sqrt(eta_sq / (1 - eta_sq))
# Using pwr.anova.test
if (requireNamespace("pwr", quietly = TRUE)) {
  library(pwr)
  power_result <- pwr.anova.test(k = length(unique(df$disturbance)),
                                  n = nrow(df) / length(unique(df$disturbance)),
                                  f = cohen_f,
                                  sig.level = 0.05)
  cat("  Post-hoc power:", round(power_result$power, 3), "\n")
} else {
  cat("  (Install 'pwr' package for power calculation)\n")
}

## ── Save model results ─────────────────────────────────────────────────────
# Extract coefficients
coef_table <- data.frame(
  term = names(fixef(m_full)$cond),
  estimate = fixef(m_full)$cond,
  std_error = summary(m_full)$coefficients$cond[, "Std. Error"],
  z_value = summary(m_full)$coefficients$cond[, "z value"],
  p_value = summary(m_full)$coefficients$cond[, "Pr(>|z|)"]
)
write_csv(coef_table, "output/glmm_erica_cover_coefficients.csv")
write_csv(model_comparison, "output/model_comparison_erica_cover.csv")
write_csv(model_comparison_germ, "output/model_comparison_erica_germination.csv")

cat("\n=== GLMM ANALYSIS COMPLETE ===\n")
cat("Results saved to output/ directory\n")
