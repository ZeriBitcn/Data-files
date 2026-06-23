###############################################################################
# 01_diversity_ordination.R
#
# Diversity analysis and ordination for:
#   Asrat et al. "Distribution and Regeneration Potential of the Soil Seed Bank
#   of Ericaceous Vegetation in the Bale Mountains"
#   Journal of Vegetation Science
#
# This script performs:
#   1. Alpha diversity (species richness S, Shannon H', Pielou's evenness J')
#   2. Sørensen similarity between standing vegetation and seed bank
#   3. Non-metric multidimensional scaling (NMDS) with Bray-Curtis dissimilarity
#   4. Canonical correspondence analysis (CCA) with environmental variables
#   5. PERMANOVA and ANOSIM tests for compositional differences among PES
#
# Required R packages: vegan, dplyr, tidyr, readxl
# R version: 4.3.0+
###############################################################################

## ── Load packages ──────────────────────────────────────────────────────────
library(vegan)
library(dplyr)
library(tidyr)
library(readr)

set.seed(42)

## ── Load data ──────────────────────────────────────────────────────────────
data_dir <- "data"

# Standing vegetation (aboveground) — plot × species abundance
veg_raw <- read_csv(file.path(data_dir, "Aboveground(stand)_3_column_dat.csv"))
# Soil seed bank (germination) — plot × species abundance
seed_raw <- read_csv(file.path(data_dir, "Soil-Seedbank_3_column_data.csv"))
# Environmental data
env <- read_csv(file.path(data_dir, "Environmental_data.csv"))
# Plot metadata (elevation, site)
plots <- read_csv(file.path(data_dir, "Plot_data.csv"))

# Clean column names
names(veg_raw) <- c("row", "plot", "species", "abundance")
names(seed_raw) <- c("row", "plot", "species", "abundance")

# Rename environmental columns
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

# Rename plot data columns
plots <- plots %>%
  rename(plot = Plot,
         elevation_plot = Elevetion,
         site = Sites,
         latitude = `Latitude_(N)`,
         longitude = `Longitude_(E)`,
         locality = `Locality/relife`)

# Merge env with plot site info
env_merged <- env %>%
  left_join(plots %>% select(plot, site), by = "plot")

# Define PES (principal experimental site) from site column
# South = PES-S, North = PES-N, Plateau = PES-P (inferred from locality/elevation)
env_merged <- env_merged %>%
  mutate(PES = case_when(
    grepl("South", site, ignore.case = TRUE) ~ "PES-S",
    grepl("North", site, ignore.case = TRUE) ~ "PES-N",
    grepl("Plateau|Sanetti", site, ignore.case = TRUE) ~ "PES-P",
    TRUE ~ NA_character_
  ))

cat("Environmental data loaded:", nrow(env_merged), "plots\n")
cat("PES distribution:\n")
print(table(env_merged$PES, useNA = "ifany"))

## ── Build community matrices ────────────────────────────────────────────────

# Standing vegetation: plot × species matrix
veg_wide <- veg_raw %>%
  select(plot, species, abundance) %>%
  pivot_wider(names_from = species, values_from = abundance, values_fill = 0) %>%
  arrange(plot)
veg_matrix <- as.data.frame(veg_wide[, -1])
rownames(veg_matrix) <- veg_wide$plot

# Seed bank: plot × species matrix
seed_wide <- seed_raw %>%
  select(plot, species, abundance) %>%
  pivot_wider(names_from = species, values_from = abundance, values_fill = 0) %>%
  arrange(plot)
seed_matrix <- as.data.frame(seed_wide[, -1])
rownames(seed_matrix) <- seed_wide$plot

cat("\nVegetation matrix:", nrow(veg_matrix), "plots ×", ncol(veg_matrix), "species\n")
cat("Seed bank matrix:", nrow(seed_matrix), "plots ×", ncol(seed_matrix), "species\n")

## ── 1. Alpha diversity ──────────────────────────────────────────────────────
cat("\n=== ALPHA DIVERSITY ===\n")

# Standing vegetation
veg_richness <- specnumber(veg_matrix)
veg_shannon <- diversity(veg_matrix, index = "shannon")
veg_evenness <- veg_shannon / log(veg_richness)

# Seed bank
seed_richness <- specnumber(seed_matrix)
seed_shannon <- diversity(seed_matrix, index = "shannon")
seed_evenness <- seed_shannon / log(seed_richness)

alpha_diversity <- data.frame(
  plot = as.integer(rownames(veg_matrix)),
  veg_richness = veg_richness,
  veg_shannon = veg_shannon,
  veg_evenness = veg_evenness,
  seed_richness = seed_richness[match(rownames(veg_matrix), rownames(seed_matrix))],
  seed_shannon = seed_shannon[match(rownames(veg_matrix), rownames(seed_matrix))],
  seed_evenness = seed_evenness[match(rownames(veg_matrix), rownames(seed_matrix))]
)

# Add PES info
alpha_diversity <- alpha_diversity %>%
  left_join(env_merged %>% select(plot, PES, elevation, disturbance), by = "plot")

# Summary by PES
cat("\nAlpha diversity by PES (mean ± SD):\n")
alpha_summary <- alpha_diversity %>%
  group_by(PES) %>%
  summarise(
    veg_richness_mean = mean(veg_richness, na.rm = TRUE),
    veg_richness_sd = sd(veg_richness, na.rm = TRUE),
    veg_shannon_mean = mean(veg_shannon, na.rm = TRUE),
    veg_shannon_sd = sd(veg_shannon, na.rm = TRUE),
    veg_evenness_mean = mean(veg_evenness, na.rm = TRUE),
    seed_richness_mean = mean(seed_richness, na.rm = TRUE),
    seed_shannon_mean = mean(seed_shannon, na.rm = TRUE),
    n = n()
  )
print(alpha_summary)

# Welch's t-tests (standing vegetation vs seed bank)
t_rich <- t.test(alpha_diversity$veg_richness, alpha_diversity$seed_richness)
t_shan <- t.test(alpha_diversity$veg_shannon, alpha_diversity$seed_shannon)
cat("\nWelch's t-test — Richness: t =", round(t_rich$statistic, 2), ", p =", t_rich$p.value, "\n")
cat("Welch's t-test — Shannon: t =", round(t_shan$statistic, 2), ", p =", t_shan$p.value, "\n")

## ── 2. Sørensen similarity ──────────────────────────────────────────────────
cat("\n=== SØRENSEN SIMILARITY ===\n")

# Species lists
veg_spp <- colnames(veg_matrix)
seed_spp <- colnames(seed_matrix)
shared <- intersect(veg_spp, seed_spp)
sorensen <- 2 * length(shared) / (length(veg_spp) + length(seed_spp))
cat("Standing vegetation species:", length(veg_spp), "\n")
cat("Seed bank species:", length(seed_spp), "\n")
cat("Shared species:", length(shared), "\n")
cat("Sørensen similarity:", round(sorensen, 3), "\n")

# Per-plot Sørensen
sorensen_plots <- sapply(1:nrow(veg_matrix), function(i) {
  v_pres <- names(which(veg_matrix[i, ] > 0))
  s_pres <- names(which(seed_matrix[match(rownames(veg_matrix)[i], rownames(seed_matrix)), ] > 0))
  c_shared <- length(intersect(v_pres, s_pres))
  if (length(v_pres) + length(s_pres) == 0) return(NA)
  2 * c_shared / (length(v_pres) + length(s_pres))
})
cat("Mean per-plot Sørensen:", round(mean(sorensen_plots, na.rm = TRUE), 3),
    "±", round(sd(sorensen_plots, na.rm = TRUE), 3), "\n")

## ── 3. NMDS ordination ──────────────────────────────────────────────────────
cat("\n=== NMDS ORDINATION ===\n")

# NMDS on standing vegetation (Bray-Curtis)
nmds <- metaMDS(veg_matrix, distance = "bray", k = 2, trymax = 100, trace = 0)
cat("NMDS stress:", round(nmds$stress, 3), "\n")

# ANOSIM test among PES
anosim_result <- anosim(veg_matrix, env_merged$PES, distance = "bray")
cat("ANOSIM: R =", round(anosim_result$statistic, 3), ", p =", anosim_result$signif, "\n")

## ── 4. CCA ordination ───────────────────────────────────────────────────────
cat("\n=== CCA ORDINATION ===\n")

# CCA with environmental variables
env_vars <- env_merged %>%
  select(elevation, slope, rock, moss, litter, disturbance) %>%
  mutate_all(as.numeric)

cca_result <- cca(veg_matrix ~ ., data = env_vars)
cat("CCA summary:\n")
print(summary(cca_result))
cat("\nCCA eigenvalues (first 2 axes):", round(cca_result$CCA$eig[1:2], 4), "\n")
cat("Constrained inertia proportion:", round(cca_result$CCA$tot.chi / cca_result$tot.chi, 3), "\n")

## ── 5. PERMANOVA ────────────────────────────────────────────────────────────
cat("\n=== PERMANOVA ===\n")

# PERMANOVA on Sørensen distances among PES
permanova_result <- adonis2(veg_matrix ~ PES, data = env_merged,
                            method = "bray", permutations = 999)
print(permanova_result)

## ── Save results ────────────────────────────────────────────────────────────
write_csv(alpha_diversity, "output/alpha_diversity.csv")
write_csv(alpha_summary, "output/alpha_diversity_summary_by_PES.csv")

cat("\n=== DIVERSITY AND ORDINATION ANALYSIS COMPLETE ===\n")
cat("Results saved to output/ directory\n")
