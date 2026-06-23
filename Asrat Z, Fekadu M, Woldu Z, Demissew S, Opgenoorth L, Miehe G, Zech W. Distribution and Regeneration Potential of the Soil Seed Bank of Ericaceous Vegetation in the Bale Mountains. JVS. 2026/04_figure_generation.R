###############################################################################
# 04_figure_generation.R
#
# Figure generation for:
#   Asrat et al. "Distribution and Regeneration Potential of the Soil Seed Bank
#   of Ericaceous Vegetation in the Bale Mountains"
#   Journal of Vegetation Science
#
# This script generates:
#   Fig. 2 — Temporal emergence patterns (180-day germination trial)
#   Fig. 4 — Environmental drivers of Erica regeneration (4-panel: moss/disturbance)
#   Fig. 5 — Elevational patterns in diversity and Erica regeneration
#   Fig. 6 — Seed density vs germination and disturbance effects
#
# Required R packages: ggplot2, dplyr, readr, vegan, tidyr
# R version: 4.3.0+
###############################################################################

## ── Load packages ──────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(vegan)

set.seed(42)

## ── Load data ──────────────────────────────────────────────────────────────
data_dir <- "data"
output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE)

env <- read_csv(file.path(data_dir, "Environmental_data.csv"))
germ <- read_csv(file.path(data_dir, "SOIL_seed-germ_dat.csv"))
germ_counts <- read_csv(file.path(data_dir, "Germination_count_per_10_days_.csv"))

# Clean column names
env <- env %>%
  rename(plot = `Unnamed:_0`, elevation = Altitude, slope = Slop,
         soil = Soil, rock = Rock, deadwood = DW, litter = Litter,
         erica_tree = Ericatree, erica_seedling = `Erika_kindergarten`,
         shrub_cover = Shrubcov, moss = Bryo, lichen = Lich,
         disturbance = Disturbance)

germ <- germ %>%
  rename(plot = Plot, moss_sum = `Moss_sum`, lichen_cover = `Lycon_cover`,
         germ_erica = `Germ_Erica`, germ_total = `Ger_cout`)

df <- env %>%
  left_join(germ, by = "plot") %>%
  mutate(erica_cover = erica_tree + erica_seedling,
         PES = case_when(
           elevation >= 3800 ~ "PES-P",
           elevation < 3700 ~ "PES-S",
           TRUE ~ "PES-N"
         ),
         PES = factor(PES, levels = c("PES-N", "PES-S", "PES-P")),
         disturbance_label = paste0("D", disturbance))

# Colour palette
col_moss   <- "#2A9D8F"  # teal-green (facilitation)
col_dist   <- "#E63946"  # red (suppression)
col_deco   <- "#6C757D"  # gray (decoupling)
col_n      <- "#0077BB"  # blue (PES-N)
col_s      <- "#EE7733"  # orange (PES-S)
col_p      <- "#CC3311"  # red (PES-P)

## ── Figure 2: Temporal emergence patterns ───────────────────────────────────
cat("=== FIGURE 2: Temporal emergence patterns ===\n")

# Reshape germination counts to long format
germ_long <- germ_counts %>%
  select(-`Ger_cout`, -`Avarage`) %>%
  pivot_longer(cols = -plot, names_to = "interval", values_to = "count") %>%
  mutate(interval = as.numeric(interval),
         midpoint_days = interval * 10) %>%
  left_join(df %>% select(plot, PES, disturbance, disturbance_label), by = "plot")

p2 <- ggplot(germ_long, aes(x = midpoint_days, y = count, color = disturbance_label)) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.15, linewidth = 1.2) +
  facet_wrap(~ PES, ncol = 3) +
  scale_color_manual(values = c("D0" = "#2A9D8F", "D1" = "#93C5AE",
                                 "D2" = "#F4A261", "D3" = "#E76F51",
                                 "D4" = "#CC3311"),
                     name = "Disturbance") +
  labs(x = "Days since trial start", y = "Seedlings emerged (per 10-day interval)",
       title = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(output_dir, "Figure2_temporal_emergence.png"),
       p2, width = 12, height = 5, dpi = 600, bg = "white")
cat("  ✓ Figure 2 saved\n")

## ── Figure 4: Environmental drivers (4-panel) ───────────────────────────────
cat("\n=== FIGURE 4: Environmental drivers of Erica regeneration ===\n")

# Panel (a): Erica cover vs Moss cover
p4a <- ggplot(df, aes(x = moss, y = erica_cover)) +
  geom_point(color = col_moss, size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = col_moss, linewidth = 1.2,
              se = FALSE, linetype = ifelse(summary(lm(erica_cover ~ moss, df))$coef[2,4] < 0.05, "solid", "dashed")) +
  annotate("label", x = min(df$moss) + 5, y = max(df$erica_cover, na.rm = TRUE) * 0.95,
           label = paste0("R² = ", round(summary(lm(erica_cover ~ moss, df))$r.squared, 3),
                          "\np = ", formatC(summary(lm(erica_cover ~ moss, df))$coef[2,4], format = "f", digits = 3), " ***"),
           hjust = 0, size = 3.5, color = col_moss,
           fill = "white", label.size = 0.3) +
  labs(x = "Moss cover (%)", y = "Erica cover (%)", tag = "(a)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14))

# Panel (b): Erica cover vs Disturbance
p4b <- ggplot(df, aes(x = disturbance, y = erica_cover)) +
  geom_point(color = col_dist, size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = col_dist, linewidth = 1.2, se = FALSE) +
  annotate("label", x = min(df$disturbance) + 0.2, y = max(df$erica_cover, na.rm = TRUE) * 0.95,
           label = paste0("R² = ", round(summary(lm(erica_cover ~ disturbance, df))$r.squared, 3),
                          "\np < 0.001 ***"),
           hjust = 0, size = 3.5, color = col_dist,
           fill = "white", label.size = 0.3) +
  labs(x = "Disturbance intensity (D₀–D₄ scale)", y = "Erica cover (%)", tag = "(b)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14))

# Panel (c): Erica seedling vs Disturbance
p4c <- ggplot(df, aes(x = disturbance, y = erica_seedling)) +
  geom_point(color = col_dist, size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = col_dist, linewidth = 1.2, se = FALSE) +
  annotate("label", x = min(df$disturbance) + 0.2, y = max(df$erica_seedling, na.rm = TRUE) * 0.95,
           label = paste0("R² = ", round(summary(lm(erica_seedling ~ disturbance, df))$r.squared, 3),
                          "\np = ", formatC(summary(lm(erica_seedling ~ disturbance, df))$coef[2,4], format = "f", digits = 3), " *"),
           hjust = 0, size = 3.5, color = col_dist,
           fill = "white", label.size = 0.3) +
  labs(x = "Disturbance intensity (D₀–D₄ scale)", y = "Erica seedlings (count per plot)", tag = "(c)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14))

# Panel (d): Erica germination vs Disturbance
p4d <- ggplot(df, aes(x = disturbance, y = germ_erica)) +
  geom_point(color = col_deco, size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = col_deco, linewidth = 1.2, se = FALSE, linetype = "dashed") +
  annotate("label", x = min(df$disturbance) + 0.2, y = max(df$germ_erica, na.rm = TRUE) * 0.95,
           label = paste0("R² = ", round(summary(lm(germ_erica ~ disturbance, df))$r.squared, 3),
                          "\np = 0.84 n.s."),
           hjust = 0, size = 3.5, color = col_deco,
           fill = "white", label.size = 0.3) +
  labs(x = "Disturbance intensity (D₀–D₄ scale)", y = "Erica germination (seedlings m⁻²)", tag = "(d)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14))

# Combine using patchwork (or cowplot)
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p4 <- (p4a + p4b) / (p4c + p4d) +
    plot_annotation(
      title = "Environmental drivers of Erica regeneration",
      subtitle = "Facilitation by moss, suppression by disturbance, and decoupling of the seed bank",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  ggsave(file.path(output_dir, "Figure4_environmental_drivers.png"),
         p4, width = 14, height = 11, dpi = 600, bg = "white")
  cat("  ✓ Figure 4 saved (patchwork)\n")
} else {
  ggsave(file.path(output_dir, "Figure4a_erica_vs_moss.png"),
         p4a, width = 6, height = 5, dpi = 600, bg = "white")
  ggsave(file.path(output_dir, "Figure4b_erica_cover_vs_disturbance.png"),
         p4b, width = 6, height = 5, dpi = 600, bg = "white")
  ggsave(file.path(output_dir, "Figure4c_erica_seedling_vs_disturbance.png"),
         p4c, width = 6, height = 5, dpi = 600, bg = "white")
  ggsave(file.path(output_dir, "Figure4d_erica_germination_vs_disturbance.png"),
         p4d, width = 6, height = 5, dpi = 600, bg = "white")
  cat("  ✓ Figure 4 panels saved individually (install 'patchwork' for combined figure)\n")
}

## ── Figure 5: Elevational patterns ──────────────────────────────────────────
cat("\n=== FIGURE 5: Elevational patterns ===\n")

# Compute family richness per plot (need to reload species data for this)
veg_raw <- read_csv(file.path(data_dir, "Aboveground(stand)_3_column_dat.csv"),
                    col_types = cols(.default = "c"))
names(veg_raw) <- c("row", "plot", "species", "abundance")
veg_raw$abundance <- as.numeric(veg_raw$abundance)

# Load species list with family info
spp_list <- read_csv(file.path(data_dir, "TOTAL_SPECIES_LIST.csv"))
names(spp_list) <- c("row", "family", "species", "habit", "life_style", "endemicity")

# Merge to get family info
veg_fam <- veg_raw %>%
  left_join(spp_list %>% select(species, family), by = "species") %>%
  filter(!is.na(family))

family_richness <- veg_fam %>%
  group_by(plot) %>%
  summarise(family_richness = n_distinct(family)) %>%
  mutate(plot = as.numeric(plot))

df_fig5 <- df %>%
  left_join(family_richness, by = "plot")

# Species evenness from vegan
veg_wide <- veg_raw %>%
  mutate(plot = as.numeric(plot)) %>%
  select(plot, species, abundance) %>%
  pivot_wider(names_from = species, values_from = abundance, values_fill = 0) %>%
  arrange(plot)
veg_matrix <- as.data.frame(veg_wide[, -1])
rownames(veg_matrix) <- veg_wide$plot

df_fig5$evenness <- diversity(veg_matrix, "shannon") / log(specnumber(veg_matrix))

p5a <- ggplot(df_fig5, aes(x = elevation, y = family_richness)) +
  geom_point(color = col_moss, size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = col_moss, se = TRUE, alpha = 0.15) +
  labs(x = "Elevation (m a.s.l.)", y = "Family richness", tag = "(a)") +
  theme_bw(base_size = 12) +
  theme(plot.tag = element_text(face = "bold", size = 14))

p5b <- ggplot(df_fig5, aes(x = elevation, y = evenness)) +
  geom_point(color = col_moss, size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = col_moss, se = TRUE, alpha = 0.15) +
  labs(x = "Elevation (m a.s.l.)", y = "Species evenness (J')", tag = "(b)") +
  theme_bw(base_size = 12) +
  theme(plot.tag = element_text(face = "bold", size = 14))

p5c <- ggplot(df_fig5, aes(x = elevation, y = erica_cover)) +
  geom_point(color = col_n, size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), color = col_n, se = TRUE, alpha = 0.15) +
  labs(x = "Elevation (m a.s.l.)", y = "Standing Erica cover (%)", tag = "(c)") +
  theme_bw(base_size = 12) +
  theme(plot.tag = element_text(face = "bold", size = 14))

p5d <- ggplot(df_fig5, aes(x = elevation, y = germ_erica)) +
  geom_point(color = col_s, size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), color = col_s, se = TRUE, alpha = 0.15) +
  labs(x = "Elevation (m a.s.l.)", y = "Erica germination (seedlings m⁻²)", tag = "(d)") +
  theme_bw(base_size = 12) +
  theme(plot.tag = element_text(face = "bold", size = 14))

if (requireNamespace("patchwork", quietly = TRUE)) {
  p5 <- (p5a + p5b) / (p5c + p5d) +
    plot_annotation(title = "Elevational patterns in diversity and Erica regeneration",
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))
  ggsave(file.path(output_dir, "Figure5_elevational_patterns.png"),
         p5, width = 14, height = 11, dpi = 600, bg = "white")
  cat("  ✓ Figure 5 saved\n")
} else {
  ggsave(file.path(output_dir, "Figure5a_family_richness.png"), p5a, width = 6, height = 5, dpi = 600, bg = "white")
  ggsave(file.path(output_dir, "Figure5b_evenness.png"), p5b, width = 6, height = 5, dpi = 600, bg = "white")
  ggsave(file.path(output_dir, "Figure5c_erica_cover.png"), p5c, width = 6, height = 5, dpi = 600, bg = "white")
  ggsave(file.path(output_dir, "Figure5d_germination.png"), p5d, width = 6, height = 5, dpi = 600, bg = "white")
  cat("  ✓ Figure 5 panels saved individually\n")
}

## ── Figure 6: Seed density vs germination + disturbance effects ─────────────
cat("\n=== FIGURE 6: Seed density vs germination and disturbance ===\n")

p6a <- ggplot(df, aes(x = elevation, y = erica_cover)) +
  geom_point(color = col_n, size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = col_deco, linetype = "dashed", se = FALSE) +
  annotate("label", x = min(df$elevation) + 50, y = max(df$erica_cover, na.rm = TRUE) * 0.9,
           label = paste0("R² = ", round(summary(lm(erica_cover ~ elevation, df))$r.squared, 3),
                          "\np = ", formatC(summary(lm(erica_cover ~ elevation, df))$coef[2,4], format = "f", digits = 3), " n.s."),
           hjust = 0, size = 3.5, color = col_deco, fill = "white", label.size = 0.3) +
  labs(x = "Elevation (m a.s.l.)", y = "Erica seed density (seeds m⁻², proxy)",
       tag = "(a)") +
  theme_bw(base_size = 12) +
  theme(plot.tag = element_text(face = "bold", size = 14))

p6b <- ggplot(df, aes(x = disturbance, y = germ_erica)) +
  geom_point(color = col_dist, size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = col_dist, se = TRUE, alpha = 0.15) +
  annotate("label", x = min(df$disturbance) + 0.2, y = max(df$germ_erica, na.rm = TRUE) * 0.9,
           label = paste0("R² = ", round(summary(lm(germ_erica ~ disturbance, df))$r.squared, 3),
                          "\np = ", formatC(summary(lm(germ_erica ~ disturbance, df))$coef[2,4], format = "f", digits = 3)),
           hjust = 0, size = 3.5, color = col_dist, fill = "white", label.size = 0.3) +
  labs(x = "Disturbance intensity (D₀–D₄ scale)",
       y = "Erica germination (seedlings m⁻²)", tag = "(b)") +
  theme_bw(base_size = 12) +
  theme(plot.tag = element_text(face = "bold", size = 14))

if (requireNamespace("patchwork", quietly = TRUE)) {
  p6 <- p6a + p6b +
    plot_annotation(title = "Erica seed density, germination, and environmental gradients",
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))
  ggsave(file.path(output_dir, "Figure6_seed_density_germination.png"),
         p6, width = 12, height = 5.5, dpi = 600, bg = "white")
  cat("  ✓ Figure 6 saved\n")
} else {
  ggsave(file.path(output_dir, "Figure6a_seed_density_elevation.png"), p6a, width = 6, height = 5, dpi = 600, bg = "white")
  ggsave(file.path(output_dir, "Figure6b_germination_disturbance.png"), p6b, width = 6, height = 5, dpi = 600, bg = "white")
  cat("  ✓ Figure 6 panels saved individually\n")
}

cat("\n=== FIGURE GENERATION COMPLETE ===\n")
cat("Figures saved to", output_dir, "/ directory\n")
cat("Note: Figure 1 (study area map) and Figure 3 (CCA/NMDS ordination)\n")
cat("are generated separately using GIS software and R ordination plots.\n")
