# This script compares the four target urine metabolites by sex, age, BMI and fluid intake
# Builds on the long-format pipeline from the wide-pivot script (samples_cr),
# adding `sex`, `age`, `bmi_calculated` and period-matched fluid intake before pivoting/plotting.
# All four plots are combined into one patchwork figure.

## Setup
# Clear environment
rm(list = ls())

# Install required packages
library(readxl)
library(writexl)
library(tidyr)
library(rstatix)
library(purrr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(this.path)

# Set working directory and upload data
setwd(this.dir())
data <- read_excel("urine_targeted_gcms_data_processed.xlsx", sheet = "Sheet1")

dir.create("results", recursive = TRUE)

## Preparation
# Extract samples
samples <- data[data$screen_group %in% c(1, 2, 3) & !is.na(data$screen_group), ]
dim(samples)

# Check sex is present and inspect its levels/labels
table(samples$sex, useNA = "ifany")

# Add timepoint column (w/o at-home samples)
samples <- samples %>%
  mutate(timepoint = case_when(
    grepl("BL", Name) ~ "bl2",
    grepl("0-3", Name) ~ "p1",
    grepl("3-6", Name) ~ "p2",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(timepoint))

# Group and correct compounds of interest
vars <- c("galactose_ug", "lactose_ug", "galactitol_ug", "galactonate_ug")
samples[vars] <- lapply(samples[vars], function(x) as.numeric(as.character(x)))

# Add fluid intake matching each urine collection period
# (fluid columns are per participant, so each sample gets the intake of its own collection period)
fluid_vars <- c("water_consumed_b2_ml", "water_consumed_p1_ml", "water_consumed_p2_ml")
samples[fluid_vars] <- lapply(samples[fluid_vars], function(x) as.numeric(as.character(x)))

samples <- samples %>%
  mutate(fluid_ml = case_when(
    timepoint == "bl2" ~ water_consumed_b2_ml,
    timepoint == "p1"  ~ water_consumed_p1_ml,
    timepoint == "p2"  ~ water_consumed_p2_ml))

# Pivot compounds to rows, keeping sex/age/bmi/fluid alongside the usual grouping columns
samples_cr <- samples %>%
  dplyr::select(study_id, screen_group, sex, age, bmi_calculated, fluid_ml, timepoint, all_of(vars)) %>%
  pivot_longer(cols = all_of(vars), names_to = "compound", values_to = "value")
dim(samples_cr)

# Remove incomplete measurements (participants missing any timepoint)
all_timepoints <- unique(samples_cr$timepoint)
samples_cr <- samples_cr %>%
  group_by(compound, study_id) %>%
  filter(all(all_timepoints %in% timepoint)) %>%
  ungroup()
dim(samples_cr)

# Check group sizes by sex (overall, and per timepoint/compound if needed)
samples_cr %>%
  distinct(study_id, sex) %>%
  count(sex)


## ============================================================
## SEX — Statistical comparison (Wilcoxon test per compound x timepoint)
## ============================================================
# (non-parametric, since metabolite concentrations are typically right-skewed;
# switch to t_test() below if you'd rather assume normality)
sex_stats <- samples_cr %>%
  group_by(compound, timepoint) %>%
  wilcox_test(value ~ sex) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

print(sex_stats, n = Inf)

write_xlsx(sex_stats, "results/metabolites_by_sex_wilcoxon.xlsx")


## ============================================================
## AGE, BMI & FLUID INTAKE — Statistical comparison (Spearman correlation per compound x timepoint)
## ============================================================
age_stats <- samples_cr %>%
  group_by(compound, timepoint) %>%
  summarise(
    rho     = cor.test(value, age, method = "spearman", exact = FALSE)$estimate,
    p_value = cor.test(value, age, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  add_significance("p_adj")

print(age_stats, n = Inf)
write_xlsx(age_stats, "results/metabolites_by_age_spearman.xlsx")

bmi_stats <- samples_cr %>%
  group_by(compound, timepoint) %>%
  summarise(
    rho     = cor.test(value, bmi_calculated, method = "spearman", exact = FALSE)$estimate,
    p_value = cor.test(value, bmi_calculated, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  add_significance("p_adj")

print(bmi_stats, n = Inf)
write_xlsx(bmi_stats, "results/metabolites_by_bmi_spearman.xlsx")

# Fluid intake: baseline excluded - intake was fixed at 650 mL for everyone (no variance to correlate)
fluid_stats <- samples_cr %>%
  filter(timepoint != "bl2") %>%
  group_by(compound, timepoint) %>%
  summarise(
    rho     = cor.test(value, fluid_ml, method = "spearman", exact = FALSE)$estimate,
    p_value = cor.test(value, fluid_ml, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  add_significance("p_adj")

print(fluid_stats, n = Inf)
write_xlsx(fluid_stats, "results/metabolites_by_fluid_spearman.xlsx")


## ============================================================
## PLOT SETUP (shared)
## ============================================================

compound_order <- c("lactose_ug", "galactose_ug", "galactonate_ug", "galactitol_ug")
compound_labels <- c("lactose_ug"     = "Lactose",
                     "galactose_ug"   = "Galactose",
                     "galactonate_ug" = "Galactonate",
                     "galactitol_ug"  = "Galactitol")

tp_levels <- c("bl2", "p1", "p2")
tp_labels <- c("BL", "0-3h", "3-6h")

sex_colors <- c("Female" = "#377EB8", "Male" = "#A6761D")  # colorblind-friendly, distinct from screen_group palette

samples_plot <- samples_cr %>%
  mutate(
    sex       = factor(sex, levels = c(1, 2), labels = c("Female", "Male")),
    timepoint = factor(timepoint, levels = tp_levels, labels = tp_labels),
    compound  = factor(compound, levels = compound_order, labels = compound_labels[compound_order]))

# Same factor labels for the stats tables (so rho labels land in the right facet)
to_plot_factors <- function(df) {
  df %>%
    mutate(
      timepoint = factor(timepoint, levels = tp_levels, labels = tp_labels),
      compound  = factor(compound, levels = compound_order, labels = compound_labels[compound_order]))
}

# Shared theme for the individual plots (matches the finalized box/line figures)
theme_pub <- theme_light(base_size = 15) +
  theme(panel.grid        = element_blank(),
        strip.background  = element_blank(),
        strip.text        = element_text(face = "bold", size = 20, color = "black"),
        axis.title        = element_text(size = 22.5),
        axis.text         = element_text(size = 21),
        plot.title        = element_text(face = "bold", size = 25, hjust = 0.5))


## ============================================================
## SEX PLOT — box plot grid (compound rows x timepoint columns)
## ============================================================

stat_test_sex <- samples_plot %>%
  group_by(compound, timepoint) %>%
  wilcox_test(value ~ sex) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "sex")

p_sex <- ggplot(samples_plot, aes(x = sex, y = value, fill = sex,
                                  group = interaction(timepoint, sex))) +
  geom_boxplot(outlier.size = 0.8, linewidth = 0.6) +
  stat_pvalue_manual(stat_test_sex, label = "p.adj.signif", tip.length = 0.01, size = 6.8) +  # geom text size is in mm: +5 pt = +1.8 mm
  facet_grid(compound ~ timepoint, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +  # headroom so the stars above the brackets aren't cut off
  labs(x = "Sex", y = "Amount / µg", fill = NULL) +
  scale_fill_manual(values = sex_colors) +
  theme_pub


## ============================================================
## AGE, BMI & FLUID INTAKE PLOTS — scatter plot grid (compound rows x timepoint columns)
## ============================================================

make_covariate_scatter <- function(data, covariate, covariate_label, stats_table, panel_title) {
  
  stats_labels <- stats_table %>%
    mutate(label = paste0("rho = ", round(rho, 2), ", p.adj ", p_adj.signif))
  
  ggplot(data, aes(x = .data[[covariate]], y = value)) +
    geom_point(alpha = 0.6, size = 1.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
    geom_text(data = stats_labels,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.3, size = 5.8, inherit.aes = FALSE) +  # +5 pt = +1.8 mm
    facet_grid(compound ~ timepoint, scales = "free") +
    labs(x = covariate_label, y = "Amount / µg") +
    theme_pub
}

p_age <- make_covariate_scatter(samples_plot, "age", "Age (years)",
                                to_plot_factors(age_stats), "Age")

p_bmi <- make_covariate_scatter(samples_plot, "bmi_calculated", "BMI (kg/m²)",
                                to_plot_factors(bmi_stats), "BMI")

p_fluid <- make_covariate_scatter(samples_plot %>% filter(timepoint != "BL"),
                                  "fluid_ml", "Fluid intake (mL)",
                                  to_plot_factors(fluid_stats), "Fluid intake")


## ============================================================
## COMBINED FIGURE
## ============================================================

p_final_cov <- (p_sex | p_age) / (p_bmi | p_fluid) +
  plot_layout(guides = "collect") &
  theme(axis.text.y     = element_text(size = 19),
        axis.title.x    = element_text(margin = margin(t = 9)),
        plot.margin     = margin(5, 10, 5, 10),
        legend.position = "none")

p_final_cov_with_title <- p_final_cov + plot_annotation(
  #title = "Lactose, galactose, galactonate and galactitol profiles by sex, age, BMI and fluid intake",
  theme = theme(plot.title = element_text(size = 35)))

print(p_final_cov_with_title)

ggsave("results/pub_plots/metabolites_by_sex_age_bmi_fluid.png", p_final_cov_with_title,
       width = 24, height = 20, dpi = 600, device = "png", create.dir = TRUE)
