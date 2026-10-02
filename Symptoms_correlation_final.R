# This script correlates targeted urine GC-MS data with symptoms

### Setup
## Clear environment
rm(list = ls())


## Install required packages
library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(writexl)
library(tidyverse)
library(this.path)


## Set working directories and upload datasets
wdir <- this.dir()

setwd(wdir)
gcms_data <- read_excel("urine_targeted_gcms_data_processed_wide.xlsx", sheet = "Sheet1")
dim(gcms_data) #111 18
names(gcms_data) # including some metadata, participants encoded as study_id
head(gcms_data) # study_id without underscore

symp_data <- read_excel("gi_symptoms_data_processed.xlsx", sheet = "Sheet1")
dim(symp_data) #112 63
names(symp_data) # no metadata, participants encoded as study_id
head(symp_data) # study_id without underscore


## Merge datasets
# Check for mismatches
ids_gcms <- unique(gcms_data$study_id)
ids_symp <- unique(symp_data$study_id)

only_in_gcms <- setdiff(ids_gcms, ids_symp)
only_in_symp <- setdiff(ids_symp, ids_gcms)

if (length(only_in_gcms) > 0) warning(paste("IDs in GCMs data but not in symptoms data:", paste(only_in_gcms, collapse = ", ")))
if (length(only_in_symp) > 0) warning(paste("IDs in symptoms data but not in GCMs data:", paste(only_in_symp, collapse = ", ")))
# LB122 was excluded because unable to provide second baseline urine sample despite multiple attempts

# Merge
data <- inner_join(gcms_data, symp_data, by = "study_id")
dim(data) #111 81

# Filter by study and intervention groups
data <- data %>% filter(treatment == "Treatment A")
#data <- data %>% filter(screen_group == 3)
dim(data) #91 81



### Subset by timepoint
## Check timepoints
names(symp_data)

# bl_mean <-- bl2 urine pool (bl1 excluded from targeted analysis)
# +30, 60, 120 min <-- urine pool 1
# +360 min <-- urine pool 2

# not straightforward symptom timepoints assignments: +180 min, +540 min
# + 180 min <-- assign to urine pool 1, as symptoms lag behind metabolism
# + 540 min <-- exclude from analysis, as 3 hours apart from urine sampling


## Average symptoms within urine pool windows
symptoms <- c("abd_pain", "nausea", "bloating", "flatulence", "diarrhoea")

for (sym in symptoms) {
  # p1 average (30-180min)
  p1_cols <- paste0(sym, c("_30min", "_60min", "_120min", "_180min"))
  data[[paste0(sym, "_p1_mean")]] <- rowMeans(data[, p1_cols], na.rm = TRUE)
  
  # p2 average (180-360min)
  p2_cols <- paste0(sym, c("_180min", "_360min"))
  data[[paste0(sym, "_p2_mean")]] <- rowMeans(data[, p2_cols], na.rm = TRUE)
}

# Check
names(data)[grepl("p1_mean|p2_mean", names(data))]



### Spearman correlation — stratified by group
## Prep
metabolites <- c("galactose_ug", "lactose_ug", "galactitol_ug", "galactonate_ug")
symptoms    <- c("abd_pain", "nausea", "bloating", "flatulence", "diarrhoea")
timepoints  <- c("p1", "p2")

cor_results_group <- list()


## Run
for (grp in c(1, 2, 3)) {
  for (tp in timepoints) {
    for (met in metabolites) {
      met_col <- paste0(met, "_", tp)
      for (sym in symptoms) {
        sym_col <- paste0(sym, "_", tp, "_mean")
        
        if (met_col %in% names(data) & sym_col %in% names(data)) {
          df <- data %>% filter(screen_group == grp)
          
          test <- tryCatch(
            cor.test(df[[met_col]], df[[sym_col]], method = "spearman", exact = FALSE),
            error = function(e) NULL)
          
          if (is.null(test)) next
          
          cor_results_group[[length(cor_results_group) + 1]] <- data.frame(
            screen_group = grp,
            timepoint    = tp,
            metabolite   = met,
            parameter    = sym,
            rho          = test$estimate,
            p_value      = test$p.value)
        }
      }
    }
  }
}


## Merge
results_spearman_group <- do.call(rbind, cor_results_group)


## FDR correction — within each group
results_spearman_group <- results_spearman_group %>%
  group_by(screen_group) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(significant = p_adj < 0.05)


## Clean up labels
results_spearman_group <- results_spearman_group %>%
  mutate(
    screen_group = factor(screen_group,
                          levels = c(1, 2, 3),
                          labels = c("Controls", "LM", "LI")),
    timepoint = recode(timepoint,
                       "p1" = "0-3h",
                       "p2" = "3-6h"),
    timepoint  = factor(timepoint, levels = c("0-3h", "3-6h")),
    metabolite = recode(metabolite,
                        "galactose_ug"   = "Galactose / µg",
                        "lactose_ug"     = "Lactose / µg",
                        "galactitol_ug"  = "Galactitol / µg",
                        "galactonate_ug" = "Galactonate / µg"),
    metabolite = factor(metabolite,
                        levels = c("Lactose / µg", "Galactose / µg",
                                   "Galactonate / µg", "Galactitol / µg")),
    parameter = recode(parameter,
                       "abd_pain"   = "Abdominal pain",
                       "nausea"     = "Nausea",
                       "bloating"   = "Bloating",
                       "flatulence" = "Flatulence",
                       "diarrhoea"  = "Diarrhoea"),
    parameter = factor(parameter,
                       levels = c("Nausea", "Flatulence", "Diarrhoea", "Bloating", "Abdominal pain")))



### Visualization
## Main plot — all screen groups, faceted
plot_group <- ggplot(results_spearman_group, aes(x = metabolite, y = parameter, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(significant, "*", "")),
            size = 5, vjust = 0.5, hjust = 0.5) +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  facet_grid(screen_group ~ timepoint) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text       = element_text(size = 11, face = "bold"),
        plot.title       = element_text(size = 14),
        axis.title.y     = element_text(margin = margin(r = 15))) +
  labs(title   = "Spearman correlations between targeted urine metabolites and GI symptoms",
       x       = "Urine metabolite",
       y       = "Symptom",
       caption = "* Adjusted p < 0.05")

plot_group


## LI-only plot
plot_li <- results_spearman_group %>%
  filter(screen_group == "LI") %>%
  ggplot(aes(x = metabolite, y = parameter, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(significant, "*", "")),
            size = 5, vjust = 0.5, hjust = 0.5) +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  facet_wrap(~ timepoint, ncol = 2) +
  theme_bw() +
  theme(legend.title     = element_text(size = 13),
        #legend.text      = element_text(size = 12),
        axis.text        = element_text(size = 12),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        axis.title       = element_text(size = 13),
        axis.title.y     = element_text(margin = margin(r = 15)),
        axis.title.x     = element_text(margin = margin(t = 19)),
        strip.background = element_blank(),
        plot.background  = element_blank(),
        panel.background = element_blank(),
        strip.text       = element_text(size = 11, face = "bold"),
        plot.title       = element_text(size = 14),
        plot.caption = element_text(size = 12)) +
  labs(x       = "Urine metabolite",
       y       = "Symptom",
       caption = "LI group only (n = 51)")

plot_li


## Export
ggsave("results/pub_plots/heatmap_sc_stratified.pdf", plot = plot_group, width = 10, height = 9, create.dir = TRUE)
ggsave("results/pub_plots/heatmap_sc_LI_only.pdf",    plot = plot_li,    width = 10, height = 5, create.dir = TRUE)


## Read out rho values
# All rho values
results_spearman_group %>%
  dplyr::select(screen_group, timepoint, metabolite, parameter, rho, p_adj, significant) %>%
  arrange(screen_group, timepoint, p_adj) %>%
  print()

# Only significant
results_spearman_group %>%
  filter(significant == TRUE) %>%
  dplyr::select(screen_group, timepoint, metabolite, parameter, rho, p_adj) %>%
  arrange(screen_group, timepoint, p_adj) %>%
  print()

# LI group only
results_spearman_group %>%
  filter(screen_group == "LI") %>%
  dplyr::select(timepoint, metabolite, parameter, rho, p_adj, significant) %>%
  arrange(timepoint, p_adj) %>%
  print()
