## ============================================================
## ROC curves: targeted urine metabolites AND breath hydrogen
## Combined: 0-3h, 3-6h, and pooled 0-6h - sex-adjusted
## ============================================================

## Setup
rm(list = ls())

library(readxl)
library(readr)
library(writexl)
library(purrr)
library(tidyverse)
library(pROC)
library(ggplot2)
library(caret)
library(this.path)

# ---- EDIT THESE PATHS IF NEEDED ----
dir_gcms <- this.dir()

setwd(dir_gcms)
data <- read_excel("urine_targeted_gcms_data_processed.xlsx", sheet = "Sheet1")

bh2_data <- read_csv("Quintron_H2_CH4_20260407.csv")


## --- Process breath hydrogen data ---
bh2_data <- bh2_data %>%
  pivot_wider(names_from = Timepoint, values_from = c(H2, CH4, CO2, Corr))

bh2_data <- bh2_data %>%
  mutate(
    H2_p1_mean = rowMeans(dplyr::select(., paste0("H2_", c("15", "30", "60", "90", "120"))), na.rm = TRUE),
    H2_p2_mean = rowMeans(dplyr::select(., paste0("H2_", c("180", "240", "300", "360"))), na.rm = TRUE)
  ) %>%
  dplyr::select(Participant_ID, H2_p1_mean, H2_p2_mean)


## --- Merge breath hydrogen onto GC-MS data ---
data <- data %>%
  left_join(bh2_data, by = c("study_id" = "Participant_ID"))

cat("Rows with no matching breath hydrogen data:",
    sum(is.na(data$H2_p1_mean) & is.na(data$H2_p2_mean)), "out of", nrow(data), "\n")


## Prep
comparisons <- list(
  list(name = "Controls vs LM", groups = c(1, 2)),
  list(name = "Controls vs LI", groups = c(1, 3)),
  list(name = "LM vs LI",       groups = c(2, 3))
)

comparison_levels <- map_chr(comparisons, "name")

timepoints <- c("0-3h", "3-6h")

## Combined model list: metabolite models + breath hydrogen model
model_versions <- list(
  list(label = "All four target compounds", vars = c("lactose_ug", "galactose_ug", "galactitol_ug", "galactonate_ug")),
  list(label = "Without lactose",           vars = c("galactose_ug", "galactitol_ug", "galactonate_ug")),
  list(label = "Galactonate & galactitol",  vars = c("galactitol_ug", "galactonate_ug")),
  list(label = "Galactitol",                vars = c("galactitol_ug")),
  list(label = "Galactonate",               vars = c("galactonate_ug")),
  list(label = "Breath hydrogen",           vars = c("H2_val"))
)

## --- Pooling function for 0-6h (averages 0-3h + 3-6h per participant)
prep_pooled_data <- function(data, comp, vars) {
  
  data %>%
    dplyr::filter(treatment == "Treatment A") %>%
    dplyr::filter(Type == "Sample") %>%
    dplyr::filter(!grepl("_AH", Name)) %>%
    mutate(
      timepoint = case_when(
        grepl("0-3", Name) ~ "0-3h",
        grepl("3-6", Name) ~ "3-6h",
        TRUE ~ NA_character_),
      screen_group = factor(screen_group),
      sex = factor(sex),
      H2_val = case_when(
        timepoint == "0-3h" ~ H2_p1_mean,
        timepoint == "3-6h" ~ H2_p2_mean,
        TRUE ~ NA_real_)
    ) %>%
    dplyr::filter(timepoint %in% c("0-3h", "3-6h")) %>%
    dplyr::filter(screen_group %in% comp$groups) %>%
    mutate(screen_group = droplevels(screen_group)) %>%
    group_by(study_id, screen_group, sex) %>%
    summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), n_tp = n(), .groups = "drop") %>%
    dplyr::filter(n_tp == 2) %>%
    dplyr::select(-n_tp) %>%
    dplyr::filter(complete.cases(pick(all_of(c(vars, "sex")))))
}

## CV parameters
n_repeats <- 25
n_folds   <- 5

roc_results <- list()

## --- Shared CV/ROC runner, so both loops below use identical logic ---
run_cv_roc <- function(df, vars, comp) {
  
  fold_aucs         <- c()
  all_sensitivities <- list()
  fold_idx          <- 1
  
  set.seed(123)
  multi_folds <- createMultiFolds(df$screen_group, k = n_folds, times = n_repeats)
  
  for (fold_name in names(multi_folds)) {
    
    train_idx <- multi_folds[[fold_name]]
    test_idx  <- setdiff(seq_len(nrow(df)), train_idx)
    train_df  <- df[train_idx, ]
    test_df   <- df[test_idx, ]
    
    train_counts <- table(train_df$screen_group)
    train_w      <- 1 / train_counts[as.character(train_df$screen_group)]
    train_w      <- as.numeric(train_w / sum(train_w) * nrow(train_df))
    
    formula_str <- paste("screen_group ~", paste(c(vars, "sex"), collapse = " + "))
    
    glm_fit <- tryCatch(
      glm(as.formula(formula_str), data = train_df,
          family = binomial(), weights = train_w),
      error = function(e) NULL
    )
    
    if (is.null(glm_fit)) next
    
    test_scores <- predict(glm_fit, newdata = test_df, type = "response")
    
    roc_obj <- tryCatch(
      roc(test_df$screen_group, test_scores,
          levels = comp$groups,
          direction = "<",
          quiet = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(roc_obj)) next
    
    fold_aucs <- c(fold_aucs, as.numeric(auc(roc_obj)))
    
    spec_points <- seq(0, 1, by = 0.01)
    roc_coords  <- coords(roc_obj, "all", ret = c("specificity", "sensitivity")) %>%
      as.data.frame() %>%
      arrange(specificity)
    
    sens_interp <- approx(x      = roc_coords$specificity,
                          y      = roc_coords$sensitivity,
                          xout   = spec_points,
                          method = "constant",
                          rule   = 2)$y
    
    all_sensitivities[[fold_idx]] <- data.frame(
      specificity = spec_points,
      sensitivity = sens_interp
    )
    
    fold_idx <- fold_idx + 1
  }
  
  if (length(all_sensitivities) == 0) return(NULL)
  
  mean_roc <- bind_rows(all_sensitivities) %>%
    group_by(specificity) %>%
    summarise(
      mean_sensitivity = mean(sensitivity, na.rm = TRUE),
      sd_sensitivity   = sd(sensitivity, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(mean_auc = mean(fold_aucs), sd_auc = sd(fold_aucs), mean_roc = mean_roc, n = nrow(df))
}

## --- Loop 1: 0-3h and 3-6h (per-timepoint, unpooled) ---
for (mv in model_versions) {
  for (comp in comparisons) {
    for (tp in timepoints) {
      
      df <- data %>%
        dplyr::filter(treatment == "Treatment A") %>%
        dplyr::filter(Type == "Sample") %>%
        dplyr::filter(!grepl("_AH", Name)) %>%
        mutate(
          timepoint = case_when(
            grepl("0-3", Name) ~ "0-3h",
            grepl("3-6", Name) ~ "3-6h",
            TRUE ~ NA_character_),
          screen_group = factor(screen_group),
          sex = factor(sex),
          H2_val = case_when(
            timepoint == "0-3h" ~ H2_p1_mean,
            timepoint == "3-6h" ~ H2_p2_mean,
            TRUE ~ NA_real_)
        ) %>%
        dplyr::filter(timepoint == tp) %>%
        dplyr::filter(screen_group %in% comp$groups) %>%
        mutate(screen_group = droplevels(screen_group)) %>%
        dplyr::filter(complete.cases(pick(all_of(c(mv$vars, "sex")))))
      
      if (nrow(df) < 10) next
      
      res <- run_cv_roc(df, mv$vars, comp)
      if (is.null(res)) next
      
      roc_results[[length(roc_results) + 1]] <- list(
        comparison = comp$name,
        timepoint  = tp,
        model      = mv$label,
        mean_auc   = res$mean_auc,
        sd_auc     = res$sd_auc,
        mean_roc   = res$mean_roc
      )
    }
  }
}

## --- Loop 2: pooled 0-6h ---
for (mv in model_versions) {
  for (comp in comparisons) {
    
    df <- prep_pooled_data(data, comp, mv$vars)
    
    if (nrow(df) < 10) next
    
    res <- run_cv_roc(df, mv$vars, comp)
    if (is.null(res)) next
    
    roc_results[[length(roc_results) + 1]] <- list(
      comparison = comp$name,
      timepoint  = "0-6h",
      model      = mv$label,
      mean_auc   = res$mean_auc,
      sd_auc     = res$sd_auc,
      mean_roc   = res$mean_roc
    )
  }
}

## Summary table
roc_summary <- map_dfr(roc_results, function(r) {
  data.frame(
    comparison = r$comparison,
    timepoint  = r$timepoint,
    model      = r$model,
    AUC        = round(r$mean_auc, 3),
    AUC_SD     = round(r$sd_auc, 3)
  )
}) %>%
  mutate(
    timepoint  = factor(timepoint, levels = c("0-3h", "3-6h", "0-6h")),
    comparison = factor(comparison, levels = comparison_levels)
  ) %>%
  arrange(comparison, timepoint, desc(AUC))

print(roc_summary)


## ============================================================
## Plot — 3 rows (0-3h, 3-6h, 0-6h) x comparisons
## ============================================================

models_to_plot <- c("All four target compounds", "Without lactose",
                    "Galactonate & galactitol", "Breath hydrogen")

model_colors <- c(
  "All four target compounds" = "#D95F02",
  "Without lactose"           = "#E7298A",
  "Galactonate & galactitol"  = "#666666",
  "Breath hydrogen"           = "#1B9E77"
)

model_linetypes <- c(
  "All four target compounds" = "solid",
  "Without lactose"           = "dashed",
  "Galactonate & galactitol"  = "dotdash",
  "Breath hydrogen"           = "longdash"
)

label_y <- c(
  "All four target compounds" = 0.25,
  "Without lactose"           = 0.18,
  "Galactonate & galactitol"  = 0.11,
  "Breath hydrogen"           = 0.04
)

roc_plot_data <- map_dfr(roc_results, function(r) {
  if (!r$model %in% models_to_plot) return(NULL)
  r$mean_roc %>%
    mutate(
      comparison = factor(r$comparison, levels = comparison_levels),
      timepoint  = factor(r$timepoint, levels = c("0-3h", "3-6h", "0-6h")),
      model      = r$model
    )
})

auc_labels <- roc_summary %>%
  dplyr::filter(model %in% models_to_plot) %>%
  mutate(
    label = paste0(model, ": ", round(AUC, 2), " \u00b1 ", round(AUC_SD, 2)),
    y_pos = label_y[model]
  )

p <- ggplot(roc_plot_data, aes(x = 1 - specificity, y = mean_sensitivity,
                               color = model, linetype = model)) +
  geom_segment(x = 0, xend = 1, y = 0, yend = 1,
               linetype = "dashed", color = "black", linewidth = 0.4,
               inherit.aes = FALSE) +
  geom_line(linewidth = 0.9) +
  geom_text(
    data = auc_labels,
    aes(label = label, color = model, y = y_pos),
    x = 0.29, hjust = 0, size = 4.5, show.legend = FALSE
  ) +
  scale_color_manual(values = model_colors, breaks = names(model_colors)) +
  scale_linetype_manual(values = model_linetypes, breaks = names(model_linetypes)) +
  facet_grid(timepoint ~ comparison) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(
    x        = "1 - Specificity (FPR)",
    y        = "Sensitivity (TPR)",
    color    = "Model",
    linetype = "Model"
  ) +
  theme_bw(base_size = 18) +  
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 18),  
    axis.title       = element_text(size = 18),                
    axis.text        = element_text(size = 15),                
    legend.title     = element_text(face = "bold", size = 18),  
    legend.position  = "bottom",
    panel.grid       = element_blank(),
    legend.text      = element_text(size = 19)                  
  )
print(p)

## Save
ggsave("results/pub_plots/roc_curves_metabolites_h2_all_timepoints_with_sex.png", plot = p,
       width = 15, height = 14, dpi = 300, create.dir = TRUE)
write_xlsx(roc_summary, paste0("results/roc_summary_metabolites_h2_all_timepoints_", n_repeats, ".xlsx"))

cat("Min AUC:", min(roc_summary$AUC), "\n")
cat("Max AUC:", max(roc_summary$AUC), "\n")
cat("Min SD: ", min(roc_summary$AUC_SD), "\n")
cat("Max SD: ", max(roc_summary$AUC_SD), "\n")

