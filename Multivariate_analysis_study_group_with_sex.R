# This script performs multivariate analyses on processed targeted urine GC-MS 

## Setup
# Clear environment
rm(list = ls())

# Install required packages
library(readxl)
library(writexl)
library(purrr)
library(tidyverse)
library(MASS) # <-- for LDA
library(this.path)


# Set working directory and upload data
setwd(this.dir())
data <- read_excel("urine_targeted_gcms_data_processed.xlsx", sheet = "Sheet1")
dim(data)
names(data)
length(unique(data$study_id[!is.na(data$study_id)])) #111

# Set timepoint - change this to switch between "0-3h", "3-6h", "0-6h"
TIMEPOINT <- "0-6h"

# Create function for BER calculation
calc_ber <- function(predicted, actual) {
  classes <- levels(factor(actual))
  class_errors <- sapply(classes, function(cls) {
    idx <- actual == cls
    1 - mean(predicted[idx] == actual[idx])
  })
  mean(class_errors)
}

# Create function for LDA data prep (handles "0-3h", "3-6h", and pooled "0-6h" the same way)
prep_lda_data <- function(data, timepoint_choice, groups, vars) {
  
  df <- data %>%
    filter(treatment == "Treatment A") %>%
    filter(Type == "Sample") %>%
    filter(!grepl("_AH", Name)) %>%
    mutate(
      timepoint = case_when(
        grepl("BL", Name) ~ "Baseline",
        grepl("0-3", Name) ~ "0-3h",
        grepl("3-6", Name) ~ "3-6h",
        TRUE ~ NA_character_),
      timepoint = factor(timepoint, levels = c("Baseline", "0-3h", "3-6h")),
      screen_group = factor(screen_group),
      sex = factor(sex)
    )
  
  if (timepoint_choice %in% c("0-3h", "3-6h")) {
    df <- df %>%
      filter(timepoint == timepoint_choice)
  } else if (timepoint_choice == "0-6h") {
    df <- df %>%
      filter(timepoint %in% c("0-3h", "3-6h")) %>%
      group_by(study_id, screen_group, sex) %>%
      summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), n_tp = n(), .groups = "drop") %>%
      filter(n_tp == 2) %>% # only pool participants with both windows present
      dplyr::select(-n_tp)
  } else {
    stop('timepoint_choice must be one of "0-3h", "3-6h", "0-6h"')
  }
  
  df %>%
    filter(screen_group %in% groups) %>%
    mutate(screen_group = droplevels(screen_group)) %>%
    dplyr::select(screen_group, sex, all_of(vars))
}


## Linear discriminant analyses (LDA)
# Across all groups
lda_data_w = prep_lda_data(data, TIMEPOINT, groups = c(1, 2, 3),
                           vars = c("galactitol_ug", "galactonate_ug"))

lda_data_w %>%
  pivot_longer(!c(screen_group, sex), names_to = "metabolite", values_to = "amount") %>%
  ggplot(aes(x = screen_group, y = amount, fill = screen_group)) + 
  geom_boxplot() +
  facet_grid(metabolite ~ .) +
  theme_light()

w = lda(screen_group~., lda_data_w)
w_cv = lda(screen_group~., lda_data_w, CV = TRUE)
w_acc <- mean(w_cv$class == lda_data_w$screen_group) * 100
w_acc
ber_w <- calc_ber(w_cv$class, lda_data_w$screen_group)
ber_w
plot(w)
w

table(predicted = predict(w)$class, actual = lda_data_w$screen_group) # Group 2 is never correctly classified
sum(diag(table(predicted = predict(w)$class, actual = lda_data_w$screen_group))) / nrow(lda_data_w) # Accuracy of 0.73 at 0-3h (0.73 at 3-6h)

# As expected, the model cannot distinguish between groups 2 and 3

# Group 1 Vs. 2
lda_data_x = prep_lda_data(data, TIMEPOINT, groups = c(1, 2),
                           vars = c("lactose_ug", "galactose_ug", "galactitol_ug", "galactonate_ug"))

lda_data_x %>%
  pivot_longer(!c(screen_group, sex), names_to = "metabolite", values_to = "amount") %>%
  ggplot(aes(x = screen_group, y = amount, fill = screen_group)) + 
  geom_boxplot() +
  facet_grid(metabolite ~ .) +
  theme_light()

x = lda(screen_group~., lda_data_x)
x_cv = lda(screen_group~., lda_data_x, CV = TRUE)
x_acc <- mean(x_cv$class == lda_data_x$screen_group) * 100
x_acc
ber_x <- calc_ber(x_cv$class, lda_data_x$screen_group)
ber_x
plot(x)
x

table(predicted = predict(x)$class, actual = lda_data_x$screen_group) # Only missclassifies one person
sum(diag(table(predicted = predict(x)$class, actual = lda_data_x$screen_group))) / nrow(lda_data_x) # Accuracy of 0.93 at 0-3h (0.90 at 3-6h)

# Group 1 Vs. 3
lda_data_y = prep_lda_data(data, TIMEPOINT, groups = c(1, 3),
                           vars = c("lactose_ug", "galactose_ug", "galactitol_ug", "galactonate_ug"))

lda_data_y  %>%
  pivot_longer(!c(screen_group, sex), names_to = "metabolite", values_to = "amount") %>%
  ggplot(aes(x = screen_group, y = amount, fill = screen_group)) + 
  geom_boxplot() +
  facet_grid(metabolite ~ .) +
  theme_light()

y = lda(screen_group~., lda_data_y)
y_cv = lda(screen_group~., lda_data_y, CV = TRUE)
y_acc <- mean(y_cv$class == lda_data_y$screen_group) * 100
y_acc
ber_y <- calc_ber(y_cv$class, lda_data_y$screen_group)
ber_y
plot(y)
y

table(predicted = predict(y)$class, actual = lda_data_y$screen_group) # Missclassifies two people
sum(diag(table(predicted = predict(y)$class, actual = lda_data_y$screen_group))) / nrow(lda_data_y) # Accuracy of 0.92 at 0-3h (0.94 at 3-6)

# Group 2 Vs. 3
lda_data_z = prep_lda_data(data, TIMEPOINT, groups = c(2, 3),
                           vars = c("lactose_ug", "galactose_ug", "galactitol_ug", "galactonate_ug"))

lda_data_z %>%
  pivot_longer(!c(screen_group, sex), names_to = "metabolite", values_to = "amount") %>%
  ggplot(aes(x = screen_group, y = amount, fill = screen_group)) + 
  geom_boxplot() +
  facet_grid(metabolite ~ .) +
  theme_light()

z = lda(screen_group~., lda_data_z)
z_cv = lda(screen_group~., lda_data_z, CV = TRUE)
z_acc <- mean(z_cv$class == lda_data_z$screen_group) * 100
z_acc
ber_z <- calc_ber(z_cv$class, lda_data_z$screen_group)
ber_z
plot(z)
z

table(predicted = predict(z)$class, actual = lda_data_z$screen_group) # Missclassifies two people
sum(diag(table(predicted = predict(z)$class, actual = lda_data_z$screen_group))) / nrow(lda_data_z) # Accuracy of 0.74 at 0-3h (0.70 at 3-6h)

# Surprisingly we do see accuracy above 50% when all metabolites are considered


## Summary table 
lda_summary <- data.frame(
  timepoint        = TIMEPOINT,
  comparison       = c("All groups", "Controls vs LM", "Controls vs LI", "LM vs LI"),
  LD1_variance_pct = c(round(w$svd^2 / sum(w$svd^2) * 100, 1)[1], NA, NA, NA),
  CV_accuracy      = c(w_acc, x_acc, y_acc, z_acc),
  BER              = round(c(ber_w, ber_x, ber_y, ber_z), 3)
)

print(lda_summary)


##Visualize
models <- list(
  list(model = w, model_cv = w_cv, data = lda_data_w, name = "All groups"),
  list(model = x, model_cv = x_cv, data = lda_data_x, name = "Controls vs LM"),
  list(model = y, model_cv = y_cv, data = lda_data_y, name = "Controls vs LI"),
  list(model = z, model_cv = z_cv, data = lda_data_z, name = "LM vs LI"))

# Build plot data
plot_data <- purrr::map_dfr(models, function(m) {
  pred   <- predict(m$model)
  actual <- m$data$screen_group
  cv_acc <- mean(m$model_cv$class == actual) * 100  # <- NEW
  classes     <- levels(factor(actual))             # <- NEW
  class_errors <- sapply(classes, function(cls) {
    idx <- actual == cls
    1 - mean(m$model_cv$class[idx] == actual[idx])
  })
  ber <- mean(class_errors)
  
  data.frame(
    LD1 = pred$x[, 1],
    group = factor(actual,
                   levels = c(1, 2, 3),
                   labels = c("Controls", "LM", "LI")),
    comparison = paste0(m$name, "\n Accuracy: ", round(cv_acc, 1), "%  |  BER: ", round(ber, 2)))})

# Plot
library(ggtext)

ggplot(plot_data, aes(x = LD1, y = group, color = group)) +
  geom_jitter(size = 3, height = 0.15, alpha = 0.8) +
  facet_wrap(~ comparison, ncol = 2, scales = "free_y",
             labeller = as_labeller(function(x) paste0(
               sub("\n", "**<br><span style='font-size:13pt'>", paste0("**", x)),
               "</span>"))) +
  scale_color_manual(values = c("Controls" = "#8DA0CB",
                                "LM"       = "#66C2A5",
                                "LI"       = "#E78AC3")) +
  theme_bw() +
  theme(panel.spacing    = unit(12, "pt"),
        legend.position  = "none",
        panel.grid       = element_blank(),
        strip.background = element_blank(),
        strip.text       =  element_markdown(size = 14, lineheight = 1.2),
        axis.text        = element_text(size = 11.8),
        axis.title       = element_text(size = 14),
        plot.title       = element_text(size = 16)) +
  labs(x = "LD1 score", y = "Study group",
       title = paste0("Study group separation by linear discriminant analysis after a lactose challenge (all target compounds, ", TIMEPOINT, ")"))

ggsave(paste0("results/pub_plots/LDA_group_separation_all_", TIMEPOINT, "_with_sex.png"), width = 12, height = 10, dpi = 300, create.dir = TRUE)
