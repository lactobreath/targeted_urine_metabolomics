### Preparation ----------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(arsenal)
library(tableHTML)
library(writexl)
library(this.path)

## Upload data
setwd(this.dir())
data <- read_excel("urine_targeted_gcms_data_processed.xlsx", sheet = "Sheet1")
dim(data)
names(data)


## Add timepoint column (w/o at-home samples)
data <- data %>%
  mutate(timepoint = case_when(
    grepl("BL", Name) ~ "bl2",
    grepl("0-3", Name) ~ "p1",
    grepl("3-6", Name) ~ "p2",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(timepoint))


## Extract samples, one row per participant x timepoint, keeping only what's needed
samples_fluid <- data %>%
  filter(screen_group %in% c(1, 2, 3) & !is.na(screen_group)) %>%
  distinct(study_id, screen_group, timepoint, urine_b2_ml, urine_p1_ml, urine_p2_ml,
           water_consumed_b2_ml, water_consumed_p1_ml, water_consumed_p2_ml)


## Reshape to one row per participant (urine/water values are already timepoint-specific columns,
## so we just need one row per study_id/screen_group - de-duplicate across timepoints)
data_fluid <- samples_fluid %>%
  dplyr::select(-timepoint) %>%
  distinct() %>%
  mutate(
    screen_group = factor(screen_group,
                          levels = c(1, 2, 3),
                          labels = c("Controls", "LM", "LI"))
  ) %>%
  dplyr::select(
    study_id,
    screen_group,
    urine_b2_ml, urine_p1_ml, urine_p2_ml,
    water_consumed_p1_ml, water_consumed_p2_ml   # water_consumed_b2_ml excluded - constant 650ml, no variance
  ) %>%
  rename(
    "Urine volume, Baseline (mL)" = urine_b2_ml,
    "Urine volume, 0-3h (mL)"     = urine_p1_ml,
    "Urine volume, 3-6h (mL)"     = urine_p2_ml,
    "Fluid intake, 0-3h (mL)"     = water_consumed_p1_ml,
    "Fluid intake, 3-6h (mL)"     = water_consumed_p2_ml
  )

dim(data_fluid)
table(data_fluid$screen_group)


### Set up table -----------------------------------------------------------------
mycontrols <- tableby.control(test = TRUE, total = TRUE,
                              numeric.test = "kwt", cat.test = "chisq", ordered.test = "trend",
                              numeric.stats = c("medianq1q3"),
                              cat.stats = c("countpct"),
                              stats.labels = list(N = 'N', medianq1q3 = 'Median (Q1,Q2)'),
                              digits = 1, digits.p = 4, format.p = FALSE)

make_table <- function(tableby_object, footer_text, col_names) {
  summary(tableby_object, text = TRUE) %>%
    data.frame() %>%
    rename(all_of(col_names)) %>%
    tableHTML(rownames = FALSE,
              footer = footer_text,
              widths = c(200, 150, 150, 150, 150, 100)) %>%
    add_theme(theme = 'scientific') %>%
    add_css_header(css = list(c('height', 'background-color', 'font-size', 'padding'),
                              c('20px', 'white', '13px', '10px')),
                   headers = 1:6) %>%
    add_css_tbody(css = list(c('height', 'background-color', 'font-size', 'padding'),
                             c('18px', 'white', '12px', '10px'))) %>%
    add_css_row(css = list(c('background-color', 'font-size'), c('white', '12px')),
                rows = 1) %>%
    add_css_column(css = list(c('text-align', 'background-color'), c('center', 'white')),
                   columns = c(1:6)) %>%
    add_css_footer(css = list(c('align-text', 'font-weight', 'color', 'font-size'),
                              c('center', 'bold', 'black', 11)))
}

## Build table
fluid_table <- tableby(screen_group ~ .,
                       data = dplyr::select(data_fluid, -study_id),
                       control = mycontrols)

summary(fluid_table, text = TRUE)   # <-- run this first to see actual auto-generated column names before renaming


## Export as Excel
setwd(file.path(this.dir(), "results"))

table_2_df <- summary(fluid_table, text = TRUE) %>% data.frame()
write_xlsx(table_2_df, "table_2a.xlsx")
