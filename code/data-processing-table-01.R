# Cleans raw data and 

#===================================================================
# Workspace prep.
#===================================================================
# to prevent scientific notation
options(scipen=999)

# Load necessary packages
library(readr)
library(tidyverse)
library(kableExtra)
library(knitr)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtsummary)
library(gt)
library(patchwork)
library(knitcitations)
library(MASS)
library(leaps)
library(cowplot)
library(tsbart)
library(bartMachine)    
library(randomForest)
library(glmnet)
library(caret)

#===================================================================
# Read both files
#===================================================================

natality2023us <- read_csv("data/natality2023us.csv") %>%
  mutate(territory = "U.S.")
dim(natality2023us)

natality2023ps <- read_csv("data/natality2023ps.csv") %>%
  mutate(territory = "U.S. Territory")
head(natality2023ps)

natality2023 = rbind(natality2023us, natality2023ps)

#===================================================================
# Clean and create analysis sample
#===================================================================

train_full_df= natality2023us %>%
  mutate(
    # Recode variables with labels
    preterm = factor(gestrec3,
                     levels = c(1, 2, 3),
                     labels = c("Yes", "No", "Unknown")),
    LBW = ifelse(dbwt <= 2500, 1, 0),
    LBW = factor(LBW, labels = c("No", "Yes")),
    territory = factor(territory,
                       levels = c("U.S.", "U.S. Territory")),
    MRACE6 = factor(mrace6, 
                    levels = 1:6,
                    labels = c("White", "Black", "AIAN", "Asian", "NHOPI", "Multiracial")),
    MAGER9 = factor(mager9,
                    levels = 1:8,
                    labels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54")),
    SEX = factor(sex, 
                 levels = c("M", "F"),
                 labels = c("Male", "Female")),
    RF_PPTERM = factor(rf_ppterm,
                       levels = c("Y", "N", "U"),
                       labels = c("Yes", "No", "Unknown")),
    RF_GDIAB = factor(rf_gdiab,
                      levels = c("Y", "N", "U"),
                      labels = c("Yes", "No", "Unknown")),
    RF_GHYPE = factor(rf_ghype,
                      levels = c("Y", "N", "U"),
                      labels = c("Yes", "No", "Unknown")),
    CIG_REC = factor(cig_rec,
                     levels = c("Y", "N", "U"),
                     labels = c("Yes", "No", "Unknown"))
  ) %>%
  dplyr::select(
    territory,
    MRACE6,
    mager,
    RF_PPTERM,
    RF_GDIAB,
    RF_GHYPE,
    previs,
    wtgain,
    CIG_REC,
    combgest,
    LBW,
    preterm,
    SEX
  ) %>%
  mutate(preterm = as.factor(ifelse(preterm == "Yes", 1, 0)))


test_full_df= natality2023ps %>%
  mutate(
    # Recode variables with labels
    preterm = factor(gestrec3,
                     levels = c(1, 2, 3),
                     labels = c("Yes", "No", "Unknown")),
    LBW = ifelse(dbwt <= 2500, 1, 0),
    LBW = factor(LBW, labels = c("No", "Yes")),
    territory = factor(territory,
                       levels = c("U.S.", "U.S. Territory")),
    MRACE6 = factor(mrace6, 
                    levels = 1:6,
                    labels = c("White", "Black", "AIAN", "Asian", "NHOPI", "Multiracial")),
    MAGER9 = factor(mager9,
                    levels = 1:8,
                    labels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54")),
    SEX = factor(sex, 
                 levels = c("M", "F"),
                 labels = c("Male", "Female")),
    RF_PPTERM = factor(rf_ppterm,
                       levels = c("Y", "N", "U"),
                       labels = c("Yes", "No", "Unknown")),
    RF_GDIAB = factor(rf_gdiab,
                      levels = c("Y", "N", "U"),
                      labels = c("Yes", "No", "Unknown")),
    RF_GHYPE = factor(rf_ghype,
                      levels = c("Y", "N", "U"),
                      labels = c("Yes", "No", "Unknown")),
    CIG_REC = factor(cig_rec,
                     levels = c("Y", "N", "U"),
                     labels = c("Yes", "No", "Unknown"))
  ) %>%
  dplyr::select(
    territory,
    MRACE6,
    mager,
    RF_PPTERM,
    RF_GDIAB,
    RF_GHYPE,
    previs,
    wtgain,
    CIG_REC,
    combgest,
    preterm,
    LBW,
    SEX
  ) %>%
  mutate(preterm = as.factor(ifelse(preterm == "Yes", 1, 0)))

train_full_df <- na.omit(train_full_df) 
saveRDS(train_full_df, "./data/full_mainland.rds")
test_full_df <- na.omit(test_full_df)
saveRDS(test_full_df, "./data/full_territory.rds")

set.seed(64)
train_idx <- sample(nrow(train_full_df), size = 2000)
mainland_sample <- train_full_df[train_idx, ]
mainland_sample$id <- 1:2000
saveRDS(mainland_sample, "./data/mainland_sample.rds")
test_idx = sample(nrow(test_full_df), size = 200)
territory_sample <- test_full_df[test_idx, ]
territory_sample$id <- 1:200
saveRDS(territory_sample, "./data/territory_sample.rds")


#####################################################################
# Assemble table.
#####################################################################

# Prepare data with subcategories and labels
table1_data <- natality2023 %>%
  mutate(
    # Recode variables with labels
    preterm = factor(gestrec3,
                     levels = c(1, 2, 3),
                     labels = c("Yes", "No", "Unknown")),
    LBW = ifelse(dbwt <= 2500, 1, 0),
    LBW = factor(LBW, labels = c("No", "Yes")),
    territory = factor(territory,
                       levels = c("U.S.", "U.S. Territory")),
    MRACE6 = factor(mrace6, 
                    levels = 1:6,
                    labels = c("White", "Black", "AIAN", "Asian", "NHOPI", "Multiracial")),
    MAGER9 = factor(mager9,
                    levels = 1:8,
                    labels = c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54")),
    SEX = factor(sex, 
                 levels = c("M", "F"),
                 labels = c("Male", "Female")),
    RF_PPTERM = factor(rf_ppterm,
                       levels = c("Y", "N", "U"),
                       labels = c("Yes", "No", "Unknown")),
    RF_GDIAB = factor(rf_gdiab,
                      levels = c("Y", "N", "U"),
                      labels = c("Yes", "No", "Unknown")),
    RF_GHYPE = factor(rf_ghype,
                      levels = c("Y", "N", "U"),
                      labels = c("Yes", "No", "Unknown")),
    CIG_REC = factor(cig_rec,
                     levels = c("Y", "N", "U"),
                     labels = c("Yes", "No", "Unknown")),
    # Create subcategory columns
    `Maternal Demographics` = NA,
    `Maternal Risk Factors` = NA,
    `Infant Characteristics` = NA
  ) %>%
  dplyr::select(
    territory,
    `Maternal Demographics`,
    MRACE6,
    MAGER9,
    `Maternal Risk Factors`,
    RF_PPTERM,
    RF_GDIAB,
    RF_GHYPE,
    previs,
    wtgain,
    CIG_REC,
    `Infant Characteristics`,
    LBW,
    preterm,
    SEX
  )

# Generate Table 1 with subcategories
table1 <- table1_data  %>%
  tbl_summary(
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    by = territory,
    digits = all_continuous() ~ 1,
    missing = "no",
    type = list(
      previs ~ "continuous",
      wtgain ~ "continuous"
    ),
    label = list(
      MRACE6 ~ "Mother's race",
      MAGER9 ~ "Mother's age (years)",
      RF_PPTERM ~ "Previous preterm birth",
      RF_GDIAB ~ "Gestational diabetes",
      RF_GHYPE ~ "Gestational hypertension",
      previs ~ "Prenatal visits (count)",
      wtgain ~ "Weight gain (lbs)",
      CIG_REC ~ "Smoked during pregnancy",
      preterm ~ "Preterm birth (<37 weeks)",
      LBW ~ "Low birth weight (<2500 gram)",
      SEX ~ "Infant sex (female)"
    )
  ) %>%
  add_overall() %>%
  #add_p() %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_caption("Table 1. Maternal and Infant Characteristics (CDC NCHS Natality 2023)") %>%
  # for sub-tab purpose
  modify_table_body(
    ~ .x %>%
      mutate(across(everything(), ~ ifelse(. == "0 (NA%)", "", .)))) %>%
  # Add bold headers for subcategories
  modify_table_styling(
    rows = label %in% c("Maternal Demographics", "Maternal Risk Factors", "Infant Characteristics"),
    columns = label,
    text_format = "bold"
  ) 
# modify_table_styling(
#   columns = label,
#   rows = variable %in% c("MRACE6", "MAGER9"),
#   text_format = "bold",
#   footnote = NA_character_  # Explicitly set footnote to character NA
# ) %>%
# modify_table_styling(
#   columns = label,
#   rows = variable %in% c("RF_PPTERM", "RF_GDIAB", "RF_GHYPE", "PREVIS", "WTGAIN", "CIG_REC"),
#   text_format = "bold",
#   footnote = NA_character_
# ) %>%
# modify_table_styling(
#   columns = label,
#   rows = variable == "SEX",
#   text_format = "bold",
#   footnote = NA_character_
# )

# Print as LaTeX
table1 %>%
  as_kable_extra(
    booktabs = TRUE,
    longtable = TRUE,
    linesep = "",
    format = "latex"
  ) %>%
  kableExtra::kable_styling(
    position = "center",
    latex_options = c("striped", "repeat_header", "hold_position", "scale_down"),
    stripe_color = "gray!15",
    font_size = 8
  )

table1 %>%
  as_gt() %>%  # Convert gtsummary table to gt object
  gt::gtsave(
    filename = "table1.png",
    path = getwd(),  # Saves to current working directory
    vwidth = 1200,   # Width in pixels
    vheight = 900,   # Height in pixels
    zoom = 2         # Zoom factor for better resolution
  )

table1 %>%
  as_kable_extra(format = "latex") %>%
  writeLines("table1.tex")
