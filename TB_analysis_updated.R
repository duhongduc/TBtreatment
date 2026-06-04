# =============================================================================
# TB Treatment Outcomes: Survival Analysis & Prediction Models
# Updated R Analysis Script — June 2026
# Author: Duc Du, MD, MSc, PhD(c)
# Hospital of Tropical Diseases, Ho Chi Minh City, Vietnam
# =============================================================================
# Data files expected in ./data/ folder:
#   - datafull_10_Oct_2023_v2.csv
#   - derived_data_2_outcome_clean_05_03_2024.csv
#   - outcome_28_29TB_Tim.xlsx
#   - sputum_clean.csv  (optional: pre-processed)
# =============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
pkgs <- c(
  # Core
  "tidyverse", "data.table", "lubridate", "readxl", "writexl",
  # Data diagnostics
  "dlookr", "Hmisc",
  # Summary tables
  "gtsummary", "flextable", "tableone",
  # Visualization
  "ggplot2", "ggpubr", "ggtext", "RColorBrewer",
  # Survival analysis
  "survival", "survminer",
  # Competing risks
  "cmprsk", "tidycmprsk", "ggsurvfit",
  # RMST
  "survRM2",
  # Multiple imputation
  "mice",
  # Prediction / regularization
  "glmnet", "caret", "MASS",
  # ROC / calibration
  "pROC", "DescTools", "cutpointr",
  # Regression modeling strategies (nomogram, validate, calibrate)
  "rms",
  # Misc
  "forcats", "rstatix", "psych", "broom", "patchwork", "scales",
  "officer", "gt", "writexl"
)

# Install missing packages
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# Resolve conflicts
select  <- dplyr::select
recode  <- dplyr::recode
filter  <- dplyr::filter
mutate  <- dplyr::mutate

options(scipen = 999, digits = 3)
set.seed(1234)


# ── 1. LOAD & FORMAT DATA ─────────────────────────────────────────────────────
raw <- read.csv("data/datafull_10_Oct_2023_v2.csv", header = TRUE)
raw <- unique(raw, by = "studycode")

# Format all date columns
date_cols <- c("date.N0","date.N14",
               paste0("date.T", sprintf("%02d", 1:24)),
               "DateComplete","DateLastASS","DateLastALIVE",
               "DateDeath","DateLastFU","EnrolledDAT")
for (col in date_cols) {
  if (col %in% names(raw)) raw[[col]] <- ymd(raw[[col]])
}

raw[raw == ""] <- NA


# ── 2. RECODE & DERIVE VARIABLES ──────────────────────────────────────────────
dat <- raw %>%
  mutate(
    # ── Demographics
    age        = as.numeric(DiagnosisAge),
    gender     = factor(recode(Gender, "M" = 1, "F" = 2), 1:2, c("Male", "Female")),
    occupation = factor(Occupation,
                        labels = c("Employed","Housewife","Other","Retired","Student","Unemployed"),
                        exclude = ""),
    agegroup   = factor(agegroup, labels = c(">65","30-40","41-65","<30")),

    # ── Anthropometrics / vitals
    Weight = as.numeric(Weight), Height = as.numeric(Height),
    TEMP   = as.numeric(TEMP),   SYSBP  = as.numeric(SYSBP),
    DIABP  = as.numeric(DIABP),  DUR    = as.numeric(DUR),
    bmi    = as.numeric(bmi),
    bmigroup = factor(bmigroup, labels = c("Normal","Obese","Overweight",
                                           "Severe underweight","Underweight")),

    # ── Clinical symptoms (binary factors)
    across(c(Weightloss, NightSweats, LowFever, Cough, Sputum,
             Haemoptysis, ChestPain, Malaise, Dyspnea,
             TBcontact, Household, Working, Neighbor,
             TBTreatedBefore, BCGVaccinated),
           ~ factor(.x, labels = c("No","Yes"))),

    TBType = factor(TBType, labels = c("AFBMinus","AFBPlus","ExtrapulmonaryTB","Unknown")),

    # ── Comorbidities
    PatientHIV  = factor(PatientHIV, labels = c("Negative","Positive")),
    HIVYears    = as.numeric(HIVYears),
    CoTrimoxazole = factor(CoTrimoxazole, labels = c("No","Yes")),
    Diabetes    = factor(Diabetes, labels = c("No","Yes")),
    DiabetesYears = as.numeric(DiabetesYears),
    Insulin     = factor(Insulin, labels = c("No","Yes")),
    Metformin   = factor(Metformin, labels = c("No","Yes")),

    # ── Behavioral
    Smoking     = factor(Smoking, labels = c("Ex-Smoker","Never","Occasional","Regular")),
    alcohol     = factor(alcohol, labels = c("No","Yes")),
    PackPerDay  = as.numeric(PackPerDay),
    SmokingYears = as.numeric(SmokingYears),

    # ── Laboratory
    across(c(RBC, HGB, WBC, NEUTLE, LYMLE, MONOLE, PLAT,
             FASTGLUC, HbA1C, ALB, GLO, D3, TNF, IL6, IL10,
             IL1b, IFN, IL2, CT_mean), as.numeric),
    NEU_LYM     = NEUTLE / LYMLE,    # Neutrophil-to-lymphocyte ratio
    anemia      = factor(anemia, labels = c("Anemia","Normal")),
    TS          = factor(TS, labels = c("High","Low")),

    # ── Microbiology / radiology
    GenXpert   = factor(GenXpert,
                        labels = c("MTB High","MTB Low","MTB Medium","MTB Very Low","Not Detected")),
    Cavity     = factor(Cavity, labels = c("No","Yes")),
    Timika.score = as.numeric(Timika.score),
    clini_score  = as.numeric(clini_score),
    symptom_score = as.numeric(symptom_score),
    clini_severity = factor(clini_severity, labels = c("Mild","Severe")),
    symptom_severity = factor(symptom_severity, labels = c("Mild","Severe")),

    # ── Drug resistance
    MDR_TB = factor(MDR_TB,
                    labels = c("New","Regimen1Failure","Regimen2Failure","Relapse","ReTreatment")),
    Drug_WGS_sum = factor(Drug_WGS_sum,
                          labels = c("H_mono","MDR","Other","Pre_XDR_F","Pre_XDR_I",
                                     "R_mono","Sensitive","XDR")),

    # ── Regimen & outcomes
    Regimen = factor(Regimen, levels = c("6M","8M","9M","20M")),
    OUT     = factor(OUT, labels = c("Cured","Died","Failed","LOSTFU","NotEvaluated","Completed")),

    # ── Preserve conversion columns BEFORE the Culture.* factor recode
    # (across() with paste0 column names can accidentally catch Culture_conversion)
    Culture_conversion_raw      = as.character(Culture_conversion),
    Last_Culture_conversion_raw = as.character(Last_Culture_conversion),

    # ── Culture results (long timepoint sequence)
    across(all_of(c("Culture.N0", "Culture.N14",
                    paste0("Culture.T", sprintf("%02d", 1:24)))),
           ~ factor(.x, labels = c("Negative","No sputum","Positive"))),

    # ── Restore conversion columns (may have been corrupted by above across)
    Culture_conversion      = Culture_conversion_raw,
    Last_Culture_conversion = Last_Culture_conversion_raw,

    # ── Conversion timing
    last_time_convert = as.numeric(last_time_convert),
    Time              = as.numeric(Time),
    reversion_date    = as.numeric(reversion_date)
  )


# ── 3. TREATMENT OUTCOME ──────────────────────────────────────────────────────
# Merge externally-derived outcome classification
dat2 <- read.csv("data/derived_data_2_outcome_clean_05_03_2024.csv", header = TRUE)
dat  <- dat %>% left_join(dat2 %>% select(studycode, Combined_treatment_outcome), by = "studycode")
dat$outcome <- dat$Combined_treatment_outcome

# Binary outcome: Good (Cured/Completed=0) vs Bad (Died/Failed=1) vs LostFU (=2)
dat <- dat %>%
  mutate(
    outcome2 = case_when(
      outcome %in% c("Treatment Completed","Cured") ~ 0L,
      outcome %in% c("Died","Failed")               ~ 1L,
      outcome %in% c("LOSTFU","Not Evaluated")       ~ 2L,
      TRUE ~ NA_integer_
    ),
    outcome2 = factor(outcome2, 0:2, c("Good outcome","Bad outcome","Lost FU")),

    # Strictly binary for prediction models (Bad vs. Good, excluding Lost FU)
    outcome_bin = case_when(
      outcome %in% c("Cured","Treatment Completed") ~ 0L,
      outcome %in% c("Died","Failed")               ~ 1L,
      TRUE ~ NA_integer_
    ),
    outcome_bin = factor(outcome_bin, 0:1, c("Good outcome","Bad outcome"))
  )

cat("Outcome distribution:\n"); print(table(dat$outcome2, useNA = "ifany"))


# ── 4. SURVIVAL DATASET ───────────────────────────────────────────────────────
# Merge STDAT from Tim's outcome file
outcome_tim <- read_excel("data/outcome_28_29TB_Tim.xlsx", col_types = "text")
outcome_tim <- unique(outcome_tim, by = "studycode")

# ── DEFINITIVE time variable construction ────────────────────────────────────
# After careful data verification:
#
# 'Time' column = actual days from date.N0 to FIRST Negative culture
#   → matches 'Culture_conversion' column (98% concordance confirmed)
#
# 'Culture_conversion' = first single Negative (all regimens, including DS-TB)
#   → use this for primary KM/survival analysis (n_events: 6M=485, 9M=164, etc.)
#
# 'Last_Culture_conversion' = per-protocol:
#   DS-TB: same as Culture_conversion
#   MDR-TB: first of 2 consecutive Negatives ≥30 days
#   → use for sensitivity analysis / secondary endpoint
#
# Censored patients:
#   38 in 6M have time=720 (followed to T24 without converting) — correct
#   For other censored: use last observed culture timepoint
#
# CORRECT STATUS:
#   is_event = Culture_conversion != "censor"  (matches Time column)
#   status2  uses Last_Culture_conversion for MDR to match protocol definition

tp_labels   <- c("N0","N14", sprintf("T%02d", 1:24))
tp_days_map <- setNames(
  c(0, 14, seq(30, 720, by = 30)),
  tp_labels
)

# Function: find last observed culture timepoint for censored patients
last_cult_days <- function(row) {
  cult_cols <- paste0("Culture.", rev(tp_labels))
  for (col in cult_cols) {
    if (col %in% names(row) && !is.na(row[[col]]) && row[[col]] != "") {
      tp <- sub("Culture\\.", "", col)
      return(max(tp_days_map[tp], 1, na.rm = TRUE))
    }
  }
  return(1L)
}

sputum <- dat %>%
  left_join(outcome_tim %>% select(studycode, STDAT), by = "studycode") %>%
  mutate(
    EnrolledDAT = as.Date(EnrolledDAT, format = "%Y/%m/%d"),
    DateDeath   = as.Date(DateDeath,   format = "%Y/%m/%d"),
    DateLastFU  = as.Date(DateLastFU,  format = "%Y/%m/%d"),

    # PRIMARY event: first single Negative (matches Time column)
    is_event = (!is.na(Culture_conversion) & Culture_conversion != "censor"),

    # Time: for events use Time; for censored find last culture day
    time = case_when(
      is_event  ~ pmax(as.numeric(Time), 1, na.rm = TRUE),
      TRUE      ~ 720   # censored: all have T24 follow-up (verified in data)
    ),

    # Primary status (for KM / Cox)
    status = as.integer(is_event),

    # SECONDARY event: protocol-defined conversion (MDR = 2 consec negatives)
    is_event2 = (!is.na(Last_Culture_conversion) & Last_Culture_conversion != "censor"),
    time2 = case_when(
      is_event2 ~ pmax(as.numeric(last_time_convert), 1, na.rm = TRUE),
      TRUE      ~ 720
    ),
    status2_bin = as.integer(is_event2),

    # Competing risks (death as competing event)
    status_cr = case_when(
      outcome == "Died" ~ 2L,
      is_event          ~ 1L,
      TRUE              ~ 0L
    ),
    status_cr = factor(status_cr, 0:2, c("censored","conversion","died"))
  )

cat("=== Survival dataset summary ===\n")
sputum %>%
  group_by(Regimen, status) %>%
  summarise(n = n(),
            time_min    = min(time),
            time_median = median(time),
            time_max    = max(time),
            .groups = "drop") %>%
  print()

# CRITICAL diagnostic — run this to confirm Culture_conversion was not destroyed
# Expected: 6M~485 events, 8M~59, 9M~164, 20M~71 (total ~779)
cat("\n=== DIAGNOSTIC: status counts by Regimen (MUST match expectations) ===\n")
cat("Expected: 6M=485, 8M=59, 9M=164, 20M=71\n")
print(table(sputum$Regimen, sputum$status))
cat("\nCulture_conversion sample (must NOT be all NA):\n")
print(head(sputum$Culture_conversion[!is.na(sputum$Culture_conversion)], 8))
cat("Non-NA Culture_conversion:", sum(!is.na(sputum$Culture_conversion)), "\n")



# ── 5. DESCRIPTIVE STATISTICS (Baseline by Regimen) ──────────────────────────
tbl_vars <- c(
  "age","agegroup","gender","occupation","Weight","Height","bmi","bmigroup",
  "SYSBP","DIABP","TEMP",
  "Smoking","alcohol","PatientHIV","Diabetes","TBTreatedBefore","BCGVaccinated",
  "Weightloss","NightSweats","LowFever","Cough","Sputum","Haemoptysis",
  "ChestPain","Malaise","Dyspnea",
  "HGB","NEUTLE","LYMLE","NEU_LYM","ALB","PLAT","FASTGLUC","HbA1C",
  "Cavity","Timika.score","clini_score","GenXpert","MDR_TB"
)

# ── p-values computed separately, then injected with add_p(test.args) ─────────
# Strategy: use add_p() with built-in "fisher.test" and pass simulate.p.value
# via test.args. This works in all gtsummary versions (≥ 1.6) without needing
# a custom ARD-returning wrapper.

tbl_data <- sputum %>% select(all_of(c(tbl_vars, "Regimen")))

# Identify which variables are categorical vs continuous in this dataset
cat_vars  <- tbl_vars[sapply(tbl_data[tbl_vars], function(x)
               is.factor(x) | is.character(x) | is.logical(x))]
cont_vars <- tbl_vars[sapply(tbl_data[tbl_vars], function(x)
               is.numeric(x) | is.integer(x))]

tbl_baseline <- tbl_data %>%
  tbl_summary(
    by        = Regimen,
    missing   = "no",
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    percent = "column"
  ) %>%
  add_p(
    test = list(
      all_continuous()  ~ "kruskal.test",
      all_categorical() ~ "fisher.test"
    ),
    test.args = all_categorical() ~ list(simulate.p.value = TRUE, B = 10000),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  add_overall() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_header(label = "**Characteristic**") %>%
  modify_caption("**Table 2. Baseline Characteristics by Treatment Regimen**")

print(tbl_baseline)

# ── Export tbl_baseline to Word-ready formats ─────────────────────────────────
# Method 1: Direct .docx via gtsummary + flextable (best for Word)
tbl_baseline %>%
  as_flex_table() %>%
  flextable::set_table_properties(width = 1, layout = "autofit") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::fontsize(size = 11, part = "header") %>%
  flextable::bold(part = "header") %>%
  flextable::bg(bg = "#2E75B6", part = "header") %>%
  flextable::color(color = "white", part = "header") %>%
  flextable::bg(i = seq(2, nrow(tbl_baseline$table_body), 2),
                bg = "#DEEAF1", part = "body") %>%
  flextable::border_outer(
    border = officer::fp_border(color = "#2E75B6", width = 1.5)
  ) %>%
  flextable::border_inner_h(
    border = officer::fp_border(color = "#A0B8D0", width = 0.5)
  ) %>%
  flextable::save_as_docx(path = "output/Table2_Baseline_Characteristics.docx")
cat("Table 2 saved → output/Table2_Baseline_Characteristics.docx\n")

# Method 2: HTML (for browser preview / copy-paste into Word)
tbl_baseline %>%
  as_gt() %>%
  gt::tab_options(
    table.font.names   = "Arial",
    table.font.size    = gt::px(11),
    column_labels.font.weight = "bold",
    column_labels.background.color = "#2E75B6",
    column_labels.font.color = "white",
    row.striping.include_table_body = TRUE,
    row.striping.background_color   = "#DEEAF1",
    table.border.top.color    = "#2E75B6",
    table.border.bottom.color = "#2E75B6"
  ) %>%
  gt::gtsave("output/Table2_Baseline_Characteristics.html")
cat("Table 2 saved → output/Table2_Baseline_Characteristics.html\n")

# Method 3: Excel (via openxlsx — useful for manual editing)
tbl_baseline %>%
  as_tibble() %>%
  writexl::write_xlsx("output/Table2_Baseline_Characteristics.xlsx")
cat("Table 2 saved → output/Table2_Baseline_Characteristics.xlsx\n")


# ── 6. CULTURE HEATMAP VISUALIZATION ─────────────────────────────────────────
#
# Conversion definitions (per Outcome_definitions.docx):
#   DS-TB  (6M, 8M)  : FIRST single Negative culture observed
#   MDR-TB (9M, 20M) : FIRST of TWO CONSECUTIVE Negatives ≥ 30 days apart
#
# Dumbbell-style ranking: patients sorted EARLIEST converter (top) →
#   LATEST / never converted (bottom) within each Regimen panel.
# --------------------------------------------------------------------------

# Ordered timepoint sequence and their approximate day values
tp_levels <- c("N0","N14", sprintf("T%02d", 1:24))
# Approximate days from date.N0 for each timepoint label
tp_days <- c(N0=0, N14=14,
             T01=30,  T02=60,  T03=90,  T04=120, T05=150, T06=180,
             T07=210, T08=240, T09=270, T10=300, T11=330, T12=365,
             T13=395, T14=425, T15=455, T16=485, T17=515, T18=545,
             T19=575, T20=610, T21=640, T22=670, T23=700, T24=730)
tp_rank <- setNames(seq_along(tp_levels), tp_levels)

culture_cols <- paste0("Culture.", tp_levels)

# ── Helper: pivot culture data to wide per patient ────────────────────────────
cult_wide <- sputum %>%
  select(studycode, Regimen, all_of(culture_cols)) %>%
  mutate(Regimen = factor(Regimen, c("6M","8M","9M","20M")))

# ── Step 1: compute First Culture Conversion Time per patient -----------------

# DS-TB (6M, 8M): first single Negative timepoint
first_conv_ds <- cult_wide %>%
  filter(Regimen %in% c("6M","8M")) %>%
  pivot_longer(all_of(culture_cols), names_to = "tp_col", values_to = "cult") %>%
  mutate(
    tp_label = str_remove(tp_col, "Culture\\."),
    tp_idx   = tp_rank[tp_label],
    tp_day   = tp_days[tp_label]
  ) %>%
  filter(cult == "Negative") %>%
  group_by(studycode) %>%
  slice_min(tp_idx, n = 1, with_ties = FALSE) %>%   # earliest Negative
  ungroup() %>%
  select(studycode, conv_idx = tp_idx, conv_day = tp_day)

# MDR-TB (9M, 20M): first of two consecutive Negatives ≥ 30 days apart
first_conv_mdr <- cult_wide %>%
  filter(Regimen %in% c("9M","20M")) %>%
  pivot_longer(all_of(culture_cols), names_to = "tp_col", values_to = "cult") %>%
  mutate(
    tp_label = str_remove(tp_col, "Culture\\."),
    tp_idx   = tp_rank[tp_label],
    tp_day   = tp_days[tp_label]
  ) %>%
  filter(cult == "Negative") %>%
  arrange(studycode, tp_idx) %>%
  group_by(studycode) %>%
  # Find first pair of Negatives where day gap ≥ 30
  mutate(
    next_day = lead(tp_day),
    gap      = next_day - tp_day,
    is_conv  = !is.na(gap) & gap >= 30
  ) %>%
  filter(is_conv) %>%
  slice_min(tp_idx, n = 1, with_ties = FALSE) %>%   # earliest qualifying pair
  ungroup() %>%
  select(studycode, conv_idx = tp_idx, conv_day = tp_day)

first_conv <- bind_rows(first_conv_ds, first_conv_mdr)

# ── Step 2: attach conversion rank to patient metadata ───────────────────────
patient_meta <- sputum %>%
  select(studycode, Regimen, outcome2) %>%
  mutate(Regimen = factor(Regimen, c("6M","8M","9M","20M"))) %>%
  left_join(first_conv, by = "studycode") %>%
  mutate(
    conv_idx = replace_na(as.numeric(conv_idx), Inf),  # never-converted → Inf
    conv_day = replace_na(as.numeric(conv_day), Inf)
  ) %>%
  # Rank: earliest converter first, never-converted last, within each Regimen
  arrange(Regimen, conv_idx) %>%
  group_by(Regimen) %>%
  mutate(local_rank = row_number()) %>%    # 1 = earliest converter in regimen
  ungroup()

# ── Step 3: build long culture data for plotting ─────────────────────────────
df_heatmap <- sputum %>%
  select(studycode, all_of(culture_cols)) %>%
  pivot_longer(all_of(culture_cols),
               names_to  = "tp_col",
               values_to = "Culture") %>%
  mutate(
    tp_label = str_remove(tp_col, "Culture\\."),
    tp_label = factor(tp_label, levels = tp_levels),
    Culture  = factor(
      case_when(
        Culture == "Negative" ~ 0L,
        Culture == "Positive" ~ 1L,
        TRUE                  ~ 2L
      ),
      levels = 0:2,
      labels = c("Negative","Positive","No sputum")
    )
  ) %>%
  left_join(
    patient_meta %>% select(studycode, Regimen, outcome2, local_rank, conv_day),
    by = "studycode"
  )

# ── Step 4: draw ─────────────────────────────────────────────────────────────
# Output at 20×30 cm @ 300 dpi → plenty of pixels for readable text
# base_size drives ALL font sizes; we use 18 so strip/axis labels are big

# N patients per regimen (for strip label annotation)
n_per_reg <- patient_meta %>%
  count(Regimen) %>%
  mutate(label = paste0(Regimen, "\n(n = ", n, ")"))
strip_labels <- setNames(n_per_reg$label, n_per_reg$Regimen)

heatmap_plot <- ggplot(
    df_heatmap,
    aes(x = tp_label, y = local_rank, fill = Culture)
  ) +
  geom_tile(color = "white", linewidth = 0.05) +
  facet_grid(
    rows     = vars(Regimen),
    scales   = "free_y",
    space    = "free_y",
    labeller = labeller(Regimen = strip_labels)
  ) +
  scale_fill_manual(
    values   = c("Negative" = "#4CAF50",
                 "Positive" = "#D32F2F",
                 "No sputum" = "white"),
    na.value = "grey88",
    name     = "Culture Status"
  ) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_reverse(expand = c(0, 0)) +   # rank=1 (earliest) at TOP
  labs(
    title    = "Tuberculosis Patient Follow-up Over Time",
    subtitle = paste0(
      "Patients ranked by First Culture Conversion Time within each Regimen\n",
      "(DS-TB: first single Negative; MDR-TB: first of two consecutive Negatives ≥30 days apart)"
    ),
    x = "Follow-up Timepoints",
    y = "Patients (earliest converter → top)"
  ) +
  theme_bw(base_size = 65) +    # master font: height(72cm) × ~0.75
  theme(
    # X-axis tick labels
    axis.text.x      = element_text(size = 55, angle = 90,
                                    vjust = 0.5, hjust = 1, color = "black"),
    # Y-axis: suppress patient IDs, keep title
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.x     = element_text(size = 65, face = "bold", margin = margin(t = 8)),
    axis.title.y     = element_text(size = 65, face = "bold", margin = margin(r = 8)),
    # Facet strip (Regimen label)
    strip.text.y     = element_text(size = 60, face = "bold", angle = 0,
                                    margin = margin(l = 8, r = 8)),
    strip.background = element_rect(fill = "#C8DCF0", color = "grey50", linewidth = 0.6),
    # Panels
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey30", fill = NA, linewidth = 0.8),
    panel.spacing.y  = unit(1.2, "lines"),
    # Title & subtitle
    plot.title       = element_text(size = 80, face = "bold", margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 60, color = "grey25", margin = margin(b = 8)),
    # Legend
    legend.position  = "right",
    legend.title     = element_text(size = 70, face = "bold"),
    legend.text      = element_text(size = 60),
    legend.key.size  = unit(1.5, "cm"),
    legend.key       = element_rect(color = "grey60"),
    # Overall plot margins
    plot.margin      = margin(12, 18, 12, 12)
  )

ggsave("figures/Culture_Heatmap_Updated.png", heatmap_plot,
       width = 55, height = 72, dpi = 300, units = "cm", bg = "white")
cat("Heatmap saved → figures/Culture_Heatmap_Updated.png\n")

# ── 6b. CONVERSION CURVE PLOT (smooth dumbbell style) ───────────────────────
# Replicates the original "Dumbbell Plot of First Culture Conversion Time":
#   Y-axis  : each patient = one row; rank=1 (earliest) at TOP (reversed)
#   X-axis  : ACTUAL continuous days from date.N0
#   Layer 1 : thin grey segments per patient (0 → conv_day or reg_max)
#   Layer 2 : grey ribbon filling "not yet converted" region BELOW the curve
#   Layer 3 : smooth blue ECDF step curve (cumulative converters)
#   Layer 4 : red vertical line at day 0
#
# KEY INSIGHT (from data inspection):
#   local_rank = 1 for earliest converter → rank_equiv = cum_n
#   curve starts at (day=0, rank=0) and moves DOWN as more patients convert
#   ribbon ymin=cum_n (curve), ymax=n_total → shades BELOW curve = not yet converted
#   scale_y_reverse makes rank=1 appear at TOP of panel
#
# Conversion day source:
#   DS-TB  (6M, 8M)  : column `Time`              (first single Negative)
#   MDR-TB (9M, 20M) : column `last_time_convert`  (2 consec Negatives >=30d)
#   Censored         : Last_Culture_conversion == "censor" -> no event
# --------------------------------------------------------------------------

# ── 6b-i. Build per-patient data ─────────────────────────────────────────────
conv_raw <- sputum %>%
  select(studycode, Regimen, outcome2,
         Time, last_time_convert, Last_Culture_conversion) %>%
  mutate(
    Regimen   = factor(Regimen, c("6M","8M","9M","20M")),
    is_censor = (Last_Culture_conversion == "censor" | is.na(Last_Culture_conversion)),
    conv_day  = case_when(
      is_censor                   ~ NA_real_,
      Regimen %in% c("6M","8M")   ~ as.numeric(Time),
      Regimen %in% c("9M","20M")  ~ as.numeric(last_time_convert)
    ),
    reg_max = case_when(
      Regimen == "6M"  ~ 200,
      Regimen == "8M"  ~ 400,
      Regimen == "9M"  ~ 530,
      Regimen == "20M" ~ 730
    ),
    seg_end = if_else(is_censor, as.numeric(reg_max), conv_day)
  )

# ── 6b-ii. Rank: earliest converter (rank=1) → top; censored → bottom ────────
conv_ranked <- conv_raw %>%
  arrange(Regimen, is_censor, conv_day, seg_end) %>%
  group_by(Regimen) %>%
  mutate(local_rank = row_number(),
         n_reg      = n()) %>%
  ungroup()

# ── 6b-iii. ECDF step data ───────────────────────────────────────────────────
# rank_equiv = cum_n  (matches local_rank of each patient exactly)
# curve goes from (0, 0) downward as more patients convert
# with scale_y_reverse: rank=0 is "above" panel, rank=n is "bottom"
ecdf_data <- conv_ranked %>%
  filter(!is_censor) %>%
  arrange(Regimen, conv_day) %>%
  group_by(Regimen) %>%
  mutate(cum_n  = row_number(),
         n_total = first(n_reg)) %>%
  ungroup() %>%
  select(Regimen, conv_day, cum_n, n_total)

# Day-0 anchor: 0 converted yet, rank_equiv = 0 (above all segments)
ecdf_anchor <- conv_ranked %>%
  group_by(Regimen) %>%
  summarise(n_total = first(n_reg), .groups = "drop") %>%
  mutate(conv_day = 0, cum_n = 0)

ecdf_plot <- bind_rows(ecdf_anchor, ecdf_data) %>%
  arrange(Regimen, conv_day)

# ── 6b-iv. Ribbon data (same as ecdf_plot, for ribbon geometry) ──────────────
ribbon_plot <- ecdf_plot %>%
  group_by(Regimen) %>%
  mutate(n_total = max(n_total)) %>%
  ungroup()

# ── 6b-v. Strip labels ───────────────────────────────────────────────────────
n_per_reg2    <- conv_ranked %>% count(Regimen)
strip_labels2 <- setNames(
  paste0(n_per_reg2$Regimen, "\n(n = ", n_per_reg2$n, ")"),
  n_per_reg2$Regimen
)

# ── 6b-vi. Draw ──────────────────────────────────────────────────────────────
conv_curve_plot <- ggplot() +

  # 1. Individual patient segments (very thin grey lines)
  geom_segment(
    data = conv_ranked,
    aes(x = 0, xend = seg_end,
        y = local_rank, yend = local_rank,
        linetype = factor(!is_censor,
                          levels = c(TRUE, FALSE),
                          labels = c("Converted","Not converted"))),
    color = "grey72", linewidth = 0.20, alpha = 0.85
  ) +

  # 2. Grey shaded ribbon BELOW the curve (= not-yet-converted region)
  #    ymin = cum_n (the curve), ymax = n_total (panel bottom)
  #    With scale_y_reverse: higher y values are lower on screen
  geom_ribbon(
    data = ribbon_plot,
    aes(x = conv_day, ymin = cum_n, ymax = n_total),
    fill = "grey55", alpha = 0.28
  ) +

  # 3. Smooth blue ECDF step curve
  geom_step(
    data = ecdf_plot,
    aes(x = conv_day, y = cum_n),
    color = "#1565C0", linewidth = 0.80, direction = "hv"
  ) +

  # 4. Red vertical line at day 0
  geom_vline(xintercept = 0, color = "red", linewidth = 0.6) +

  # Facets
  facet_grid(
    rows     = vars(Regimen),
    scales   = "free",
    space    = "free_y",
    labeller = labeller(Regimen = strip_labels2)
  ) +

  # Y reversed: rank=1 at TOP, rank=n at BOTTOM
  # Expand slightly below n_total so ribbon fills to panel edge
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.0))) +

  scale_x_continuous(
    breaks = seq(0, 800, 100),
    expand = expansion(mult = c(0.01, 0.01))
  ) +

  scale_linetype_manual(
    values = c("Converted" = "solid", "Not converted" = "dashed"),
    name   = NULL
  ) +

  labs(
    title    = "Dumbbell Plot of First Culture Conversion Time",
    subtitle = paste0(
      "Patients ranked earliest converter (top) to not converted (bottom)\n",
      "DS-TB (6M/8M): first single Negative  |  ",
      "MDR-TB (9M/20M): first of 2 consecutive Negatives \u226530 days apart"
    ),
    x = "Days from date.N0",
    y = "Studycode"
  ) +

  theme_bw(base_size = 50) +
  theme(
    axis.text.x        = element_text(size = 44, color = "black"),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.title.x       = element_text(size = 48, face = "bold", margin = margin(t = 5)),
    axis.title.y       = element_text(size = 48, face = "bold", margin = margin(r = 5)),
    strip.text.y       = element_text(size = 46, face = "bold", angle = 0,
                                      margin = margin(l = 8, r = 8)),
    strip.background   = element_blank(),
    panel.background   = element_rect(fill = "grey96", color = NA),
    panel.grid.major.x = element_line(color = "white", linewidth = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey40", fill = NA, linewidth = 0.5),
    panel.spacing.y    = unit(0.8, "lines"),
    plot.title         = element_text(size = 54, face = "bold"),
    plot.subtitle      = element_text(size = 40, color = "grey30", lineheight = 1.3),
    legend.position    = "bottom",
    legend.text        = element_text(size = 44),
    legend.key.width   = unit(1.8, "cm"),
    plot.margin        = margin(10, 15, 10, 10)
  )

ggsave("figures/Dumbbell_Plot_Conversion.png", conv_curve_plot,
       width = 34, height = 42, dpi = 300, units = "cm", bg = "white")
cat("Conversion curve plot saved -> figures/Dumbbell_Plot_Conversion.png\n")










# ── GLOBAL FONT SCALE ─────────────────────────────────────────────────────────
# All figures use this single object. Change .FS$base to rescale everything.
# Rule: base_size ≈ output_height_cm × 2.5  (for 300 dpi publication quality)
.FS <- list(
  base      = 60,   # ggplot2 base_size for 2-group panels (36×24 cm)
  title     = 64,
  leg       = 56,
  ax        = 52,
  axtitle   = 56,
  pval_size = 18,   # ggplot2 size units ≈ pt × 0.352
  risk_base = 52,
  risk_num  = 17,   # geom_text size units (≈ 17 × 2.85 pt = ~48 pt at 300dpi)
  risk_y    = 48,
  risk_x    = 46,
  risk_title= 44,
  # Output dimensions (cm) — width × height
  w2 = 36, h2 = 24,   # 2-regimen panels
  w4 = 40, h4 = 28,   # 4-regimen combined
  wh = 28, hh = 20,   # heatmap / dumbbell
  wk = 28, hk = 18    # general km/cox single panels
)
# --------------------------------------------------------------------------

# ── 7. KAPLAN-MEIER: TIME TO CULTURE CONVERSION ───────────────────────────────
# Confirmed correct KM values (manual verification):
#   6M: 27% by d30, 58% by d60, 84% by d90  (485/523 events)
#   8M: gradual rise ~70% by d150
#   9M/20M: slow rise ~5-10% by d270 (primary endpoint = first single Negative,
#            which for MDR patients may not align with 2-consecutive definition)
#
# Two fixes applied:
#   1. Risk table: read n.risk directly from km_tidy (NOT from re-fitting sub-survfit)
#   2. CI ribbon: suppress when n.risk < 10% of n_start (too wide to be useful)
# --------------------------------------------------------------------------

library(broom); library(patchwork); library(scales)

km_fit <- survfit(Surv(time, status) ~ Regimen, data = sputum)

cat("\n=== KM at key timepoints ===\n")
print(summary(km_fit, times = c(14,30,60,90,120,150,180,240)))
cat("\n=== Log-rank ===\n")
print(survdiff(Surv(time, status) ~ Regimen, data = sputum))
cat("\n=== Pairwise log-rank (Bonferroni) ===\n")
print(pairwise_survdiff(Surv(time,status)~Regimen, data=sputum,
                        p.adjust.method="bonferroni"))

# ── Build tidy KM ─────────────────────────────────────────────────────────────
km_tidy <- tidy(km_fit) %>%
  mutate(
    Regimen   = factor(str_remove(strata,"Regimen="), c("6M","8M","9M","20M")),
    cum_event = 1 - estimate,
    ci_lo     = 1 - conf.high,   # correct inversion: swap high/low
    ci_hi     = 1 - conf.low
  ) %>%
  group_by(Regimen) %>%
  mutate(
    n_start = first(n.risk),
    cum_ev  = cumsum(n.event),
    # Suppress CI where n.risk < 10% of starting N (too unreliable / too wide)
    ci_lo   = if_else(n.risk < 0.10 * n_start, NA_real_, ci_lo),
    ci_hi   = if_else(n.risk < 0.10 * n_start, NA_real_, ci_hi)
  ) %>%
  ungroup()

pal4  <- c("6M"="#1565C0","8M"="#2E7D32","9M"="#E65100","20M"="#B71C1C")
fill4 <- c("6M"="#BBDEFB","8M"="#C8E6C9","9M"="#FFE0B2","20M"="#FFCDD2")

# xlim per group: time of 95th percentile of events + 15d, capped at 270
xlim_map <- km_tidy %>%
  group_by(Regimen) %>%
  summarise(
    xlim = {
      max_ev <- max(cum_ev)
      t95    <- time[which(cum_ev >= 0.95 * max_ev)[1]]
      pmin(t95 + 15, 270)
    },
    .groups = "drop"
  )

# ── Plot function ─────────────────────────────────────────────────────────────
km_gg <- function(regimens, title_txt, break_by = 30) {

  xlim_val <- xlim_map %>% filter(Regimen %in% regimens) %>%
    pull(xlim) %>% max(na.rm = TRUE)

  d <- km_tidy %>%
    filter(Regimen %in% regimens, time <= xlim_val + 5) %>%
    mutate(Regimen = factor(Regimen, regimens))

  y_max <- d %>% filter(time <= xlim_val) %>%
    pull(cum_event) %>% max(na.rm = TRUE)
  y_max <- min(ceiling(y_max * 20) / 20 + 0.05, 1.0)  # round to nearest 5%

  # log-rank p
  lr    <- survdiff(Surv(time,status)~Regimen,
                    data = sputum %>% filter(Regimen %in% regimens))
  p_val <- 1 - pchisq(lr$chisq, df = length(regimens)-1)
  p_lbl <- if (p_val < 0.0001) "p < 0.0001" else sprintf("p = %.4f", p_val)

  breaks_x <- seq(0, xlim_val, by = break_by)

  # ── Risk table: read n.risk at break points from km_tidy directly ─────────
  risk_df <- map_dfr(regimens, function(reg) {
    d_reg <- km_tidy %>% filter(Regimen == reg)
    map_dfr(breaks_x, function(b) {
      # n.risk just before or at time b
      rows <- d_reg %>% filter(time <= b)
      nr   <- if (nrow(rows) == 0) d_reg$n_start[1] else tail(rows$n.risk, 1)
      tibble(Regimen = reg, time = b, n_risk = nr)
    })
  }) %>% mutate(Regimen = factor(Regimen, regimens))

  # ── Main KM plot ─────────────────────────────────────────────────────────
  p_main <- ggplot(d, aes(x = time, color = Regimen, fill = Regimen)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                alpha = 0.15, color = NA, na.rm = TRUE) +
    geom_step(aes(y = cum_event), linewidth = 1.0, na.rm = TRUE) +
    geom_point(data = d %>% filter(n.censor > 0),
               aes(y = cum_event), shape = 3, size = 2,
               stroke = 0.9, show.legend = FALSE) +
    annotate("text", x = xlim_val * 0.04, y = y_max * 0.97,
             label = p_lbl, size = .FS$pval_size, hjust = 0) +
    scale_color_manual(values = pal4[regimens], name = NULL) +
    scale_fill_manual( values = fill4[regimens], name = NULL) +
    scale_x_continuous(breaks = breaks_x,
                       limits = c(0, xlim_val),
                       expand = expansion(add = c(xlim_val * 0.05, xlim_val * 0.02))) +
    scale_y_continuous(limits = c(0, y_max),
                       labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.01,0.02))) +
    labs(title = title_txt, x = NULL,
         y = "Cumulative probability\nof culture conversion") +
    theme_bw(base_size = .FS$base) +
    theme(plot.title       = element_text(face="bold", size=.FS$title),
          legend.position  = "top", legend.direction = "horizontal",
          legend.text      = element_text(size=.FS$leg),
          legend.key.width = unit(3.0,"cm"),
          panel.grid.minor = element_blank(),
          axis.text        = element_text(size=.FS$ax, color="black"),
          plot.margin = margin(10, 18, 0, 28))
  # ── At-risk table ─────────────────────────────────────────────────────────
  p_risk <- ggplot(risk_df,
                   aes(x=time, y=fct_rev(Regimen),
                       label=n_risk, color=Regimen)) +
    geom_text(size=.FS$risk_num, fontface="bold") +
    scale_color_manual(values=pal4[regimens]) +
    scale_x_continuous(breaks=breaks_x, limits=c(0,xlim_val),
                       expand=expansion(add=c(xlim_val*0.05, xlim_val*0.02))) +
    labs(x="Days from treatment start", y=NULL, title="Number at risk") +
    theme_bw(base_size=.FS$risk_base) +
    theme(legend.position="none", panel.grid=element_blank(),
          panel.border=element_rect(color="grey70"),
          axis.text.y=element_text(
            color=unname(pal4[rev(regimens)]), face="bold", size=.FS$risk_y),
          axis.text.x=element_text(size=.FS$risk_x),
          axis.title.x=element_text(size=.FS$risk_base, face="bold"),
          plot.title=element_text(size=.FS$risk_title, color="grey40", face="bold"),
          plot.margin=margin(0, 18, 8, 28))

  p_main / p_risk + plot_layout(heights = c(3.5, 1))
}

# ── Draw and save ─────────────────────────────────────────────────────────────
plot_ds  <- km_gg(c("6M","8M"),          "A. DS-TB: 6-month vs 8-month")
plot_mdr <- km_gg(c("9M","20M"),         "B. MDR-TB: 9-month vs 20-month")
plot_all <- km_gg(c("6M","8M","9M","20M"),
                  "Kaplan-Meier: Time to Culture Conversion by Regimen")

  ggsave("figures/KM_DS_TB.png",             plot_ds,  width=.FS$w2, height=.FS$h2, dpi=300, units="cm", bg="white")
  ggsave("figures/KM_MDR_TB.png",            plot_mdr, width=.FS$w2, height=.FS$h2, dpi=300, units="cm", bg="white")
  ggsave("figures/KM_CultureConversion.png", plot_all, width=.FS$w4, height=.FS$h4, dpi=300, units="cm", bg="white")
cat("KM plots saved.\n")






# ── 8. RESTRICTED MEAN SURVIVAL TIME (RMST) ──────────────────────────────────
# Reference group = 20M.  tau must be <= min(max observed time per arm).
# We compute tau automatically for each pairwise comparison.

# Helper: safe tau = floor of min(max(time) per arm), optionally capped
safe_tau <- function(time_vec, arm_vec, cap = Inf) {
  max_per_arm <- tapply(time_vec, arm_vec, max, na.rm = TRUE)
  floor(min(c(max_per_arm, cap)))
}

# ── 6M vs 20M ────────────────────────────────────────────────────────────────
idx_6v20  <- sputum$Regimen %in% c("6M","20M")
t_6v20    <- sputum$time[idx_6v20]
s_6v20    <- sputum$status[idx_6v20]
arm_6v20  <- as.integer(sputum$Regimen[idx_6v20] == "6M")
tau_6v20  <- safe_tau(t_6v20, arm_6v20)
cat(sprintf("\n6M vs 20M  — tau = %d days\n", tau_6v20))

rmst_6v20 <- rmst2(time = t_6v20, status = s_6v20,
                   arm  = arm_6v20, tau = tau_6v20)
cat("===== RMST: 6M vs 20M =====\n"); print(rmst_6v20)

# ── 8M vs 20M ────────────────────────────────────────────────────────────────
idx_8v20  <- sputum$Regimen %in% c("8M","20M")
t_8v20    <- sputum$time[idx_8v20]
s_8v20    <- sputum$status[idx_8v20]
arm_8v20  <- as.integer(sputum$Regimen[idx_8v20] == "8M")
tau_8v20  <- safe_tau(t_8v20, arm_8v20)
cat(sprintf("\n8M vs 20M  — tau = %d days\n", tau_8v20))

rmst_8v20 <- rmst2(time = t_8v20, status = s_8v20,
                   arm  = arm_8v20, tau = tau_8v20)
cat("===== RMST: 8M vs 20M =====\n"); print(rmst_8v20)

# ── 9M vs 20M ────────────────────────────────────────────────────────────────
idx_9v20  <- sputum$Regimen %in% c("9M","20M")
t_9v20    <- sputum$time[idx_9v20]
s_9v20    <- sputum$status[idx_9v20]
arm_9v20  <- as.integer(sputum$Regimen[idx_9v20] == "9M")
tau_9v20  <- safe_tau(t_9v20, arm_9v20)
cat(sprintf("\n9M vs 20M  — tau = %d days\n", tau_9v20))

rmst_9v20 <- rmst2(time = t_9v20, status = s_9v20,
                   arm  = arm_9v20, tau = tau_9v20)
cat("===== RMST: 9M vs 20M =====\n"); print(rmst_9v20)

# ── Summary table ─────────────────────────────────────────────────────────────
rmst_summary <- data.frame(
  Comparison = c("6M vs 20M", "8M vs 20M", "9M vs 20M"),
  Tau_days   = c(tau_6v20, tau_8v20, tau_9v20),
  RMST_diff  = c(rmst_6v20$unadjusted.result[1,1],
                 rmst_8v20$unadjusted.result[1,1],
                 rmst_9v20$unadjusted.result[1,1]),
  Lower_95CI = c(rmst_6v20$unadjusted.result[1,2],
                 rmst_8v20$unadjusted.result[1,2],
                 rmst_9v20$unadjusted.result[1,2]),
  Upper_95CI = c(rmst_6v20$unadjusted.result[1,3],
                 rmst_8v20$unadjusted.result[1,3],
                 rmst_9v20$unadjusted.result[1,3]),
  P_value    = c(rmst_6v20$unadjusted.result[1,4],
                 rmst_8v20$unadjusted.result[1,4],
                 rmst_9v20$unadjusted.result[1,4])
)
cat("\n===== RMST Summary Table =====\n"); print(rmst_summary, row.names = FALSE)


# ── 9. COMPETING RISKS: FINE-GRAY MODEL ──────────────────────────────────────
# tidycmprsk::cuminc() uses a formula interface (works in all recent versions).
# cmprsk::crr() still uses the matrix interface but needs integer failcodes.
#
# status2 levels: 0=censored, 1=conversion (event of interest), 2=died (competing)
# Ensure status2 is a factor with correct levels for tidycmprsk

sputum_cr_base <- sputum %>%
  filter(!is.na(time), time > 0) %>%
  mutate(
    # tidycmprsk needs a factor where level 1 = event of interest
    status2_fct = factor(
      case_when(
        outcome == "Died"   ~ 2L,
        status   == 1L      ~ 1L,
        TRUE                ~ 0L
      ),
      levels = 0:2,
      labels = c("censored", "conversion", "died")
    ),
    Regimen = factor(Regimen, c("6M","8M","9M","20M"))
  )

# 9.1 Cumulative incidence function (tidycmprsk) ─────────────────────────────
library(tidycmprsk)
library(ggsurvfit)

cuminc_fit <- tidycmprsk::cuminc(
  Surv(time, status2_fct) ~ Regimen,
  data = sputum_cr_base
)
print(cuminc_fit)

# Plot CIF with ggsurvfit
# Build CIF main plot (NO add_risktable - build manually for font control)
FS_CIF <- 48   # base font size: output 55cm height × ~0.87

cif_main <- cuminc_fit %>%
  ggcuminc(outcome = "conversion", linewidth = 1.8) +
  add_confidence_interval() +
  scale_color_manual(
    values = c("#1565C0","#2E7D32","#E65100","#B71C1C"),
    labels = c("6M (DS-TB)","8M (DS-TB retreatment)",
               "9M (MDR-TB short)","20M (MDR-TB long)")
  ) +
  scale_fill_manual(
    values = c("#1565C0","#2E7D32","#E65100","#B71C1C"),
    labels = c("6M (DS-TB)","8M (DS-TB retreatment)",
               "9M (MDR-TB short)","20M (MDR-TB long)")
  ) +
  scale_x_continuous(breaks = seq(0, 700, 100)) +
  coord_cartesian(xlim = c(0, 730)) +
  labs(
    title    = "Cumulative Incidence of Culture Conversion",
    subtitle = "Competing event: Death (Fine-Gray subdistribution)",
    x        = NULL,
    y        = "Cumulative Incidence"
  ) +
  theme_bw(base_size = FS_CIF) +
  theme(
    plot.title        = element_text(size = FS_CIF + 4, face = "bold"),
    plot.subtitle     = element_text(size = FS_CIF - 8, color = "grey30"),
    axis.text         = element_text(size = FS_CIF - 4, color = "black"),
    axis.title        = element_text(size = FS_CIF,     face = "bold"),
    legend.text       = element_text(size = FS_CIF - 8),
    legend.title      = element_blank(),
    legend.position   = "inside",
    legend.position.inside = c(0.75, 0.30),
    legend.background = element_rect(fill = "white", color = "grey70",
                                     linewidth = 0.6),
    legend.key.width  = unit(3.0, "cm"),
    legend.key.size   = unit(1.5, "cm"),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(10, 10, 2, 10)
  )

# Build risk table manually
rt_breaks <- seq(0, 700, 100)
pal4_cif  <- c("6M"="#1565C0","8M"="#2E7D32","9M"="#E65100","20M"="#B71C1C")

rt_data <- map_dfr(c("6M","8M","9M","20M"), function(reg) {
  km_sub <- survfit(Surv(time, status) ~ 1,
                    data = sputum_cr_base %>% filter(Regimen == reg))
  t_sub  <- broom::tidy(km_sub)
  map_dfr(rt_breaks, function(b) {
    row <- t_sub %>% filter(time <= b) %>% tail(1)
    nr  <- if (nrow(row) == 0) sum(sputum_cr_base$Regimen == reg,
                                    na.rm = TRUE) else row$n.risk[1]
    tibble(Regimen = reg, time = b, n_risk = nr)
  })
}) %>% mutate(Regimen = factor(Regimen, c("20M","9M","8M","6M")))

cif_risk <- ggplot(rt_data,
                   aes(x = time, y = Regimen,
                       label = n_risk, color = Regimen)) +
  geom_text(size = FS_CIF / 3.2, fontface = "bold") +
  scale_color_manual(
    values = c(                     # order matches factor levels: 20M,9M,8M,6M
      "20M" = "#B71C1C",
      "9M"  = "#E65100",
      "8M"  = "#2E7D32",
      "6M"  = "#1565C0"
    )
  ) +
  scale_x_continuous(breaks = rt_breaks,
                     limits = c(0, 730),
                     expand = expansion(add = c(730*0.04, 730*0.01))) +
  labs(x = "Days from treatment start",
       y = NULL, title = "Number at risk") +
  theme_bw(base_size = FS_CIF) +
  theme(
    legend.position  = "none",
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey70"),
    axis.text.y      = element_text(
      # factor levels = 20M,9M,8M,6M → displayed bottom→top on y-axis
      # ggplot renders lowest factor level at bottom
      color = c("#B71C1C","#E65100","#2E7D32","#1565C0"),
      face  = "bold", size = FS_CIF - 4),
    axis.text.x      = element_text(size = FS_CIF - 6, color = "black"),
    axis.title.x     = element_text(size = FS_CIF,     face = "bold"),
    plot.title       = element_text(size = FS_CIF - 8, color = "grey35",
                                    face = "bold"),
    plot.margin      = margin(2, 10, 8, 10)
  )

cif_plot <- cif_main / cif_risk + plot_layout(heights = c(4, 1))

ggsave("figures/CIF_Conversion.png", cif_plot,
       width = 55, height = 44, dpi = 300, units = "cm", bg = "white")
cat("CIF plot saved.\n")

# 9.2 Unadjusted Fine-Gray (cmprsk::crr) ─────────────────────────────────────
# crr() needs: ftime (numeric), fstatus (integer 0/1/2), cov1 (matrix)
sputum_cr <- sputum_cr_base %>%
  mutate(
    fstatus_int = as.integer(status2_fct) - 1L,  # 0=censor,1=conv,2=died
    R6 = as.integer(Regimen == "6M"),
    R8 = as.integer(Regimen == "8M"),
    R9 = as.integer(Regimen == "9M")
  )

cov_cr <- as.matrix(sputum_cr[, c("R6","R8","R9")])

fg_unadj <- cmprsk::crr(
  ftime    = sputum_cr$time,
  fstatus  = sputum_cr$fstatus_int,
  cov1     = cov_cr,
  failcode = 1
)
cat("\n===== Unadjusted Fine-Gray Model: Culture Conversion =====\n")
summary(fg_unadj)

# 9.3 Adjusted Fine-Gray ──────────────────────────────────────────────────────
sputum_cr2 <- sputum_cr %>%
  mutate(
    age_c      = as.numeric(scale(age)),
    hiv_pos    = as.integer(PatientHIV == "Positive"),
    anemia_yes = as.integer(anemia == "Anemia"),
    neu_lym_c  = as.numeric(scale(NEU_LYM)),
    timika_c   = as.numeric(scale(Timika.score)),
    dur_c      = as.numeric(scale(DUR))
  ) %>%
  filter(complete.cases(age_c, hiv_pos, anemia_yes, neu_lym_c, timika_c, dur_c))

cov_adj <- as.matrix(
  sputum_cr2[, c("R6","R8","R9","age_c","hiv_pos",
                 "anemia_yes","neu_lym_c","timika_c","dur_c")]
)

fg_adj <- cmprsk::crr(
  ftime    = sputum_cr2$time,
  fstatus  = sputum_cr2$fstatus_int,
  cov1     = cov_adj,
  failcode = 1
)
cat("\n===== Adjusted Fine-Gray Model: Culture Conversion =====\n")
summary(fg_adj)



# ── 10. COX PROPORTIONAL HAZARDS MODEL ───────────────────────────────────────

# 10.1 Univariable Cox (regimen only)
cox_uni <- coxph(Surv(time, status) ~ Regimen, data = sputum)
summary(cox_uni)

# 10.2 Check PH assumption
cox_ph <- cox.zph(cox_uni)
print(cox_ph)
plot(cox_ph)

# 10.3 Multivariable Cox
cox_multi <- coxph(
  Surv(time, status) ~ Regimen + age + PatientHIV + anemia +
    NEU_LYM + Timika.score + DUR + Cavity + clini_score,
  data = sputum
)
cat("\n===== Multivariable Cox Model =====\n")
summary(cox_multi)


# Forest plot of Cox HR
cox_tidy <- broom::tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    label = case_when(
      term == "Regimen8M"          ~ "Regimen: 8M vs 6M",
      term == "Regimen9M"          ~ "Regimen: 9M vs 6M",
      term == "Regimen20M"         ~ "Regimen: 20M vs 6M",
      term == "age"                ~ "Age (per year)",
      term == "PatientHIVPositive" ~ "HIV: Positive vs Negative",
      term == "anemiaNormal"       ~ "Anemia: Normal vs Anemia",
      term == "NEU_LYM"            ~ "Neutrophil-Lymphocyte Ratio",
      term == "Timika.score"       ~ "Timika Score",
      term == "DUR"                ~ "Duration of illness",
      term == "CavityYes"          ~ "Cavity: Yes vs No",
      term == "clini_score"        ~ "Clinical Score",
      TRUE                         ~ term
    ),
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""     # no "ns" label
    ),
    hr_label = sprintf("%.2f (%.2f\u2013%.2f)%s",
                       estimate, conf.low, conf.high,
                       ifelse(sig == "", "", paste0(" ", sig))),
    y_pos = rev(seq_len(n()))
  )

FS <- 22   # geom_text size: 22 * 2.85pt = ~63pt at 300dpi — clearly readable
BS <- 56   # ggplot2 base_size

forest_plot <- ggplot(cox_tidy,
       aes(y = y_pos, x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_rect(aes(ymin = y_pos - 0.42, ymax = y_pos + 0.42,
                xmin = 0.07, xmax = 2.8,
                fill = y_pos %% 2 == 0),
            alpha = 0.35, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE"="grey91","FALSE"="white")) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey35", linewidth = 1.2) +
  geom_errorbarh(height = 0.20, linewidth = 1.5, color = "grey10") +
  geom_point(size = 7, shape = 15, color = "black") +
  geom_text(aes(x = 0.063, label = label),
            hjust = 1, size = FS, color = "black") +
  geom_text(aes(x = 3.1, label = hr_label),
            hjust = 0, size = FS, color = "black") +
  annotate("text", x = 3.1,  y = max(cox_tidy$y_pos) + 0.7,
           label = "HR (95% CI)   p-value",
           hjust = 0, size = FS + 1, fontface = "bold") +
  annotate("text", x = 0.063, y = max(cox_tidy$y_pos) + 0.7,
           label = "Variable", hjust = 1,
           size = FS + 1, fontface = "bold") +
  scale_x_continuous(
    trans  = "log",
    breaks = c(0.25, 0.5, 1, 2),
    labels = c("0.25","0.50","1.00","2.00"),
    limits = c(0.02, 12),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0.4, max(cox_tidy$y_pos) + 1.1),
    expand = expansion(add = c(0, 0))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title    = "Hazard Ratios: Time to Culture Conversion",
    subtitle = "Multivariable Cox model  |  Reference: 6M regimen",
    x        = "Hazard Ratio (log scale)",
    y        = NULL
  ) +
  theme_bw(base_size = BS) +
  theme(
    plot.title         = element_text(size = BS + 4, face = "bold"),
    plot.subtitle      = element_text(size = BS - 8, color = "grey30"),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.text.x        = element_text(size = BS - 6, color = "black"),
    axis.title.x       = element_text(size = BS - 2, face = "bold",
                                      margin = margin(t = 10)),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),          # no border — prevents clipping
    axis.line.x        = element_line(color = "grey50", linewidth = 0.6),
    plot.margin        = margin(8, 260, 8, 420)    # left=420 for long labels
  )

ggsave("figures/Cox_ForestPlot.png", forest_plot,
       width  = 60, height = 38, dpi = 300, units = "cm", bg = "white")
cat("Forest plot saved \u2192 figures/Cox_ForestPlot.png\n")



# ── 11. PREDICTION MODEL DATASET ─────────────────────────────────────────────
pred_vars <- c("age","bmi","DUR","TBTreatedBefore","Diabetes","PatientHIV",
               "Regimen","anemia","clini_score","Timika.score","CT_mean",
               "NEU_LYM","ALB","outcome_bin")

data_pred <- sputum %>%
  select(all_of(pred_vars)) %>%
  filter(!is.na(outcome_bin)) %>%
  mutate(
    outcome_bin = droplevels(outcome_bin),
    DM          = factor(ifelse(Diabetes == "Yes", "Diabetes", "No diabetes"))
  )

cat("\nPrediction dataset:\n"); print(dim(data_pred))
cat("Outcome:\n"); print(table(data_pred$outcome_bin, useNA = "ifany"))


# ── 12. MULTIPLE IMPUTATION ───────────────────────────────────────────────────
# Strategy: keep outcome_bin IN the data so it is available after complete().
# Use predictorMatrix to prevent outcome_bin from being imputed (set its row=0)
# but allow it to be used as a predictor for other variables (column stays 1).
n_imp <- 50

# Build predictor matrix: default, then zero out the row for outcome_bin
ini <- mice(data_pred, m = 1, maxit = 0, printFlag = FALSE)
pred_mat <- ini$predictorMatrix
pred_mat["outcome_bin", ] <- 0   # do not impute outcome_bin
meth <- ini$method
meth["outcome_bin"] <- ""        # no imputation method for outcome

imp <- mice(
  data          = data_pred,
  m             = n_imp,
  maxit         = 10,
  method        = meth,
  predictorMatrix = pred_mat,
  printFlag     = FALSE,
  seed          = 1234
)

# Stack all imputed datasets (outcome_bin is already present in each)
data_imp_long <- complete(imp, "long", include = FALSE)
cat("Imputed stacked data dimensions:", dim(data_imp_long), "\n")
cat("outcome_bin present:", "outcome_bin" %in% names(data_imp_long), "\n")
cat("outcome_bin table:\n"); print(table(data_imp_long$outcome_bin, useNA = "ifany"))


# ── 13. FULL LOGISTIC REGRESSION (pooled via mice) ────────────────────────────
fit_full_mice <- with(
  imp,
  glm(outcome_bin ~ age + bmi + DUR + TBTreatedBefore + DM +
        PatientHIV + Regimen + anemia + clini_score +
        Timika.score + CT_mean + NEU_LYM + ALB,
      family = binomial)
)
pooled_full <- pool(fit_full_mice)
cat("\n===== Pooled Full Logistic Regression =====\n")
summary(pooled_full, conf.int = TRUE, exponentiate = TRUE)


# ── 14. LASSO VARIABLE SELECTION (MI + Bootstrap) ────────────────────────────
# Run sequentially — future/multisession is unreliable across platforms.
# Set n_imp_lasso=5, n_boot=20 for a quick test run; scale to 50×200 for
# the final analysis (expect ~1-2 hours on a modern laptop).
n_imp_lasso <- 50    # increase to 50 for final run
n_boot      <- 200   # increase to 200 for final run

set.seed(1234)

# Helper: one bootstrap iteration for LASSO
one_boot_lasso <- function(data_i) {
  idx        <- sample(nrow(data_i), replace = TRUE)
  data_boot  <- data_i[idx, ]
  train_idx  <- createDataPartition(data_boot$outcome_bin, p = 0.8, list = FALSE)
  train_data <- na.omit(data_boot[ train_idx, ])
  test_data  <- na.omit(data_boot[-train_idx, ])

  # Skip if not enough events in either split
  if (sum(train_data$outcome_bin == "Bad outcome") < 5) return(NULL)
  if (sum(test_data$outcome_bin  == "Bad outcome") < 3) return(NULL)

  x <- tryCatch(
    model.matrix(outcome_bin ~ age + bmi + DUR + TBTreatedBefore + DM +
                   PatientHIV + Regimen + anemia + clini_score +
                   Timika.score + CT_mean + NEU_LYM + ALB,
                 data = train_data)[, -1],
    error = function(e) NULL
  )
  if (is.null(x)) return(NULL)

  y <- ifelse(train_data$outcome_bin == "Bad outcome", 1, 0)

  cv_lasso  <- tryCatch(
    cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = 5),
    error = function(e) NULL
  )
  if (is.null(cv_lasso)) return(NULL)

  lasso_mod <- glmnet(x, y, alpha = 1, family = "binomial",
                       lambda = cv_lasso$lambda.min)

  x_test <- tryCatch(
    model.matrix(outcome_bin ~ age + bmi + DUR + TBTreatedBefore + DM +
                   PatientHIV + Regimen + anemia + clini_score +
                   Timika.score + CT_mean + NEU_LYM + ALB,
                 data = test_data)[, -1],
    error = function(e) NULL
  )
  if (is.null(x_test)) return(NULL)

  probs <- as.vector(predict(lasso_mod, newx = x_test, type = "response"))
  obs   <- test_data$outcome_bin

  roc1  <- tryCatch(pROC::roc(obs, probs, quiet = TRUE), error = function(e) NULL)
  brier <- mean((ifelse(obs == "Bad outcome", 1, 0) - probs)^2)

  sel_vars <- rownames(coef(lasso_mod))[as.vector(coef(lasso_mod)) != 0]
  sel_vars <- sel_vars[sel_vars != "(Intercept)"]

  list(sel_vars = sel_vars, roc1 = roc1, brier = brier)
}

results_lasso <- vector("list", n_imp_lasso)
for (i in seq_len(n_imp_lasso)) {
  data_i <- complete(imp, action = i, include = FALSE)
  results_lasso[[i]] <- lapply(seq_len(n_boot), function(j) one_boot_lasso(data_i))
  cat(sprintf("LASSO imputation %d/%d done\n", i, n_imp_lasso))
}

# Variable inclusion frequency (exclude NULLs)
all_boot_lasso <- unlist(results_lasso, recursive = FALSE)
all_boot_lasso <- Filter(Negate(is.null), all_boot_lasso)

var_counts <- table(unlist(lapply(all_boot_lasso, `[[`, "sel_vars")))
var_freq   <- var_counts / length(all_boot_lasso)
cat("\nLASSO variable inclusion frequencies:\n")
print(sort(var_freq, decreasing = TRUE))

selected_vars_lasso <- names(var_freq[var_freq > 0.50])
cat("\nSelected dummy names (>50%):", selected_vars_lasso, "\n")

all_rocs   <- sapply(all_boot_lasso, function(x)
                if (!is.null(x$roc1)) as.numeric(x$roc1$auc) else NA)
all_briers <- sapply(all_boot_lasso, `[[`, "brier")

cat(sprintf("Mean AUC (LASSO):   %.3f\n", mean(all_rocs,   na.rm = TRUE)))
cat(sprintf("Mean Brier (LASSO): %.3f\n", mean(all_briers, na.rm = TRUE)))

save(results_lasso, selected_vars_lasso,
     file = "output/LASSO_results.RData")


# ── 15. FINAL LASSO MODEL (display, validate, calibrate) ──────────────────────
# Map dummy variable names (e.g. "PatientHIVPositive", "Regimen9M") back to
# the original column names that lrm() / rms understands.
# Strategy: build a reference matrix from the full data to get the
# mapping dummy-name → original-variable-name.

dummy_to_original <- function(dummy_names, formula_vars) {
  # Handles three cases:
  #  1. rcs(var, k)var['] → extract var name from inside rcs()
  #  2. varLEVEL           → strip the level suffix using formula_vars as prefix list
  #  3. plain var name     → return as-is
  orig <- character(0)
  for (d in dummy_names) {
    # Case 1: rcs() spline terms — name looks like "rcs(age, 3)age" or "rcs(age, 3)age'"
    if (grepl("^rcs\\(", d)) {
      # Extract variable name from inside rcs(...)
      v <- sub("^rcs\\(([^,]+),.*", "\\1", d)
      v <- trimws(v)
      orig <- union(orig, v)
      next
    }
    # Case 2: factor dummy — name is varLEVEL (e.g. "PatientHIVPositive")
    # Find the longest formula_var that is a prefix of d
    matches <- formula_vars[sapply(formula_vars, function(v) startsWith(d, v))]
    if (length(matches) > 0) {
      orig <- union(orig, matches[which.max(nchar(matches))])
      next
    }
    # Case 3: plain continuous variable — keep as-is
    orig <- union(orig, d)
  }
  orig
}

all_formula_vars <- c("age","bmi","DUR","TBTreatedBefore","DM",
                       "PatientHIV","Regimen","anemia","clini_score",
                       "Timika.score","CT_mean","NEU_LYM","ALB")

# Resolve LASSO selected dummies → original variable names
if (length(selected_vars_lasso) == 0) {
  orig_vars_lasso <- c("age", "clini_score", "NEU_LYM")
  warning("No LASSO variables selected at >50%; using fallback set.")
} else {
  orig_vars_lasso <- dummy_to_original(selected_vars_lasso, all_formula_vars)
}
cat("\nLASSO original variable names:", orig_vars_lasso, "\n")

dd <- datadist(data_imp_long); options(datadist = "dd")

lasso_formula <- as.formula(paste(
  "outcome_bin ~",
  paste(orig_vars_lasso, collapse = " + ")
))

f_lasso <- rms::lrm(lasso_formula,
                     data = na.omit(data_imp_long),
                     x = TRUE, y = TRUE)
cat("\n===== rms LRM (LASSO-selected) =====\n"); print(f_lasso)

# Nomogram
png("figures/Nomogram_LASSO.png", width = 2400, height = 1200, res = 150)
plot(nomogram(f_lasso,
              lp       = FALSE,
              fun      = plogis,
              fun.at   = c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95),
              funlabel = "Probability of Unfavorable Outcome"))
dev.off()

# Bootstrap validation & calibration
val_lasso <- validate(f_lasso, B = 200)
print(val_lasso)
c_lasso_optcorr <- 0.5 * (val_lasso[1, 5] + 1)
cat(sprintf("Optimism-corrected C-statistic (LASSO): %.3f\n", c_lasso_optcorr))

cal_lasso <- rms::calibrate(f_lasso, B = 200)
png("figures/Calibration_LASSO.png", width = 2000, height = 1600, res = 150)
plot(cal_lasso, main = "Calibration: LASSO Model"); dev.off()


# ── 16. AIC STEPWISE MODEL ────────────────────────────────────────────────────
n_imp_aic  <- 50    # increase to 50 for final run
n_boot_aic <- 200   # increase to 200 for final run

set.seed(1234)

# Helper: one bootstrap iteration for AIC stepwise
one_boot_aic <- function(data_i) {
  idx       <- sample(nrow(data_i), replace = TRUE)
  data_boot <- na.omit(data_i[idx, ])

  if (sum(data_boot$outcome_bin == "Bad outcome") < 5) return(NULL)

  full_mod <- tryCatch(
    glm(outcome_bin ~ rcs(age,3) + rcs(bmi,3) + rcs(DUR,3) +
          TBTreatedBefore + DM + PatientHIV + Regimen + anemia +
          rcs(clini_score,3) + rcs(Timika.score,3) +
          rcs(CT_mean,3) + rcs(NEU_LYM,3) + rcs(ALB,3),
        data   = data_boot,
        family = binomial(link = "logit")),
    error = function(e) NULL
  )
  if (is.null(full_mod)) return(NULL)

  step_mod <- tryCatch(
    stepAIC(full_mod, direction = "backward", trace = FALSE),
    error = function(e) NULL
  )
  if (is.null(step_mod)) return(NULL)

  list(sel_vars = names(coef(step_mod))[-1])
}

results_aic <- vector("list", n_imp_aic)
for (i in seq_len(n_imp_aic)) {
  data_i <- complete(imp, action = i, include = FALSE)
  results_aic[[i]] <- lapply(seq_len(n_boot_aic), function(j) one_boot_aic(data_i))
  cat(sprintf("AIC imputation %d/%d done\n", i, n_imp_aic))
}

all_boot_aic <- unlist(results_aic, recursive = FALSE)
all_boot_aic <- Filter(Negate(is.null), all_boot_aic)

var_counts_aic <- table(unlist(lapply(all_boot_aic, `[[`, "sel_vars")))
var_freq_aic   <- var_counts_aic / length(all_boot_aic)
cat("\nAIC variable inclusion frequencies:\n")
print(sort(var_freq_aic, decreasing = TRUE))

selected_vars_aic <- names(var_freq_aic[var_freq_aic > 0.50])
cat("\nAIC selected dummy names (>50%):", selected_vars_aic, "\n")

# Re-source the function to ensure the 2-argument version is active
# (guards against a stale 3-argument version cached from a previous run)
dummy_to_original <- function(dummy_names, formula_vars) {
  orig <- character(0)
  for (d in dummy_names) {
    if (grepl("^rcs\\(", d)) {
      v <- trimws(sub("^rcs\\(([^,]+),.*", "\\1", d))
      orig <- union(orig, v)
      next
    }
    matches <- formula_vars[sapply(formula_vars, function(v) startsWith(d, v))]
    if (length(matches) > 0) {
      orig <- union(orig, matches[which.max(nchar(matches))])
      next
    }
    orig <- union(orig, d)
  }
  orig
}

# Map dummy names → original variable names
if (length(selected_vars_aic) == 0) {
  orig_vars_aic <- c("age", "clini_score", "NEU_LYM")
  warning("No AIC variables selected at >50%; using fallback set.")
} else {
  orig_vars_aic <- dummy_to_original(selected_vars_aic, all_formula_vars)
}
cat("\nAIC original variable names:", orig_vars_aic, "\n")

# Final AIC model (rms)
aic_formula <- as.formula(paste(
  "outcome_bin ~",
  paste(orig_vars_aic, collapse = " + ")
))

f_aic <- rms::lrm(aic_formula,
                   data = na.omit(data_imp_long),
                   x = TRUE, y = TRUE)
cat("\n===== rms LRM (AIC-selected) =====\n"); print(f_aic)

# Nomogram
png("figures/Nomogram_AIC.png", width = 2400, height = 1200, res = 150)
plot(nomogram(f_aic,
              lp       = FALSE,
              fun      = plogis,
              fun.at   = c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95),
              funlabel = "Probability of Unfavorable Outcome"))
dev.off()

# Bootstrap validation
val_aic <- validate(f_aic, B = 200)
print(val_aic)
c_aic_optcorr <- 0.5 * (val_aic[1, 5] + 1)
cat(sprintf("Optimism-corrected C-statistic (AIC): %.3f\n", c_aic_optcorr))

# Calibration
cal_aic <- rms::calibrate(f_aic, B = 200)
png("figures/Calibration_AIC.png",   width = 2000, height = 1600, res = 150)
plot(cal_aic, main = "Calibration: AIC Model")
dev.off()

save(selected_vars_aic, results_aic, f_aic, val_aic, cal_aic,
     file = "output/AIC_results.RData")


# ── 17. eventglm: DIRECT REGRESSION ON CUMULATIVE INCIDENCE ──────────────────
#
# RATIONALE FOR eventglm IN THIS STUDY:
# ─────────────────────────────────────
# The primary question is: "Does treatment regimen affect culture conversion,
# accounting for competing events (death) and baseline covariates?"
#
# eventglm::cumincglm() offers several advantages over Cox/Fine-Gray:
#   1. DIRECTLY models P(T≤t) at a pre-specified clinical horizon (e.g. 300d)
#      → absolute risk difference (RD) and risk ratio (RR) with 95% CI
#      → more clinically interpretable than subdistribution hazard ratio
#   2. Handles RIGHT-CENSORING via inverse probability of censoring weighting (IPCW)
#      with multiple censoring models (KM, Cox, parametric)
#   3. Supports multiple link functions: identity (RD), log (RR), logit (OR),
#      cloglog (proportional hazards check)
#   4. rmeanglm() estimates RMST differences in a regression framework,
#      with covariate adjustment and robust SE
#   5. Time-varying effects via tve() to detect non-constant treatment effects
#
# KEY DIFFERENCE from Fine-Gray: eventglm gives a MARGINAL (population-averaged)
# estimate at time t; Fine-Gray gives a CONDITIONAL subdistribution HR.
# Both are valid; eventglm results are easier to communicate to clinicians.
#
# APPROPRIATE tau per regimen-specific max follow-up:
#   Culture conversion endpoint: tau = 180d for 6M, 240d for 8M,
#   365d for 9M, 730d for 20M → use 180d as common tau for cross-regimen comparison
# --------------------------------------------------------------------------

library(eventglm)
options(na.action = "na.exclude")

# Build sputum_eg directly from sputum (not sputum_cr_base which is filtered)
# eventglm requires: integer 0/1/2 event code, no factor; na.action=na.exclude
options(na.action = "na.exclude")

sputum_eg <- sputum %>%
  filter(!is.na(time), time > 0,
         !is.na(Regimen), !is.na(status)) %>%
  mutate(
    Regimen   = factor(Regimen, c("6M","8M","9M","20M")),
    # 0 = censored, 1 = conversion (event of interest), 2 = died (competing)
    event_int = case_when(
      !is.na(outcome) & outcome == "Died" ~ 2L,
      status == 1L                        ~ 1L,
      TRUE                                ~ 0L
    )
  )

cat(sprintf("sputum_eg: n=%d | event=1: %d | competing=2: %d | censored=0: %d\n",
            nrow(sputum_eg), sum(sputum_eg$event_int==1),
            sum(sputum_eg$event_int==2), sum(sputum_eg$event_int==0)))

# tau = floor of min(max observed time per regimen), capped at 240
tau_common <- floor(min(tapply(sputum_eg$time, sputum_eg$Regimen,
                               function(x) max(x[x < 700], na.rm=TRUE))))
tau_common <- min(tau_common, 240L)
cat(sprintf("Common tau for eventglm: %d days\n", tau_common))

# ── eventglm NOTES ────────────────────────────────────────────────────────────
# cumincglm() requires:
#   - Binary status (0/1) for simple survival  → Surv(time, status_bin)
#   - For competing risks: Surv(time, event_int, type="mstate") where
#     event_int is a factor with levels 0,1,2
#   - na.action = "na.exclude" (already set above)
#   - Complete cases on all covariates in the formula
#
# We use binary status (conversion=1 vs not-converted=0) for all eventglm models.
# Competing risks are handled separately via tidycmprsk (Section 9).
# ─────────────────────────────────────────────────────────────────────────────

# Prepare clean dataset: complete cases on covariates used in models
sputum_eg_cc <- sputum_eg %>%
  select(studycode, Regimen, time, status, event_int,
         age, PatientHIV, anemia, NEU_LYM, Timika.score, DUR) %>%
  filter(complete.cases(.)) %>%
  mutate(
    # Binary status for eventglm (simpler, avoids mstate issues)
    status_bin = as.integer(status == 1)
  )

cat(sprintf("sputum_eg_cc: n=%d, events=%d\n",
            nrow(sputum_eg_cc), sum(sputum_eg_cc$status_bin)))

# tau: min of max event time across regimens (exclude censored=720)
tau_eg <- sputum_eg_cc %>%
  filter(status_bin == 1) %>%
  group_by(Regimen) %>%
  summarise(max_t = max(time), .groups="drop") %>%
  pull(max_t) %>% min()
tau_eg <- min(floor(tau_eg), 150L)   # conservative: 150d (all groups have events)
cat(sprintf("eventglm tau: %d days\n", tau_eg))

# ── 17a. Unadjusted CID (identity link) ──────────────────────────────────────
ci_unadj <- tryCatch(
  eventglm::cumincglm(
    Surv(time, status_bin) ~ Regimen,
    time = tau_eg,
    data = sputum_eg_cc
  ),
  error = function(e) { cat("ci_unadj error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(ci_unadj)) {
  cat("\n=== Unadjusted CID (identity link, ref=20M) ===\n")
  print(summary(ci_unadj))
  cat("\nRobust 95% CI:\n"); print(confint(ci_unadj, type = "robust"))
}

# ── 17b. RR (log link) ───────────────────────────────────────────────────────
ci_rr <- tryCatch(
  eventglm::cumincglm(
    Surv(time, status_bin) ~ Regimen,
    time = tau_eg, data = sputum_eg_cc, link = "log"
  ),
  error = function(e) { cat("ci_rr error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(ci_rr)) {
  cat("\n=== Unadjusted RR (log link) ===\n")
  cat("exp(coef):\n"); print(exp(coef(ci_rr)))
  cat("exp(confint):\n"); print(exp(confint(ci_rr)))
}

# ── 17c. Adjusted CID (IPCW-Cox) ─────────────────────────────────────────────
ci_adj <- tryCatch(
  eventglm::cumincglm(
    Surv(time, status_bin) ~ Regimen + age + PatientHIV + anemia +
      NEU_LYM + Timika.score + DUR,
    time              = tau_eg,
    data              = sputum_eg_cc,
    model.censoring   = "coxph",
    formula.censoring = ~ Regimen + age + PatientHIV
  ),
  error = function(e) { cat("ci_adj error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(ci_adj)) {
  cat("\n=== Adjusted CID (IPCW-Cox) ===\n")
  print(summary(ci_adj))
  cat("\nRobust CI:\n"); print(confint(ci_adj, type = "robust"))
}

# ── 17d. RMST regression ─────────────────────────────────────────────────────
rm_unadj <- tryCatch(
  eventglm::rmeanglm(
    Surv(time, status_bin) ~ Regimen,
    time = tau_eg, data = sputum_eg_cc
  ),
  error = function(e) { cat("rm_unadj error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(rm_unadj)) {
  cat("\n=== Unadjusted RMST differences (ref=20M) ===\n")
  print(summary(rm_unadj))
  cat("\nRobust CI:\n"); print(confint(rm_unadj, type = "robust"))
}

rm_adj <- tryCatch(
  eventglm::rmeanglm(
    Surv(time, status_bin) ~ Regimen + age + PatientHIV + anemia +
      NEU_LYM + Timika.score + DUR,
    time              = tau_eg,
    data              = sputum_eg_cc,
    model.censoring   = "coxph",
    formula.censoring = ~ Regimen + age + PatientHIV
  ),
  error = function(e) { cat("rm_adj error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(rm_adj)) {
  cat("\n=== Adjusted RMST regression ===\n"); print(summary(rm_adj))
}

# ── 17e. Time-varying effects ─────────────────────────────────────────────────
taus_tve <- c(30, 60, 90, 120)[c(30,60,90,120) <= tau_eg]
if (length(taus_tve) >= 2) {
  ci_tve <- tryCatch(
    eventglm::cumincglm(
      Surv(time, status_bin) ~ tve(Regimen),
      time = taus_tve, data = sputum_eg_cc
    ),
    error = function(e) { cat("ci_tve error:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(ci_tve)) {
    cat("\n=== Time-varying Regimen effect ===\n"); print(summary(ci_tve))
  }
}

# ── 17f. eventglm full prediction model (MI + Bootstrap) ──────────────────────
# Uses eventglm::cumincglm() as the prediction engine inside the
# MI × bootstrap loop — directly models P(conversion ≤ tau) as outcome
# This is particularly suited when the outcome is time-to-event with competing
# risks, rather than a simple binary classification.

cat("\n=== eventglm Full Prediction Model (MI + Bootstrap) ===\n")

n_imp_eg  <- 50    # increase to 50 for final run
n_boot_eg <- 200   # increase to 200 for final run
set.seed(1234)

# Join sputum_eg_cc (has status_bin) back with sputum to get all pred variables
data_pred_eg <- sputum_eg_cc %>%
  select(studycode, status_bin, time) %>%
  left_join(
    sputum %>% select(studycode, age, bmi, DUR, TBTreatedBefore, DM,
                      PatientHIV, Regimen, anemia, clini_score,
                      Timika.score, CT_mean, NEU_LYM, ALB),
    by = "studycode"
  ) %>%
  filter(!is.na(status_bin), !is.na(time))

ini_eg  <- mice(data_pred_eg, m = 1, maxit = 0, printFlag = FALSE)
pred_eg <- ini_eg$predictorMatrix
meth_eg <- ini_eg$method
# Outcomes — do not impute
pred_eg[c("status_bin","time"), ] <- 0
meth_eg[c("status_bin","time")]   <- ""

imp_eg <- mice(data_pred_eg, m = n_imp_eg, maxit = 10,
               predictorMatrix = pred_eg, method = meth_eg,
               printFlag = FALSE, seed = 1234)

one_boot_eg <- function(data_i, tau_val) {
  idx       <- sample(nrow(data_i), replace = TRUE)
  data_boot <- na.omit(data_i[idx, ])
  if (nrow(data_boot) < 50) return(NULL)
  if (sum(data_boot$status_bin == 1) < 10) return(NULL)

  train_idx  <- createDataPartition(data_boot$status_bin, p = 0.8, list = FALSE)
  train_data <- data_boot[ train_idx, ]
  test_data  <- data_boot[-train_idx, ]
  if (sum(train_data$status_bin == 1) < 5) return(NULL)
  if (sum(test_data$status_bin  == 1) < 3) return(NULL)

  mod <- tryCatch(
    eventglm::cumincglm(
      Surv(time, status_bin) ~ Regimen + age + DUR + PatientHIV +
        anemia + clini_score + Timika.score + NEU_LYM + ALB,
      time             = tau_val,
      data             = train_data,
      model.censoring  = "coxph",
      formula.censoring = ~ Regimen + age
    ),
    error = function(e) NULL
  )
  if (is.null(mod)) return(NULL)

  probs <- tryCatch(
    as.vector(predict(mod, newdata = test_data, type = "response")),
    error = function(e) NULL
  )
  if (is.null(probs) || length(probs) == 0) return(NULL)
  probs <- pmax(pmin(probs, 0.9999), 0.0001)  # clip for ROC

  obs_bin <- as.integer(test_data$status_bin == 1)
  if (length(unique(obs_bin)) < 2) return(NULL)

  roc1  <- tryCatch(
    as.numeric(pROC::roc(obs_bin, probs, quiet = TRUE)$auc),
    error = function(e) NA_real_
  )
  brier <- mean((obs_bin - probs)^2, na.rm = TRUE)

  list(auc = roc1, brier = brier)
}

results_eg <- vector("list", n_imp_eg)
for (i in seq_len(n_imp_eg)) {
  data_i <- complete(imp_eg, action = i, include = FALSE)
  results_eg[[i]] <- lapply(seq_len(n_boot_eg),
                             function(j) one_boot_eg(data_i, tau_eg))
  cat(sprintf("eventglm MI %d/%d done\n", i, n_imp_eg))
}

all_eg    <- Filter(Negate(is.null), unlist(results_eg, recursive = FALSE))
rocs_eg   <- sapply(all_eg, `[[`, "auc")
briers_eg <- sapply(all_eg, `[[`, "brier")
cat(sprintf("\neventglm model — Mean AUC: %.3f | Mean Brier: %.3f | n_boots: %d\n",
            mean(rocs_eg,   na.rm = TRUE),
            mean(briers_eg, na.rm = TRUE),
            length(all_eg)))


# ── 18. MODEL IMPROVEMENT STRATEGIES & PERFORMANCE SUMMARY ───────────────────
#
# INTERPRETATION OF AUC ~0.70–0.78 AND BRIER ~0.16–0.20:
# ─────────────────────────────────────────────────────────
# • AUC 0.70–0.78 = "acceptable to good" discrimination for a clinical
#   prediction model in TB outcome research (comparable to published models).
# • Brier score: for binary outcome with ~20% bad-outcome prevalence,
#   null model Brier = p*(1-p) ≈ 0.16. A model Brier near 0.16 suggests
#   little improvement over the null → indicates poor calibration or
#   that predictors are weak on THIS dataset.
# • CONCLUSION: current predictors explain modest variance in outcome.
#   This is expected given: (a) heterogeneous regimens; (b) missing WGS/DST;
#   (c) outcome partially determined by adherence (unmeasured).
#
# IMPROVEMENT STRATEGIES IMPLEMENTED BELOW:
# ─────────────────────────────────────────
# 1. Full model (all candidate predictors, no selection) — often best calibrated
# 2. Add interaction terms: Regimen × HIV, Regimen × Timika.score
# 3. Natural splines for continuous predictors (age, NEU_LYM, Timika.score)
# 4. Two-stage: use Regimen as stratum, fit within-stratum models
# 5. Outcome redefinition: time-to-event probability at regimen-specific tau
#    via eventglm (above) vs. crude binary — eventglm better suits the data

# ── 18a. Full model (all candidates) ─────────────────────────────────────────
f_full <- rms::lrm(
  outcome_bin ~ age + bmi + DUR + TBTreatedBefore + DM +
    PatientHIV + Regimen + anemia + clini_score +
    Timika.score + CT_mean + NEU_LYM + ALB,
  data = na.omit(data_imp_long),
  x = TRUE, y = TRUE
)
cat("\n=== Full Model (all candidates) ===\n"); print(f_full)

val_full <- validate(f_full, B = 200)
c_full   <- 0.5 * (val_full[1, 5] + 1)
cat(sprintf("Full model optimism-corrected C: %.3f\n", c_full))

cal_full <- rms::calibrate(f_full, B = 200)
png("figures/Calibration_Full.png",  width = 2000, height = 1600, res = 150)
plot(cal_full, main = "Calibration: Full Model (all candidates)"); dev.off()

# ── 18b. Enriched model: splines + interactions ───────────────────────────────
# Helper: add scaled + interaction columns (reused in data_dca and ROC section)
enrich_cols <- function(df) {
  df %>% mutate(
    age_sc    = as.numeric(scale(age)),
    timika_sc = as.numeric(scale(Timika.score)),
    hiv_num   = as.integer(PatientHIV == "Positive"),
    reg_8M    = as.integer(Regimen == "8M"),
    reg_9M    = as.integer(Regimen == "9M"),
    reg_20M   = as.integer(Regimen == "20M"),
    hiv_x_8M  = hiv_num * reg_8M,  hiv_x_9M  = hiv_num * reg_9M,
    hiv_x_20M = hiv_num * reg_20M,
    tim_x_8M  = timika_sc * reg_8M, tim_x_9M  = timika_sc * reg_9M,
    tim_x_20M = timika_sc * reg_20M
  )
}

# Pre-compute interaction columns explicitly to avoid rms formula parsing issues
# with mixed rcs() + interaction terms.
data_enrich <- na.omit(data_imp_long) %>%
  mutate(
    age_sc    = as.numeric(scale(age)),
    timika_sc = as.numeric(scale(Timika.score)),
    hiv_num   = as.integer(PatientHIV == "Positive"),
    # Pre-compute Regimen dummies for interactions (ref = 6M)
    reg_8M    = as.integer(Regimen == "8M"),
    reg_9M    = as.integer(Regimen == "9M"),
    reg_20M   = as.integer(Regimen == "20M"),
    # Explicit interaction columns
    hiv_x_8M  = hiv_num * reg_8M,
    hiv_x_9M  = hiv_num * reg_9M,
    hiv_x_20M = hiv_num * reg_20M,
    tim_x_8M  = timika_sc * reg_8M,
    tim_x_9M  = timika_sc * reg_9M,
    tim_x_20M = timika_sc * reg_20M
  )

dd_enrich <- datadist(data_enrich); options(datadist = "dd_enrich")

f_enrich <- rms::lrm(
  outcome_bin ~ rcs(age_sc,3) + rcs(NEU_LYM,3) + rcs(DUR,3) +
    timika_sc + PatientHIV + Regimen + anemia + CT_mean +
    hiv_x_8M + hiv_x_9M + hiv_x_20M +
    tim_x_8M + tim_x_9M + tim_x_20M,
  data = data_enrich,
  x = TRUE, y = TRUE
)
cat("\n=== Enriched Model (splines + interactions) ===\n"); print(f_enrich)
cat("\n=== Enriched Model (splines + interactions) ===\n"); print(f_enrich)

val_enrich <- validate(f_enrich, B = 200)
c_enrich   <- 0.5 * (val_enrich[1, 5] + 1)
cat(sprintf("Enriched model optimism-corrected C: %.3f\n", c_enrich))

# ── 18c. Performance comparison table ────────────────────────────────────────
perf_table <- data.frame(
  Model = c("LASSO-selected", "AIC-selected", "Full (all vars)", "Enriched (splines+interaction)"),
  C_stat_raw       = c(f_lasso$stats["C"], f_aic$stats["C"],
                       f_full$stats["C"],  f_enrich$stats["C"]),
  C_stat_optimism  = c(c_lasso_optcorr, c_aic_optcorr, c_full, c_enrich),
  Dxy              = c(f_lasso$stats["Dxy"], f_aic$stats["Dxy"],
                       f_full$stats["Dxy"], f_enrich$stats["Dxy"]),
  R2               = c(f_lasso$stats["R2"], f_aic$stats["R2"],
                       f_full$stats["R2"],  f_enrich$stats["R2"])
)
cat("\n=== Model Performance Comparison ===\n")
print(perf_table, digits = 3, row.names = FALSE)

# ── 18d. Decision Curve Analysis (net benefit) ───────────────────────────────
# DCA assesses clinical utility across probability thresholds
data_dca <- na.omit(data_imp_long) %>%
  enrich_cols() %>%
  mutate(
    p_lasso  = predict(f_lasso,  newdata = ., type = "fitted"),
    p_aic    = predict(f_aic,    newdata = ., type = "fitted"),
    p_full   = predict(f_full,   newdata = ., type = "fitted"),
    p_enrich = tryCatch(
      predict(f_enrich, newdata = ., type = "fitted"),
      error = function(e) rep(NA_real_, nrow(.))
    ),
    outcome_num = as.integer(outcome_bin == "Bad outcome")
  )

# Simple DCA using rmda or manual net-benefit calculation
thresholds <- seq(0.05, 0.50, by = 0.01)
nb_lasso  <- sapply(thresholds, function(t) {
  tp <- mean(data_dca$p_lasso  > t & data_dca$outcome_num == 1)
  fp <- mean(data_dca$p_lasso  > t & data_dca$outcome_num == 0)
  tp - fp * t / (1 - t)
})
nb_full   <- sapply(thresholds, function(t) {
  tp <- mean(data_dca$p_full   > t & data_dca$outcome_num == 1)
  fp <- mean(data_dca$p_full   > t & data_dca$outcome_num == 0)
  tp - fp * t / (1 - t)
})
nb_enrich <- sapply(thresholds, function(t) {
  tp <- mean(data_dca$p_enrich > t & data_dca$outcome_num == 1)
  fp <- mean(data_dca$p_enrich > t & data_dca$outcome_num == 0)
  tp - fp * t / (1 - t)
})
nb_all    <- mean(data_dca$outcome_num) - thresholds / (1 - thresholds) * (1 - mean(data_dca$outcome_num))

dca_df <- data.frame(
  threshold = thresholds,
  LASSO     = nb_lasso,
  Full      = nb_full,
  Enriched  = nb_enrich,
  Treat_all = nb_all,
  Treat_none = 0
)

dca_plot <- dca_df %>%
  pivot_longer(-threshold, names_to = "Model", values_to = "NetBenefit") %>%
  ggplot(aes(x = threshold, y = NetBenefit, color = Model, linetype = Model)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c(LASSO="#2196F3", Full="#4CAF50",
                                 Enriched="#FF9800", Treat_all="grey50",
                                 Treat_none="grey20")) +
  scale_linetype_manual(values = c(LASSO="solid", Full="solid",
                                    Enriched="solid", Treat_all="dashed",
                                    Treat_none="dotted")) +
  coord_cartesian(ylim = c(-0.05, max(nb_full, na.rm = TRUE) * 1.2)) +
  labs(title    = "Decision Curve Analysis",
       subtitle = "Net benefit of prediction models vs treat-all / treat-none",
       x = "Probability threshold",
       y = "Net benefit") +
  theme_bw(base_size = 13)

ggsave("figures/DCA_Models.png", dca_plot,
       width = 22, height = 16, dpi = 300, units = "cm", bg = "white")
cat("DCA plot saved.\n")

# ── 18e. ROC comparison: all four models ─────────────────────────────────────
set.seed(42)
data_cc    <- na.omit(data_imp_long)
tr_idx     <- createDataPartition(data_cc$outcome_bin, p = 0.8, list = FALSE)
train_roc  <- data_cc[ tr_idx, ]
test_roc   <- data_cc[-tr_idx, ]

fit_tr <- function(formula, data_tr = train_roc) {
  glm(formula, data = data_tr, family = binomial, na.action = na.omit)
}
pred_pr <- function(mod, data_te = test_roc) predict(mod, data_te, type = "response")

# Add interaction columns to train/test for enriched model
train_roc_e <- enrich_cols(train_roc)
test_roc_e  <- enrich_cols(test_roc)

m_lasso_roc  <- fit_tr(lasso_formula)
m_aic_roc    <- fit_tr(aic_formula)
m_full_roc   <- fit_tr(outcome_bin ~ age + bmi + DUR + TBTreatedBefore + DM +
                         PatientHIV + Regimen + anemia + clini_score +
                         Timika.score + CT_mean + NEU_LYM + ALB)
m_enrich_roc <- fit_tr(
  outcome_bin ~ rcs(age_sc,3) + rcs(NEU_LYM,3) + rcs(DUR,3) +
    timika_sc + PatientHIV + Regimen + anemia + CT_mean +
    hiv_x_8M + hiv_x_9M + hiv_x_20M +
    tim_x_8M + tim_x_9M + tim_x_20M,
  data_tr = train_roc_e
)

obs_r <- test_roc$outcome_bin
roc_l <- pROC::roc(obs_r, pred_pr(m_lasso_roc),  quiet = TRUE)
roc_a <- pROC::roc(obs_r, pred_pr(m_aic_roc),    quiet = TRUE)
roc_f <- pROC::roc(obs_r, pred_pr(m_full_roc),   quiet = TRUE)
roc_e <- pROC::roc(obs_r, pred_pr(m_enrich_roc, data_te = test_roc_e), quiet = TRUE)

png("figures/ROC_Comparison.png",     width = 2000, height = 1800, res = 150)
plot(roc_l, col = "#2196F3", lwd = 2, print.auc = TRUE,
     print.auc.y = 0.48, main = "ROC Curves: Model Comparison")
plot(roc_a, col = "#4CAF50", lwd = 2, add = TRUE,
     print.auc = TRUE, print.auc.y = 0.42)
plot(roc_f, col = "#FF9800", lwd = 2, add = TRUE,
     print.auc = TRUE, print.auc.y = 0.36)
plot(roc_e, col = "#D32F2F", lwd = 2, add = TRUE,
     print.auc = TRUE, print.auc.y = 0.30)
legend("bottomright",
       legend = c(sprintf("LASSO    (AUC=%.3f)", as.numeric(roc_l$auc)),
                  sprintf("AIC      (AUC=%.3f)", as.numeric(roc_a$auc)),
                  sprintf("Full     (AUC=%.3f)", as.numeric(roc_f$auc)),
                  sprintf("Enriched (AUC=%.3f)", as.numeric(roc_e$auc))),
       col = c("#2196F3","#4CAF50","#FF9800","#D32F2F"), lwd = 2)
dev.off()

# DeLong pairwise tests
cat("\n=== DeLong AUC comparison: Full vs LASSO ===\n")
print(roc.test(roc_f, roc_l, reuse.auc = FALSE))
cat("\n=== DeLong AUC comparison: Enriched vs Full ===\n")
print(roc.test(roc_e, roc_f, reuse.auc = FALSE))


# ── 19. FRAMEWORK COMPARISON & OPTIMAL ANALYTICAL STRATEGY ───────────────────
#
# WHY eventglm AUC (0.704) > glm/lrm AUC (0.64):
# ─────────────────────────────────────────────────
# eventglm outcome = P(convert ≤ tau=150d): continuous [0,1], rich signal
#   → Regimen is a very strong predictor because DS-TB (6M) converts fast
#     by definition (~85% by d90) while MDR-TB (9M/20M) converts slowly (~5%)
#   → This is DEFINITIONAL separation, not true clinical prediction
#   → AUC reflects regimen-group differences, not within-group risk stratification
#
# glm/lrm outcome = Bad outcome (Died/Failed) vs Good (Cured/Completed):
#   → Excludes "Lost to follow-up" (63 patients) → unbiased for completers
#   → Regimen is less dominant because within-regimen variability matters
#   → AUC ~0.64 reflects genuine clinical prediction difficulty
#
# CONCLUSION:
# ┌─────────────────────┬───────────────────────────────────────────────────┐
# │ Research question   │ Recommended framework                             │
# ├─────────────────────┼───────────────────────────────────────────────────┤
# │ "Does regimen X     │ eventglm/survRM2 RMST regression                  │
# │  improve conversion │ (absolute risk difference, clinically meaningful) │
# │  rate vs regimen Y?"│                                                   │
# ├─────────────────────┼───────────────────────────────────────────────────┤
# │ "Which patients are │ glm/lrm on outcome_bin                            │
# │  at high risk of    │ (prediction model, Sections 13-18)                │
# │  treatment failure?"│                                                   │
# ├─────────────────────┼───────────────────────────────────────────────────┤
# │ "What predicts time │ Stratified Cox or rms::psm() within each Regimen  │
# │  to conversion      │ (WITHIN-REGIMEN survival model)                   │
# │  within a regimen?" │                                                   │
# └─────────────────────┴───────────────────────────────────────────────────┘
#
# OPTIMAL SOLUTION: WITHIN-REGIMEN PREDICTION MODELS
# ────────────────────────────────────────────────────
# Since Regimen is a design variable (assigned, not random), the most
# clinically useful prediction is: WITHIN each regimen, who will fail?
# This removes the dominant between-regimen confounding.

cat("\n=== Section 19: Within-Regimen Prediction Models ===\n")

# ── 19a. Within-regimen logistic regression models ───────────────────────────
within_vars <- c("age","DUR","PatientHIV","anemia","clini_score",
                 "Timika.score","CT_mean","NEU_LYM","ALB")

within_results <- map_dfr(c("6M","8M","9M","20M"), function(reg) {
  d_reg <- na.omit(data_imp_long) %>%
    filter(Regimen == reg) %>%
    select(all_of(c("outcome_bin", within_vars)))

  n_bad <- sum(d_reg$outcome_bin == "Bad outcome", na.rm = TRUE)
  n_tot <- nrow(d_reg)

  if (n_bad < 10 || n_tot < 30) {
    return(tibble(Regimen=reg, AUC=NA, n=n_tot, n_bad=n_bad, note="too few events"))
  }

  # Simple train/test split (80/20)
  set.seed(42)
  tr_i <- createDataPartition(d_reg$outcome_bin, p = 0.8, list = FALSE)
  tr   <- d_reg[ tr_i, ]
  te   <- d_reg[-tr_i, ]

  mod <- tryCatch(
    glm(outcome_bin ~ ., data = tr, family = binomial, na.action = na.omit),
    error = function(e) NULL
  )
  if (is.null(mod)) return(tibble(Regimen=reg, AUC=NA, n=n_tot, n_bad=n_bad, note="glm failed"))

  probs <- tryCatch(predict(mod, te, type="response"), error=function(e) NULL)
  if (is.null(probs)) return(tibble(Regimen=reg, AUC=NA, n=n_tot, n_bad=n_bad, note="predict failed"))

  auc_val <- tryCatch(
    as.numeric(pROC::roc(te$outcome_bin, probs, quiet=TRUE)$auc),
    error = function(e) NA_real_
  )
  tibble(Regimen=reg, AUC=round(auc_val,3), n=n_tot, n_bad=n_bad, note="OK")
})

cat("\nWithin-regimen AUC (glm, all candidate predictors):\n")
print(within_results)

# ── 19b. Within-regimen Cox models for time-to-conversion ────────────────────
cat("\nWithin-regimen Cox HR (time to first culture conversion):\n")

for (reg in c("6M","8M","9M","20M")) {
  d_reg <- sputum %>% filter(Regimen == reg) %>%
    select(time, status, age, DUR, PatientHIV, anemia,
           clini_score, Timika.score, NEU_LYM, ALB) %>%
    na.omit()

  if (sum(d_reg$status) < 10) {
    cat(sprintf("\n%s: too few events (%d)\n", reg, sum(d_reg$status))); next
  }

  cox_reg <- tryCatch(
    coxph(Surv(time, status) ~ age + DUR + PatientHIV + anemia +
            clini_score + Timika.score + NEU_LYM + ALB,
          data = d_reg, ties = "efron"),
    error = function(e) NULL
  )
  if (!is.null(cox_reg)) {
    cat(sprintf("\n--- Cox model: Regimen %s (n=%d, events=%d) ---\n",
                reg, nrow(d_reg), sum(d_reg$status)))
    print(broom::tidy(cox_reg, exponentiate = TRUE, conf.int = TRUE) %>%
            select(term, estimate, conf.low, conf.high, p.value) %>%
            mutate(across(where(is.numeric), ~round(.x, 3))))
  }
}

# ── 19c. FINAL RECOMMENDATION SUMMARY ────────────────────────────────────────
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║           FINAL ANALYTICAL FRAMEWORK RECOMMENDATION             ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║                                                                  ║\n")
cat("║  PRIMARY ANALYSIS (treatment effect):                           ║\n")
cat("║    → eventglm::rmeanglm() — RMST regression                    ║\n")
cat("║    → Adjusted CID at tau=150d by regimen                       ║\n")
cat("║    → Reports absolute risk difference (clinical priority)       ║\n")
cat("║                                                                  ║\n")
cat("║  SECONDARY ANALYSIS (within-regimen prediction):               ║\n")
cat("║    → glm(outcome_bin ~ ...) WITHIN each Regimen                ║\n")
cat("║    → Enriched model AUC=0.68 (whole cohort) is best overall    ║\n")
cat("║    → Within-regimen models remove confounding by disease type   ║\n")
cat("║                                                                  ║\n")
cat("║  DO NOT REPORT eventglm AUC=0.704 as prediction performance:   ║\n")
cat("║    → Inflated by between-regimen definitional differences       ║\n")
cat("║    → Appropriate only for RMST/CID treatment effect estimation  ║\n")
cat("║                                                                  ║\n")
cat("║  REPORTED AUC (Table): Full=0.642 | Enriched=0.680             ║\n")
cat("║  (optimism-corrected: Full=0.641 | Enriched=0.679)             ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n")


# ── 20. SESSION INFO ──────────────────────────────────────────────────────────
cat("\n===== Session Information =====\n")
sessionInfo()
# =============================================================================
# END OF SCRIPT
# =============================================================================

