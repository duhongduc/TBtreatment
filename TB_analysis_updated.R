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
  "cmprsk", "tidycmprsk",
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
  # Parallel processing
  "future", "future.apply",
  # Misc
  "forcats", "rstatix", "psych"
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

    # ── Culture results (long timepoint sequence)
    across(c(Culture.N0, Culture.N14,
             paste0("Culture.T", sprintf("%02d", 1:24))),
           ~ factor(.x, labels = c("Negative","No sputum","Positive"))),

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

sputum <- dat %>%
  left_join(outcome_tim %>% select(studycode, STDAT), by = "studycode") %>%
  mutate(
    EnrolledDAT = as.Date(EnrolledDAT, format = "%Y/%m/%d"),
    DateDeath   = as.Date(DateDeath,   format = "%Y/%m/%d"),
    DateLastFU  = as.Date(DateLastFU,  format = "%Y/%m/%d"),

    # End date = min(regimen-end, last FU, death)
    regimen_days = case_when(
      Regimen == "6M"  ~ 180L,
      Regimen == "8M"  ~ 240L,
      Regimen == "9M"  ~ 365L,
      Regimen == "20M" ~ 730L
    ),
    EndDate = pmin(EnrolledDAT + days(regimen_days), DateLastFU, DateDeath, na.rm = TRUE),

    # Time-to-event in days
    time = as.numeric(difftime(EndDate, EnrolledDAT, units = "days")),
    time = pmax(time, 1),   # avoid zero/negative times

    # Conversion status (1 = converted, 0 = censored/no conversion)
    status = case_when(
      Culture_conversion != "censor" & !is.na(Culture_conversion) ~ 1L,
      TRUE ~ 0L
    ),

    # Competing risks: 0=censored, 1=conversion, 2=died
    status2 = case_when(
      outcome == "Died"                        ~ 2L,
      status == 1L                             ~ 1L,
      TRUE                                     ~ 0L
    ),
    status2 = factor(status2, 0:2, c("censored","conversion","died"))
  )


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

tbl_baseline <- sputum %>%
  select(all_of(c(tbl_vars,"Regimen"))) %>%
  tbl_summary(
    by      = Regimen,
    missing = "no",
    statistic = list(all_continuous()  ~ "{median} ({p25}, {p75})",
                     all_categorical() ~ "{n} ({p}%)"),
    percent = "column"
  ) %>%
  add_p(
    test = list(
      all_continuous()  ~ "kruskal.test",
      all_categorical() ~ "fisher.test.simulate.p.value"
    )
  ) %>%
  add_overall() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_header(label = "**Characteristic**") %>%
  modify_caption("**Table 2. Baseline Characteristics by Treatment Regimen**")

print(tbl_baseline)


# ── 6. CULTURE HEATMAP VISUALIZATION ─────────────────────────────────────────
culture_cols <- c("Culture.N0","Culture.N14", paste0("Culture.T",sprintf("%02d",1:24)))

df_heatmap <- sputum %>%
  select(studycode, Regimen, outcome2, all_of(culture_cols)) %>%
  pivot_longer(cols = all_of(culture_cols),
               names_to = "timepoint", values_to = "Culture") %>%
  mutate(
    timepoint = str_remove(timepoint, "Culture\\."),
    Culture   = factor(
      case_when(Culture == "Negative" ~ 0, Culture == "Positive" ~ 1, TRUE ~ 2),
      0:2, c("Negative","Positive","No sputum")
    ),
    Regimen = factor(Regimen, c("6M","8M","9M","20M"))
  )

# Ordered fakeid for Y-axis
uid_order <- sputum %>%
  arrange(Regimen, outcome2) %>%
  pull(studycode)

df_heatmap <- df_heatmap %>%
  mutate(studycode = factor(studycode, levels = uid_order))

fakeids <- data.frame(
  studycode = levels(df_heatmap$studycode),
  fakeid    = seq_along(levels(df_heatmap$studycode))
)
df_heatmap <- left_join(df_heatmap, fakeids, by = "studycode")

heatmap_plot <- ggplot(df_heatmap, aes(x = timepoint, y = desc(fakeid), fill = Culture)) +
  geom_tile(color = "white", linewidth = 0.05) +
  facet_grid(rows = vars(Regimen, outcome2), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Negative" = "#4CAF50", "Positive" = "#F44336",
                                "No sputum" = "white"),
                    na.value = "grey90") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(title    = "Tuberculosis Patient Follow-up Over Time",
       subtitle = "Heatmap of Culture Status Ordered by Regimen & Outcome",
       x = "Follow-up Timepoints", y = "Number of patients",
       fill = "Culture Status") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid    = element_blank(),
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    strip.text.y  = element_text(angle = 0, face = "bold", size = 7),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing.y = unit(0.5, "lines"),
    legend.position = "right"
  )

ggsave("figures/Culture_Heatmap_Updated.png", heatmap_plot,
       width = 30, height = 45, dpi = 300, units = "cm", bg = "white")
cat("Heatmap saved.\n")


# ── 7. KAPLAN-MEIER: TIME TO CULTURE CONVERSION ───────────────────────────────

# 7.1 KM fit
km_fit <- survfit(Surv(time, status) ~ Regimen, data = sputum)
summary(km_fit, times = c(30, 60, 120, 180, 300))

# 7.2 KM plot (inverted: probability of culture conversion = 1 - survival)
km_plot <- ggsurvplot(
  km_fit,
  fun         = "event",       # plot cumulative incidence (1 - S(t))
  data        = sputum,
  palette     = c("#2196F3","#4CAF50","#FF9800","#F44336"),
  conf.int    = TRUE,
  pval        = TRUE,
  pval.method = TRUE,
  risk.table  = TRUE,
  risk.table.col = "strata",
  xlab        = "Time (days)",
  ylab        = "Cumulative Probability of Culture Conversion",
  title       = "Kaplan-Meier: Time to Sputum Culture Conversion by Regimen",
  legend.title = "Regimen",
  break.time.by = 30,
  xlim        = c(0, 300),
  ggtheme     = theme_bw()
)

# Save KM plot
png("figures/KM_CultureConversion.png", width = 25, height = 20, res = 300, units = "cm")
print(km_plot)
dev.off()

# 7.3 Log-rank test
survdiff(Surv(time, status) ~ Regimen, data = sputum)

# 7.4 Pairwise log-rank (Bonferroni-adjusted)
pairwise_survdiff(Surv(time, status) ~ Regimen, data = sputum, p.adjust.method = "bonferroni")


# ── 8. RESTRICTED MEAN SURVIVAL TIME (RMST) ──────────────────────────────────
# RMST at tau = 300 days, reference = 20M

# Using survRM2 package
tau <- 300

rmst_result <- rmst2(
  time   = sputum$time,
  status = sputum$status,
  arm    = as.integer(sputum$Regimen == "6M"),   # 1 = 6M vs 0 = 20M
  tau    = tau
)
cat("\n===== RMST: 6M vs 20M =====\n"); print(rmst_result)

rmst_8vs20 <- rmst2(
  time   = sputum$time[sputum$Regimen %in% c("8M","20M")],
  status = sputum$status[sputum$Regimen %in% c("8M","20M")],
  arm    = as.integer(sputum$Regimen[sputum$Regimen %in% c("8M","20M")] == "8M"),
  tau    = tau
)
cat("\n===== RMST: 8M vs 20M =====\n"); print(rmst_8vs20)

rmst_9vs20 <- rmst2(
  time   = sputum$time[sputum$Regimen %in% c("9M","20M")],
  status = sputum$status[sputum$Regimen %in% c("9M","20M")],
  arm    = as.integer(sputum$Regimen[sputum$Regimen %in% c("9M","20M")] == "9M"),
  tau    = tau
)
cat("\n===== RMST: 9M vs 20M =====\n"); print(rmst_9vs20)


# ── 9. COMPETING RISKS: FINE-GRAY MODEL ──────────────────────────────────────

# 9.1 Cumulative incidence (conversion, death as competing event)
cuminc_fit <- cuminc(
  ftime   = sputum$time,
  fstatus = as.numeric(sputum$status2) - 1,  # 0=censor,1=conv,2=died
  group   = sputum$Regimen
)

# Plot CIF
plot(cuminc_fit,
     curvlab = levels(sputum$Regimen),
     main = "Cumulative Incidence of Culture Conversion (Competing: Death)",
     xlab = "Days", ylab = "Cumulative Incidence",
     col  = c("#2196F3","#4CAF50","#FF9800","#F44336",
              "#9C27B0","#795548","#009688","#607D8B"))

# 9.2 Unadjusted Fine-Gray model
# Note: crr() uses factor coding; create dummy for Regimen (reference = 20M)
sputum_cr <- sputum %>%
  filter(!is.na(time) & !is.na(status2)) %>%
  mutate(
    status_cr = as.integer(status2) - 1,  # 0=censor,1=conv,2=died
    R6  = as.integer(Regimen == "6M"),
    R8  = as.integer(Regimen == "8M"),
    R9  = as.integer(Regimen == "9M")
  )

cov_cr <- model.matrix(~ R6 + R8 + R9, data = sputum_cr)[, -1]

fg_unadj <- crr(
  ftime   = sputum_cr$time,
  fstatus = sputum_cr$status_cr,
  cov1    = cov_cr,
  failcode = 1
)
summary(fg_unadj)

# 9.3 Adjusted Fine-Gray model (add key covariates)
sputum_cr2 <- sputum_cr %>%
  mutate(
    age_c         = scale(age),
    hiv_pos       = as.integer(PatientHIV == "Positive"),
    anemia_yes    = as.integer(anemia == "Anemia"),
    neu_lym_c     = scale(NEU_LYM),
    timika_c      = scale(Timika.score),
    dur_c         = scale(DUR)
  ) %>%
  filter(complete.cases(age_c, hiv_pos, anemia_yes, neu_lym_c, timika_c, dur_c))

cov_adj <- model.matrix(
  ~ R6 + R8 + R9 + age_c + hiv_pos + anemia_yes + neu_lym_c + timika_c + dur_c,
  data = sputum_cr2
)[, -1]

fg_adj <- crr(
  ftime    = sputum_cr2$time,
  fstatus  = sputum_cr2$status_cr,
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
ggforest(cox_multi, data = sputum,
         main = "Hazard Ratios: Time to Culture Conversion")


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
n_imp <- 50

imp <- mice(
  data    = data_pred %>% select(-outcome_bin),
  m       = n_imp,
  maxit   = 10,
  method  = "pmm",
  printFlag = FALSE,
  seed    = 1234
)

# Stack all imputed datasets (for rms display)
data_imp_long <- complete(imp, "long", include = FALSE)
data_imp_long$outcome_bin <- rep(data_pred$outcome_bin[complete.cases(data_pred)],
                                  n_imp)


# ── 13. FULL LOGISTIC REGRESSION (pooled via mice) ────────────────────────────
fit_full_mice <- with(
  data = mice::complete(imp, "long") %>%
    mutate(.imp = .imp),
  expr = glm(outcome_bin ~ age + bmi + DUR + TBTreatedBefore + DM +
               PatientHIV + Regimen + anemia + clini_score +
               Timika.score + CT_mean + NEU_LYM + ALB,
             family = binomial)
)
pooled_full <- pool(fit_full_mice)
cat("\n===== Pooled Full Logistic Regression =====\n")
summary(pooled_full, conf.int = TRUE, exponentiate = TRUE)


# ── 14. LASSO VARIABLE SELECTION (MI + Bootstrap) ────────────────────────────
n_imp_lasso <- 50
n_boot      <- 200
plan(multisession, workers = min(parallel::detectCores() - 1, 6))

results_lasso <- vector("list", n_imp_lasso)

for (i in seq_len(n_imp_lasso)) {
  data_i <- complete(imp, action = i, include = FALSE)

  results_boot <- vector("list", n_boot)
  for (j in seq_len(n_boot)) {
    results_boot[[j]] <- future({
      idx        <- sample(nrow(data_i), replace = TRUE)
      data_boot  <- data_i[idx, ]
      train_idx  <- createDataPartition(data_boot$outcome_bin, p = 0.8, list = FALSE)
      train_data <- na.omit(data_boot[train_idx, ])
      test_data  <- na.omit(data_boot[-train_idx, ])

      x <- model.matrix(outcome_bin ~ age + bmi + DUR + TBTreatedBefore + DM +
                           PatientHIV + Regimen + anemia + clini_score +
                           Timika.score + CT_mean + NEU_LYM + ALB,
                         data = train_data)[, -1]
      y <- ifelse(train_data$outcome_bin == "Bad outcome", 1, 0)

      cv_lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial", nfolds = 10)
      lasso_mod <- glmnet(x, y, alpha = 1, family = "binomial",
                           lambda = cv_lasso$lambda.min)

      x_test  <- model.matrix(outcome_bin ~ age + bmi + DUR + TBTreatedBefore + DM +
                                 PatientHIV + Regimen + anemia + clini_score +
                                 Timika.score + CT_mean + NEU_LYM + ALB,
                               data = test_data)[, -1]
      probs   <- predict(lasso_mod, newx = x_test, type = "response")
      obs     <- test_data$outcome_bin

      roc1  <- pROC::roc(obs, as.vector(probs), quiet = TRUE)
      brier <- mean((ifelse(obs == "Bad outcome", 1, 0) - as.vector(probs))^2)

      sel_vars <- rownames(coef(lasso_mod))[as.vector(coef(lasso_mod)) != 0]
      sel_vars <- sel_vars[sel_vars != "(Intercept)"]

      list(sel_vars = sel_vars, roc1 = roc1, brier = brier)
    }, gc = TRUE, seed = TRUE)
  }
  results_lasso[[i]] <- lapply(results_boot, value)
  cat(sprintf("Imputation %d/%d done\n", i, n_imp_lasso))
}

plan(sequential)

# Variable inclusion frequency
var_counts <- table(unlist(lapply(results_lasso, function(x)
  unlist(lapply(x, `[[`, "sel_vars")))))
var_freq   <- var_counts / (n_imp_lasso * n_boot)
cat("\nLASSO variable inclusion frequencies:\n"); print(sort(var_freq, decreasing = TRUE))

selected_vars_lasso <- names(var_freq[var_freq > 0.50])
cat("\nSelected (>50%):", selected_vars_lasso, "\n")

# Pool ROC & Brier
all_rocs   <- unlist(lapply(unlist(results_lasso, recursive = FALSE),
                             function(x) as.numeric(x$roc1$auc)))
all_briers <- unlist(lapply(unlist(results_lasso, recursive = FALSE), `[[`, "brier"))

cat(sprintf("Mean AUC (LASSO):   %.3f\n", mean(all_rocs, na.rm = TRUE)))
cat(sprintf("Mean Brier (LASSO): %.3f\n", mean(all_briers, na.rm = TRUE)))

save(results_lasso, selected_vars_lasso,
     file = "output/LASSO_results.RData")


# ── 15. FINAL LASSO MODEL (display, validate, calibrate) ──────────────────────
dd <- datadist(data_imp_long); options(datadist = "dd")

lasso_formula <- as.formula(paste(
  "outcome_bin ~",
  paste(selected_vars_lasso, collapse = " + ")
))

f_lasso <- rms::lrm(lasso_formula,
                     data  = na.omit(data_imp_long),
                     x = TRUE, y = TRUE)

cat("\n===== rms LRM (LASSO-selected) =====\n"); print(f_lasso)

# Nomogram
png("figures/Nomogram_LASSO.png", width = 2400, height = 1200, res = 150)
plot(nomogram(f_lasso,
              lp      = FALSE,
              fun     = plogis,
              fun.at  = c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95),
              funlabel = "Probability of Unfavorable Outcome"))
dev.off()

# Bootstrap validation
val_lasso <- validate(f_lasso, B = 200)
print(val_lasso)
c_lasso_optcorr <- 0.5 * (val_lasso[1, 5] + 1)
cat(sprintf("Optimism-corrected C-statistic (LASSO): %.3f\n", c_lasso_optcorr))

# Calibration
cal_lasso <- rms::calibrate(f_lasso, B = 200)
png("figures/Calibration_LASSO.png", width = 800, height = 600)
plot(cal_lasso, main = "Calibration: LASSO Model")
dev.off()


# ── 16. AIC STEPWISE MODEL ────────────────────────────────────────────────────
n_imp_aic <- 50
n_boot_aic <- 200
plan(multisession, workers = min(parallel::detectCores() - 1, 6))

results_aic <- vector("list", n_imp_aic)
for (i in seq_len(n_imp_aic)) {
  data_i <- complete(imp, action = i, include = FALSE)

  results_boot_aic <- vector("list", n_boot_aic)
  for (j in seq_len(n_boot_aic)) {
    results_boot_aic[[j]] <- future({
      idx       <- sample(nrow(data_i), replace = TRUE)
      data_boot <- data_i[idx, ]

      full_mod <- tryCatch(
        glm(outcome_bin ~ rcs(age, 3) + rcs(bmi, 3) + rcs(DUR, 3) +
              TBTreatedBefore + DM + PatientHIV + Regimen + anemia +
              rcs(clini_score, 3) + rcs(Timika.score, 3) +
              rcs(CT_mean, 3) + rcs(NEU_LYM, 3) + rcs(ALB, 3),
            data   = na.omit(data_boot),
            family = binomial(link = "logit")),
        error = function(e) NULL
      )
      if (is.null(full_mod)) return(NULL)

      step_mod  <- stepAIC(full_mod, direction = "backward", trace = FALSE)
      sel_vars  <- names(coef(step_mod))[-1]
      list(sel_vars = sel_vars)
    }, gc = TRUE, seed = TRUE)
  }
  results_aic[[i]] <- lapply(results_boot_aic, value)
  cat(sprintf("AIC imputation %d/%d done\n", i, n_imp_aic))
}
plan(sequential)

var_counts_aic <- table(unlist(lapply(results_aic, function(x)
  unlist(lapply(x, function(y) if (!is.null(y)) y$sel_vars else NULL)))))
var_freq_aic   <- var_counts_aic / (n_imp_aic * n_boot_aic)
cat("\nAIC variable inclusion frequencies:\n"); print(sort(var_freq_aic, decreasing = TRUE))
selected_vars_aic <- names(var_freq_aic[var_freq_aic > 0.50])
cat("\nAIC selected (>50%):", selected_vars_aic, "\n")

# Final AIC model (rms)
aic_formula <- as.formula(paste(
  "outcome_bin ~",
  paste(selected_vars_aic, collapse = " + ")
))

f_aic <- rms::lrm(aic_formula,
                   data  = na.omit(data_imp_long),
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
png("figures/Calibration_AIC.png", width = 800, height = 600)
plot(cal_aic, main = "Calibration: AIC Model")
dev.off()

save(selected_vars_aic, results_aic, f_aic, val_aic, cal_aic,
     file = "output/AIC_results.RData")


# ── 17. ROC CURVE COMPARISON ──────────────────────────────────────────────────
# Quick single-split comparison for illustration
set.seed(123)
data_cc <- na.omit(data_pred)
train_idx_roc <- createDataPartition(data_cc$outcome_bin, p = 0.8, list = FALSE)
train_roc     <- data_cc[train_idx_roc, ]
test_roc      <- data_cc[-train_idx_roc, ]

# Full model
full_roc_mod  <- glm(outcome_bin ~ age + bmi + DUR + TBTreatedBefore + DM +
                       PatientHIV + Regimen + anemia + clini_score +
                       Timika.score + CT_mean + NEU_LYM + ALB,
                     data = train_roc, family = binomial, na.action = na.omit)
probs_full <- predict(full_roc_mod, test_roc, type = "response")
obs_roc    <- test_roc$outcome_bin

roc_full  <- pROC::roc(obs_roc, probs_full, quiet = TRUE)

# LASSO model (lambda.min)
x_train_r <- model.matrix(outcome_bin ~ age + bmi + DUR + TBTreatedBefore + DM +
                             PatientHIV + Regimen + anemia + clini_score +
                             Timika.score + CT_mean + NEU_LYM + ALB,
                           data = train_roc)[, -1]
y_train_r <- ifelse(train_roc$outcome_bin == "Bad outcome", 1, 0)
cv_r <- cv.glmnet(x_train_r, y_train_r, alpha = 1, family = "binomial")
lasso_r <- glmnet(x_train_r, y_train_r, alpha = 1, family = "binomial",
                   lambda = cv_r$lambda.min)
x_test_r  <- model.matrix(outcome_bin ~ age + bmi + DUR + TBTreatedBefore + DM +
                             PatientHIV + Regimen + anemia + clini_score +
                             Timika.score + CT_mean + NEU_LYM + ALB,
                           data = test_roc)[, -1]
probs_lasso <- as.vector(predict(lasso_r, newx = x_test_r, type = "response"))
roc_lasso   <- pROC::roc(obs_roc, probs_lasso, quiet = TRUE)

# Plot
png("figures/ROC_Comparison.png", width = 800, height = 700)
plot(roc_full,  col = "#2196F3", lwd = 2, print.auc = TRUE, print.auc.y = 0.45,
     main = "ROC Curves: Full vs. LASSO Logistic Regression")
plot(roc_lasso, col = "#F44336", lwd = 2, add = TRUE, print.auc = TRUE, print.auc.y = 0.38)
legend("bottomright",
       legend = c(sprintf("Full Model (AUC=%.3f)", as.numeric(roc_full$auc)),
                  sprintf("LASSO     (AUC=%.3f)", as.numeric(roc_lasso$auc))),
       col = c("#2196F3","#F44336"), lwd = 2)
dev.off()

# DeLong test for AUC comparison
roc.test(roc_full, roc_lasso, reuse.auc = FALSE)


# ── 18. SESSION INFO ──────────────────────────────────────────────────────────
cat("\n===== Session Information =====\n")
sessionInfo()
# =============================================================================
# END OF SCRIPT
# =============================================================================
