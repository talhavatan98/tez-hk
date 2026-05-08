#!/usr/bin/env Rscript
# ════════════════════════════════════════════════════════════════════════════
#  COMPLETE CLINICAL ANALYSIS
#  Using DIABLO samples (n=561) + Clinical data merge
#  This is the CORRECT approach for publication
# ════════════════════════════════════════════════════════════════════════════

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║  CLINICAL ANALYSIS - COMPLETE (DIABLO Samples)                ║\n")
cat("║  Table 1 for Publication                                      ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Disable OpenMP
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)

library(survival)
set.seed(42)

# Output directory
OUTPUT_DIR <- "Section_4_1_Clinical_Cohort"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ═══════════════════════════════════════════════════════════════════════════
# 1. LOAD IDH STATUS (FROM DIABLO - OUR ANALYSIS SAMPLES)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 1: Loading IDH Status from DIABLO Analysis\n")
cat("─────────────────────────────────────────────────────────────────\n")

load("DIABLO_ready_data.RData")
Y_idh <- Y # IDH status
sample_ids <- names(Y_idh)

cat(sprintf("✓ Loaded DIABLO samples: %d\n", length(sample_ids)))
cat("  IDH distribution:\n")
print(table(Y_idh))
cat("\n")

# Create IDH dataframe
idh_data <- data.frame(
  sampleID = sample_ids,
  IDH_status = as.character(Y_idh),
  stringsAsFactors = FALSE
)

# ═══════════════════════════════════════════════════════════════════════════
# 2. LOAD CLINICAL DATA
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 2: Loading Clinical Matrix\n")
cat("─────────────────────────────────────────────────────────────────\n")

clinical <- read.delim("TCGA.GBMLGG.sampleMap_GBMLGG_clinicalMatrix",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

cat(sprintf(
  "✓ Loaded clinical data: %d samples, %d variables\n\n",
  nrow(clinical), ncol(clinical)
))

# ═══════════════════════════════════════════════════════════════════════════
# 3. MERGE DIABLO SAMPLES WITH CLINICAL DATA
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 3: Merging DIABLO Samples with Clinical Data\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Merge by sampleID
merged <- merge(
  idh_data,
  clinical[, c(
    "sampleID", "age_at_initial_pathologic_diagnosis", "gender",
    "histological_type", "neoplasm_histologic_grade", "vital_status",
    "days_to_death", "days_to_last_followup"
  )],
  by = "sampleID",
  all.x = TRUE # Keep all DIABLO samples
)

cat(sprintf("✓ Merged %d DIABLO samples with clinical data\n", nrow(merged)))
cat(sprintf(
  "  • Samples with age data: %d\n",
  sum(!is.na(as.numeric(merged$age_at_initial_pathologic_diagnosis)))
))
cat(sprintf(
  "  • Samples with gender data: %d\n",
  sum(trimws(merged$gender) != "", na.rm = TRUE)
))
cat(sprintf(
  "  • Samples with grade data: %d\n",
  sum(trimws(merged$neoplasm_histologic_grade) != "", na.rm = TRUE)
))
cat("\n")

# ═══════════════════════════════════════════════════════════════════════════
# 4. CLEAN AND PREPARE DATA
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 4: Cleaning and Preparing Data\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Age
merged$age <- as.numeric(merged$age_at_initial_pathologic_diagnosis)
merged$age_group <- cut(merged$age,
  breaks = c(0, 40, 60, 100),
  labels = c("<40", "40-60", ">60"),
  include.lowest = TRUE
)

# Gender
merged$gender <- toupper(trimws(merged$gender))
merged$gender[merged$gender == ""] <- NA

# Tumor type
merged$tumor_type <- ifelse(
  grepl("glioblastoma|GBM", merged$histological_type, ignore.case = TRUE),
  "GBM", "LGG"
)

# Grade
merged$grade <- trimws(merged$neoplasm_histologic_grade)
merged$grade[merged$grade == ""] <- NA
merged$grade[merged$tumor_type == "GBM" & is.na(merged$grade)] <- "G4"
merged$grade_num <- gsub("^G", "", merged$grade)

# Survival
merged$vital_status <- toupper(trimws(merged$vital_status))
merged$survival_days <- ifelse(
  merged$vital_status == "DECEASED",
  as.numeric(merged$days_to_death),
  as.numeric(merged$days_to_last_followup)
)
merged$survival_months <- merged$survival_days / 30.44
merged$survival_status <- ifelse(merged$vital_status == "DECEASED", 1, 0)

cat("✓ Data cleaned\n")
cat(sprintf(
  "  • Valid age: %d (%.1f%%)\n",
  sum(!is.na(merged$age)), 100 * sum(!is.na(merged$age)) / nrow(merged)
))
cat(sprintf(
  "  • Valid gender: %d (%.1f%%)\n",
  sum(!is.na(merged$gender)), 100 * sum(!is.na(merged$gender)) / nrow(merged)
))
cat(sprintf(
  "  • Valid grade: %d (%.1f%%)\n",
  sum(!is.na(merged$grade_num)), 100 * sum(!is.na(merged$grade_num)) / nrow(merged)
))
cat(sprintf(
  "  • Valid survival: %d (%.1f%%)\n\n",
  sum(!is.na(merged$survival_days) & merged$survival_days > 0),
  100 * sum(!is.na(merged$survival_days) & merged$survival_days > 0) / nrow(merged)
))

# ═══════════════════════════════════════════════════════════════════════════
# 5. CREATE TABLE 1
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 5: Creating Table 1 (Demographics by IDH Status)\n")
cat("─────────────────────────────────────────────────────────────────\n\n")

mutant_data <- merged[merged$IDH_status == "Mutant", ]
wt_data <- merged[merged$IDH_status == "WT", ]

# Helper function for percentages
calc_pct <- function(x, total) {
  sprintf("%d (%.1f%%)", x, 100 * x / total)
}

# Create table
table1 <- data.frame(
  Variable = character(),
  Overall = character(),
  IDH_Mutant = character(),
  IDH_WT = character(),
  stringsAsFactors = FALSE
)

# Sample size
table1 <- rbind(table1, data.frame(
  Variable = "N",
  Overall = as.character(nrow(merged)),
  IDH_Mutant = as.character(nrow(mutant_data)),
  IDH_WT = as.character(nrow(wt_data))
))

# Age
n_age_all <- sum(!is.na(merged$age))
n_age_mut <- sum(!is.na(mutant_data$age))
n_age_wt <- sum(!is.na(wt_data$age))

table1 <- rbind(table1, data.frame(
  Variable = sprintf("Age (mean ± SD), n=%d", n_age_all),
  Overall = sprintf("%.1f ± %.1f", mean(merged$age, na.rm = TRUE), sd(merged$age, na.rm = TRUE)),
  IDH_Mutant = sprintf("%.1f ± %.1f", mean(mutant_data$age, na.rm = TRUE), sd(mutant_data$age, na.rm = TRUE)),
  IDH_WT = sprintf("%.1f ± %.1f", mean(wt_data$age, na.rm = TRUE), sd(wt_data$age, na.rm = TRUE))
))

# Age groups
for (ag in c("<40", "40-60", ">60")) {
  n_ag_all <- sum(!is.na(merged$age_group))
  table1 <- rbind(table1, data.frame(
    Variable = sprintf("  Age %s", ag),
    Overall = calc_pct(sum(merged$age_group == ag, na.rm = TRUE), n_ag_all),
    IDH_Mutant = calc_pct(
      sum(mutant_data$age_group == ag, na.rm = TRUE),
      sum(!is.na(mutant_data$age_group))
    ),
    IDH_WT = calc_pct(
      sum(wt_data$age_group == ag, na.rm = TRUE),
      sum(!is.na(wt_data$age_group))
    )
  ))
}

# Gender
n_gender_all <- sum(!is.na(merged$gender))
table1 <- rbind(table1, data.frame(
  Variable = sprintf("Gender, n=%d", n_gender_all),
  Overall = "",
  IDH_Mutant = "",
  IDH_WT = ""
))

for (gend in c("FEMALE", "MALE")) {
  table1 <- rbind(table1, data.frame(
    Variable = sprintf("  %s", gend),
    Overall = calc_pct(sum(merged$gender == gend, na.rm = TRUE), n_gender_all),
    IDH_Mutant = calc_pct(
      sum(mutant_data$gender == gend, na.rm = TRUE),
      sum(!is.na(mutant_data$gender))
    ),
    IDH_WT = calc_pct(
      sum(wt_data$gender == gend, na.rm = TRUE),
      sum(!is.na(wt_data$gender))
    )
  ))
}

# Tumor type
table1 <- rbind(table1, data.frame(
  Variable = "Tumor Type",
  Overall = "",
  IDH_Mutant = "",
  IDH_WT = ""
))

for (tt in c("GBM", "LGG")) {
  table1 <- rbind(table1, data.frame(
    Variable = sprintf("  %s", tt),
    Overall = calc_pct(sum(merged$tumor_type == tt, na.rm = TRUE), nrow(merged)),
    IDH_Mutant = calc_pct(sum(mutant_data$tumor_type == tt, na.rm = TRUE), nrow(mutant_data)),
    IDH_WT = calc_pct(sum(wt_data$tumor_type == tt, na.rm = TRUE), nrow(wt_data))
  ))
}

# Grade
n_grade_all <- sum(!is.na(merged$grade_num))
table1 <- rbind(table1, data.frame(
  Variable = sprintf("Grade, n=%d", n_grade_all),
  Overall = "",
  IDH_Mutant = "",
  IDH_WT = ""
))

for (gr in c("2", "3", "4")) {
  table1 <- rbind(table1, data.frame(
    Variable = sprintf("  Grade %s", gr),
    Overall = calc_pct(sum(merged$grade_num == gr, na.rm = TRUE), n_grade_all),
    IDH_Mutant = calc_pct(
      sum(mutant_data$grade_num == gr, na.rm = TRUE),
      sum(!is.na(mutant_data$grade_num))
    ),
    IDH_WT = calc_pct(
      sum(wt_data$grade_num == gr, na.rm = TRUE),
      sum(!is.na(wt_data$grade_num))
    )
  ))
}

# Vital status
n_vital_all <- sum(!is.na(merged$vital_status))
table1 <- rbind(table1, data.frame(
  Variable = sprintf("Vital Status, n=%d", n_vital_all),
  Overall = "",
  IDH_Mutant = "",
  IDH_WT = ""
))

table1 <- rbind(table1, data.frame(
  Variable = "  Deceased",
  Overall = calc_pct(sum(merged$vital_status == "DECEASED", na.rm = TRUE), n_vital_all),
  IDH_Mutant = calc_pct(
    sum(mutant_data$vital_status == "DECEASED", na.rm = TRUE),
    sum(!is.na(mutant_data$vital_status))
  ),
  IDH_WT = calc_pct(
    sum(wt_data$vital_status == "DECEASED", na.rm = TRUE),
    sum(!is.na(wt_data$vital_status))
  )
))

# Survival (deceased only)
deceased_all <- merged[merged$vital_status == "DECEASED" & !is.na(merged$survival_months), ]
deceased_mut <- mutant_data[mutant_data$vital_status == "DECEASED" & !is.na(mutant_data$survival_months), ]
deceased_wt <- wt_data[wt_data$vital_status == "DECEASED" & !is.na(wt_data$survival_months), ]

if (nrow(deceased_all) > 0) {
  table1 <- rbind(table1, data.frame(
    Variable = sprintf("Survival (months), median [IQR], n=%d deceased", nrow(deceased_all)),
    Overall = sprintf(
      "%.1f [%.1f-%.1f]",
      median(deceased_all$survival_months, na.rm = TRUE),
      quantile(deceased_all$survival_months, 0.25, na.rm = TRUE),
      quantile(deceased_all$survival_months, 0.75, na.rm = TRUE)
    ),
    IDH_Mutant = sprintf(
      "%.1f [%.1f-%.1f]",
      median(deceased_mut$survival_months, na.rm = TRUE),
      quantile(deceased_mut$survival_months, 0.25, na.rm = TRUE),
      quantile(deceased_mut$survival_months, 0.75, na.rm = TRUE)
    ),
    IDH_WT = sprintf(
      "%.1f [%.1f-%.1f]",
      median(deceased_wt$survival_months, na.rm = TRUE),
      quantile(deceased_wt$survival_months, 0.25, na.rm = TRUE),
      quantile(deceased_wt$survival_months, 0.75, na.rm = TRUE)
    )
  ))
}

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                    TABLE 1: DEMOGRAPHICS                       ║\n")
cat("║                  (DIABLO Analysis Cohort)                      ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

print(table1, row.names = FALSE)

write.csv(table1, file.path(OUTPUT_DIR, "Table1_Demographics_DIABLO_Cohort.csv"), row.names = FALSE)
cat("\n✓ Table 1 saved\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# 6. STATISTICAL TESTS
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 6: Statistical Tests (Clinical Variables vs IDH Status)\n")
cat("─────────────────────────────────────────────────────────────────\n")

stat_tests <- data.frame(
  Variable = character(),
  Test = character(),
  N_analyzed = character(),
  Statistic = character(),
  P_value = character(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Age: t-test
age_data <- merged[!is.na(merged$age), ]
if (nrow(age_data) > 0) {
  age_test <- t.test(age ~ IDH_status, data = age_data)
  stat_tests <- rbind(stat_tests, data.frame(
    Variable = "Age",
    Test = "Two-sample t-test",
    N_analyzed = as.character(nrow(age_data)),
    Statistic = sprintf("t = %.3f", age_test$statistic),
    P_value = sprintf("%.4e", age_test$p.value),
    Significant = ifelse(age_test$p.value < 0.05, "Yes", "No")
  ))
}

# Gender: Chi-square
gender_data <- merged[!is.na(merged$gender), ]
if (nrow(gender_data) > 0) {
  gender_tab <- table(gender_data$gender, gender_data$IDH_status)
  if (min(gender_tab) >= 5) {
    gender_test <- chisq.test(gender_tab)
    stat_tests <- rbind(stat_tests, data.frame(
      Variable = "Gender",
      Test = "Chi-square test",
      N_analyzed = as.character(nrow(gender_data)),
      Statistic = sprintf("χ² = %.3f", gender_test$statistic),
      P_value = sprintf("%.4e", gender_test$p.value),
      Significant = ifelse(gender_test$p.value < 0.05, "Yes", "No")
    ))
  }
}

# Tumor type: Chi-square
tumor_tab <- table(merged$tumor_type, merged$IDH_status)
if (min(tumor_tab) >= 5) {
  tumor_test <- chisq.test(tumor_tab)
  stat_tests <- rbind(stat_tests, data.frame(
    Variable = "Tumor_Type",
    Test = "Chi-square test",
    N_analyzed = as.character(nrow(merged)),
    Statistic = sprintf("χ² = %.3f", tumor_test$statistic),
    P_value = sprintf("%.4e", tumor_test$p.value),
    Significant = ifelse(tumor_test$p.value < 0.05, "Yes", "No")
  ))
}

# Grade: Chi-square
grade_data <- merged[!is.na(merged$grade_num), ]
if (nrow(grade_data) > 0) {
  grade_tab <- table(grade_data$grade_num, grade_data$IDH_status)
  if (min(dim(grade_tab)) > 0 && min(grade_tab) >= 5) {
    grade_test <- chisq.test(grade_tab)
    stat_tests <- rbind(stat_tests, data.frame(
      Variable = "Grade",
      Test = "Chi-square test",
      N_analyzed = as.character(nrow(grade_data)),
      Statistic = sprintf("χ² = %.3f", grade_test$statistic),
      P_value = sprintf("%.4e", grade_test$p.value),
      Significant = ifelse(grade_test$p.value < 0.05, "Yes", "No")
    ))
  }
}

cat("\nStatistical Test Results:\n")
print(stat_tests, row.names = FALSE)

write.csv(stat_tests, file.path(OUTPUT_DIR, "Statistical_Tests_DIABLO_Cohort.csv"), row.names = FALSE)
cat("\n✓ Statistical tests saved\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# 7. SURVIVAL ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 7: Survival Analysis (Kaplan-Meier)\n")
cat("─────────────────────────────────────────────────────────────────\n")

surv_data <- merged[
  !is.na(merged$survival_days) &
    merged$survival_days > 0 &
    !is.na(merged$survival_status),
]

cat(sprintf(
  "Survival data available: %d samples (%.1f%% of cohort)\n",
  nrow(surv_data), 100 * nrow(surv_data) / nrow(merged)
))

if (nrow(surv_data) > 10) {
  surv_obj <- Surv(time = surv_data$survival_months, event = surv_data$survival_status)
  km_fit <- survfit(surv_obj ~ IDH_status, data = surv_data)

  cat("\nMedian Survival by IDH Status:\n")
  print(km_fit)

  # Log-rank test
  logrank_test <- survdiff(surv_obj ~ IDH_status, data = surv_data)
  logrank_pval <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)

  if (logrank_pval < 2e-16) {
    p_display <- "< 2e-16"
  } else {
    p_display <- sprintf("%.4e", logrank_pval)
  }

  cat(sprintf("\n\nLog-rank test p-value: %s\n", p_display))
  if (logrank_pval < 0.001) {
    cat("★★★ HIGHLY SIGNIFICANT survival difference (p < 0.001)\n")
  } else if (logrank_pval < 0.05) {
    cat("★★ SIGNIFICANT survival difference (p < 0.05)\n")
  } else {
    cat("Not significant (p ≥ 0.05)\n")
  }

  # Plot
  pdf(file.path(OUTPUT_DIR, "Kaplan_Meier_Survival_Curves.pdf"), width = 10, height = 7)
  par(mar = c(5, 5, 4, 2))
  plot(km_fit,
    col = c("#E63946", "#457B9D"),
    lwd = 3,
    xlab = "Time (months)",
    ylab = "Survival Probability",
    main = sprintf(
      "Kaplan-Meier Survival Curves (DIABLO Cohort)\n(Log-rank p %s)",
      ifelse(grepl("<", p_display), p_display, paste("=", p_display))
    ),
    cex.main = 1.5,
    cex.lab = 1.3,
    cex.axis = 1.1
  )

  legend("topright",
    legend = c(
      sprintf("IDH Mutant (n=%d)", sum(surv_data$IDH_status == "Mutant")),
      sprintf("IDH WT (n=%d)", sum(surv_data$IDH_status == "WT"))
    ),
    col = c("#E63946", "#457B9D"),
    lwd = 3,
    bty = "n",
    cex = 1.2
  )

  dev.off()
  cat("\n✓ Kaplan-Meier curves saved\n")

  # Cox regression
  cox_uni <- coxph(surv_obj ~ IDH_status, data = surv_data)
  cox_summary <- summary(cox_uni)

  cox_results <- data.frame(
    Variable = "IDH_status_WT_vs_Mutant",
    N_events = sum(surv_data$survival_status),
    HR = exp(coef(cox_uni)),
    CI_lower = cox_summary$conf.int[, "lower .95"],
    CI_upper = cox_summary$conf.int[, "upper .95"],
    P_value = cox_summary$coefficients[, "Pr(>|z|)"]
  )

  cat("\nCox Regression (Univariate):\n")
  cat(sprintf(
    "  Hazard Ratio (WT vs Mutant): %.3f [95%% CI: %.3f-%.3f]\n",
    cox_results$HR, cox_results$CI_lower, cox_results$CI_upper
  ))
  cat(sprintf("  P-value: %.4e\n", cox_results$P_value))
  cat(sprintf("  Events: %d\n", cox_results$N_events))

  write.csv(cox_results, file.path(OUTPUT_DIR, "Cox_Regression_DIABLO_Cohort.csv"), row.names = FALSE)
  cat("\n✓ Cox regression saved\n")
} else {
  cat("⚠ Insufficient survival data\n")
}

cat("\n")

# ═══════════════════════════════════════════════════════════════════════════
# 8. DATA COVERAGE SUMMARY
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 8: Data Coverage Summary\n")
cat("─────────────────────────────────────────────────────────────────\n")

coverage <- data.frame(
  Variable = c("IDH Status", "Age", "Gender", "Tumor Type", "Grade", "Vital Status", "Survival Data"),
  N_Available = c(
    nrow(merged),
    sum(!is.na(merged$age)),
    sum(!is.na(merged$gender)),
    nrow(merged),
    sum(!is.na(merged$grade_num)),
    sum(!is.na(merged$vital_status)),
    sum(!is.na(merged$survival_days) & merged$survival_days > 0)
  ),
  Percent_Coverage = c(
    100.0,
    100 * sum(!is.na(merged$age)) / nrow(merged),
    100 * sum(!is.na(merged$gender)) / nrow(merged),
    100.0,
    100 * sum(!is.na(merged$grade_num)) / nrow(merged),
    100 * sum(!is.na(merged$vital_status)) / nrow(merged),
    100 * sum(!is.na(merged$survival_days) & merged$survival_days > 0) / nrow(merged)
  )
)

coverage$Percent_Coverage <- sprintf("%.1f%%", coverage$Percent_Coverage)

cat("\nData Coverage (n=561 DIABLO samples):\n")
print(coverage, row.names = FALSE)

write.csv(coverage, file.path(OUTPUT_DIR, "Data_Coverage_Summary.csv"), row.names = FALSE)
cat("\n✓ Coverage summary saved\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# 9. FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                  ANALYSIS COMPLETE                             ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("COHORT SUMMARY:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  • Total DIABLO samples: %d\n", nrow(merged)))
cat(sprintf(
  "  • IDH Mutant: %d (%.1f%%)\n",
  sum(merged$IDH_status == "Mutant"),
  100 * sum(merged$IDH_status == "Mutant") / nrow(merged)
))
cat(sprintf(
  "  • IDH WT: %d (%.1f%%)\n\n",
  sum(merged$IDH_status == "WT"),
  100 * sum(merged$IDH_status == "WT") / nrow(merged)
))

cat("FILES SAVED:\n")
cat(sprintf("  📁 %s/\n", OUTPUT_DIR))
cat("    • Table1_Demographics_DIABLO_Cohort.csv\n")
cat("    • Statistical_Tests_DIABLO_Cohort.csv\n")
cat("    • Kaplan_Meier_Survival_Curves.pdf\n")
cat("    • Cox_Regression_DIABLO_Cohort.csv\n")
cat("    • Data_Coverage_Summary.csv\n\n")

cat("✓ READY FOR MANUSCRIPT\n")
cat("  This Table 1 represents the DIABLO analysis cohort (n=561)\n")
cat("  Use this for consistency with all other analyses\n\n")
