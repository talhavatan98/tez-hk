#!/usr/bin/env Rscript
# ════════════════════════════════════════════════════════════════════════════
#  STATISTICAL METHOD COMPARISON
#  Bootstrap-based pairwise comparisons between all methods
#  Generates p-values for publication
# ════════════════════════════════════════════════════════════════════════════

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║  STATISTICAL METHOD COMPARISON (Bootstrap)                     ║\n")
cat("║  Pairwise comparisons for publication                          ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

suppressPackageStartupMessages({
  library(mixOmics)  # For predict() on DIABLO and sPLS-DA models
  library(pROC)
  library(caret)  # For RF predict()
})

set.seed(42)  # Reproducible bootstrap

# ═══════════════════════════════════════════════════════════════════════════
# 1. LOAD TEST SET PREDICTIONS
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 1: Loading Test Set Predictions\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Load workspace files that contain test predictions
load("DIABLO_Final_Publication_Fixed/Model_Workspace.RData")
Y_test_diablo <- Y_test
# DIABLO uses weighted average of block predictions
diablo_pred <- predict(final_model, newdata = X_test)
pred_scores_mrna <- diablo_pred$predict$mRNA[, "Mutant", 2]  # component 2
pred_scores_meth <- diablo_pred$predict$Methylation[, "Mutant", 2]
diablo_scores <- (pred_scores_mrna + pred_scores_meth) / 2

# Single-omics mRNA
load("Single_Omics_mRNA_Results/Model_Workspace.RData")
Y_test_mrna <- Y_test
mrna_pred <- predict(final_model, newdata = X_test)
mrna_scores <- mrna_pred$predict[, "Mutant", 2]

# Single-omics Methylation
load("Single_Omics_Methylation_Results/Model_Workspace.RData")
Y_test_meth <- Y_test
meth_pred <- predict(final_model, newdata = X_test)
meth_scores <- meth_pred$predict[, "Mutant", 2]

# Random Forest
load("Random_Forest_Concatenated_Results/Model_Workspace.RData")
Y_test_rf <- Y_test
rf_scores <- predict(rf_model, newdata = X_test, type = "prob")[, "Mutant"]

# Verify all have same test samples
if (!all(names(Y_test_diablo) == names(Y_test_mrna)) ||
    !all(names(Y_test_diablo) == names(Y_test_meth)) ||
    !all(names(Y_test_diablo) == names(Y_test_rf))) {
  stop("Test sets are not identical across methods!")
}

cat(sprintf("✓ Loaded predictions for %d test samples\n", length(Y_test_diablo)))
cat("✓ All methods have identical test sets\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# 2. BOOTSTRAP RESAMPLING
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 2: Bootstrap Resampling (B=1000)\n")
cat("─────────────────────────────────────────────────────────────────\n")

n_bootstrap <- 1000
n_samples <- length(Y_test_diablo)

# Storage for bootstrap AUCs
boot_auc_diablo <- numeric(n_bootstrap)
boot_auc_mrna <- numeric(n_bootstrap)
boot_auc_meth <- numeric(n_bootstrap)
boot_auc_rf <- numeric(n_bootstrap)

cat("Running bootstrap iterations...\n")
pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3, width = 50)

for (i in 1:n_bootstrap) {
  # Resample indices with replacement
  idx <- sample(1:n_samples, n_samples, replace = TRUE)
  
  # Calculate AUC for each method on bootstrap sample
  boot_auc_diablo[i] <- as.numeric(auc(roc(Y_test_diablo[idx], diablo_scores[idx], 
                                          levels = c("WT", "Mutant"), direction = "<", quiet = TRUE)))
  boot_auc_mrna[i] <- as.numeric(auc(roc(Y_test_mrna[idx], mrna_scores[idx], 
                                        levels = c("WT", "Mutant"), direction = "<", quiet = TRUE)))
  boot_auc_meth[i] <- as.numeric(auc(roc(Y_test_meth[idx], meth_scores[idx], 
                                        levels = c("WT", "Mutant"), direction = "<", quiet = TRUE)))
  boot_auc_rf[i] <- as.numeric(auc(roc(Y_test_rf[idx], rf_scores[idx], 
                                      levels = c("WT", "Mutant"), direction = "<", quiet = TRUE)))
  
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n\n✓ Bootstrap complete\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# 3. DELONG TEST - ROC CURVE COMPARISON (GOLD STANDARD)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 3: DeLong Test for ROC Curve Comparison\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("  DeLong test: Gold standard for comparing correlated ROC curves\n")
cat("  (Same test set, different models)\n\n")

# Create ROC objects for each method
roc_diablo <- roc(Y_test_diablo, diablo_scores, levels = c("WT", "Mutant"), 
                  direction = "<", quiet = TRUE)
roc_mrna <- roc(Y_test_mrna, mrna_scores, levels = c("WT", "Mutant"), 
                direction = "<", quiet = TRUE)
roc_meth <- roc(Y_test_meth, meth_scores, levels = c("WT", "Mutant"), 
                direction = "<", quiet = TRUE)
roc_rf <- roc(Y_test_rf, rf_scores, levels = c("WT", "Mutant"), 
              direction = "<", quiet = TRUE)

# Function to perform DeLong test
delong_compare <- function(roc1, roc2, name1, name2) {
  # DeLong test (default method in roc.test for paired/correlated ROCs)
  test_result <- roc.test(roc1, roc2, method = "delong", paired = TRUE)
  
  data.frame(
    Comparison = paste(name1, "vs", name2),
    Method_1 = name1,
    Method_2 = name2,
    AUC_1 = as.numeric(auc(roc1)),
    AUC_2 = as.numeric(auc(roc2)),
    AUC_Difference = as.numeric(auc(roc1)) - as.numeric(auc(roc2)),
    DeLong_Z = test_result$statistic,
    DeLong_p_value = test_result$p.value,
    Significant = ifelse(test_result$p.value < 0.05, "Yes", "No"),
    CI_95_lower = NA,  # Will be filled from bootstrap
    CI_95_upper = NA,
    stringsAsFactors = FALSE
  )
}

# Perform all pairwise DeLong tests
delong_results <- rbind(
  delong_compare(roc_diablo, roc_mrna, "DIABLO", "mRNA"),
  delong_compare(roc_diablo, roc_meth, "DIABLO", "Methylation"),
  delong_compare(roc_diablo, roc_rf, "DIABLO", "Random_Forest"),
  delong_compare(roc_mrna, roc_meth, "mRNA", "Methylation"),
  delong_compare(roc_mrna, roc_rf, "mRNA", "Random_Forest"),
  delong_compare(roc_meth, roc_rf, "Methylation", "Random_Forest")
)

cat("✓ DeLong tests complete\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# 4. PAIRWISE BOOTSTRAP COMPARISONS (FOR 95% CI)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 4: Bootstrap Comparisons for 95% Confidence Intervals\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Function to compare two methods via bootstrap
compare_methods_bootstrap <- function(auc1, auc2, name1, name2) {
  # Wilcoxon paired test (bootstrap distributions are paired)
  test_result <- wilcox.test(auc1, auc2, paired = TRUE, alternative = "two.sided")
  
  # Calculate 95% CI of the difference
  diff <- auc1 - auc2
  ci_lower <- quantile(diff, 0.025)
  ci_upper <- quantile(diff, 0.975)
  
  data.frame(
    Comparison = paste(name1, "vs", name2),
    Method_1 = name1,
    Method_2 = name2,
    Mean_AUC_1 = mean(auc1),
    Mean_AUC_2 = mean(auc2),
    Median_AUC_1 = median(auc1),
    Median_AUC_2 = median(auc2),
    Mean_Difference = mean(auc1) - mean(auc2),
    CI_95_lower = ci_lower,
    CI_95_upper = ci_upper,
    Bootstrap_p_value = test_result$p.value,
    Significant = ifelse(test_result$p.value < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

# Perform all pairwise bootstrap comparisons
bootstrap_comparisons <- rbind(
  compare_methods_bootstrap(boot_auc_diablo, boot_auc_mrna, "DIABLO", "mRNA"),
  compare_methods_bootstrap(boot_auc_diablo, boot_auc_meth, "DIABLO", "Methylation"),
  compare_methods_bootstrap(boot_auc_diablo, boot_auc_rf, "DIABLO", "Random_Forest"),
  compare_methods_bootstrap(boot_auc_mrna, boot_auc_meth, "mRNA", "Methylation"),
  compare_methods_bootstrap(boot_auc_mrna, boot_auc_rf, "mRNA", "Random_Forest"),
  compare_methods_bootstrap(boot_auc_meth, boot_auc_rf, "Methylation", "Random_Forest")
)

# Add 95% CI to DeLong results from bootstrap
for (i in seq_len(nrow(delong_results))) {
  delong_results$CI_95_lower[i] <- bootstrap_comparisons$CI_95_lower[i]
  delong_results$CI_95_upper[i] <- bootstrap_comparisons$CI_95_upper[i]
}

cat("✓ Bootstrap 95% CI calculated\n\n")

# Display results
cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║             METHOD COMPARISON RESULTS (DeLong Test)           ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

for (i in seq_len(nrow(delong_results))) {
  cat(sprintf("%s:\n", delong_results$Comparison[i]))
  cat(sprintf("  • AUC: %.4f vs %.4f (Δ = %.4f)\n", 
              delong_results$AUC_1[i], delong_results$AUC_2[i], 
              delong_results$AUC_Difference[i]))
  cat(sprintf("  • DeLong Z-statistic: %.4f\n", delong_results$DeLong_Z[i]))
  cat(sprintf("  • DeLong p-value: %.4e\n", delong_results$DeLong_p_value[i]))
  cat(sprintf("  • 95%% CI (Bootstrap): [%.4f, %.4f]\n", 
              delong_results$CI_95_lower[i], delong_results$CI_95_upper[i]))
  cat(sprintf("  %s\n", 
              ifelse(delong_results$Significant[i] == "Yes", 
                     "★ SIGNIFICANT (p < 0.05)", 
                     "Not significant (p ≥ 0.05)")))
  cat("\n")
}

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║         BOOTSTRAP DISTRIBUTION SUMMARY (B=1000)                ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

for (i in seq_len(nrow(bootstrap_comparisons))) {
  cat(sprintf("%s:\n", bootstrap_comparisons$Comparison[i]))
  cat(sprintf("  • Mean AUC: %.4f vs %.4f\n", 
              bootstrap_comparisons$Mean_AUC_1[i], bootstrap_comparisons$Mean_AUC_2[i]))
  cat(sprintf("  • Mean Difference: %.4f [95%% CI: %.4f, %.4f]\n", 
              bootstrap_comparisons$Mean_Difference[i],
              bootstrap_comparisons$CI_95_lower[i],
              bootstrap_comparisons$CI_95_upper[i]))
  cat(sprintf("  • Bootstrap p-value (Wilcoxon): %.4e\n", 
              bootstrap_comparisons$Bootstrap_p_value[i]))
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# 5. SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 5: Saving Results\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Save DeLong test results (PRIMARY - Gold Standard for ROC comparison)
write.csv(delong_results, "Method_Comparison_DeLong_Tests.csv", row.names = FALSE)

# Save bootstrap comparison table (SUPPLEMENTARY - for 95% CI)
write.csv(bootstrap_comparisons, "Method_Comparison_Bootstrap_Pvalues.csv", row.names = FALSE)

# Save bootstrap distributions (RAW DATA)
boot_data <- data.frame(
  DIABLO = boot_auc_diablo,
  mRNA = boot_auc_mrna,
  Methylation = boot_auc_meth,
  Random_Forest = boot_auc_rf
)
write.csv(boot_data, "Method_Comparison_Bootstrap_Distributions.csv", row.names = FALSE)

# Create boxplot with DeLong results annotated
pdf("Method_Comparison_Boxplot.pdf", width = 12, height = 8)
par(mar = c(5, 5, 5, 2))

boxplot(boot_auc_diablo, boot_auc_mrna, boot_auc_meth, boot_auc_rf,
        names = c("DIABLO", "mRNA", "Methylation", "RF"),
        col = c("#2E86AB", "#A23B72", "#F18F01", "#06A77D"),
        main = "Method Comparison: Bootstrap AUC Distributions (B=1000)\nStatistical Significance via DeLong Test",
        ylab = "AUC (Bootstrap Resampling)",
        cex.main = 1.3, cex.lab = 1.3, cex.axis = 1.1,
        ylim = c(0.85, 1.0))

# Add horizontal line at 0.9
abline(h = 0.9, lty = 2, col = "gray50", lwd = 2)

# Add median AUC values as text
medians <- c(median(boot_auc_diablo), median(boot_auc_mrna), 
             median(boot_auc_meth), median(boot_auc_rf))
text(1:4, medians + 0.01, sprintf("%.3f", medians), 
     col = "black", cex = 1.1, font = 2)

dev.off()

cat("✓ Method_Comparison_DeLong_Tests.csv (PRIMARY - Gold Standard)\n")
cat("✓ Method_Comparison_Bootstrap_Pvalues.csv (SUPPLEMENTARY - 95% CI)\n")
cat("✓ Method_Comparison_Bootstrap_Distributions.csv (RAW DATA)\n")
cat("✓ Method_Comparison_Boxplot.pdf\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# 6. FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                    ANALYSIS COMPLETE                           ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

n_sig_delong <- sum(delong_results$Significant == "Yes")
n_sig_bootstrap <- sum(bootstrap_comparisons$Significant == "Yes")

cat(sprintf("STATISTICAL COMPARISON SUMMARY:\n"))
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  • Total comparisons: 6 (all pairwise)\n"))
cat(sprintf("  • Test set size: %d samples (identical for all methods)\n", n_samples))
cat(sprintf("  • Bootstrap iterations: %d\n", n_bootstrap))
cat(sprintf("  • Seed: 42 (fully reproducible)\n\n"))

cat("DELONG TEST RESULTS (PRIMARY - Gold Standard for ROC):\n")
cat(sprintf("  • Significant differences (p < 0.05): %d / 6\n", n_sig_delong))
if (n_sig_delong > 0) {
  cat("  ★ Significant comparisons:\n")
  sig_delong <- delong_results[delong_results$Significant == "Yes", ]
  for (i in seq_len(nrow(sig_delong))) {
    cat(sprintf("    • %s (DeLong p = %.4e, Z = %.3f)\n", 
                sig_delong$Comparison[i], 
                sig_delong$DeLong_p_value[i],
                sig_delong$DeLong_Z[i]))
  }
}
cat("\n")

cat("BOOTSTRAP ANALYSIS (SUPPLEMENTARY - 95% Confidence Intervals):\n")
cat(sprintf("  • Significant differences (Wilcoxon p < 0.05): %d / 6\n", n_sig_bootstrap))
cat("  • Bootstrap distributions provide 95% CI for AUC differences\n\n")

if (n_sig_delong == 0) {
  cat("INTERPRETATION:\n")
  cat("  ✓ No significant differences detected (DeLong test)\n")
  cat("  → All methods show statistically comparable performance\n")
  cat("  → Multi-omics integration does not significantly outperform single-omics\n\n")
} else {
  cat("INTERPRETATION:\n")
  cat("  ★ Statistically significant performance differences detected\n")
  cat("  → Some methods perform significantly better than others\n")
  cat("  → See DeLong test results for specific comparisons\n\n")
}

cat("FOR MANUSCRIPT:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("  PRIMARY ANALYSIS:\n")
cat("    • Use 'Method_Comparison_DeLong_Tests.csv'\n")
cat("    • Report: \"DeLong test for correlated ROC curves\"\n")
cat("    • Citation: DeLong et al. (1988) Biometrics\n\n")
cat("  SUPPLEMENTARY ANALYSIS:\n")
cat("    • Use 'Method_Comparison_Bootstrap_Pvalues.csv'\n")
cat("    • Report: \"Bootstrap resampling (B=1000) for 95% CI\"\n")
cat("    • Include 'Method_Comparison_Boxplot.pdf'\n\n")

cat("METHODOLOGY STATEMENT (for Methods section):\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("\"Statistical comparison of model performance was performed using\n")
cat("DeLong's test for correlated ROC curves (DeLong et al., 1988),\n")
cat("implemented in the pROC R package. Bootstrap resampling (B=1000)\n")
cat("was used to estimate 95% confidence intervals for AUC differences.\n")
cat("All methods were evaluated on an identical test set (n=%d) using\n", n_samples)
cat("seed=42 for full reproducibility. P-values < 0.05 were considered\n")
cat("statistically significant.\"\n\n")

cat("✓ READY FOR PUBLICATION\n\n")
