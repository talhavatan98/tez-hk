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
# 3. PAIRWISE COMPARISONS
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 3: Pairwise Statistical Comparisons\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Function to compare two methods
compare_methods <- function(auc1, auc2, name1, name2) {
  # Wilcoxon paired test (bootstrap distributions are paired)
  test_result <- wilcox.test(auc1, auc2, paired = TRUE, alternative = "two.sided")
  
  data.frame(
    Comparison = paste(name1, "vs", name2),
    Method_1 = name1,
    Method_2 = name2,
    Mean_AUC_1 = mean(auc1),
    Mean_AUC_2 = mean(auc2),
    Median_AUC_1 = median(auc1),
    Median_AUC_2 = median(auc2),
    Difference = mean(auc1) - mean(auc2),
    p_value = test_result$p.value,
    Significant = ifelse(test_result$p.value < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

# Perform all pairwise comparisons
comparisons <- rbind(
  compare_methods(boot_auc_diablo, boot_auc_mrna, "DIABLO", "mRNA"),
  compare_methods(boot_auc_diablo, boot_auc_meth, "DIABLO", "Methylation"),
  compare_methods(boot_auc_diablo, boot_auc_rf, "DIABLO", "Random_Forest"),
  compare_methods(boot_auc_mrna, boot_auc_meth, "mRNA", "Methylation"),
  compare_methods(boot_auc_mrna, boot_auc_rf, "mRNA", "Random_Forest"),
  compare_methods(boot_auc_meth, boot_auc_rf, "Methylation", "Random_Forest")
)

# Display results
cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║               METHOD COMPARISON RESULTS                        ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

for (i in 1:nrow(comparisons)) {
  cat(sprintf("%s:\n", comparisons$Comparison[i]))
  cat(sprintf("  • Mean AUC: %.4f vs %.4f (Δ = %.4f)\n", 
              comparisons$Mean_AUC_1[i], comparisons$Mean_AUC_2[i], 
              comparisons$Difference[i]))
  cat(sprintf("  • P-value: %.4f %s\n", 
              comparisons$p_value[i], 
              ifelse(comparisons$Significant[i] == "Yes", "★ SIGNIFICANT", "")))
  cat("\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# 4. SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 4: Saving Results\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Save comparison table
write.csv(comparisons, "Method_Comparison_Pvalues.csv", row.names = FALSE)

# Save bootstrap distributions
boot_data <- data.frame(
  DIABLO = boot_auc_diablo,
  mRNA = boot_auc_mrna,
  Methylation = boot_auc_meth,
  Random_Forest = boot_auc_rf
)
write.csv(boot_data, "Method_Comparison_Bootstrap_Distributions.csv", row.names = FALSE)

# Create boxplot
pdf("Method_Comparison_Boxplot.pdf", width = 10, height = 7)
par(mar = c(5, 5, 4, 2))

boxplot(boot_auc_diablo, boot_auc_mrna, boot_auc_meth, boot_auc_rf,
        names = c("DIABLO", "mRNA", "Methylation", "RF"),
        col = c("#2E86AB", "#A23B72", "#F18F01", "#06A77D"),
        main = "Method Comparison (Bootstrap AUC Distributions)",
        ylab = "AUC (Bootstrap Resampling, B=1000)",
        cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.1,
        ylim = c(0.7, 1.0))

# Add horizontal line at 0.9
abline(h = 0.9, lty = 2, col = "gray50", lwd = 2)

dev.off()

cat("✓ Method_Comparison_Pvalues.csv\n")
cat("✓ Method_Comparison_Bootstrap_Distributions.csv\n")
cat("✓ Method_Comparison_Boxplot.pdf\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# 5. SUMMARY
# ═══════════════════════════════════════════════════════════════════════════

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                    ANALYSIS COMPLETE                           ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

n_sig <- sum(comparisons$Significant == "Yes")
cat(sprintf("KEY FINDINGS:\n"))
cat(sprintf("  • %d out of 6 comparisons show significant differences (p < 0.05)\n", n_sig))
cat(sprintf("  • Bootstrap iterations: %d\n", n_bootstrap))
cat(sprintf("  • Test samples: %d\n\n", n_samples))

if (n_sig == 0) {
  cat("INTERPRETATION:\n")
  cat("  ✓ No significant differences detected between methods\n")
  cat("  → All approaches show comparable performance\n")
  cat("  → Multi-omics integration does not significantly outperform single-omics\n\n")
} else {
  cat("INTERPRETATION:\n")
  cat("  ★ Significant differences detected:\n")
  sig_comps <- comparisons[comparisons$Significant == "Yes", ]
  for (i in 1:nrow(sig_comps)) {
    cat(sprintf("    • %s (p = %.4f)\n", sig_comps$Comparison[i], sig_comps$p_value[i]))
  }
  cat("\n")
}

cat("FOR MANUSCRIPT:\n")
cat("  Use 'Method_Comparison_Pvalues.csv' for Results/Methods sections\n")
cat("  Include 'Method_Comparison_Boxplot.pdf' as Supplementary Figure\n\n")
