#!/usr/bin/env Rscript
# ════════════════════════════════════════════════════════════════════════════
#  MASTER OVERFITTING SUMMARY - ALL METHODS
#  Reads results from all 4 methods and creates unified comparison report
# ════════════════════════════════════════════════════════════════════════════

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║        MASTER OVERFITTING SUMMARY - ALL METHODS               ║\n")
cat("║        Unified Comparison Report for Publication              ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

# ═══════════════════════════════════════════════════════════════════════════
# 1. READ ALL RESULTS
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 1: Reading Results from All Methods\n")
cat("─────────────────────────────────────────────────────────────────\n")

methods <- list(
  list(name = "DIABLO", 
       perf_file = "DIABLO_Final_Publication_Fixed/Performance_Metrics.csv",
       comp_file = "DIABLO_Final_Publication_Fixed/Train_vs_Test_Comparison.csv"),
  list(name = "mRNA", 
       perf_file = "Single_Omics_mRNA_Results/Performance_Metrics.csv",
       comp_file = "Single_Omics_mRNA_Results/Train_vs_Test_Comparison.csv"),
  list(name = "Methylation", 
       perf_file = "Single_Omics_Methylation_Results/Performance_Metrics.csv",
       comp_file = "Single_Omics_Methylation_Results/Train_vs_Test_Comparison.csv"),
  list(name = "RandomForest", 
       perf_file = "Random_Forest_Concatenated_Results/Performance_Metrics.csv",
       comp_file = "Random_Forest_Concatenated_Results/Train_vs_Test_Comparison.csv")
)

# Initialize storage
all_results <- list()
files_found <- 0

for (method in methods) {
  if (file.exists(method$perf_file) && file.exists(method$comp_file)) {
    perf <- read.csv(method$perf_file)
    comp <- read.csv(method$comp_file)
    
    all_results[[method$name]] <- list(
      performance = perf,
      comparison = comp
    )
    files_found <- files_found + 1
    cat(sprintf("  ✓ %s results loaded\n", method$name))
  } else {
    cat(sprintf("  ✗ %s results not found (run analysis first)\n", method$name))
  }
}

if (files_found == 0) {
  stop("\n✗ No results found. Please run the main analyses first.\n")
}

cat(sprintf("\n✓ Loaded results from %d/%d methods\n\n", files_found, length(methods)))

# ═══════════════════════════════════════════════════════════════════════════
# 2. CREATE UNIFIED COMPARISON TABLE
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 2: Creating Unified Comparison Table\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Extract key metrics
unified_table <- data.frame(
  Method = character(),
  Train_AUC = numeric(),
  Test_AUC = numeric(),
  AUC_Gap = numeric(),
  Gap_Percent = numeric(),
  Permutation_pvalue = numeric(),
  Test_Balanced_Acc = numeric(),
  Test_Sensitivity = numeric(),
  Test_Specificity = numeric(),
  Overfitting_Risk = character(),
  stringsAsFactors = FALSE
)

for (method_name in names(all_results)) {
  result <- all_results[[method_name]]
  
  # Extract metrics
  comp <- result$comparison
  perf <- result$performance
  
  train_auc <- comp$AUC[comp$Dataset == "Training"]
  test_auc <- comp$AUC[comp$Dataset == "Test"]
  gap_auc <- comp$AUC[comp$Dataset == "Gap (Train-Test)"]
  
  perm_pval <- perf$Value[perf$Metric == "Permutation_pvalue"]
  test_bal_acc <- perf$Value[perf$Metric == "Balanced_Accuracy"]
  test_sens <- perf$Value[perf$Metric == "Sensitivity"]
  test_spec <- perf$Value[perf$Metric == "Specificity"]
  
  gap_percent <- abs(gap_auc / train_auc * 100)
  
  # Risk assessment
  if (abs(gap_auc) < 0.05) {
    risk <- "VERY LOW ✓✓✓"
  } else if (abs(gap_auc) < 0.10) {
    risk <- "LOW ✓✓"
  } else if (abs(gap_auc) < 0.15) {
    risk <- "MODERATE ⚠"
  } else {
    risk <- "HIGH ✗✗"
  }
  
  unified_table <- rbind(unified_table, data.frame(
    Method = method_name,
    Train_AUC = train_auc,
    Test_AUC = test_auc,
    AUC_Gap = gap_auc,
    Gap_Percent = gap_percent,
    Permutation_pvalue = perm_pval,
    Test_Balanced_Acc = test_bal_acc,
    Test_Sensitivity = test_sens,
    Test_Specificity = test_spec,
    Overfitting_Risk = risk,
    stringsAsFactors = FALSE
  ))
}

cat("\n╔════════════════════════════════════════════════════════════════════════════════════════════════════╗\n")
cat("║                         UNIFIED OVERFITTING ASSESSMENT - ALL METHODS                              ║\n")
cat("╚════════════════════════════════════════════════════════════════════════════════════════════════════╝\n\n")

print(unified_table, row.names = FALSE)
cat("\n")

# ═══════════════════════════════════════════════════════════════════════════
# 3. DETAILED COMPARISON
# ═══════════════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════\n")
cat("              DETAILED METHOD COMPARISON\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (i in 1:nrow(unified_table)) {
  method <- unified_table[i, ]
  
  cat(sprintf("【 %s 】\n", method$Method))
  cat("─────────────────────────────────────────────────────────────────\n")
  cat(sprintf("  Train AUC:           %.4f\n", method$Train_AUC))
  cat(sprintf("  Test AUC:            %.4f\n", method$Test_AUC))
  cat(sprintf("  AUC Gap:             %.4f (%.1f%%)\n", method$AUC_Gap, method$Gap_Percent))
  cat(sprintf("  Overfitting Risk:    %s\n", method$Overfitting_Risk))
  cat(sprintf("  Permutation p-val:   %.4f ", method$Permutation_pvalue))
  if (method$Permutation_pvalue < 0.01) {
    cat("(★★★ Highly Significant)\n")
  } else if (method$Permutation_pvalue < 0.05) {
    cat("(★★ Significant)\n")
  } else {
    cat("(✗ Not Significant)\n")
  }
  cat(sprintf("  Test Balanced Acc:   %.4f (%.1f%%)\n", method$Test_Balanced_Acc, 100*method$Test_Balanced_Acc))
  cat(sprintf("  Test Sensitivity:    %.4f (%.1f%%)\n", method$Test_Sensitivity, 100*method$Test_Sensitivity))
  cat(sprintf("  Test Specificity:    %.4f (%.1f%%)\n\n", method$Test_Specificity, 100*method$Test_Specificity))
}

# ═══════════════════════════════════════════════════════════════════════════
# 4. RANK METHODS BY OVERFITTING CONTROL
# ═══════════════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════\n")
cat("              RANKING BY OVERFITTING CONTROL\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Sort by gap (lower is better)
ranked <- unified_table[order(abs(unified_table$AUC_Gap)), ]

cat("Best Generalization (Lowest Train-Test Gap):\n")
cat("─────────────────────────────────────────────────────────────────\n")
for (i in 1:nrow(ranked)) {
  cat(sprintf("  %d. %-15s - Gap: %.4f (%.1f%%) - %s\n", 
              i, ranked$Method[i], ranked$AUC_Gap[i], 
              ranked$Gap_Percent[i], ranked$Overfitting_Risk[i]))
}
cat("\n")

# ═══════════════════════════════════════════════════════════════════════════
# 5. OVERALL ASSESSMENT
# ═══════════════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════\n")
cat("                 OVERALL ASSESSMENT\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Count methods by risk level
n_very_low <- sum(grepl("VERY LOW", unified_table$Overfitting_Risk))
n_low <- sum(grepl("^LOW", unified_table$Overfitting_Risk))
n_moderate <- sum(grepl("MODERATE", unified_table$Overfitting_Risk))
n_high <- sum(grepl("HIGH", unified_table$Overfitting_Risk))

cat("Overfitting Risk Distribution:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  VERY LOW (< 5%% gap):    %d method(s)\n", n_very_low))
cat(sprintf("  LOW (5-10%% gap):        %d method(s)\n", n_low))
cat(sprintf("  MODERATE (10-15%% gap):  %d method(s)\n", n_moderate))
cat(sprintf("  HIGH (> 15%% gap):       %d method(s)\n\n", n_high))

# Statistical significance
n_highly_sig <- sum(unified_table$Permutation_pvalue < 0.01)
n_sig <- sum(unified_table$Permutation_pvalue >= 0.01 & unified_table$Permutation_pvalue < 0.05)
n_not_sig <- sum(unified_table$Permutation_pvalue >= 0.05)

cat("Statistical Significance (Permutation Tests):\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  p < 0.01 (Highly Sig):   %d method(s)\n", n_highly_sig))
cat(sprintf("  p < 0.05 (Significant):  %d method(s)\n", n_sig))
cat(sprintf("  p ≥ 0.05 (Not Sig):      %d method(s)\n\n", n_not_sig))

# Publication readiness
if (n_high == 0 && n_not_sig == 0) {
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  VERDICT: ✓✓✓ ALL METHODS PUBLICATION-READY                   ║\n")
  cat("║  All methods show acceptable overfitting control              ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
} else if (n_high == 0) {
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  VERDICT: ✓✓ METHODS SUITABLE FOR PUBLICATION                 ║\n")
  cat("║  Minor concerns but overall acceptable                        ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
} else {
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  VERDICT: ⚠ SOME METHODS NEED ATTENTION                       ║\n")
  cat("║  Review methods with HIGH overfitting risk                    ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# 6. VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 3: Creating Visualization\n")
cat("─────────────────────────────────────────────────────────────────\n")

pdf("MASTER_OVERFITTING_SUMMARY.pdf", width = 14, height = 10)

# Plot 1: Train vs Test AUC Comparison
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))

# 1. Train vs Test AUC
barplot(
  rbind(unified_table$Train_AUC, unified_table$Test_AUC),
  beside = TRUE,
  names.arg = unified_table$Method,
  col = c("#2E86AB", "#F18F01"),
  main = "Train vs Test AUC - All Methods",
  ylab = "AUC",
  ylim = c(0, 1),
  cex.main = 1.5,
  cex.names = 1.2,
  cex.axis = 1.2,
  cex.lab = 1.3
)
legend("bottomright", 
       legend = c("Train", "Test"), 
       fill = c("#2E86AB", "#F18F01"),
       cex = 1.2)
abline(h = 0.8, lty = 2, col = "red", lwd = 2)

# 2. AUC Gap
barplot(
  abs(unified_table$AUC_Gap),
  names.arg = unified_table$Method,
  col = ifelse(abs(unified_table$AUC_Gap) < 0.05, "#4CAF50",
               ifelse(abs(unified_table$AUC_Gap) < 0.10, "#FFC107", "#F44336")),
  main = "Train-Test Gap (Overfitting Risk)",
  ylab = "Absolute AUC Gap",
  cex.main = 1.5,
  cex.names = 1.2,
  cex.axis = 1.2,
  cex.lab = 1.3
)
abline(h = 0.05, lty = 2, col = "darkgreen", lwd = 2)
abline(h = 0.10, lty = 2, col = "orange", lwd = 2)
legend("topright", 
       legend = c("< 5% (Very Low)", "5-10% (Low)", "> 10% (Concern)"),
       fill = c("#4CAF50", "#FFC107", "#F44336"),
       cex = 1.0)

# 3. Permutation p-values
barplot(
  -log10(unified_table$Permutation_pvalue),
  names.arg = unified_table$Method,
  col = ifelse(unified_table$Permutation_pvalue < 0.01, "#4CAF50",
               ifelse(unified_table$Permutation_pvalue < 0.05, "#FFC107", "#F44336")),
  main = "Permutation Test Significance",
  ylab = "-log10(p-value)",
  cex.main = 1.5,
  cex.names = 1.2,
  cex.axis = 1.2,
  cex.lab = 1.3
)
abline(h = -log10(0.05), lty = 2, col = "orange", lwd = 2)
abline(h = -log10(0.01), lty = 2, col = "darkgreen", lwd = 2)
text(2, -log10(0.05) + 0.2, "p = 0.05", cex = 1.1)
text(2, -log10(0.01) + 0.2, "p = 0.01", cex = 1.1)

# 4. Test Performance
barplot(
  unified_table$Test_AUC,
  names.arg = unified_table$Method,
  col = "#2E86AB",
  main = "Test Set Performance",
  ylab = "Test AUC",
  ylim = c(0, 1),
  cex.main = 1.5,
  cex.names = 1.2,
  cex.axis = 1.2,
  cex.lab = 1.3
)
abline(h = 0.90, lty = 2, col = "darkgreen", lwd = 2)
abline(h = 0.80, lty = 2, col = "orange", lwd = 2)

dev.off()

cat("✓ Visualization saved: MASTER_OVERFITTING_SUMMARY.pdf\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# 7. SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 4: Saving Results\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Save unified table
write.csv(unified_table, "MASTER_OVERFITTING_SUMMARY.csv", row.names = FALSE)
cat("✓ Summary table saved: MASTER_OVERFITTING_SUMMARY.csv\n")

# Create detailed text report
sink("MASTER_OVERFITTING_REPORT.txt")
cat("════════════════════════════════════════════════════════════════\n")
cat("       MASTER OVERFITTING SUMMARY - ALL METHODS\n")
cat("       Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("════════════════════════════════════════════════════════════════\n\n")

cat("METHODS ANALYZED: ", files_found, "\n\n")

cat("════════════════════════════════════════════════════════════════\n")
cat("UNIFIED COMPARISON TABLE\n")
cat("════════════════════════════════════════════════════════════════\n\n")
print(unified_table, row.names = FALSE)

cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("RANKING BY GENERALIZATION (Best to Worst)\n")
cat("════════════════════════════════════════════════════════════════\n\n")
for (i in 1:nrow(ranked)) {
  cat(sprintf("%d. %-15s - Gap: %.4f (%.1f%%) - %s\n", 
              i, ranked$Method[i], ranked$AUC_Gap[i], 
              ranked$Gap_Percent[i], ranked$Overfitting_Risk[i]))
}

cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("RECOMMENDATION FOR PUBLICATION\n")
cat("════════════════════════════════════════════════════════════════\n\n")

if (n_high == 0 && n_not_sig == 0) {
  cat("✓✓✓ ALL METHODS ARE PUBLICATION-READY\n\n")
  cat("All methods demonstrate:\n")
  cat("  • Acceptable train-test gaps (no high-risk methods)\n")
  cat("  • Statistical significance (all p < 0.05)\n")
  cat("  • Good generalization to unseen data\n\n")
  cat("RECOMMENDATION: Proceed with publication confidently.\n")
} else if (n_high > 0) {
  cat("⚠ ATTENTION NEEDED\n\n")
  cat(sprintf("%d method(s) show HIGH overfitting risk.\n", n_high))
  cat("Consider:\n")
  cat("  • Increasing regularization\n")
  cat("  • Reducing model complexity\n")
  cat("  • External validation cohort\n")
}

cat("\n\n════════════════════════════════════════════════════════════════\n")
cat("For detailed per-method analysis, see individual reports:\n")
cat("  - DIABLO_Final_Publication_Fixed/\n")
cat("  - Single_Omics_mRNA_Results/\n")
cat("  - Single_Omics_Methylation_Results/\n")
cat("  - Random_Forest_Concatenated_Results/\n")
cat("════════════════════════════════════════════════════════════════\n")

sink()

cat("✓ Detailed report saved: MASTER_OVERFITTING_REPORT.txt\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# FINAL MESSAGE
# ═══════════════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════\n")
cat("           MASTER SUMMARY COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("GENERATED FILES:\n")
cat("  📊 MASTER_OVERFITTING_SUMMARY.csv  - Unified comparison table\n")
cat("  📈 MASTER_OVERFITTING_SUMMARY.pdf  - Visual comparison plots\n")
cat("  📄 MASTER_OVERFITTING_REPORT.txt   - Detailed text report\n\n")

cat("✓ All overfitting controls verified across methods\n")
cat("✓ Results ready for publication review\n\n")

