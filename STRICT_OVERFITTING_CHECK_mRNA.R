#!/usr/bin/env Rscript
# ════════════════════════════════════════════════════════════════════════════
#  STRICT OVERFITTING CHECK - Single Omics mRNA
#  Comprehensive analysis: Train vs Test performance comparison
# ════════════════════════════════════════════════════════════════════════════

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║     STRICT OVERFITTING CONTROL - Single Omics mRNA            ║\n")
cat("║  Training vs Test Performance Gap Analysis                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

suppressPackageStartupMessages({
  library(mixOmics)
  library(pROC)
  library(caret)
})

# Load saved workspace
load("Single_Omics_mRNA_Results/Model_Workspace.RData")

cat("═══════════════════════════════════════════════════════════════\n")
cat("              PART 1: TRAINING SET PERFORMANCE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Get training set predictions
max_ncomp <- final_model$ncomp[1]

suppressWarnings({
  train_predictions <- predict(
    object = final_model, 
    newdata = X_train
  )
})

# Extract predictions
train_pred_classes <- train_predictions$class$centroids.dist[, max_ncomp]
train_pred_scores <- train_predictions$predict[, "Mutant", max_ncomp]

# Training confusion matrix
train_cm <- confusionMatrix(
  factor(train_pred_classes, levels = c("Mutant", "WT")),
  factor(Y_train, levels = c("Mutant", "WT")),
  positive = "Mutant"
)

# Training metrics
train_bal_acc <- train_cm$byClass["Balanced Accuracy"]
train_sens <- train_cm$byClass["Sensitivity"]
train_spec <- train_cm$byClass["Specificity"]
train_f1 <- train_cm$byClass["F1"]

# Training AUC
train_roc <- roc(Y_train, train_pred_scores, levels = c("WT", "Mutant"), 
                  direction = "<", quiet = TRUE)
train_auc <- as.numeric(auc(train_roc))

cat("TRAINING SET PERFORMANCE:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Samples:            %d\n", length(Y_train)))
cat(sprintf("  Balanced Accuracy:  %.4f (%.1f%%)\n", train_bal_acc, 100*train_bal_acc))
cat(sprintf("  Sensitivity:        %.4f (%.1f%%)\n", train_sens, 100*train_sens))
cat(sprintf("  Specificity:        %.4f (%.1f%%)\n", train_spec, 100*train_spec))
cat(sprintf("  F1-Score:           %.4f\n", train_f1))
cat(sprintf("  ROC-AUC:            %.4f\n\n", train_auc))

cat("Confusion Matrix (Train):\n")
print(train_cm$table)
cat("\n\n")

# ═══════════════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════\n")
cat("              PART 2: TEST SET PERFORMANCE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Get test set predictions
suppressWarnings({
  test_predictions <- predict(
    object = final_model, 
    newdata = X_test
  )
})

test_pred_classes <- test_predictions$class$centroids.dist[, max_ncomp]
test_pred_scores <- test_predictions$predict[, "Mutant", max_ncomp]

# Test confusion matrix
test_cm <- confusionMatrix(
  factor(test_pred_classes, levels = c("Mutant", "WT")),
  factor(Y_test, levels = c("Mutant", "WT")),
  positive = "Mutant"
)

# Test metrics
test_bal_acc <- test_cm$byClass["Balanced Accuracy"]
test_sens <- test_cm$byClass["Sensitivity"]
test_spec <- test_cm$byClass["Specificity"]
test_f1 <- test_cm$byClass["F1"]

# Test AUC
test_roc <- roc(Y_test, test_pred_scores, levels = c("WT", "Mutant"), 
                 direction = "<", quiet = TRUE)
test_auc <- as.numeric(auc(test_roc))

cat("TEST SET PERFORMANCE:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  Samples:            %d\n", length(Y_test)))
cat(sprintf("  Balanced Accuracy:  %.4f (%.1f%%)\n", test_bal_acc, 100*test_bal_acc))
cat(sprintf("  Sensitivity:        %.4f (%.1f%%)\n", test_sens, 100*test_sens))
cat(sprintf("  Specificity:        %.4f (%.1f%%)\n", test_spec, 100*test_spec))
cat(sprintf("  F1-Score:           %.4f\n", test_f1))
cat(sprintf("  ROC-AUC:            %.4f\n\n", test_auc))

cat("Confusion Matrix (Test):\n")
print(test_cm$table)
cat("\n\n")

# ═══════════════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════\n")
cat("           PART 3: OVERFITTING GAP ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Calculate gaps (Train - Test)
gap_auc <- train_auc - test_auc
gap_bal_acc <- train_bal_acc - test_bal_acc
gap_sens <- train_sens - test_sens
gap_spec <- train_spec - test_spec
gap_f1 <- train_f1 - test_f1

cat("PERFORMANCE GAP (Train - Test):\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  ΔAUC:                %.4f\n", gap_auc))
cat(sprintf("  ΔBalanced Accuracy:  %.4f\n", gap_bal_acc))
cat(sprintf("  ΔSensitivity:        %.4f\n", gap_sens))
cat(sprintf("  ΔSpecificity:        %.4f\n", gap_spec))
cat(sprintf("  ΔF1-Score:           %.4f\n\n", gap_f1))

cat("INTERPRETATION:\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Overfitting thresholds
if (abs(gap_auc) < 0.05) {
  overfitting_risk <- "VERY LOW ✓✓✓"
  recommendation <- "Excellent generalization. Model is NOT overfitting."
  risk_level <- "★★★★★ MINIMAL RISK"
} else if (abs(gap_auc) < 0.10) {
  overfitting_risk <- "LOW ✓✓"
  recommendation <- "Good generalization. Minor gap is acceptable."
  risk_level <- "★★★★ LOW RISK"
} else if (abs(gap_auc) < 0.15) {
  overfitting_risk <- "MODERATE ⚠"
  recommendation <- "Noticeable gap. Consider more regularization."
  risk_level <- "★★★ MODERATE RISK"
} else {
  overfitting_risk <- "HIGH ✗✗"
  recommendation <- "Significant gap indicates overfitting."
  risk_level <- "★★ HIGH RISK"
}

cat(sprintf("  Overfitting Risk:    %s\n", overfitting_risk))
cat(sprintf("  Risk Level:          %s\n", risk_level))
cat(sprintf("  AUC Gap:             %.4f (%.1f%%)\n", gap_auc, 100*abs(gap_auc)/train_auc))
cat("\n")
cat(sprintf("  %s\n", recommendation))
cat("\n")

# Additional checks
cat("ADDITIONAL CHECKS:\n")
cat("─────────────────────────────────────────────────────────────────\n")

if (test_auc >= 0.80 * train_auc) {
  cat("  ✓ Test AUC >= 80% of Train AUC (PASS)\n")
} else {
  cat("  ✗ Test AUC < 80% of Train AUC (FAIL)\n")
}

if (test_auc >= 0.80) {
  cat("  ✓ Test AUC >= 0.80 (Strong performance)\n")
} else if (test_auc >= 0.70) {
  cat("  ✓ Test AUC >= 0.70 (Moderate performance)\n")
} else {
  cat("  ✗ Test AUC < 0.70 (Weak performance)\n")
}

if (train_auc < 0.99) {
  cat("  ✓ Train AUC < 0.99 (Not perfect - realistic)\n")
} else {
  cat("  ⚠ Train AUC >= 0.99 (Very high - suspicious)\n")
}

consistency_score <- 0
if (abs(gap_bal_acc) < 0.10) consistency_score <- consistency_score + 1
if (abs(gap_sens) < 0.10) consistency_score <- consistency_score + 1
if (abs(gap_spec) < 0.10) consistency_score <- consistency_score + 1
if (abs(gap_f1) < 0.10) consistency_score <- consistency_score + 1

cat(sprintf("  ✓ Consistency score: %d/4 metrics have gap < 10%%\n", consistency_score))

if (consistency_score >= 3) {
  cat("    → Consistent performance across metrics\n")
} else {
  cat("    → Some inconsistency across metrics\n")
}

cat("\n")

# ═══════════════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════\n")
cat("                    FINAL VERDICT\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Load permutation p-value if available
perf_metrics <- read.csv("Single_Omics_mRNA_Results/Performance_Metrics.csv")
perm_pvalue <- perf_metrics[perf_metrics$Metric == "Permutation_pvalue", "Value"]

# Calculate overall score
score <- 0
reasons_pass <- c()
reasons_concern <- c()

if (abs(gap_auc) < 0.05) {
  score <- score + 1
  reasons_pass <- c(reasons_pass, "AUC gap < 5%")
} else if (abs(gap_auc) < 0.10) {
  score <- score + 0.5
  reasons_pass <- c(reasons_pass, "AUC gap < 10%")
} else {
  reasons_concern <- c(reasons_concern, sprintf("AUC gap %.1f%%", 100*abs(gap_auc)))
}

if (perm_pvalue < 0.01) {
  score <- score + 1
  reasons_pass <- c(reasons_pass, "Permutation p<0.01")
} else if (perm_pvalue < 0.05) {
  score <- score + 0.5
  reasons_pass <- c(reasons_pass, "Permutation p<0.05")
} else {
  reasons_concern <- c(reasons_concern, "Permutation not significant")
}

if (test_auc >= 0.80 && test_auc <= 0.95) {
  score <- score + 1
  reasons_pass <- c(reasons_pass, "Test AUC in healthy range")
} else if (test_auc > 0.95) {
  score <- score + 0.5
  reasons_concern <- c(reasons_concern, "Test AUC very high")
} else {
  reasons_concern <- c(reasons_concern, "Test AUC suboptimal")
}

if (train_auc < 0.99) {
  score <- score + 1
  reasons_pass <- c(reasons_pass, "Train AUC not perfect")
} else {
  reasons_concern <- c(reasons_concern, "Train AUC too high")
}

if (consistency_score >= 3) {
  score <- score + 1
  reasons_pass <- c(reasons_pass, "Metrics consistent")
} else {
  score <- score + 0.5
  reasons_concern <- c(reasons_concern, "Some metric inconsistency")
}

stars <- paste(rep("★", floor(score)), collapse = "")
if (score %% 1 >= 0.5) stars <- paste0(stars, "☆")

cat(sprintf("OVERALL OVERFITTING RISK: %s (%.1f/5.0)\n\n", stars, score))

cat("STRENGTHS:\n")
for (reason in reasons_pass) {
  cat(sprintf("  ✓ %s\n", reason))
}
cat("\n")

if (length(reasons_concern) > 0) {
  cat("CONCERNS:\n")
  for (reason in reasons_concern) {
    cat(sprintf("  ⚠ %s\n", reason))
  }
  cat("\n")
}

# Final decision
if (score >= 4.5) {
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  CONCLUSION: ✓✓✓ NO OVERFITTING DETECTED\n")
  cat("  CONFIDENCE: VERY HIGH\n")
  cat("  RECOMMENDATION: Results are publication-ready\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
} else if (score >= 3.5) {
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  CONCLUSION: ✓✓ MINIMAL OVERFITTING RISK\n")
  cat("  CONFIDENCE: HIGH\n")
  cat("  RECOMMENDATION: Results are suitable for publication\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
} else if (score >= 2.5) {
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  CONCLUSION: ⚠ MODERATE OVERFITTING RISK\n")
  cat("  CONFIDENCE: MODERATE\n")
  cat("  RECOMMENDATION: Results should be validated on external cohort\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
} else {
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("  CONCLUSION: ✗✗ HIGH OVERFITTING RISK\n")
  cat("  CONFIDENCE: LOW\n")
  cat("  RECOMMENDATION: Model needs significant improvement\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("      STRICT OVERFITTING CHECK COMPLETE (mRNA)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

