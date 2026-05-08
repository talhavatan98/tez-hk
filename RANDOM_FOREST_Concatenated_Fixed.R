#!/usr/bin/env Rscript
# ════════════════════════════════════════════════════════════════════════════
#  RANDOM FOREST ANALYSIS - Concatenated Data
#  Benchmark comparison for DIABLO
# ════════════════════════════════════════════════════════════════════════════

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║  RANDOM FOREST ANALYSIS - Concatenated Data                    ║\n")
cat("║  Benchmark for Multi-Omics Integration                        ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

suppressPackageStartupMessages({
  if (!require("randomForest")) install.packages("randomForest", repos="https://cloud.r-project.org")
  library(randomForest)
  library(caret)
  library(pROC)
  library(PRROC)
  library(ggplot2)
})

# CRITICAL: Fixed seed for full reproducibility (Same as DIABLO)
set.seed(42)
options(width = 120)

# Configuration (Matches DIABLO_Final_Honest.R)
CONFIG <- list(
  test_size = 0.30,
  cv_folds = 10,
  cv_repeats = 5,
  max_features_prefilter = 2000,
  output_dir = "Random_Forest_Concatenated_Results"
)

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

# ═══════════════════════════════════════════════════════════════════════════
# 1. DATA LOADING & ALIGNMENT
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 1: Loading Data\n")
cat("─────────────────────────────────────────────────────────────────\n")

load("DIABLO_ready_data.RData")
Y_all <- Y
rm(Y, X)

temp_env <- new.env()
load("mRNA_processed_for_DIABLO.RData", envir = temp_env)
X_mrna <- temp_env$mrna_final

temp_env <- new.env()
load("methylation_processed_for_DIABLO.RData", envir = temp_env)
X_meth <- temp_env$meth_final

# Align samples
common_samples <- Reduce(intersect, list(names(Y_all), rownames(X_mrna), rownames(X_meth)))
Y_all <- Y_all[common_samples]
X_mrna <- X_mrna[common_samples, ]
X_meth <- X_meth[common_samples, ]

cat(sprintf("✓ Loaded: %d samples, %d mRNA, %d Methylation features\n\n",
            length(Y_all), ncol(X_mrna), ncol(X_meth)))

# ═══════════════════════════════════════════════════════════════════════════
# 2. STRATIFIED TRAIN-TEST SPLIT (FIRST - TO PREVENT LEAKAGE)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 2: Stratified Train-Test Split (Seed=42) - BEFORE FILTERING\n")
cat("─────────────────────────────────────────────────────────────────\n")

set.seed(42)
train_idx <- createDataPartition(Y_all, p = 1 - CONFIG$test_size, list = FALSE)[,1]

Y_train <- Y_all[train_idx]
Y_test <- Y_all[-train_idx]

# Split RAW data
X_mrna_train_raw <- X_mrna[train_idx, ]
X_mrna_test_raw <- X_mrna[-train_idx, ]

X_meth_train_raw <- X_meth[train_idx, ]
X_meth_test_raw <- X_meth[-train_idx, ]

cat(sprintf("✓ Training: %d samples\n", length(Y_train)))
cat(sprintf("✓ Testing: %d samples\n\n", length(Y_test)))
# ═══════════════════════════════════════════════════════════════════════════
# 3. FEATURE PRE-FILTERING (TRAIN ONLY)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 3: Feature Pre-Filtering (TRAIN ONLY - No Leakage)\n")
cat("─────────────────────────────────────────────────────────────────\n")

filter_features <- function(X, Y, max_features = 2000) {
  feat_vars <- apply(X, 2, var, na.rm = TRUE)
  
  if (ncol(X) > max_features) {
    # T-test filter removed to prevent dominance. Only using unsupervised variance filter.
    top_idx <- order(feat_vars, decreasing = TRUE)[1:max_features]
    X <- X[, top_idx, drop = FALSE]
  } else {
    keep_var <- feat_vars > quantile(feat_vars, 0.20)
    X <- X[, keep_var, drop = FALSE]
  }
  return(X)
}

cat("Filtering mRNA (Train Only):\n")
X_mrna_train_filtered <- filter_features(X_mrna_train_raw, Y_train, CONFIG$max_features_prefilter)
# Apply to test
X_mrna_test_filtered <- X_mrna_test_raw[, colnames(X_mrna_train_filtered)]

cat("\nFiltering Methylation (Train Only):\n")
X_meth_train_filtered <- filter_features(X_meth_train_raw, Y_train, CONFIG$max_features_prefilter)
# Apply to test
X_meth_test_filtered <- X_meth_test_raw[, colnames(X_meth_train_filtered)]

cat("\n")

# ═══════════════════════════════════════════════════════════════════════════
# 4. SCALING & CONCATENATION
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 4: Scaling & Concatenating Data\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Scale based on TRAIN statistics (avoid leakage)
scale_mrna <- preProcess(X_mrna_train_filtered, method = c("center", "scale"))
X_mrna_train_scaled <- predict(scale_mrna, X_mrna_train_filtered)
X_mrna_test_scaled <- predict(scale_mrna, X_mrna_test_filtered)

scale_meth <- preProcess(X_meth_train_filtered, method = c("center", "scale"))
X_meth_train_scaled <- predict(scale_meth, X_meth_train_filtered)
X_meth_test_scaled <- predict(scale_meth, X_meth_test_filtered)

# Ensure unique column names
colnames(X_mrna_train_scaled) <- paste0("mRNA_", colnames(X_mrna_train_scaled))
colnames(X_mrna_test_scaled) <- paste0("mRNA_", colnames(X_mrna_test_scaled))
colnames(X_meth_train_scaled) <- paste0("Meth_", colnames(X_meth_train_scaled))
colnames(X_meth_test_scaled) <- paste0("Meth_", colnames(X_meth_test_scaled))

# Concatenate
X_train <- cbind(X_mrna_train_scaled, X_meth_train_scaled)
X_test <- cbind(X_mrna_test_scaled, X_meth_test_scaled)

cat(sprintf("✓ Concatenated Train Matrix: %d samples x %d features\n", nrow(X_train), ncol(X_train)))
cat(sprintf("✓ Concatenated Test Matrix: %d samples x %d features\n\n", nrow(X_test), ncol(X_test)))



# ═══════════════════════════════════════════════════════════════════════════
# 5. MODEL TRAINING (Random Forest)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 5: Training Random Forest\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Using caret for consistent CV
ctrl <- trainControl(
  method = "repeatedcv",
  number = CONFIG$cv_folds,
  repeats = CONFIG$cv_repeats,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

set.seed(42)
rf_model <- train(
  x = X_train,
  y = Y_train,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5, # Tune mtry
  ntree = 500
)

cat("✓ Model trained. Best tuning parameters:\n")
print(rf_model$bestTune)

# Get OOB error (internal RF validation - more realistic than train performance)
final_rf <- rf_model$finalModel
oob_error <- final_rf$err.rate[final_rf$ntree, "OOB"]
oob_accuracy <- 1 - oob_error

cat(sprintf("\n  OOB Error Rate:     %.4f\n", oob_error))
cat(sprintf("  OOB Accuracy:       %.4f (%.1f%%)\n", oob_accuracy, 100*oob_accuracy))
cat("  Note: OOB is RF's internal validation (more realistic than train set)\n")
cat("\n")

# ═══════════════════════════════════════════════════════════════════════════
# 6. TRAINING SET EVALUATION (Overfitting Control)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 5.5: Training Set Performance (Overfitting Control)\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Get training predictions
train_pred_classes <- predict(rf_model, newdata = X_train)
train_pred_probs <- predict(rf_model, newdata = X_train, type = "prob")
train_pred_scores <- train_pred_probs[, "Mutant"]

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

cat(sprintf("  Training Samples:       %d\n", length(Y_train)))
cat(sprintf("  Training Balanced Acc:  %.4f (%.1f%%)\n", train_bal_acc, 100*train_bal_acc))
cat(sprintf("  Training ROC-AUC:       %.4f\n", train_auc))

if (train_auc >= 0.99) {
  cat("\n  ⚠ NOTE: Train AUC ≈ 1.0 is NORMAL for Random Forest\n")
  cat("    RF can memorize training data with deep trees\n")
  cat("    This is NOT overfitting if test performance is good\n")
  cat(sprintf("    OOB Accuracy (more realistic): %.1f%%\n", 100*oob_accuracy))
}
cat("\n")

# ═══════════════════════════════════════════════════════════════════════════
# 7. TEST SET EVALUATION (EXTENDED METRICS)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 6: Test Set Evaluation\n")
cat("─────────────────────────────────────────────────────────────────\n")

pred_classes <- predict(rf_model, newdata = X_test)
pred_probs <- predict(rf_model, newdata = X_test, type = "prob")
pred_scores <- pred_probs[, "Mutant"]

# Confusion Matrix
cm <- confusionMatrix(
  factor(pred_classes, levels = c("Mutant", "WT")),
  factor(Y_test, levels = c("Mutant", "WT")),
  positive = "Mutant"
)

# Metrics (Match DIABLO)
balanced_accuracy <- cm$byClass["Balanced Accuracy"]
sensitivity <- cm$byClass["Sensitivity"]
specificity <- cm$byClass["Specificity"]
precision <- cm$byClass["Precision"]
f1_score <- cm$byClass["F1"]

# MCC
TP <- cm$table["Mutant", "Mutant"]
TN <- cm$table["WT", "WT"]
FP <- cm$table["Mutant", "WT"]
FN <- cm$table["WT", "Mutant"]
mcc <- (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

# ROC & PR AUC
roc_obj <- roc(Y_test, pred_scores, levels = c("WT", "Mutant"), direction = "<", quiet = TRUE)
roc_auc <- as.numeric(auc(roc_obj))

binary_labels <- ifelse(Y_test == "Mutant", 1, 0)
pr_obj <- pr.curve(scores.class0 = pred_scores[binary_labels == 1],
                   scores.class1 = pred_scores[binary_labels == 0],
                   curve = TRUE)
pr_auc <- pr_obj$auc.integral

cat(sprintf("  Test Samples:           %d\n", length(Y_test)))
cat(sprintf("  Balanced Accuracy:      %.4f (%.1f%%)\n", balanced_accuracy, 100*balanced_accuracy))
cat(sprintf("  Sensitivity:            %.4f (%.1f%%)\n", sensitivity, 100*sensitivity))
cat(sprintf("  Specificity:            %.4f (%.1f%%)\n", specificity, 100*specificity))
cat(sprintf("  F1-Score:               %.4f\n", f1_score))
cat(sprintf("  MCC:                    %.4f\n", mcc))
cat(sprintf("  ROC-AUC:                %.4f\n", roc_auc))
cat(sprintf("  PR-AUC:                 %.4f\n\n", pr_auc))

# ═══════════════════════════════════════════════════════════════════════════
# 7.5. TRAIN-TEST GAP ANALYSIS (Overfitting Assessment)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 6.5: Train-Test Gap Analysis (Overfitting Check)\n")
cat("─────────────────────────────────────────────────────────────────\n")

gap_auc <- train_auc - roc_auc
gap_bal_acc <- train_bal_acc - balanced_accuracy

cat(sprintf("  Train AUC:              %.4f\n", train_auc))
cat(sprintf("  Test AUC:               %.4f\n", roc_auc))
cat(sprintf("  ΔAUC (Train-Test):      %.4f\n", gap_auc))
cat(sprintf("  ΔBalanced Acc:          %.4f\n\n", gap_bal_acc))

if (abs(gap_auc) < 0.05) {
  cat("  ✓✓✓ OVERFITTING RISK: VERY LOW (gap < 5%%)\n")
} else if (abs(gap_auc) < 0.10) {
  cat("  ✓✓ OVERFITTING RISK: LOW (gap < 10%%)\n")
} else if (abs(gap_auc) < 0.15) {
  cat("  ⚠ OVERFITTING RISK: MODERATE (gap < 15%%)\n")
} else {
  cat("  ✗✗ OVERFITTING RISK: HIGH (gap >= 15%%)\n")
}
cat("\n")

# ═══════════════════════════════════════════════════════════════════════════
# 8. PERMUTATION TEST (Statistical Validation)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 7: Permutation Test (200 permutations)\n")
cat("─────────────────────────────────────────────────────────────────\n")

n_perm <- 200
perm_aucs <- numeric(n_perm)

cat("Running permutations (this may take a few minutes)...\n")
pb <- txtProgressBar(min = 0, max = n_perm, style = 3, width = 50)

set.seed(42)  # Reproducible permutations
for (i in 1:n_perm) {
  Y_perm <- sample(Y_train)
  suppressWarnings({
    perm_model <- train(
      x = X_train,
      y = Y_perm,
      method = "rf",
      trControl = trainControl(method = "none"),  # No CV for permutation
      tuneGrid = rf_model$bestTune,  # Use same tuning parameters
      ntree = 500
    )
  })
  perm_probs <- predict(perm_model, newdata = X_test, type = "prob")
  perm_scores <- perm_probs[, "Mutant"]
  perm_roc <- roc(Y_test, perm_scores, levels = c("WT", "Mutant"), direction = "<", quiet = TRUE)
  perm_aucs[i] <- as.numeric(auc(perm_roc))
  setTxtProgressBar(pb, i)
}
close(pb)

# Calculate p-value
n_exceed <- sum(perm_aucs >= roc_auc)
perm_pvalue <- (n_exceed + 1) / (n_perm + 1)

cat("\n\n✓ Permutation Test Results:\n")
cat(sprintf("  • Real AUC:         %.4f\n", roc_auc))
cat(sprintf("  • Mean Perm AUC:    %.4f ± %.4f\n", mean(perm_aucs), sd(perm_aucs)))
cat(sprintf("  • N exceeding real: %d / %d\n", n_exceed, n_perm))
cat(sprintf("  • P-value:          %.4f\n", perm_pvalue))

if (perm_pvalue < 0.01) {
  cat("  ★★★ HIGHLY SIGNIFICANT (p < 0.01)\n")
} else if (perm_pvalue < 0.05) {
  cat("  ★★ SIGNIFICANT (p < 0.05)\n")
} else {
  cat("  ✗ NOT SIGNIFICANT (p ≥ 0.05)\n")
}
cat("\n")

# ═══════════════════════════════════════════════════════════════════════════
# 9. SAVE RESULTS (Including Overfitting Controls)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 8: Saving Results\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Save Plots
pdf(file.path(CONFIG$output_dir, "ROC_Curve.pdf"), width = 7, height = 6)
plot(roc_obj, main = sprintf("ROC Curve (AUC = %.4f)", roc_auc), col = "#2E86AB", lwd = 3)
abline(a = 0, b = 1, lty = 2, col = "gray50", lwd = 2)
dev.off()

# Save Permutation Test Plot
pdf(file.path(CONFIG$output_dir, "Permutation_Test.pdf"), width = 10, height = 7)
hist(perm_aucs, breaks = 30, col = "#E8E8E8", border = "gray30",
     main = sprintf("Permutation Test (p = %.4f)", perm_pvalue),
     xlab = "AUC from Permuted Labels", cex.main = 1.8, cex.lab = 1.4)
abline(v = roc_auc, col = "#2E86AB", lwd = 4)
abline(v = mean(perm_aucs), col = "#F18F01", lwd = 3, lty = 2)
legend("topright", legend = c(sprintf("Real (%.4f)", roc_auc),
                              sprintf("Perm Mean (%.4f)", mean(perm_aucs))),
       col = c("#2E86AB", "#F18F01"), lty = c(1, 2), lwd = c(4, 3))
dev.off()

# Save Train vs Test Comparison Plot
pdf(file.path(CONFIG$output_dir, "Train_vs_Test_Gap_Visual.pdf"), width = 10, height = 7)
par(mar = c(5, 5, 4, 2))
metrics_names <- c("Balanced\nAccuracy", "Sensitivity", "Specificity", "F1-Score", "ROC-AUC")
train_vals <- c(train_bal_acc, train_sens, train_spec, train_f1, train_auc)
test_vals <- c(balanced_accuracy, sensitivity, specificity, f1_score, roc_auc)
x_pos <- barplot(rbind(train_vals, test_vals), beside = TRUE, 
                 names.arg = metrics_names,
                 col = c("#2E86AB", "#F18F01"), ylim = c(0, 1),
                 main = "Train vs Test Performance", 
                 ylab = "Score", cex.main = 1.8, cex.lab = 1.4, cex.names = 1.2)
legend("bottomright", legend = c("Train", "Test"), 
       fill = c("#2E86AB", "#F18F01"), cex = 1.3)
dev.off()

# Save Confusion Matrix Plot
pdf(file.path(CONFIG$output_dir, "Confusion_Matrix.pdf"), width = 7, height = 6)
cm_df <- as.data.frame(cm$table)
cm_df$Prediction <- factor(cm_df$Prediction, levels = c("WT", "Mutant"))
cm_df$Reference <- factor(cm_df$Reference, levels = c("Mutant", "WT"))
midpoint <- mean(range(cm_df$Freq))
cm_df$TextCol <- ifelse(cm_df$Freq > midpoint, "white", "black")

p_cm <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Freq, color = TextCol), size = 10, fontface = "bold") +
  scale_color_identity() +
  scale_fill_gradient(low = "#86BBD8", high = "#2E86AB") +
  labs(title = "Confusion Matrix (RF Concatenated)",
       subtitle = sprintf("Balanced Accuracy: %.4f", balanced_accuracy),
       x = "True Class", y = "Predicted Class") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")
print(p_cm)
dev.off()

# Save Metrics (Including Permutation p-value)
write.csv(data.frame(
  Metric = c("Balanced_Accuracy", "Sensitivity", "Specificity", 
             "Precision", "F1_Score", "MCC", "ROC_AUC", "PR_AUC", "Permutation_pvalue"),
  Value = c(balanced_accuracy, sensitivity, specificity, 
            precision, f1_score, mcc, roc_auc, pr_auc, perm_pvalue)
), file.path(CONFIG$output_dir, "Performance_Metrics.csv"), row.names = FALSE)

# Save Train vs Test Comparison
comparison_df <- data.frame(
  Dataset = c("Training", "Test", "Gap (Train-Test)"),
  Samples = c(length(Y_train), length(Y_test), NA),
  AUC = c(train_auc, roc_auc, gap_auc),
  Balanced_Accuracy = c(train_bal_acc, balanced_accuracy, gap_bal_acc),
  Sensitivity = c(train_sens, sensitivity, train_sens - sensitivity),
  Specificity = c(train_spec, specificity, train_spec - specificity),
  F1_Score = c(train_f1, f1_score, train_f1 - f1_score)
)
write.csv(comparison_df, file.path(CONFIG$output_dir, "Train_vs_Test_Comparison.csv"), row.names = FALSE)

# Save Variable Importance
var_imp <- varImp(rf_model)$importance
var_imp$Feature <- rownames(var_imp)
var_imp <- var_imp[order(var_imp$Overall, decreasing = TRUE), ]
write.csv(var_imp, file.path(CONFIG$output_dir, "Variable_Importance.csv"), row.names = FALSE)

# Save workspace for statistical comparison
save(rf_model, X_train, X_test, Y_train, Y_test, pred_classes, pred_probs,
     file = file.path(CONFIG$output_dir, "Model_Workspace.RData"))

cat(sprintf("✓ Results saved to %s/\n", CONFIG$output_dir))
cat("  • Performance_Metrics.csv (including permutation p-value)\n")
cat("  • Train_vs_Test_Comparison.csv (overfitting assessment)\n")
cat("  • Variable_Importance.csv\n")
cat("  • ROC_Curve.pdf\n")
cat("  • Confusion_Matrix.pdf\n")
cat("  • Permutation_Test.pdf\n")
cat("  • Train_vs_Test_Gap_Visual.pdf\n\n")

cat("═══════════════════════════════════════════════════════════════\n")
cat("             ANALYSIS COMPLETE - SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("PERFORMANCE SUMMARY:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  • Test AUC:             %.4f\n", roc_auc))
cat(sprintf("  • Test Balanced Acc:    %.1f%%\n", 100*balanced_accuracy))
cat(sprintf("  • Permutation p-value:  %.4f ", perm_pvalue))
if (perm_pvalue < 0.05) {
  cat("(SIGNIFICANT ★)\n")
} else {
  cat("(not significant)\n")
}
cat("\n")

cat("OVERFITTING ASSESSMENT:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat(sprintf("  • Train AUC:            %.4f", train_auc))
if (train_auc >= 0.99) {
  cat(" (Expected for RF)")
}
cat("\n")
cat(sprintf("  • OOB Accuracy:         %.4f (RF internal validation)\n", oob_accuracy))
cat(sprintf("  • Test AUC:             %.4f\n", roc_auc))
cat(sprintf("  • Gap (Train-Test):     %.4f (%.1f%%)\n", gap_auc, 100*abs(gap_auc)/train_auc))
if (abs(gap_auc) < 0.05) {
  cat("  • Risk Level:           ✓✓✓ VERY LOW\n")
} else if (abs(gap_auc) < 0.10) {
  cat("  • Risk Level:           ✓✓ LOW\n")
} else if (abs(gap_auc) < 0.15) {
  cat("  • Risk Level:           ⚠ MODERATE\n")
} else {
  cat("  • Risk Level:           ✗✗ HIGH\n")
}
cat("\n")

cat("METHODOLOGY:\n")
cat("  ✓ Fixed seed (42) - fully reproducible\n")
cat("  ✓ Stratified train-test split (30%%)\n")
cat("  ✓ Feature selection on train only (no leakage)\n")
cat("  ✓ Cross-validation (10-fold × 5 repeats)\n")
cat("  ✓ Permutation test (200 iterations)\n")
cat("  ✓ Train-test gap analysis\n")
cat("  ✓ OOB validation (RF's internal CV)\n")
cat("\n")
cat("NOTE: Random Forest can achieve near-perfect training accuracy\n")
cat("      This is NOT overfitting if test performance remains strong\n")
cat("      OOB error provides more realistic training estimate\n")
cat("\n")

cat("═══════════════════════════════════════════════════════════════\n")
cat("READY FOR PUBLICATION ✓\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
