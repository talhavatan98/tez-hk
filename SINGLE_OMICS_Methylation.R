#!/usr/bin/env Rscript
# ════════════════════════════════════════════════════════════════════════════
#  SINGLE OMICS ANALYSIS - Methylation (sPLS-DA)
#  Benchmark comparison for DIABLO
# ════════════════════════════════════════════════════════════════════════════

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║  SINGLE OMICS ANALYSIS - Methylation (sPLS-DA)                 ║\n")
cat("║  Benchmark for Multi-Omics Integration                        ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

suppressPackageStartupMessages({
  library(mixOmics)
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
  max_features_prefilter = 200,
  keepX_grid = c(3, 5, 8, 10, 15, 20, 30, 50), # Compatible with DIABLO + extended range
  max_ncomp = 2,
  output_dir = "Single_Omics_Methylation_Results"
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
load("methylation_processed_for_DIABLO.RData", envir = temp_env)
X_meth <- temp_env$meth_final

# Align samples
common_samples <- intersect(names(Y_all), rownames(X_meth))
Y_all <- Y_all[common_samples]
X_meth <- X_meth[common_samples, ]

cat(sprintf("✓ Loaded: %d samples, %d Methylation features\n\n",
            length(Y_all), ncol(X_meth)))

# ═══════════════════════════════════════════════════════════════════════════
# 2. FEATURE PRE-FILTERING (Same as DIABLO)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 2: Feature Pre-Filtering\n")
cat("─────────────────────────────────────────────────────────────────\n")

filter_features <- function(X, Y, max_features = 200) {
  feat_vars <- apply(X, 2, var, na.rm = TRUE)
  keep_var <- feat_vars > quantile(feat_vars, 0.20)
  X <- X[, keep_var, drop = FALSE]
  
  if (ncol(X) > max_features) {
    p_values <- apply(X, 2, function(feat) {
      tryCatch(t.test(feat ~ Y)$p.value, error = function(e) 1.0)
    })
    top_idx <- order(p_values)[1:min(max_features, length(p_values))]
    X <- X[, top_idx, drop = FALSE]
  }
  return(X)
}

X_meth_filtered <- filter_features(X_meth, Y_all, CONFIG$max_features_prefilter)
cat(sprintf("✓ Filtered to %d features\n\n", ncol(X_meth_filtered)))

# ═══════════════════════════════════════════════════════════════════════════
# 3. TRAIN-TEST SPLIT (Same Seed = 42)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 3: Stratified Train-Test Split (Seed=42)\n")
cat("─────────────────────────────────────────────────────────────────\n")

set.seed(42)
train_idx <- createDataPartition(Y_all, p = 1 - CONFIG$test_size, list = FALSE)[,1]

Y_train <- Y_all[train_idx]
Y_test <- Y_all[-train_idx]
X_train <- X_meth_filtered[train_idx, ]
X_test <- X_meth_filtered[-train_idx, ]

cat(sprintf("✓ Training: %d samples\n", length(Y_train)))
cat(sprintf("✓ Testing: %d samples\n\n", length(Y_test)))

# ═══════════════════════════════════════════════════════════════════════════
# 4. MODEL TUNING & TRAINING
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 4: Tuning sPLS-DA\n")
cat("─────────────────────────────────────────────────────────────────\n")

tune_splsda <- tune.splsda(
  X = X_train, Y = Y_train, ncomp = CONFIG$max_ncomp,
  test.keepX = CONFIG$keepX_grid,
  validation = "Mfold", folds = CONFIG$cv_folds,
  nrepeat = CONFIG$cv_repeats, dist = "centroids.dist",
  measure = "BER", scale = TRUE, progressBar = FALSE
)

optimal_keepX <- tune_splsda$choice.keepX
cat("✓ Optimal keepX:\n")
print(optimal_keepX)

cat("\nSTEP 5: Training Final Model\n")
cat("─────────────────────────────────────────────────────────────────\n")

final_model <- mixOmics::splsda(
  X = X_train, Y = Y_train, ncomp = CONFIG$max_ncomp,
  keepX = optimal_keepX, scale = TRUE
)

# ═══════════════════════════════════════════════════════════════════════════
# 5. EVALUATION (EXTENDED METRICS)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 6: Test Set Evaluation\n")
cat("─────────────────────────────────────────────────────────────────\n")

predictions <- predict(final_model, newdata = X_test)
pred_classes <- predictions$class$centroids.dist[, CONFIG$max_ncomp]
pred_scores <- predictions$predict[, "Mutant", CONFIG$max_ncomp]

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

cat(sprintf("  Balanced Accuracy:  %.4f (%.1f%%)\n", balanced_accuracy, 100*balanced_accuracy))
cat(sprintf("  Sensitivity:        %.4f (%.1f%%)\n", sensitivity, 100*sensitivity))
cat(sprintf("  Specificity:        %.4f (%.1f%%)\n", specificity, 100*specificity))
cat(sprintf("  F1-Score:           %.4f\n", f1_score))
cat(sprintf("  MCC:                %.4f\n", mcc))
cat(sprintf("  ROC-AUC:            %.4f\n", roc_auc))
cat(sprintf("  PR-AUC:             %.4f\n\n", pr_auc))

# ═══════════════════════════════════════════════════════════════════════════
# 6. SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════════════

# Save Plots
pdf(file.path(CONFIG$output_dir, "ROC_Curve.pdf"), width = 7, height = 6)
plot(roc_obj, main = sprintf("ROC Curve (AUC = %.4f)", roc_auc), col = "#2E86AB", lwd = 3)
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
  labs(title = "Confusion Matrix (Methylation Only)",
       subtitle = sprintf("Balanced Accuracy: %.4f", balanced_accuracy),
       x = "True Class", y = "Predicted Class") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")
print(p_cm)
dev.off()

# Save ALL Metrics
write.csv(data.frame(
  Metric = c("Balanced_Accuracy", "Sensitivity", "Specificity", 
             "Precision", "F1_Score", "MCC", "ROC_AUC", "PR_AUC"),
  Value = c(balanced_accuracy, sensitivity, specificity, 
            precision, f1_score, mcc, roc_auc, pr_auc)
), file.path(CONFIG$output_dir, "Performance_Metrics.csv"), row.names = FALSE)

cat(sprintf("✓ Results saved to %s/\n", CONFIG$output_dir))
