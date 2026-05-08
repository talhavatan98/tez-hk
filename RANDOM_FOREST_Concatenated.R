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
  max_features_prefilter = 200,
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

cat("Filtering mRNA:\n")
X_mrna_filtered <- filter_features(X_mrna, Y_all, CONFIG$max_features_prefilter)
cat(sprintf("  ✓ mRNA features: %d\n", ncol(X_mrna_filtered)))

cat("Filtering Methylation:\n")
X_meth_filtered <- filter_features(X_meth, Y_all, CONFIG$max_features_prefilter)
cat(sprintf("  ✓ Methylation features: %d\n\n", ncol(X_meth_filtered)))

# ═══════════════════════════════════════════════════════════════════════════
# 3. SCALING & CONCATENATION (CRITICAL!)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 3: Scaling & Concatenating Data\n")
cat("─────────────────────────────────────────────────────────────────\n")

# CRITICAL: Scale each omics separately BEFORE concatenation
# This ensures that mRNA (high values) doesn't dominate Methylation (0-1)
X_mrna_scaled <- scale(X_mrna_filtered)
X_meth_scaled <- scale(X_meth_filtered)

# Ensure unique column names
colnames(X_mrna_scaled) <- paste0("mRNA_", colnames(X_mrna_scaled))
colnames(X_meth_scaled) <- paste0("Meth_", colnames(X_meth_scaled))

X_concat <- cbind(X_mrna_scaled, X_meth_scaled)

cat(sprintf("✓ Scaled mRNA: mean=%.2f, sd=%.2f\n", 
            mean(X_mrna_scaled, na.rm=TRUE), sd(X_mrna_scaled, na.rm=TRUE)))
cat(sprintf("✓ Scaled Meth: mean=%.2f, sd=%.2f\n", 
            mean(X_meth_scaled, na.rm=TRUE), sd(X_meth_scaled, na.rm=TRUE)))
cat(sprintf("✓ Concatenated Matrix: %d samples x %d features\n\n",
            nrow(X_concat), ncol(X_concat)))

# ═══════════════════════════════════════════════════════════════════════════
# 4. TRAIN-TEST SPLIT (Same Seed = 42)
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 4: Stratified Train-Test Split (Seed=42)\n")
cat("─────────────────────────────────────────────────────────────────\n")

set.seed(42)
train_idx <- createDataPartition(Y_all, p = 1 - CONFIG$test_size, list = FALSE)[,1]

Y_train <- Y_all[train_idx]
Y_test <- Y_all[-train_idx]
X_train <- X_concat[train_idx, ]
X_test <- X_concat[-train_idx, ]

cat(sprintf("✓ Training: %d samples\n", length(Y_train)))
cat(sprintf("✓ Testing: %d samples\n\n", length(Y_test)))

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
cat("\n")

# ═══════════════════════════════════════════════════════════════════════════
# 6. EVALUATION (EXTENDED METRICS)
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

cat(sprintf("  Balanced Accuracy:  %.4f (%.1f%%)\n", balanced_accuracy, 100*balanced_accuracy))
cat(sprintf("  Sensitivity:        %.4f (%.1f%%)\n", sensitivity, 100*sensitivity))
cat(sprintf("  Specificity:        %.4f (%.1f%%)\n", specificity, 100*specificity))
cat(sprintf("  F1-Score:           %.4f\n", f1_score))
cat(sprintf("  MCC:                %.4f\n", mcc))
cat(sprintf("  ROC-AUC:            %.4f\n", roc_auc))
cat(sprintf("  PR-AUC:             %.4f\n\n", pr_auc))

# ═══════════════════════════════════════════════════════════════════════════
# 7. SAVE RESULTS
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
  labs(title = "Confusion Matrix (RF Concatenated)",
       subtitle = sprintf("Balanced Accuracy: %.4f", balanced_accuracy),
       x = "True Class", y = "Predicted Class") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")
print(p_cm)
dev.off()

# Save Metrics
write.csv(data.frame(
  Metric = c("Balanced_Accuracy", "Sensitivity", "Specificity", 
             "Precision", "F1_Score", "MCC", "ROC_AUC", "PR_AUC"),
  Value = c(balanced_accuracy, sensitivity, specificity, 
            precision, f1_score, mcc, roc_auc, pr_auc)
), file.path(CONFIG$output_dir, "Performance_Metrics.csv"), row.names = FALSE)

# Save Variable Importance
var_imp <- varImp(rf_model)$importance
var_imp$Feature <- rownames(var_imp)
var_imp <- var_imp[order(var_imp$Overall, decreasing = TRUE), ]
write.csv(var_imp, file.path(CONFIG$output_dir, "Variable_Importance.csv"), row.names = FALSE)

cat(sprintf("✓ Results saved to %s/\n", CONFIG$output_dir))
