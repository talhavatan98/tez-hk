#!/usr/bin/env Rscript
# ════════════════════════════════════════════════════════════════════════════
#  FEATURE OVERLAP ANALYSIS
#  Compare features selected by DIABLO vs Random Forest
# ════════════════════════════════════════════════════════════════════════════

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║  FEATURE OVERLAP ANALYSIS - DIABLO vs Random Forest            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Load DIABLO selected features
diablo_mrna <- read.csv("DIABLO_Final_Publication_Fixed/Selected_mRNA.csv")$Feature
diablo_meth <- read.csv("DIABLO_Final_Publication_Fixed/Selected_Methylation.csv")$Feature

# Load RF variable importance (top 50)
rf_importance <- read.csv("Random_Forest_Concatenated_Results/Variable_Importance.csv")
rf_top50 <- rf_importance[1:50, ]

# Extract feature names (remove mRNA_ and Meth_ prefixes)
rf_mrna <- gsub("mRNA_", "", rf_top50$Feature[grepl("mRNA_", rf_top50$Feature)])
rf_meth <- gsub("Meth_", "", rf_top50$Feature[grepl("Meth_", rf_top50$Feature)])

# Calculate overlaps
overlap_mrna <- intersect(diablo_mrna, rf_mrna)
overlap_meth <- intersect(diablo_meth, rf_meth)

# Print results
cat("═══════════════════════════════════════════════════════════════\n")
cat("mRNA FEATURES\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat(sprintf("DIABLO selected: %d mRNA features\n", length(diablo_mrna)))
cat(sprintf("RF top 50 contains: %d mRNA features\n", length(rf_mrna)))
cat(sprintf("Overlap: %d features (%.1f%% of DIABLO)\n\n", 
            length(overlap_mrna), 
            100*length(overlap_mrna)/length(diablo_mrna)))

cat("Overlapping mRNA features:\n")
if (length(overlap_mrna) > 0) {
  for (feat in overlap_mrna) {
    rf_rank <- which(rf_importance$Feature == paste0("mRNA_", feat))
    rf_score <- rf_importance$Overall[rf_rank]
    cat(sprintf("  • %s (RF rank: %d, importance: %.1f)\n", feat, rf_rank, rf_score))
  }
} else {
  cat("  (none)\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("METHYLATION FEATURES\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat(sprintf("DIABLO selected: %d methylation features\n", length(diablo_meth)))
cat(sprintf("RF top 50 contains: %d methylation features\n", length(rf_meth)))
cat(sprintf("Overlap: %d features (%.1f%% of DIABLO)\n\n", 
            length(overlap_meth), 
            100*length(overlap_meth)/length(diablo_meth)))

cat("Overlapping methylation features:\n")
if (length(overlap_meth) > 0) {
  for (feat in overlap_meth) {
    rf_rank <- which(rf_importance$Feature == paste0("Meth_", feat))
    rf_score <- rf_importance$Overall[rf_rank]
    cat(sprintf("  • %s (RF rank: %d, importance: %.1f)\n", feat, rf_rank, rf_score))
  }
} else {
  cat("(none)\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

total_overlap <- length(overlap_mrna) + length(overlap_meth)
total_diablo <- length(diablo_mrna) + length(diablo_meth)

cat(sprintf("Total DIABLO features: %d\n", total_diablo))
cat(sprintf("Total overlap with RF top 50: %d (%.1f%%)\n\n", 
            total_overlap, 
            100*total_overlap/total_diablo))

if (total_overlap/total_diablo > 0.5) {
  cat("✅ STRONG AGREEMENT: >50% overlap indicates both methods\n")
  cat("   identify similar important features.\n\n")
} else if (total_overlap/total_diablo > 0.3) {
  cat("✅ MODERATE AGREEMENT: 30-50% overlap is reasonable.\n")
  cat("   Different methods capture different aspects.\n\n")
} else {
  cat("⚠️  LOW AGREEMENT: <30% overlap suggests methods capture\n")
  cat("   different feature sets (expected for sparse vs dense models).\n\n")
}

cat("INTERPRETATION:\n")
cat("───────────────────────────────────────────────────────────────\n")
cat("• DIABLO: Sparse model (41 features), optimized for correlation\n")
cat("• Random Forest: Dense model (400 features), optimized for accuracy\n")
cat("• Overlap indicates ROBUST features important across methods\n")
cat("• Non-overlap is EXPECTED and shows complementary information\n\n")

# Save overlap results
write.csv(data.frame(
  Feature = c(paste0("mRNA_", overlap_mrna), paste0("Meth_", overlap_meth)),
  Type = c(rep("mRNA", length(overlap_mrna)), rep("Methylation", length(overlap_meth)))
), "Feature_Overlap_DIABLO_RF.csv", row.names = FALSE)

cat("✓ Overlap results saved to Feature_Overlap_DIABLO_RF.csv\n\n")
