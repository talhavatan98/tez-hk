#!/usr/bin/env Rscript
# ════════════════════════════════════════════════════════════════════════════
#  COMPREHENSIVE FEATURE OVERLAP ANALYSIS
#  Compare features selected across ALL methods:
#  - DIABLO (multi-omics)
#  - Single-omics mRNA (sPLS-DA)
#  - Single-omics Methylation (sPLS-DA)
#  - Random Forest (concatenated)
# ════════════════════════════════════════════════════════════════════════════

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║  COMPREHENSIVE FEATURE OVERLAP ANALYSIS                        ║\n")
cat("║  All Methods: DIABLO, Single-Omics, Random Forest              ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# 1. LOAD ALL SELECTED FEATURES
# ═══════════════════════════════════════════════════════════════════════════

cat("STEP 1: Loading Selected Features from All Methods\n")
cat("─────────────────────────────────────────────────────────────────\n")

# DIABLO
diablo_mrna <- read.csv("DIABLO_Final_Publication_Fixed/Selected_mRNA.csv")$Feature
diablo_meth <- read.csv("DIABLO_Final_Publication_Fixed/Selected_Methylation.csv")$Feature

# DIABLO Feature Importance (NEW)
diablo_mrna_imp <- read.csv("DIABLO_Final_Publication_Fixed/Feature_Importance_mRNA.csv")
diablo_meth_imp <- read.csv("DIABLO_Final_Publication_Fixed/Feature_Importance_Methylation.csv")

# Single-omics
single_mrna <- read.csv("Single_Omics_mRNA_Results/Selected_Features.csv")$Feature
single_meth <- read.csv("Single_Omics_Methylation_Results/Selected_Features.csv")$Feature

# Single-omics Feature Importance (NEW)
single_mrna_imp <- read.csv("Single_Omics_mRNA_Results/Feature_Importance.csv")
single_meth_imp <- read.csv("Single_Omics_Methylation_Results/Feature_Importance.csv")

# Random Forest (top 50)
rf_importance <- read.csv("Random_Forest_Concatenated_Results/Variable_Importance.csv")
rf_top50 <- rf_importance[1:50, ]
rf_mrna <- gsub("mRNA_", "", rf_top50$Feature[grepl("mRNA_", rf_top50$Feature)])
rf_meth <- gsub("Meth_", "", rf_top50$Feature[grepl("Meth_", rf_top50$Feature)])

cat(sprintf("✓ DIABLO: %d mRNA + %d Methylation = %d total\n", 
            length(diablo_mrna), length(diablo_meth), 
            length(diablo_mrna) + length(diablo_meth)))
cat(sprintf("✓ Single-omics mRNA: %d features\n", length(single_mrna)))
cat(sprintf("✓ Single-omics Methylation: %d features\n", length(single_meth)))
cat(sprintf("✓ Random Forest (top 50): %d mRNA + %d Methylation\n\n", 
            length(rf_mrna), length(rf_meth)))

# ═══════════════════════════════════════════════════════════════════════════
# 2. mRNA FEATURE OVERLAP ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════\n")
cat("mRNA FEATURES - PAIRWISE OVERLAP\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Pairwise overlaps
overlap_diablo_single_mrna <- intersect(diablo_mrna, single_mrna)
overlap_diablo_rf_mrna <- intersect(diablo_mrna, rf_mrna)
overlap_single_rf_mrna <- intersect(single_mrna, rf_mrna)

# Three-way overlap
overlap_all_mrna <- Reduce(intersect, list(diablo_mrna, single_mrna, rf_mrna))

cat("Pairwise Overlaps:\n")
cat(sprintf("  • DIABLO ∩ Single-omics: %d features (%.1f%% of DIABLO)\n", 
            length(overlap_diablo_single_mrna),
            100*length(overlap_diablo_single_mrna)/length(diablo_mrna)))
cat(sprintf("  • DIABLO ∩ RF top 50:    %d features (%.1f%% of DIABLO)\n", 
            length(overlap_diablo_rf_mrna),
            100*length(overlap_diablo_rf_mrna)/length(diablo_mrna)))
cat(sprintf("  • Single ∩ RF top 50:    %d features (%.1f%% of Single)\n\n", 
            length(overlap_single_rf_mrna),
            100*length(overlap_single_rf_mrna)/length(single_mrna)))

cat(sprintf("Three-way Overlap (DIABLO ∩ Single ∩ RF): %d features\n", 
            length(overlap_all_mrna)))

if (length(overlap_all_mrna) > 0) {
  cat("\n✅ ROBUST mRNA BIOMARKERS (validated across ALL methods):\n")
  for (feat in overlap_all_mrna) {
    rf_rank <- which(rf_importance$Feature == paste0("mRNA_", feat))
    rf_score <- rf_importance$Overall[rf_rank]
    cat(sprintf("  ★ %s (RF rank: %d, importance: %.1f)\n", feat, rf_rank, rf_score))
  }
} else {
  cat("  (No three-way overlap)\n")
}

cat("\nDIABLO ∩ Single-omics overlap:\n")
if (length(overlap_diablo_single_mrna) > 0) {
  for (feat in overlap_diablo_single_mrna) {
    # Get DIABLO loading
    diablo_loading <- diablo_mrna_imp$Max_Abs_Loading[diablo_mrna_imp$Feature == feat]
    if (length(diablo_loading) == 0) diablo_loading <- NA
    
    # Get Single-omics loading
    single_loading <- single_mrna_imp$Max_Abs_Loading[single_mrna_imp$Feature == feat]
    if (length(single_loading) == 0) single_loading <- NA
    
    cat(sprintf("  • %s (DIABLO: %.3f, Single: %.3f)\n", 
                feat, 
                ifelse(is.na(diablo_loading), 0, diablo_loading),
                ifelse(is.na(single_loading), 0, single_loading)))
  }
} else {
  cat("  (none)\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# 3. METHYLATION FEATURE OVERLAP ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("METHYLATION FEATURES - PAIRWISE OVERLAP\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Pairwise overlaps
overlap_diablo_single_meth <- intersect(diablo_meth, single_meth)
overlap_diablo_rf_meth <- intersect(diablo_meth, rf_meth)
overlap_single_rf_meth <- intersect(single_meth, rf_meth)

# Three-way overlap
overlap_all_meth <- Reduce(intersect, list(diablo_meth, single_meth, rf_meth))

cat("Pairwise Overlaps:\n")
cat(sprintf("  • DIABLO ∩ Single-omics: %d features (%.1f%% of DIABLO)\n", 
            length(overlap_diablo_single_meth),
            100*length(overlap_diablo_single_meth)/length(diablo_meth)))
cat(sprintf("  • DIABLO ∩ RF top 50:    %d features (%.1f%% of DIABLO)\n", 
            length(overlap_diablo_rf_meth),
            100*length(overlap_diablo_rf_meth)/length(diablo_meth)))
cat(sprintf("  • Single ∩ RF top 50:    %d features (%.1f%% of Single)\n\n", 
            length(overlap_single_rf_meth),
            100*length(overlap_single_rf_meth)/length(single_meth)))

cat(sprintf("Three-way Overlap (DIABLO ∩ Single ∩ RF): %d features\n", 
            length(overlap_all_meth)))

if (length(overlap_all_meth) > 0) {
  cat("\n✅ ROBUST METHYLATION BIOMARKERS (validated across ALL methods):\n")
  for (feat in overlap_all_meth) {
    rf_rank <- which(rf_importance$Feature == paste0("Meth_", feat))
    rf_score <- rf_importance$Overall[rf_rank]
    cat(sprintf("  ★ %s (RF rank: %d, importance: %.1f)\n", feat, rf_rank, rf_score))
  }
} else {
  cat("  (No three-way overlap)\n")
}

cat("\nDIABLO ∩ Single-omics overlap:\n")
if (length(overlap_diablo_single_meth) > 0) {
  for (feat in overlap_diablo_single_meth) {
    # Get DIABLO loading
    diablo_loading <- diablo_meth_imp$Max_Abs_Loading[diablo_meth_imp$Feature == feat]
    if (length(diablo_loading) == 0) diablo_loading <- NA
    
    # Get Single-omics loading
    single_loading <- single_meth_imp$Max_Abs_Loading[single_meth_imp$Feature == feat]
    if (length(single_loading) == 0) single_loading <- NA
    
    cat(sprintf("  • %s (DIABLO: %.3f, Single: %.3f)\n", 
                feat, 
                ifelse(is.na(diablo_loading), 0, diablo_loading),
                ifelse(is.na(single_loading), 0, single_loading)))
  }
} else {
  cat("(none)\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# 4. SUMMARY STATISTICS
# ═══════════════════════════════════════════════════════════════════════════

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("OVERALL SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

total_diablo <- length(diablo_mrna) + length(diablo_meth)
total_single <- length(single_mrna) + length(single_meth)
total_overlap_diablo_single <- length(overlap_diablo_single_mrna) + length(overlap_diablo_single_meth)
total_overlap_all <- length(overlap_all_mrna) + length(overlap_all_meth)

cat("Feature Counts:\n")
cat(sprintf("  • DIABLO total:        %d features\n", total_diablo))
cat(sprintf("  • Single-omics total:  %d features\n", total_single))
cat(sprintf("  • RF top 50:           50 features\n\n"))

cat("Cross-Method Validation:\n")
cat(sprintf("  • DIABLO ∩ Single:     %d features (%.1f%% of DIABLO)\n", 
            total_overlap_diablo_single,
            100*total_overlap_diablo_single/total_diablo))
cat(sprintf("  • Three-way overlap:   %d features (ROBUST BIOMARKERS)\n\n", 
            total_overlap_all))

# Interpretation
if (total_overlap_diablo_single/total_diablo > 0.5) {
  cat("✅ STRONG AGREEMENT: >50% overlap between DIABLO and single-omics\n")
  cat("   → Both methods identify similar important features\n\n")
} else if (total_overlap_diablo_single/total_diablo > 0.3) {
  cat("✅ MODERATE AGREEMENT: 30-50% overlap is reasonable\n")
  cat("   → Different methods capture different aspects\n\n")
} else {
  cat("⚠️  LOW AGREEMENT: <30% overlap\n")
  cat("   → Expected for sparse (DIABLO) vs dense (single-omics) models\n\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# 5. SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════\n")
cat("SAVING RESULTS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Save all overlaps
if (length(overlap_all_mrna) > 0 || length(overlap_all_meth) > 0) {
  all_features <- c()
  all_types <- c()
  
  if (length(overlap_all_mrna) > 0) {
    all_features <- c(all_features, paste0("mRNA_", overlap_all_mrna))
    all_types <- c(all_types, rep("mRNA", length(overlap_all_mrna)))
  }
  
  if (length(overlap_all_meth) > 0) {
    all_features <- c(all_features, paste0("Meth_", overlap_all_meth))
    all_types <- c(all_types, rep("Methylation", length(overlap_all_meth)))
  }
  
  write.csv(data.frame(
    Feature = all_features,
    Type = all_types,
    Validated_By = rep("DIABLO + Single-omics + RF", length(all_features))
  ), "Feature_Overlap_All_Methods_ThreeWay.csv", row.names = FALSE)
} else {
  write.csv(data.frame(
    Feature = character(0),
    Type = character(0),
    Validated_By = character(0)
  ), "Feature_Overlap_All_Methods_ThreeWay.csv", row.names = FALSE)
}

write.csv(data.frame(
  Feature = c(paste0("mRNA_", overlap_diablo_single_mrna), 
              paste0("Meth_", overlap_diablo_single_meth)),
  Type = c(rep("mRNA", length(overlap_diablo_single_mrna)), 
           rep("Methylation", length(overlap_diablo_single_meth))),
  Validated_By = "DIABLO + Single-omics"
), "Feature_Overlap_DIABLO_SingleOmics.csv", row.names = FALSE)

# Create summary table
summary_table <- data.frame(
  Method = c("DIABLO", "Single-omics mRNA", "Single-omics Methylation", "Random Forest (top 50)"),
  mRNA_Features = c(length(diablo_mrna), length(single_mrna), 0, length(rf_mrna)),
  Meth_Features = c(length(diablo_meth), 0, length(single_meth), length(rf_meth)),
  Total_Features = c(total_diablo, length(single_mrna), length(single_meth), 50)
)

write.csv(summary_table, "Feature_Selection_Summary.csv", row.names = FALSE)

cat("✓ Feature_Overlap_All_Methods_ThreeWay.csv (robust biomarkers)\n")
cat("✓ Feature_Overlap_DIABLO_SingleOmics.csv (DIABLO vs single-omics)\n")
cat("✓ Feature_Selection_Summary.csv (method comparison)\n\n")

# ═══════════════════════════════════════════════════════════════════════════
# 6. VISUALIZATION DATA
# ═══════════════════════════════════════════════════════════════════════════

cat("═══════════════════════════════════════════════════════════════\n")
cat("VENN DIAGRAM DATA\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("For mRNA features:\n")
cat(sprintf("  • DIABLO only:           %d\n", 
            length(setdiff(diablo_mrna, union(single_mrna, rf_mrna)))))
cat(sprintf("  • Single-omics only:     %d\n", 
            length(setdiff(single_mrna, union(diablo_mrna, rf_mrna)))))
cat(sprintf("  • RF only:               %d\n", 
            length(setdiff(rf_mrna, union(diablo_mrna, single_mrna)))))
cat(sprintf("  • DIABLO ∩ Single:       %d\n", 
            length(setdiff(overlap_diablo_single_mrna, rf_mrna))))
cat(sprintf("  • DIABLO ∩ RF:           %d\n", 
            length(setdiff(overlap_diablo_rf_mrna, single_mrna))))
cat(sprintf("  • Single ∩ RF:           %d\n", 
            length(setdiff(overlap_single_rf_mrna, diablo_mrna))))
cat(sprintf("  • All three:             %d\n\n", length(overlap_all_mrna)))

cat("For Methylation features:\n")
cat(sprintf("  • DIABLO only:           %d\n", 
            length(setdiff(diablo_meth, union(single_meth, rf_meth)))))
cat(sprintf("  • Single-omics only:     %d\n", 
            length(setdiff(single_meth, union(diablo_meth, rf_meth)))))
cat(sprintf("  • RF only:               %d\n", 
            length(setdiff(rf_meth, union(diablo_meth, single_meth)))))
cat(sprintf("  • DIABLO ∩ Single:       %d\n", 
            length(setdiff(overlap_diablo_single_meth, rf_meth))))
cat(sprintf("  • DIABLO ∩ RF:           %d\n", 
            length(setdiff(overlap_diablo_rf_meth, single_meth))))
cat(sprintf("  • Single ∩ RF:           %d\n", 
            length(setdiff(overlap_single_rf_meth, diablo_meth))))
cat(sprintf("  • All three:             %d\n\n", length(overlap_all_meth)))

# ═══════════════════════════════════════════════════════════════════════════
# 7. BIOLOGICAL VALIDATION PREPARATION (ENRICHR READY)
# ═══════════════════════════════════════════════════════════════════════════

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("BIOLOGICAL VALIDATION PREPARATION\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Function to create consensus master list
create_consensus_list <- function(diablo_feats, single_feats, rf_feats, 
                                  diablo_imp, single_imp, rf_imp, type) {
  
  all_feats <- unique(c(diablo_feats, single_feats, rf_feats))
  
  consensus_df <- data.frame(
    Feature = all_feats,
    Type = type,
    Found_By_DIABLO = all_feats %in% diablo_feats,
    Found_By_Single = all_feats %in% single_feats,
    Found_By_RF = all_feats %in% rf_feats,
    DIABLO_Load = NA,
    Single_Load = NA,
    RF_Imp = NA,
    stringsAsFactors = FALSE
  )
  
  # Fill importance scores
  for(i in 1:nrow(consensus_df)) {
    f <- consensus_df$Feature[i]
    
    # DIABLO
    if(consensus_df$Found_By_DIABLO[i]) {
      val <- diablo_imp$Max_Abs_Loading[diablo_imp$Feature == f]
      if(length(val)>0) consensus_df$DIABLO_Load[i] <- val
    }
    
    # Single
    if(consensus_df$Found_By_Single[i]) {
      val <- single_imp$Max_Abs_Loading[single_imp$Feature == f]
      if(length(val)>0) consensus_df$Single_Load[i] <- val
    }
    
    # RF
    if(consensus_df$Found_By_RF[i]) {
      # RF names have prefix
      rf_name <- paste0(ifelse(type=="mRNA", "mRNA_", "Meth_"), f)
      val <- rf_imp$Overall[rf_imp$Feature == rf_name]
      if(length(val)>0) consensus_df$RF_Imp[i] <- val
    }
  }
  
  # Calculate Consensus Score (How many methods found it?)
  consensus_df$Method_Count <- rowSums(consensus_df[,c("Found_By_DIABLO", "Found_By_Single", "Found_By_RF")])
  
  # Assign Tiers
  # Tier 1 (Gold): Found by 2+ methods OR (1 method + very high importance)
  # Tier 2 (Silver): Found by 1 method
  consensus_df$Tier <- ifelse(consensus_df$Method_Count >= 2, "Tier 1 (Robust)", "Tier 2 (Method-Specific)")
  
  # Sort by Method Count then by max importance (normalized roughly)
  # Simple sorting: Method Count descending
  consensus_df <- consensus_df[order(consensus_df$Method_Count, decreasing = TRUE), ]
  
  return(consensus_df)
}

# Create Master Lists
mrna_consensus <- create_consensus_list(diablo_mrna, single_mrna, rf_mrna,
                                        diablo_mrna_imp, single_mrna_imp, rf_importance, "mRNA")

meth_consensus <- create_consensus_list(diablo_meth, single_meth, rf_meth,
                                        diablo_meth_imp, single_meth_imp, rf_importance, "Methylation")

# Save Master Tables
write.csv(mrna_consensus, "Biological_Validation_mRNA_MasterList.csv", row.names = FALSE)
write.csv(meth_consensus, "Biological_Validation_Methylation_MasterList.csv", row.names = FALSE)

# Export Clean Lists for Enrichr (Text Files)
# 1. mRNA - Tier 1 (Robust)
tier1_mrna <- mrna_consensus$Feature[mrna_consensus$Tier == "Tier 1 (Robust)"]
write.table(tier1_mrna, "Enrichr_mRNA_Tier1_Robust.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# 2. mRNA - All (Flexible)
write.table(mrna_consensus$Feature, "Enrichr_mRNA_All_Flexible.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# 3. Methylation - Tier 1 (Robust)
tier1_meth <- meth_consensus$Feature[meth_consensus$Tier == "Tier 1 (Robust)"]
write.table(tier1_meth, "Enrichr_Methylation_Tier1_Robust.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# 4. Methylation - All (Flexible)
write.table(meth_consensus$Feature, "Enrichr_Methylation_All_Flexible.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("✓ Generated Master Lists (CSV) with consensus scoring\n")
cat("✓ Generated Enrichr-ready text files (.txt):\n")
cat(sprintf("  • mRNA Tier 1 (Robust): %d genes\n", length(tier1_mrna)))
cat(sprintf("  • mRNA All (Flexible):  %d genes\n", nrow(mrna_consensus)))
cat(sprintf("  • Meth Tier 1 (Robust): %d CpGs\n", length(tier1_meth)))
cat(sprintf("  • Meth All (Flexible):  %d CpGs\n\n", nrow(meth_consensus)))

cat("═══════════════════════════════════════════════════════════════\n")
cat("ANALYSIS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("KEY FINDINGS:\n")
cat(sprintf("  • %d ROBUST mRNA biomarkers (validated across all methods)\n", 
            length(overlap_all_mrna)))
cat(sprintf("  • %d ROBUST methylation biomarkers (validated across all methods)\n", 
            length(overlap_all_meth)))
cat(sprintf("  • %.1f%% agreement between DIABLO and single-omics methods\n", 
            100*total_overlap_diablo_single/total_diablo))
cat("  • Enrichr lists are ready for biological validation step.\n\n")
