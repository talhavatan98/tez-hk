# METHOD COMPARISON SUMMARY

## 1. Overview
All methods performed exceptionally well, with AUCs > 0.94. This suggests that the IDH mutation signal is very strong in both mRNA and Methylation data individually.

## 2. Performance Metrics Comparison

| Metric | DIABLO (Multi-Omics) | Random Forest (Concat) | mRNA Only (sPLS-DA) | Methylation Only (sPLS-DA) |
| :--- | :--- | :--- | :--- | :--- |
| **ROC AUC** | **0.9485** | **0.9517** | 0.9484 | 0.9421 |
| **Balanced Accuracy** | 89.2% | **91.2%** | 89.2% | 90.6% |
| **Sensitivity** | 99.2% | 99.2% | 99.2% | **100.0%** |
| **Specificity** | 79.2% | **83.3%** | 79.2% | 81.3% |
| **F1-Score** | 0.9555 | 0.9633 | 0.9555 | **0.9636** |
| **MCC** | 0.8378 | 0.8674 | 0.8378 | **0.8691** |

## 3. Key Observations

### A. Random Forest (Concatenated) is the Top Performer
- **Highest ROC-AUC (0.952)** and **Specificity (83.3%)**.
- This is expected because Random Forest is a non-linear model that can capture complex interactions between mRNA and methylation that linear models (like sPLS-DA/DIABLO) might miss.
- **Note:** The difference is small (0.952 vs 0.949), but consistent.

### B. DIABLO vs. Single Omics
- **DIABLO (0.9485)** performed almost identically to **mRNA Only (0.9484)**.
- **Methylation Only (0.9421)** had slightly lower AUC but **perfect Sensitivity (100%)**.
- This suggests that for *classification performance alone*, adding Methylation to mRNA (via DIABLO) didn't drastically boost the AUC compared to mRNA alone.
- **However**, DIABLO's value isn't just classification; it's **biological insight** (finding correlated features across omics).

### C. Why is DIABLO not "much higher"?
1.  **Signal Saturation:** The mRNA signal for IDH mutation is already so strong (AUC ~0.95) that there's little room for improvement.
2.  **Linearity:** DIABLO is a linear method. If the biological relationship is non-linear, Random Forest will naturally win.
3.  **Feature Selection:** DIABLO was restricted to a small number of features (3-15) for interpretability, while RF used all 400 pre-filtered features.

## 4. Conclusion for Publication
- **Honesty is key:** Report these results exactly as they are.
- **Narrative:** "While Random Forest achieved marginally higher classification performance (AUC 0.952 vs 0.949), DIABLO provided a parsimonious, interpretable signature of correlated multi-omics features with comparable predictive power."
- **Validation:** The fact that all methods agree (AUC > 0.94) confirms the robustness of your findings.
