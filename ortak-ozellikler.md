# Multi-Omics Analysis: Common Features and Methodological Differences

This document outlines the shared methodological foundations and specific characteristics of the four machine learning approaches used in this study to predict IDH mutation status.

## 1. Common Methodological Foundations (Shared by All Methods)

All four methods adhere to a strict, publication-grade statistical framework to ensure fair comparison, reproducibility, and prevention of overfitting.

### A. Data Integrity & Preprocessing
*   **Identical Dataset:** All methods use the exact same set of samples (n=561) derived from the intersection of mRNA, Methylation, and Clinical data.
*   **Sample Alignment:** Samples are rigorously aligned across all omics layers before any analysis begins.
*   **Reproducibility:** A fixed random seed (`set.seed(42)`) is used in every script to ensure that data splits and stochastic processes are 100% reproducible.

### B. Overfitting Control & Validation Strategy
*   **Stratified Train-Test Split (70/30):**
    *   The data is split into a **Training Set (70%, n=394)** and an independent **Test Set (30%, n=167)**.
    *   **Crucial:** This split is performed *before* any feature filtering, scaling, or model training to prevent "data leakage" (information from the test set influencing the model).
    *   **Stratification:** The split maintains the same ratio of Mutant/WT classes in both sets.
*   **Feature Pre-filtering (Train-Only):**
    *   Feature selection (e.g., keeping top 200 variable features) is calculated *only* using the Training data.
    *   The Test data is then subsetted to match these features, ensuring the test set remains truly "unseen."
*   **Cross-Validation (CV):**
    *   All models use **10-Fold Cross-Validation** repeated **5 times** on the training set for hyperparameter tuning.

### C. Performance Evaluation
*   **Metrics:** All methods are evaluated using the same set of robust metrics:
    *   **ROC-AUC:** Area Under the Receiver Operating Characteristic Curve (primary metric).
    *   **Balanced Accuracy:** Average of Sensitivity and Specificity (handles class imbalance).
    *   **Sensitivity & Specificity:** To assess false negatives and false positives.
    *   **F1-Score:** Harmonic mean of Precision and Recall.
*   **Permutation Testing:** Statistical significance of the models is validated using permutation tests (shuffling labels 200-1000 times) to generate empirical p-values.

---

## 2. Method-Specific Characteristics

### A. DIABLO (Multi-Omics Integration)
*   **Type:** Intermediate Fusion (Joint Latent Space Learning).
*   **Algorithm:** `block.splsda` (Sparse Partial Least Squares Discriminant Analysis).
*   **Key Feature:** Maximizes the covariance (correlation) between mRNA, Methylation, and the Outcome (Y) simultaneously.
*   **Design Matrix:** Uses a "design matrix" (set to 0.5) to balance the trade-off between maximizing correlation *between* omics layers vs. maximizing discrimination of the outcome.
*   **Output:** Identifies a "molecular signature" of highly correlated mRNA and Methylation features that jointly predict IDH status.

### B. Single-Omics mRNA
*   **Type:** Single View Learning.
*   **Algorithm:** `splsda` (Sparse PLS-DA).
*   **Key Feature:** Focuses solely on Transcriptomic data.
*   **Purpose:** Serves as a benchmark to see if adding methylation data (via DIABLO) improves performance over gene expression alone.
*   **Tuning:** Optimizes the number of mRNA features (`keepX`) to maximize classification accuracy.

### C. Single-Omics Methylation
*   **Type:** Single View Learning.
*   **Algorithm:** `splsda` (Sparse PLS-DA).
*   **Key Feature:** Focuses solely on Epigenetic (Methylation) data.
*   **Purpose:** Serves as a benchmark to evaluate the predictive power of CpG methylation sites independently.
*   **Tuning:** Optimizes the number of CpG probes (`keepX`) to maximize classification accuracy.

### D. Random Forest (Concatenated)
*   **Type:** Early Fusion (Concatenation).
*   **Algorithm:** `randomForest` (Ensemble of Decision Trees).
*   **Key Feature:** Concatenates mRNA and Methylation matrices into one large matrix ($X_{all} = [X_{mRNA}, X_{meth}]$) before training.
*   **Mechanism:** Builds hundreds of decision trees and aggregates their votes. Handles non-linear relationships and interactions between features well.
*   **Special Note:** Often achieves Train AUC = 1.0 (perfect fit) due to its nature, so "Out-of-Bag" (OOB) error is used as the internal validation metric instead of standard CV error.

---

## 3. Summary Comparison Table

| Feature | DIABLO | Single-Omics mRNA | Single-Omics Meth | Random Forest |
| :--- | :--- | :--- | :--- | :--- |
| **Integration Type** | Intermediate Fusion | None | None | Early Fusion |
| **Algorithm** | block.splsda | splsda | splsda | Random Forest |
| **Input Data** | List(mRNA, Meth) | Matrix(mRNA) | Matrix(Meth) | Matrix(mRNA + Meth) |
| **Feature Selection** | Built-in (L1 penalty) | Built-in (L1 penalty) | Built-in (L1 penalty) | Gini Importance |
| **Tuning Method** | 10-Fold CV (5x) | 10-Fold CV (5x) | 10-Fold CV (5x) | 10-Fold CV (5x) |
| **Linearity** | Linear Combinations | Linear Combinations | Linear Combinations | Non-Linear |
| **Goal** | Correlation + Prediction | Prediction | Prediction | Prediction |

---

## 4. Why This Setup is Robust
By using **identical data splits** and **evaluation metrics** across all four methods, we ensure that any difference in performance is due to the **methodology itself** (e.g., integration strategy, algorithm) and not due to random chance or data artifacts. This provides a rigorous scientific basis for comparing Multi-Omics integration against Single-Omics baselines.
