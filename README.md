
### **PART 1: The Logical Strategy (What to compare)**

To publish your article, you must run these **three** new models to compare against your existing DIABLO result.

1.  **sPLS-DA on mRNA only:**
    *   *Purpose:* The "Control Group." Proves that adding Methylation actually improves the score.
    *   *Why sPLS-DA?* It is the exact same math as DIABLO, just with one dataset. This isolates "Data" as the only variable.
2.  **sPLS-DA on Methylation only:**
    *   *Purpose:* Proves that Methylation alone isn't enough.
3.  **Random Forest on Concatenated Data (mRNA + Meth):**
    *   *Purpose:* The "Dumb Integration" test. Proves that DIABLO’s "smart" correlation search is better than just pasting tables together.

---

### **PART 2: Technical "Must-Dos" Before Coding**

*   **Row Matching is Non-Negotiable:**
    *   Ensure `rownames(mRNA)` and `rownames(Methylation)` are identical.
    *   Ensure the `Y` vector (Cancer Type) matches the rows exactly.
*   **Handle Missing Data:**
    *   Remove any patient that has `NA` in *either* mRNA or Methylation before you start.
*   **Scaling for Concatenation:**
    *   For the Random Forest model, you **must** scale the data before combining.
    *   *Code:* `X_concat <- cbind(scale(mRNA), scale(Methylation))`
    *   *Reason:* mRNA values are huge (e.g., 0–10,000); Methylation is tiny (0–1). If not scaled, Random Forest will ignore Methylation.

---

### **PART 3: Coding the Models (Technical Checklist)**

#### **A. For the Single-Omics sPLS-DA Models (mRNA & Meth)**
*   **Tuning is Mandatory:** You tuned DIABLO, so you must tune these too.
    *   Use `tune.splsda()` to find the best `ncomp` and `keepX`.
    *   *Fairness:* If you use default parameters here but optimized parameters for DIABLO, the reviewer will reject the comparison.
*   **Code Setup:**
    *   Use `mixOmics` package.
    *   Syntax: `splsda(X = mRNA, Y = Y, ncomp = tuned_ncomp, keepX = tuned_keepX)`

#### **B. For the Concatenated Random Forest Model**
*   **Column Name Conflict:**
    *   Check if a gene and a methylation probe have the same name (e.g., "TP53").
    *   *Fix:* `colnames(Methylation) <- paste0("Meth_", colnames(Methylation))` before binding.
*   **Feature Selection:**
    *   Concatenating creates a huge matrix. It is fair to select only the top variance features (e.g., top 1000 mRNA + top 1000 Meth) if the model takes too long to run.

---

### **PART 4: The "Golden Rule" (Validation Fairness)**

**This is the single most important technical point.**

*   **Identical Cross-Validation (CV) Folds:**
    *   To compare F1 scores fairly, Model A and Model B must be tested on the **exact same patients**.
*   **How to achieve this in R:**
    *   **Method 1 (The Seed):** Run `set.seed(123)` immediately before running the validation function (`perf` or `randomForest`).
    *   **Method 2 (The Manual Folds - Better):** Create a list of folds using `caret::createFolds` and pass this exact list to every model's validation step (if the package allows) or write a `for` loop to enforce it.
*   **Consistency:**
    *   If DIABLO used **5-fold, 50 repeats**, the sPLS-DA models must also use **5-fold, 50 repeats**.

# COMPLETE SUMMARY: Technical Notes for Fair Comparison

## 🎯 **THE GOLDEN RULE (Non-Negotiable)**

```r
# Create CV folds ONCE, use for ALL 4 methods
set.seed(123)
folds <- createFolds(Y, k = 5, returnTrain = FALSE)
saveRDS(folds, "cv_folds.rds")
```

**Why:** Using different test samples = invalid comparison. This is the #1 reason papers get rejected.


# **CRITICAL CHECKLIST (Before Running Analysis)**

## **Setup Phase:**
- [ ] ✅ Create CV folds ONCE with `set.seed(123)`
- [ ] ✅ Save folds: `saveRDS(folds, "cv_folds.rds")`
- [ ] ✅ Check `rownames(mRNA) == rownames(methylation)`
- [ ] ✅ Check `length(Y) == nrow(mRNA)`
- [ ] ✅ Remove incomplete cases consistently
- [ ] ✅ Check class balance: `table(Y)`

## **For sPLS-DA Methods:**
- [ ] ✅ Run `tune.splsda()` first (don't use defaults)
- [ ] ✅ Use `nrepeat = 50` (match DIABLO)
- [ ] ✅ Use same `dist = 'max.dist'` as DIABLO
- [ ] ✅ Apply same preprocessing as DIABLO

## **For Random Forest Methods:**
- [ ] ✅ Use manual CV loop with saved folds
- [ ] ✅ Don't rely on OOB error alone
- [ ] ✅ Consider feature reduction (match DIABLO's feature count)
- [ ] ✅ Calculate F1 from confusion matrix

## **For Concatenation Method:**
- [ ] ✅ Scale each omics BEFORE combining
- [ ] ✅ Check for duplicate column names
- [ ] ✅ Reduce features (prevent p >> n problem)
- [ ] ✅ Use same CV folds

---

# **FINAL WORKFLOW (Step-by-Step)**

```r
# ==========================================
# STEP 1: SETUP (DO ONCE)
# ==========================================
library(mixOmics)
library(randomForest)
library(caret)
library(pROC)

set.seed(123)
folds <- createFolds(Y, k = 5, returnTrain = FALSE)
saveRDS(folds, "cv_folds.rds")

# Check alignment
stopifnot(rownames(mRNA) == rownames(methylation))

# ==========================================
# STEP 2: METHOD 1 - sPLS-DA mRNA
# ==========================================
tune.mRNA <- tune.splsda(mRNA, Y, ncomp = 5, 
                         validation = 'Mfold', 
                         folds = 5, nrepeat = 50,
                         test.keepX = c(5,10,15,20,30,50,100))

splsda.mRNA <- splsda(mRNA, Y, 
                      ncomp = tune.mRNA$choice.ncomp$ncomp,
                      keepX = tune.mRNA$choice.keepX)

perf.mRNA <- perf(splsda.mRNA, validation = "Mfold", 
                  folds = 5, nrepeat = 50)

f1.mRNA <- # Extract from perf.mRNA

# ==========================================
# STEP 3: METHOD 2 - sPLS-DA Methylation
# ==========================================
# Repeat Step 2 with methylation data

# ==========================================
# STEP 4: METHOD 3 - RF Single-omics
# ==========================================
folds <- readRDS("cv_folds.rds")
results.rf <- data.frame()

for(i in 1:length(folds)) {
  test.idx <- folds[[i]]
  train.idx <- setdiff(1:nrow(mRNA), test.idx)
  
  rf <- randomForest(mRNA[train.idx,], Y[train.idx], ntree=500)
  pred <- predict(rf, mRNA[test.idx,])
  
  f1 <- confusionMatrix(pred, Y[test.idx])$byClass["F1"]
  results.rf <- rbind(results.rf, data.frame(fold=i, f1=f1))
}

f1.rf <- mean(results.rf$f1)

# ==========================================
# STEP 5: METHOD 4 - RF Concatenated
# ==========================================
# Scale each omics
mRNA.scaled <- scale(mRNA)
meth.scaled <- scale(methylation)

# Reduce features
mRNA.top <- mRNA.scaled[, order(apply(mRNA.scaled,2,var), 
                                decreasing=TRUE)[1:500]]
meth.top <- meth.scaled[, order(apply(meth.scaled,2,var), 
                                decreasing=TRUE)[1:500]]

# Concatenate
concat <- cbind(mRNA.top, meth.top)
colnames(concat) <- make.unique(colnames(concat))

# Use same CV loop as Method 3
# Replace mRNA with concat

# ==========================================
# STEP 6: COMPARE RESULTS
# ==========================================
results.table <- data.frame(
  Method = c("sPLS-DA mRNA", "sPLS-DA Meth", 
             "RF mRNA", "RF Concat", "DIABLO"),
  F1 = c(f1.mRNA, f1.meth, f1.rf, f1.concat, f1.diablo),
  AUC = c(auc.mRNA, auc.meth, auc.rf, auc.concat, auc.diablo)
)

print(results.table)
```

---

# **TOP 5 MOST COMMON MISTAKES (AVOID THESE!)**

| Mistake | Why It's Wrong | Correct Approach |
|---------|----------------|------------------|
| **Different CV folds** | Comparing different test sets | Create folds ONCE, save, reuse |
| **No parameter tuning for sPLS-DA** | Unfair vs optimized DIABLO | Always run `tune.splsda()` |
| **Concatenating without scaling** | Magnitude bias | Scale EACH omics before `cbind()` |
| **Using all features for concatenation** | Overfitting (p >> n) | Reduce to top variable features |
| **Different preprocessing** | Inconsistent comparison | Apply SAME steps to all methods |

---

# **ONE-SENTENCE SUMMARY FOR EACH METHOD**

1. **sPLS-DA mRNA:** Tune parameters, use same CV as DIABLO, apply same preprocessing
2. **sPLS-DA Methylation:** Identical to sPLS-DA mRNA, just different data
3. **RF Single-omics:** Manual CV loop with saved folds, consider feature reduction
4. **RF Concatenated:** Scale each omics first, reduce features, use same CV folds

---

# **TIME ESTIMATE**

- Setup (CV folds, data checking): **15 minutes**
- sPLS-DA mRNA (tuning + eval): **30 minutes**
- sPLS-DA Methylation: **30 minutes**
- RF single-omics: **20 minutes**
- RF concatenated: **25 minutes**


---

# **My Actual Recommendation:**

## **For YOUR paper (applied cancer prediction), I recommend:**

### **Primary Comparison: sPLS-DA (same as DIABLO)**
**Why:** Clean, defensible, easy to explain

```r
# Fair comparison - same method family
splsda.mRNA <- splsda(mRNA, Y)
splsda.meth <- splsda(methylation, Y)
diablo.multi <- block.splsda(list(mRNA, meth), Y)
```

### **Secondary Comparison: Also show Random Forest**
**Why:** Shows you beat "real world" methods too

```r
# Practical comparison
rf.mRNA <- randomForest(mRNA, Y)
# Compare this to DIABLO too
```

---


## **Bottom Line:**

**Both AIs recommended sPLS-DA because it's the "cleanest" scientific comparison, but Random Forest is also perfectly valid for an applied paper.**

**The honest truth:** Most reviewers will accept either, as long as you:
1. Use the same CV splits
2. Use the same preprocessing  
3. Report statistical significance
4. Explain your choice

----

PART 1 — Core Strategy: What You Must Compare

You need three new baseline models to compare against your DIABLO results:
	1.	sPLS-DA on mRNA only
	•	Serves as the control model.
	•	Shows whether adding methylation data improves performance.
	•	Uses the same mathematical foundation as DIABLO, ensuring the method is not confounded.
	2.	sPLS-DA on Methylation only
	•	Tests whether methylation alone carries enough predictive signal.
	3.	Random Forest on Concatenated mRNA + Methylation
	•	Represents a naive integration method.
	•	Tests whether simply “gluing” omics datasets together can match DIABLO.

⸻

PART 2 — Technical Requirements Before Coding

Row Matching

Sample alignment is mandatory:
	•	rownames(mRNA) must equal rownames(Methylation)
	•	names(Y) must match the same samples
Any mismatch invalidates all results.

Missing Data Handling

Remove any patient with missing entries in either omics dataset.

Scaling Before Concatenation (for Random Forest)

Because mRNA and methylation are on different scales (large vs tiny values), scale each dataset before combining:

X_concat <- cbind(scale(mRNA), scale(Methylation))


⸻

PART 3 — Model Implementation: Technical Checklist

A. sPLS-DA (Single Omics)
	•	Must tune ncomp and keepX using tune.splsda()
	•	Must use same CV settings as DIABLO
	•	Must use the same preprocessing pipeline

Code example:

splsda(X = mRNA, Y = Y, ncomp = tuned_ncomp, keepX = tuned_keepX)

B. Random Forest (Concatenated Data)
	•	Add prefixes to avoid column name collisions (e.g., “Meth_”)
	•	Scale OMICS datasets before concatenation
	•	Reduce to top-variable features if needed (e.g., 1000 per omics)

⸻

PART 4 — The Golden Rule: FAIR Validation

To compare models validly:

Use Identical Cross-Validation Splits

This is the most important rule.

You can:
	1.	Use the same seed before running CV, or
	2.	Even better: pre-create CV folds using:

set.seed(123)
folds <- createFolds(Y, k = 5, returnTrain = FALSE)
saveRDS(folds, "cv_folds.rds")

All models (sPLS-DA, DIABLO, RF) must use exactly these folds.

CV Settings Must Match DIABLO

If DIABLO used 5-fold, 50 repeats, then all baselines must also use 5-fold, 50 repeats.

⸻

PART 5 — Final Results Table (To Put in Your Article)

You will ultimately fill a table like this:

Model Comparison	Algorithm	Data Used	F1 Score	Conclusion
Baseline 1	sPLS-DA	mRNA	0.xx	mRNA signal present
Baseline 2	sPLS-DA	Methylation	0.xx	Methylation weaker alone
Naive Integration	Random Forest	cbind(mRNA, Meth)	0.xx	Naive integration underperforms
Proposed	DIABLO	mRNA + Meth	0.xx	DIABLO integrates effectively


⸻

Setup Requirements (Essential Checklist)

Data Alignment
	•	Check row consistency
	•	Remove incomplete cases
	•	Check class balance

Consistent Preprocessing

Apply exactly the same preprocessing as DIABLO to all methods.

⸻

sPLS-DA (mRNA & Methylation)

Parameter Tuning

Must run tune.splsda() with:
	•	same number of CV folds
	•	same repeats
	•	same test.keepX grid

Evaluation

Use perf() with the same settings as DIABLO (M-fold, nrepeat, etc.)

⸻

Random Forest (Single-Omics and Concatenated)

Use Manual CV Loop

Do not use OOB error — it is not comparable to DIABLO.

Use the saved folds:

folds <- readRDS("cv_folds.rds")

Feature Reduction

Random Forest may require using top-variable features for fairness and speed.

⸻

Concatenated RF Method

Required:
	•	scale each omics dataset before combining
	•	avoid duplicate names
	•	reduce features
	•	use identical CV folds

⸻

Critical Checklist Before Running Analysis

ALL Must Be True:
	•	Same CV folds
	•	Same preprocessing
	•	Tune sPLS-DA
	•	Scale before concatenation
	•	Feature reduction for RF
	•	Use same number of repeats/folds

These ensure a fair comparison.



⸻
