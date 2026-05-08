#!/usr/bin/env Rscript
# ════════════════════════════════════════════════════════════════════════════
#  METHYLATION PROBE TO GENE CONVERTER
#  Converts CpG probe IDs (cg...) to Gene Symbols for Enrichr Analysis
# ════════════════════════════════════════════════════════════════════════════

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║  METHYLATION PROBE -> GENE CONVERTER                           ║\n")
cat("║  Preparing list for Biological Validation (Enrichr)            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# 1. Install/Load Required Package
# ═══════════════════════════════════════════════════════════════════════════
cat("STEP 1: Loading Annotation Package\n")
cat("─────────────────────────────────────────────────────────────────\n")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Check if annotation package is installed
if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  cat("  • Installing IlluminaHumanMethylation450kanno.ilmn12.hg19 (this may take a moment)...\n")
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", update = FALSE, ask = FALSE)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
} else {
  cat("  ✓ Package already installed.\n")
}

# 2. Load Probe List
# ═══════════════════════════════════════════════════════════════════════════
cat("\nSTEP 2: Loading Probe List\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Read the "Flexible" list (contains all selected features)
input_file <- "Enrichr_Methylation_All_Flexible.txt"

if (!file.exists(input_file)) {
  stop(sprintf("Error: Input file '%s' not found. Run Feature_Overlap_Analysis_Comprehensive.R first.", input_file))
}

meth_probes <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)$V1
cat(sprintf("✓ Loaded %d CpG probes from %s\n", length(meth_probes), input_file))

# 3. Convert to Genes
# ═══════════════════════════════════════════════════════════════════════════
cat("\nSTEP 3: Converting Probes to Genes\n")
cat("─────────────────────────────────────────────────────────────────\n")

# Get annotation data
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Match probes
matched_genes <- ann[meth_probes, "UCSC_RefGene_Name"]

# Handle missing matches
missing_count <- sum(is.na(matched_genes) | matched_genes == "")
if (missing_count > 0) {
  cat(sprintf("  ⚠ %d probes did not map to any gene (intergenic regions)\n", missing_count))
}

# Split multiple genes (some probes map to multiple isoforms/genes, separated by ;)
meth_gene_list <- unique(unlist(strsplit(matched_genes, ";")))

# Clean empty strings and NAs
meth_gene_list <- meth_gene_list[!is.na(meth_gene_list) & meth_gene_list != ""]

# Sort alphabetically
meth_gene_list <- sort(meth_gene_list)

cat(sprintf("✓ Successfully mapped to %d unique genes\n", length(meth_gene_list)))

# 4. Save Results
# ═══════════════════════════════════════════════════════════════════════════
cat("\nSTEP 4: Saving Results\n")
cat("─────────────────────────────────────────────────────────────────\n")

output_file <- "Enrichr_Methylation_Genes_Ready.txt"
write.table(meth_gene_list, output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(sprintf("✓ Saved gene list to: %s\n\n", output_file))

cat("INSTRUCTIONS FOR ENRICHR:\n")
cat("1. Open 'Enrichr_Methylation_Genes_Ready.txt'\n")
cat("2. Copy the gene list\n")
cat("3. Open 'Enrichr_mRNA_All_Flexible.txt'\n")
cat("4. Copy that gene list too\n")
cat("5. Paste BOTH lists into https://maayanlab.cloud/Enrichr/\n")
cat("6. Analyze!\n\n")
