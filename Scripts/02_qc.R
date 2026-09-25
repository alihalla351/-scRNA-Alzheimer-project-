###############################################################################
# GSE141044 — IMPORT + QC
# APP23 mouse hippocampal single-nucleus RNA-seq
#
# This script:
# 1. Loads the raw expression table (mislabeled as .mtx, actually a
#    dense whitespace-delimited table)
# 2. Extracts nucleus IDs, gene names, and sample metadata (age/genotype)
# 3. Builds the Seurat object
# 4. Applies QC filtering (thresholds follow the original paper)
#
# Expects the raw data file at: data/GSE141044_matrix.mtx
# (relative to the repo's root folder — run this script from there)
###############################################################################

# ---- 0. Clean session ----
rm(list = ls())
gc()
set.seed(100)
options(stringsAsFactors = FALSE)

# ---- 1. Load packages ----
packages <- c("Matrix", "data.table", "Seurat")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(Matrix)
library(data.table)
library(Seurat)

# ---- 2. Locate the matrix file ----
matrix_file <- "data/GSE141044_matrix.mtx"

if (!file.exists(matrix_file)) {
  stop("Could not find data/GSE141044_matrix.mtx. Place the raw data file in the data/ folder.")
}

cat("\nMatrix file:\n")
print(matrix_file)

# ---- 3. Read first line separately (nucleus IDs only) ----
first_line <- readLines(matrix_file, n = 1)
cell_ids <- strsplit(first_line, "\\s+")[[1]]
cell_ids <- gsub('^"|"$', "", cell_ids)
cell_ids <- cell_ids[cell_ids != ""]

cat("\nNumber of nucleus IDs read:\n")
print(length(cell_ids))

# ---- 4. Read expression values (skip first line, no header) ----
expr_dt <- data.table::fread(
  matrix_file, skip = 1, header = FALSE, sep = " ",
  data.table = TRUE, check.names = FALSE, showProgress = TRUE
)

expected_columns <- length(cell_ids) + 1
if (ncol(expr_dt) != expected_columns) {
  stop("Column mismatch between nucleus IDs and expression data. Stop and check the file.")
}

# ---- 5. Extract and clean gene names ----
genes <- as.character(expr_dt[[1]])
genes <- trimws(gsub('^"|"$', "", genes))
genes[is.na(genes) | genes == ""] <- paste0("UnknownGene_", which(is.na(genes) | genes == ""))
genes <- make.unique(genes, sep = "_dup")

expr_dt[, 1 := NULL]

# ---- 6. Convert to matrix, assign names ----
counts_dense <- as.matrix(expr_dt)
storage.mode(counts_dense) <- "numeric"
rownames(counts_dense) <- genes
colnames(counts_dense) <- cell_ids

if (!all(colnames(counts_dense) == cell_ids)) {
  stop("Matrix cell IDs were not assigned correctly.")
}

# ---- 7. Convert to sparse matrix ----
counts <- Matrix::Matrix(counts_dense, sparse = TRUE)
counts <- as(counts, "dgCMatrix")
rm(counts_dense, expr_dt)
gc()

if (anyDuplicated(rownames(counts)) > 0) {
  rownames(counts) <- make.unique(rownames(counts), sep = "_dup")
}
if (anyDuplicated(colnames(counts)) > 0) {
  colnames(counts) <- make.unique(colnames(counts), sep = "_dup")
}

# ---- 8. Extract sample metadata from nucleus IDs ----
# Example ID: AACGTGAT_M6_APP1_L1
cell_ids <- colnames(counts)
sample_full <- sub("^[^_]+_", "", cell_ids)

age <- ifelse(grepl("^M6_", sample_full, ignore.case = TRUE), "6_month",
       ifelse(grepl("^M24_", sample_full, ignore.case = TRUE), "24_month", NA))

genotype <- ifelse(grepl("APP", sample_full, ignore.case = TRUE), "APP23",
            ifelse(grepl("WT", sample_full, ignore.case = TRUE), "WT", NA))

sample_id <- sub("_L[0-9]+$", "", sample_full)

if (anyNA(age)) stop("Some nuclei could not be assigned an age.")
if (anyNA(genotype)) stop("Some nuclei could not be assigned a genotype.")

# ---- 9. Create Seurat object ----
app23 <- CreateSeuratObject(counts = counts, project = "GSE141044_APP23",
                             min.cells = 0, min.features = 0)

app23$sample_full <- sample_full
app23$sample_id   <- sample_id
app23$age         <- factor(age, levels = c("6_month", "24_month"))
app23$genotype    <- factor(genotype, levels = c("WT", "APP23"))
app23$condition   <- factor(paste(app23$age, app23$genotype, sep = "_"),
                             levels = c("6_month_WT", "6_month_APP23",
                                        "24_month_WT", "24_month_APP23"))

# ---- 10. Remove known artifact feature "a" if present ----
if ("a" %in% rownames(app23)) {
  app23 <- subset(app23, features = setdiff(rownames(app23), "a"))
}

app23$orig.ident <- app23$sample_id

# ---- 11. QC metrics ----
app23[["percent.mt"]]   <- PercentageFeatureSet(app23, pattern = "^[mM][tT]-")
app23[["percent.ribo"]] <- PercentageFeatureSet(app23, pattern = "^[Rr][Pp][LlSs]")

cat("\nMitochondrial genes detected:\n")
print(sum(grepl("^[mM][tT]-", rownames(app23))))

cat("\nQC summary BEFORE filtering:\n")
print(summary(app23@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.ribo")]))

# ---- 12. QC filtering ----
# Thresholds follow the original paper: nFeature_RNA > 500, nCount_RNA > 4000
# percent.mt is NOT used as a filter — mitochondrial genes are absent
# from this processed matrix.
n_before <- ncol(app23)

app23 <- subset(app23, subset = nFeature_RNA > 500 & nCount_RNA > 4000)

n_after <- ncol(app23)

cat("\nNuclei before QC:", n_before, "\n")
cat("Nuclei after QC:", n_after, "\n")
cat("Nuclei removed:", n_before - n_after, "\n")
cat("Percentage retained:", round(100 * n_after / n_before, 2), "\n")

cat("\nNuclei per biological sample AFTER QC:\n")
print(table(app23$sample_id))

cat("\nAge x genotype AFTER QC:\n")
print(table(app23$age, app23$genotype))

# ---- 13. QC figures ----
dir.create("FIGURES/QC", recursive = TRUE, showWarnings = FALSE)

pdf("FIGURES/QC/QC_violin_plots.pdf", width = 11, height = 6)
print(VlnPlot(app23, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"),
              group.by = "sample_id", pt.size = 0, ncol = 3))
dev.off()

pdf("FIGURES/QC/QC_nCount_vs_nFeature.pdf", width = 7, height = 6)
print(FeatureScatter(app23, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

# ---- 14. Save post-QC object ----
dir.create("RESULTS/OBJECTS", recursive = TRUE, showWarnings = FALSE)
saveRDS(app23, "RESULTS/OBJECTS/GSE141044_post_QC.rds")

cat("\nImport + QC completed. Saved to RESULTS/OBJECTS/GSE141044_post_QC.rds\n")