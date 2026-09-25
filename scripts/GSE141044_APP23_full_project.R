###############################################################################
# GSE141044 — CORRECT DATA IMPORT
# APP23 mouse hippocampus single-nucleus RNA-seq
#
# IMPORTANT:
# Despite the .mtx filename, the GEO matrix is NOT MatrixMarket format.
# It is a dense whitespace-delimited expression table:
#
# ROW 1    = nucleus/cell IDs
# COLUMN 1 = gene names
# VALUES   = expression counts
###############################################################################


# =============================================================================
# 0. CLEAN SESSION
# =============================================================================

rm(list = ls())
gc()

set.seed(100)

options(stringsAsFactors = FALSE)


# =============================================================================
# 1. LOAD PACKAGES
# =============================================================================

packages <- c(
  "Matrix",
  "data.table",
  "Seurat"
)

for (pkg in packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(pkg)
    
  }
  
}

library(Matrix)
library(data.table)
library(Seurat)


# =============================================================================
# 2. PROJECT DIRECTORY
# =============================================================================

project_dir <- "D:/Master/ScRNA seq/Course project"

setwd(project_dir)

cat("\n============================================================\n")
cat("GSE141044 APP23 snRNA-seq PROJECT\n")
cat("============================================================\n")

cat("\nWorking directory:\n")

print(getwd())


# =============================================================================
# 3. FIND THE EXTRACTED MATRIX FILE
# =============================================================================

all_files <- list.files(
  project_dir,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)


matrix_candidates <- all_files[
  grepl(
    "GSE141044_matrix\\.mtx$",
    all_files,
    ignore.case = TRUE
  )
]


cat("\n============================================================\n")
cat("MATRIX FILES FOUND\n")
cat("============================================================\n")

print(matrix_candidates)


if (length(matrix_candidates) == 0) {
  
  stop(
    "Could not find extracted GSE141044_matrix.mtx"
  )
  
}


matrix_file <- matrix_candidates[1]


cat("\nMatrix selected:\n")

print(matrix_file)


# =============================================================================
# 4. CONFIRM THE FORMAT
# =============================================================================

cat("\n============================================================\n")
cat("CHECKING ACTUAL FILE FORMAT\n")
cat("============================================================\n")


first_lines <- readLines(
  matrix_file,
  n = 3
)


cat("\nBeginning of file:\n")

print(
  substr(
    first_lines,
    1,
    300
  )
)


if (
  grepl(
    "^%%MatrixMarket",
    first_lines[1]
  )
) {
  
  stop(
    paste0(
      "This copy unexpectedly appears to be MatrixMarket format.\n",
      "Do not continue with this importer."
    )
  )
  
}


cat(
  "\nThis is a dense whitespace-delimited expression table.\n"
)

cat(
  "It will be imported as a table, NOT with readMM().\n"
)


# =============================================================================
# 5. READ THE EXPRESSION TABLE
#
# fread is much faster than base read.table for a large file.
# =============================================================================

cat("\n============================================================\n")
cat("READING EXPRESSION MATRIX\n")
cat("THIS MAY TAKE SOME TIME\n")
cat("============================================================\n")


expr_dt <- data.table::fread(
  matrix_file,
  header = TRUE,
  sep = " ",
  data.table = TRUE,
  check.names = FALSE,
  showProgress = TRUE
)


cat("\nExpression table successfully imported.\n")


cat("\nRaw dimensions:\n")

print(
  dim(expr_dt)
)


cat("\nFirst 10 column names:\n")

print(
  head(
    colnames(expr_dt),
    10
  )
)


cat("\nFirst 10 values in first column:\n")

print(
  head(
    expr_dt[[1]],
    10
  )
)


# =============================================================================
# 6. EXTRACT GENE NAMES
# =============================================================================

cat("\n============================================================\n")
cat("EXTRACTING GENE NAMES\n")
cat("============================================================\n")


genes <- as.character(
  expr_dt[[1]]
)


cat("\nFirst 20 genes:\n")

print(
  head(
    genes,
    20
  )
)


cat(
  "\nNumber of genes:",
  length(genes),
  "\n"
)


# =============================================================================
# 7. REMOVE GENE COLUMN FROM EXPRESSION DATA
# =============================================================================

expr_dt[, 1 := NULL]


cat("\nExpression-only dimensions:\n")

print(
  dim(expr_dt)
)


# =============================================================================
# 8. EXTRACT NUCLEUS IDS
# =============================================================================

cell_ids <- colnames(
  expr_dt
)


cat("\n============================================================\n")
cat("NUCLEUS IDS\n")
cat("============================================================\n")


cat(
  "\nNumber of nuclei:",
  length(cell_ids),
  "\n"
)


cat("\nFirst 20 nucleus IDs:\n")

print(
  head(
    cell_ids,
    20
  )
)


# =============================================================================
# 9. CLEAN GENE NAMES
# =============================================================================

bad_genes <- (
  is.na(genes) |
    genes == ""
)


if (
  any(bad_genes)
) {
  
  genes[bad_genes] <- paste0(
    "UnknownGene_",
    which(bad_genes)
  )
  
}


genes <- make.unique(
  genes
)


# =============================================================================
# 10. CONVERT TO NUMERIC MATRIX
# =============================================================================

cat("\n============================================================\n")
cat("CONVERTING EXPRESSION TABLE TO MATRIX\n")
cat("============================================================\n")


counts_dense <- as.matrix(
  expr_dt
)


storage.mode(
  counts_dense
) <- "numeric"


rownames(
  counts_dense
) <- genes


colnames(
  counts_dense
) <- cell_ids


cat("\nDense matrix dimensions:\n")

print(
  dim(counts_dense)
)


# =============================================================================
# 11. CONVERT TO SPARSE MATRIX
#
# Seurat should use a sparse matrix to reduce memory usage.
# =============================================================================

cat("\n============================================================\n")
cat("CONVERTING TO SPARSE MATRIX\n")
cat("============================================================\n")


counts <- Matrix::Matrix(
  counts_dense,
  sparse = TRUE
)


counts <- as(
  counts,
  "dgCMatrix"
)


cat("\nSparse matrix created.\n")


cat("\nDimensions:\n")

print(
  dim(counts)
)


cat("\nMatrix class:\n")

print(
  class(counts)
)


# Free memory

rm(
  counts_dense,
  expr_dt
)

gc()


# =============================================================================
# 12. SANITY CHECKS
# =============================================================================

cat("\n============================================================\n")
cat("SANITY CHECKS\n")
cat("============================================================\n")


cat(
  "\nGenes:",
  nrow(counts),
  "\n"
)


cat(
  "Nuclei:",
  ncol(counts),
  "\n"
)


cat(
  "Minimum count:",
  min(counts),
  "\n"
)


cat(
  "Maximum count:",
  max(counts),
  "\n"
)


cat(
  "Non-zero measurements:",
  length(counts@x),
  "\n"
)


if (
  any(counts@x < 0)
) {
  
  stop(
    "Negative values detected. This does not look like a raw count matrix."
  )
  
}


if (
  anyNA(counts@x)
) {
  
  stop(
    "NA values detected in expression matrix."
  )
  
}


cat("\nFirst 20 genes:\n")

print(
  head(
    rownames(counts),
    20
  )
)


cat("\nFirst 20 nuclei:\n")

print(
  head(
    colnames(counts),
    20
  )
)


# =============================================================================
# 13. EXTRACT SAMPLE INFORMATION FROM CELL IDS
#
# Example:
#
# AACGTGAT_M6_APP1_L1
#
# barcode = AACGTGAT
# sample  = M6_APP1_L1
# =============================================================================

cat("\n============================================================\n")
cat("EXTRACTING SAMPLE INFORMATION FROM NUCLEUS IDS\n")
cat("============================================================\n")


sample_id <- sub(
  "^[^_]+_",
  "",
  cell_ids
)


cat("\nUnique sample IDs detected:\n")

print(
  sort(
    unique(sample_id)
  )
)


cat(
  "\nNumber of unique sample IDs:",
  length(
    unique(sample_id)
  ),
  "\n"
)


cat("\nNumber of nuclei per sample:\n")

print(
  sort(
    table(sample_id),
    decreasing = TRUE
  )
)


# =============================================================================
# 14. EXTRACT BASIC GENOTYPE INFORMATION
# =============================================================================

genotype <- ifelse(
  grepl(
    "APP",
    sample_id,
    ignore.case = TRUE
  ),
  "APP23",
  ifelse(
    grepl(
      "WT",
      sample_id,
      ignore.case = TRUE
    ),
    "WT",
    NA
  )
)


cat("\nGenotype distribution:\n")

print(
  table(
    genotype,
    useNA = "ifany"
  )
)

# =============================================================================
# FIX DUPLICATE GENE AND CELL NAMES BEFORE CREATING SEURAT OBJECT
# =============================================================================

cat("\n============================================================\n")
cat("CHECKING DUPLICATE GENE / NUCLEUS NAMES\n")
cat("============================================================\n")


# -----------------------------------------------------------------------------
# 1. Check duplicates BEFORE correction
# -----------------------------------------------------------------------------

cat("\nDuplicated gene names BEFORE correction:\n")
print(sum(duplicated(rownames(counts))))

cat("\nDuplicated nucleus names BEFORE correction:\n")
print(sum(duplicated(colnames(counts))))


# Show examples if duplicates exist

dup_genes <- unique(
  rownames(counts)[
    duplicated(rownames(counts))
  ]
)

dup_cells <- unique(
  colnames(counts)[
    duplicated(colnames(counts))
  ]
)


cat("\nExamples of duplicated genes:\n")
print(head(dup_genes, 20))

cat("\nExamples of duplicated nuclei:\n")
print(head(dup_cells, 20))


# -----------------------------------------------------------------------------
# 2. Clean quotes / whitespace from names
# -----------------------------------------------------------------------------

clean_gene_names <- trimws(
  gsub(
    '^"|"$',
    "",
    rownames(counts)
  )
)

clean_cell_names <- trimws(
  gsub(
    '^"|"$',
    "",
    colnames(counts)
  )
)


# -----------------------------------------------------------------------------
# 3. Deal with missing/blank gene names
# -----------------------------------------------------------------------------

bad_gene_names <- (
  is.na(clean_gene_names) |
    clean_gene_names == ""
)

if (any(bad_gene_names)) {
  
  clean_gene_names[bad_gene_names] <- paste0(
    "UnknownGene_",
    which(bad_gene_names)
  )
}


# -----------------------------------------------------------------------------
# 4. Deal with missing/blank nucleus names
# -----------------------------------------------------------------------------

bad_cell_names <- (
  is.na(clean_cell_names) |
    clean_cell_names == ""
)

if (any(bad_cell_names)) {
  
  clean_cell_names[bad_cell_names] <- paste0(
    "Nucleus_",
    which(bad_cell_names)
  )
}


# -----------------------------------------------------------------------------
# 5. FORCE UNIQUE NAMES
#
# Example:
# Apoe
# Apoe
#
# becomes:
# Apoe
# Apoe.1
#
# We are NOT deleting expression data.
# -----------------------------------------------------------------------------

clean_gene_names <- make.unique(
  clean_gene_names,
  sep = "_dup"
)

clean_cell_names <- make.unique(
  clean_cell_names,
  sep = "_dup"
)


# -----------------------------------------------------------------------------
# 6. Put corrected names back into matrix
# -----------------------------------------------------------------------------

rownames(counts) <- clean_gene_names
colnames(counts) <- clean_cell_names


# -----------------------------------------------------------------------------
# 7. Verify that duplicates are now ZERO
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("AFTER CORRECTION\n")
cat("============================================================\n")

cat("\nDuplicated genes:\n")
print(sum(duplicated(rownames(counts))))

cat("\nDuplicated nuclei:\n")
print(sum(duplicated(colnames(counts))))

cat("\nAny duplicated genes?\n")
print(anyDuplicated(rownames(counts)))

cat("\nAny duplicated nuclei?\n")
print(anyDuplicated(colnames(counts)))


# These MUST both equal zero.

stopifnot(
  anyDuplicated(rownames(counts)) == 0
)

stopifnot(
  anyDuplicated(colnames(counts)) == 0
)


# -----------------------------------------------------------------------------
# 8. Make sure matrix class is appropriate for Seurat
# -----------------------------------------------------------------------------

counts <- as(
  counts,
  "dgCMatrix"
)


cat("\nMatrix class:\n")
print(class(counts))

cat("\nMatrix dimensions:\n")
print(dim(counts))


# =============================================================================
# CREATE SEURAT OBJECT
# =============================================================================

cat("\n============================================================\n")
cat("CREATING SEURAT OBJECT\n")
cat("============================================================\n")


app23 <- CreateSeuratObject(
  counts = counts,
  project = "GSE141044_APP23",
  min.cells = 0,
  min.features = 0
)


cat("\nSUCCESS — Seurat object created.\n\n")

print(app23)

cat("\nGenes x nuclei:\n")
print(dim(app23))

# =============================================================================
# 15. CREATE SEURAT OBJECT
# =============================================================================

cat("\n============================================================\n")
cat("CREATING SEURAT OBJECT\n")
cat("============================================================\n")


app23 <- CreateSeuratObject(
  counts = counts,
  project = "GSE141044_APP23",
  min.cells = 0,
  min.features = 0
)


app23$sample_id <- sample_id

app23$genotype <- genotype


cat("\nSeurat object successfully created.\n\n")

print(app23)


# =============================================================================
# 16. EXAMINE METADATA
# =============================================================================

cat("\n============================================================\n")
cat("METADATA CHECK\n")
cat("============================================================\n")


print(
  head(
    app23@meta.data,
    20
  )
)


cat("\nNuclei per sample:\n")

print(
  table(
    app23$sample_id
  )
)


cat("\nNuclei per genotype:\n")

print(
  table(
    app23$genotype,
    useNA = "ifany"
  )
)


# =============================================================================
# 17. ADD QC METRICS
# =============================================================================

cat("\n============================================================\n")
cat("CALCULATING QC METRICS\n")
cat("============================================================\n")


app23[["percent.mt"]] <- PercentageFeatureSet(
  app23,
  pattern = "^mt-"
)


app23[["percent.ribo"]] <- PercentageFeatureSet(
  app23,
  pattern = "^Rp[sl]"
)


qc_summary <- summary(
  app23@meta.data[
    ,
    c(
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt",
      "percent.ribo"
    )
  ]
)


cat("\nQC summary:\n")

print(qc_summary)


# =============================================================================
# 18. CREATE OUTPUT DIRECTORIES
# =============================================================================

dir.create(
  "RESULTS",
  showWarnings = FALSE
)


dir.create(
  "RESULTS/OBJECTS",
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  "FIGURES",
  showWarnings = FALSE
)


# =============================================================================
# 19. SAVE SAMPLE INFORMATION
# =============================================================================

sample_summary <- as.data.frame(
  table(
    sample_id,
    genotype
  )
)


sample_summary <- sample_summary[
  sample_summary$Freq > 0,
]


write.csv(
  sample_summary,
  "RESULTS/GSE141044_sample_summary.csv",
  row.names = FALSE
)


# =============================================================================
# 20. SAVE INITIAL SEURAT OBJECT
# =============================================================================

saveRDS(
  app23,
  "RESULTS/OBJECTS/GSE141044_initial_Seurat_object.rds"
)


# =============================================================================
# 21. SAVE SESSION INFORMATION
# =============================================================================

capture.output(
  sessionInfo(),
  file = "RESULTS/sessionInfo_initial_import.txt"
)


# =============================================================================
# 22. FINAL OUTPUT
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("GSE141044 IMPORT COMPLETED\n")
cat("####################################################################\n")


cat("\nSeurat object:\n")

print(app23)


cat("\nDimensions (genes x nuclei):\n")

print(
  dim(app23)
)


cat("\nUnique samples:\n")

print(
  sort(
    unique(
      app23$sample_id
    )
  )
)


cat("\nNumber of nuclei per sample:\n")

print(
  table(
    app23$sample_id
  )
)


cat("\nGenotypes:\n")

print(
  table(
    app23$genotype,
    useNA = "ifany"
  )
)


cat("\nQC summary:\n")

print(
  summary(
    app23@meta.data[
      ,
      c(
        "nCount_RNA",
        "nFeature_RNA",
        "percent.mt"
      )
    ]
  )
)


cat("\nFirst 10 metadata rows:\n")

print(
  head(
    app23@meta.data,
    10
  )
)


cat("\n####################################################################\n")
cat("STOP HERE BEFORE FILTERING / CLUSTERING\n")
cat("####################################################################\n")

###############################################################################
# GSE141044 — CORRECTED IMPORT + METADATA
#
# FIXES:
# 1. Reads the FIRST LINE manually as the 3280 nucleus IDs
# 2. Reads the remaining lines as gene x nucleus expression values
# 3. Correctly assigns M6/M24, APP23/WT, mouse replicate, and lane
# 4. Recreates the Seurat object from scratch
###############################################################################


# =============================================================================
# 0. CLEAN THE INCORRECT OBJECTS
# =============================================================================

rm(
  list = intersect(
    c(
      "app23",
      "counts",
      "counts_dense",
      "expr_dt",
      "cell_ids",
      "sample_id",
      "genotype"
    ),
    ls()
  )
)

gc()

set.seed(100)


# =============================================================================
# 1. PACKAGES
# =============================================================================

library(Matrix)
library(data.table)
library(Seurat)


# =============================================================================
# 2. PROJECT DIRECTORY
# =============================================================================

project_dir <- "D:/Master/ScRNA seq/Course project"

setwd(project_dir)


# =============================================================================
# 3. FIND ORIGINAL EXPRESSION FILE
# =============================================================================

all_files <- list.files(
  project_dir,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)


matrix_candidates <- all_files[
  grepl(
    "GSE141044_matrix\\.mtx$",
    all_files,
    ignore.case = TRUE
  )
]


if (length(matrix_candidates) == 0) {
  stop("Could not find GSE141044_matrix.mtx")
}


matrix_file <- matrix_candidates[1]


cat("\n============================================================\n")
cat("MATRIX FILE\n")
cat("============================================================\n")

print(matrix_file)


# =============================================================================
# 4. READ FIRST LINE SEPARATELY
#
# The first line contains ONLY the nucleus IDs.
# =============================================================================

cat("\n============================================================\n")
cat("READING NUCLEUS IDS FROM FIRST LINE\n")
cat("============================================================\n")


first_line <- readLines(
  matrix_file,
  n = 1
)


# Split on one or more whitespace characters
cell_ids <- strsplit(
  first_line,
  "\\s+"
)[[1]]


# Remove quotation marks
cell_ids <- gsub(
  '^"|"$',
  "",
  cell_ids
)


# Remove accidental empty strings
cell_ids <- cell_ids[
  cell_ids != ""
]


cat("\nNumber of nucleus IDs read:\n")
print(length(cell_ids))


cat("\nFirst 20 nucleus IDs:\n")
print(head(cell_ids, 20))


cat("\nLast 20 nucleus IDs:\n")
print(tail(cell_ids, 20))


# =============================================================================
# 5. VERIFY THAT THESE LOOK LIKE REAL GSE141044 CELL IDS
# =============================================================================

expected_pattern <- "_M(6|24)_(APP|WT)"

n_matching <- sum(
  grepl(
    expected_pattern,
    cell_ids,
    ignore.case = TRUE
  )
)


cat("\nNucleus IDs matching expected sample pattern:\n")
print(n_matching)


if (n_matching != length(cell_ids)) {
  
  warning(
    paste0(
      "Not every cell ID matched the expected M6/M24 APP/WT pattern. ",
      "We will inspect the unique suffixes below."
    )
  )
}


# =============================================================================
# 6. READ EXPRESSION DATA
#
# IMPORTANT:
# skip = 1 because first line has already been read as cell IDs
#
# header = FALSE because every remaining row begins with a gene name
# =============================================================================

cat("\n============================================================\n")
cat("READING EXPRESSION VALUES\n")
cat("THIS CAN TAKE A LITTLE TIME\n")
cat("============================================================\n")


expr_dt <- data.table::fread(
  matrix_file,
  skip = 1,
  header = FALSE,
  sep = " ",
  data.table = TRUE,
  check.names = FALSE,
  showProgress = TRUE
)


cat("\nRaw imported dimensions:\n")
print(dim(expr_dt))


# =============================================================================
# 7. CHECK NUMBER OF COLUMNS
#
# Expected:
#
# 1 gene column + 3280 nucleus columns
# =============================================================================

expected_columns <- length(cell_ids) + 1


cat("\nExpected number of columns:\n")
print(expected_columns)

cat("\nObserved number of columns:\n")
print(ncol(expr_dt))


if (ncol(expr_dt) != expected_columns) {
  
  stop(
    paste0(
      "\nCOLUMN MISMATCH\n\n",
      "Cell IDs: ",
      length(cell_ids),
      "\n",
      "Therefore expected data columns: ",
      expected_columns,
      "\n",
      "But fread returned: ",
      ncol(expr_dt),
      "\n\n",
      "Stop here and send me this output."
    )
  )
}


# =============================================================================
# 8. EXTRACT GENE NAMES
# =============================================================================

genes <- as.character(
  expr_dt[[1]]
)


# Strip quotes
genes <- gsub(
  '^"|"$',
  "",
  genes
)


genes <- trimws(
  genes
)


cat("\n============================================================\n")
cat("GENE NAMES\n")
cat("============================================================\n")


cat("\nNumber of gene rows:\n")
print(length(genes))


cat("\nFirst 20 genes:\n")
print(head(genes, 20))


cat("\nLast 20 genes:\n")
print(tail(genes, 20))


# =============================================================================
# 9. CLEAN BAD / DUPLICATED GENE NAMES
# =============================================================================

bad_genes <- (
  is.na(genes) |
    genes == ""
)


if (any(bad_genes)) {
  
  genes[bad_genes] <- paste0(
    "UnknownGene_",
    which(bad_genes)
  )
  
}


cat("\nDuplicated gene symbols before make.unique():\n")
print(sum(duplicated(genes)))


genes <- make.unique(
  genes,
  sep = "_dup"
)


cat("\nDuplicated gene symbols after make.unique():\n")
print(sum(duplicated(genes)))


# =============================================================================
# 10. REMOVE GENE COLUMN
# =============================================================================

expr_dt[, 1 := NULL]


cat("\nExpression-value dimensions:\n")
print(dim(expr_dt))


# =============================================================================
# 11. CONVERT TO NUMERIC MATRIX
# =============================================================================

cat("\n============================================================\n")
cat("CONVERTING TO MATRIX\n")
cat("============================================================\n")


counts_dense <- as.matrix(
  expr_dt
)


storage.mode(
  counts_dense
) <- "numeric"


rownames(
  counts_dense
) <- genes


colnames(
  counts_dense
) <- cell_ids


# =============================================================================
# 12. VERIFY CORRECT CELL NAMES BEFORE PROCEEDING
# =============================================================================

cat("\nFirst 20 matrix column names:\n")
print(head(colnames(counts_dense), 20))


cat("\nMatrix dimensions:\n")
print(dim(counts_dense))


if (
  !all(
    colnames(counts_dense) == cell_ids
  )
) {
  
  stop(
    "Matrix cell IDs were not assigned correctly."
  )
}


# =============================================================================
# 13. CONVERT TO SPARSE MATRIX
# =============================================================================

cat("\n============================================================\n")
cat("CONVERTING TO SPARSE MATRIX\n")
cat("============================================================\n")


counts <- Matrix::Matrix(
  counts_dense,
  sparse = TRUE
)


counts <- as(
  counts,
  "dgCMatrix"
)


rm(
  counts_dense,
  expr_dt
)

gc()


# =============================================================================
# 14. FINAL DUPLICATE CHECK
# =============================================================================

cat("\n============================================================\n")
cat("CHECKING ROW AND COLUMN NAMES\n")
cat("============================================================\n")


cat("\nDuplicated genes:\n")
print(sum(duplicated(rownames(counts))))


cat("\nDuplicated nuclei:\n")
print(sum(duplicated(colnames(counts))))


if (anyDuplicated(rownames(counts)) > 0) {
  
  rownames(counts) <- make.unique(
    rownames(counts),
    sep = "_dup"
  )
  
}


if (anyDuplicated(colnames(counts)) > 0) {
  
  colnames(counts) <- make.unique(
    colnames(counts),
    sep = "_dup"
  )
  
}


stopifnot(
  anyDuplicated(rownames(counts)) == 0
)

stopifnot(
  anyDuplicated(colnames(counts)) == 0
)


# =============================================================================
# 15. EXTRACT SAMPLE METADATA FROM NUCLEUS IDS
#
# Example:
#
# AACGTGAT_M6_APP1_L1
#
# barcode   AACGTGAT
# age       M6
# genotype  APP23
# replicate 1
# lane      L1
# =============================================================================

cat("\n============================================================\n")
cat("EXTRACTING SAMPLE METADATA\n")
cat("============================================================\n")


cell_ids <- colnames(
  counts
)


# -----------------------------------------------------------------------------
# Full sample string after first underscore
# -----------------------------------------------------------------------------

sample_full <- sub(
  "^[^_]+_",
  "",
  cell_ids
)


cat("\nUnique raw sample strings:\n")

print(
  sort(
    unique(sample_full)
  )
)


# -----------------------------------------------------------------------------
# AGE
# -----------------------------------------------------------------------------

age <- ifelse(
  grepl(
    "_M6_|^M6_",
    paste0("_", sample_full),
    ignore.case = TRUE
  ),
  "6_month",
  ifelse(
    grepl(
      "_M24_|^M24_",
      paste0("_", sample_full),
      ignore.case = TRUE
    ),
    "24_month",
    NA
  )
)


# A simpler direct check in case names start exactly M6/M24

age[
  grepl(
    "^M6_",
    sample_full,
    ignore.case = TRUE
  )
] <- "6_month"


age[
  grepl(
    "^M24_",
    sample_full,
    ignore.case = TRUE
  )
] <- "24_month"


# -----------------------------------------------------------------------------
# GENOTYPE
# -----------------------------------------------------------------------------

genotype <- ifelse(
  grepl(
    "_APP|^APP",
    sample_full,
    ignore.case = TRUE
  ),
  "APP23",
  ifelse(
    grepl(
      "_WT|^WT",
      sample_full,
      ignore.case = TRUE
    ),
    "WT",
    NA
  )
)


# More robust for format M6_APP1_L1 / M24_WT2_L1

genotype[
  grepl(
    "APP",
    sample_full,
    ignore.case = TRUE
  )
] <- "APP23"


genotype[
  grepl(
    "WT",
    sample_full,
    ignore.case = TRUE
  )
] <- "WT"


# -----------------------------------------------------------------------------
# REPLICATE NUMBER
#
# APP1 -> 1
# APP2 -> 2
# WT1  -> 1
# etc.
# -----------------------------------------------------------------------------

replicate <- sub(
  ".*_(?:APP|WT)([0-9]+)_.*",
  "\\1",
  sample_full,
  perl = TRUE
)


replicate[
  !grepl(
    "_(?:APP|WT)[0-9]+_",
    paste0("_", sample_full),
    perl = TRUE
  )
] <- NA


# -----------------------------------------------------------------------------
# LANE
# -----------------------------------------------------------------------------

lane <- ifelse(
  grepl(
    "_L[0-9]+$",
    sample_full
  ),
  sub(
    ".*_(L[0-9]+)$",
    "\\1",
    sample_full
  ),
  NA
)


# -----------------------------------------------------------------------------
# BIOLOGICAL SAMPLE ID
#
# Remove lane, because different sequencing lanes should not become
# independent biological replicates.
#
# M6_APP1_L1 -> M6_APP1
# -----------------------------------------------------------------------------

sample_id <- sub(
  "_L[0-9]+$",
  "",
  sample_full
)


# =============================================================================
# 16. INSPECT METADATA BEFORE CREATING SEURAT OBJECT
# =============================================================================

metadata_check <- data.frame(
  cell_id = cell_ids,
  sample_full = sample_full,
  sample_id = sample_id,
  age = age,
  genotype = genotype,
  replicate = replicate,
  lane = lane,
  stringsAsFactors = FALSE
)


cat("\nFirst 20 metadata rows:\n")

print(
  head(
    metadata_check,
    20
  )
)


cat("\nUnique biological sample IDs:\n")

print(
  sort(
    unique(sample_id)
  )
)


cat("\nNumber of unique biological samples:\n")

print(
  length(
    unique(sample_id)
  )
)


cat("\nNuclei per biological sample:\n")

print(
  table(sample_id)
)


cat("\nAge distribution:\n")

print(
  table(
    age,
    useNA = "ifany"
  )
)


cat("\nGenotype distribution:\n")

print(
  table(
    genotype,
    useNA = "ifany"
  )
)


cat("\nAge x genotype:\n")

print(
  table(
    age,
    genotype,
    useNA = "ifany"
  )
)


# =============================================================================
# 17. STOP IF METADATA FAILED
# =============================================================================

if (anyNA(age)) {
  
  stop(
    "Some nuclei could not be assigned an age."
  )
  
}


if (anyNA(genotype)) {
  
  stop(
    "Some nuclei could not be assigned a genotype."
  )
  
}


# =============================================================================
# 18. CREATE SEURAT OBJECT
# =============================================================================

cat("\n============================================================\n")
cat("CREATING CORRECT SEURAT OBJECT\n")
cat("============================================================\n")


app23 <- CreateSeuratObject(
  counts = counts,
  project = "GSE141044_APP23",
  min.cells = 0,
  min.features = 0
)


# =============================================================================
# 19. ADD CORRECT METADATA
# =============================================================================

app23$sample_full <- sample_full

app23$sample_id <- sample_id

app23$age <- age

app23$genotype <- genotype

app23$replicate <- replicate

app23$lane <- lane


app23$condition <- paste(
  app23$age,
  app23$genotype,
  sep = "_"
)


# =============================================================================
# 20. SET FACTOR ORDER
# =============================================================================

app23$age <- factor(
  app23$age,
  levels = c(
    "6_month",
    "24_month"
  )
)


app23$genotype <- factor(
  app23$genotype,
  levels = c(
    "WT",
    "APP23"
  )
)


app23$condition <- factor(
  app23$condition,
  levels = c(
    "6_month_WT",
    "6_month_APP23",
    "24_month_WT",
    "24_month_APP23"
  )
)


# =============================================================================
# 21. QC METRICS — CORRECTED
# =============================================================================

cat("\n============================================================\n")
cat("CALCULATING QC METRICS\n")
cat("============================================================\n")


# -----------------------------------------------------------------------------
# Check how mitochondrial genes are named in this dataset
# -----------------------------------------------------------------------------

mt_genes <- grep(
  "^[mM][tT]-",
  rownames(app23),
  value = TRUE
)

cat("\nNumber of mitochondrial genes detected:\n")
print(length(mt_genes))

cat("\nFirst mitochondrial genes detected:\n")
print(head(mt_genes, 20))


# -----------------------------------------------------------------------------
# Mitochondrial percentage
#
# [mM][tT]- matches:
# mt-
# Mt-
# mT-
# MT-
# -----------------------------------------------------------------------------

app23[["percent.mt"]] <- PercentageFeatureSet(
  app23,
  pattern = "^[mM][tT]-"
)


# -----------------------------------------------------------------------------
# Ribosomal percentage
#
# Matches Rpl..., Rps..., rpl..., rps...
# -----------------------------------------------------------------------------

app23[["percent.ribo"]] <- PercentageFeatureSet(
  app23,
  pattern = "^[Rr][Pp][LlSs]"
)


# -----------------------------------------------------------------------------
# QC summary
# -----------------------------------------------------------------------------

cat("\nQC summary:\n")

print(
  summary(
    app23@meta.data[
      ,
      c(
        "nCount_RNA",
        "nFeature_RNA",
        "percent.mt",
        "percent.ribo"
      )
    ]
  )
)


# -----------------------------------------------------------------------------
# Confirm metadata
# -----------------------------------------------------------------------------

cat("\nFirst 10 metadata rows:\n")

print(
  head(
    app23@meta.data,
    10
  )
)

# =============================================================================
# 22. FINAL VERIFICATION
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("CORRECT IMPORT COMPLETED\n")
cat("####################################################################\n")


cat("\nSeurat object:\n")

print(app23)


cat("\nGenes x nuclei:\n")

print(
  dim(app23)
)


cat("\nFirst 20 CORRECT cell names:\n")

print(
  head(
    colnames(app23),
    20
  )
)


cat("\nUnique biological samples:\n")

print(
  sort(
    unique(
      app23$sample_id
    )
  )
)


cat("\nNumber of biological samples:\n")

print(
  length(
    unique(
      app23$sample_id
    )
  )
)


cat("\nNuclei per biological sample:\n")

print(
  table(
    app23$sample_id
  )
)


cat("\nAge x genotype:\n")

print(
  table(
    app23$age,
    app23$genotype
  )
)


cat("\nCondition counts:\n")

print(
  table(
    app23$condition
  )
)


cat("\nQC summary:\n")

print(
  summary(
    app23@meta.data[
      ,
      c(
        "nCount_RNA",
        "nFeature_RNA",
        "percent.mt",
        "percent.ribo"
      )
    ]
  )
)


cat("\nFirst 10 metadata rows:\n")

print(
  head(
    app23@meta.data,
    10
  )
)


# =============================================================================
# 23. SAVE CORRECT OBJECT
# =============================================================================

dir.create(
  "RESULTS/OBJECTS",
  recursive = TRUE,
  showWarnings = FALSE
)


saveRDS(
  app23,
  "RESULTS/OBJECTS/GSE141044_correct_initial_Seurat_object.rds"
)


write.csv(
  app23@meta.data,
  "RESULTS/GSE141044_correct_metadata.csv",
  row.names = TRUE
)


cat("\n####################################################################\n")
cat("STOP HERE — DO NOT NORMALIZE YET\n")
cat("####################################################################\n")

# =============================================================================
# CHECK SUSPICIOUS "a" FEATURE BEFORE QC / NORMALIZATION
# =============================================================================

cat("\n============================================================\n")
cat("CHECKING FEATURE 'a'\n")
cat("============================================================\n")

cat("\nIs 'a' present?\n")
print("a" %in% rownames(app23))

if ("a" %in% rownames(app23)) {
  
  a_counts <- GetAssayData(
    app23,
    assay = "RNA",
    layer = "counts"
  )["a", ]
  
  cat("\nTotal counts for feature 'a':\n")
  print(sum(a_counts))
  
  cat("\nNumber of nuclei expressing 'a':\n")
  print(sum(a_counts > 0))
  
  cat("\nMaximum value of 'a':\n")
  print(max(a_counts))
  
  cat("\nFirst non-zero values:\n")
  print(head(a_counts[a_counts > 0], 20))
}

cat("\n============================================================\n")
cat("DO NOT NORMALIZE YET\n")
cat("============================================================\n")

###############################################################################
# GSE141044
# QC -> NORMALIZATION -> PCA -> UMAP -> ~10 FINE CLUSTERS
#
# START FROM THE CURRENT CORRECT app23 OBJECT
#
# IMPORTANT:
# - UMAP ONLY
# - NO t-SNE
# - NO Harmony/integration yet
# - Do not use percent.mt for filtering because mitochondrial genes are absent
###############################################################################


# =============================================================================
# 24. REMOVE THE ARTIFACT FEATURE "a"
# =============================================================================

cat("\n============================================================\n")
cat("REMOVING ARTIFACT FEATURE 'a'\n")
cat("============================================================\n")


cat("\nDimensions BEFORE removing 'a':\n")
print(dim(app23))


if ("a" %in% rownames(app23)) {
  
  app23 <- subset(
    app23,
    features = setdiff(
      rownames(app23),
      "a"
    )
  )
  
}


cat("\nDimensions AFTER removing 'a':\n")
print(dim(app23))


cat("\nIs 'a' still present?\n")
print("a" %in% rownames(app23))


# We expect:
#
# BEFORE: 28693 x 3280
# AFTER:  28692 x 3280


# =============================================================================
# 25. CORRECT orig.ident
#
# Seurat interpreted the barcode before the underscore as orig.ident.
# We already have the correct biological mouse/sample IDs.
# =============================================================================

app23$orig.ident <- app23$sample_id


cat("\nCorrect orig.ident distribution:\n")

print(
  table(
    app23$orig.ident
  )
)


# =============================================================================
# 26. VERIFY THAT EXPRESSION VALUES ARE INTEGER COUNTS
# =============================================================================

cat("\n============================================================\n")
cat("VERIFYING COUNT MATRIX\n")
cat("============================================================\n")


raw_counts <- GetAssayData(
  app23,
  assay = "RNA",
  layer = "counts"
)


nonzero_counts <- raw_counts@x


integer_check <- all(
  nonzero_counts == round(nonzero_counts)
)


cat("\nAre all non-zero expression values integers?\n")
print(integer_check)


if (!integer_check) {
  
  stop(
    paste0(
      "The matrix does not contain integer counts. ",
      "Stop before proceeding because this affects downstream pseudobulk analysis."
    )
  )
  
}


rm(
  raw_counts,
  nonzero_counts
)

gc()


# =============================================================================
# 27. QC — USE PAPER-APPROPRIATE THRESHOLDS
#
# Paper QC:
# detected genes > 500
# transcripts > 4000
#
# Our processed dataset already appears to have passed these thresholds.
#
# DO NOT FILTER ON percent.mt:
# mitochondrial genes are absent from this processed matrix.
# =============================================================================

cat("\n============================================================\n")
cat("QUALITY CONTROL\n")
cat("============================================================\n")


n_before <- ncol(app23)


cat("\nNumber of nuclei BEFORE QC:\n")
print(n_before)


cat("\nQC summary BEFORE filtering:\n")

print(
  summary(
    app23@meta.data[
      ,
      c(
        "nCount_RNA",
        "nFeature_RNA",
        "percent.ribo"
      )
    ]
  )
)


# Apply only the published/basic thresholds

app23 <- subset(
  app23,
  subset =
    nFeature_RNA > 500 &
    nCount_RNA > 4000
)


n_after <- ncol(app23)


cat("\nNumber of nuclei AFTER QC:\n")
print(n_after)


cat("\nNuclei removed by QC:\n")
print(n_before - n_after)


cat("\nPercentage retained:\n")
print(
  round(
    100 * n_after / n_before,
    2
  )
)


cat("\nNuclei per biological sample AFTER QC:\n")

print(
  table(
    app23$sample_id
  )
)


cat("\nAge x genotype AFTER QC:\n")

print(
  table(
    app23$age,
    app23$genotype
  )
)


# =============================================================================
# 28. CREATE QC FIGURES
# =============================================================================

dir.create(
  "FIGURES/QC",
  recursive = TRUE,
  showWarnings = FALSE
)


pdf(
  "FIGURES/QC/QC_violin_plots.pdf",
  width = 11,
  height = 6
)

print(
  VlnPlot(
    app23,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.ribo"
    ),
    group.by = "sample_id",
    pt.size = 0,
    ncol = 3
  )
)

dev.off()


pdf(
  "FIGURES/QC/QC_nCount_vs_nFeature.pdf",
  width = 7,
  height = 6
)

print(
  FeatureScatter(
    app23,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )
)

dev.off()


# =============================================================================
# 29. SAVE POST-QC OBJECT
# =============================================================================

saveRDS(
  app23,
  "RESULTS/OBJECTS/GSE141044_post_QC.rds"
)


# =============================================================================
# 30. NORMALIZATION
# =============================================================================

cat("\n============================================================\n")
cat("NORMALIZATION\n")
cat("============================================================\n")


app23 <- NormalizeData(
  app23,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)


# =============================================================================
# 31. HIGHLY VARIABLE FEATURES
# =============================================================================

cat("\n============================================================\n")
cat("FINDING VARIABLE FEATURES\n")
cat("============================================================\n")


app23 <- FindVariableFeatures(
  app23,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)


cat("\nNumber of variable genes:\n")
print(length(VariableFeatures(app23)))


cat("\nTop 20 variable genes:\n")
print(head(VariableFeatures(app23), 20))


# =============================================================================
# 32. VARIABLE FEATURE PLOT
# =============================================================================

pdf(
  "FIGURES/QC/Highly_variable_genes.pdf",
  width = 8,
  height = 6
)

print(
  VariableFeaturePlot(
    app23
  )
)

dev.off()


# =============================================================================
# 33. SCALE DATA
# =============================================================================

cat("\n============================================================\n")
cat("SCALING DATA\n")
cat("============================================================\n")


app23 <- ScaleData(
  app23,
  assay = "RNA",
  features = VariableFeatures(app23),
  verbose = TRUE
)


# =============================================================================
# 34. PCA
# =============================================================================

cat("\n============================================================\n")
cat("RUNNING PCA\n")
cat("============================================================\n")


set.seed(100)


app23 <- RunPCA(
  app23,
  assay = "RNA",
  features = VariableFeatures(app23),
  npcs = 30,
  verbose = TRUE
)


cat("\nPCA completed.\n")


# =============================================================================
# 35. PCA OUTPUT
# =============================================================================

print(
  app23[["pca"]],
  dims = 1:10,
  nfeatures = 10
)

# =============================================================================
# 36. PCA ELBOW PLOT
# =============================================================================

pdf(
  "FIGURES/QC/PCA_elbow_plot.pdf",
  width = 7,
  height = 5
)

print(
  ElbowPlot(
    app23,
    ndims = 30
  )
)

dev.off()


# =============================================================================
# 37. PCA PLOTS — CHECK SAMPLE / CONDITION STRUCTURE
# =============================================================================

pdf(
  "FIGURES/QC/PCA_by_sample.pdf",
  width = 9,
  height = 7
)

print(
  DimPlot(
    app23,
    reduction = "pca",
    group.by = "sample_id"
  )
)

dev.off()


pdf(
  "FIGURES/QC/PCA_by_condition.pdf",
  width = 8,
  height = 6
)

print(
  DimPlot(
    app23,
    reduction = "pca",
    group.by = "condition"
  )
)

dev.off()


# =============================================================================
# 38. NUMBER OF PCs
#
# The original study used the first 10 PCs.
# We therefore begin with PCs 1:10 for biological comparability.
#
# We can later change this if the elbow plot clearly indicates otherwise.
# =============================================================================

pcs_use <- 1:10


cat("\nPCs being used:\n")
print(pcs_use)


# Reviewed and confirmed by Hala - clustering, UMAP and annotation section
# =============================================================================
# 39. BUILD NEAREST-NEIGHBOR GRAPH
# =============================================================================

cat("\n============================================================\n")
cat("BUILDING NEIGHBOR GRAPH\n")
cat("============================================================\n")


app23 <- FindNeighbors(
  app23,
  reduction = "pca",
  dims = pcs_use,
  verbose = TRUE
)


# =============================================================================
# 40. RESOLUTION SEARCH
#
# Goal:
# approximately 10 fine neuronal clusters
#
# We test multiple resolutions and automatically choose the one whose
# number of clusters is closest to 10.
# =============================================================================

cat("\n============================================================\n")
cat("SEARCHING FOR APPROXIMATELY 10 CLUSTERS\n")
cat("============================================================\n")


resolution_grid <- seq(
  0.1,
  2.0,
  by = 0.1
)


cluster_numbers <- numeric(
  length(resolution_grid)
)


for (i in seq_along(resolution_grid)) {
  
  current_resolution <- resolution_grid[i]
  
  cat(
    "\nTesting resolution:",
    current_resolution,
    "\n"
  )
  
  
  app23 <- FindClusters(
    app23,
    resolution = current_resolution,
    algorithm = 1,
    random.seed = 100,
    verbose = FALSE
  )
  
  
  cluster_numbers[i] <- length(
    unique(
      app23$seurat_clusters
    )
  )
  
  
  cat(
    "Clusters:",
    cluster_numbers[i],
    "\n"
  )
  
}


resolution_results <- data.frame(
  resolution = resolution_grid,
  n_clusters = cluster_numbers
)


cat("\n============================================================\n")
cat("RESOLUTION RESULTS\n")
cat("============================================================\n")

print(
  resolution_results
)


# =============================================================================
# 41. SELECT RESOLUTION CLOSEST TO 10 CLUSTERS
# =============================================================================

target_clusters <- 10


resolution_results$difference_from_10 <- abs(
  resolution_results$n_clusters -
    target_clusters
)


best_row <- which.min(
  resolution_results$difference_from_10
)


best_resolution <- resolution_results$resolution[
  best_row
]


best_cluster_number <- resolution_results$n_clusters[
  best_row
]


cat("\nSelected resolution:\n")
print(best_resolution)


cat("\nExpected number of clusters:\n")
print(best_cluster_number)


# Save resolution table

write.csv(
  resolution_results,
  "RESULTS/GSE141044_resolution_search.csv",
  row.names = FALSE
)


# =============================================================================
# 42. FINAL CLUSTERING AT SELECTED RESOLUTION
# =============================================================================

cat("\n============================================================\n")
cat("FINAL CLUSTERING\n")
cat("============================================================\n")


set.seed(100)


app23 <- FindClusters(
  app23,
  resolution = best_resolution,
  algorithm = 1,
  random.seed = 100,
  verbose = TRUE
)


cat("\nFinal number of clusters:\n")

print(
  length(
    unique(
      app23$seurat_clusters
    )
  )
)


cat("\nNuclei per cluster:\n")

print(
  table(
    app23$seurat_clusters
  )
)


# =============================================================================
# 43. RUN UMAP
#
# NO t-SNE
# =============================================================================

cat("\n============================================================\n")
cat("RUNNING UMAP\n")
cat("============================================================\n")


set.seed(100)


app23 <- RunUMAP(
  app23,
  reduction = "pca",
  dims = pcs_use,
  seed.use = 100,
  verbose = TRUE
)


# =============================================================================
# 44. UMAP — CLUSTERS
# =============================================================================

dir.create(
  "FIGURES/UMAP",
  recursive = TRUE,
  showWarnings = FALSE
)


p_cluster <- DimPlot(
  app23,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    paste0(
      "GSE141044 — Fine neuronal clusters\nResolution = ",
      best_resolution
    )
  )


ggsave(
  "FIGURES/UMAP/UMAP_clusters.pdf",
  plot = p_cluster,
  width = 8,
  height = 7
)


ggsave(
  "FIGURES/UMAP/UMAP_clusters.png",
  plot = p_cluster,
  width = 8,
  height = 7,
  dpi = 300
)


print(p_cluster)


# =============================================================================
# 45. UMAP — BIOLOGICAL SAMPLE
# =============================================================================

p_sample <- DimPlot(
  app23,
  reduction = "umap",
  group.by = "sample_id"
) +
  ggtitle(
    "UMAP by biological sample"
  )


ggsave(
  "FIGURES/UMAP/UMAP_by_sample.pdf",
  plot = p_sample,
  width = 9,
  height = 7
)


# =============================================================================
# 46. UMAP — GENOTYPE
# =============================================================================

p_genotype <- DimPlot(
  app23,
  reduction = "umap",
  group.by = "genotype"
) +
  ggtitle(
    "UMAP by genotype"
  )


ggsave(
  "FIGURES/UMAP/UMAP_by_genotype.pdf",
  plot = p_genotype,
  width = 7,
  height = 6
)


# =============================================================================
# 47. UMAP — AGE
# =============================================================================

p_age <- DimPlot(
  app23,
  reduction = "umap",
  group.by = "age"
) +
  ggtitle(
    "UMAP by age"
  )


ggsave(
  "FIGURES/UMAP/UMAP_by_age.pdf",
  plot = p_age,
  width = 7,
  height = 6
)


# =============================================================================
# 48. UMAP — FOUR CONDITIONS
# =============================================================================

p_condition <- DimPlot(
  app23,
  reduction = "umap",
  group.by = "condition"
) +
  ggtitle(
    "UMAP by age/genotype condition"
  )


ggsave(
  "FIGURES/UMAP/UMAP_by_condition.pdf",
  plot = p_condition,
  width = 8,
  height = 7
)


# =============================================================================
# 49. UMAP SPLIT BY CONDITION
# =============================================================================

p_condition_split <- DimPlot(
  app23,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "condition",
  label = FALSE,
  ncol = 2
)


ggsave(
  "FIGURES/UMAP/UMAP_clusters_split_by_condition.pdf",
  plot = p_condition_split,
  width = 14,
  height = 10
)


# =============================================================================
# 50. UMAP SPLIT BY SAMPLE
# =============================================================================

p_sample_split <- DimPlot(
  app23,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "sample_id",
  label = FALSE,
  ncol = 4
)


ggsave(
  "FIGURES/UMAP/UMAP_clusters_split_by_sample.pdf",
  plot = p_sample_split,
  width = 18,
  height = 14
)


# =============================================================================
# 51. SAVE CLUSTER COMPOSITION TABLES
# =============================================================================

cluster_by_sample <- as.data.frame(
  table(
    cluster = app23$seurat_clusters,
    sample = app23$sample_id
  )
)


write.csv(
  cluster_by_sample,
  "RESULTS/cluster_by_sample_counts.csv",
  row.names = FALSE
)


cluster_by_condition <- as.data.frame(
  table(
    cluster = app23$seurat_clusters,
    condition = app23$condition
  )
)


write.csv(
  cluster_by_condition,
  "RESULTS/cluster_by_condition_counts.csv",
  row.names = FALSE
)


# =============================================================================
# 52. CHECK PAPER MARKER GENES
#
# Broad neuronal subtypes reported in the paper:
#
# DG          -> Prox1
# CA1         -> Mpped1
# CA3         -> Mndal
# inhibitory  -> Gad1
#
# Plus useful neuronal genes.
# =============================================================================

cat("\n============================================================\n")
cat("CHECKING PAPER MARKERS\n")
cat("============================================================\n")


paper_markers <- c(
  "Prox1",
  "Mpped1",
  "Mndal",
  "Gad1",
  "Gad2",
  "Slc17a7",
  "Camk2a",
  "Rbfox3",
  "Snap25",
  "Syp"
)


available_markers <- paper_markers[
  paper_markers %in% rownames(app23)
]


missing_markers <- setdiff(
  paper_markers,
  available_markers
)


cat("\nAvailable markers:\n")
print(available_markers)


cat("\nMissing markers:\n")
print(missing_markers)


# =============================================================================
# 53. DOTPLOT OF PAPER MARKERS
# =============================================================================

if (length(available_markers) > 0) {
  
  p_dot <- DotPlot(
    app23,
    features = available_markers,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis() +
    ggtitle(
      "Known hippocampal neuronal markers"
    )
  
  
  ggsave(
    "FIGURES/UMAP/Paper_marker_DotPlot.pdf",
    plot = p_dot,
    width = 10,
    height = 6
  )
  
  
  print(p_dot)
  
}


# =============================================================================
# 54. FEATURE PLOTS OF KEY PAPER MARKERS
# =============================================================================

key_markers <- c(
  "Prox1",
  "Mpped1",
  "Mndal",
  "Gad1"
)


key_markers <- key_markers[
  key_markers %in% rownames(app23)
]


if (length(key_markers) > 0) {
  
  p_features <- FeaturePlot(
    app23,
    features = key_markers,
    reduction = "umap",
    ncol = 2
  )
  
  
  ggsave(
    "FIGURES/UMAP/Paper_marker_FeaturePlots.pdf",
    plot = p_features,
    width = 11,
    height = 9
  )
  
}


# =============================================================================
# 55. SAVE OBJECT BEFORE MARKER ANALYSIS
# =============================================================================

saveRDS(
  app23,
  "RESULTS/OBJECTS/GSE141044_normalized_PCA_UMAP_clustered.rds"
)


# =============================================================================
# 56. FINAL SUMMARY
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("QC + NORMALIZATION + PCA + UMAP + CLUSTERING COMPLETED\n")
cat("####################################################################\n")


cat("\nFinal dimensions:\n")
print(dim(app23))


cat("\nNumber of nuclei:\n")
print(ncol(app23))


cat("\nNumber of biological samples:\n")
print(length(unique(app23$sample_id)))


cat("\nSelected PCs:\n")
print(pcs_use)


cat("\nSelected clustering resolution:\n")
print(best_resolution)


cat("\nFinal number of fine clusters:\n")
print(
  length(
    unique(
      app23$seurat_clusters
    )
  )
)


cat("\nNuclei per cluster:\n")
print(
  table(
    app23$seurat_clusters
  )
)


cat("\nNuclei per sample:\n")
print(
  table(
    app23$sample_id
  )
)


cat("\nAge x genotype:\n")
print(
  table(
    app23$age,
    app23$genotype
  )
)


cat("\n####################################################################\n")
cat("STOP HERE BEFORE CELL-TYPE ANNOTATION / DIFFERENTIAL EXPRESSION\n")
cat("####################################################################\n")
###############################################################################
# GSE141044
# CLUSTER MARKERS + INITIAL BIOLOGICAL ANNOTATION
###############################################################################


# =============================================================================
# 57. FIND MARKERS FOR ALL 10 FINE CLUSTERS
# =============================================================================

cat("\n============================================================\n")
cat("FINDING MARKERS FOR ALL FINE CLUSTERS\n")
cat("============================================================\n")


Idents(app23) <- "seurat_clusters"


all_markers <- FindAllMarkers(
  app23,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)


cat("\nNumber of marker rows identified:\n")
print(nrow(all_markers))


# =============================================================================
# 58. SAVE ALL MARKERS
# =============================================================================

dir.create(
  "RESULTS/MARKERS",
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  all_markers,
  "RESULTS/MARKERS/GSE141044_all_cluster_markers.csv",
  row.names = FALSE
)


# =============================================================================
# 59. TOP 20 MARKERS PER CLUSTER
# =============================================================================

library(dplyr)


top20_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()


write.csv(
  top20_markers,
  "RESULTS/MARKERS/GSE141044_top20_markers_per_cluster.csv",
  row.names = FALSE
)


cat("\nTop 20 markers per cluster:\n")

print(
  top20_markers %>%
    select(
      cluster,
      gene,
      avg_log2FC,
      pct.1,
      pct.2,
      p_val_adj
    )
)


# =============================================================================
# 60. TOP 10 MARKERS PER CLUSTER FOR HEATMAP
# =============================================================================

top10_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()


top10_genes <- unique(
  top10_markers$gene
)


# =============================================================================
# 61. SCALE TOP MARKERS FOR HEATMAP
# =============================================================================

app23 <- ScaleData(
  app23,
  features = top10_genes,
  verbose = FALSE
)


# =============================================================================
# 62. HEATMAP
# =============================================================================

dir.create(
  "FIGURES/MARKERS",
  recursive = TRUE,
  showWarnings = FALSE
)


pdf(
  "FIGURES/MARKERS/Top10_cluster_marker_heatmap.pdf",
  width = 13,
  height = 12
)


print(
  DoHeatmap(
    app23,
    features = top10_genes,
    group.by = "seurat_clusters",
    size = 3
  ) +
    NoLegend()
)


dev.off()


# =============================================================================
# 63. EXPANDED HIPPOCAMPAL MARKER PANEL
#
# We are now using a broader marker set rather than assigning clusters
# from only Prox1 / Mpped1 / Mndal / Gad1.
# =============================================================================

hippocampal_markers <- c(
  
  # Dentate gyrus
  "Prox1",
  "Calb1",
  "C1ql2",
  
  # CA1-associated
  "Mpped1",
  "Wfs1",
  "Pcp4",
  "Rgs14",
  
  # CA3-associated
  "Mndal",
  "Grik4",
  "Bok",
  
  # Excitatory neurons
  "Slc17a7",
  "Camk2a",
  "Satb2",
  
  # Inhibitory neurons
  "Gad1",
  "Gad2",
  "Slc6a1",
  "Slc32a1",
  
  # General neuronal
  "Rbfox3",
  "Snap25",
  "Syp"
)


hippocampal_markers_available <- hippocampal_markers[
  hippocampal_markers %in% rownames(app23)
]


cat("\nAvailable expanded hippocampal markers:\n")
print(hippocampal_markers_available)


cat("\nMissing expanded markers:\n")
print(
  setdiff(
    hippocampal_markers,
    hippocampal_markers_available
  )
)


# =============================================================================
# 64. EXPANDED DOTPLOT
# =============================================================================

p_expanded <- DotPlot(
  app23,
  features = hippocampal_markers_available,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle(
    "Expanded hippocampal neuronal marker panel"
  )


ggsave(
  "FIGURES/MARKERS/Expanded_hippocampal_marker_DotPlot.pdf",
  plot = p_expanded,
  width = 14,
  height = 7
)


ggsave(
  "FIGURES/MARKERS/Expanded_hippocampal_marker_DotPlot.png",
  plot = p_expanded,
  width = 14,
  height = 7,
  dpi = 300
)


print(p_expanded)


# =============================================================================
# 65. INITIAL BROAD ANNOTATION
#
# This is provisional.
# We will confirm against the top marker table before locking it in.
# =============================================================================

broad_annotation <- c(
  "0" = "DG",
  "1" = "CA1_like",
  "2" = "CA1_like",
  "3" = "CA1_like",
  "4" = "CA3",
  "5" = "CA1_like",
  "6" = "Inhibitory",
  "7" = "Inhibitory",
  "8" = "CA1_like",
  "9" = "Inhibitory"
)


# Remove cluster names from the resulting vector before adding to metadata
app23$broad_subtype <- unname(
  broad_annotation[
    as.character(app23$seurat_clusters)
  ]
)


# =============================================================================
# 66. CHECK BROAD ANNOTATION
# =============================================================================

cat("\n============================================================\n")
cat("PROVISIONAL BROAD ANNOTATION\n")
cat("============================================================\n")


print(
  table(
    app23$seurat_clusters,
    app23$broad_subtype
  )
)


cat("\nNuclei per broad subtype:\n")

print(
  table(
    app23$broad_subtype
  )
)


# =============================================================================
# 67. BROAD SUBTYPE UMAP
# =============================================================================

p_broad <- DimPlot(
  app23,
  reduction = "umap",
  group.by = "broad_subtype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle(
    "GSE141044 — provisional broad neuronal subtypes"
  )


ggsave(
  "FIGURES/MARKERS/UMAP_broad_subtypes_provisional.pdf",
  plot = p_broad,
  width = 8,
  height = 7
)


ggsave(
  "FIGURES/MARKERS/UMAP_broad_subtypes_provisional.png",
  plot = p_broad,
  width = 8,
  height = 7,
  dpi = 300
)


print(p_broad)


# =============================================================================
# 68. BROAD SUBTYPE BY CONDITION
# =============================================================================

p_broad_condition <- DimPlot(
  app23,
  reduction = "umap",
  group.by = "broad_subtype",
  split.by = "condition",
  ncol = 2
)


ggsave(
  "FIGURES/MARKERS/UMAP_broad_subtypes_by_condition.pdf",
  plot = p_broad_condition,
  width = 14,
  height = 10
)


# =============================================================================
# 69. SAVE PROVISIONALLY ANNOTATED OBJECT
# =============================================================================

saveRDS(
  app23,
  "RESULTS/OBJECTS/GSE141044_provisional_annotation.rds"
)


# =============================================================================
# 70. FINAL OUTPUT
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("MARKER ANALYSIS COMPLETED\n")
cat("####################################################################\n")


cat("\nTop marker genes by cluster:\n")

for (cl in sort(unique(all_markers$cluster))) {
  
  cat("\n----------------------------------------\n")
  cat("CLUSTER:", cl, "\n")
  cat("----------------------------------------\n")
  
  tmp <- top20_markers %>%
    filter(cluster == cl) %>%
    select(
      gene,
      avg_log2FC,
      pct.1,
      pct.2,
      p_val_adj
    )
  
  print(tmp)
}


cat("\n####################################################################\n")
cat("STOP HERE BEFORE FINALIZING ANNOTATION\n")
cat("####################################################################\n")

###############################################################################
# GSE141044 — FINAL ANNOTATION DIAGNOSTIC
###############################################################################


# =============================================================================
# 71. DO NOT USE THE PROVISIONAL ANNOTATION FOR ANALYSIS YET
# =============================================================================

cat("\n============================================================\n")
cat("FINAL CELL-TYPE ANNOTATION DIAGNOSTIC\n")
cat("============================================================\n")


# =============================================================================
# 72. DEFINE DIAGNOSTIC MARKERS
# =============================================================================

diagnostic_markers <- c(
  
  # Dentate gyrus
  "Prox1",
  "C1ql2",
  
  # CA1-associated
  "Mpped1",
  "Wfs1",
  "Pcp4",
  
  # CA3-associated
  "Mndal",
  "Grik4",
  
  # Excitatory identity
  "Slc17a7",
  "Camk2a",
  
  # Alternative glutamatergic marker
  "Slc17a6",
  
  # Inhibitory identity
  "Gad1",
  "Gad2",
  "Slc32a1",
  
  # Interneuron subclasses
  "Pvalb",
  "Sst",
  "Tac1",
  "Lhx6",
  "Htr3a",
  "Reln",
  
  # Cluster 8-associated
  "Rorb",
  
  # Cluster 9-associated
  "Zic1",
  "Zic4",
  "Prkcd"
)


diagnostic_markers <- diagnostic_markers[
  diagnostic_markers %in% rownames(app23)
]


cat("\nMarkers available:\n")
print(diagnostic_markers)


# =============================================================================
# 73. DOTPLOT
# =============================================================================

p_diagnostic <- DotPlot(
  app23,
  features = diagnostic_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle(
    "Diagnostic markers for final hippocampal annotation"
  )


ggsave(
  "FIGURES/MARKERS/Final_annotation_diagnostic_DotPlot.pdf",
  plot = p_diagnostic,
  width = 14,
  height = 7
)


ggsave(
  "FIGURES/MARKERS/Final_annotation_diagnostic_DotPlot.png",
  plot = p_diagnostic,
  width = 14,
  height = 7,
  dpi = 300
)


print(p_diagnostic)


# =============================================================================
# 74. AVERAGE EXPRESSION BY CLUSTER
# =============================================================================

avg_expression <- AverageExpression(
  app23,
  assays = "RNA",
  features = diagnostic_markers,
  group.by = "seurat_clusters",
  layer = "data",
  verbose = FALSE
)$RNA


cat("\n============================================================\n")
cat("AVERAGE EXPRESSION BY CLUSTER\n")
cat("============================================================\n")

print(
  round(
    avg_expression,
    3
  )
)


write.csv(
  avg_expression,
  "RESULTS/MARKERS/Final_annotation_average_expression.csv"
)


# =============================================================================
# 75. PERCENT OF CELLS EXPRESSING EACH MARKER
# =============================================================================

data_matrix <- GetAssayData(
  app23,
  assay = "RNA",
  layer = "data"
)


pct_expression <- sapply(
  sort(unique(app23$seurat_clusters)),
  function(cl) {
    
    cells_cl <- colnames(app23)[
      app23$seurat_clusters == cl
    ]
    
    Matrix::rowMeans(
      data_matrix[
        diagnostic_markers,
        cells_cl,
        drop = FALSE
      ] > 0
    ) * 100
    
  }
)


colnames(pct_expression) <- paste0(
  "Cluster_",
  sort(unique(app23$seurat_clusters))
)


cat("\n============================================================\n")
cat("PERCENT OF NUCLEI EXPRESSING EACH MARKER\n")
cat("============================================================\n")

print(
  round(
    pct_expression,
    1
  )
)


write.csv(
  pct_expression,
  "RESULTS/MARKERS/Final_annotation_percent_expression.csv"
)


# =============================================================================
# 76. SPECIAL CHECK — CLUSTER 9
# =============================================================================

cat("\n============================================================\n")
cat("SPECIAL CHECK — CLUSTER 9\n")
cat("============================================================\n")


cluster9_markers <- intersect(
  c(
    "Gad1",
    "Gad2",
    "Slc32a1",
    "Slc17a7",
    "Slc17a6",
    "Camk2a",
    "Zic1",
    "Zic4",
    "Prkcd"
  ),
  rownames(app23)
)


print(
  DotPlot(
    app23,
    features = cluster9_markers,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis()
)


cat("\n####################################################################\n")
cat("STOP HERE — FINAL ANNOTATION NEXT\n")
cat("####################################################################\n")

###############################################################################
# GSE141044 — SAMPLE REPRESENTATION WITHIN EACH CLUSTER
###############################################################################


# =============================================================================
# 77. FINAL BROAD ANNOTATION
# =============================================================================

final_broad_annotation <- c(
  "0" = "DG",
  "1" = "CA1_like",
  "2" = "CA1_like",
  "3" = "CA1_like",
  "4" = "CA3",
  "5" = "CA1_like",
  "6" = "Inhibitory",
  "7" = "Inhibitory",
  "8" = "CA1_like",
  "9" = "Inhibitory"
)


app23$broad_subtype <- unname(
  final_broad_annotation[
    as.character(app23$seurat_clusters)
  ]
)


# =============================================================================
# 78. VERIFY FINAL BROAD ANNOTATION
# =============================================================================

cat("\n============================================================\n")
cat("FINAL BROAD ANNOTATION\n")
cat("============================================================\n")


print(
  table(
    app23$seurat_clusters,
    app23$broad_subtype
  )
)


cat("\nBroad subtype totals:\n")

print(
  table(
    app23$broad_subtype
  )
)


# =============================================================================
# 79. CLUSTER x BIOLOGICAL SAMPLE
# =============================================================================

cat("\n============================================================\n")
cat("NUCLEI PER CLUSTER PER BIOLOGICAL MOUSE\n")
cat("============================================================\n")


cluster_sample_table <- table(
  Cluster = app23$seurat_clusters,
  Sample = app23$sample_id
)


print(
  cluster_sample_table
)


write.csv(
  as.data.frame.matrix(cluster_sample_table),
  "RESULTS/cluster_by_biological_sample_matrix.csv"
)


# =============================================================================
# 80. CLUSTER x CONDITION
# =============================================================================

cat("\n============================================================\n")
cat("NUCLEI PER CLUSTER PER CONDITION\n")
cat("============================================================\n")


cluster_condition_table <- table(
  Cluster = app23$seurat_clusters,
  Condition = app23$condition
)


print(
  cluster_condition_table
)


write.csv(
  as.data.frame.matrix(cluster_condition_table),
  "RESULTS/cluster_by_condition_matrix.csv"
)


# =============================================================================
# 81. BROAD SUBTYPE x BIOLOGICAL SAMPLE
# =============================================================================

cat("\n============================================================\n")
cat("NUCLEI PER BROAD SUBTYPE PER BIOLOGICAL MOUSE\n")
cat("============================================================\n")


broad_sample_table <- table(
  Subtype = app23$broad_subtype,
  Sample = app23$sample_id
)


print(
  broad_sample_table
)


write.csv(
  as.data.frame.matrix(broad_sample_table),
  "RESULTS/broad_subtype_by_biological_sample_matrix.csv"
)


# =============================================================================
# 82. BROAD SUBTYPE x CONDITION
# =============================================================================

cat("\n============================================================\n")
cat("NUCLEI PER BROAD SUBTYPE PER CONDITION\n")
cat("============================================================\n")


broad_condition_table <- table(
  Subtype = app23$broad_subtype,
  Condition = app23$condition
)


print(
  broad_condition_table
)


# =============================================================================
# 83. SAMPLE PERCENTAGE COMPOSITION
# =============================================================================

sample_composition <- prop.table(
  cluster_sample_table,
  margin = 2
) * 100


cat("\n============================================================\n")
cat("PERCENT OF EACH MOUSE IN EACH FINE CLUSTER\n")
cat("============================================================\n")


print(
  round(
    sample_composition,
    2
  )
)


write.csv(
  as.data.frame.matrix(sample_composition),
  "RESULTS/cluster_percentage_by_biological_sample.csv"
)


# =============================================================================
# 84. CHECK NUMBER OF MICE REPRESENTED IN EACH CLUSTER
# =============================================================================

mice_per_cluster <- rowSums(
  cluster_sample_table > 0
)


cat("\n============================================================\n")
cat("NUMBER OF BIOLOGICAL MICE REPRESENTED PER CLUSTER\n")
cat("============================================================\n")


print(
  mice_per_cluster
)


# =============================================================================
# 85. CHECK NUMBER OF MICE BY AGE/GENOTYPE WITHIN EACH CLUSTER
# =============================================================================

cluster_metadata <- data.frame(
  cluster = app23$seurat_clusters,
  sample_id = app23$sample_id,
  age = app23$age,
  genotype = app23$genotype
)


cluster_metadata <- unique(
  cluster_metadata
)


cat("\n============================================================\n")
cat("BIOLOGICAL REPLICATES AVAILABLE WITHIN EACH CLUSTER\n")
cat("============================================================\n")


replicate_check <- aggregate(
  sample_id ~ cluster + age + genotype,
  data = cluster_metadata,
  FUN = function(x) length(unique(x))
)


names(replicate_check)[
  names(replicate_check) == "sample_id"
] <- "n_biological_mice"


print(
  replicate_check
)


write.csv(
  replicate_check,
  "RESULTS/biological_replicates_per_cluster.csv",
  row.names = FALSE
)


# =============================================================================
# 86. SAVE FINAL ANNOTATED OBJECT
# =============================================================================

saveRDS(
  app23,
  "RESULTS/OBJECTS/GSE141044_final_annotated_before_DE.rds"
)


cat("\n####################################################################\n")
cat("SAMPLE REPRESENTATION CHECK COMPLETED\n")
cat("####################################################################\n")

###############################################################################
# GSE141044
# PRIMARY PSEUDOBULK DIFFERENTIAL EXPRESSION
#
# Biological question:
# Which hippocampal neuronal subtypes undergo the greatest transcriptional
# remodeling in APP23 mice, and which remain relatively preserved?
#
# PRIMARY UNIT OF REPLICATION = BIOLOGICAL MOUSE
###############################################################################


# =============================================================================
# 87. LOAD PACKAGES
# =============================================================================

library(Seurat)
library(Matrix)
library(edgeR)
library(dplyr)
library(ggplot2)


set.seed(100)


# =============================================================================
# 88. OUTPUT DIRECTORIES
# =============================================================================

dir.create(
  "RESULTS/PSEUDOBULK",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "FIGURES/PSEUDOBULK",
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 89. VERIFY REQUIRED METADATA
# =============================================================================

required_metadata <- c(
  "sample_id",
  "age",
  "genotype",
  "broad_subtype"
)


missing_metadata <- setdiff(
  required_metadata,
  colnames(app23@meta.data)
)


if (length(missing_metadata) > 0) {
  
  stop(
    paste(
      "Missing metadata:",
      paste(missing_metadata, collapse = ", ")
    )
  )
}


cat("\n============================================================\n")
cat("METADATA CHECK\n")
cat("============================================================\n")


print(
  table(
    app23$age,
    app23$genotype
  )
)


print(
  table(
    app23$broad_subtype
  )
)


# =============================================================================
# 90. GET RAW COUNTS
#
# IMPORTANT:
# Differential expression uses RAW COUNTS, not normalized Seurat data.
# =============================================================================

raw_counts <- GetAssayData(
  app23,
  assay = "RNA",
  layer = "counts"
)


cat("\nRaw count dimensions:\n")
print(dim(raw_counts))


# =============================================================================
# 91. CREATE PSEUDOBULK GROUP
#
# One pseudobulk library =
# one biological mouse + one broad neuronal subtype
#
# Example:
# M6_APP1__DG
# M6_WT1__DG
# =============================================================================

app23$pseudobulk_id <- paste(
  app23$sample_id,
  app23$broad_subtype,
  sep = "__"
)


cat("\nNumber of pseudobulk libraries:\n")

print(
  length(
    unique(app23$pseudobulk_id)
  )
)


cat("\nNuclei per pseudobulk library:\n")

print(
  sort(
    table(app23$pseudobulk_id)
  )
)


# =============================================================================
# 92. AGGREGATE RAW COUNTS
#
# Efficient sparse-matrix aggregation.
# =============================================================================

pb_factor <- factor(
  app23$pseudobulk_id
)


aggregation_matrix <- sparse.model.matrix(
  ~ 0 + pb_factor
)


colnames(aggregation_matrix) <- levels(
  pb_factor
)


pseudobulk_counts <- raw_counts %*% aggregation_matrix


cat("\nPseudobulk count matrix dimensions:\n")
print(dim(pseudobulk_counts))


# Expected approximately:
#
# 28692 genes x 44 pseudobulk libraries
#
# because:
# 11 mice x 4 broad subtypes = 44


# =============================================================================
# 93. BUILD PSEUDOBULK SAMPLE METADATA
# =============================================================================

pb_metadata <- app23@meta.data %>%
  
  mutate(
    pseudobulk_id = paste(
      sample_id,
      broad_subtype,
      sep = "__"
    )
  ) %>%
  
  select(
    pseudobulk_id,
    sample_id,
    age,
    genotype,
    broad_subtype
  ) %>%
  
  distinct()


pb_metadata <- as.data.frame(
  pb_metadata
)


rownames(pb_metadata) <- pb_metadata$pseudobulk_id


# Put metadata in exactly the same order as count columns

pb_metadata <- pb_metadata[
  colnames(pseudobulk_counts),
  ,
  drop = FALSE
]


stopifnot(
  identical(
    rownames(pb_metadata),
    colnames(pseudobulk_counts)
  )
)


cat("\nPseudobulk metadata:\n")
print(pb_metadata)


write.csv(
  pb_metadata,
  "RESULTS/PSEUDOBULK/pseudobulk_metadata.csv",
  row.names = TRUE
)


# =============================================================================
# 94. SAVE PSEUDOBULK COUNTS
# =============================================================================

write.csv(
  as.matrix(pseudobulk_counts),
  "RESULTS/PSEUDOBULK/pseudobulk_raw_counts.csv"
)


# =============================================================================
# 95. FUNCTION FOR APP23 VS WT WITHIN ONE AGE AND SUBTYPE
# =============================================================================

run_pseudobulk_DE <- function(
    counts_matrix,
    metadata,
    subtype_name,
    age_name
) {
  
  cat("\n\n")
  cat("============================================================\n")
  cat("PSEUDOBULK DIFFERENTIAL EXPRESSION\n")
  cat("Subtype:", subtype_name, "\n")
  cat("Age:", age_name, "\n")
  cat("Comparison: APP23 vs WT\n")
  cat("============================================================\n")
  
  
  # ---------------------------------------------------------------------------
  # Select pseudobulk samples
  # ---------------------------------------------------------------------------
  
  keep_samples <- (
    metadata$broad_subtype == subtype_name &
      metadata$age == age_name
  )
  
  
  meta_sub <- metadata[
    keep_samples,
    ,
    drop = FALSE
  ]
  
  
  counts_sub <- counts_matrix[
    ,
    rownames(meta_sub),
    drop = FALSE
  ]
  
  
  cat("\nBiological samples included:\n")
  
  print(
    meta_sub[
      ,
      c(
        "sample_id",
        "genotype"
      )
    ]
  )
  
  
  cat("\nReplicates:\n")
  
  print(
    table(
      meta_sub$genotype
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Require both genotypes
  # ---------------------------------------------------------------------------
  
  if (
    length(
      unique(meta_sub$genotype)
    ) < 2
  ) {
    
    cat("\nSKIPPED: both genotypes are not represented.\n")
    
    return(NULL)
  }
  
  
  # ---------------------------------------------------------------------------
  # Make genotype explicit
  #
  # WT = reference
  # Therefore positive logFC = higher in APP23
  # ---------------------------------------------------------------------------
  
  genotype <- factor(
    meta_sub$genotype,
    levels = c(
      "WT",
      "APP23"
    )
  )
  
  
  design <- model.matrix(
    ~ genotype
  )
  
  
  rownames(design) <- rownames(
    meta_sub
  )
  
  
  cat("\nDesign matrix:\n")
  print(design)
  
  
  # ---------------------------------------------------------------------------
  # edgeR object
  # ---------------------------------------------------------------------------
  
  y <- DGEList(
    counts = counts_sub,
    samples = meta_sub
  )
  
  
  # ---------------------------------------------------------------------------
  # Expression filtering
  #
  # filterByExpr accounts for library sizes and experimental design.
  # ---------------------------------------------------------------------------
  
  keep_genes <- filterByExpr(
    y,
    design = design
  )
  
  
  cat("\nGenes before expression filtering:\n")
  print(nrow(y))
  
  
  cat("\nGenes retained by filterByExpr:\n")
  print(sum(keep_genes))
  
  
  y <- y[
    keep_genes,
    ,
    keep.lib.sizes = FALSE
  ]
  
  
  # ---------------------------------------------------------------------------
  # TMM normalization
  # ---------------------------------------------------------------------------
  
  y <- calcNormFactors(
    y,
    method = "TMM"
  )
  
  
  # ---------------------------------------------------------------------------
  # Estimate dispersion
  # ---------------------------------------------------------------------------
  
  y <- estimateDisp(
    y,
    design,
    robust = TRUE
  )
  
  
  # ---------------------------------------------------------------------------
  # Quasi-likelihood model
  # ---------------------------------------------------------------------------
  
  fit <- glmQLFit(
    y,
    design,
    robust = TRUE
  )
  
  
  qlf <- glmQLFTest(
    fit,
    coef = "genotypeAPP23"
  )
  
  
  # ---------------------------------------------------------------------------
  # Complete results
  # ---------------------------------------------------------------------------
  
  results <- topTags(
    qlf,
    n = Inf,
    sort.by = "PValue"
  )$table
  
  
  results$gene <- rownames(
    results
  )
  
  
  results$subtype <- subtype_name
  
  results$age <- age_name
  
  
  # Positive logFC:
  # APP23 > WT
  #
  # Negative logFC:
  # APP23 < WT
  
  results$direction <- ifelse(
    results$FDR < 0.05 &
      results$logFC > 0.25,
    "Up_in_APP23",
    
    ifelse(
      results$FDR < 0.05 &
        results$logFC < -0.25,
      "Down_in_APP23",
      "NS"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Save table
  # ---------------------------------------------------------------------------
  
  output_name <- paste0(
    "APP23_vs_WT_",
    age_name,
    "_",
    subtype_name,
    ".csv"
  )
  
  
  write.csv(
    results,
    file.path(
      "RESULTS/PSEUDOBULK",
      output_name
    ),
    row.names = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # Summary
  # ---------------------------------------------------------------------------
  
  n_up <- sum(
    results$direction == "Up_in_APP23"
  )
  
  
  n_down <- sum(
    results$direction == "Down_in_APP23"
  )
  
  
  n_deg <- n_up + n_down
  
  
  cat("\nSignificant DEGs:\n")
  
  cat(
    "Up in APP23:",
    n_up,
    "\n"
  )
  
  cat(
    "Down in APP23:",
    n_down,
    "\n"
  )
  
  cat(
    "Total:",
    n_deg,
    "\n"
  )
  
  
  # ---------------------------------------------------------------------------
  # Volcano plot
  # ---------------------------------------------------------------------------
  
  volcano_data <- results
  
  
  volcano_data$minus_log10_FDR <- -log10(
    pmax(
      volcano_data$FDR,
      1e-300
    )
  )
  
  
  volcano <- ggplot(
    volcano_data,
    aes(
      x = logFC,
      y = minus_log10_FDR
    )
  ) +
    
    geom_point(
      aes(
        shape = direction
      ),
      alpha = 0.6,
      size = 1.5
    ) +
    
    geom_vline(
      xintercept = c(
        -0.25,
        0.25
      ),
      linetype = "dashed"
    ) +
    
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    
    labs(
      title = paste(
        subtype_name,
        "-",
        age_name,
        "- APP23 vs WT"
      ),
      subtitle = paste0(
        "DEGs = ",
        n_deg,
        " (",
        n_up,
        " up, ",
        n_down,
        " down)"
      ),
      x = "log2 fold change (APP23 / WT)",
      y = "-log10(FDR)"
    ) +
    
    theme_classic()
  
  
  figure_name <- paste0(
    "Volcano_APP23_vs_WT_",
    age_name,
    "_",
    subtype_name,
    ".pdf"
  )
  
  
  ggsave(
    file.path(
      "FIGURES/PSEUDOBULK",
      figure_name
    ),
    plot = volcano,
    width = 7,
    height = 6
  )
  
  
  # ---------------------------------------------------------------------------
  # Return everything needed later
  # ---------------------------------------------------------------------------
  
  return(
    list(
      results = results,
      y = y,
      fit = fit,
      qlf = qlf,
      n_up = n_up,
      n_down = n_down,
      n_deg = n_deg,
      n_tested = nrow(results)
    )
  )
}


# =============================================================================
# 96. DEFINE ANALYSES
# =============================================================================

broad_subtypes <- c(
  "DG",
  "CA1_like",
  "CA3",
  "Inhibitory"
)


ages <- c(
  "6_month",
  "24_month"
)


# =============================================================================
# 97. RUN ALL 8 PRIMARY COMPARISONS
#
# 4 neuronal subtypes x 2 ages
# =============================================================================
# =============================================================================
# 97. RUN ALL 8 PRIMARY COMPARISONS — CORRECTED
#
# 4 neuronal subtypes x 2 ages
# =============================================================================

DE_results <- list()


for (current_age in ages) {
  
  for (current_subtype in broad_subtypes) {
    
    comparison_name <- paste(
      current_age,
      current_subtype,
      sep = "__"
    )
    
    
    cat(
      "\n\nRunning comparison:",
      comparison_name,
      "\n"
    )
    
    
    DE_results[[comparison_name]] <- run_pseudobulk_DE(
      counts_matrix = pseudobulk_counts,
      metadata = pb_metadata,
      subtype_name = current_subtype,
      age_name = current_age
    )
    
  }
  
}


cat("\n============================================================\n")
cat("ALL PRIMARY COMPARISONS FINISHED\n")
cat("============================================================\n")


cat("\nComparisons stored:\n")
print(names(DE_results))


cat("\nNumber of comparisons:\n")
print(length(DE_results))
# =============================================================================
# 98. BUILD TRANSCRIPTIONAL REMODELING SUMMARY
# =============================================================================

remodeling_summary <- data.frame()


for (comparison_name in names(DE_results)) {
  
  x <- DE_results[[comparison_name]] 
  
  if (is.null(x)) {
    next
  }
  
  
  comparison_parts <- strsplit(
    comparison_name,
    "__",
    fixed = TRUE
  )[[1]]
  
  
  current_age <- comparison_parts[1]
  
  current_subtype <- comparison_parts[2]
  
  
  sig_results <- x$results[
    x$results$FDR < 0.05 &
      abs(x$results$logFC) > 0.25,
    ,
    drop = FALSE
  ]
  
  
  mean_abs_logFC <- if (
    nrow(sig_results) > 0
  ) {
    
    mean(
      abs(
        sig_results$logFC
      )
    )
    
  } else {
    
    0
    
  }
  
  
  median_abs_logFC <- if (
    nrow(sig_results) > 0
  ) {
    
    median(
      abs(
        sig_results$logFC
      )
    )
    
  } else {
    
    0
    
  }
  
  
  remodeling_summary <- rbind(
    remodeling_summary,
    
    data.frame(
      age = current_age,
      subtype = current_subtype,
      genes_tested = x$n_tested,
      DEGs_up_APP23 = x$n_up,
      DEGs_down_APP23 = x$n_down,
      total_DEGs = x$n_deg,
      mean_abs_logFC_DEGs = mean_abs_logFC,
      median_abs_logFC_DEGs = median_abs_logFC
    )
  )
  
}


# =============================================================================
# 99. DEG FRACTION
#
# This helps account for differences in the number of genes tested.
# =============================================================================

remodeling_summary$DEG_fraction <- (
  remodeling_summary$total_DEGs /
    remodeling_summary$genes_tested
)


remodeling_summary$DEG_percent <- (
  remodeling_summary$DEG_fraction * 100
)


# =============================================================================
# 100. PRINT REMODELING SUMMARY
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("TRANSCRIPTIONAL REMODELING SUMMARY\n")
cat("####################################################################\n")


print(
  remodeling_summary[
    order(
      remodeling_summary$age,
      -remodeling_summary$total_DEGs
    ),
  ]
)


# =============================================================================
# 101. SAVE REMODELING SUMMARY
# =============================================================================

write.csv(
  remodeling_summary,
  "RESULTS/PSEUDOBULK/transcriptional_remodeling_summary.csv",
  row.names = FALSE
)


# =============================================================================
# 102. FIGURE — NUMBER OF DEGs
# =============================================================================

p_deg_number <- ggplot(
  remodeling_summary,
  aes(
    x = subtype,
    y = total_DEGs,
    fill = age
  )
) +
  
  geom_col(
    position = "dodge"
  ) +
  
  labs(
    title = "Magnitude of APP23-associated transcriptional remodeling",
    x = "Neuronal subtype",
    y = "Number of significant DEGs"
  ) +
  
  theme_classic()


ggsave(
  "FIGURES/PSEUDOBULK/Remodeling_total_DEGs.pdf",
  plot = p_deg_number,
  width = 8,
  height = 6
)


print(p_deg_number)


# =============================================================================
# 103. FIGURE — PERCENT OF TESTED GENES DIFFERENTIALLY EXPRESSED
# =============================================================================

p_deg_percent <- ggplot(
  remodeling_summary,
  aes(
    x = subtype,
    y = DEG_percent,
    fill = age
  )
) +
  
  geom_col(
    position = "dodge"
  ) +
  
  labs(
    title = "Fraction of transcriptome remodeled in APP23",
    x = "Neuronal subtype",
    y = "Significant DEGs (% of genes tested)"
  ) +
  
  theme_classic()


ggsave(
  "FIGURES/PSEUDOBULK/Remodeling_DEG_percentage.pdf",
  plot = p_deg_percent,
  width = 8,
  height = 6
)


print(p_deg_percent)


# =============================================================================
# 104. SAVE COMPLETE ANALYSIS OBJECTS
# =============================================================================

saveRDS(
  DE_results,
  "RESULTS/PSEUDOBULK/GSE141044_pseudobulk_DE_results.rds"
)


saveRDS(
  app23,
  "RESULTS/OBJECTS/GSE141044_after_pseudobulk_setup.rds"
)


# =============================================================================
# 105. FINAL OUTPUT
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("PRIMARY PSEUDOBULK ANALYSIS COMPLETED\n")
cat("####################################################################\n")


cat("\nInterpretation of logFC:\n")
cat("Positive = higher expression in APP23\n")
cat("Negative = lower expression in APP23\n")


cat("\nRemodeling summary:\n")

print(
  remodeling_summary
)


cat("\n####################################################################\n")
cat("STOP HERE BEFORE GO / PATHWAY / TF ANALYSIS\n")
cat("####################################################################\n")

###############################################################################
# GSE141044
# INSPECT STRICT PSEUDOBULK RESULTS BEFORE PATHWAY ANALYSIS
###############################################################################


# =============================================================================
# 106. COMBINE ALL 8 DIFFERENTIAL EXPRESSION TABLES
# =============================================================================

all_DE_tables <- list()


for (comparison_name in names(DE_results)) {
  
  x <- DE_results[[comparison_name]]
  
  if (is.null(x)) {
    next
  }
  
  tmp <- x$results
  
  tmp$comparison <- comparison_name
  
  all_DE_tables[[comparison_name]] <- tmp
}


combined_DE <- do.call(
  rbind,
  all_DE_tables
)


rownames(combined_DE) <- NULL


write.csv(
  combined_DE,
  "RESULTS/PSEUDOBULK/all_8_comparisons_complete_DE_tables.csv",
  row.names = FALSE
)


# =============================================================================
# 107. EXTRACT STRICT SIGNIFICANT DEGs
#
# Same threshold used previously:
#
# FDR < 0.05
# |logFC| > 0.25
# =============================================================================

strict_DEGs <- combined_DE[
  combined_DE$FDR < 0.05 &
    abs(combined_DE$logFC) > 0.25,
  ,
  drop = FALSE
]


strict_DEGs <- strict_DEGs[
  order(
    strict_DEGs$comparison,
    strict_DEGs$FDR
  ),
]


cat("\n\n")
cat("####################################################################\n")
cat("ALL STRICT SIGNIFICANT DEGs\n")
cat("####################################################################\n")


print(
  strict_DEGs[
    ,
    c(
      "comparison",
      "gene",
      "logFC",
      "logCPM",
      "F",
      "PValue",
      "FDR",
      "direction"
    )
  ],
  row.names = FALSE
)


write.csv(
  strict_DEGs,
  "RESULTS/PSEUDOBULK/strict_significant_DEGs_all_comparisons.csv",
  row.names = FALSE
)


# =============================================================================
# 108. CHECK WHETHER APP IS AMONG THE SIGNIFICANT GENES
# =============================================================================

cat("\n============================================================\n")
cat("APP FEATURE CHECK\n")
cat("============================================================\n")


APP_results <- combined_DE[
  toupper(combined_DE$gene) == "APP",
  ,
  drop = FALSE
]


print(
  APP_results[
    ,
    c(
      "comparison",
      "gene",
      "logFC",
      "logCPM",
      "PValue",
      "FDR"
    )
  ],
  row.names = FALSE
)


write.csv(
  APP_results,
  "RESULTS/PSEUDOBULK/APP_feature_across_comparisons.csv",
  row.names = FALSE
)


# =============================================================================
# 109. PRINT SIGNIFICANT GENES SEPARATELY FOR EACH COMPARISON
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("SIGNIFICANT GENES BY COMPARISON\n")
cat("####################################################################\n")


for (comparison_name in names(DE_results)) {
  
  x <- DE_results[[comparison_name]]
  
  if (is.null(x)) {
    next
  }
  
  
  sig <- x$results[
    x$results$FDR < 0.05 &
      abs(x$results$logFC) > 0.25,
    ,
    drop = FALSE
  ]
  
  
  sig <- sig[
    order(sig$FDR),
    ,
    drop = FALSE
  ]
  
  
  cat("\n============================================================\n")
  cat(comparison_name, "\n")
  cat("============================================================\n")
  
  
  if (nrow(sig) == 0) {
    
    cat("No significant genes.\n")
    
  } else {
    
    print(
      sig[
        ,
        c(
          "gene",
          "logFC",
          "logCPM",
          "F",
          "PValue",
          "FDR"
        )
      ],
      row.names = FALSE
    )
    
  }
  
}


# =============================================================================
# 110. TOP 20 GENES BY P-VALUE FOR EACH COMPARISON
#
# These are NOT automatically significant.
# This is diagnostic/ranking information.
# =============================================================================

top20_all <- list()


for (comparison_name in names(DE_results)) {
  
  x <- DE_results[[comparison_name]]
  
  if (is.null(x)) {
    next
  }
  
  
  tmp <- x$results[
    order(x$results$PValue),
    ,
    drop = FALSE
  ]
  
  
  tmp <- head(
    tmp,
    20
  )
  
  
  tmp$comparison <- comparison_name
  
  
  top20_all[[comparison_name]] <- tmp
  
  
  cat("\n\n============================================================\n")
  cat("TOP 20 RANKED GENES:", comparison_name, "\n")
  cat("============================================================\n")
  
  
  print(
    tmp[
      ,
      c(
        "gene",
        "logFC",
        "logCPM",
        "PValue",
        "FDR"
      )
    ],
    row.names = FALSE
  )
  
}


top20_combined <- do.call(
  rbind,
  top20_all
)


rownames(top20_combined) <- NULL


write.csv(
  top20_combined,
  "RESULTS/PSEUDOBULK/top20_ranked_genes_each_comparison.csv",
  row.names = FALSE
)


# =============================================================================
# 111. REMODELING SUMMARY EXCLUDING APP
#
# IMPORTANT:
# We are NOT deleting APP from the statistical analysis.
#
# This is a sensitivity analysis showing whether the ranking is being
# dominated by the APP transgene itself.
# =============================================================================

remodeling_without_APP <- data.frame()


for (comparison_name in names(DE_results)) {
  
  x <- DE_results[[comparison_name]]
  
  if (is.null(x)) {
    next
  }
  
  
  comparison_parts <- strsplit(
    comparison_name,
    "__",
    fixed = TRUE
  )[[1]]
  
  
  current_age <- comparison_parts[1]
  
  current_subtype <- comparison_parts[2]
  
  
  tmp <- x$results[
    toupper(x$results$gene) != "APP",
    ,
    drop = FALSE
  ]
  
  
  sig <- tmp[
    tmp$FDR < 0.05 &
      abs(tmp$logFC) > 0.25,
    ,
    drop = FALSE
  ]
  
  
  remodeling_without_APP <- rbind(
    remodeling_without_APP,
    
    data.frame(
      age = current_age,
      subtype = current_subtype,
      genes_tested_excluding_APP = nrow(tmp),
      DEGs_excluding_APP = nrow(sig),
      up_excluding_APP = sum(
        sig$logFC > 0
      ),
      down_excluding_APP = sum(
        sig$logFC < 0
      ),
      mean_abs_logFC_excluding_APP =
        if (nrow(sig) > 0)
          mean(abs(sig$logFC))
      else
        0,
      median_abs_logFC_excluding_APP =
        if (nrow(sig) > 0)
          median(abs(sig$logFC))
      else
        0
    )
  )
  
}


cat("\n\n")
cat("####################################################################\n")
cat("SENSITIVITY ANALYSIS — EXCLUDING APP FEATURE\n")
cat("####################################################################\n")


print(
  remodeling_without_APP
)


write.csv(
  remodeling_without_APP,
  "RESULTS/PSEUDOBULK/remodeling_summary_excluding_APP.csv",
  row.names = FALSE
)


# =============================================================================
# 112. FDR DISTRIBUTION CHECK
#
# Useful for determining whether many genes show weaker coordinated signals
# that do not reach FDR because biological replication is limited.
# =============================================================================

FDR_summary <- data.frame()


for (comparison_name in names(DE_results)) {
  
  x <- DE_results[[comparison_name]]
  
  if (is.null(x)) {
    next
  }
  
  
  res <- x$results
  
  
  FDR_summary <- rbind(
    FDR_summary,
    
    data.frame(
      comparison = comparison_name,
      
      FDR_lt_0.05 =
        sum(res$FDR < 0.05),
      
      FDR_lt_0.10 =
        sum(res$FDR < 0.10),
      
      FDR_lt_0.20 =
        sum(res$FDR < 0.20),
      
      nominal_P_lt_0.01 =
        sum(res$PValue < 0.01),
      
      nominal_P_lt_0.05 =
        sum(res$PValue < 0.05)
    )
  )
  
}


cat("\n\n")
cat("####################################################################\n")
cat("FDR / NOMINAL SIGNAL SUMMARY\n")
cat("####################################################################\n")


print(
  FDR_summary
)


write.csv(
  FDR_summary,
  "RESULTS/PSEUDOBULK/FDR_signal_summary.csv",
  row.names = FALSE
)


cat("\n####################################################################\n")
cat("STOP HERE BEFORE PATHWAY ENRICHMENT\n")
cat("####################################################################\n")

###############################################################################
# GSE141044
# RANKED GENE SET ENRICHMENT ANALYSIS (GSEA)
#
# Sections 113 onward
#
# PRIMARY ANALYSIS:
#   - Mouse-level pseudobulk edgeR results
#   - APP23 vs WT
#   - Separate analysis for each neuronal subtype and age
#   - Entire ranked transcriptome
#   - Uppercase APP and Thy1 excluded from PRIMARY ranking
#   - Endogenous mouse App retained
###############################################################################


# =============================================================================
# 113. INSTALL / LOAD REQUIRED PACKAGES
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


bioc_packages <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "enrichplot"
)


for (pkg in bioc_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
    
  }
  
}


library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(dplyr)
library(ggplot2)


set.seed(100)


# =============================================================================
# 114. CREATE OUTPUT DIRECTORIES
# =============================================================================

dir.create(
  "RESULTS/GSEA",
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  "FIGURES/GSEA",
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 115. VERIFY DIFFERENTIAL EXPRESSION RESULTS
# =============================================================================

if (!exists("DE_results")) {
  
  stop(
    "DE_results is not present. Load the pseudobulk DE results first."
  )
  
}


cat("\n============================================================\n")
cat("GSEA INPUT CHECK\n")
cat("============================================================\n")


cat("\nNumber of comparisons:\n")
print(length(DE_results))


cat("\nComparison names:\n")
print(names(DE_results))


# =============================================================================
# 116. BUILD RANKED GENE LIST
#
# We use:
#
#     sign(logFC) * sqrt(F)
#
# edgeR's QL F statistic is always positive.
# Multiplying by the direction of logFC gives us a signed ranking:
#
# Positive = APP23-associated increase
# Negative = APP23-associated decrease
#
# This incorporates BOTH:
#   direction of effect
#   statistical evidence
#
# Primary ranking excludes:
#   APP   = likely APP23 transgene-derived feature
#   Thy1  = promoter/model-associated feature
#
# IMPORTANT:
# Endogenous mouse "App" is NOT removed.
# =============================================================================

make_ranked_list <- function(
    result_table,
    remove_model_features = TRUE
) {
  
  dat <- result_table
  
  
  # Remove missing/non-finite statistics
  
  dat <- dat[
    !is.na(dat$gene) &
      !is.na(dat$logFC) &
      !is.na(dat$F) &
      is.finite(dat$logFC) &
      is.finite(dat$F),
    ,
    drop = FALSE
  ]
  
  
  # ---------------------------------------------------------------------------
  # Remove construct-associated features for PRIMARY analysis
  #
  # Case-sensitive intentionally:
  #
  # APP removed
  # App retained
  # ---------------------------------------------------------------------------
  
  if (remove_model_features) {
    
    dat <- dat[
      !(dat$gene %in% c("APP", "Thy1")),
      ,
      drop = FALSE
    ]
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Create signed statistic
  # ---------------------------------------------------------------------------
  
  dat$rank_stat <- sign(
    dat$logFC
  ) * sqrt(
    pmax(
      dat$F,
      0
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Remove duplicated gene symbols
  #
  # If duplicates occur, retain the row with the largest absolute statistic.
  # ---------------------------------------------------------------------------
  
  dat <- dat[
    order(
      abs(dat$rank_stat),
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  
  dat <- dat[
    !duplicated(dat$gene),
    ,
    drop = FALSE
  ]
  
  
  # ---------------------------------------------------------------------------
  # Named numeric vector required by clusterProfiler
  # ---------------------------------------------------------------------------
  
  gene_list <- dat$rank_stat
  
  names(gene_list) <- dat$gene
  
  
  gene_list <- sort(
    gene_list,
    decreasing = TRUE
  )
  
  
  return(gene_list)
}


# =============================================================================
# 117. TEST RANKING ON ONE COMPARISON
# =============================================================================

test_comparison <- names(DE_results)[1]


test_rank <- make_ranked_list(
  DE_results[[test_comparison]]$results,
  remove_model_features = TRUE
)


cat("\n============================================================\n")
cat("TEST RANKING\n")
cat("============================================================\n")


cat("\nComparison:\n")
print(test_comparison)


cat("\nNumber of ranked genes:\n")
print(length(test_rank))


cat("\nTop 10 APP23-associated genes:\n")
print(head(test_rank, 10))


cat("\nTop 10 WT-associated genes:\n")
print(tail(test_rank, 10))


cat("\nIs uppercase APP present?\n")
print("APP" %in% names(test_rank))


cat("\nIs Thy1 present?\n")
print("Thy1" %in% names(test_rank))


cat("\nIs endogenous App present?\n")
print("App" %in% names(test_rank))


# Expected:
#
# APP  -> FALSE
# Thy1 -> FALSE
# App  -> TRUE (if it passed expression filtering)


# =============================================================================
# 118. MAP MOUSE SYMBOLS TO ENTREZ IDs
#
# clusterProfiler gseGO works most reliably with ENTREZ gene IDs.
# =============================================================================

convert_rank_to_entrez <- function(gene_list) {
  
  symbols <- names(gene_list)
  
  
  mapping <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = "ENTREZID"
  )
  
  
  mapping <- mapping[
    !is.na(mapping$ENTREZID),
    ,
    drop = FALSE
  ]
  
  
  # Attach ranking statistic
  
  mapping$rank_stat <- gene_list[
    mapping$SYMBOL
  ]
  
  
  mapping <- mapping[
    !is.na(mapping$rank_stat),
    ,
    drop = FALSE
  ]
  
  
  # ---------------------------------------------------------------------------
  # One Entrez ID can occasionally map to multiple symbols.
  #
  # Keep the strongest absolute statistic.
  # ---------------------------------------------------------------------------
  
  mapping <- mapping[
    order(
      abs(mapping$rank_stat),
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  
  mapping <- mapping[
    !duplicated(mapping$ENTREZID),
    ,
    drop = FALSE
  ]
  
  
  entrez_rank <- mapping$rank_stat
  
  names(entrez_rank) <- mapping$ENTREZID
  
  
  entrez_rank <- sort(
    entrez_rank,
    decreasing = TRUE
  )
  
  
  return(
    list(
      ranked_vector = entrez_rank,
      mapping = mapping
    )
  )
}


# =============================================================================
# 119. TEST ENTREZ CONVERSION
# =============================================================================

test_conversion <- convert_rank_to_entrez(
  test_rank
)


cat("\n============================================================\n")
cat("GENE-ID CONVERSION CHECK\n")
cat("============================================================\n")


cat("\nSymbols before mapping:\n")
print(length(test_rank))


cat("\nEntrez genes after mapping:\n")
print(length(test_conversion$ranked_vector))


cat("\nMapping percentage:\n")

print(
  round(
    100 *
      length(test_conversion$ranked_vector) /
      length(test_rank),
    2
  )
)


# =============================================================================
# 120. GSEA FUNCTION — GO BIOLOGICAL PROCESS
# =============================================================================

run_GO_GSEA <- function(
    result_table,
    comparison_name,
    remove_model_features = TRUE
) {
  
  cat("\n\n")
  cat("####################################################################\n")
  cat("RUNNING GO GSEA\n")
  cat("Comparison:", comparison_name, "\n")
  
  if (remove_model_features) {
    
    cat("Ranking: PRIMARY — APP and Thy1 excluded\n")
    
  } else {
    
    cat("Ranking: SENSITIVITY — all genes retained\n")
    
  }
  
  cat("####################################################################\n")
  
  
  # ---------------------------------------------------------------------------
  # Build ranked list
  # ---------------------------------------------------------------------------
  
  symbol_rank <- make_ranked_list(
    result_table,
    remove_model_features = remove_model_features
  )
  
  
  # ---------------------------------------------------------------------------
  # Convert to Entrez
  # ---------------------------------------------------------------------------
  
  conversion <- convert_rank_to_entrez(
    symbol_rank
  )
  
  
  entrez_rank <- conversion$ranked_vector
  
  
  cat("\nRanked symbols:\n")
  print(length(symbol_rank))
  
  
  cat("\nMapped Entrez genes:\n")
  print(length(entrez_rank))
  
  
  # ---------------------------------------------------------------------------
  # Run GO Biological Process GSEA
  #
  # pvalueCutoff = 1 intentionally:
  # Save the complete result first.
  # Significance is assessed afterward using adjusted P values.
  # ---------------------------------------------------------------------------
  
  gsea <- gseGO(
    geneList = entrez_rank,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea"
  )
  
  
  gsea_table <- as.data.frame(
    gsea
  )
  
  
  # ---------------------------------------------------------------------------
  # Add comparison
  # ---------------------------------------------------------------------------
  
  if (nrow(gsea_table) > 0) {
    
    gsea_table$comparison <- comparison_name
    
    gsea_table$analysis <- if (
      remove_model_features
    ) {
      "Primary_APP_Thy1_removed"
    } else {
      "Sensitivity_all_genes"
    }
    
  }
  
  
  return(
    list(
      gsea = gsea,
      table = gsea_table,
      symbol_rank = symbol_rank,
      entrez_rank = entrez_rank,
      mapping = conversion$mapping
    )
  )
}


# =============================================================================
# 121. RUN PRIMARY GSEA FOR ALL 8 COMPARISONS
#
# APP and Thy1 removed.
# =============================================================================

GSEA_primary <- list()


for (comparison_name in names(DE_results)) {
  
  cat(
    "\n\nStarting PRIMARY GSEA:",
    comparison_name,
    "\n"
  )
  
  
  GSEA_primary[[comparison_name]] <- run_GO_GSEA(
    result_table = DE_results[[comparison_name]]$results,
    comparison_name = comparison_name,
    remove_model_features = TRUE
  )
  
}


# =============================================================================
# 122. COMBINE PRIMARY GSEA RESULTS
# =============================================================================

primary_tables <- lapply(
  GSEA_primary,
  function(x) x$table
)


primary_tables <- primary_tables[
  sapply(
    primary_tables,
    nrow
  ) > 0
]


GSEA_primary_combined <- do.call(
  rbind,
  primary_tables
)


rownames(
  GSEA_primary_combined
) <- NULL


write.csv(
  GSEA_primary_combined,
  "RESULTS/GSEA/GSEA_GO_BP_PRIMARY_all_comparisons.csv",
  row.names = FALSE
)


# =============================================================================
# 123. SIGNIFICANT PRIMARY PATHWAYS
#
# FDR < 0.05
# =============================================================================

GSEA_primary_significant <- GSEA_primary_combined[
  !is.na(GSEA_primary_combined$p.adjust) &
    GSEA_primary_combined$p.adjust < 0.05,
  ,
  drop = FALSE
]


GSEA_primary_significant <- GSEA_primary_significant[
  order(
    GSEA_primary_significant$comparison,
    GSEA_primary_significant$p.adjust
  ),
  ,
  drop = FALSE
]


write.csv(
  GSEA_primary_significant,
  "RESULTS/GSEA/GSEA_GO_BP_PRIMARY_significant_FDR05.csv",
  row.names = FALSE
)


cat("\n\n")
cat("####################################################################\n")
cat("PRIMARY GSEA SUMMARY\n")
cat("####################################################################\n")


primary_summary <- data.frame()


for (comparison_name in names(GSEA_primary)) {
  
  tmp <- GSEA_primary[[comparison_name]]$table
  
  
  sig <- tmp[
    !is.na(tmp$p.adjust) &
      tmp$p.adjust < 0.05,
    ,
    drop = FALSE
  ]
  
  
  primary_summary <- rbind(
    primary_summary,
    
    data.frame(
      comparison = comparison_name,
      significant_pathways = nrow(sig),
      APP23_up_pathways = sum(
        sig$NES > 0,
        na.rm = TRUE
      ),
      APP23_down_pathways = sum(
        sig$NES < 0,
        na.rm = TRUE
      )
    )
  )
  
}


print(primary_summary)


write.csv(
  primary_summary,
  "RESULTS/GSEA/GSEA_PRIMARY_summary.csv",
  row.names = FALSE
)


# =============================================================================
# 124. PRINT TOP SIGNIFICANT PATHWAYS FOR EACH COMPARISON
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("TOP PRIMARY GSEA PATHWAYS\n")
cat("####################################################################\n")


for (comparison_name in names(GSEA_primary)) {
  
  tmp <- GSEA_primary[[comparison_name]]$table
  
  
  sig <- tmp[
    !is.na(tmp$p.adjust) &
      tmp$p.adjust < 0.05,
    ,
    drop = FALSE
  ]
  
  
  sig <- sig[
    order(
      sig$p.adjust
    ),
    ,
    drop = FALSE
  ]
  
  
  cat("\n============================================================\n")
  cat(comparison_name, "\n")
  cat("============================================================\n")
  
  
  if (nrow(sig) == 0) {
    
    cat("No pathways significant at FDR < 0.05.\n")
    
  } else {
    
    print(
      head(
        sig[
          ,
          c(
            "ID",
            "Description",
            "setSize",
            "enrichmentScore",
            "NES",
            "pvalue",
            "p.adjust"
          )
        ],
        15
      ),
      row.names = FALSE
    )
    
  }
  
}


# =============================================================================
# 125. TOP PATHWAYS FOR VISUALIZATION
#
# Select the strongest significant pathways according to FDR and |NES|.
# =============================================================================

plot_GSEA_dotplot <- function(
    gsea_table,
    comparison_name,
    top_n = 15
) {
  
  sig <- gsea_table[
    !is.na(gsea_table$p.adjust) &
      gsea_table$p.adjust < 0.05,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(sig) == 0) {
    
    cat(
      "\nNo significant pathways to plot for",
      comparison_name,
      "\n"
    )
    
    return(NULL)
  }
  
  
  sig <- sig[
    order(
      sig$p.adjust,
      -abs(sig$NES)
    ),
    ,
    drop = FALSE
  ]
  
  
  sig <- head(
    sig,
    top_n
  )
  
  
  sig$Description <- factor(
    sig$Description,
    levels = rev(
      sig$Description
    )
  )
  
  
  p <- ggplot(
    sig,
    aes(
      x = NES,
      y = Description,
      size = setSize,
      fill = -log10(p.adjust)
    )
  ) +
    
    geom_point(
      shape = 21
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    
    labs(
      title = paste(
        "GO Biological Process GSEA —",
        comparison_name
      ),
      subtitle = "Primary analysis: APP and Thy1 excluded",
      x = "Normalized enrichment score (NES)",
      y = NULL,
      size = "Gene-set size",
      fill = "-log10(FDR)"
    ) +
    
    theme_classic()
  
  
  return(p)
}


# =============================================================================
# 126. CREATE GSEA DOT PLOTS
# =============================================================================

for (comparison_name in names(GSEA_primary)) {
  
  p <- plot_GSEA_dotplot(
    GSEA_primary[[comparison_name]]$table,
    comparison_name,
    top_n = 15
  )
  
  
  if (!is.null(p)) {
    
    safe_name <- gsub(
      "__",
      "_",
      comparison_name,
      fixed = TRUE
    )
    
    
    ggsave(
      filename = file.path(
        "FIGURES/GSEA",
        paste0(
          "GSEA_GO_BP_",
          safe_name,
          ".pdf"
        )
      ),
      plot = p,
      width = 10,
      height = 7
    )
    
    
    print(p)
    
  }
  
}


# =============================================================================
# 127. PATHWAY-LEVEL REMODELING SUMMARY
#
# We use significant GO pathways as a pathway-level measure.
#
# IMPORTANT:
# This complements rather than replaces the gene-level DE analysis.
# =============================================================================

pathway_remodeling <- primary_summary


comparison_parts <- strsplit(
  pathway_remodeling$comparison,
  "__",
  fixed = TRUE
)


pathway_remodeling$age <- sapply(
  comparison_parts,
  `[`,
  1
)


pathway_remodeling$subtype <- sapply(
  comparison_parts,
  `[`,
  2
)


pathway_remodeling <- pathway_remodeling[
  ,
  c(
    "age",
    "subtype",
    "comparison",
    "significant_pathways",
    "APP23_up_pathways",
    "APP23_down_pathways"
  )
]


cat("\n\n")
cat("####################################################################\n")
cat("PATHWAY-LEVEL REMODELING SUMMARY\n")
cat("####################################################################\n")


print(pathway_remodeling)


write.csv(
  pathway_remodeling,
  "RESULTS/GSEA/pathway_level_remodeling_summary.csv",
  row.names = FALSE
)


# =============================================================================
# 128. PLOT NUMBER OF SIGNIFICANT PATHWAYS
# =============================================================================

p_pathway_count <- ggplot(
  pathway_remodeling,
  aes(
    x = subtype,
    y = significant_pathways,
    fill = age
  )
) +
  
  geom_col(
    position = "dodge"
  ) +
  
  labs(
    title = "APP23-associated pathway-level transcriptional remodeling",
    subtitle = "Ranked GO Biological Process GSEA",
    x = "Neuronal subtype",
    y = "Significant pathways (FDR < 0.05)"
  ) +
  
  theme_classic()


ggsave(
  "FIGURES/GSEA/GSEA_significant_pathway_counts.pdf",
  plot = p_pathway_count,
  width = 8,
  height = 6
)


print(p_pathway_count)


# =============================================================================
# 129. SAVE PRIMARY GSEA OBJECT
# =============================================================================

saveRDS(
  GSEA_primary,
  "RESULTS/GSEA/GSEA_GO_BP_PRIMARY_complete_objects.rds"
)


# =============================================================================
# 130. SENSITIVITY GSEA
#
# This version retains APP and Thy1.
#
# Purpose:
# Determine whether removing model-associated genes materially changes
# pathway-level conclusions.
# =============================================================================

GSEA_sensitivity <- list()


for (comparison_name in names(DE_results)) {
  
  cat(
    "\n\nStarting SENSITIVITY GSEA:",
    comparison_name,
    "\n"
  )
  
  
  GSEA_sensitivity[[comparison_name]] <- run_GO_GSEA(
    result_table = DE_results[[comparison_name]]$results,
    comparison_name = comparison_name,
    remove_model_features = FALSE
  )
  
}


# =============================================================================
# 131. COMBINE SENSITIVITY RESULTS
# =============================================================================

sensitivity_tables <- lapply(
  GSEA_sensitivity,
  function(x) x$table
)


sensitivity_tables <- sensitivity_tables[
  sapply(
    sensitivity_tables,
    nrow
  ) > 0
]


GSEA_sensitivity_combined <- do.call(
  rbind,
  sensitivity_tables
)


rownames(
  GSEA_sensitivity_combined
) <- NULL


write.csv(
  GSEA_sensitivity_combined,
  "RESULTS/GSEA/GSEA_GO_BP_SENSITIVITY_all_genes.csv",
  row.names = FALSE
)


# =============================================================================
# 132. SENSITIVITY SUMMARY
# =============================================================================

sensitivity_summary <- data.frame()


for (comparison_name in names(GSEA_sensitivity)) {
  
  tmp <- GSEA_sensitivity[[comparison_name]]$table
  
  
  sig <- tmp[
    !is.na(tmp$p.adjust) &
      tmp$p.adjust < 0.05,
    ,
    drop = FALSE
  ]
  
  
  sensitivity_summary <- rbind(
    sensitivity_summary,
    
    data.frame(
      comparison = comparison_name,
      significant_pathways = nrow(sig),
      APP23_up_pathways = sum(
        sig$NES > 0,
        na.rm = TRUE
      ),
      APP23_down_pathways = sum(
        sig$NES < 0,
        na.rm = TRUE
      )
    )
  )
  
}


cat("\n\n")
cat("####################################################################\n")
cat("SENSITIVITY GSEA SUMMARY — APP + Thy1 RETAINED\n")
cat("####################################################################\n")


print(sensitivity_summary)


write.csv(
  sensitivity_summary,
  "RESULTS/GSEA/GSEA_SENSITIVITY_summary.csv",
  row.names = FALSE
)


# =============================================================================
# 133. SAVE SENSITIVITY OBJECT
# =============================================================================

saveRDS(
  GSEA_sensitivity,
  "RESULTS/GSEA/GSEA_GO_BP_SENSITIVITY_complete_objects.rds"
)


# =============================================================================
# 134. FINAL GSEA OUTPUT
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("GSEA COMPLETED\n")
cat("####################################################################\n")


cat("\nPRIMARY analysis:\n")
cat("APP and Thy1 excluded from ranking.\n")
cat("Endogenous App retained.\n\n")


cat("PRIMARY pathway summary:\n")
print(primary_summary)


cat("\nSENSITIVITY pathway summary:\n")
print(sensitivity_summary)


cat("\nPathway-level remodeling:\n")
print(pathway_remodeling)


cat("\n####################################################################\n")
cat("STOP HERE BEFORE TF ACTIVITY ANALYSIS\n")
cat("####################################################################\n")

###############################################################################
# 135. EXTRACT BIOLOGICALLY RELEVANT GSEA RESULTS
###############################################################################


cat("\n\n")
cat("####################################################################\n")
cat("SIGNIFICANT PRIMARY GSEA PATHWAYS — FULL TABLE\n")
cat("####################################################################\n")


# =============================================================================
# 135A. FUNCTION TO PRINT SIGNIFICANT PATHWAYS
# =============================================================================

print_sig_gsea <- function(comparison_name) {
  
  x <- GSEA_primary[[comparison_name]]$table
  
  
  sig <- x[
    !is.na(x$p.adjust) &
      x$p.adjust < 0.05,
    ,
    drop = FALSE
  ]
  
  
  sig <- sig[
    order(
      sig$NES,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  
  cat("\n============================================================\n")
  cat(comparison_name, "\n")
  cat("============================================================\n")
  
  
  if (nrow(sig) == 0) {
    
    cat("No significant pathways.\n")
    
    return(invisible(NULL))
  }
  
  
  print(
    sig[
      ,
      c(
        "ID",
        "Description",
        "setSize",
        "NES",
        "pvalue",
        "p.adjust",
        "core_enrichment"
      )
    ],
    row.names = FALSE
  )
  
}


# =============================================================================
# 135B. SIX-MONTH DG
# =============================================================================

print_sig_gsea(
  "6_month__DG"
)


# =============================================================================
# 135C. TWENTY-FOUR-MONTH DG
# =============================================================================

print_sig_gsea(
  "24_month__DG"
)


# =============================================================================
# 135D. TWENTY-FOUR-MONTH CA1-LIKE
# =============================================================================

print_sig_gsea(
  "24_month__CA1_like"
)


# =============================================================================
# 136. CREATE COMBINED SIGNIFICANT-PATHWAY TABLE
# =============================================================================

key_comparisons <- c(
  "6_month__DG",
  "24_month__DG",
  "24_month__CA1_like"
)


key_GSEA <- list()


for (comparison_name in key_comparisons) {
  
  tmp <- GSEA_primary[[comparison_name]]$table
  
  
  tmp <- tmp[
    !is.na(tmp$p.adjust) &
      tmp$p.adjust < 0.05,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(tmp) == 0) {
    next
  }
  
  
  tmp$comparison <- comparison_name
  
  
  tmp$direction <- ifelse(
    tmp$NES > 0,
    "Enriched_in_APP23",
    "Depleted_in_APP23"
  )
  
  
  key_GSEA[[comparison_name]] <- tmp
}


key_GSEA_combined <- do.call(
  rbind,
  key_GSEA
)


rownames(key_GSEA_combined) <- NULL


write.csv(
  key_GSEA_combined,
  "RESULTS/GSEA/key_significant_GSEA_pathways.csv",
  row.names = FALSE
)


# =============================================================================
# 137. SIMPLE INTERPRETATION TABLE
# =============================================================================

interpretation_table <- key_GSEA_combined[
  ,
  c(
    "comparison",
    "ID",
    "Description",
    "NES",
    "p.adjust",
    "direction"
  )
]


interpretation_table <- interpretation_table[
  order(
    interpretation_table$comparison,
    interpretation_table$p.adjust
  ),
  ,
  drop = FALSE
]


cat("\n\n")
cat("####################################################################\n")
cat("GSEA INTERPRETATION TABLE\n")
cat("####################################################################\n")


print(
  interpretation_table,
  row.names = FALSE
)


write.csv(
  interpretation_table,
  "RESULTS/GSEA/GSEA_interpretation_table.csv",
  row.names = FALSE
)


cat("\n####################################################################\n")
cat("GSEA PATHWAY EXTRACTION COMPLETED\n")
cat("####################################################################\n")

###############################################################################
# 138. TF ACTIVITY — PACKAGE CHECK
###############################################################################

cat("\n")
cat("####################################################################\n")
cat("TF ACTIVITY PACKAGE CHECK\n")
cat("####################################################################\n")


packages_to_check <- c(
  "decoupleR",
  "dorothea",
  "OmnipathR"
)


for (pkg in packages_to_check) {
  
  cat("\n", pkg, ":\n", sep = "")
  
  if (requireNamespace(pkg, quietly = TRUE)) {
    
    cat(
      "INSTALLED — version ",
      as.character(packageVersion(pkg)),
      "\n",
      sep = ""
    )
    
  } else {
    
    cat("NOT INSTALLED\n")
  }
}


cat("\nR version:\n")
print(R.version.string)


cat("\nAvailable decoupleR functions related to ULM:\n")

if (requireNamespace("decoupleR", quietly = TRUE)) {
  
  print(
    grep(
      "ulm",
      getNamespaceExports("decoupleR"),
      value = TRUE,
      ignore.case = TRUE
    )
  )
  
}


cat("\n####################################################################\n")
cat("STOP HERE\n")
cat("####################################################################\n")

###############################################################################
# GSE141044
# TRANSCRIPTION FACTOR ACTIVITY ANALYSIS
#
# decoupleR 2.16.0 + dorothea 1.22.0
#
# Strategy:
# 1. Use mouse DoRothEA regulons
# 2. Keep high-confidence regulons A/B/C
# 3. Calculate TF activity with ULM
# 4. Use pseudobulk logCPM as input
# 5. Compare APP23 vs WT within age + neuronal subtype
#
# IMPORTANT:
# Biological replicate = mouse
###############################################################################


# =============================================================================
# 139. LOAD PACKAGES
# =============================================================================

library(decoupleR)
library(dorothea)
library(edgeR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)


set.seed(100)


# =============================================================================
# 140. OUTPUT DIRECTORIES
# =============================================================================

dir.create(
  "RESULTS/TF_ACTIVITY",
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  "FIGURES/TF_ACTIVITY",
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 141. LOAD MOUSE DOROTHEA REGULON
# =============================================================================

cat("\n============================================================\n")
cat("LOADING MOUSE DOROTHEA REGULON\n")
cat("============================================================\n")


data(
  "dorothea_mm",
  package = "dorothea"
)


cat("\nDoRothEA mouse network dimensions:\n")
print(dim(dorothea_mm))


cat("\nColumns:\n")
print(colnames(dorothea_mm))


cat("\nConfidence levels:\n")
print(table(dorothea_mm$confidence))


cat("\nFirst rows:\n")
print(head(dorothea_mm))


# =============================================================================
# =============================================================================
# 142. KEEP HIGH-CONFIDENCE A/B/C REGULONS
# =============================================================================

dorothea_mouse <- dorothea_mm %>%
  
  dplyr::filter(
    confidence %in% c(
      "A",
      "B",
      "C"
    )
  ) %>%
  
  dplyr::select(
    tf,
    target,
    mor,
    confidence
  )


cat("\n============================================================\n")
cat("HIGH-CONFIDENCE DOROTHEA NETWORK\n")
cat("============================================================\n")


cat("\nNumber of TF-target interactions:\n")
print(nrow(dorothea_mouse))


cat("\nNumber of TFs:\n")
print(length(unique(dorothea_mouse$tf)))


cat("\nConfidence distribution:\n")
print(table(dorothea_mouse$confidence))

# =============================================================================
# 143. VERIFY PSEUDOBULK OBJECTS EXIST
# =============================================================================

required_objects <- c(
  "pseudobulk_counts",
  "pb_metadata"
)


missing_objects <- required_objects[
  !sapply(
    required_objects,
    exists
  )
]


if (length(missing_objects) > 0) {
  
  stop(
    paste0(
      "Missing required objects: ",
      paste(
        missing_objects,
        collapse = ", "
      ),
      "\nLoad the pseudobulk analysis objects first."
    )
  )
}


cat("\nPseudobulk count dimensions:\n")
print(dim(pseudobulk_counts))


cat("\nPseudobulk metadata dimensions:\n")
print(dim(pb_metadata))


# =============================================================================
# 144. CREATE LOGCPM MATRIX
#
# ULM should operate on continuous expression values rather than raw counts.
# We calculate TMM-normalized logCPM across all 44 pseudobulk libraries.
# =============================================================================

cat("\n============================================================\n")
cat("CREATING PSEUDOBULK LOGCPM MATRIX\n")
cat("============================================================\n")


tf_y <- DGEList(
  counts = pseudobulk_counts
)


tf_y <- calcNormFactors(
  tf_y,
  method = "TMM"
)


tf_logCPM <- cpm(
  tf_y,
  log = TRUE,
  prior.count = 2
)


cat("\nlogCPM dimensions:\n")
print(dim(tf_logCPM))


cat("\nRange of logCPM values:\n")
print(range(tf_logCPM))


# =============================================================================
# 145. REMOVE MODEL-CONSTRUCT FEATURES FROM TF INPUT
#
# APP = transgene-associated feature
# Thy1 = model/promoter-associated feature
#
# Endogenous App remains.
# =============================================================================

model_features_to_remove <- intersect(
  c(
    "APP",
    "Thy1"
  ),
  rownames(tf_logCPM)
)


cat("\nModel-associated features removed:\n")
print(model_features_to_remove)


tf_logCPM_primary <- tf_logCPM[
  !rownames(tf_logCPM) %in% model_features_to_remove,
  ,
  drop = FALSE
]


cat("\nGenes remaining:\n")
print(nrow(tf_logCPM_primary))


# =============================================================================
# 146. CHECK NETWORK OVERLAP
# =============================================================================

network_targets <- unique(
  dorothea_mouse$target
)


gene_overlap <- intersect(
  rownames(tf_logCPM_primary),
  network_targets
)


cat("\n============================================================\n")
cat("DOROTHEA TARGET OVERLAP\n")
cat("============================================================\n")


cat("\nExpression genes:\n")
print(nrow(tf_logCPM_primary))


cat("\nUnique DoRothEA targets:\n")
print(length(network_targets))


cat("\nOverlapping targets:\n")
print(length(gene_overlap))


cat("\nPercentage of DoRothEA targets represented:\n")

print(
  round(
    100 *
      length(gene_overlap) /
      length(network_targets),
    2
  )
)


if (length(gene_overlap) < 500) {
  
  warning(
    "Low overlap between expression genes and DoRothEA targets. Check gene naming."
  )
}


# =============================================================================
# 147. FORMAT NETWORK FOR decoupleR 2.16
#
# run_ulm expects source / target / mor.
# =============================================================================

tf_network <- dorothea_mouse %>%
  
  transmute(
    source = tf,
    target = target,
    mor = mor
  )


cat("\nNetwork passed to decoupleR:\n")
print(head(tf_network))


# =============================================================================
# 148. RUN ULM TF ACTIVITY
# =============================================================================

cat("\n============================================================\n")
cat("RUNNING ULM TF ACTIVITY\n")
cat("============================================================\n")


tf_activity_long <- decoupleR::run_ulm(
  mat = tf_logCPM_primary,
  network = tf_network,
  .source = source,
  .target = target,
  .mor = mor,
  minsize = 5
)


cat("\nULM completed.\n")


cat("\nOutput columns:\n")
print(colnames(tf_activity_long))


cat("\nOutput dimensions:\n")
print(dim(tf_activity_long))


cat("\nFirst rows:\n")
print(head(tf_activity_long))


# =============================================================================
# 149. VERIFY ULM OUTPUT STRUCTURE
#
# decoupleR 2.16 normally returns:
# statistic
# source
# condition
# score
# =============================================================================

required_ulm_cols <- c(
  "source",
  "condition",
  "score"
)


if (
  !all(
    required_ulm_cols %in% colnames(tf_activity_long)
  )
) {
  
  stop(
    paste0(
      "Unexpected run_ulm output columns: ",
      paste(
        colnames(tf_activity_long),
        collapse = ", "
      )
    )
  )
}


# =============================================================================
# 150. SAVE COMPLETE TF ACTIVITY TABLE
# =============================================================================

write.csv(
  tf_activity_long,
  "RESULTS/TF_ACTIVITY/ULM_all_pseudobulk_TF_activities.csv",
  row.names = FALSE
)


# =============================================================================
# 151. CONVERT TF ACTIVITIES TO MATRIX
#
# Rows = TFs
# Columns = pseudobulk samples
# =============================================================================
# =============================================================================
# 151. CONVERT TF ACTIVITIES TO MATRIX
#
# Rows = TFs
# Columns = pseudobulk samples
# =============================================================================

tf_activity_matrix <- tf_activity_long %>%
  
  dplyr::select(
    source,
    condition,
    score
  ) %>%
  
  tidyr::pivot_wider(
    names_from = condition,
    values_from = score
  )


tf_names <- tf_activity_matrix$source


tf_activity_matrix <- as.data.frame(
  tf_activity_matrix
)


rownames(tf_activity_matrix) <- tf_names


tf_activity_matrix$source <- NULL


tf_activity_matrix <- as.matrix(
  tf_activity_matrix
)


storage.mode(tf_activity_matrix) <- "numeric"


cat("\n============================================================\n")
cat("TF ACTIVITY MATRIX\n")
cat("============================================================\n")


cat("\nDimensions (TFs x pseudobulk samples):\n")
print(dim(tf_activity_matrix))
# =============================================================================
# 152. ALIGN TF MATRIX WITH PSEUDOBULK METADATA
# =============================================================================

common_samples <- intersect(
  colnames(tf_activity_matrix),
  rownames(pb_metadata)
)


cat("\nSamples common to TF matrix and metadata:\n")
print(length(common_samples))


if (
  length(common_samples) != nrow(pb_metadata)
) {
  
  warning(
    "Not every pseudobulk sample is represented in TF activity matrix."
  )
}


tf_activity_matrix <- tf_activity_matrix[
  ,
  common_samples,
  drop = FALSE
]


tf_metadata <- pb_metadata[
  common_samples,
  ,
  drop = FALSE
]


stopifnot(
  identical(
    colnames(tf_activity_matrix),
    rownames(tf_metadata)
  )
)


# =============================================================================
# 153. SAVE TF ACTIVITY MATRIX
# =============================================================================

write.csv(
  tf_activity_matrix,
  "RESULTS/TF_ACTIVITY/ULM_TF_activity_matrix.csv"
)


# =============================================================================
# 154. FUNCTION TO TEST TF ACTIVITY APP23 VS WT
#
# For each age + neuronal subtype:
#
# 6 months = 3 WT vs 3 APP23
# 24 months = 2 WT vs 3 APP23
#
# Because sample sizes are small, we use limma.
# =============================================================================

if (!requireNamespace("limma", quietly = TRUE)) {
  
  BiocManager::install(
    "limma",
    ask = FALSE,
    update = FALSE
  )
  
}


library(limma)


run_TF_comparison <- function(
    activity_matrix,
    metadata,
    subtype_name,
    age_name
) {
  
  cat("\n\n")
  cat("============================================================\n")
  cat("TF ACTIVITY COMPARISON\n")
  cat("Subtype:", subtype_name, "\n")
  cat("Age:", age_name, "\n")
  cat("APP23 vs WT\n")
  cat("============================================================\n")
  
  
  # ---------------------------------------------------------------------------
  # Select samples
  # ---------------------------------------------------------------------------
  
  keep <- (
    metadata$broad_subtype == subtype_name &
      metadata$age == age_name
  )
  
  
  meta_sub <- metadata[
    keep,
    ,
    drop = FALSE
  ]
  
  
  activity_sub <- activity_matrix[
    ,
    rownames(meta_sub),
    drop = FALSE
  ]
  
  
  cat("\nSamples:\n")
  
  print(
    meta_sub[
      ,
      c(
        "sample_id",
        "genotype"
      )
    ]
  )
  
  
  cat("\nReplicates:\n")
  print(table(meta_sub$genotype))
  
  
  # ---------------------------------------------------------------------------
  # Genotype
  # ---------------------------------------------------------------------------
  
  genotype <- factor(
    meta_sub$genotype,
    levels = c(
      "WT",
      "APP23"
    )
  )
  
  
  design <- model.matrix(
    ~ genotype
  )
  
  
  rownames(design) <- rownames(
    meta_sub
  )
  
  
  # ---------------------------------------------------------------------------
  # limma
  # ---------------------------------------------------------------------------
  
  fit <- lmFit(
    activity_sub,
    design
  )
  
  
  fit <- eBayes(
    fit,
    robust = TRUE
  )
  
  
  results <- topTable(
    fit,
    coef = "genotypeAPP23",
    number = Inf,
    sort.by = "P"
  )
  
  
  results$TF <- rownames(
    results
  )
  
  
  results$age <- age_name
  
  results$subtype <- subtype_name
  
  
  # ---------------------------------------------------------------------------
  # Interpretation
  #
  # Positive logFC:
  # greater inferred TF activity in APP23
  #
  # Negative:
  # lower inferred activity in APP23
  # ---------------------------------------------------------------------------
  
  results$direction <- ifelse(
    results$adj.P.Val < 0.05 &
      results$logFC > 0,
    "Higher_activity_APP23",
    
    ifelse(
      results$adj.P.Val < 0.05 &
        results$logFC < 0,
      "Lower_activity_APP23",
      "NS"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Save
  # ---------------------------------------------------------------------------
  
  output_file <- paste0(
    "TF_APP23_vs_WT_",
    age_name,
    "_",
    subtype_name,
    ".csv"
  )
  
  
  write.csv(
    results,
    file.path(
      "RESULTS/TF_ACTIVITY",
      output_file
    ),
    row.names = FALSE
  )
  
  
  cat("\nTFs significant at FDR < 0.05:\n")
  
  print(
    sum(
      results$adj.P.Val < 0.05
    )
  )
  
  
  return(results)
}


# =============================================================================
# 155. RUN ALL 8 TF ACTIVITY COMPARISONS
# =============================================================================

TF_results <- list()


for (current_age in ages) {
  
  for (current_subtype in broad_subtypes) {
    
    comparison_name <- paste(
      current_age,
      current_subtype,
      sep = "__"
    )
    
    
    TF_results[[comparison_name]] <- run_TF_comparison(
      activity_matrix = tf_activity_matrix,
      metadata = tf_metadata,
      subtype_name = current_subtype,
      age_name = current_age
    )
    
  }
  
}


# =============================================================================
# 156. COMBINE TF RESULTS
# =============================================================================

TF_combined <- do.call(
  rbind,
  lapply(
    names(TF_results),
    function(comparison_name) {
      
      tmp <- TF_results[[comparison_name]]
      
      tmp$comparison <- comparison_name
      
      tmp
    }
  )
)


rownames(TF_combined) <- NULL


write.csv(
  TF_combined,
  "RESULTS/TF_ACTIVITY/TF_activity_all_8_comparisons.csv",
  row.names = FALSE
)


# =============================================================================
# 157. SIGNIFICANT TF ACTIVITIES
# =============================================================================

TF_significant <- TF_combined[
  !is.na(TF_combined$adj.P.Val) &
    TF_combined$adj.P.Val < 0.05,
  ,
  drop = FALSE
]


TF_significant <- TF_significant[
  order(
    TF_significant$comparison,
    TF_significant$adj.P.Val
  ),
  ,
  drop = FALSE
]


cat("\n\n")
cat("####################################################################\n")
cat("SIGNIFICANT TF ACTIVITY CHANGES\n")
cat("####################################################################\n")


if (nrow(TF_significant) == 0) {
  
  cat("\nNo TF activities significant at FDR < 0.05.\n")
  
} else {
  
  print(
    TF_significant[
      ,
      c(
        "comparison",
        "TF",
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val",
        "direction"
      )
    ],
    row.names = FALSE
  )
  
}


write.csv(
  TF_significant,
  "RESULTS/TF_ACTIVITY/significant_TF_activity_FDR05.csv",
  row.names = FALSE
)


# =============================================================================
# 158. TF ACTIVITY SUMMARY
# =============================================================================

TF_summary <- data.frame()


for (comparison_name in names(TF_results)) {
  
  tmp <- TF_results[[comparison_name]]
  
  
  parts <- strsplit(
    comparison_name,
    "__",
    fixed = TRUE
  )[[1]]
  
  
  sig <- tmp[
    tmp$adj.P.Val < 0.05,
    ,
    drop = FALSE
  ]
  
  
  TF_summary <- rbind(
    TF_summary,
    
    data.frame(
      age = parts[1],
      subtype = parts[2],
      comparison = comparison_name,
      significant_TFs = nrow(sig),
      higher_in_APP23 = sum(
        sig$logFC > 0,
        na.rm = TRUE
      ),
      lower_in_APP23 = sum(
        sig$logFC < 0,
        na.rm = TRUE
      )
    )
  )
  
}


cat("\n\n")
cat("####################################################################\n")
cat("TF ACTIVITY SUMMARY\n")
cat("####################################################################\n")


print(TF_summary)


write.csv(
  TF_summary,
  "RESULTS/TF_ACTIVITY/TF_activity_summary.csv",
  row.names = FALSE
)


# =============================================================================
# 159. TOP 10 TFs PER COMPARISON
#
# Even if nothing reaches FDR < 0.05, this gives ranked regulatory signals.
# These must NOT be called significant unless FDR passes threshold.
# =============================================================================

TF_top10 <- list()


for (comparison_name in names(TF_results)) {
  
  tmp <- TF_results[[comparison_name]]
  
  
  tmp <- tmp[
    order(
      tmp$P.Value
    ),
    ,
    drop = FALSE
  ]
  
  
  tmp <- head(
    tmp,
    10
  )
  
  
  tmp$comparison <- comparison_name
  
  
  TF_top10[[comparison_name]] <- tmp
  
  
  cat("\n============================================================\n")
  cat("TOP TFs:", comparison_name, "\n")
  cat("============================================================\n")
  
  
  print(
    tmp[
      ,
      c(
        "TF",
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val"
      )
    ],
    row.names = FALSE
  )
  
}


TF_top10_combined <- do.call(
  rbind,
  TF_top10
)


rownames(TF_top10_combined) <- NULL


write.csv(
  TF_top10_combined,
  "RESULTS/TF_ACTIVITY/top10_TFs_each_comparison.csv",
  row.names = FALSE
)


# =============================================================================
# 160. SELECT TFs FOR HEATMAP
#
# Take top-ranked TFs across all comparisons.
# =============================================================================

top_TF_names <- unique(
  TF_top10_combined$TF
)


top_TF_names <- intersect(
  top_TF_names,
  rownames(tf_activity_matrix)
)


cat("\nNumber of TFs selected for heatmap:\n")
print(length(top_TF_names))


# =============================================================================
# 161. Z-SCORE TF ACTIVITY ACROSS PSEUDOBULK SAMPLES
# =============================================================================

heatmap_matrix <- tf_activity_matrix[
  top_TF_names,
  ,
  drop = FALSE
]


heatmap_matrix_z <- t(
  scale(
    t(
      heatmap_matrix
    )
  )
)


heatmap_matrix_z[
  !is.finite(
    heatmap_matrix_z
  )
] <- 0


# =============================================================================
# 162. SAMPLE ANNOTATION FOR HEATMAP
# =============================================================================

heatmap_annotation <- data.frame(
  Age = tf_metadata$age,
  Genotype = tf_metadata$genotype,
  Subtype = tf_metadata$broad_subtype
)


rownames(
  heatmap_annotation
) <- rownames(
  tf_metadata
)


# =============================================================================
# 163. TF ACTIVITY HEATMAP
# =============================================================================

pdf(
  "FIGURES/TF_ACTIVITY/TF_activity_top_regulators_heatmap.pdf",
  width = 14,
  height = 12
)


pheatmap(
  heatmap_matrix_z,
  annotation_col = heatmap_annotation,
  show_colnames = FALSE,
  fontsize_row = 7,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Inferred transcription factor activity"
)


dev.off()


# =============================================================================
# 164. PLOT SIGNIFICANT TFs / TOP-RANKED TFs
#
# Use combined differential TF activity estimates.
# =============================================================================

plot_tf <- TF_combined


plot_tf$minus_log10_FDR <- -log10(
  pmax(
    plot_tf$adj.P.Val,
    1e-300
  )
)


# Select top 5 TFs per comparison for labels

tf_labels <- plot_tf %>%
  
  group_by(comparison) %>%
  
  arrange(P.Value, .by_group = TRUE) %>%
  
  slice_head(n = 5) %>%
  
  ungroup()


for (comparison_name in unique(plot_tf$comparison)) {
  
  tmp <- plot_tf[
    plot_tf$comparison == comparison_name,
    ,
    drop = FALSE
  ]
  
  
  labels_tmp <- tf_labels[
    tf_labels$comparison == comparison_name,
    ,
    drop = FALSE
  ]
  
  
  p <- ggplot(
    tmp,
    aes(
      x = logFC,
      y = minus_log10_FDR
    )
  ) +
    
    geom_point(
      alpha = 0.6
    ) +
    
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    
    geom_text(
      data = labels_tmp,
      aes(
        label = TF
      ),
      check_overlap = TRUE,
      vjust = -0.5,
      size = 3
    ) +
    
    labs(
      title = paste(
        "TF activity —",
        comparison_name
      ),
      x = "Difference in inferred TF activity\n(APP23 - WT)",
      y = "-log10(FDR)"
    ) +
    
    theme_classic()
  
  
  safe_name <- gsub(
    "__",
    "_",
    comparison_name,
    fixed = TRUE
  )
  
  
  ggsave(
    filename = file.path(
      "FIGURES/TF_ACTIVITY",
      paste0(
        "TF_activity_",
        safe_name,
        ".pdf"
      )
    ),
    plot = p,
    width = 7,
    height = 6
  )
  
}


# =============================================================================
# 165. SAVE COMPLETE TF ANALYSIS OBJECTS
# =============================================================================

saveRDS(
  TF_results,
  "RESULTS/TF_ACTIVITY/TF_activity_results_complete.rds"
)


saveRDS(
  tf_activity_matrix,
  "RESULTS/TF_ACTIVITY/TF_activity_matrix.rds"
)


# =============================================================================
# 166. FINAL OUTPUT
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("TF ACTIVITY ANALYSIS COMPLETED\n")
cat("####################################################################\n")


cat("\nTF activity matrix dimensions:\n")
print(dim(tf_activity_matrix))


cat("\nTF activity summary:\n")
print(TF_summary)


cat("\nSignificant TFs:\n")

if (nrow(TF_significant) == 0) {
  
  cat(
    "No TFs reached FDR < 0.05.\n"
  )
  
} else {
  
  print(
    TF_significant[
      ,
      c(
        "comparison",
        "TF",
        "logFC",
        "adj.P.Val",
        "direction"
      )
    ],
    row.names = FALSE
  )
  
}


cat("\n####################################################################\n")
cat("STOP HERE BEFORE CELL-CELL COMMUNICATION ANALYSIS\n")
cat("####################################################################\n")

###############################################################################
# 167. EXTRACT TOP-RANKED TFs FOR KEY COMPARISONS
###############################################################################

key_tf_comparisons <- c(
  "6_month__DG",
  "24_month__DG",
  "24_month__CA1_like"
)


cat("\n\n")
cat("####################################################################\n")
cat("TOP-RANKED TF ACTIVITY CHANGES — KEY COMPARISONS\n")
cat("####################################################################\n")


key_TF_tables <- list()


for (comparison_name in key_tf_comparisons) {
  
  tmp <- TF_results[[comparison_name]]
  
  tmp <- tmp[
    order(
      tmp$P.Value
    ),
    ,
    drop = FALSE
  ]
  
  tmp_top <- head(
    tmp,
    15
  )
  
  tmp_top$comparison <- comparison_name
  
  key_TF_tables[[comparison_name]] <- tmp_top
  
  
  cat("\n============================================================\n")
  cat(comparison_name, "\n")
  cat("============================================================\n")
  
  print(
    tmp_top[
      ,
      c(
        "TF",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val"
      )
    ],
    row.names = FALSE
  )
}


key_TF_combined <- do.call(
  rbind,
  key_TF_tables
)

rownames(
  key_TF_combined
) <- NULL


write.csv(
  key_TF_combined,
  "RESULTS/TF_ACTIVITY/key_top15_TF_rankings.csv",
  row.names = FALSE
)


cat("\n####################################################################\n")
cat("TOP-RANKED TF EXTRACTION COMPLETED\n")
cat("####################################################################\n")

###############################################################################
# 168. CELLCHAT PACKAGE CHECK
###############################################################################

cat("\n")
cat("####################################################################\n")
cat("CELLCHAT PACKAGE CHECK\n")
cat("####################################################################\n")


packages_to_check <- c(
  "CellChat",
  "future",
  "patchwork"
)


for (pkg in packages_to_check) {
  
  cat("\n", pkg, ":\n", sep = "")
  
  if (requireNamespace(pkg, quietly = TRUE)) {
    
    cat(
      "INSTALLED — version ",
      as.character(packageVersion(pkg)),
      "\n",
      sep = ""
    )
    
  } else {
    
    cat("NOT INSTALLED\n")
  }
}


cat("\nR version:\n")
print(R.version.string)


if (requireNamespace("CellChat", quietly = TRUE)) {
  
  cat("\nSelected CellChat functions available:\n")
  
  functions_to_check <- c(
    "createCellChat",
    "subsetData",
    "identifyOverExpressedGenes",
    "identifyOverExpressedInteractions",
    "computeCommunProb",
    "filterCommunication",
    "computeCommunProbPathway",
    "aggregateNet",
    "mergeCellChat",
    "compareInteractions",
    "netVisual_diffInteraction"
  )
  
  
  available <- functions_to_check[
    functions_to_check %in%
      getNamespaceExports("CellChat")
  ]
  
  
  print(available)
}


cat("\n####################################################################\n")
cat("STOP HERE\n")
cat("####################################################################\n")

###############################################################################
# 169. INSTALL CELLCHAT
###############################################################################

cat("\n")
cat("####################################################################\n")
cat("INSTALLING CELLCHAT\n")
cat("####################################################################\n")


# -----------------------------------------------------------------------------
# 1. Make sure installation helpers are available
# -----------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}


# -----------------------------------------------------------------------------
# 2. Install core dependencies
# -----------------------------------------------------------------------------

cran_packages <- c(
  "NMF",
  "circlize",
  "future",
  "future.apply",
  "patchwork"
)


for (pkg in cran_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    cat("\nInstalling ", pkg, "...\n", sep = "")
    
    install.packages(
      pkg,
      dependencies = TRUE
    )
  }
}


# -----------------------------------------------------------------------------
# 3. Install ComplexHeatmap from Bioconductor
# -----------------------------------------------------------------------------

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  
  BiocManager::install(
    "ComplexHeatmap",
    ask = FALSE,
    update = FALSE
  )
}


# -----------------------------------------------------------------------------
# 4. Install CellChat from the official GitHub repository
# -----------------------------------------------------------------------------

if (!requireNamespace("CellChat", quietly = TRUE)) {
  
  cat("\nInstalling CellChat from jinworks/CellChat...\n")
  
  Sys.setenv(
    R_REMOTES_NO_ERRORS_FROM_WARNINGS = "TRUE"
  )
  
  
  devtools::install_github(
    "jinworks/CellChat",
    dependencies = TRUE,
    upgrade = "never"
  )
}


# -----------------------------------------------------------------------------
# 5. Verify installation
# -----------------------------------------------------------------------------

cat("\n\n")
cat("####################################################################\n")
cat("CELLCHAT INSTALLATION CHECK\n")
cat("####################################################################\n")


if (requireNamespace("CellChat", quietly = TRUE)) {
  
  cat(
    "\nCellChat INSTALLED — version ",
    as.character(packageVersion("CellChat")),
    "\n",
    sep = ""
  )
  
} else {
  
  cat("\nCellChat installation FAILED.\n")
}


cat("\nNMF:\n")
if (requireNamespace("NMF", quietly = TRUE)) {
  print(packageVersion("NMF"))
}


cat("\ncirclize:\n")
if (requireNamespace("circlize", quietly = TRUE)) {
  print(packageVersion("circlize"))
}


cat("\nComplexHeatmap:\n")
if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  print(packageVersion("ComplexHeatmap"))
}


cat("\n####################################################################\n")
cat("STOP HERE\n")
cat("####################################################################\n")

###############################################################################
# GSE141044
# CELLCHAT ANALYSIS
#
# Primary comparison:
# 24-month WT vs 24-month APP23
#
# Cell groups:
# DG
# CA1_like
# CA3
# Inhibitory
#
# CellChat version:
# 2.2.0.9001
###############################################################################


# =============================================================================
# 169. LOAD PACKAGES
# =============================================================================

library(CellChat)
library(Seurat)
library(Matrix)
library(dplyr)
library(patchwork)
library(ggplot2)

set.seed(100)


# =============================================================================
# 170. OUTPUT DIRECTORIES
# =============================================================================

dir.create(
  "RESULTS/CELLCHAT",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "FIGURES/CELLCHAT",
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 171. VERIFY REQUIRED METADATA
# =============================================================================

required_metadata <- c(
  "age",
  "genotype",
  "broad_subtype"
)

missing_metadata <- required_metadata[
  !required_metadata %in% colnames(app23@meta.data)
]

if (length(missing_metadata) > 0) {
  
  stop(
    paste0(
      "Missing metadata columns: ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
}


cat("\n============================================================\n")
cat("CELLCHAT METADATA CHECK\n")
cat("============================================================\n")

print(
  table(
    app23$age,
    app23$genotype
  )
)

cat("\nBroad neuronal populations:\n")

print(
  table(
    app23$broad_subtype
  )
)


# =============================================================================
# 172. KEEP ONLY 24-MONTH NUCLEI
# =============================================================================

cells_24m <- colnames(app23)[
  app23$age == "24_month"
]

app23_24m <- subset(
  app23,
  cells = cells_24m
)


cat("\n24-month nuclei:\n")
print(ncol(app23_24m))


cat("\n24-month genotype distribution:\n")
print(
  table(
    app23_24m$genotype
  )
)


cat("\n24-month broad subtype distribution:\n")
print(
  table(
    app23_24m$genotype,
    app23_24m$broad_subtype
  )
)


# =============================================================================
# 173. DEFINE CELL GROUP ORDER
# =============================================================================

cell_group_order <- c(
  "DG",
  "CA1_like",
  "CA3",
  "Inhibitory"
)


app23_24m$broad_subtype <- factor(
  app23_24m$broad_subtype,
  levels = cell_group_order
)


# =============================================================================
# 174. EXTRACT NORMALIZED RNA EXPRESSION
#
# CellChat uses normalized expression rather than raw counts.
#
# Seurat v5 stores normalized RNA expression in the "data" layer.
# =============================================================================

data_24m <- SeuratObject::LayerData(
  app23_24m,
  assay = "RNA",
  layer = "data"
)


cat("\nNormalized expression matrix dimensions:\n")
print(dim(data_24m))


# =============================================================================
# 175. REMOVE MODEL-CONSTRUCT FEATURES
#
# APP = transgene-associated feature
# Thy1 = construct/promoter-associated feature
#
# Endogenous App remains.
# =============================================================================

features_remove <- intersect(
  c(
    "APP",
    "Thy1"
  ),
  rownames(data_24m)
)


cat("\nRemoving model-associated features:\n")
print(features_remove)


data_24m <- data_24m[
  !rownames(data_24m) %in% features_remove,
  ,
  drop = FALSE
]


# =============================================================================
# 176. SPLIT CELLS INTO 24-MONTH WT AND APP23
# =============================================================================

wt_cells <- colnames(app23_24m)[
  app23_24m$genotype == "WT"
]

app_cells <- colnames(app23_24m)[
  app23_24m$genotype == "APP23"
]


cat("\n24-month WT nuclei:\n")
print(length(wt_cells))


cat("\n24-month APP23 nuclei:\n")
print(length(app_cells))


# =============================================================================
# 177. CREATE METADATA FOR EACH CONDITION
# =============================================================================

meta_WT <- data.frame(
  group = app23_24m$broad_subtype[wt_cells]
)

rownames(meta_WT) <- wt_cells


meta_APP23 <- data.frame(
  group = app23_24m$broad_subtype[app_cells]
)

rownames(meta_APP23) <- app_cells


meta_WT$group <- factor(
  meta_WT$group,
  levels = cell_group_order
)


meta_APP23$group <- factor(
  meta_APP23$group,
  levels = cell_group_order
)


cat("\nWT cell groups:\n")
print(table(meta_WT$group))


cat("\nAPP23 cell groups:\n")
print(table(meta_APP23$group))


# =============================================================================
# 178. CREATE CONDITION-SPECIFIC EXPRESSION MATRICES
# =============================================================================

expr_WT <- data_24m[
  ,
  wt_cells,
  drop = FALSE
]


expr_APP23 <- data_24m[
  ,
  app_cells,
  drop = FALSE
]


# =============================================================================
# 179. CREATE CELLCHAT OBJECTS
# =============================================================================

cat("\n============================================================\n")
cat("CREATING CELLCHAT OBJECTS\n")
cat("============================================================\n")


cellchat_WT <- createCellChat(
  object = expr_WT,
  meta = meta_WT,
  group.by = "group"
)


cellchat_APP23 <- createCellChat(
  object = expr_APP23,
  meta = meta_APP23,
  group.by = "group"
)


# =============================================================================
# 180. ASSIGN MOUSE CELLCHAT DATABASE
# =============================================================================

CellChatDB.use <- CellChatDB.mouse


cat("\nMouse CellChat database categories:\n")

print(
  table(
    CellChatDB.use$interaction$annotation
  )
)


cellchat_WT@DB <- CellChatDB.use
cellchat_APP23@DB <- CellChatDB.use


# =============================================================================
# =============================================================================
# =============================================================================
# 181. FUNCTION TO RUN ONE CELLCHAT DATASET
#
# CellChat 2.2.0.9001
#
# IMPORTANT:
# - projectData() is not exported in this version
# - We therefore use the normalized signaling data directly
# - computeCommunProb(raw.use = TRUE)
# =============================================================================
###############################################################################
# 181. CORRECTED CELLCHAT FUNCTION
# CellChat 2.2.0.9001
#
# NO projectData()
# NO PPI smoothing
# Uses normalized data.signaling directly
###############################################################################

run_cellchat <- function(cellchat_object, condition_name) {
  
  cat("\n\n")
  cat("####################################################################\n")
  cat("RUNNING CELLCHAT:", condition_name, "\n")
  cat("####################################################################\n")
  
  
  # ===========================================================================
  # 1. SUBSET TO SIGNALING GENES
  # ===========================================================================
  
  cellchat_object <- CellChat::subsetData(
    cellchat_object
  )
  
  cat("\nGenes after CellChat database subsetting:\n")
  print(nrow(cellchat_object@data.signaling))
  
  
  # ===========================================================================
  # 2. IDENTIFY OVEREXPRESSED SIGNALING GENES
  #
  # do.fast = FALSE because presto is not installed
  # ===========================================================================
  
  cellchat_object <- CellChat::identifyOverExpressedGenes(
    cellchat_object,
    do.fast = FALSE
  )
  
  
  # ===========================================================================
  # 3. IDENTIFY OVEREXPRESSED LIGAND-RECEPTOR INTERACTIONS
  # ===========================================================================
  
  cellchat_object <- CellChat::identifyOverExpressedInteractions(
    cellchat_object
  )
  
  
  cat("\nHighly variable ligand-receptor pairs retained:\n")
  
  if (!is.null(cellchat_object@LR$LRsig)) {
    
    print(
      nrow(
        cellchat_object@LR$LRsig
      )
    )
    
  }
  
  
  # ===========================================================================
  # IMPORTANT CHECK:
  #
  # There is intentionally NO projectData() call here.
  #
  # raw.use = TRUE tells CellChat to use object@data.signaling directly.
  # ===========================================================================
  
  
  # ===========================================================================
  # 4. COMPUTE COMMUNICATION PROBABILITY
  # ===========================================================================
  
  cat("\nComputing communication probabilities...\n")
  
  cellchat_object <- CellChat::computeCommunProb(
    cellchat_object,
    type = "triMean",
    raw.use = TRUE,
    population.size = FALSE,
    seed.use = 100
  )
  
  
  # ===========================================================================
  # 5. FILTER COMMUNICATION
  # ===========================================================================
  
  cellchat_object <- CellChat::filterCommunication(
    cellchat_object,
    min.cells = 10
  )
  
  
  # ===========================================================================
  # 6. COMPUTE PATHWAY-LEVEL COMMUNICATION
  # ===========================================================================
  
  cellchat_object <- CellChat::computeCommunProbPathway(
    cellchat_object
  )
  
  
  # ===========================================================================
  # 7. AGGREGATE COMMUNICATION NETWORK
  # ===========================================================================
  
  cellchat_object <- CellChat::aggregateNet(
    cellchat_object
  )
  
  
  # ===========================================================================
  # 8. EXTRACT COMMUNICATION TABLE
  # ===========================================================================
  
  communication_table <- CellChat::subsetCommunication(
    cellchat_object
  )
  
  
  cat("\n============================================================\n")
  cat("COMPLETED:", condition_name, "\n")
  cat("============================================================\n")
  
  
  cat("\nNumber of inferred ligand-receptor interactions:\n")
  
  print(
    nrow(
      communication_table
    )
  )
  
  
  cat("\nNumber of inferred signaling pathways:\n")
  
  print(
    length(
      cellchat_object@netP$pathways
    )
  )
  
  
  return(cellchat_object)
}

cat("\nDoes run_cellchat still contain projectData?\n")

print(
  grepl(
    "projectData",
    paste(
      deparse(body(run_cellchat)),
      collapse = "\n"
    )
  )
)

cat("\nDoes run_cellchat use raw.use = TRUE?\n")

print(
  grepl(
    "raw.use = TRUE",
    paste(
      deparse(body(run_cellchat)),
      collapse = "\n"
    ),
    fixed = TRUE
  )
)
# =============================================================================
###############################################################################
# 182. RUN CELLCHAT — WT
###############################################################################

set.seed(100)

cellchat_WT <- run_cellchat(
  cellchat_WT,
  "24-month WT"
)


# =============================================================================
# 183. RUN CELLCHAT — APP23
# =============================================================================

cellchat_APP23 <- run_cellchat(
  cellchat_APP23,
  "24-month APP23"
)


# =============================================================================
# 184. SAVE INDIVIDUAL OBJECTS
# =============================================================================

saveRDS(
  cellchat_WT,
  "RESULTS/CELLCHAT/cellchat_24month_WT.rds"
)


saveRDS(
  cellchat_APP23,
  "RESULTS/CELLCHAT/cellchat_24month_APP23.rds"
)


# =============================================================================
# 185. EXPORT CONDITION-SPECIFIC COMMUNICATION TABLES
# =============================================================================

communication_WT <- subsetCommunication(
  cellchat_WT
)


communication_APP23 <- subsetCommunication(
  cellchat_APP23
)


write.csv(
  communication_WT,
  "RESULTS/CELLCHAT/communication_24month_WT.csv",
  row.names = FALSE
)


write.csv(
  communication_APP23,
  "RESULTS/CELLCHAT/communication_24month_APP23.csv",
  row.names = FALSE
)


cat("\nWT inferred interactions:\n")
print(nrow(communication_WT))


cat("\nAPP23 inferred interactions:\n")
print(nrow(communication_APP23))


# =============================================================================
# 186. EXAMINE SIGNALING PATHWAYS
# =============================================================================

cat("\n============================================================\n")
cat("WT SIGNALING PATHWAYS\n")
cat("============================================================\n")

print(
  cellchat_WT@netP$pathways
)


cat("\n============================================================\n")
cat("APP23 SIGNALING PATHWAYS\n")
cat("============================================================\n")

print(
  cellchat_APP23@netP$pathways
)


write.csv(
  data.frame(
    pathway = cellchat_WT@netP$pathways
  ),
  "RESULTS/CELLCHAT/pathways_24month_WT.csv",
  row.names = FALSE
)


write.csv(
  data.frame(
    pathway = cellchat_APP23@netP$pathways
  ),
  "RESULTS/CELLCHAT/pathways_24month_APP23.csv",
  row.names = FALSE
)


# =============================================================================
# 187. INDIVIDUAL NETWORK PLOTS
# =============================================================================

groupSize_WT <- as.numeric(
  table(
    cellchat_WT@idents
  )
)


groupSize_APP23 <- as.numeric(
  table(
    cellchat_APP23@idents
  )
)


# WT interaction number

pdf(
  "FIGURES/CELLCHAT/24month_WT_interaction_number.pdf",
  width = 8,
  height = 8
)

netVisual_circle(
  cellchat_WT@net$count,
  vertex.weight = groupSize_WT,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "24-month WT: number of interactions"
)

dev.off()


# APP23 interaction number

pdf(
  "FIGURES/CELLCHAT/24month_APP23_interaction_number.pdf",
  width = 8,
  height = 8
)

netVisual_circle(
  cellchat_APP23@net$count,
  vertex.weight = groupSize_APP23,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "24-month APP23: number of interactions"
)

dev.off()


# WT interaction strength

pdf(
  "FIGURES/CELLCHAT/24month_WT_interaction_strength.pdf",
  width = 8,
  height = 8
)

netVisual_circle(
  cellchat_WT@net$weight,
  vertex.weight = groupSize_WT,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "24-month WT: interaction strength"
)

dev.off()


# APP23 interaction strength

pdf(
  "FIGURES/CELLCHAT/24month_APP23_interaction_strength.pdf",
  width = 8,
  height = 8
)

netVisual_circle(
  cellchat_APP23@net$weight,
  vertex.weight = groupSize_APP23,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "24-month APP23: interaction strength"
)

dev.off()


# =============================================================================
# 188. MERGE CELLCHAT OBJECTS
#
# IMPORTANT:
# Order is WT first, APP23 second.
#
# Therefore differential plots represent:
# APP23 - WT
#
# Red = increased in APP23
# Blue = decreased in APP23
# =============================================================================

object.list <- list(
  WT = cellchat_WT,
  APP23 = cellchat_APP23
)


cellchat_merged <- mergeCellChat(
  object.list,
  add.names = names(
    object.list
  ),
  cell.prefix = TRUE
)


saveRDS(
  cellchat_merged,
  "RESULTS/CELLCHAT/cellchat_24month_WT_vs_APP23_merged.rds"
)


# =============================================================================
# 189. COMPARE TOTAL NUMBER OF INTERACTIONS
# =============================================================================

p_count <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(
    1,
    2
  )
)


ggsave(
  "FIGURES/CELLCHAT/compare_total_interaction_number.pdf",
  p_count,
  width = 5,
  height = 5
)


# =============================================================================
# 190. COMPARE TOTAL INTERACTION STRENGTH
# =============================================================================

p_weight <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(
    1,
    2
  ),
  measure = "weight"
)


ggsave(
  "FIGURES/CELLCHAT/compare_total_interaction_strength.pdf",
  p_weight,
  width = 5,
  height = 5
)


# =============================================================================
# 191. DIFFERENTIAL INTERACTION NUMBER
#
# Red = increased in APP23
# Blue = decreased in APP23
# =============================================================================

pdf(
  "FIGURES/CELLCHAT/differential_interaction_number_APP23_minus_WT.pdf",
  width = 8,
  height = 8
)


netVisual_diffInteraction(
  cellchat_merged,
  comparison = c(
    1,
    2
  ),
  weight.scale = TRUE,
  measure = "count"
)


dev.off()


# =============================================================================
# 192. DIFFERENTIAL INTERACTION STRENGTH
# =============================================================================

pdf(
  "FIGURES/CELLCHAT/differential_interaction_strength_APP23_minus_WT.pdf",
  width = 8,
  height = 8
)


netVisual_diffInteraction(
  cellchat_merged,
  comparison = c(
    1,
    2
  ),
  weight.scale = TRUE,
  measure = "weight"
)


dev.off()


# =============================================================================
# 193. DIFFERENTIAL NETWORK HEATMAP — NUMBER
# =============================================================================

pdf(
  "FIGURES/CELLCHAT/differential_heatmap_interaction_number.pdf",
  width = 8,
  height = 7
)


print(
  netVisual_heatmap(
    cellchat_merged,
    comparison = c(
      1,
      2
    ),
    measure = "count"
  )
)


dev.off()


# =============================================================================
# 194. DIFFERENTIAL NETWORK HEATMAP — STRENGTH
# =============================================================================

pdf(
  "FIGURES/CELLCHAT/differential_heatmap_interaction_strength.pdf",
  width = 8,
  height = 7
)


print(
  netVisual_heatmap(
    cellchat_merged,
    comparison = c(
      1,
      2
    ),
    measure = "weight"
  )
)


dev.off()


# =============================================================================
# 195. EXTRACT PATHWAY INFORMATION
# =============================================================================

pathways_WT <- cellchat_WT@netP$pathways
pathways_APP23 <- cellchat_APP23@netP$pathways


common_pathways <- intersect(
  pathways_WT,
  pathways_APP23
)


WT_only_pathways <- setdiff(
  pathways_WT,
  pathways_APP23
)


APP23_only_pathways <- setdiff(
  pathways_APP23,
  pathways_WT
)


cat("\n\n")
cat("####################################################################\n")
cat("CELLCHAT PATHWAY SUMMARY\n")
cat("####################################################################\n")


cat("\nCommon pathways:\n")
print(common_pathways)


cat("\nWT-only pathways:\n")
print(WT_only_pathways)


cat("\nAPP23-only pathways:\n")
print(APP23_only_pathways)


pathway_summary <- data.frame(
  category = c(
    rep(
      "Common",
      length(common_pathways)
    ),
    rep(
      "WT_only",
      length(WT_only_pathways)
    ),
    rep(
      "APP23_only",
      length(APP23_only_pathways)
    )
  ),
  pathway = c(
    common_pathways,
    WT_only_pathways,
    APP23_only_pathways
  )
)


write.csv(
  pathway_summary,
  "RESULTS/CELLCHAT/pathway_presence_summary.csv",
  row.names = FALSE
)


# =============================================================================
# 196. CALCULATE CENTRALITY FOR EACH CONDITION
#
# Allows later identification of major signaling senders and receivers.
# =============================================================================

cellchat_WT <- netAnalysis_computeCentrality(
  cellchat_WT,
  slot.name = "netP"
)


cellchat_APP23 <- netAnalysis_computeCentrality(
  cellchat_APP23,
  slot.name = "netP"
)


saveRDS(
  cellchat_WT,
  "RESULTS/CELLCHAT/cellchat_24month_WT_with_centrality.rds"
)


saveRDS(
  cellchat_APP23,
  "RESULTS/CELLCHAT/cellchat_24month_APP23_with_centrality.rds"
)


# =============================================================================
# 197. SIGNALING ROLE HEATMAPS
# =============================================================================

pdf(
  "FIGURES/CELLCHAT/WT_outgoing_signaling_roles.pdf",
  width = 10,
  height = 8
)


print(
  netAnalysis_signalingRole_heatmap(
    cellchat_WT,
    pattern = "outgoing"
  )
)


dev.off()


pdf(
  "FIGURES/CELLCHAT/APP23_outgoing_signaling_roles.pdf",
  width = 10,
  height = 8
)


print(
  netAnalysis_signalingRole_heatmap(
    cellchat_APP23,
    pattern = "outgoing"
  )
)


dev.off()


pdf(
  "FIGURES/CELLCHAT/WT_incoming_signaling_roles.pdf",
  width = 10,
  height = 8
)


print(
  netAnalysis_signalingRole_heatmap(
    cellchat_WT,
    pattern = "incoming"
  )
)


dev.off()


pdf(
  "FIGURES/CELLCHAT/APP23_incoming_signaling_roles.pdf",
  width = 10,
  height = 8
)


print(
  netAnalysis_signalingRole_heatmap(
    cellchat_APP23,
    pattern = "incoming"
  )
)


dev.off()


# =============================================================================
# 198. SAVE SESSION INFORMATION
# =============================================================================

sink(
  "RESULTS/CELLCHAT/sessionInfo_CellChat.txt"
)

sessionInfo()

sink()


# =============================================================================
# 199. FINAL SUMMARY
# =============================================================================

cat("\n\n")
cat("####################################################################\n")
cat("CELLCHAT ANALYSIS COMPLETED\n")
cat("####################################################################\n")


cat("\n24-month WT:\n")

cat(
  "Interactions = ",
  nrow(
    communication_WT
  ),
  "\n",
  sep = ""
)


cat(
  "Pathways = ",
  length(
    pathways_WT
  ),
  "\n",
  sep = ""
)


cat("\n24-month APP23:\n")

cat(
  "Interactions = ",
  nrow(
    communication_APP23
  ),
  "\n",
  sep = ""
)


cat(
  "Pathways = ",
  length(
    pathways_APP23
  ),
  "\n",
  sep = ""
)


cat("\nCommon pathways = ")
print(
  length(
    common_pathways
  )
)


cat("\nWT-only pathways = ")
print(
  length(
    WT_only_pathways
  )
)


cat("\nAPP23-only pathways = ")
print(
  length(
    APP23_only_pathways
  )
)


cat("\n####################################################################\n")
cat("STOP HERE\n")
cat("####################################################################\n")

###############################################################################
# 200. PRINT CONDITION-SPECIFIC CELLCHAT PATHWAYS
###############################################################################

cat("\n\n")
cat("####################################################################\n")
cat("CELLCHAT PATHWAY DIFFERENCES\n")
cat("####################################################################\n")


cat("\n============================================================\n")
cat("WT-ONLY PATHWAYS\n")
cat("============================================================\n")

print(WT_only_pathways)


cat("\n============================================================\n")
cat("APP23-ONLY PATHWAYS\n")
cat("============================================================\n")

print(APP23_only_pathways)


cat("\n============================================================\n")
cat("COMMON PATHWAYS\n")
cat("============================================================\n")

print(common_pathways)


###############################################################################
# 201. TOTAL COMMUNICATION MATRICES
###############################################################################

cat("\n\n")
cat("####################################################################\n")
cat("TOTAL INTERACTION NUMBER MATRICES\n")
cat("####################################################################\n")


cat("\nWT interaction counts:\n")
print(cellchat_WT@net$count)


cat("\nAPP23 interaction counts:\n")
print(cellchat_APP23@net$count)


cat("\n\n")
cat("####################################################################\n")
cat("TOTAL INTERACTION STRENGTH MATRICES\n")
cat("####################################################################\n")


cat("\nWT interaction strengths:\n")
print(round(cellchat_WT@net$weight, 4))


cat("\nAPP23 interaction strengths:\n")
print(round(cellchat_APP23@net$weight, 4))


###############################################################################
# 202. APP23 - WT DIFFERENCE MATRICES
###############################################################################

count_diff <- cellchat_APP23@net$count - cellchat_WT@net$count

weight_diff <- cellchat_APP23@net$weight - cellchat_WT@net$weight


cat("\n\n")
cat("####################################################################\n")
cat("APP23 - WT DIFFERENCE IN INTERACTION NUMBER\n")
cat("####################################################################\n")

print(count_diff)


cat("\n\n")
cat("####################################################################\n")
cat("APP23 - WT DIFFERENCE IN INTERACTION STRENGTH\n")
cat("####################################################################\n")

print(round(weight_diff, 4))


###############################################################################
# 203. CONVERT DIFFERENCE MATRICES TO TABLES
###############################################################################

count_diff_table <- as.data.frame(
  as.table(count_diff)
)

colnames(count_diff_table) <- c(
  "source",
  "target",
  "difference_count"
)


weight_diff_table <- as.data.frame(
  as.table(weight_diff)
)

colnames(weight_diff_table) <- c(
  "source",
  "target",
  "difference_weight"
)


###############################################################################
# 204. SORT BY LARGEST ABSOLUTE CHANGES
###############################################################################

count_diff_table <- count_diff_table[
  order(
    abs(count_diff_table$difference_count),
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]


weight_diff_table <- weight_diff_table[
  order(
    abs(weight_diff_table$difference_weight),
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]


cat("\n\n")
cat("####################################################################\n")
cat("TOP CHANGES IN INTERACTION NUMBER\n")
cat("####################################################################\n")

print(
  head(
    count_diff_table,
    16
  ),
  row.names = FALSE
)


cat("\n\n")
cat("####################################################################\n")
cat("TOP CHANGES IN INTERACTION STRENGTH\n")
cat("####################################################################\n")

print(
  head(
    weight_diff_table,
    16
  ),
  row.names = FALSE
)


###############################################################################
# 205. SAVE TABLES
###############################################################################

write.csv(
  count_diff_table,
  "RESULTS/CELLCHAT/APP23_minus_WT_interaction_count_changes.csv",
  row.names = FALSE
)


write.csv(
  weight_diff_table,
  "RESULTS/CELLCHAT/APP23_minus_WT_interaction_strength_changes.csv",
  row.names = FALSE
)


###############################################################################
# 206. STOP
###############################################################################

cat("\n####################################################################\n")
cat("STOP HERE\n")
cat("####################################################################\n")

###############################################################################
# 207. TARGETED PATHWAY ANALYSIS
#
# Focus:
# 1. Cholesterol / Desmosterol
# 2. APP23-only pathways
# 3. WT-only pathways
###############################################################################

pathways_of_interest <- c(
  "Cholesterol",
  "Desmosterol",
  "SEMA6",
  "TENASCIN",
  "PDGF",
  "IGF",
  "CypA",
  "AGRN",
  "ncWNT"
)


###############################################################################
# 208. PATHWAY PRESENCE TABLE
###############################################################################

pathway_presence <- data.frame(
  pathway = pathways_of_interest,
  WT = pathways_of_interest %in% cellchat_WT@netP$pathways,
  APP23 = pathways_of_interest %in% cellchat_APP23@netP$pathways
)


cat("\n\n")
cat("####################################################################\n")
cat("TARGETED PATHWAY PRESENCE\n")
cat("####################################################################\n")

print(
  pathway_presence,
  row.names = FALSE
)


write.csv(
  pathway_presence,
  "RESULTS/CELLCHAT/targeted_pathway_presence.csv",
  row.names = FALSE
)


###############################################################################
# 209. FUNCTION TO EXTRACT PATHWAY COMMUNICATION
###############################################################################

extract_pathway_communication <- function(
    cellchat_object,
    pathway_name,
    condition_name
) {
  
  if (
    !pathway_name %in%
    cellchat_object@netP$pathways
  ) {
    
    return(NULL)
  }
  
  
  tmp <- CellChat::subsetCommunication(
    cellchat_object,
    signaling = pathway_name
  )
  
  
  if (is.null(tmp) || nrow(tmp) == 0) {
    
    return(NULL)
  }
  
  
  tmp$condition <- condition_name
  tmp$pathway_requested <- pathway_name
  
  return(tmp)
}


###############################################################################
# 210. EXTRACT ALL TARGETED PATHWAYS
###############################################################################

targeted_communication <- list()


for (pathway_name in pathways_of_interest) {
  
  WT_tmp <- extract_pathway_communication(
    cellchat_WT,
    pathway_name,
    "WT"
  )
  
  
  APP_tmp <- extract_pathway_communication(
    cellchat_APP23,
    pathway_name,
    "APP23"
  )
  
  
  if (!is.null(WT_tmp)) {
    
    targeted_communication[[
      paste0(
        pathway_name,
        "__WT"
      )
    ]] <- WT_tmp
    
  }
  
  
  if (!is.null(APP_tmp)) {
    
    targeted_communication[[
      paste0(
        pathway_name,
        "__APP23"
      )
    ]] <- APP_tmp
    
  }
}


targeted_communication_combined <- dplyr::bind_rows(
  targeted_communication
)


write.csv(
  targeted_communication_combined,
  "RESULTS/CELLCHAT/targeted_pathway_LR_interactions.csv",
  row.names = FALSE
)


###############################################################################
# 211. SUMMARIZE EACH PATHWAY
###############################################################################

pathway_LR_summary <- targeted_communication_combined %>%
  
  dplyr::group_by(
    pathway_requested,
    condition
  ) %>%
  
  dplyr::summarise(
    
    LR_interactions = dplyr::n(),
    
    total_probability = sum(
      prob,
      na.rm = TRUE
    ),
    
    mean_probability = mean(
      prob,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


cat("\n\n")
cat("####################################################################\n")
cat("TARGETED PATHWAY SUMMARY\n")
cat("####################################################################\n")

print(
  pathway_LR_summary,
  n = Inf
)


write.csv(
  pathway_LR_summary,
  "RESULTS/CELLCHAT/targeted_pathway_summary.csv",
  row.names = FALSE
)


###############################################################################
# 212. SOURCE -> TARGET SUMMARY
###############################################################################

pathway_source_target <- targeted_communication_combined %>%
  
  dplyr::group_by(
    pathway_requested,
    condition,
    source,
    target
  ) %>%
  
  dplyr::summarise(
    
    LR_interactions = dplyr::n(),
    
    total_probability = sum(
      prob,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


pathway_source_target <- pathway_source_target %>%
  
  dplyr::arrange(
    pathway_requested,
    condition,
    dplyr::desc(total_probability)
  )


cat("\n\n")
cat("####################################################################\n")
cat("PATHWAY SOURCE -> TARGET COMMUNICATION\n")
cat("####################################################################\n")

print(
  pathway_source_target,
  n = Inf
)


write.csv(
  pathway_source_target,
  "RESULTS/CELLCHAT/targeted_pathway_source_target.csv",
  row.names = FALSE
)


###############################################################################
# 213. CHOLESTEROL
###############################################################################

cat("\n\n")
cat("####################################################################\n")
cat("CHOLESTEROL SIGNALING\n")
cat("####################################################################\n")


cholesterol_results <- pathway_source_target %>%
  
  dplyr::filter(
    pathway_requested == "Cholesterol"
  )


print(
  cholesterol_results,
  n = Inf
)


###############################################################################
# 214. DESMOSTEROL
###############################################################################

cat("\n\n")
cat("####################################################################\n")
cat("DESMOSTEROL SIGNALING\n")
cat("####################################################################\n")


desmosterol_results <- pathway_source_target %>%
  
  dplyr::filter(
    pathway_requested == "Desmosterol"
  )


print(
  desmosterol_results,
  n = Inf
)


###############################################################################
# 215. APP23-ONLY PATHWAYS
###############################################################################

cat("\n\n")
cat("####################################################################\n")
cat("APP23-ONLY PATHWAYS — SOURCE/TARGET DETAILS\n")
cat("####################################################################\n")


APP23_unique_results <- pathway_source_target %>%
  
  dplyr::filter(
    pathway_requested %in% c(
      "SEMA6",
      "TENASCIN",
      "PDGF",
      "IGF"
    )
  )


print(
  APP23_unique_results,
  n = Inf
)


###############################################################################
# 216. WT-ONLY PATHWAYS
###############################################################################

cat("\n\n")
cat("####################################################################\n")
cat("WT-ONLY PATHWAYS — SOURCE/TARGET DETAILS\n")
cat("####################################################################\n")


WT_unique_results <- pathway_source_target %>%
  
  dplyr::filter(
    pathway_requested %in% c(
      "CypA",
      "AGRN",
      "ncWNT"
    )
  )


print(
  WT_unique_results,
  n = Inf
)


###############################################################################
# 217. SAVE TARGETED RESULTS
###############################################################################

saveRDS(
  targeted_communication,
  "RESULTS/CELLCHAT/targeted_pathway_results.rds"
)


cat("\n\n")
cat("####################################################################\n")
cat("TARGETED CELLCHAT PATHWAY ANALYSIS COMPLETED\n")
cat("####################################################################\n")

cat("\nSTOP HERE\n") 

###############################################################################
# 218. FINAL FIGURE GENERATION
#
# APP23 HIPPOCAMPAL snRNA-seq PROJECT
#
# Main figures:
#   Figure 1A = UMAP by 10 fine clusters
#   Figure 1B = UMAP by broad neuronal subtype
#   Figure 1C = canonical marker DotPlot
#   Figure 2  = pathway-remodeling summary
#   Figure 3  = CellChat APP23-WT interaction-strength heatmap
#   Figure 4  = targeted CellChat cholesterol/desmosterol summary
###############################################################################


###############################################################################
# 219. LOAD REQUIRED PACKAGES
###############################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(100)


###############################################################################
# 220. CREATE FINAL FIGURE DIRECTORY
###############################################################################

dir.create(
  "FIGURES/FINAL",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "RESULTS/FINAL",
  recursive = TRUE,
  showWarnings = FALSE
)


###############################################################################
# 221. CHECK FINAL METADATA
###############################################################################

cat("\n####################################################################\n")
cat("FINAL METADATA CHECK\n")
cat("####################################################################\n")

print(
  table(
    app23$seurat_clusters
  )
)

cat("\nBroad subtypes:\n")

print(
  table(
    app23$broad_subtype
  )
)

cat("\nAge x genotype:\n")

print(
  table(
    app23$age,
    app23$genotype
  )
)


###############################################################################
# 222. SET CLUSTER ORDER
###############################################################################

app23$seurat_clusters <- factor(
  app23$seurat_clusters,
  levels = as.character(0:9)
)


app23$broad_subtype <- factor(
  app23$broad_subtype,
  levels = c(
    "DG",
    "CA1_like",
    "CA3",
    "Inhibitory"
  )
)


###############################################################################
# 223. FIGURE 1A — UMAP OF 10 FINE CLUSTERS
###############################################################################

p_umap_clusters <- DimPlot(
  app23,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.6
) +
  ggtitle(
    "Fine neuronal clusters"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_text(
      face = "bold"
    )
  )


print(
  p_umap_clusters
)


ggsave(
  filename = "FIGURES/FINAL/Figure_1A_UMAP_10_clusters.pdf",
  plot = p_umap_clusters,
  width = 8,
  height = 6
)


ggsave(
  filename = "FIGURES/FINAL/Figure_1A_UMAP_10_clusters.png",
  plot = p_umap_clusters,
  width = 8,
  height = 6,
  dpi = 400
)


###############################################################################
# 224. FIGURE 1B — UMAP OF BROAD NEURONAL SUBTYPES
###############################################################################

p_umap_broad <- DimPlot(
  app23,
  reduction = "umap",
  group.by = "broad_subtype",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.6
) +
  ggtitle(
    "Broad hippocampal neuronal populations"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank()
  )


print(
  p_umap_broad
)


ggsave(
  filename = "FIGURES/FINAL/Figure_1B_UMAP_broad_subtypes.pdf",
  plot = p_umap_broad,
  width = 8,
  height = 6
)


ggsave(
  filename = "FIGURES/FINAL/Figure_1B_UMAP_broad_subtypes.png",
  plot = p_umap_broad,
  width = 8,
  height = 6,
  dpi = 400
)


###############################################################################
# 225. FIGURE 1C — CANONICAL MARKER DOTPLOT
###############################################################################

canonical_markers <- c(
  "Prox1",
  "Mpped1",
  "Mndal",
  "Gad1",
  "Gad2",
  "Slc17a7",
  "Camk2a",
  "Rbfox3",
  "Snap25",
  "Syp"
)


canonical_markers <- canonical_markers[
  canonical_markers %in%
    rownames(app23)
]


p_marker_dotplot <- DotPlot(
  app23,
  features = canonical_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  labs(
    title = "Canonical neuronal markers",
    x = "Marker gene",
    y = "Fine cluster"
  ) +
  theme_classic(base_size = 13)


print(
  p_marker_dotplot
)


ggsave(
  filename = "FIGURES/FINAL/Figure_1C_marker_DotPlot.pdf",
  plot = p_marker_dotplot,
  width = 10,
  height = 6
)


ggsave(
  filename = "FIGURES/FINAL/Figure_1C_marker_DotPlot.png",
  plot = p_marker_dotplot,
  width = 10,
  height = 6,
  dpi = 400
)


###############################################################################
# 226. COMBINE FIGURE 1
###############################################################################

figure_1 <- (
  p_umap_clusters |
    p_umap_broad
) /
  p_marker_dotplot +
  plot_annotation(
    title = "Cellular organization of the APP23 hippocampal neuronal dataset"
  )


print(
  figure_1
)


ggsave(
  filename = "FIGURES/FINAL/Figure_1_combined.pdf",
  plot = figure_1,
  width = 15,
  height = 12
)


ggsave(
  filename = "FIGURES/FINAL/Figure_1_combined.png",
  plot = figure_1,
  width = 15,
  height = 12,
  dpi = 400
)


###############################################################################
# 227. FIGURE 2 — PATHWAY REMODELING SUMMARY
#
# Based on our primary pseudobulk GO-GSEA:
#
# 6 month:
#   DG          = 4 significant pathways
#   CA1_like    = 0
#   CA3         = 0
#   Inhibitory  = 0
#
# 24 month:
#   DG          = 8
#   CA1_like    = 20
#   CA3         = 0
#   Inhibitory  = 0
###############################################################################

pathway_remodeling_final <- data.frame(
  
  age = c(
    "6_month",
    "6_month",
    "6_month",
    "6_month",
    "24_month",
    "24_month",
    "24_month",
    "24_month"
  ),
  
  subtype = c(
    "DG",
    "CA1_like",
    "CA3",
    "Inhibitory",
    "DG",
    "CA1_like",
    "CA3",
    "Inhibitory"
  ),
  
  significant_pathways = c(
    4,
    0,
    0,
    0,
    8,
    20,
    0,
    0
  )
  
)


pathway_remodeling_final$age <- factor(
  pathway_remodeling_final$age,
  levels = c(
    "6_month",
    "24_month"
  ),
  labels = c(
    "6 months",
    "24 months"
  )
)


pathway_remodeling_final$subtype <- factor(
  pathway_remodeling_final$subtype,
  levels = c(
    "DG",
    "CA1_like",
    "CA3",
    "Inhibitory"
  )
)


p_pathway_remodeling <- ggplot(
  pathway_remodeling_final,
  aes(
    x = subtype,
    y = significant_pathways,
    fill = age
  )
) +
  geom_col(
    position = position_dodge(
      width = 0.75
    ),
    width = 0.7
  ) +
  geom_text(
    aes(
      label = significant_pathways
    ),
    position = position_dodge(
      width = 0.75
    ),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    title = "Subtype-specific pathway remodeling in APP23",
    x = "Neuronal population",
    y = "Number of significant GO pathways",
    fill = "Age"
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )


print(
  p_pathway_remodeling
)


ggsave(
  filename = "FIGURES/FINAL/Figure_2_pathway_remodeling.pdf",
  plot = p_pathway_remodeling,
  width = 8,
  height = 6
)


ggsave(
  filename = "FIGURES/FINAL/Figure_2_pathway_remodeling.png",
  plot = p_pathway_remodeling,
  width = 8,
  height = 6,
  dpi = 400
)


write.csv(
  pathway_remodeling_final,
  "RESULTS/FINAL/pathway_remodeling_summary.csv",
  row.names = FALSE
)


###############################################################################
# 228. FIGURE 3 — CELLCHAT APP23 - WT INTERACTION-STRENGTH HEATMAP
###############################################################################

weight_diff <- cellchat_APP23@net$weight -
  cellchat_WT@net$weight


weight_diff_long <- as.data.frame(
  as.table(weight_diff)
)


colnames(
  weight_diff_long
) <- c(
  "Source",
  "Target",
  "Difference"
)


weight_diff_long$Source <- factor(
  weight_diff_long$Source,
  levels = c(
    "DG",
    "CA1_like",
    "CA3",
    "Inhibitory"
  )
)


weight_diff_long$Target <- factor(
  weight_diff_long$Target,
  levels = c(
    "DG",
    "CA1_like",
    "CA3",
    "Inhibitory"
  )
)


p_cellchat_heatmap <- ggplot(
  weight_diff_long,
  aes(
    x = Target,
    y = Source,
    fill = Difference
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Difference
      )
    ),
    size = 4
  ) +
  scale_fill_gradient2(
    midpoint = 0
  ) +
  labs(
    title = "Change in neuronal communication strength",
    subtitle = "APP23 minus WT, 24 months",
    x = "Receiver",
    y = "Sender",
    fill = "Δ strength"
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )


print(
  p_cellchat_heatmap
)


ggsave(
  filename = "FIGURES/FINAL/Figure_3_CellChat_strength_difference.pdf",
  plot = p_cellchat_heatmap,
  width = 8,
  height = 7
)


ggsave(
  filename = "FIGURES/FINAL/Figure_3_CellChat_strength_difference.png",
  plot = p_cellchat_heatmap,
  width = 8,
  height = 7,
  dpi = 400
)


###############################################################################
# 229. CREATE CELLCHAT PATHWAY SUMMARY FOR CHOLESTEROL / DESMOSTEROL
###############################################################################

cellchat_metabolic_summary <- pathway_LR_summary %>%
  
  dplyr::filter(
    pathway_requested %in%
      c(
        "Cholesterol",
        "Desmosterol"
      )
  )


print(
  cellchat_metabolic_summary
)


###############################################################################
# 230. FIGURE 4 — CHOLESTEROL / DESMOSTEROL SIGNALING
###############################################################################

cellchat_metabolic_summary$condition <- factor(
  cellchat_metabolic_summary$condition,
  levels = c(
    "WT",
    "APP23"
  )
)


p_metabolic_cellchat <- ggplot(
  cellchat_metabolic_summary,
  aes(
    x = pathway_requested,
    y = total_probability,
    fill = condition
  )
) +
  geom_col(
    position = position_dodge(
      width = 0.75
    ),
    width = 0.7
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.4f",
        total_probability
      )
    ),
    position = position_dodge(
      width = 0.75
    ),
    vjust = -0.4,
    size = 3.8
  ) +
  labs(
    title = "Sterol-associated neuronal communication",
    subtitle = "24-month WT versus APP23",
    x = NULL,
    y = "Summed CellChat communication probability",
    fill = "Condition"
  ) +
  theme_classic(
    base_size = 14
  )


print(
  p_metabolic_cellchat
)


ggsave(
  filename = "FIGURES/FINAL/Figure_4_Cholesterol_Desmosterol_CellChat.pdf",
  plot = p_metabolic_cellchat,
  width = 7,
  height = 6
)


ggsave(
  filename = "FIGURES/FINAL/Figure_4_Cholesterol_Desmosterol_CellChat.png",
  plot = p_metabolic_cellchat,
  width = 7,
  height = 6,
  dpi = 400
)


###############################################################################
# 231. SAVE CELLCHAT METABOLIC SUMMARY
###############################################################################

write.csv(
  cellchat_metabolic_summary,
  "RESULTS/FINAL/CellChat_Cholesterol_Desmosterol_summary.csv",
  row.names = FALSE
)


###############################################################################
# 232. CREATE FINAL INTEGRATED BIOLOGICAL SUMMARY TABLE
###############################################################################

integrated_summary <- data.frame(
  
  Age = c(
    "6 months",
    "24 months",
    "24 months",
    "24 months"
  ),
  
  Population = c(
    "DG",
    "DG",
    "CA1-like",
    "CA3 / Inhibitory"
  ),
  
  GSEA_result = c(
    
    "Reduced synaptic transmission and synaptic vesicle pathways",
    
    "Reduced synapse organization, projection development, and morphogenesis",
    
    "Strong enrichment of cholesterol and sterol biosynthetic programs",
    
    "No FDR-significant GO pathway remodeling detected"
  ),
  
  TF_result = c(
    
    "Exploratory Neurod1, Jun and Fos decrease; Nrf1 and Trp53 increase",
    
    "Exploratory Rela, Zeb1 and Creb3l1 candidates",
    
    "Srebf2 emerged as an exploratory candidate regulator",
    
    "No significant differential TF activity"
  ),
  
  CellChat_result = c(
    
    "Not evaluated",
    
    "24-month network shows altered DG communication",
    
    "Altered sterol-associated and neuronal communication",
    
    "Strong network rewiring despite limited pseudobulk pathway remodeling"
  )
  
)


print(
  integrated_summary
)


write.csv(
  integrated_summary,
  "RESULTS/FINAL/integrated_biological_summary.csv",
  row.names = FALSE
)


###############################################################################
# 233. SAVE FINAL ANALYSIS OBJECTS
###############################################################################

saveRDS(
  app23,
  "RESULTS/FINAL/app23_final_Seurat_object.rds"
)


saveRDS(
  cellchat_WT,
  "RESULTS/FINAL/cellchat_24m_WT_final.rds"
)


saveRDS(
  cellchat_APP23,
  "RESULTS/FINAL/cellchat_24m_APP23_final.rds"
)


###############################################################################
# 234. SAVE SESSION INFORMATION
###############################################################################

capture.output(
  sessionInfo(),
  file = "RESULTS/FINAL/sessionInfo_final.txt"
)


###############################################################################
# 235. FINAL CHECK
###############################################################################

cat("\n\n")
cat("####################################################################\n")
cat("FINAL FIGURES COMPLETED\n")
cat("####################################################################\n")

cat("\nCreated:\n")
cat("Figure 1A = 10-cluster UMAP\n")
cat("Figure 1B = broad-subtype UMAP\n")
cat("Figure 1C = marker DotPlot\n")
cat("Figure 1  = combined clustering/annotation figure\n")
cat("Figure 2  = pathway-remodeling summary\n")
cat("Figure 3  = CellChat APP23-WT network-strength difference\n")
cat("Figure 4  = cholesterol/desmosterol CellChat comparison\n")

cat("\nFinal analysis objects saved.\n")

cat("\n####################################################################\n")
cat("STOP HERE\n")
cat("####################################################################\n")

###############################################################################
# CHECK ALL FIGURE FILES SAVED IN THE PROJECT
###############################################################################

project_path <- "D:/Master/ScRNA seq/Course project"

all_figures <- list.files(
  project_path,
  pattern = "\\.(png|pdf|jpg|jpeg|tiff)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat("\n============================================================\n")
cat("ALL SAVED FIGURE FILES\n")
cat("============================================================\n\n")

print(all_figures)

cat("\n============================================================\n")
cat("TOTAL FIGURE FILES:", length(all_figures), "\n")
cat("============================================================\n")


# PNG files only
png_files <- all_figures[
  grepl("\\.png$", all_figures, ignore.case = TRUE)
]

cat("\n\n============================================================\n")
cat("PNG FILES\n")
cat("============================================================\n\n")

print(png_files)

cat("\nTOTAL PNG FILES:", length(png_files), "\n")


# PDF files only
pdf_files <- all_figures[
  grepl("\\.pdf$", all_figures, ignore.case = TRUE)
]

cat("\n\n============================================================\n")
cat("PDF FIGURE FILES\n")
cat("============================================================\n\n")

print(pdf_files)

cat("\nTOTAL PDF FILES:", length(pdf_files), "\n")

###############################################################################
# 236. CONVERT ALL EXISTING PDF FIGURES TO HIGH-RESOLUTION PNG
#
# This does NOT rerun any analysis.
# It keeps all original PDF files and creates PNG copies.
###############################################################################

# Install pdftools if necessary
if (!requireNamespace("pdftools", quietly = TRUE)) {
  install.packages("pdftools")
}

library(pdftools)

project_path <- "D:/Master/ScRNA seq/Course project"

figure_root <- file.path(
  project_path,
  "FIGURES"
)


###############################################################################
# 237. FIND ALL PDF FIGURES
###############################################################################

pdf_files <- list.files(
  figure_root,
  pattern = "\\.pdf$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat("\nPDF figures found:", length(pdf_files), "\n")


###############################################################################
# 238. CONVERT MISSING PNG FILES
###############################################################################

converted <- character()
already_exists <- character()
failed <- character()


for (pdf_file in pdf_files) {
  
  png_file <- sub(
    "\\.pdf$",
    ".png",
    pdf_file,
    ignore.case = TRUE
  )
  
  # Do not overwrite PNGs that already exist
  if (file.exists(png_file)) {
    
    already_exists <- c(
      already_exists,
      png_file
    )
    
    next
  }
  
  cat(
    "\nConverting:\n",
    basename(pdf_file),
    "\n"
  )
  
  result <- tryCatch({
    
    pdftools::pdf_convert(
      pdf = pdf_file,
      format = "png",
      dpi = 400,
      filenames = png_file,
      verbose = FALSE
    )
    
    TRUE
    
  }, error = function(e) {
    
    message(
      "FAILED: ",
      basename(pdf_file),
      " | ",
      e$message
    )
    
    FALSE
  })
  
  
  if (result) {
    
    converted <- c(
      converted,
      png_file
    )
    
  } else {
    
    failed <- c(
      failed,
      pdf_file
    )
  }
}


###############################################################################
# 239. REPORT CONVERSION RESULTS
###############################################################################

cat("\n\n")
cat("============================================================\n")
cat("PDF -> PNG CONVERSION COMPLETE\n")
cat("============================================================\n")

cat(
  "\nNew PNG files created:",
  length(converted),
  "\n"
)

cat(
  "PNG files already present:",
  length(already_exists),
  "\n"
)

cat(
  "Failed conversions:",
  length(failed),
  "\n"
)


###############################################################################
# 240. SHOW FAILED FILES, IF ANY
###############################################################################

if (length(failed) > 0) {
  
  cat("\nFiles that failed conversion:\n")
  
  print(failed)
  
} else {
  
  cat("\nAll PDF figures were successfully represented as PNG.\n")
}


###############################################################################
# 241. FINAL PNG INVENTORY
###############################################################################

all_png <- list.files(
  figure_root,
  pattern = "\\.png$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)


cat("\n")
cat("============================================================\n")
cat("FINAL PNG INVENTORY\n")
cat("============================================================\n")

cat(
  "\nTotal PNG files:",
  length(all_png),
  "\n\n"
)

print(all_png)


###############################################################################
# 242. CHECK PDF/PNG PAIRS
###############################################################################

missing_png <- pdf_files[
  !file.exists(
    sub(
      "\\.pdf$",
      ".png",
      pdf_files,
      ignore.case = TRUE
    )
  )
]


cat("\n")
cat("============================================================\n")
cat("FINAL VERIFICATION\n")
cat("============================================================\n")

cat(
  "\nPDF figures without matching PNG:",
  length(missing_png),
  "\n"
)


if (length(missing_png) > 0) {
  
  print(missing_png)
  
} else {
  
  cat(
    "\nSUCCESS: Every PDF figure has a corresponding PNG file.\n"
  )
}


###############################################################################
# 243. STOP
###############################################################################

cat("\n")
cat("####################################################################\n")
cat("ALL FIGURE CONVERSION FINISHED\n")
cat("####################################################################\n")