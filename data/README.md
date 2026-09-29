# Data acquisition

The primary dataset is GEO accession **GSE141044**. Large source files are not committed to GitHub.

## GEO supplementary files used

- `GSE141044_Genes.csv.gz`
- `GSE141044_barcode.csv.gz`
- `GSE141044_barcode_sequence.txt.gz`
- `GSE141044_matrix.mtx.gz`

The processed expression file named `GSE141044_matrix.mtx.gz` is not a conventional Matrix Market sparse matrix. After decompression it is a dense whitespace-delimited table in which the first row contains nucleus identifiers, the first column of subsequent rows contains gene symbols, and the remaining entries contain expression counts. The analysis import code must therefore parse this layout explicitly rather than using `Matrix::readMM()`.

## Biological samples

| GEO sample | Sample label |
|---|---|
| GSM4194799 | 6-month_WT_1 |
| GSM4194800 | 6-month_WT_2 |
| GSM4194801 | 6-month_WT_3 |
| GSM4194802 | 6-month_APP23_1 |
| GSM4194803 | 6-month_APP23_2 |
| GSM4194804 | 6-month_APP23_3 |
| GSM4194805 | 24-month_WT_1 |
| GSM4194806 | 24-month_WT_2 |
| GSM4194807 | 24-month_APP23_1 |
| GSM4194808 | 24-month_APP23_2 |
| GSM4194809 | 24-month_APP23_3 |

The original study began with 12 mice; one 24-month WT sample was excluded for low quality, leaving 11 biological samples in the processed dataset.

## Processed matrix used here

Initial dimensions after import were 28,693 features × 3,280 nuclei. A feature named exactly `a` contained only eight counts across two nuclei and was removed as an import/artifactual feature, leaving **28,692 genes × 3,280 nuclei** for analysis.

Do not commit downloaded GEO files to this repository. Place them locally under `data/raw/` or another local path ignored by Git.
