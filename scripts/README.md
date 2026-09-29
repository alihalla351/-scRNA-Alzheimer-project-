# Analysis scripts

## `GSE141044_APP23_full_project.R`

This is the exact cumulative R project script used during the APP23 hippocampal single-nucleus RNA-seq analysis. It is retained as the primary executable/provenance record rather than reconstructing code from narrative summaries.

The script includes:

- corrected import of the non-standard GSE141044 dense expression matrix
- metadata construction for the 11 biological samples
- QC verification using the study-compatible gene/transcript thresholds
- Seurat normalization and 2,000 highly variable genes
- PCA and PC selection
- graph-based clustering with the selected resolution 0.5
- UMAP visualization only; no t-SNE is used
- marker analysis and neuronal subtype annotation
- mouse/sample-level pseudobulk differential expression with edgeR
- ranked GO Biological Process GSEA with clusterProfiler
- DoRothEA/decoupleR TF activity inference and limma testing
- descriptive 24-month CellChat WT-versus-APP23 analysis
- targeted cholesterol/desmosterol CellChat summaries
- final figure generation and PDF-to-PNG figure conversion

### Important provenance note

The file preserves the cumulative analysis history and therefore contains earlier exploratory/troubleshooting blocks as well as the corrected final workflow. When an earlier block conflicts with a later explicitly corrected block, the later corrected workflow represents the analysis used for final interpretation.

Large raw GEO files and serialized R objects are intentionally excluded from version control.