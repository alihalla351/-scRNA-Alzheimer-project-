# Methods summary

## Dataset and experimental design

The analysis uses the processed neuronal single-nucleus expression dataset deposited under GEO accession GSE141044. The dataset comprises 3,280 hippocampal neuronal nuclei from 11 mice distributed across 6-month WT, 6-month APP23, 24-month WT, and 24-month APP23 groups. The source study used a modified Smart-seq2 workflow on neuronal nuclei.

## Quality control and preprocessing

The processed matrix was imported into Seurat. A low-information feature named exactly `a` (eight total counts across two nuclei) was removed. Study-compatible QC thresholds of >500 detected genes and >4,000 transcripts per nucleus were applied; all 3,280 nuclei passed. Mitochondrial-percentage filtering was not performed because mitochondrial genes were absent from the processed matrix. Expression was log-normalized with a scale factor of 10,000, 2,000 highly variable features were selected using the variance-stabilizing transformation method, and variable features were scaled. PCA was performed and PCs 1–10 were used for graph construction, clustering, and UMAP. No batch integration was applied in the primary workflow. Graph-based clustering at resolution 0.5 yielded 10 fine neuronal clusters. Random seed 100 was used for reproducible dimensional reduction where applicable.

## Neuronal annotation

Fine clusters were characterized with Seurat marker analysis and canonical hippocampal neuronal markers. The primary biological analysis collapsed fine clusters into DG, CA1-like, CA3, and inhibitory populations. The broad mapping was cluster 0 to DG; clusters 1, 2, 3, 5, and 8 to CA1-like; cluster 4 to CA3; and clusters 6, 7, and 9 to inhibitory neurons.

## Pseudobulk differential expression

Raw counts were aggregated for each biological mouse within each broad neuronal subtype. Differential expression between APP23 and WT was tested separately within each age and subtype using edgeR quasi-likelihood models. Lowly expressed genes were filtered with `filterByExpr`, libraries were TMM normalized, dispersion was estimated robustly, and quasi-likelihood models were fitted with `glmQLFit` followed by `glmQLFTest`. Positive log fold change indicates higher expression in APP23. Strict gene-level significance was defined as FDR < 0.05 and absolute log2 fold change > 0.25.

## Gene-set enrichment analysis

Because strict pseudobulk DEG counts were sparse, pathway remodeling was assessed using ranked GO Biological Process GSEA. The primary ranking excluded exact uppercase `APP` and `Thy1` as a sensitivity strategy for APP23 construct-associated signals while retaining endogenous title-case `App`. Genes were ranked as sign(logFC) × sqrt(F-statistic). Mouse gene symbols were mapped to Entrez identifiers, duplicated mappings were reduced to the strongest absolute rank, and `clusterProfiler::gseGO` was run using BH correction, minimum gene-set size 10, and maximum gene-set size 500.

## Transcription-factor activity

TF activities were inferred from pseudobulk logCPM expression using mouse DoRothEA A/B/C regulons and the decoupleR univariate linear model method. APP23-versus-WT activity differences were tested within each age and broad neuronal subtype using limma. TF results were interpreted as exploratory because no tested TF remained significant at FDR < 0.05.

## Cell-cell communication

CellChat was applied descriptively to 24-month WT and APP23 neuronal nuclei, the age comparison showing the strongest late pathway phenotype. Analyses used the four broad neuronal populations and normalized RNA expression. Exact uppercase `APP` and `Thy1` were removed for the primary sensitivity analysis while endogenous `App` was retained. Communication probabilities were estimated with the trimean approach without population-size weighting, followed by pathway-level probability calculation and network aggregation. Because nuclei were pooled within condition and the dataset contains neurons only, CellChat results are treated as descriptive inferred neuron-to-neuron signaling rather than biological-replicate statistical inference or a complete hippocampal signaling map.
