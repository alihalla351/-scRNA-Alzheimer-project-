# Analysis workflow

## 1. Import and metadata

Import the processed GSE141044 dense expression table, construct a Seurat object, and parse each nucleus identifier into biological sample, age, and genotype metadata. Set `orig.ident` to the biological `sample_id`.

## 2. QC verification

The authors' processed dataset contains 3,280 neuronal nuclei. Apply the study-compatible thresholds `nFeature_RNA > 500` and `nCount_RNA > 4000`; all 3,280 nuclei pass. No mitochondrial filtering is applied because mitochondrial features are absent from the processed matrix. Ribosomal percentage can be retained as a diagnostic metric. No additional computational doublet removal is part of the primary workflow.

## 3. Seurat preprocessing

- `NormalizeData(method = "LogNormalize", scale.factor = 10000)`
- `FindVariableFeatures(selection.method = "vst", nfeatures = 2000)`
- Scale variable features
- PCA with 30 PCs
- Use PCs 1–10 for neighborhood graph and UMAP
- No Harmony/integration in the primary workflow
- Graph-based clustering at resolution 0.5
- Fixed random seed 100
- UMAP only; t-SNE is not used

The selected solution contains 10 fine neuronal clusters.

## 4. Annotation

Canonical markers include `Prox1` (DG), `Mpped1` (CA1-associated), `Mndal` (CA3), and `Gad1/Gad2` (inhibitory neurons), supported by broader neuronal/excitatory markers such as `Slc17a7`, `Camk2a`, `Rbfox3`, `Snap25`, and `Syp`.

Primary broad mapping:

- cluster 0 → DG
- clusters 1, 2, 3, 5, 8 → CA1_like
- cluster 4 → CA3
- clusters 6, 7, 9 → Inhibitory

Broad totals: DG 1,119; CA1_like 1,493; CA3 328; Inhibitory 340.

## 5. Mouse-level pseudobulk differential expression

Aggregate raw counts by `sample_id × broad_subtype`. Within each age and broad subtype, use edgeR with `filterByExpr`, TMM normalization, robust dispersion estimation, `glmQLFit`, and `glmQLFTest` for APP23 versus WT. Positive logFC denotes higher expression in APP23.

A strict gene-level summary uses FDR < 0.05 and |logFC| > 0.25. Because strict DEG counts are sparse, biological interpretation emphasizes ranked pathway analysis rather than DEG counts alone.

## 6. Ranked GO Biological Process GSEA

For the primary pathway ranking, exclude exact uppercase `APP` and `Thy1`; retain endogenous title-case `App`. Rank genes using `sign(logFC) * sqrt(F)`, map mouse symbols to Entrez identifiers, retain the strongest absolute rank for duplicated mappings, and run `clusterProfiler::gseGO` with GO Biological Process, BH correction, minimum gene-set size 10 and maximum 500.

Positive NES indicates APP23-upregulated enrichment; negative NES indicates depletion in APP23 relative to WT.

Key FDR-significant programs:

- 6m DG: depletion of synaptic transmission/vesicle programs.
- 24m DG: depletion of synapse organization, assembly, neuronal projection, and morphogenesis programs.
- 24m CA1_like: strong positive cholesterol/sterol biosynthetic enrichment.
- CA3 and Inhibitory: no FDR-significant GO-BP pathways under the primary criterion.

## 7. TF activity inference

Use mouse DoRothEA A/B/C regulons with decoupleR ULM on pseudobulk logCPM expression. Test APP23 versus WT within age × broad subtype with limma. Primary input excludes exact uppercase `APP` and `Thy1` and retains endogenous `App`.

No TF passes FDR < 0.05 across the primary comparisons. TF findings are therefore exploratory. `Srebf2` is an exploratory candidate consistent with the 24m CA1-like sterol program but is not significant after multiple-testing correction.

## 8. CellChat

Primary descriptive CellChat comparison: 24-month WT versus APP23. Use the four broad neuronal populations, normalized RNA, and remove exact uppercase `APP` and `Thy1` while retaining endogenous `App`. Run overexpressed-gene/interactions identification, communication probability with `type = "triMean"`, `raw.use = TRUE`, `population.size = FALSE`, pathway probability, and network aggregation.

Interpret CellChat as descriptive condition-level inference rather than mouse-level statistical evidence. The dataset contains neuronal nuclei only, so signaling represents inferred neuron-to-neuron communication and excludes glia and other hippocampal cell types.

## 9. Integrated interpretation

The primary interpretation is subtype- and age-dependent remodeling: DG shows early synaptic dysfunction followed by later synaptic/structural remodeling, whereas CA1-like neurons show pronounced late cholesterol/sterol-associated remodeling. CA3 and inhibitory neurons are comparatively preserved by the FDR-defined pathway criterion. CellChat suggests communication rewiring rather than a uniform increase in signaling.
