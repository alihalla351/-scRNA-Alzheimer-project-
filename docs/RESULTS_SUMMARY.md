# Results summary

## Biological question

**Which hippocampal neuronal subtypes exhibit the strongest pathway-level transcriptional remodeling in APP23 mice, and which biological programs distinguish susceptible from relatively preserved neuronal populations?**

## Dataset and neuronal populations

The analysis uses the processed GSE141044 neuronal single-nucleus expression dataset: 3,280 nuclei from 11 mice spanning 6- and 24-month WT and APP23 groups. Primary inference uses four broad neuronal populations: DG, CA1-like, CA3, and inhibitory neurons.

## Pseudobulk differential expression

Differential expression was performed at the biological-sample level with edgeR rather than treating individual nuclei as independent replicates. Under the strict criterion of FDR < 0.05 and |logFC| > 0.25, relatively few individual genes were significant in each age-by-subtype comparison. Therefore, the biological interpretation emphasizes ranked pathway analysis rather than claiming extensive remodeling from DEG counts alone.

Exact uppercase `APP` was strongly increased across comparisons and is interpreted as transgene-associated. `Thy1` was also strongly associated with APP23 status. Primary pathway/TF sensitivity analyses exclude exact uppercase `APP` and `Thy1`, while endogenous title-case `App` is retained.

## Ranked GO Biological Process GSEA

### 6-month DG

Four FDR-significant pathways were depleted in APP23 DG neurons, centered on synaptic function:

- regulation of trans-synaptic signaling
- modulation of chemical synaptic transmission
- synaptic vesicle cycle
- vesicle-mediated transport in synapse

This supports an early DG phenotype involving reduced synaptic transmission/vesicle-associated programs.

### 24-month DG

The later DG phenotype is dominated by depleted pathways related to neuronal structure and synaptic organization, including regulation of neuron projection development, postsynapse organization, regulation of synapse organization, regulation of cell morphogenesis, and regulation of synapse assembly.

Together with the 6-month result, this is consistent with a progression from early synaptic functional disruption toward later synaptic/structural remodeling.

### 24-month CA1-like neurons

CA1-like neurons show the clearest late metabolic pathway phenotype. Strong positive enrichment is observed for cholesterol/sterol biosynthetic programs, including cholesterol biosynthesis via lathosterol/desmosterol and zymosterol metabolism/biosynthesis.

GO terms are redundant and are therefore interpreted collectively as a **sterol/cholesterol biosynthetic program**, rather than as independent biological discoveries.

### CA3 and inhibitory neurons

No GO-BP pathways reached the primary FDR threshold in the analyzed 6- or 24-month CA3 and inhibitory comparisons. These populations are therefore described as relatively preserved under the FDR-defined pathway-remodeling criterion, not as biologically unaffected.

## Transcription-factor activity

DoRothEA/decoupleR activity inference followed by sample-level limma testing identified no TFs significant at FDR < 0.05 in the age-by-subtype comparisons. TF results are therefore exploratory.

In 24-month CA1-like neurons, Srebf2 showed a positive APP23-associated activity estimate consistent with the independent sterol/cholesterol GSEA signal, but it was not significant after multiple-testing correction. It is treated as a candidate regulator rather than a confirmed mechanism.

## CellChat

CellChat was used descriptively to compare pooled 24-month WT and APP23 neuronal populations. Because the dataset contains neuronal nuclei only, the analysis can infer only neuron-to-neuron signaling represented by the measured transcriptome; it does not reconstruct the full hippocampal signaling environment.

The comparison indicates **communication rewiring rather than a uniform global increase**. Many source-target interaction strengths decrease in APP23, whereas selected interactions increase, including DG-to-CA1-like signaling. Targeted analysis also identifies changes in cholesterol/desmosterol-associated inferred signaling.

These CellChat findings are descriptive condition-level observations and should not be interpreted as mouse-level statistical evidence or as proof that the intracellular CA1-like sterol program causes the signaling changes.

## Integrated interpretation

APP23-associated hippocampal remodeling is subtype- and age-dependent. DG neurons exhibit early disruption of synaptic transmission/vesicle programs and later structural/synaptic remodeling. CA1-like neurons develop a pronounced cholesterol/sterol-associated transcriptional phenotype at 24 months. CA3 and inhibitory populations show comparatively limited FDR-significant pathway changes under the primary analysis. TF inference nominates Srebf2 only as an exploratory candidate, while CellChat suggests broader reorganization of inferred neuronal communication.

## Key limitations

- The starting dataset contains neuronal nuclei only; non-neuronal hippocampal cell types are absent.
- The processed Smart-seq2-style dataset was already QC-filtered by the original study; no additional computational doublet removal was performed.
- Mitochondrial features were absent from the processed expression matrix, preventing mitochondrial-percentage filtering.
- Sample sizes are small, particularly at 24 months (2 WT versus 3 APP23 mice).
- Strict pseudobulk DE yielded few FDR-significant genes, motivating ranked pathway analysis.
- TF activity results do not survive FDR correction and remain exploratory.
- CellChat pools nuclei by condition and is descriptive rather than replicate-level inference.
- Potential sex-associated signals should not be interpreted without verified sample sex metadata.

## Overall conclusion

The strongest supported project-level conclusion is that APP23-associated neuronal transcriptional remodeling is heterogeneous across hippocampal neuronal populations and changes with age, with early/later synaptic remodeling in DG and a prominent late sterol/cholesterol program in CA1-like neurons.