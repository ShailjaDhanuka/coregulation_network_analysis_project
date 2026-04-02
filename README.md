# coregulation_network_analysis_project

Multi-ome Data Reference: https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-x-1-standard-2-0-0

Workflow Sections
1. WGCNA
2. Cytoscape
3. ARACNe

Repo structure: all data files are centralized in `data/` (for example `data/raw`, `data/WGCNA`, `data/Cytoscape`), while workflow code and plots are organized under `workflow/WGCNA`, `workflow/Cytoscape`, and `workflow/ARACNe`.

## Project Description

This project applies a multi-tool computational pipeline to reconstruct gene co-expression networks and infer transcription factor (TF) regulatory programs from single-cell RNA-seq (scRNA-seq) data.

Because WGCNA was originally designed for bulk transcriptomics, scRNA-seq data presents limitations for WGCNA due to the inherent sparsity of single-cell data (MDPI). To address this, the data is first pseudobulked by aggregating cell-level counts into sample-level profiles so it becomes compatible with bulk-oriented methods.

WGCNA then builds a weighted co-expression network and groups genes into modules based on shared expression patterns. Module activity can be measured by average gene expression level, and modules may contain key regulators, often those with high network centrality (Nature).

ARACNe is subsequently applied to infer direct TF-to-target relationships within those modules using mutual information, while pruning indirect interactions via the data processing inequality. Finally, Cytoscape is used to visualize the resulting networks and identify hub genes and high-centrality TF nodes as candidate master regulators.

## Gaps In The Field This Addresses

1. WGCNA was not built for single-cell data.
WGCNA reconstructs undirected gene co-expression networks from bulk transcriptomics data, and the resulting network can contain many false-positive associations and reduced interpretability (Oxford Academic) when applied naively to sparse scRNA-seq data. Pseudobulking is a pragmatic solution, but it can collapse cellular heterogeneity.

2. TF inference from expression alone is fundamentally noisy.
Network inference, especially at the individual interaction level, is prone to false positives and false negatives. Predicted regulatory interactions (Oxford Academic) from transcriptomics alone remain difficult to validate. GRN inference from transcriptome data alone can produce many false positives because it does not directly account for mechanisms like TF binding (Oxford Academic).

3. ARACNe has sample-size requirements that scRNA-seq often struggles to meet.
ARACNe requires N >= 100 tissue-specific gene expression profiles representing statistically independent samples (Nature), which is difficult to achieve from pseudobulked scRNA-seq without many biological replicates.

4. GRN benchmarking is unreliable.
Regulatory interactions in databases are often aggregated from broad datasets and may lack specificity to particular biological systems, making evaluation of GRN inference performance less reliable (arXiv).

5. The field is moving toward multimodal integration, which this pipeline does not yet use.
Many newer methods, including LINGER and SCENIC+, combine expression with chromatin accessibility (scATAC-seq) to anchor TF inference to likely binding sites. By contrast, an expression-only WGCNA + ARACNe pipeline does not directly capture this layer of regulatory evidence (Nature).

## Why This Pipeline Is Still Valuable

Despite these gaps, this approach is interpretable, computationally accessible, and well-established.

TF-target interaction inference generates regulons (gene sets regulated by each TF), and comparison of regulon activity between states can suggest key regulators (Nature). Combining module-level co-expression (WGCNA), regulon-level TF inference (ARACNe), and network topology analysis (Cytoscape) provides three orthogonal lines of evidence converging on candidate master regulators, which is often more compelling than any single method alone.

This also serves as a strong foundation before adding multimodal data layers.

## Contributors

## Shailja Dhanuka

-  **Data acquisition**, sourcing and organizing relevant biological datasets for analysis  
- Performed **data preprocessing**, including cleaning, normalization, and quality control  
-  **WGCNA (Weighted Gene Co-expression Network Analysis)** to identify gene modules and key clusters  
-  **Module gene identification** to extract biologically significant gene groups  
- Designed the **end-to-end workflow pipeline**, selecting appropriate tools and methodologies for each stage  
- Managed **tool selectionn** to ensure efficiency and reproducibility  
- Handled **GitHub repository management**, including structure, version control, and documentation organization  

## Ushta Samal
- Performed Cytoscape and identified hub genes
- Designed the **end-to-end workflow pipeline**, selecting appropriate tools and methodologies for each stage  

## Cindy 
- Performed ARACNe to infer direct TF-to-target relationships from the gene modules found from WGCNA
## Mahita 

