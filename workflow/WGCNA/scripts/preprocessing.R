

library(Seurat)
library(hdf5r)
# library(SingleR)
# library(celldex)
library(dplyr)
# library(Azimuth)


input_h5 <- "data/raw/10k_PBMC_Multiome_nextgem_Chromium_X_filtered_feature_bc_matrix.h5"
plot_dir <- "workflow/WGCNA/plots"
output_dir <- "data/WGCNA"

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- Read10X_h5(input_h5)
names(counts)

rna_counts <- counts$`Gene Expression`

obj <- CreateSeuratObject(
  counts = rna_counts,
  project = "PBMC_multiome"
)

#removing cells with more than 10% mitochondria percentage
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

obj_no_mt <- subset(
  obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)



# obj@meta.data$percent.mt

obj <- NormalizeData(obj_no_mt)

normalized_counts<- GetAssayData(obj,assay="RNA",layer="data")

obj <- FindVariableFeatures(obj)

obj <- ScaleData(obj)

obj <- RunPCA(obj)

obj <- FindNeighbors(obj)

obj <- FindClusters(obj, resolution = 0.5)

obj <- RunUMAP(obj, dims = 1:30)

# At this point we have:
# 
# cells × genes matrix
# clusters
# UMAP

# some plots for analysis - data visualization

p<-DimPlot(obj, label = TRUE)
png(file.path(plot_dir, "umap_plot.png"), width = 2000, height = 1500, res = 300)  # high-res
print(p)   # you must explicitly print the ggplot object
dev.off()


RunTSNE(obj)

##########################
# save object
#####################333
saveRDS(obj,file = file.path(output_dir, "seu_object_preprocessed.rds"))


# this aggregates the counts not the normalizef counts
pseudobulk_seu <- AggregateExpression(
  obj,
  group.by = "seurat_clusters",
  return.seurat = TRUE  # <-- keep as Seurat
)

saveRDS(pseudobulk_seu,file=file.path(output_dir, "pseudobulked_seu_obj.rds"))

###########




######################################33
# Preprocess fro WGCNA
##################################

setwd("./coregulation_network_analysis_project")
# input_raw <- "data/WGCNA/seu_object_preprocessed.rds"
input_pb <-"data/WGCNA/pseudobulked_seu_obj.rds"
plot_dir <- "workflow/WGCNA/plots"
output_dir <- "data/WGCNA"

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# obj_raw<-readRDS(input_raw)
obj_pb<-readRDS(input_pb)

pb_counts <- GetAssayData(obj_pb, assay = "RNA", layer = "data")
pb_matrix <- as.matrix(pb_counts)

pb_matrix_t<- t(pb_matrix)  # WGCNA: rows = samples, cols = genes

dim(pb_matrix_t)  # n_clusters x n_genes — confirm orientation   (15 36601)


# Remove genes with zero or near-zero variance ─────────────────────
# Genes with no variation across pseudobulk samples carry no coexpression signal
gene_vars <- apply(pb_matrix_t, 2, var)
pb <- pb_matrix_t[, gene_vars > 0]
cat("Genes after zero-variance filter:", ncol(pb), "\n")   # 28745


# Remove genes not meaningfully expressed across pseudobulk samples.
# Using log-normalized data: threshold of 0.1 means expressed at low but
# detectable levels. Adjust based on your data distribution.
gene_means <- colMeans(pb)
hist(gene_means, breaks = 100, main = "Gene mean expression distribution")
abline(v = 0.1, col = "red")  # visualize your cutoff

pb_mean_rm <- pb[, gene_means > 0.1]
cat("Genes after mean expression filter:", ncol(pb_mean_rm), "\n")



# select highly variable genes
# WGCNA on all remaining genes is noisy and slow.
# Top 5000 by coefficient of variation (CV = sd/mean) is standard for bulk-like data.
# CV is preferred over raw variance because it is mean-independent.
gene_means_filt <- colMeans(pb_mean_rm)
gene_sds_filt   <- apply(pb_mean_rm, 2, sd)
gene_cv         <- gene_sds_filt / gene_means_filt

hist(gene_cv, breaks = 100, main = "Coefficient of Variation distribution")
dim(pb_mean_rm)



# Select top 5000 most variable genes by CV
# For smaller pseudobulk matrices (few clusters), 2000-3000 may be more appropriate
top_n <- 5000
top_genes <- names(sort(gene_cv, decreasing = TRUE))[1:min(top_n, length(gene_cv))]
pb_5000 <- pb_mean_rm[, top_genes]
cat("Genes after CV-based HVG selection:", ncol(pb_5000), "\n")



# Outlier sample detection ─────────────────────────────────────────
# Check for pseudobulk samples (clusters) that are transcriptionally aberrant.
# These can distort module detection.
sample_tree <- hclust(dist(pb_5000), method = "average")

png(file.path(plot_dir, "sample_clustering_pretom.png"), 
    width = 2000, height = 1200, res = 300)
plot(sample_tree,
     main = "Pseudobulk sample clustering (check for outliers)",
     sub  = "",
     xlab = "")
dev.off()


# WGCNA goodSamplesGenes check ─────────────────────────────────────
# This checks for genes with too many missing values and samples with
# too many missing genes. Should pass cleanly after the above filters.
gsg <- goodSamplesGenes(pb_5000, verbose = 3)

if (!gsg$allOK) {
  cat("Removing", sum(!gsg$goodGenes), "bad genes and",
      sum(!gsg$goodSamples), "bad samples\n")
  pb_5000 <- pb_5000[gsg$goodSamples, gsg$goodGenes]  # overwrite in place
} else {
  cat("All genes and samples passed QC\n")
}

# Now always print the correct final matrix
cat("Final matrix dimensions:", dim(pb_5000), "\n")

# rows = pseudobulk samples (clusters), cols = genes

# save the filtered file for WGCNA
saveRDS(pb_5000, file = file.path(output_dir, "wgcna_input_matrix.rds"))  # final matrix has 15 clusters and 5000 genes






