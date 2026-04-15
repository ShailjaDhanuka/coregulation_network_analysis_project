library(Seurat)
library(hdf5r)
library(dplyr)
library(zellkonverter)
library(WGCNA)

setwd("Desktop/Shailja_everything/CMU_courses/FundamentalsOfBioinformatics/Bioinfo_Project/")

input_h5   <- "GSE260657_filtered.h5ad"
plot_dir   <- "plots"
output_dir <- "./"

dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load and convert ──────────────────────────────────────────────────────────
sce        <- readH5AD(input_h5)
seurat_obj <- as.Seurat(sce, counts = "X", data = NULL)
seurat_obj[["RNA"]] <- seurat_obj[["originalexp"]]
DefaultAssay(seurat_obj) <- "RNA"
obj <- seurat_obj

# ── Preprocess ────────────────────────────────────────────────────────────────
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30)
obj <- RunTSNE(obj, dims = 1:30)

# ── UMAP plot ─────────────────────────────────────────────────────────────────
png(file.path(plot_dir, "umap_plot.png"), width = 2000, height = 1500, res = 300)
print(DimPlot(obj, label = TRUE))
dev.off()

# ── Add symptomatic label ─────────────────────────────────────────────────────
symptom_map <- data.frame(
  human       = factor(1:15),
  symptomatic = c(0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1)
  #               1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
)

obj@meta.data$symptomatic <- symptom_map$symptomatic[
  match(obj@meta.data$human, symptom_map$human)
]

# Verify
table(obj@meta.data$symptomatic)
head(obj@meta.data[, c("human", "symptomatic")])

# ── Save Seurat object ────────────────────────────────────────────────────────
saveRDS(obj, file = file.path(output_dir, "seu_object_preprocessed.rds"))

# ── Pseudobulk by donor ───────────────────────────────────────────────────────
pb_by_donor <- AggregateExpression(
  obj,
  group.by      = "human",
  assay         = "RNA",
  slot          = "counts",
  return.seurat = FALSE
)$RNA

pb_matrix <- t(as.matrix(pb_by_donor))

# ── Build metadata ────────────────────────────────────────────────────────────
pb_meta <- obj@meta.data %>%
  group_by(human) %>%
  summarise(
    orig.ident   = first(orig.ident),
    symptomatic  = first(symptomatic),
    n_cells      = n(),
    mean_n_genes = mean(n_genes),
    mean_counts  = mean(total_counts)
  ) %>%
  arrange(as.numeric(as.character(human))) %>%
  as.data.frame()

rownames(pb_meta) <- pb_meta$human

# ── Normalize ─────────────────────────────────────────────────────────────────
pb_norm <- log1p(sweep(pb_matrix, 1, rowSums(pb_matrix), "/") * 1e6)

# ── Zero variance filter ──────────────────────────────────────────────────────
gene_vars <- apply(pb_norm, 2, var)
pb_filt   <- pb_norm[, gene_vars > 0]
cat("Genes after zero-variance filter:", ncol(pb_filt), "\n")

# ── Mean expression filter ────────────────────────────────────────────────────
gene_means <- colMeans(pb_filt)
pb_filt    <- pb_filt[, gene_means > 0.1]
cat("Genes after mean expression filter:", ncol(pb_filt), "\n")

# ── CV-based HVG selection ────────────────────────────────────────────────────
gene_means_filt <- colMeans(pb_filt)
gene_sds_filt   <- apply(pb_filt, 2, sd)
gene_cv         <- gene_sds_filt / gene_means_filt

top_n     <- 3000
top_genes <- names(sort(gene_cv, decreasing = TRUE))[1:min(top_n, length(gene_cv))]
pb_final  <- pb_filt[, top_genes]
cat("Genes after CV-based HVG selection:", ncol(pb_final), "\n")

# ── Outlier sample detection ──────────────────────────────────────────────────
sample_tree <- hclust(dist(pb_final), method = "average")
png(file.path(plot_dir, "sample_clustering_pretom.png"), width = 2000, height = 1200, res = 300)
plot(sample_tree, main = "Pseudobulk sample clustering", sub = "", xlab = "")
dev.off()

# ── WGCNA QC check ────────────────────────────────────────────────────────────
gsg <- goodSamplesGenes(pb_final, verbose = 3)
if (!gsg$allOK) {
  cat("Removing", sum(!gsg$goodGenes), "bad genes and",
      sum(!gsg$goodSamples), "bad samples\n")
  pb_final <- pb_final[gsg$goodSamples, gsg$goodGenes]
} else {
  cat("All genes and samples passed QC\n")
}
cat("Final matrix dimensions:", dim(pb_final), "\n")  # should be 15 x ~3000

# ── Attach metadata and save ──────────────────────────────────────────────────
attr(pb_final, "metadata") <- pb_meta
saveRDS(pb_final, file = file.path(output_dir, "wgcna_input_matrix.rds"))

# ── Verify ────────────────────────────────────────────────────────────────────
test <- readRDS("wgcna_input_matrix.rds")
cat("Metadata present:", !is.null(attr(test, "metadata")), "\n")  # must be TRUE
cat("Metadata rows:",    nrow(attr(test, "metadata")), "\n")       # must be 15
cat("Symptomatic counts:\n")
print(table(attr(test, "metadata")$symptomatic))                   # must be 7 zeros, 8 ones

