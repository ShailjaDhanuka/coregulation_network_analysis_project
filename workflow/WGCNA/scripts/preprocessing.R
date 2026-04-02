

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

pseudobulk_seu <- AggregateExpression(
  obj,
  group.by = "seurat_clusters",
  return.seurat = TRUE  # <-- keep as Seurat
)

saveRDS(pseudobulk_seu,file=file.path(output_dir, "pseudobulked_seu_obj.rds"))

###########


