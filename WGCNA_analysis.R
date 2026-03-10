




p<-DimPlot(obj, label = TRUE)

print(p)



#we still dont have any metadata, so, 
# either manual annotations or using azimuth (automatic annotations). We will do the automatic annotations using azimuth

# manual annotations for practice
# markers <- FindAllMarkers(obj,
#                           only.pos = TRUE,
#                           min.pct = 0.25,
#                           logfc.threshold = 0.25)
# markers %>%
#   group_by(cluster) %>%
#   top_n(10, avg_log2FC)

# after this, manually identify markers and then label clusters into celtype_sum(

# automatic annotation(cell type identification using azimuth)

# obj <- RunAzimuth(obj, reference = "pbmc")

# to do in preprocessing - remove some cells violin plot, umap and other visualization things

pseudobulk <- AggregateExpression(
  obj,
  group.by = "seurat_clusters",
  return.seurat = FALSE
)

expr <- pseudobulk$RNA
expr <- log1p(expr)
# filter low variance 
gene_variance <- apply(expr, 1, var)
# keep most variable genes
expr_filtered <- expr[gene_variance > quantile(gene_variance, 0.75), ]



# Transpose for wgcna

datExpr <- t(expr_filtered)

library(WGCNA)
#Pick soft threshold:
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers)

# this plot is to pick a power for the actual wgcna thign
plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  xlab="Soft Threshold (power)",
  ylab="Scale Free Topology Model Fit (R^2)",
  type="n"
)

text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  labels=powers,
  col="red"
)

abline(h=0.8, col="red") # looks like 6 is the first number above 0.8

power_picked=6
net <- blockwiseModules(
  datExpr,
  power = power_picked,
  TOMType = "unsigned",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  verbose = 3
)
# head(net)

# WGCNA returns gene modules.
# 
# Example:
#   
#   Module	Genes
# blue	immune activation genes
# brown	metabolism genes
# yellow	cell cycle genes
# 
# These represent co-expression networks.

moduleColors <- labels2colors(net$colors)
# length(net$blockGenes)

# [1] 2 we have 2 blocks, In WGCNA, blocks are simply subsets of genes that are analyzed separately to reduce memory and computation time.

# They do NOT represent biological groups — they are just technical partitions of the gene list.
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)


length(net$colors)
length(net$blockGenes[[1]])


for(i in 1:length(net$dendrograms)){
  plotDendroAndColors(
    net$dendrograms[[i]],
    moduleColors[net$blockGenes[[i]]],
    paste("Module colors - block", i),
    dendroLabels = FALSE,
    hang = 0.03
  )
}


# 1 module eigengene = representative expression of all genes in the module

datExpr <- as.data.frame(datExpr)
datExpr <- as.matrix(datExpr)
moduleColors <- labels2colors(net$colors)
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,2,0))
#This shows:

# relationships between modules
# Interpretation:
#   
#   Modules close together → similar biological programs
# 
# Modules far apart → different biological functions


# 2 hub gene identification # these are for single colours/,modules
geneNames <- colnames(datExpr)

module <- "turquoise"

moduleGenes <- moduleColors == module

kMEturq <- kME[moduleGenes, paste0("kME", module)]

names(kMEturq) <- geneNames[moduleGenes]

hubGenes <- names(sort(kMEturq, decreasing = TRUE))[1:10]

hubGenes # genes with strongest connectivity within a module


# 3 3️⃣ Functional enrichment analysis
# what biological processes does this module represent?
library(clusterProfiler)
library(org.Hs.eg.db)


moduleGenes <- colnames(datExpr)[moduleColors == "turquoise"]
head(moduleGenes)
length(moduleGenes)

ego <- enrichGO(
  gene = moduleGenes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)

dotplot(ego)

# GO enrichment heatmap across modules lets you see which biological processes correspond to which WGCNA modules in a single figure. This is very common in gene network papers.

# modules <- unique(moduleColors)
# modules <- modules[modules != "grey"] #Usually we exclude the grey module (unassigned genes).
# modules
# GO_results <- list()
# 
# for(m in modules){
#   
#   genes <- colnames(datExpr)[moduleColors == m]
#   
#   ego <- enrichGO(
#     gene = genes,
#     OrgDb = org.Hs.eg.db,
#     keyType = "SYMBOL",
#     ont = "BP",
#     readable = TRUE
#   )
#   
#   GO_results[[m]] <- ego
# }
# 
# top_terms <- lapply(GO_results, function(x){
#   head(x@result$Description, 5)
# })
# 
# # matrix for the heatmap
# library(dplyr)
# 
# heatmap_data <- data.frame()
# 
# for(m in modules){
#   
#   res <- GO_results[[m]]@result
#   
#   if(!is.null(res)){
#     
#     tmp <- res %>%
#       dplyr::select(Description, p.adjust) %>%
#       head(10)
#     
#     tmp$module <- m
#     
#     heatmap_data <- rbind(heatmap_data, tmp)
#   }
# }
# 
# library(tidyr)
# 
# heatmap_matrix <- heatmap_data %>%
#   mutate(score = -log10(p.adjust)) %>%
#   select(module, Description, score) %>%
#   pivot_wider(names_from = module, values_from = score)
# 
# #plotheatmap
# library(pheatmap)
# 
# mat <- as.matrix(heatmap_matrix[,-1])
# rownames(mat) <- heatmap_matrix$Description
# mat[is.na(mat)] <- 0
# pheatmap(
#   mat,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE
# )

# save into a file that can be used in cytoscape
modules_of_interest <- modules[modules != "grey"] #Usually we exclude the grey module (unassigned genes).
module_df <- data.frame(
  gene_id = names(net$colors),
  colors = labels2colors(net$colors)
)

write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

head(genes_of_interest)
sum(genes_of_interest$gene_id %in% rownames(normalized_counts))
length(genes_of_interest$gene_id)
expr_of_interest = normalized_counts[genes_of_interest$gene_id,]

TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = power_picked)


row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)

write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")