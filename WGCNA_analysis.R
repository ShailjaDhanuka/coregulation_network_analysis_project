library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)



obj<-readRDS("seu_object_preprocessed.rds")

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


obj@meta.data
pseudobulk <- AggregateExpression(
  obj,
  group.by = "seurat_clusters",
  return.seurat = FALSE
)

pseudobulk$RNA

expr <- pseudobulk$RNA
expr <- log1p(expr)
# filter low variance 
gene_variance <- apply(expr, 1, var)
# keep most variable genes
expr_filtered <- expr[gene_variance > quantile(gene_variance, 0.75), ]



# Transpose for wgcna

datExpr <- t(expr_filtered)


#Pick soft threshold:
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers)

# this plot is to pick a power for the actual wgcna thign

png("picking_power_for_wgcna.png", width = 2000, height = 1500, res = 300)  # high-res
# you must explicitly print the ggplot object

p<-plot(
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
print(p)
dev.off()



# abline(h=0.8, col="red") # looks like 6 is the first number above 0.8

# looks like 6 is the first number above 0.8 - dendrograms looked plateaued, not much to infer, increase threshold
power_picked=18 # pick between 16 to 20 for low sampels  

net <- blockwiseModules(
  datExpr,
  power = power_picked,
  TOMType = "unsigned",
  minModuleSize = 10,  # keep this low for small samples
  mergeCutHeight = 0.30, # slightly higher
  deepSplit      = 2, # 0-4, higher = more smaller modules
  numericLabels = TRUE,
  verbose = 3
)

#################################
# wgcna done!!
# now analysis part
##################################

  

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
# plotDendroAndColors(
#   net$dendrograms[[1]],
#   moduleColors[net$blockGenes[[1]]],
#   "Module colors",
#   dendroLabels = FALSE,
#   hang = 0.03,
#   addGuide = TRUE,
#   guideHang = 0.05
# )


length(net$colors)
length(net$blockGenes[[1]])

# plotting dendrograms ###########################3

for(i in 1:length(net$dendrograms)){
  save_name<-paste0("block",i, "_Cluster_Dendrogram.png")
  png(save_name, width = 2000, height = 1500, res = 300)  # high-res
  
  p<- plotDendroAndColors(
    title=save_name,
    net$dendrograms[[i]],
    moduleColors[net$blockGenes[[i]]],
    paste("Module colors - block", i),
    dendroLabels = FALSE,
    hang = 0.03
  )
  
  
  print(p)
  dev.off()
}

# plotting heatmaps


datExpr <- as.data.frame(datExpr)
datExpr <- as.matrix(datExpr)
moduleColors <- labels2colors(net$colors)
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,2,0))

# moduleTraitCor <- cor(MEs, datTraits, use = "p")
# moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
# 
# labeledHeatmap(
#   Matrix    = moduleTraitCor,
#   xLabels   = colnames(datTraits),
#   yLabels   = names(MEs),
#   ySymbols  = names(MEs),
#   colorLabels = FALSE,
#   colors    = blueWhiteRed(50),
#   textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
#                      signif(moduleTraitPvalue, 1), ")", sep = ""),
#   setStdMargins = FALSE,
#   cex.text  = 0.5,
#   zlim      = c(-1, 1),
#   main      = "Module-trait relationships"
# )

# 1 module eigengene = representative expression of all genes in the module

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
kME <- signedKME(datExpr, MEs)
kMEturq <- kME[moduleGenes, paste0("kME", module)]

names(kMEturq) <- geneNames[moduleGenes]

hubGenes <- names(sort(kMEturq, decreasing = TRUE))[1:10]

hubGenes # genes with strongest connectivity within a module


# 3 3️⃣ Functional enrichment analysis
# what biological processes does this module represent?



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
modules <- unique(moduleColors)
modules <- modules[modules != "grey"] #Usually we exclude the grey module (unassigned genes).
modules
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


# to be able to use it for cytoscape:so w
module_sizes <- table(moduleColors)
module_sizes <- module_sizes[names(module_sizes) != "grey"]

#  pick the best modules from here or use clusters to get the best clusters - biologically more meaningful

# # 1. See all module sizes (bigger = more robust)
# sort(table(moduleColors), decreasing = TRUE)
# 
# # 2. Get internal connectivity (intramodular connectivity) for each module
# # Higher = genes in that module are more tightly connected to each other
# connectivity <- intramodularConnectivity(adjacency(datExpr, power = 18), moduleColors)
# 
# # Average connectivity per module — higher is better
# module_connectivity <- tapply(connectivity$kWithin, moduleColors, mean)
# sort(module_connectivity, decreasing = TRUE)
# 
# # 3. Module eigengene correlation — pick modules that are DISTINCT from each other
# MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
# ME_cor <- cor(MEs)
# print(round(ME_cor, 2))
# # Avoid picking two modules with cor > 0.8 — they're redundant
# 


# these below streps are to make the traits from clusters

rownames(datExpr)
clusters <- rownames(datExpr)

datTraits <- as.data.frame(
  model.matrix(~ 0 + factor(clusters))  # creates one binary column per cluster
)
colnames(datTraits) <- gsub("factor\\(clusters\\)", "Cluster_", colnames(datTraits))
rownames(datTraits) <- clusters

datTraits
# Get modules with strongest trait correlations
# (assumes you have moduleTraitCor from your labeled heatmap step)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

labeledHeatmap(
  Matrix      = moduleTraitCor,
  xLabels     = colnames(datTraits),
  yLabels     = names(MEs),
  ySymbols    = names(MEs),
  colorLabels = FALSE,
  colors      = blueWhiteRed(50),
  textMatrix  = paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = ""),
  setStdMargins = FALSE,
  cex.text    = 0.5,
  zlim        = c(-1, 1),
  main        = "Module-cluster relationships"
)

top_modules_by_cor <- apply(abs(moduleTraitCor), 1, max) 
top_modules_by_cor <- sort(top_modules_by_cor, decreasing = TRUE)
top_modules_by_cor
# Pick modules that are: large enough + strongly correlated with a trait
meaningful_modules <- names(top_modules_by_cor[top_modules_by_cor > 0.5])
meaningful_modules <- meaningful_modules[meaningful_modules != "MEgrey"]
# strip the "ME" prefix that moduleEigengenes adds
meaningful_modules <- gsub("^ME", "", meaningful_modules)

print(meaningful_modules)
print(module_sizes[meaningful_modules])

# 
# print(meaningful_modules)
# [1] "yellow"    "green"     "blue"      "brown"     "turquoise"
# > print(module_sizes[meaningful_modules])
# moduleColors
# yellow     green      blue     brown turquoise 
# 45        41      1044        74      7880  
# turquoise is a problem, too many genes and they are similarly regulated. taking only the top 100 genes from it

turquoise_hub <- kME %>%
  filter(moduleColors == "turquoise") %>%
  arrange(desc(kMEturquoise)) %>%
  head(100) %>%
  rownames()

# Then combine with other modules
genes_of_interest <- module_df %>% 
  filter(colors %in% c("yellow", "green", "blue", "brown") |
           (colors == "turquoise" & gene_id %in% turquoise_hub))

chosen_modules<-modules

keep_genes <- genes_of_interest$gene_id

# 2. Build TOM (this may take a few minutes)
TOM <- TOMsimilarityFromExpr(datExpr, power = 18)


rownames(TOM) <- colnames(TOM) <- colnames(datExpr)

# 3. Subset TOM to only your chosen genes
TOM_sub <- TOM[keep_genes, keep_genes]

# 4. Export edge list — threshold removes weak connections
threshold <- 0.1  # adjust: higher = fewer edges = cleaner network
edge_list  <- which(TOM_sub > threshold, arr.ind = TRUE)

edges <- data.frame(
  source = rownames(TOM_sub)[edge_list[,1]],
  target = colnames(TOM_sub)[edge_list[,2]],
  weight = TOM_sub[edge_list]
) %>% 
  filter(source != target) %>%       # remove self loops
  filter(source < target)            # remove duplicate edges (A-B and B-A)

# 5. Node table with module color
nodes <- data.frame(
  id     = keep_genes,
  module = genes_of_interest$colors
)

# 6. Save
write.csv(edges, "cytoscape_edges.csv", row.names = FALSE)
write.csv(nodes, "cytoscape_nodes.csv", row.names = FALSE)

cat("Edges:", nrow(edges), "\n")
cat("Nodes:", nrow(nodes), "\n")
