#Visualizing the results of WGCNA saved in plots
# hub genes per modules saved as hub_genes_per_modules in data



library(pheatmap)



setwd("./coregulation_network_analysis_project")

plot_dir <- "workflow/WGCNA/plots"
output_dir <- "data/WGCNA"



# Module eigengene heatmap ──────────────────────────────────────────
# Shows how each module's activity varies across your pseudobulk clusters
# This is one of the most biologically interpretable WGCNA plots


# Clean up ME column names (remove "ME" prefix for readability)
ME_matrix <- readRDS("data/WGCNA/module_eigengenes.rds")
colnames(ME_matrix) <- gsub("^ME", "", colnames(ME_matrix))

# Remove grey if present
ME_matrix <- ME_matrix[, colnames(ME_matrix) != "grey"]

dev.off()  # run this once or twice until you get "null device"
graphics.off()  # closes all devices to be safe

png(file.path(plot_dir, "eigengene_heatmap.png"), width = 2400, height = 2000, res = 300)
pheatmap(t(ME_matrix),
         cluster_rows    = TRUE,
         cluster_cols    = TRUE,
         show_colnames   = TRUE,
         fontsize        = 8,
         color           = colorRampPalette(c("blue", "white", "red"))(100),
         main            = "Module eigengene expression across pseudobulk clusters")
dev.off()



# Eigengene correlation network ─────────────────────────────────────
# Shows which modules are transcriptionally related to each other
png(file.path(plot_dir, "eigengene_correlation.png"), width = 2000, height = 2000, res = 300)
plotEigengeneNetworks(
  ME_matrix,
  "Module eigengene network",
  marHeatmap      = c(3, 4, 2, 2),
  plotDendrograms = TRUE,
  xLabelsAngle    = 90
)
dev.off()



# Gene module membership vs. intramodular connectivity ───────────────
# Hub genes sit in the top right: high connectivity AND high module membership
# These are your candidate hub genes for Cytoscape

# Compute intramodular connectivity
pb   <- readRDS("data/WGCNA/wgcna_input_matrix.rds")
net        <- readRDS("data/WGCNA/wgcna_net.rds")
MEs        <- readRDS("data/WGCNA/module_eigengenes.rds")
module_genes <- readRDS("data/WGCNA/module_gene_lists.rds")

# Recreate colors_no_grey and softPower if not in environment
softPower      <- 15
connectivity <- intramodularConnectivity.fromExpr(
  pb,
  colors      = net$colors,
  networkType = "signed",
  power       = softPower,
  corFnc      = "bicor"
)

# Assign names explicitly from pb_5000 column order
rownames(connectivity) <- colnames(pb)
# Spot check: manually correlate one gene with its module eigengene
# and compare to what connectivity2 reports
# test_gene   <- colnames(pb_5000)[1]
# test_module <- net$colors[test_gene]
# cat("Test gene:", test_gene, "\n")
# cat("Its module:", test_module, "\n")
# cat("kWithin from connectivity2:", connectivity2[test_gene, "kWithin"], "\n")
# 
# # Cross-check: does this gene appear in module_genes for that module?
# cat("Gene in module list:", test_gene %in% module_genes[[test_module]], "\n")

# Save hub gene table for each module ───────────────────────────────────────
MM_list <- lapply(names(module_genes), function(mod) {
  ME_col       <- paste0("ME", mod)
  if (!ME_col %in% colnames(MEs)) return(NULL)
  
  genes        <- intersect(module_genes[[mod]], colnames(pb))
  mm_vals      <- cor(pb[, genes], MEs[, ME_col], use = "p")
  
  data.frame(
    gene             = genes,
    module           = mod,
    module_membership = as.numeric(mm_vals),
    stringsAsFactors = FALSE
  )
})

MM_df <- do.call(rbind, MM_list)

# ── Merge with connectivity ────────────────────────────────────────────────────
conn_df <- data.frame(
  gene    = rownames(connectivity),
  kWithin = connectivity$kWithin,
  kTotal  = connectivity$kTotal,
  stringsAsFactors = FALSE
)

gene_stats <- merge(MM_df, conn_df, by = "gene")

# ── Apply dual threshold to define hub genes ──────────────────────────────────
# MM > 0.8 AND kWithin in top 25% of its module
# This keeps biologically meaningful genes without arbitrary top N

gene_stats <- gene_stats %>%
  group_by(module) %>%
  mutate(
    kWithin_quantile = percent_rank(kWithin),   # 0-1, higher = more connected
    is_hub = module_membership > 0.7 & kWithin_quantile > 0.75
  ) %>%
  ungroup()

cat("Total genes:", nrow(gene_stats), "\n")
cat("Hub genes (MM>0.8 & top 25% connectivity):", sum(gene_stats$is_hub), "\n")

# Hub genes per module
gene_stats %>%
  filter(is_hub) %>%
  group_by(module) %>%
  summarise(n_hubs = n()) %>%
  arrange(desc(n_hubs)) %>%
  print()

# ── Save full gene stats table ─────────────────────────────────────────────────
write.table(gene_stats,
            file      = file.path(output_dir, "gene_module_stats.txt"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# Save hub genes separately
hub_gene_table <- gene_stats %>% filter(is_hub) %>% arrange(module, desc(kWithin))

write.table(hub_gene_table,
            file      = file.path(output_dir, "hub_genes_per_module.txt"),
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

cat("Hub gene table saved:", nrow(hub_gene_table), "genes\n")


# visualise per module hub genes
png(file.path(plot_dir, "hub_gene_scatter.png"), width = 3000, height = 3000, res = 300)
par(mfrow = c(ceiling(sqrt(length(module_genes))),
              ceiling(sqrt(length(module_genes)))))

for (mod in names(module_genes)) {
  
  ME_col <- paste0("ME", mod)
  if (!ME_col %in% colnames(MEs)) next
  
  genes_in_mod <- module_genes[[mod]]
  genes_in_mod <- intersect(genes_in_mod, colnames(pb))
  genes_in_mod <- intersect(genes_in_mod, rownames(connectivity))
  
  if (length(genes_in_mod) < 3) next
  
  mm <- cor(pb[, genes_in_mod], MEs[, ME_col], use = "p")
  kw <- connectivity[genes_in_mod, "kWithin"]
  
  if (all(!is.finite(kw)) | all(!is.finite(mm))) next
  
  gene_is_hub  <- gene_stats$is_hub[match(genes_in_mod, gene_stats$gene)]
  point_colors <- ifelse(gene_is_hub, "red", mod)
  
  plot(kw, mm,
       col  = point_colors,
       pch  = 20,
       xlab = "Intramodular connectivity",
       ylab = "Module membership",
       main = paste0("Module: ", mod, " (n=", length(genes_in_mod), ")"),
       cex  = 0.6)
  
  abline(h = 0.7, col = "red",  lty = 2, lwd = 0.8)  # MM threshold updated to 0.7
  abline(v = quantile(kw, 0.75), col = "blue", lty = 2, lwd = 0.8)
  abline(lm(mm ~ kw), col = "black", lty = 1, lwd = 0.8)
}

dev.off()
cat("Plot saved to:", file.path(plot_dir, "hub_gene_scatter.png"), "\n")

###########################################################################################
