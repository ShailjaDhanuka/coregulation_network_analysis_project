library(WGCNA)
library(pheatmap)
library(dplyr)
options(stringsAsFactors = FALSE)

setwd("./coregulation_network_analysis_project")
plot_dir   <- "plots"
output_dir <- "./"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load everything ───────────────────────────────────────────────────────────
pb           <- readRDS("wgcna_input_matrix.rds")
pb_meta      <- attr(pb, "metadata")
net          <- readRDS("wgcna_net.rds")
MEs          <- readRDS("module_eigengenes.rds")
module_genes <- readRDS("module_gene_lists.rds")
softPower    <- 15  # update if you changed this

# Rebuild traits
traits <- data.frame(
  row.names   = rownames(pb),
  symptomatic = pb_meta$symptomatic[order(as.numeric(as.character(pb_meta$human)))]
)

# ── Module eigengene heatmap ──────────────────────────────────────────────────
ME_matrix <- MEs
colnames(ME_matrix) <- gsub("^ME", "", colnames(ME_matrix))
ME_matrix <- ME_matrix[, colnames(ME_matrix) != "grey"]

# Add symptomatic annotation for columns
sample_annotation <- data.frame(
  row.names   = rownames(ME_matrix),
  symptomatic = factor(traits$symptomatic, labels = c("Asymptomatic", "Symptomatic"))
)

png(file.path(plot_dir, "eigengene_heatmap.png"), width = 2400, height = 2000, res = 300)
pheatmap(t(ME_matrix),
         cluster_rows    = TRUE,
         cluster_cols    = TRUE,
         show_colnames   = TRUE,
         annotation_col  = sample_annotation,
         fontsize        = 8,
         color           = colorRampPalette(c("blue", "white", "red"))(100),
         main            = "Module eigengene expression across donors")
dev.off()

# ── Eigengene correlation network ─────────────────────────────────────────────
png(file.path(plot_dir, "eigengene_correlation.png"), width = 2000, height = 2000, res = 300)
plotEigengeneNetworks(
  ME_matrix,
  "Module eigengene network",
  marHeatmap      = c(3, 4, 2, 2),
  plotDendrograms = TRUE,
  xLabelsAngle    = 90
)
dev.off()

# ── Intramodular connectivity ─────────────────────────────────────────────────
connectivity <- intramodularConnectivity.fromExpr(
  pb,
  colors      = net$colors,
  networkType = "signed",
  power       = softPower,
  corFnc      = "bicor"
)
rownames(connectivity) <- colnames(pb)

# ── Module membership ─────────────────────────────────────────────────────────
MM_list <- lapply(names(module_genes), function(mod) {
  ME_col <- paste0("ME", mod)
  if (!ME_col %in% colnames(MEs)) return(NULL)
  
  genes  <- intersect(module_genes[[mod]], colnames(pb))
  if (length(genes) == 0) return(NULL)
  
  mm_vals <- cor(pb[, genes, drop = FALSE], MEs[, ME_col], use = "pairwise.complete.obs")
  
  data.frame(
    gene              = genes,
    module            = mod,
    module_membership = as.numeric(mm_vals),
    stringsAsFactors  = FALSE
  )
})

MM_df <- do.call(rbind, Filter(Negate(is.null), MM_list))

# ── Merge with connectivity ───────────────────────────────────────────────────
conn_df <- data.frame(
  gene    = rownames(connectivity),
  kWithin = connectivity$kWithin,
  kTotal  = connectivity$kTotal,
  stringsAsFactors = FALSE
)

gene_stats <- merge(MM_df, conn_df, by = "gene")

# ── Define hub genes ──────────────────────────────────────────────────────────
gene_stats <- gene_stats %>%
  group_by(module) %>%
  mutate(
    kWithin_quantile = percent_rank(kWithin),
    is_hub           = module_membership > 0.7 & kWithin_quantile > 0.75
  ) %>%
  ungroup()

cat("Total genes:", nrow(gene_stats), "\n")
cat("Hub genes (MM>0.7 & top 25% connectivity):", sum(gene_stats$is_hub), "\n")

gene_stats %>%
  filter(is_hub) %>%
  group_by(module) %>%
  summarise(n_hubs = n()) %>%
  arrange(desc(n_hubs)) %>%
  print()

# ── Save gene stats and hub genes ────────────────────────────────────────────
write.table(gene_stats,
            file      = file.path(output_dir, "gene_module_stats.txt"),
            sep       = "\t", quote = FALSE, row.names = FALSE)

hub_gene_table <- gene_stats %>%
  filter(is_hub) %>%
  arrange(module, desc(kWithin))

write.table(hub_gene_table,
            file      = file.path(output_dir, "hub_genes_per_module.txt"),
            sep       = "\t", quote = FALSE, row.names = FALSE)

cat("Hub genes saved:", nrow(hub_gene_table), "genes\n")

# ── Hub gene scatter plots per module ────────────────────────────────────────
n_mods <- length(module_genes)
n_cols  <- ceiling(sqrt(n_mods))
n_rows  <- ceiling(n_mods / n_cols)

png(file.path(plot_dir, "hub_gene_scatter.png"), width = 3000, height = 3000, res = 300)
par(mfrow = c(n_rows, n_cols))

for (mod in names(module_genes)) {
  ME_col <- paste0("ME", mod)
  if (!ME_col %in% colnames(MEs)) next
  
  genes_in_mod <- intersect(module_genes[[mod]], colnames(pb))
  genes_in_mod <- intersect(genes_in_mod, rownames(connectivity))
  if (length(genes_in_mod) < 3) next
  
  mm <- cor(pb[, genes_in_mod, drop = FALSE], MEs[, ME_col], use = "pairwise.complete.obs")
  kw <- connectivity[genes_in_mod, "kWithin"]
  
  if (!any(is.finite(kw)) || !any(is.finite(mm))) next
  
  gene_is_hub  <- gene_stats$is_hub[match(genes_in_mod, gene_stats$gene)]
  point_colors <- ifelse(gene_is_hub, "red", mod)
  
  plot(kw, mm,
       col  = point_colors,
       pch  = 20,
       xlab = "Intramodular connectivity",
       ylab = "Module membership",
       main = paste0(mod, " (n=", length(genes_in_mod), ")"),
       cex  = 0.6)
  
  abline(h = 0.7,                    col = "red",   lty = 2, lwd = 0.8)
  abline(v = quantile(kw, 0.75),     col = "blue",  lty = 2, lwd = 0.8)
  abline(lm(mm ~ kw),                col = "black", lty = 1, lwd = 0.8)
}

dev.off()
cat("Hub gene scatter saved to:", file.path(plot_dir, "hub_gene_scatter.png"), "\n")

