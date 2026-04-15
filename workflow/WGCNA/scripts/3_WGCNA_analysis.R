library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)
library(dplyr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

setwd("./coregulation_network_analysis_project")
plot_dir   <- "plots"
output_dir <- "./"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load ──────────────────────────────────────────────────────────────────────
pb      <- readRDS("wgcna_input_matrix.rds")
pb_meta <- attr(pb, "metadata")

# Sanity checks
dim(pb)                          # must be 15 x 3000
table(pb_meta$symptomatic)       # must be 7 zeros, 8 ones

# ── Traits ────────────────────────────────────────────────────────────────────
traits <- data.frame(
  row.names   = rownames(pb),
  symptomatic = pb_meta$symptomatic[order(as.numeric(as.character(pb_meta$human)))]
)
table(traits$symptomatic)        # must be 7 zeros, 8 ones
all(!is.na(traits$symptomatic))  # must be TRUE

# ── Soft threshold ────────────────────────────────────────────────────────────
powers <- c(1:20)
sft <- pickSoftThreshold(
  pb,
  powerVector = powers,
  networkType = "signed",
  corFnc      = "bicor",
  verbose     = 5
)

png(file.path(plot_dir, "soft_threshold.png"), width = 2400, height = 1200, res = 300)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)", ylab = "Scale-free topology fit R²",
     main = "Scale independence", type = "n")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.8, col = "red", lty = 2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold (power)", ylab = "Mean connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, col = "red")
dev.off()

print(sft$fitIndices)
# !! STOP HERE - inspect soft_threshold.png and pick softPower before continuing !!

# ── Set soft power (update after inspecting plot) ─────────────────────────────
softPower <- 15  # update this after inspecting the plot

# ── Build network and detect modules ─────────────────────────────────────────
net <- blockwiseModules(
  pb,
  power             = softPower,
  networkType       = "signed",
  corType           = "bicor",
  maxPOutliers      = 0.1,
  TOMType           = "signed",
  saveTOMs          = TRUE,
  saveTOMFileBase   = file.path(output_dir, "TOM_block"),
  minModuleSize     = 30,
  deepSplit         = 2,
  mergeCutHeight    = 0.25,
  pamRespectsDendro = FALSE,
  maxBlockSize      = 10000,
  numericLabels     = FALSE,
  verbose           = 3,
  seed              = 42
)

# ── Inspect modules ───────────────────────────────────────────────────────────
table(net$colors)  # grey = unassigned, not a module

# ── Dendrogram plot ───────────────────────────────────────────────────────────
png(file.path(plot_dir, "module_dendrogram.png"), width = 3000, height = 1500, res = 300)
plotDendroAndColors(
  net$dendrograms[[1]],
  net$colors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

# ── Extract outputs ───────────────────────────────────────────────────────────
MEs <- net$MEs

gene_modules <- data.frame(
  gene   = names(net$colors),
  module = net$colors,
  stringsAsFactors = FALSE
)

gene_modules_filtered <- gene_modules[gene_modules$module != "grey", ]
colors_no_grey        <- net$colors[net$colors != "grey"]
module_genes          <- split(names(colors_no_grey), colors_no_grey)

cat("Number of modules detected:", length(module_genes), "\n")
cat("Genes assigned to modules:", nrow(gene_modules_filtered), "\n")
cat("Unassigned (grey) genes removed:", sum(net$colors == "grey"), "\n")

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(net,                   file = file.path(output_dir, "wgcna_net.rds"))
saveRDS(MEs,                   file = file.path(output_dir, "module_eigengenes.rds"))
saveRDS(module_genes,          file = file.path(output_dir, "module_gene_lists.rds"))
write.table(gene_modules_filtered,
            file      = file.path(output_dir, "gene_module_assignments.txt"),
            sep       = "\t", row.names = FALSE, quote = FALSE)

# ── Module-trait correlation ──────────────────────────────────────────────────
MEs_ordered <- orderMEs(MEs)
MEs_aligned <- MEs_ordered[rownames(MEs_ordered) != "g0", ]

# Align traits to MEs row order
traits_aligned <- traits[rownames(MEs_aligned), , drop = FALSE]
table(traits_aligned$symptomatic)        # must be 7 zeros, 8 ones
all(!is.na(traits_aligned$symptomatic))  # must be TRUE

module_trait_cor  <- cor(MEs_aligned, traits_aligned, use = "pairwise.complete.obs")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(MEs_aligned))

# ── Heatmap ───────────────────────────────────────────────────────────────────
png(file.path(plot_dir, "module_trait_heatmap.png"), width = 1200, height = 2400, res = 300)
labeledHeatmap(
  Matrix        = module_trait_cor,
  xLabels       = "Symptomatic",
  yLabels       = rownames(module_trait_cor),
  colorLabels   = FALSE,
  colors        = blueWhiteRed(50),
  textMatrix    = paste0(signif(module_trait_cor, 2),
                         "\n(p=", signif(module_trait_pval, 1), ")"),
  setStdMargins = FALSE,
  main          = "Module-trait correlation: Symptomatic vs Asymptomatic"
)
dev.off()

# ── Ranked summary table ──────────────────────────────────────────────────────
module_summary <- data.frame(
  module  = rownames(module_trait_cor),
  cor     = round(module_trait_cor[, 1], 3),
  pval    = round(module_trait_pval[, 1], 4),
  n_genes = as.numeric(table(net$colors)[gsub("ME", "", rownames(module_trait_cor))])
) %>%
  arrange(pval)

print(module_summary)

# ── Top 5 module sizes ────────────────────────────────────────────────────────
table(net$colors)[c("cyan", "pink", "red", "blue", "turquoise")]
