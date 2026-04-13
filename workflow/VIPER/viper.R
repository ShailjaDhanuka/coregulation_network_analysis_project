library(Seurat)
library(viper)
#Read TF table
df <- read.table("black_aracne.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

df <- df[df$MI > 1.0, ]

#Read expression matrix
seurat_obj <- readRDS("pseudobulked_seu_obj.rds")
expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

expr <- as.matrix(expr)
mode(expr) <- "numeric"

#Convert Arcane network → regulon
regulon_list <- lapply(split(df, df$Regulator), function(x) {
  targets <- x$Target
  weights <- x$MI
  list(
    tfmode = setNames(rep(1, length(targets)), targets),
    likelihood = setNames(weights, targets)
  )
})

#Filter weak edges
regulon_list <- lapply(regulon_list, function(reg) {
  top <- order(reg$likelihood, decreasing = TRUE)[seq_len(min(50, length(reg$likelihood)))]  # Keep top 50 targets
  list(
    tfmode = reg$tfmode[top],
    likelihood = reg$likelihood[top]
  )
})

#Run VIPER analysis
vpres <- viper(expr, regulon_list, method = "scale")