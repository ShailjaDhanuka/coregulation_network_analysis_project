library(Seurat)

# Load Seurat object
pb_obj <- readRDS("pseudobulked_seu_obj.rds")

expr <- as.matrix(
  GetAssayData(readRDS("pseudobulked_seu_obj.rds"), assay = "RNA", layer = "data")
)

# Ensure rownames exist
if(is.null(rownames(expr_mat))) {
  rownames(expr_mat) <- paste0("Gene", 1:nrow(expr_mat))
}

# Filter genes with variance > 0
expr <- expr[apply(expr, 1, sd) > 0, ]

write.table(expr, "pseudobulk_expression_filtered.tsv",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)