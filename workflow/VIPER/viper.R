library(viper)
library(dplyr)

outdir <- "data/VIPER"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", outdir, "\n\n")

cat("=== Loading expression matrix ===\n")

mat <- readRDS("data/WGCNA/wgcna_input_matrix.rds")
cat("Raw matrix dimensions:", dim(mat), "\n")
cat("Raw rownames (samples):", head(rownames(mat)), "\n")
cat("Raw colnames (genes):", head(colnames(mat)), "\n")

expr_mat <- t(mat)
cat("Transposed matrix dimensions:", dim(expr_mat), "\n")

cat("\n=== Loading ARACNe network ===\n")

#Put aracne txt here!!! YOU CAN CHANGE IT TO YOUR ARACNE FILE NAME!
raw <- read.delim("data/ARACNe/red_aracne.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(raw) <- c("tf", "target", "mi")
raw <- raw[complete.cases(raw), ]

cat("Total edges in network:", nrow(raw), "\n")
cat("Unique regulators:", length(unique(raw$tf)), "\n")
cat("Unique targets:", length(unique(raw$target)), "\n")

cat("\n=== Checking gene ID overlap ===\n")

aracne_genes <- unique(c(raw$tf, raw$target))
matrix_genes <- rownames(expr_mat)

cat("ARACNe unique genes:", length(aracne_genes), "\n")
cat("Matrix unique genes:", length(matrix_genes), "\n")
cat("Overlapping genes:  ", length(intersect(aracne_genes, matrix_genes)), "\n")

missing_from_matrix <- setdiff(aracne_genes, matrix_genes)
cat("\nARACNe genes NOT found in matrix (first 20):\n")
print(head(missing_from_matrix, 20))

cat("\nMatrix gene examples (first 20):\n")
print(head(matrix_genes, 20))

overlap_n <- length(intersect(aracne_genes, matrix_genes))
if (overlap_n == 0) {
  stop(paste0(
    "\nNo gene overlap between ARACNe network and expression matrix.\n",
    "Check that both files use the same gene ID format (symbols vs Ensembl).\n"
  ))
} else if (overlap_n < 10) {
  warning(paste0("Very low overlap (", overlap_n, " genes). Results may be unreliable."))
}

cat("\n=== Building regulon object ===\n")

regulon <- lapply(split(raw, raw$tf), function(df) {
  tf_name <- df$tf[1]
  targets <- df$target

  targets <- targets[targets %in% rownames(expr_mat)]
  if (length(targets) == 0) return(NULL)

  if (tf_name %in% rownames(expr_mat)) {
    tf_expr <- as.numeric(expr_mat[tf_name, ])
    mor <- sapply(targets, function(g) {
      cor(tf_expr, as.numeric(expr_mat[g, ]), use = "pairwise.complete.obs")
    })
    mor[is.na(mor)] <- 0
  } else {
    mor <- rep(1, length(targets))
    names(mor) <- targets
    cat("  TF not in matrix, assuming activation:", tf_name, "\n")
  }

  list(
    tfmode     = setNames(mor, targets),
    likelihood = rep(1, length(targets))
  )
})

regulon <- Filter(Negate(is.null), regulon)
cat("Regulators retained after overlap filter:", length(regulon), "\n")

if (length(regulon) == 0) {
  stop("No regulators survived the overlap filter. Check gene ID consistency.")
}

cat("\n=== Running VIPER ===\n")

viper_result <- viper(
  eset       = expr_mat,
  regulon    = regulon,
  minsize    = 2,
  method     = "none",
  verbose    = TRUE
)

cat("\nVIPER result dimensions:", dim(viper_result), "\n")
cat("(rows = TFs/regulators, columns = samples)\n")

cat("\n=== Top regulators by mean absolute activity ===\n")

mean_abs_nes <- sort(rowMeans(abs(viper_result)), decreasing = TRUE)
print(head(mean_abs_nes, 20))

cat("\n=== Saving results ===\n")

write.csv(
  as.data.frame(viper_result),
  file.path("data/VIPER", "viper_activity_scores.csv"),
  quote = FALSE
)
cat("Saved to:", file.path(outdir, "viper_activity_scores.csv"), "\n")

cat("\n=== VIPER pipeline complete ===\n")
