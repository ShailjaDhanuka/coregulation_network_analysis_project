library(pheatmap)
#load viper matrix, edit viper_matrix to your csv file name for VIPER
viper_matrix <- read.csv("data/VIPER/red_viper_activity_scores.csv", row.names = 1)
viper_matrix <- as.matrix(viper_matrix)
mode(viper_matrix) <- "numeric"
viper_matrix <- viper_matrix[, paste0("g", 1:15)]

symptomatic <- c(0,0,0,1,0,1,1,1,0,1,1,1,0,0,1)
names(symptomatic) <- colnames(viper_matrix)

stopifnot(all(names(symptomatic) == colnames(viper_matrix)))

#Asymptomatic (0) vs Symptomatic (1)
group0 <- viper_matrix[, symptomatic == 0]
group1 <- viper_matrix[, symptomatic == 1]

row_var <- apply(viper_matrix, 1, var)
top_idx <- order(row_var, decreasing = TRUE)[1:30]

top_matrix <- viper_matrix[top_idx, ]
#Remove rows with zero variance
top_matrix <- top_matrix[apply(top_matrix, 1, sd) > 0, ]

scaled_matrix <- t(scale(t(top_matrix)))

#Remove the sparse rows
scaled_matrix <- scaled_matrix[
  apply(scaled_matrix, 1, function(x) all(is.finite(x))),
]

scaled_matrix[is.na(scaled_matrix)] <- 0
scaled_matrix[is.nan(scaled_matrix)] <- 0
scaled_matrix[is.infinite(scaled_matrix)] <- 0

annotation <- data.frame(
  Condition = factor(symptomatic)
)
rownames(annotation) <- colnames(viper_matrix)

ord <- order(symptomatic)

viper_matrix <- viper_matrix[, ord]
scaled_matrix <- scaled_matrix[, ord]

# Rebuild annotation to match new order
annotation <- data.frame(
  Regulator = factor(symptomatic[ord])
)
rownames(annotation) <- colnames(viper_matrix)

pheatmap(
  scaled_matrix,
  annotation_col = annotation,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_cols = FALSE,
  gaps_col = sum(symptomatic == 0),
  main = "VIPER Activity Heatmap (Asymptomatic vs Symptomatic)",
  filename = "data/VIPER/heatmap_outputs/viper_heatmap.png",
  width = 10,
  height = 8
)

#pheatmap(
#  t(scale(t(group0))),
#  color = colorRampPalette(c("blue", "white", "red"))(100),
#  main = "Asymptomatic (0)"
#)

#pheatmap(
#  t(scale(t(group1))),
#  color = colorRampPalette(c("blue", "white", "red"))(100),
#  main = "Symptomatic (1)"
#)
