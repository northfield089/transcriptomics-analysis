

expressionMatrix <- read.csv("expm.csv")
filteredExpm <- expressionMatrix[rowSums(expressionMatrix[,-1] > 10) >= 2,]
expm_t <- t(filteredExpm[,-1])
pca <- prcomp(expm_t)

pdf("PCA.pdf")
plot(
  pca$x[,1],
  pca$x[,2]
)
text(
  pca$x[,1],
  pca$x[,2],
)
dev.off()

# 
# write.csv(expressionMatrix,"expm.csv",row.names=FALSE)

# head(expressionMatrix)


# heat_data <- expressionMatrix[
#   expressionMatrix$GeneID %in% topGenes,
# ]


# write.csv(
#   heat_data,
#   "heat_data.csv"
#)