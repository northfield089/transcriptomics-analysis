counts1 <- read.table(
	"SRR390728_blast/SRR390728_counts.txt",
	header=TRUE,
	comment.char = "#"
)
counts2 <- read.table(
        "SRR390729_blast/SRR390729_counts.txt",
        header=TRUE,
        comment.char = "#"
)
counts3 <- read.table(
        "SRR390730_blast/SRR390730_counts.txt",
        header=TRUE,
        comment.char = "#"
)
counts4 <- read.table(
        "SRR390731_blast/SRR390731_counts.txt",
        header=TRUE,
        comment.char = "#"
)
expressionMatrix <- data.frame(
	GeneID = counts1$Geneid,
	SRR390728 = counts1[,7],
	SRR390729 = counts2[,7],
	SRR390730 = counts3[,7],
	SRR390731 = counts4[,7]
)

#write.csv(expm ,"expm.csv",row.names = F)

library(DESeq2)
group <- factor(c(
	"Control",
	"Control",
	"Case",
	"Case"
))

DESeq2DataSet <- DESeqDataSetFromMatrix(		 #DESeq2分析的源数据
	countData = expressionMatrix[,-1],
	colData = data.frame(group),
	design = ~ group
)

rownames(DESeq2DataSet) <- expressionMatrix$GeneID

DESeq2DataSet <- DESeq(DESeq2DataSet)

result <- results(DESeq2DataSet)	#取出dds里的差异基因结果

sigRes = subset(			#筛选padj<0.05的
	result,
	padj < 0.05,

)

write.csv(
	sigRes,
	"DEG_results.csv"
)

plot(					#火山图
	sigRes$log2FoldChange,
	-log10(sigRes$padj),
	pch = 20			#点形状第20号
)

#按padj排序
#前50个基因
topGenes <- head(
    rownames(sigRes[order(sigRes$padj), ]),
    50
)


heat_data <- expressionMatrix[
    expressionMatrix$GeneID %in% topGenes,
]
heat_data <- heat_data[,-1]	#去掉基因名一列
library(pheatmap)

pheatmap(as.maheat_data)		#画热图

help(package = "pheatmap")



