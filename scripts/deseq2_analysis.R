library(DESeq2)
library(ggplot2)

# Load combined counts and metadata
count_data <- read.delim("combined_counts.txt", row.names=1)
col_data <- read.delim("sample_metadata.txt", row.names=1)

# Ensure column order matches
col_data <- col_data[colnames(count_data), ]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

# PCA plot
vsd <- vst(dds, blind=FALSE)
pca_data <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
ggplot(pca_data, aes(x=PC1, y=PC2, color=condition)) + geom_point(size=3) + ggtitle("PCA of RNA-Seq Samples")

# Get results
res <- results(dds, alpha=0.05, lfcThreshold=1)
res <- res[order(res$padj, na.last=NA), ]
top_30_genes <- head(res, 30)

# Save results
write.table(res, "DESeq2_results.txt", sep="\t", quote=FALSE)
write.table(top_30_genes, "DESeq2_top30_genes.txt", sep="\t", quote=FALSE)


