library(edgeR)
library(ggplot2)

# Load data
count_data <- read.delim("combined_counts.txt", row.names=1)
col_data <- read.delim("sample_metadata.txt", row.names=1)

# Ensure column order matches
col_data <- col_data[colnames(count_data), ]

# Create DGEList object
group <- factor(col_data$condition)
dge <- DGEList(counts=count_data, group=group)

# Normalize data
dge <- calcNormFactors(dge)

# MDS Plot
plotMDS(dge, col=as.numeric(group), main="MDS Plot")

# Estimate dispersion and fit model
dge <- estimateDisp(dge)
fit <- glmQLFit(dge, design=model.matrix(~ group))

# Perform differential expression analysis
qlf <- glmQLFTest(fit, coef=2)
top_genes <- topTags(qlf, n=30)$table

# Save results
write.table(qlf$table, "edgeR_results.txt", sep="\t", quote=FALSE)
write.table(top_genes, "edgeR_top30_genes.txt", sep="\t", quote=FALSE)

