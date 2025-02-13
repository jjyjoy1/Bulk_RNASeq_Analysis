library(limma)
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

# Normalize and apply voom transformation
dge <- calcNormFactors(dge)
design <- model.matrix(~ group)
v <- voom(dge, design, plot=TRUE)

# Fit linear model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Get results
results <- topTable(fit, coef=2, adjust="fdr", n=Inf)
top_30_genes <- head(results, 30)

# Save results
write.table(results, "limma_voom_results.txt", sep="\t", quote=FALSE)
write.table(top_30_genes, "limma_voom_top30_genes.txt", sep="\t", quote=FALSE)



