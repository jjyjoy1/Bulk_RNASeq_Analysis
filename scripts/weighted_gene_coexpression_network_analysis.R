# Load required libraries
library(WGCNA)

# Load count data (genes as rows, samples as columns)
count_data <- read.delim("combined_counts.txt", row.names=1)
count_data <- t(count_data)  # Transpose for WGCNA format

# Sample metadata
metadata <- read.delim("sample_metadata.txt", row.names=1)

# Check for missing values
gsg <- goodSamplesGenes(count_data, verbose = 3)
if (!gsg$allOK) {
    count_data <- count_data[gsg$goodSamples, gsg$goodGenes]
}

# Soft threshold power selection
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(count_data, powerVector = powers, verbose = 5)

# Choose power = 6 (adjust based on sft results)
softPower <- 6
adjacency <- adjacency(count_data, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Hierarchical clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Dynamic Tree Cut
moduleColors <- dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE)
moduleColors <- labels2colors(moduleColors)

# Module-trait relationships
MEs <- moduleEigengenes(count_data, colors = moduleColors)$eigengenes
moduleTraitCor <- cor(MEs, metadata$condition, use = "p")

# Save Module-Trait Relationship Data
write.table(moduleTraitCor, "module_trait_relationships.txt", sep="\t", quote=FALSE)

# Cytoscape Export
cyt <- exportNetworkToCytoscape(TOM, filenameBase = "wgcna_network", threshold = 0.1)
write.table(cyt$edgeData, "wgcna_edges.txt", sep="\t", quote=FALSE)
write.table(cyt$nodeData, "wgcna_nodes.txt", sep="\t", quote=FALSE)

