# End-to-End RNAseq Analysis with WGCNA and Bayesian Networks
# This workflow assumes you have a normalized expression matrix ready
# Format: rows = genes, columns = samples

#===============================================================================
# PART 1: SETUP AND LIBRARIES
#===============================================================================

# Install required packages if not already installed
required_packages <- c("WGCNA", "bnlearn", "Rgraphviz", "clusterProfiler", 
                      "org.Hs.eg.db", "flashClust", "impute", "preprocessCore", 
                      "dynamicTreeCut", "limma", "ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

# Load libraries
library(WGCNA)
library(bnlearn)
library(Rgraphviz)
library(clusterProfiler)
library(org.Hs.eg.db)
library(flashClust)
library(preprocessCore)
library(dynamicTreeCut)
library(limma)
library(ggplot2)

# Enable multi-threading for WGCNA
enableWGCNAThreads()

#===============================================================================
# PART 2: DATA PREPARATION 
#===============================================================================

# Load your normalized expression matrix
# Replace with your actual file path
# Format should be: rows = genes, columns = samples
# Assuming data is already log-transformed and normalized

# Example loading data (replace with your actual data)
# expr_data <- read.csv("path/to/normalized_expression_matrix.csv", row.names=1)

# For this example, let's create a dummy dataset
set.seed(12345)
gene_count <- 5000
sample_count <- 30

# Create gene names
gene_names <- paste0("GENE", 1:gene_count)

# Create sample names
sample_names <- paste0("SAMPLE", 1:sample_count)

# Sample traits (e.g., disease status, treatment groups)
# 0 = control, 1 = treatment
sample_traits <- data.frame(
  treatment = c(rep(0, sample_count/2), rep(1, sample_count/2)),
  age = round(runif(sample_count, 20, 80)),
  sex = sample(c("M", "F"), sample_count, replace = TRUE)
)
rownames(sample_traits) <- sample_names

# Create random expression data
expr_data <- matrix(rnorm(gene_count * sample_count), nrow = gene_count)
rownames(expr_data) <- gene_names
colnames(expr_data) <- sample_names

# Check for missing values
if(sum(is.na(expr_data)) > 0) {
  # Impute missing values
  expr_data <- impute.knn(expr_data)$data
}

# Transpose for WGCNA (WGCNA expects samples in rows, genes in columns)
expr_data_transposed <- t(expr_data)

#===============================================================================
# PART 3: WGCNA ANALYSIS
#===============================================================================

# Step 1: Check for outlier samples
sample_tree <- flashClust(dist(expr_data_transposed), method = "average")
plot(sample_tree, main = "Sample clustering to detect outliers", 
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Optional: Remove outlier samples if necessary
# For this example, we'll assume no outliers
good_samples <- rownames(expr_data_transposed)
expr_data_cleaned <- expr_data_transposed[good_samples, ]

# Step 2: Choose soft-thresholding power (scale-free topology criterion)
powers <- c(1:20)
sft <- pickSoftThreshold(expr_data_cleaned, powerVector = powers, 
                         verbose = 5, networkType = "signed")

# Plot the results
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale independence")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     main = "Mean connectivity")

# Step 3: Choose the optimal power based on the plot
# For this example, let's say the optimal power is 6
selected_power <- 6

# Step 4: Construct the co-expression network
net <- blockwiseModules(expr_data_cleaned, 
                        power = selected_power,
                        TOMType = "signed",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM",
                        verbose = 3)

# Step 5: Convert numeric labels to colors for easier visualization
module_colors <- labels2colors(net$colors)

# Step 6: Plot the dendrogram and module colors
plotDendroAndColors(net$dendrograms[[1]], 
                    module_colors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05)

# Step 7: Get module-trait relationships
# Calculate module eigengenes
MEs <- net$MEs
# Order by module labels
MEs_ordered <- MEs[, order(colnames(MEs))]
# Add treatment information
MEs_with_traits <- cbind(MEs_ordered, sample_traits)

# Calculate correlations between module eigengenes and traits
module_trait_corr <- cor(MEs_ordered, sample_traits[, 1, drop = FALSE], use = "p")
module_trait_p <- corPvalueStudent(module_trait_corr, nrow(expr_data_cleaned))

# Visualize module-trait correlations
textMatrix <- paste(signif(module_trait_corr, 2), "\n(", signif(module_trait_p, 1), ")", sep = "")
dim(textMatrix) <- dim(module_trait_corr)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = module_trait_corr,
               xLabels = colnames(module_trait_corr),
               yLabels = colnames(MEs_ordered),
               ySymbols = colnames(MEs_ordered),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Step 8: Find genes with high module membership and trait correlation
# Choose a module of interest (e.g., the one with highest correlation to treatment)
module_of_interest <- which.max(abs(module_trait_corr[, "treatment"]))
module_color <- colnames(MEs_ordered)[module_of_interest]
column_idx <- match(module_color, colnames(MEs_ordered))
module_genes <- module_colors == gsub("ME", "", module_color)

# Calculate gene significance and module membership
gene_significance <- abs(cor(expr_data_cleaned, sample_traits[, "treatment", drop = FALSE], use = "p"))
module_membership <- abs(cor(expr_data_cleaned, MEs_ordered[, column_idx], use = "p"))

# Create data frame with gene info
gene_info <- data.frame(
  gene_name = colnames(expr_data_cleaned),
  module_color = module_colors,
  gene_significance = gene_significance,
  module_membership = module_membership[, 1]
)

# Filter for genes in the module of interest
module_gene_info <- subset(gene_info, module_color == gsub("ME", "", module_color))

# Sort by module membership and gene significance
module_gene_info <- module_gene_info[order(-module_gene_info$module_membership, 
                                           -module_gene_info$gene_significance), ]

# Top hub genes in the module (high module membership)
top_hub_genes <- head(module_gene_info, 20)

# Step 9: Export the gene list for the module of interest for enrichment analysis
write.csv(module_gene_info, file = "module_genes_for_enrichment.csv", row.names = FALSE)

# Get gene lists for all modules
module_gene_lists <- list()
for (color in unique(module_colors)) {
  module_gene_lists[[color]] <- colnames(expr_data_cleaned)[module_colors == color]
}

#===============================================================================
# PART 4: BAYESIAN NETWORK ANALYSIS ON WGCNA MODULES 
#===============================================================================

# Focus Bayesian Network analysis on genes within the module of interest
# For computational reasons, we'll use a subset of the top genes in the module

# Get top 50 genes from the module based on module membership
top_module_genes <- head(module_gene_info, 50)$gene_name

# Extract expression data for these genes
# Use original orientation (genes in rows) and transpose to samples in rows
module_expr <- t(expr_data[top_module_genes, ])

# Step 1: Discretize the data for Bayesian Network (optional but common)
# Here we use 3 levels: low, medium, high expression
discretized_data <- data.frame(apply(module_expr, 2, function(x) {
  cut(x, breaks = quantile(x, probs = c(0, 0.33, 0.66, 1)), 
      labels = c("low", "medium", "high"), include.lowest = TRUE)
}))

# Step 2: Learn the Bayesian Network structure
# Use hill-climbing algorithm with BIC score
bn_structure <- hc(discretized_data, score = "bic")

# Step 3: Learn the parameters of the network
bn_model <- bn.fit(bn_structure, discretized_data)

# Step 4: Visualize the network
# Plot the network structure
plot(bn_structure, main = "Bayesian Network of Top Module Genes")

# Step 5: Find key regulatory genes (nodes with many outgoing edges)
node_strength <- strength(bn_structure)
key_regulators <- node_strength[order(-node_strength$to), ]
head(key_regulators, 10)  # Top 10 potential regulatory genes

# Step 6: Predict effects of gene perturbations
# Example: What happens if we set a key regulator gene to "high" expression?
# Identify a key regulator gene
top_regulator <- rownames(key_regulators)[1]

# Create evidence for the prediction (set the key regulator to "high")
evidence <- list()
evidence[[top_regulator]] <- "high"

# Predict the effects on other genes
predicted_effects <- predict(bn_model, node = names(discretized_data), 
                             evidence = evidence, method = "bayes-lw")

# Compare predicted vs. observed frequencies for each gene
# This shows how forcing the key regulator to "high" affects other genes
for (gene in names(discretized_data)) {
  if (gene != top_regulator) {
    cat("Gene:", gene, "\n")
    observed_table <- table(discretized_data[[gene]])
    predicted_table <- table(predicted_effects[[gene]])
    observed_prop <- prop.table(observed_table)
    predicted_prop <- prop.table(predicted_table)
    
    cat("  Observed frequencies:", "\n")
    print(observed_prop)
    cat("  Predicted frequencies after setting", top_regulator, "to high:", "\n")
    print(predicted_prop)
    cat("\n")
  }
}

#===============================================================================
# PART 5: FUNCTIONAL ENRICHMENT ANALYSIS OF WGCNA MODULES
#===============================================================================

# Perform enrichment analysis for each module
# Requires gene IDs to be in a format recognized by clusterProfiler (e.g., Entrez, Symbol)
# For this example, we'll assume our gene names are gene symbols

# Function to perform GO enrichment on a gene list
perform_GO_enrichment <- function(gene_list, ont = "BP") {
  enrich_result <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",  # Change based on your gene ID type
    ont = ont,  # BP = Biological Process, MF = Molecular Function, CC = Cellular Component
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  return(enrich_result)
}

# For example, enrichment analysis for the module of interest
module_genes <- module_gene_lists[[gsub("ME", "", module_color)]]

# For this example, we'll assume our dummy gene names need to be converted
# In a real scenario, you would use your actual gene symbols
# Here we'll just simulate some real gene names for demonstration
set.seed(123)
real_genes <- sample(keys(org.Hs.eg.db, keytype = "SYMBOL"), length(module_genes))

# Perform GO Biological Process enrichment
bp_enrichment <- perform_GO_enrichment(real_genes, ont = "BP")

# Visualize top 10 enriched terms
barplot(bp_enrichment, showCategory = 10)
dotplot(bp_enrichment, showCategory = 10)

# KEGG pathway enrichment if gene IDs are available
# Would require converting gene symbols to Entrez IDs first
gene_entrez <- bitr(real_genes, fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

kegg_enrichment <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = 'hsa',  # human
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Visualize top 10 KEGG pathways
if (nrow(kegg_enrichment) > 0) {
  barplot(kegg_enrichment, showCategory = 10)
}

#===============================================================================
# PART 6: INTEGRATING WGCNA AND BAYESIAN NETWORK RESULTS
#===============================================================================

# Identify hub genes from WGCNA
wgcna_hub_genes <- head(module_gene_info, 10)$gene_name

# Identify key regulators from Bayesian Network
bn_regulators <- rownames(head(key_regulators, 10))

# Find overlapping genes (important in both analyses)
overlapping_genes <- intersect(wgcna_hub_genes, bn_regulators)

cat("Genes identified as important in both WGCNA and Bayesian Network analysis:\n")
print(overlapping_genes)

# Prioritize these genes for further investigation
cat("\nThese genes are highly connected hub genes in the co-expression network")
cat(" AND likely have regulatory roles according to the Bayesian Network model.\n")

# Summarize biological processes associated with these key genes
if (length(overlapping_genes) > 0) {
  # Convert to real gene symbols for demonstration (in real analysis, use actual overlapping genes)
  sample_real_genes <- sample(keys(org.Hs.eg.db, keytype = "SYMBOL"), length(overlapping_genes))
  
  key_gene_enrichment <- perform_GO_enrichment(sample_real_genes)
  
  cat("\nBiological processes associated with key regulatory hub genes:\n")
  print(head(key_gene_enrichment, 5))
}

#===============================================================================
# FINAL SUMMARY
#===============================================================================

cat("\n==== ANALYSIS SUMMARY ====\n")
cat("1. WGCNA identified", length(unique(module_colors)), "co-expression modules\n")
cat("2. Module", gsub("ME", "", module_color), "showed strongest association with treatment\n")
cat("3. Bayesian Network analysis identified potential regulatory relationships among top module genes\n")
cat("4. Top regulatory genes to investigate further:", paste(rownames(head(key_regulators, 5)), collapse=", "), "\n")
cat("5. Hub genes from WGCNA:", paste(head(wgcna_hub_genes, 5), collapse=", "), "\n")
cat("6. Enrichment analysis revealed pathways associated with the key module\n")

# In a real analysis:
# 1. You would use your actual normalized expression data
# 2. Include more comprehensive sample traits
# 3. Carefully choose parameters based on your data characteristics
# 4. Validate findings through literature or experimental validation
# 5. Consider integrating with other omics data types
# 6. For very large datasets, you might need to focus Bayesian Network analysis on smaller subsets of genes

# End of workflow
