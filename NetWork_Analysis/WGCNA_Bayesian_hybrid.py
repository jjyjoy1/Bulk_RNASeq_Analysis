# End-to-End RNAseq Analysis with PyWGCNA and Python Tools
# This workflow uses Python alternatives to WGCNA and includes similar functionalities

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pyWGCNA import WGCNA
import networkx as nx
from pgmpy.models import BayesianNetwork
from pgmpy.estimators import HillClimbSearch, BicScore, MaximumLikelihoodEstimator
from pgmpy.inference import VariableElimination
from sklearn.preprocessing import quantile_transform
import gprofiler
import warnings
warnings.filterwarnings('ignore')

#===============================================================================
# PART 1: SETUP AND DATA GENERATION (In practice, you'd load your actual data)
#===============================================================================

# Create a dummy dataset for demonstration
# Replace this with your actual normalized expression matrix
np.random.seed(42)
n_genes = 3000  # Number of genes
n_samples = 30  # Number of samples

# Create gene names
gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]

# Create sample names
sample_names = [f"SAMPLE_{i:03d}" for i in range(n_samples)]

# Create dummy expression data (normally distributed)
expr_data = np.random.randn(n_genes, n_samples)

# Create DataFrame (samples as columns, genes as rows)
expr_df = pd.DataFrame(expr_data, index=gene_names, columns=sample_names)

# Create sample information
sample_info = pd.DataFrame({
    'treatment': [0] * (n_samples//2) + [1] * (n_samples//2),
    'age': np.random.randint(20, 80, n_samples),
    'sex': np.random.choice(['M', 'F'], n_samples)
}, index=sample_names)

#===============================================================================
# PART 2: PYWGCNA ANALYSIS
#===============================================================================

# Initialize PyWGCNA
pyWGCNA = WGCNA(
    name='RNA_seq_analysis',
    species='homo sapiens',  # or your species of choice
    geneExpPath=None  # We'll use the DataFrame directly
)

# Prepare data for PyWGCNA (samples as rows, genes as columns)
pyWGCNA.geneExpr = expr_df.T  # Transpose to have samples as rows

# Set gene list
pyWGCNA.gene_list = list(pyWGCNA.geneExpr.columns)

# Pre-filtering (optional)
pyWGCNA.prefilterGenes()

# Sample clustering to detect outliers
pyWGCNA.detectOutliers(doPlot=True)

# Choose soft-thresholding power
try:
    pyWGCNA.findSoftPower(corFnc='bicor', corOptions={'power': 2, 'cutoff': 0.8}, 
                         powerVector=list(range(1, 21)), 
                         networkType='signed', 
                         RsquaredCut=0.8)
except:
    # If automatic detection fails, set a default power
    pyWGCNA.power = 6

# Construct the network and identify modules
pyWGCNA.buildWGCNAnet(power=pyWGCNA.power, 
                      TOMType='signed', 
                      networkType='signed',
                      minModuleSize=30)

# Get module colors and preservation
module_colors = pyWGCNA.datME.index

# Plot the dendrogram
pyWGCNA.plotDendro(plot_gene_dendrogram=True)

# Get module-trait relationships
# For this, we'll correlate module eigengenes with traits
module_eigengenes = pyWGCNA.datME[module_colors].values
trait_data = sample_info[['treatment']].values

# Calculate correlation between module eigengenes and traits
module_trait_corr = pd.DataFrame(
    index=module_colors,
    columns=['treatment']
)

for i, module in enumerate(module_colors):
    corr, pval = stats.pearsonr(module_eigengenes[:, i], trait_data[:, 0])
    module_trait_corr.loc[module, 'treatment'] = corr

# Plot module-trait relationships
plt.figure(figsize=(10, 8))
sns.heatmap(module_trait_corr.astype(float), 
            cmap='RdBu_r', 
            center=0, 
            annot=True)
plt.title('Module-Trait Relationships')
plt.xlabel('Traits')
plt.ylabel('Modules')
plt.tight_layout()
plt.show()

# Get genes in modules
genes_in_modules = pyWGCNA.datModules

# Identify hub genes in modules
# Find the module with highest correlation to treatment
max_corr_module = module_trait_corr['treatment'].idxmax()

# Get genes in the module of interest
module_genes = genes_in_modules[genes_in_modules['moduleColors'] == max_corr_module]

# Calculate gene significance and module membership
gene_significance = pd.Series(index=gene_names)
module_membership = pd.Series(index=gene_names)

for gene in gene_names:
    # Gene significance (correlation with trait)
    gene_expr = expr_df.loc[gene, :]
    corr, _ = stats.pearsonr(gene_expr, sample_info['treatment'])
    gene_significance[gene] = abs(corr)
    
    # Module membership (if in the module)
    if gene in module_genes.index:
        module_idx = module_colors.get_loc(max_corr_module)
        corr, _ = stats.pearsonr(gene_expr, module_eigengenes[:, module_idx])
        module_membership[gene] = abs(corr)

# Get top hub genes (high module membership and gene significance)
hub_genes = module_genes[module_genes['moduleColors'] == max_corr_module].copy()
hub_genes['gene_significance'] = gene_significance[hub_genes.index]
hub_genes['module_membership'] = module_membership[hub_genes.index]
hub_genes = hub_genes.sort_values(['module_membership', 'gene_significance'], ascending=[False, False])

# Top 20 hub genes
top_hub_genes = hub_genes.head(20)

#===============================================================================
# PART 3: BAYESIAN NETWORK ANALYSIS WITH PGMPY
#===============================================================================

# Select top genes from the module for Bayesian Network analysis
# Using top 30 genes to make computation manageable
top_module_genes = top_hub_genes.head(30).index

# Extract expression data for these genes
module_expr_data = expr_df.loc[top_module_genes, :].T

# Discretize the data for Bayesian Network
# Convert to categorical (low, medium, high)
bn_data = pd.DataFrame()
for gene in top_module_genes:
    gene_data = module_expr_data[gene]
    gene_data_disc = pd.qcut(gene_data, 3, labels=['low', 'medium', 'high'])
    bn_data[gene] = gene_data_disc

# Learn the Bayesian Network structure using Hill Climbing
hc = HillClimbSearch(bn_data)
bic = BicScore(bn_data)
best_model = hc.estimate(scoring_method=bic)

# Create the Bayesian Network
model = BayesianNetwork(best_model.edges())

# Learn the parameters (CPDs)
model.fit(bn_data, estimator=MaximumLikelihoodEstimator)

# Find key regulatory genes (nodes with many outgoing edges)
out_degrees = dict(model.out_degree())
key_regulators = sorted(out_degrees.items(), key=lambda x: x[1], reverse=True)

# Convert to NetworkX for visualization
G = nx.DiGraph()
G.add_edges_from(best_model.edges())

# Plot the Bayesian Network
plt.figure(figsize=(12, 10))
pos = nx.spring_layout(G, k=0.5, iterations=50)
nx.draw(G, pos, with_labels=True, node_size=1500, node_color='lightblue', 
        font_size=8, arrows=True, arrowsize=20)
plt.title('Bayesian Network of Top Module Genes')
plt.tight_layout()
plt.show()

# Predict effects of gene perturbations
# Set up inference
infer = VariableElimination(model)

# Simulate setting a key regulator to high expression
top_regulator = key_regulators[0][0]
evidence = {top_regulator: 'high'}

# Get predictions for other genes
predictions = {}
for gene in top_module_genes:
    if gene != top_regulator:
        result = infer.query(variables=[gene], evidence=evidence)
        predictions[gene] = result

# Display predictions
print("\nEffect of setting", top_regulator, "to high:")
for gene, result in predictions.items():
    print(f"{gene}: {result}")

#===============================================================================
# PART 4: ENRICHMENT ANALYSIS USING GPROFILER
#===============================================================================

# Get gene lists for enrichment analysis
# Convert module genes to appropriate format (gene symbols)
# Note: In real analysis, you'd use actual gene identifiers

# Use gprofiler for Gene Ontology and pathway enrichment
gp = gprofiler.GProfiler(return_dataframe=True)

# Perform enrichment for the module of interest
# Note: For this example, we're using dummy gene names
# In real analysis, you would use actual gene identifiers
try:
    enrichment_results = gp.profile(query=top_hub_genes.head(100).index.tolist(),
                                    organism='hsapiens',
                                    sources=['GO:BP', 'KEGG'],
                                    significant=True)
    
    # Display top enrichment results
    print("\nTop Enriched Terms:")
    print(enrichment_results[['name', 'p_value', 'source']].head(10))
    
    # Visualize enrichment results
    if not enrichment_results.empty:
        top_terms = enrichment_results.head(15)
        plt.figure(figsize=(10, 8))
        sns.barplot(data=top_terms, x='p_value', y='name')
        plt.xscale('log')
        plt.xlabel('-log10 p-value')
        plt.title('Top 15 Enriched Terms')
        plt.tight_layout()
        plt.show()

except Exception as e:
    print(f"Enrichment analysis failed with dummy data: {e}")
    print("In real analysis, you would use actual gene identifiers.")

#===============================================================================
# PART 5: INTEGRATING PYWGCNA AND BAYESIAN NETWORK RESULTS
#===============================================================================

# Find overlapping important genes
wgcna_hub_genes = top_hub_genes.head(10).index.tolist()
bn_regulators = [reg[0] for reg in key_regulators[:10]]

overlapping_genes = list(set(wgcna_hub_genes).intersection(set(bn_regulators)))

print("\n==== INTEGRATION SUMMARY ====")
print(f"Top WGCNA hub genes: {wgcna_hub_genes[:5]}")
print(f"Top Bayesian Network regulators: {bn_regulators[:5]}")
print(f"Overlapping key genes: {overlapping_genes}")

# Visualize gene importance from both methods
importance_df = pd.DataFrame({
    'WGCNA_membership': module_membership,
    'WGCNA_significance': gene_significance
})

for gene, score in out_degrees.items():
    if gene in importance_df.index:
        importance_df.loc[gene, 'BN_out_degree'] = score

# Create a scatter plot showing both WGCNA and BN importance
important_genes = list(set(wgcna_hub_genes + bn_regulators))
plt.figure(figsize=(10, 8))
for gene in important_genes:
    if gene in importance_df.index:
        plt.scatter(importance_df.loc[gene, 'WGCNA_membership'], 
                   importance_df.loc[gene, 'BN_out_degree'],
                   s=100, alpha=0.7)
        plt.annotate(gene, (importance_df.loc[gene, 'WGCNA_membership'], 
                           importance_df.loc[gene, 'BN_out_degree']),
                    fontsize=8)

plt.xlabel('WGCNA Module Membership')
plt.ylabel('Bayesian Network Out-degree')
plt.title('Gene Importance: WGCNA vs Bayesian Network')
plt.tight_layout()
plt.show()

#===============================================================================
# PART 6: SAVING RESULTS
#===============================================================================

# Save module assignments
genes_in_modules.to_csv('module_assignments.csv')

# Save hub genes
top_hub_genes.to_csv('top_hub_genes.csv')

# Save Bayesian Network structure
nx.write_edgelist(G, 'bayesian_network_edges.txt')

# Save module-trait correlations
module_trait_corr.to_csv('module_trait_correlations.csv')

print("\n==== ANALYSIS COMPLETE ====")
print("Results saved to files:")
print("- module_assignments.csv")
print("- top_hub_genes.csv")  
print("- bayesian_network_edges.txt")
print("- module_trait_correlations.csv")

# Additional comparison with R WGCNA results (if available)
print("\nKey Differences from R WGCNA:")
print("1. PyWGCNA uses slightly different algorithms for module detection")
print("2. Visualization options are more limited in Python")
print("3. Some advanced features may require manual implementation")
print("4. Integration with other Python tools is easier")


