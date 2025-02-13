import pandas as pd
import numpy as np
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr
from pyvis.network import Network

# Load Expression Data (Top Selected Genes)
df = pd.read_csv("selected_genes.txt", header=None)
top_genes = df[0].tolist()
expression_data = pd.read_csv("combined_counts.txt", sep="\t", index_col=0)
expression_data = expression_data.loc[top_genes].T  # Transpose for genes as columns

# Compute Pairwise Correlations (Spearman)
cor_matrix = expression_data.corr(method="spearman")
cor_matrix.to_csv("gene_correlation_matrix.txt", sep="\t")

# Heatmap of Correlations
plt.figure(figsize=(12, 8))
sns.heatmap(cor_matrix, cmap="coolwarm", annot=False)
plt.title("Gene Expression Correlation Matrix")
plt.savefig("gene_correlation_heatmap.png")
plt.show()

### GRN Inference using GENIE3 ###
def genie3_grn(data, n_trees=1000):
    """Infer a Gene Regulatory Network using Random Forest"""
    num_genes = data.shape[1]
    gene_names = data.columns
    feature_importance_matrix = np.zeros((num_genes, num_genes))

    for target_idx in range(num_genes):
        target_gene = gene_names[target_idx]
        X = np.delete(data.values, target_idx, axis=1)  # Remove target gene
        y = data[target_gene]

        model = RandomForestRegressor(n_estimators=n_trees, random_state=42)
        model.fit(X, y)

        feature_importance_matrix[target_idx, np.arange(num_genes) != target_idx] = model.feature_importances_

    return pd.DataFrame(feature_importance_matrix, index=gene_names, columns=gene_names)

grn_matrix = genie3_grn(expression_data)
grn_matrix.to_csv("gene_regulatory_network_matrix.txt", sep="\t")

# Construct Graph
G = nx.Graph()
threshold = np.percentile(grn_matrix.values, 95)  # Keep top 5% connections

for gene1 in grn_matrix.columns:
    for gene2 in grn_matrix.columns:
        if gene1 != gene2 and grn_matrix.loc[gene1, gene2] > threshold:
            G.add_edge(gene1, gene2, weight=grn_matrix.loc[gene1, gene2])

# Save Network Visualization
nt = Network(height="600px", width="100%", notebook=True)
nt.from_nx(G)
nt.show("gene_regulatory_network.html")



