import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
from sklearn.feature_selection import SelectKBest, f_classif, VarianceThreshold
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout

# Load count data
data = pd.read_csv("combined_counts.txt", sep="\t", index_col=0)
metadata = pd.read_csv("sample_metadata.txt", sep="\t", index_col=0)

# Encode labels (control=0, test=1)
le = LabelEncoder()
metadata['condition'] = le.fit_transform(metadata['condition'])

# Ensure correct sample order
X = data.T  # Transpose so samples are rows
y = metadata.loc[X.index, "condition"]

# Standardize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

### Feature Selection ###
# Remove low variance genes
var_thresh = VarianceThreshold(threshold=0.01)
X_var = var_thresh.fit_transform(X_scaled)

# Select top 500 genes using ANOVA F-test
selector = SelectKBest(score_func=f_classif, k=500)
X_selected = selector.fit_transform(X_var, y)
selected_genes = X.columns[selector.get_support()]

# Save selected genes
pd.DataFrame(selected_genes).to_csv("selected_genes.txt", index=False, header=["Gene"])

### Dimensionality Reduction ###
# PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_selected)

# t-SNE
X_tsne = TSNE(n_components=2, perplexity=5, random_state=42).fit_transform(X_selected)

# UMAP
X_umap = umap.UMAP(n_components=2).fit_transform(X_selected)

# Plot results
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
axes[0].scatter(X_pca[:, 0], X_pca[:, 1], c=y, cmap="coolwarm")
axes[0].set_title("PCA")

axes[1].scatter(X_tsne[:, 0], X_tsne[:, 1], c=y, cmap="coolwarm")
axes[1].set_title("t-SNE")

axes[2].scatter(X_umap[:, 0], X_umap[:, 1], c=y, cmap="coolwarm")
axes[2].set_title("UMAP")

plt.savefig("dimension_reduction_plots.png")
plt.show()

### ML Model Training ###
X_train, X_test, y_train, y_test = train_test_split(X_selected, y, test_size=0.2, random_state=42)

# Random Forest
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)
y_pred_rf = rf.predict(X_test)
print("Random Forest Accuracy:", accuracy_score(y_test, y_pred_rf))
print(classification_report(y_test, y_pred_rf))

# XGBoost
xgb = XGBClassifier(use_label_encoder=False, eval_metric='logloss')
xgb.fit(X_train, y_train)
y_pred_xgb = xgb.predict(X_test)
print("XGBoost Accuracy:", accuracy_score(y_test, y_pred_xgb))
print(classification_report(y_test, y_pred_xgb))

### Neural Network ###
model = Sequential([
    Dense(256, activation='relu', input_shape=(X_selected.shape[1],)),
    Dropout(0.3),
    Dense(128, activation='relu'),
    Dropout(0.3),
    Dense(1, activation='sigmoid')
])

model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.fit(X_train, y_train, epochs=20, batch_size=16, validation_data=(X_test, y_test))

# Evaluate Neural Network
y_pred_nn = (model.predict(X_test) > 0.5).astype("int32")
print("Neural Network Accuracy:", accuracy_score(y_test, y_pred_nn))

### Explainable AI (XAI) ###
explainer = shap.TreeExplainer(xgb)
shap_values = explainer.shap_values(X_test)

# Plot SHAP summary
shap.summary_plot(shap_values, features=X_test, feature_names=selected_genes)

# Identify top genes driving predictions
shap_importance = np.abs(shap_values).mean(axis=0)
top_shap_genes = pd.DataFrame({'Gene': selected_genes, 'SHAP_importance': shap_importance})
top_shap_genes = top_shap_genes.sort_values(by="SHAP_importance", ascending=False).head(30)

# Save top genes
top_shap_genes.to_csv("top_xai_genes.txt", index=False)

# Bar plot of top 30 genes
plt.figure(figsize=(10, 6))
sns.barplot(x="SHAP_importance", y="Gene", data=top_shap_genes, palette="coolwarm")
plt.title("Top 30 Genes Driving Predictions (XAI)")
plt.savefig("top_xai_genes.png")
plt.show()



