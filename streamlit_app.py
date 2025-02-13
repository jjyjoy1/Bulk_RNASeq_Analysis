import streamlit as st
import pandas as pd
import numpy as np
import shap
import joblib
import matplotlib.pyplot as plt
import seaborn as sns

# Load Data
@st.cache_data
def load_data():
    counts = pd.read_csv("combined_counts.txt", sep="\t", index_col=0)
    metadata = pd.read_csv("sample_metadata.txt", sep="\t", index_col=0)
    top_genes = pd.read_csv("top_xai_genes.txt")
    enrichment = pd.read_csv("gene_enrichment_results.txt", sep="\t")
    return counts, metadata, top_genes, enrichment

counts, metadata, top_genes, enrichment = load_data()

# Sidebar
st.sidebar.header("RNA-Seq ML Analysis")
page = st.sidebar.radio("Select a Page", ["Dimensionality Reduction", "Feature Importance (XAI)", "Gene Enrichment", "Model Performance"])

# Dimensionality Reduction
if page == "Dimensionality Reduction":
    st.title("Dimensionality Reduction Visualization")
    st.image("dimension_reduction_plots.png", caption="PCA, t-SNE, and UMAP Visualizations")

# Feature Importance
elif page == "Feature Importance (XAI)":
    st.title("Top Genes Driving Predictions (XAI)")
    st.image("top_xai_genes.png", caption="SHAP-based Feature Importance")

    st.subheader("Top 30 Genes Identified by XAI")
    st.write(top_genes)

# Gene Enrichment Analysis
elif page == "Gene Enrichment":
    st.title("Gene Enrichment & Pathway Analysis")
    st.subheader("Top Enriched Pathways")
    st.write(enrichment.head(20))  # Show top 20 pathways

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(y=enrichment['name'][:10], x=enrichment['p.value'][:10], palette="coolwarm", ax=ax)
    ax.set_title("Top 10 Enriched Pathways")
    ax.set_xlabel("p-value")
    ax.set_ylabel("Pathway")
    st.pyplot(fig)

# Model Performance
elif page == "Model Performance":
    st.title("Model Performance & Predictions")

    # Load Trained Models
    rf_model = joblib.load("rf_best_model.pkl")
    xgb_model = joblib.load("xgb_best_model.pkl")

    st.subheader("Model Accuracy")
    st.write("Random Forest Accuracy: 90% (example)")
    st.write("XGBoost Accuracy: 92% (example)")

    st.subheader("Upload a New Sample for Prediction")
    uploaded_file = st.file_uploader("Upload a CSV file with gene expression values", type=["csv"])

    if uploaded_file is not None:
        input_data = pd.read_csv(uploaded_file, index_col=0)
        input_scaled = (input_data - counts.mean()) / counts.std()

        rf_pred = rf_model.predict(input_scaled.T)
        xgb_pred = xgb_model.predict(input_scaled.T)

        st.write(f"Random Forest Prediction: {'Test' if rf_pred[0] else 'Control'}")
        st.write(f"XGBoost Prediction: {'Test' if xgb_pred[0] else 'Control'}")

