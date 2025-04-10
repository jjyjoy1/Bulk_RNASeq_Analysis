Here's the Multi-Sample Analysis with Enformer document converted to well-structured GitHub Markdown:

```markdown
# Multi-Sample Analysis with Enformer for Biomarker Extraction

This advanced pipeline demonstrates how to analyze multiple genetic samples (from VCF files) using Enformer to identify regulatory patterns as biomarkers for sample classification.

## Purpose
This approach is particularly useful for:

- **Disease subtyping**: Finding regulatory signatures that distinguish disease subtypes
- **Drug response prediction**: Identifying regulatory patterns associated with treatment response
- **Personalized medicine**: Classifying patients based on predicted gene regulation patterns

## Analysis Pipeline

### 1. Input Preparation and Processing
Processes a VCF file containing genetic variants from multiple samples:

```python
vcf_data = self.process_vcf_region(vcf_path, region, sample_limit)
```

For each sample:
- Extracts the reference DNA sequence for a genomic region
- Applies the sample's specific variants to create a personalized sequence

### 2. Enformer Predictions
Runs Enformer on both reference and sample sequences:

```python
ref_predictions = self.predict_sequence(ref_sequence)

for sample, variants in variants_by_sample.items():
    sample_sequence = self.extract_sequence_with_variants(enformer_interval, variants)
    sample_predictions[sample] = self.predict_sequence(sample_sequence)
```

Each prediction includes thousands of gene regulatory tracks across different tissues and cell types.

### 3. Comparative Analysis
Computes regulatory differences between samples and reference:

```python
diff_predictions[sample] = sample_predictions[sample] - ref_predictions
```

Builds a feature matrix capturing maximum regulatory impact:

```python
feature_matrix = np.zeros((len(variants_by_sample), len(tracks_to_analyze)))
for i, sample in enumerate(variants_by_sample.keys()):
    for j, track_idx in enumerate(tracks_to_analyze):
        feature_matrix[i, j] = np.abs(diff_predictions[sample][:, track_idx]).max()
```

### 4. Biomarker Identification
Identifies significant regulatory changes:

```python
significant_tracks = self._identify_significant_tracks(feature_matrix, tracks_to_analyze)
```

Calculates biomarker scores combining significance and cluster separation:

```python
track['biomarker_score'] = -np.log10(track['p_value']) * track['cluster_separation']
```

### 5. Sample Clustering
Clusters samples using dimensionality reduction and K-means:

```python
# PCA for visualization
pca = PCA(n_components=2)
pca_result = pca.fit_transform(normalized_features)

# K-means clustering
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
cluster_labels = kmeans.fit_predict(normalized_features)
```

## Visualization and Interpretation

The pipeline generates several key visualizations:

- **Sample clustering plots**: Shows sample grouping based on regulatory profiles
- **Biomarker importance plots**: Highlights most significant regulatory elements
- **Regulatory profiles**: Displays predicted gene regulation patterns
- **Individual sample effects**: Visualizes variant effects on regulation per sample

## Usage Example

```python
analyzer = EnformerBiomarkerAnalyzer(
    vcf_path="samples.vcf",
    reference_genome="hg38",
    region="chr1:1000000-2000000"
)
results = analyzer.run_analysis()
analyzer.visualize_results()
```

## Dependencies
- Enformer model
- NumPy
- Scikit-learn
- Matplotlib/Seaborn
- PyVCF
```




