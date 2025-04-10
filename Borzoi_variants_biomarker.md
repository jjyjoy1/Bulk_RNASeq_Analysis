Here's the Borzoi analysis documentation converted to GitHub Markdown format:

```markdown
# Multi-Sample Analysis with Borzoi

This document explains how to use the Borzoi model for analyzing multiple genomic samples, highlighting key differences from Enformer-based approaches.

## Key Differences Between Enformer and Borzoi

| Feature               | Borzoi                          | Enformer                      |
|-----------------------|---------------------------------|-------------------------------|
| **Primary Focus**     | RNA-seq coverage prediction     | Gene regulation prediction    |
| **Resolution**        | 32bp                           | 128bp                        |
| **Input Sequence**    | 524,288bp                      | 393,216bp                    |
| **Model Architecture**| Transformer + U-net connections | Transformer                  |
| **Regulatory Layers** | Transcription, splicing, polyadenylation | Gene expression       |

## Analysis Pipeline

The Borzoi multi-sample analysis implements:

1. **Sequence Extraction**
   - Reference sequences
   - Variant-containing sample sequences

2. **Borzoi Predictions**
   - RNA-seq coverage predictions for all sequences

3. **Comparative Analysis**
   - Identifies regulatory pattern differences

4. **Biomarker Identification**
   - Statistical analysis of differences across samples

5. **Sample Clustering**
   - Groups samples based on regulatory profiles

6. **Visualization**
   - Biomarker plots
   - Sample clustering diagrams
   - Individual effect visualizations

## Example Usage

### With Real Data
```bash
python borzoi_analysis.py \
  --vcf samples.vcf \
  --region chr7:55000000-55100000 \
  --genome hg38.fa \
  --model_dir ~/borzoi/models/human
```

### With Mock Data (Demonstration)
```bash
python borzoi_analysis.py --demo
```

## Practical Applications

- **Disease Subtyping**: Identify regulatory signatures distinguishing disease subtypes
- **Drug Response**: Find regulatory patterns associated with treatment response
- **RNA Splicing**: Predict variant effects on RNA processing
- **Polyadenylation**: Understand impacts on RNA stability and 3' end processing

## Installation Requirements

### 1. Clone and Install Packages
```bash
# Baskerville (Borzoi dependency)
git clone https://github.com/calico/baskerville.git
cd baskerville
pip install -e .

# Borzoi
git clone https://github.com/calico/borzoi.git
cd borzoi
pip install -e .
```

### 2. Download Model Weights
```bash
cd borzoi
./download_models.sh
```

### 3. Install Additional Dependencies
```bash
pip install tensorflow pandas scikit-learn matplotlib seaborn pyfaidx pyvcf
```

## Configuration Options

Key parameters in `borzoi_analysis.py`:

```python
DEFAULT_PARAMS = {
    'batch_size': 4,          # Prediction batch size
    'min_abs_diff': 0.5,      # Minimum absolute difference threshold
    'p_value_cutoff': 0.01,   # Significance threshold for biomarkers
    'n_clusters': 3,          # Default number of sample clusters
    'tracks_to_analyze': [...] # Specific regulatory tracks to analyze
}
```

## Output Files

The pipeline generates:
- `results/biomarkers.csv`: Significant regulatory biomarkers
- `results/cluster_assignments.tsv`: Sample cluster memberships
- `plots/`: Directory containing all visualizations
```

Would you like me to:
1. Add a troubleshooting section for common installation issues?
2. Include example output images with captions?
3. Add more detailed configuration options?
4. Include a sample VCF file format example?

The Markdown includes a comparison table, clear code blocks, and organized sections for easy navigation. I can adjust the level of technical detail or add more practical examples as needed.


