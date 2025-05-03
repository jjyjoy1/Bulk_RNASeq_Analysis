# PyWGCNA Network Analysis Toolkit

A comprehensive Python implementation of Weighted Gene Co-expression Network Analysis (WGCNA) with additional network analysis tools.

## Overview

This toolkit provides PyWGCNA features equivalent to R WGCNA, along with additional Python tools for comprehensive network analysis. The implementation focuses on maintaining the core functionality of the original R package while leveraging Python's ecosystem for improved performance and integration capabilities.

## Features

### PyWGCNA Core Features (R WGCNA Equivalents)

| Feature | Status | Description |
|---------|--------|-------------|
| Network construction | ✅ | Same algorithm as R implementation |
| Module detection | ✅ | Uses Python clustering libraries (similar results) |
| Soft thresholding | ✅ | Available with automatic power selection |
| Outlier detection | ✅ | Sample clustering and visualization |
| Module-trait relationships | ✅ | Correlation analysis |

### Key Differences from R WGCNA

#### 1. Visualization
- **PyWGCNA**: Relies on matplotlib/seaborn for visualization
- **R WGCNA**: Has more specialized plots
- Trade-off: Python offers greater flexibility for custom visualizations

#### 2. Integration
- **PyWGCNA**: 
  - Better integration with scikit-learn and pandas
  - Seamless data handling with NumPy and DataFrame structures
- **R WGCNA**: 
  - Deeper integration with Bioconductor packages
  - Better for R-centric bioinformatics workflows

#### 3. Performance
- **PyWGCNA**: Generally faster for large datasets
- **R WGCNA**: More optimization for specific biological use cases

#### 4. Documentation
- **R WGCNA**: Extensive tutorials and larger community examples
- **PyWGCNA**: Growing documentation (contributions welcome!)

## Additional Network Analysis Tools

### 1. PGMPY - Bayesian Network Analysis
Equivalent to R's `bnlearn`, PGMPY provides:
- Bayesian Network structure learning
- Parameter learning
- Inference algorithms
- Model evaluation

### 2. GProfiler - Functional Enrichment
Alternative to R's `clusterProfiler`, offering:
- Gene Ontology (GO) enrichment
- Pathway analysis
- Enrichment visualization
- Multiple organism support

### 3. NetworkX - Network Analysis & Visualization
Comprehensive network analysis toolkit providing:
- Network creation and manipulation
- Graph algorithms
- Centrality measures
- Community detection
- Network visualization

## Installation

```bash
# Core installation
pip install pywgcna

# Additional tools
pip install pgmpy gprofiler-official networkx matplotlib seaborn pandas numpy scikit-learn
```

## Quick Start

```python
import pywgcna

# Initialize WGCNA object
wgcna = pywgcna.WGCNA()

# Load your expression data
data = pywgcna.load_data('expression_data.csv')

# Perform network construction
wgcna.construct_network(data, power=12)

# Detect modules
modules = wgcna.detect_modules()

# Analyze module-trait relationships
trait_corr = wgcna.module_trait_correlation(traits)
```

## Use Cases

This toolkit is ideal for:
- Gene co-expression network analysis
- Systems biology research
- Biological pathway analysis
- Complex network analysis in bioinformatics
- Integration of Python-based machine learning with network biology


## Acknowledgments

- Original R WGCNA developers for the foundational algorithm
- Python scientific computing community
- All contributors to this project

## Roadmap

- [ ] Enhanced visualization tools
- [ ] Additional clustering algorithms
- [ ] GPU acceleration for large datasets
- [ ] Integration with deep learning frameworks
- [ ] Extended documentation and tutorials

## Other Complementary Python Tools

These additional tools could enhance the analysis:

- PhenoGraph: For cell type clustering in single-cell data
- scanpy: Comprehensive single-cell analysis toolkit
- PyComplexHeatmap: Advanced heatmap visualizations
- DiffBind: RNA-seq and ChIP-seq integration
- rPy2: Bridge to use R packages directly in Python
- pySurvival: For survival analysis with expression data
- XGBoost/LightGBM: Machine learning for predicting outcomes from modules


