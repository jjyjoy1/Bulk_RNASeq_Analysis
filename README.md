```markdown
# Comparison of RNA-seq Differential Expression Methods

Each method (DESeq2, edgeR, and limma/voom) approaches RNA-seq count data analysis with different statistical frameworks, particularly in modeling, normalization, and hypothesis testing.

## Statistical Models Overview

### DESeq2
**Statistical Model**:
- Uses negative binomial (NB) distribution
- Models gene counts with condition-dependent means and dispersion parameters

**Key Features**:
- **Normalization**: Median-of-ratios method
- **Dispersion Estimation**:
  - Gene-wise dispersions
  - Empirical Bayes shrinkage for stabilization
- **Testing**:
  - Wald test (or LRT for complex designs)
  - Log₂ fold change shrinkage

### edgeR
**Statistical Model**:
- Negative binomial distribution
- Flexible dispersion estimation options

**Key Features**:
- **Normalization**: TMM (Trimmed Mean of M-values)
- **Dispersion Options**:
  - Common (shared across genes)
  - Trended (function of mean expression)
  - Tagwise (gene-specific)
- **Testing**:
  - Exact tests (two-group comparisons)
  - GLM-based tests (complex designs)
  - Quasi-likelihood F-tests

### limma/voom
**Original Framework**:
- Developed for microarrays
- Linear models with empirical Bayes moderation

**voom Transformation**:
1. Converts counts to logCPM
2. Models mean-variance relationship
3. Assigns precision weights

**Testing**:
- Moderated t-tests
- Empirical Bayes variance shrinkage

## Key Differences

| Feature               | DESeq2               | edgeR                | limma/voom           |
|-----------------------|----------------------|----------------------|----------------------|
| **Distribution**      | Negative Binomial    | Negative Binomial    | Linear Models       |
| **Normalization**     | Median-of-ratios     | TMM                  | TMM + voom          |
| **Dispersion**        | Gene-wise + shrinkage| Multiple options     | Mean-variance trend  |
| **Testing**           | Wald/LRT             | Exact/GLM tests      | Moderated t-tests   |
| **Best For**          | Small samples        | Flexible designs     | Large samples       |

## When to Use Each Method

### DESeq2
- Preferred for:
  - Small sample sizes
  - Experiments with low counts
  - When fold change shrinkage is beneficial

### edgeR
- Preferred for:
  - Flexible experimental designs
  - When different dispersion estimation methods are needed
  - Exact tests for simple comparisons

### limma/voom
- Preferred for:
  - Larger sample sizes (>10 per group)
  - Complex experimental designs
  - When leveraging limma's mature linear modeling framework

## Code Example Comparison

```r
# DESeq2
dds <- DESeqDataSetFromMatrix(countData, colData, design = ~condition)
dds <- DESeq(dds)
res <- results(dds)

# edgeR
y <- DGEList(counts)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)

# limma/voom
v <- voom(counts, design)
fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit)
```

## Practical Considerations

1. **Sample Size**:
   - <5 samples/group: DESeq2 or edgeR
   - >10 samples/group: limma/voom often performs well

2. **Computational Resources**:
   - DESeq2 is generally more memory intensive
   - limma/voom can be faster for large datasets

3. **Complex Designs**:
   - All support complex designs, but edgeR and limma offer particularly flexible GLM frameworks

---

# Comparison of Enformer and Borzoi for Genomic Sequence Modeling

Both models take DNA sequences as input but are designed for distinct biological scenarios.

## Input Specifications

| Model       | Input Length | Resolution | Primary Outputs |
|-------------|-------------|------------|-----------------|
| Enformer    | 196,608 bp  | 128 bp     | Chromatin states, gene expression |
| Borzoi      | 524,288 bp  | 32 bp      | RNA-seq coverage (transcription, splicing) |

## Key Limitations

### Enformer
1. **Resolution Constraints**:
   - 128 bp output resolution limits fine-grained analysis
2. **RNA Processes**:
   - Cannot model splicing or polyadenylation
3. **Non-Coding Variants**:
   - Struggles with distal regulatory elements
4. **Fixed Input**:
   - Requires exact 196 kb sequences
5. **Data Bias**:
   - Trained primarily on bulk ENCODE/FANTOM5 data

### Borzoi
1. **Computational Cost**:
   - Requires extensive resources for training
2. **Splicing Limitations**:
   - Poor performance on non-canonical splicing
3. **Context Dependency**:
   - Needs full 524 kb window for accurate predictions
4. **Overfitting Risk**:
   - High model capacity may overfit small datasets
5. **Tissue Specificity**:
   - Variable performance in rare cell types

## Shared Limitations
- Ignores 3D chromatin structure
- Static predictions (no dynamics)
- Limited non-human applications

## Model Selection Guide

| Use Case                                  | Preferred Model |
|-------------------------------------------|-----------------|
| Enhancer-promoter interactions           | Enformer        |
| Splicing/isoform predictions            | Borzoi          |
| Non-coding variant effects              | Enformer        |
| Tissue-specific RNA profiling           | Borzoi          |

**Note**: Always validate predictions with orthogonal assays to address model-specific limitations.
```

Key improvements made:
1. Removed duplicate heading and merged sections logically
2. Standardized all tables with consistent formatting
3. Added clear horizontal rule (`---`) to separate the two main comparisons
4. Unified the style of limitation lists between models
5. Improved the flow from general to specific information
6. Ensured all headers follow proper hierarchy
7. Added a note about validation in the conclusion

The document now has better visual separation between topics while maintaining a cohesive structure. Each comparison stands on its own while being clearly part of the larger document.

