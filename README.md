```markdown
# Comparison of RNA-seq Differential Expression Methods
# The new code contains comparision of Enformer and Borzoi for Genomic Sequence Modeling

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
```

*Comparision of Enformer with Borzoi:
**Both Enformer and Borzoi take DNA sequences as input but are designed for distinct biological scenarios and tasks. Here’s a breakdown of their input designs and **key limitations**:

---

### **1. DNA Input Scenarios**
| Model       | Input DNA Sequence Length | Key Scenario/Purpose                                                                 |
|-------------|---------------------------|-------------------------------------------------------------------------------------|
| **Enformer**| ~200 kb (196,608 bp)      | Predicts *chromatin states* (e.g., TF binding, accessibility) and *gene expression*. |
| **Borzoi**  | ~524 kb (524,288 bp)      | Predicts *RNA-seq coverage* (transcription, splicing, polyadenylation) at base-pair resolution. |

---

### **2. Limitations of Each Model**
#### **Enformer**
1. **Resolution Constraints**:
   - Outputs predictions at **128 bp resolution**, limiting fine-grained analysis (e.g., precise exon-intron boundaries or splicing events).
2. **RNA Complexity Blindness**:
   - Cannot model RNA-level processes like splicing, polyadenylation, or isoform diversity.
3. **Non-Coding Variant Interpretation**:
   - Struggles to predict effects of variants far from promoters or in unannotated regulatory regions.
4. **Fixed Input Size**:
   - Requires sequences of **exactly 196 kb**, necessitating padding/trimming for shorter/longer regions.
5. **Data Bias**:
   - Trained primarily on bulk ENCODE/FANTOM5 data, limiting accuracy in rare cell types or tissues.

#### **Borzoi**
1. **Computational Cost**:
   - Training from scratch requires **terabytes of data** and weeks on GPU clusters (mitigated partially by Flashzoi).
2. **Splicing Limitations**:
   - Struggles with unannotated splice junctions or non-canonical splicing signals.
3. **Input Context Dependency**:
   - Predictions depend on the full 524 kb input window; truncating sequences may miss distal regulatory elements.
4. **Overfitting Risks**:
   - High model capacity may lead to overfitting on small custom datasets during fine-tuning.
5. **Tissue-Specific Gaps**:
   - Performance varies in tissues with sparse training data (e.g., brain subtypes, rare cancers).

---

### **3. Shared Limitations**
- **Sequence-Centric Focus**:
  Both models ignore 3D chromatin structure, RNA-protein interactions, or epigenetic modifications (e.g., methylation).
- **Static Predictions**:
  Predictions reflect population-level averages, not dynamic cell-state transitions (e.g., differentiation or stress responses).
- **Non-Human Applications**:
  Limited generalizability to non-model organisms without retraining.

---

### **4. When to Use Which Model?**
| Use Case                                  | Preferred Model | Reason                                                                 |
|-------------------------------------------|-----------------|-----------------------------------------------------------------------|
| Enhancer-promoter interactions           | Enformer        | Optimized for long-range chromatin interactions and TF binding.       |
| Splicing or isoform-level predictions    | Borzoi          | Higher resolution (32 bp) and RNA-seq coverage modeling.              |
| Variant effect scoring in non-coding DNA | Enformer        | Better for linking variants to chromatin states and gene expression. |
| Tissue-specific RNA expression profiling | Borzoi          | Integrates RNA-seq data across tissues with cell-type specificity.    |

---

While both models use DNA sequences as input, their limitations stem from **architectural constraints** (resolution, input size) and **training data biases**. Choose Enformer for chromatin/gene regulation tasks and Borzoi for RNA-level predictions, but validate results with orthogonal assays (e.g., CRISPR screens) to address model-specific blind spots.


