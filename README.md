Each of these popular methods DESeq2, edgeR, and limma/voom approaches the analysis of RNA-seq count data with a slightly different statistical framework, particularly in how they model the count distribution, estimate variance (or dispersion), normalize data, and perform hypothesis testing. 

Below is a detailed explanation of each method and their main differences.

1. Modeling Count Data

DESeq2
Statistical Model:
DESeq2 models RNA-seq counts using a negative binomial (NB) distribution. In this model, each gene count is assumed to be NB-distributed with a mean that depends on the experimental condition and a dispersion parameter that captures extra-Poisson variability (i.e., biological variation).

Normalization:
It uses the median-of-ratios method to normalize counts, which corrects for differences in sequencing depth and RNA composition.

Dispersion Estimation:
DESeq2 estimates gene-wise dispersions and then uses empirical Bayes shrinkage to "borrow strength" across all genes. This shrinkage stabilizes the dispersion estimates, especially when sample sizes are small.
Testing and Shrinkage:
Differential expression is typically assessed with a Wald test (or a likelihood ratio test for more complex designs). Additionally, DESeq2 applies shrinkage to the log₂ fold changes, reducing the impact of low-count variability and providing more reliable effect size estimates.

edgeR
Statistical Model:
Like DESeq2, edgeR assumes that counts follow a negative binomial distribution. However, edgeR provides several ways to estimate dispersion:
Common Dispersion: A single dispersion value shared across all genes.
Trended Dispersion: A trend of dispersion values as a function of the genes mean expression.
Tagwise (Gene-specific) Dispersion: Individual dispersion estimates for each gene.
Normalization:
edgeR uses the TMM (Trimmed Mean of M-values) normalization method to adjust for differences in library sizes and composition biases.
Dispersion Estimation and Testing:
EdgeR employs empirical Bayes methods to moderate the dispersion estimates. For testing, it offers:
Exact Tests: For simple two-group comparisons (analogous to a Fishers exact test adapted for NB data).
Generalized Linear Models (GLMs): For more complex experimental designs, using either likelihood ratio tests or quasi-likelihood F-tests, which account for variability more robustly.

limma with voom
Original Framework (limma):
Originally developed for microarray data, limma uses linear models along with empirical Bayes moderation to improve variance estimates across genes.
The voom Transformation:
RNA-seq count data are inherently heteroscedastic (i.e., variance depends on the mean). The voom method addresses this by:
Transforming Counts: Converting raw counts into log counts per million (log CPM) values.
Estimating the Mean-Variance Relationship: It models how the variance changes with the mean expression level.
Assigning Weights: Each observation is given a precision weight based on its estimated variance. This weighting makes the data amenable to standard linear modeling.
Testing:
After the voom transformation, limma applies linear modeling and computes moderated t-statistics (using empirical Bayes shrinkage), which helps stabilize variance estimates, especially in studies with small sample sizes.

2. Key Differences and Considerations

Distributional Assumptions:
DESeq2 & edgeR:
Both explicitly model the count data using the negative binomial distribution. This is well suited for RNA-seq counts, particularly when the counts are low or when there is overdispersion relative to the Poisson model.
limma/voom:
Instead of modeling the counts directly, voom transforms the data so that the variability across genes becomes approximately constant (homoscedastic) on the log scale, thereby allowing the use of linear models.

Normalization Approaches:
DESeq2: Uses median-of-ratios.
edgeR: Uses TMM normalization.
limma/voom: Often uses TMM (or other similar methods) before applying the voom transformation to log CPM.

Dispersion/Variance Estimation:
DESeq2: Estimates gene-wise dispersions and shrinks these estimates using an empirical Bayes approach.
edgeR: Offers flexibility by providing common, trended, and tagwise dispersion estimates, and then uses empirical Bayes to stabilize these estimates.
limma/voom: Instead of directly estimating dispersion for count data, voom estimates the mean-variance trend and then assigns weights to each observation for use in linear models.

Testing Procedures:
DESeq2:
Uses the Wald test (or likelihood ratio tests) on the NB model parameters. It also shrinks log fold changes to improve interpretability, especially for low-count genes.
edgeR:
Uses exact tests for simple comparisons and GLM-based tests (likelihood ratio or quasi-likelihood tests) for complex designs.
limma/voom:
After transforming the data, it uses moderated t-tests in a linear modeling framework, which benefits from the empirical Bayes shrinkage of variance estimates.

When to Use Which:
DESeq2 and edgeR are often preferred for RNA-seq data, especially when dealing with low counts or small sample sizes, because their NB modeling is directly suited to count data.
limma/voom is particularly effective when you have a larger sample size or when you prefer the linear modeling framework. The voom transformation makes it possible to apply the powerful and flexible methods from limma, especially in complex experimental designs.


