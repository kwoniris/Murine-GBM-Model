# GBM Murine Model Data Analysis

Iris Kwon
July, 2023

#### Description

Repository for differential gene analysis on the data set of 4 murine GBM samples.

### Data

Current data is derived from ORlow and ORhigh syngeneic tumor models from two murine glioma cell lines, GL261 and CT-2A. A pilot spatial transcriptomics experiment (n=2 animals per OR group) was performed to test how presence of ORlow and ORhigh tumors affect the tumor microenvironment and distribution of cells within the tissue.

Our raw dataset is a large dataframe in R storing gene expression data for 4 samples (two ORhigh, two ORlow). There are 32,272 genes (stored as columns) with 16,085 spots (stored as rows) total. 

Note that `STdeconvolve` has been already performed on the current data set with outputs as `deconProp` and `deconGexp`.

-   `STdeconvolve` (an unsupervised machine learning approach) enables reference-free cell-type deconvolution of multi-cellular pixel-resolution spatial transcriptomics data [(See original paper on Nature Communications by Miller et al. 2022).](https://www.nature.com/articles/s41467-022-30033-z) From a gene expression matrix, STdeconvolve selects genes most likely to be relevant for distinguishing between cell types by looking for highly overdispersed genes across ST pixels. It also filters out spots with too few genes and genes with too few reads. 
-   `deconProp` (data type: matrix) represents the proportion of each deconvolved cell-type across each spatially resolved pixel or spot
-   `deconGexp` (data type: matrix) represents the putative gene expression profile for each deconvolved cell-type normalized to a library size of 1

### Goal 

The main biological question in this data analysis is as follows: 

**Is there a difference in GBM tumor microenvironment depending on the OR type?** 

Using statistical analyses and data visualization, we attempted to explore this research question by looking for differential gene expression patterns in the tumor-infiltrating region, if any, to identify how cell-types and processes are differentially organized in the GBM peripheral regions based on each OR group. 

### Analyses

Multiple statistical analyses were performed, including: 

- Pearson Correlation 
- Spearman Correlation 
- 2-sample t-tests (one-sided) 
- Wilcox Rank-Sum Test 
- Fold change 

### Discussion 

After exploring the dataset with multiple analyses methods and tests, we observe that there seems to be differences in the tumor microenvironment for `ORlow` vs. `ORhigh`. 

The tumor infiltrative regions for `ORlow` tumors seem to be enriched with immunoglobulin (Ig) genes, while `ORhigh` tumors show almost negligible Ig expression. This may indicate a greater composition of immune cells in the `ORlow` tumor microenvironment compared to that of `ORhigh`. 

Additionally, gene expression tended to be higher around `ORlow` tumors compared to `ORhigh` tumors. We couldn't specify a particular family of genes that showed higher expression around `ORhigh` tumor infiltrating regions while being nearly absent in expression around `ORlow` tumor infiltrating regions.  

