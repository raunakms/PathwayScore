# PathwayScore

PathwayScore is collection of R-functions to compute pathway z-score given gene-expression data and genes mapped to corresponding pathways.

The rationale behind this concept is to find out how the pathway is behaving (either activating or repressing) given the expression pattern of expression of multiple genes.
  - For every gene, in expression data, might have different scales of expression-value.
  - So to bring down expression-values of every gene/protein in the same scale we compute Z-score (see below)
  - And the pathway-score is the average of z-score (from above) of all genes mapped to that pathway.

# Z-score
```sh
Z-score = x_ij - m_i / sd_i
```
 Here,
  - i = gene
  - j = sample
  - x_ij = expression value of gene i in sample j
  - m_i = mean of expression values in for gene i across all j
  - sd_i = standard-deviation of expression values in for gene i across all j


# Usage
```sh
get.pathway.zcore(file.expr, file.pathway)
```
# Description:
  - file.expr = path of file containing gene-expression matrix with Row as Genes and Columns as Samples.
  - file.pathway = path of file containing pathway-gene relationship
   -- each line should contain information about a unique pathway
   -- first-column must contain name of the pathwayj
   -- second-column must contain genes mapped to the corresponding pathway multiple genes must be separated by comma symbol (,)
