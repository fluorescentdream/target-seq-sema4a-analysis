##TARGET-seq SEMA4A Analysis - Day 1

objective: load the published TARGET-seq dataset into R and perform an initial exploratory analysis of SEMA4A expression across cells.

tasks completed:
- successfully imported TARGET-seq metadata.
- successfully imported processed raw count matrix.
- matched SEMA4A expression values to each cell using cell IDs.
- explored the distribution of SEMA4A expression.
- compared SEMA4A expression across genotypes.
- compared SEMA4A expression across cell clusters.
- compared SEMA4A expression between CH and control samples.
- generated supa cool figures using ggplot2!

observations:
- plasma cells exhibited the highest average SEMA4A expression, although only 18 plasma cells were present in the dataset.
- endothelial cells and plasmacytoid dendritic cells also showed relatively high expression.
- WT-CH samples demonstrated moderately higher average SEMA4A expression than WT-control samples.
- differences between genotypes may partially reflect differences in cell-type composition rather than direct effects of the mutations.

challenges: 
- importing the large count matrix needed increasing the VROOM connection buffer.
- the metadata file initially loaded weirdly due to using the wrong import function.
- matching the expression matrix with metadata needed understanding that cell IDs aligned exactly between the two files.

next steps:
- investigate genes whose expression correlates with SEMA4A.
- begin statistical testing of observed expression differences.
- potentially use ML to predict responses based on features. 

##SEMA4A Correlation Expression Analysis

objective: identify genes correlated with SEMA4A expression.

tasks completed:
- loaded TARGET-seq single-cell RNA-seq raw counts and metadata into R
- extracted SEMA4A expression across all cells
- assigned SEMA4A expression values to each cell in the metadata table
- computed gene-wise correlation between SEMA4A and all ~36000 genes across 14073 cells

observations:
- many genes showed near-zero or undefined correlation with SEMA4A
- some correlations were NA due to genes with zero variance (no expression across cells)
- final output: a correlation value for every gene vs SEMA4A
- biological inference: high positive correlation suggests co-expression across the same cell states

next steps: 
- identify top positively and negatively correlated genes
- investigate whether correlated genes are enriched in specific pathways or cell types and look at biological interpretations
