##TARGET-seq SEMA4A Analysis

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

##dplyr 

- filter() - keeps rows
- select() - keeps columns
- arrange() - sorts rows
- mutate() - creates new columns
- summarise() - calculate statistics

## top 20 gene correlation information
- RP11-550A5.2 - long non-coding RNA (lncRNA) gene. referenced in research regarding the transthyretin (TTR) gene and methylation sites associated with heart disease in certain populations
- ACBD3-AS1 - long non-coding RNA (lncRNA) gene. regulating the expression of the nearby ACBD3 gene, which is a crucial membrane adaptor involved in Golgi apparatus maintenance, lipid metabolism, and viral host-pathogen interactions. tentatively linked to inherited eye disorders such as Retinitis Pigmentosa
- RP11-817G13.3 - intergenic long non-coding RNA (lncRNA). does not have a well-defined, established role
- RP11-271M24.2 - long intergenic non-coding RNA (lincRNA). no clear role.
- SAPCD1-AS1 - SAPCD1 Antisense RNA 1. long non-coding RNA (lncRNA) gene. notably associated with Giant Cell Glioblastoma, an aggressive type of brain tumor. 
- GCK - glucokinase gene. provides instructions for making the enzyme glucokinase, which functions as the body's primary glucose sensor. signals pancreatic beta cells to release insulin when blood sugar rises. facilitates glucose absorption/storage in liver.
- CH17-343M10.2 - long non-coding RNA (lncRNA). no clue what it does
- BTC - betacellulin. encodes a protein that belongs to the Epidermal Growth Factor (EGF) family. acts as a potent mitogen (induces mitosis/cell division). highly expressed in the pancreas and is thought to play a major role in the differentiation and development of pancreatic beta cells.
- CTB-102L5.8 - human long non-coding RNA (lncRNA) classified as a sense intronic RNA. on chromosome 19q13.2. function unknown
- ARFGEF1-DT - long non-coding RNA (lncRNA) gene located on human chromosome 8.  it interacts with chromatin and proteins to regulate gene expression, such as binding to partners like ILF2. upregulated in some cancers (such as T-cell acute lymphoblastic leukemia) where it aids in tumor cell proliferation. 
- TMEM202 - encodes a protein that localizes to the plasma membrane. critical for cell-to-cell communication, maintaining membrane structure, and facilitating paracellular transport (the movement of ions and molecules between cells). 
- CBS - provides instructions for making the enzyme cystathionine β-synthase. regulates homocysteine levels and produces cysteine, which is essential for building proteins and producing antioxidants like glutathione.
- AC002064.5 - antisense long non-coding RNA (lncRNA) on human chromosome 7 (locus 7q21.13). regulating the expression of nearby or downstream genes and interacts with RNA-binding proteins or microRNAs. biomarker in several gene-expression profiles and survival risk models used to study cancers, including pancreatic cancer and clear cell renal cell carcinoma.
- AC000085.4 - encodes the CLTCL1 (Clathrin heavy chain like 1) gene, which plays a key role in intracellular trafficking. major component of coated vesicles involved in endocytosis and intracellular protein sorting. Highly expressed in skeletal muscles, the placenta, and the brain, regulating metabolic processes such as glucose transporter (GLUT4) trafficking.
- OR14L1P - Olfactory Receptor Family 14 Subfamily L Member 1 Pseudogene. inactive human genetic sequence. while the active OR14L1 gene aids in the G protein-mediated signal transduction of smell perception, the OR14L1P sequence has acquired mutations over the course of human evolution that render it a non-functional pseudogene.
- TMX4 - Thioredoxin-Related Transmembrane Protein 4. endoplasmic reticulum (ER) protein that plays a key role in protein folding and cellular redox regulation. assist in quality control during protein synthesis.
- RP11-702H23.4 - uncharacterized human Long Non-Coding RNA (lncRNA) gene. on chromosome 11q13.4. generally associated with post-transcriptional gene regulation, chromatin remodeling, and acting as competing endogenous RNAs (ceRNAs) that interact with microRNAs to regulate protein-coding genes. 
- CTB-50L17.9 - is an uncharacterized human long non-coding RNA (lncRNA) gene located on chromosome 17. links genes in its vicinity to epigenetic processes and the regulation of m⁶A RNA methylation in neurological and neurovascular conditions, such as vascular dementia and ischemic stroke. 
- RP11-356O9.2 - long non-coding RNA (lncRNA) on chromosome 14. lncRNAs from the RP11 group (including this locus) are frequently studied because they form fusion transcripts or undergo structural variations during tumorigenesis. aberrant expression or chromosomal rearrangements involving these regions can serve as potential prognostic or diagnostic biomarkers for aggressive cancers.
- OR2T1 - Olfactory Receptor Family 2 Subfamily T Member 1 gene. provides instructions for making an olfactory receptor protein responsible for detecting odor molecules and initiating the sense of smell. 

##GCK Correlation Analysis
objective: investigate the relationship between GCK and SEMA4A expression across cell populations.

tasks completed:
- refactored repeated ggsave() calls into a reusable save_plot() function
- extracted GCK expression and added it to the metadata table
- created scatterplots of GCK vs SEMA4A expression across all cell clusters
- calculated GCK-SEMA4A correlation within each cell cluster
- visualized mean SEMA4A expression and cluster-specific correlations

observations:
- GCK was the highest-ranked protein-coding gene positively correlated with SEMA4A
- overall GCK-SEMA4A correlation varied substantially between cell clusters
- strongest positive correlations were observed in HSC/MPP, Small Pre-B, and Monocyte populations
- several clusters produced NA correlations due to zero variance in gene expression

