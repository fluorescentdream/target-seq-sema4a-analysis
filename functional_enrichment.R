install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("biomaRt")

library(dplyr)
library(readr)
library(ggplot2)
library(forcats)
library(stringr)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(biomaRt)

save_plot <- function(plot, filename, width = 8, height = 6) {
  ggsave(
    paste0("figures/", filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
  
  ggsave(
    paste0("figures/", filename, ".pdf"),
    plot = plot,
    width = width,
    height = height
  )
}

Sys.setenv("VROOM_CONNECTION_SIZE" = 10000000)

TARGET_seq_raw_counts <- read_tsv("data/TARGET-seq_raw_counts_matrix.txt.gz")

sema4a <- TARGET_seq_raw_counts %>%
  filter(Gene == "SEMA4A")

sema4a_expression <- as.numeric(sema4a[1, -1])

gene_correlations <- apply(TARGET_seq_raw_counts[, -1], 1, function(x) cor(as.numeric(x), sema4a_expression))

gene_results <- data.frame(Gene = TARGET_seq_raw_counts[[1]], Correlation = gene_correlations)
# %>% is like applying the filter to gene_results based on the condition inside. take gene_results and sort the rows by the Correlation column in descending order
gene_results <- gene_results %>% na.omit() %>% filter(Gene != "SEMA4A") %>% arrange(desc(Correlation))

top_genes <- gene_results %>% slice_head(n = 100) %>% pull(Gene)
top_genes

# connecting R to a huge online biology database called Ensembl
# useMart - connects to a human gene database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# getBM = download info from ensembl
# attributes = asking for gene name and gene type (protein coding, lcrna, pseudogene, etc)
# filters = gives gene names
# values = based on top genes, get it for that list
biotype_map <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters = "hgnc_symbol",
  values = top_genes,
  mart = mart
)

protein_coding <- biotype_map %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# Gene %in% protein_coding - KEEP ONLY GENES THAT ARE IN CLEANED PROTEIN-CODING LIST. take top 100 AFTER cleaning
filtered_top <- gene_results %>% filter(Gene %in% protein_coding) %>% slice_head(n = 100)

genes<- filtered_top$Gene
gene_map <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# gene = filtered_top$Gene: cleaned gene list
# OrgDb = org.Hs.eg.db: human gene annotation database
# ont = "BP": biological processes
# pAdjustMethod = "BH": statistical correction (controls false positives)
# readable = TRUE: converts IDs into readable gene names


ego <- enrichGO(
  gene = gene_map$ENTREZID,
  universe = bitr(TARGET_seq_raw_counts$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

nrow(ego@result)                  # total terms tested
sum(ego@result$p.adjust < 0.05)   # how many pass
head(ego@result[order(ego@result$pvalue), ], 10)  # look at the best ones regardless of cutoff
length(bitr(TARGET_seq_raw_counts$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)
figure12 <- dotplot(ego)
figure12
save_plot(figure12, "Figure12_gene_enrichment_protein_coding")
