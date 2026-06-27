install.packages("readr")
library(readr)
library(ggplot2)
library(tidyverse)

Sys.setenv("VROOM_CONNECTION_SIZE" = 10000000)
TARGET_seq_raw_counts <- read_tsv("data/TARGET-seq_raw_counts_matrix.txt.gz")

TARGET_seq_metadata <- read_tsv("data/TARGET-seq_metadata.tsv")

sema4a <- TARGET_seq_raw_counts[which(TARGET_seq_raw_counts[[1]] == "SEMA4A"),]
sema4a_expression <- as.numeric(sema4a[1, -1])
TARGET_seq_metadata$SEMA4A <- sema4a_expression

# gene1 contains the expression of one gene across all 14,073 cells.
gene1 <- as.numeric(TARGET_seq_raw_counts[1, -1])
length(gene1)

# correlation between expression of first gene in dataset and sema4a. 1 is positive correlation (when gene1 is high, SEMA4A is usually high too.), -1 is perfect negative. 
cor(gene1, sema4a_expression)

#[, -1] drops first column of Gene, 1 means row by row since each row is one gene. temporary function(x) for each row; for each gene, calculate correlation with SEMA4A
gene_correlations <- apply(TARGET_seq_raw_counts[, -1], 1, function(x) cor(as.numeric(x), sema4a_expression))
length(gene_correlations)
sum(is.na(gene_correlations)) # 2445 NA correlations because gene expression never changes, so 0 standard deviation. there is nothing to correlate with SEMA4A

# R treats a data frame as a list of columns using [column] - you can just isolate one particular column
gene_results <- data.frame(Gene = TARGET_seq_raw_counts[[1]], Correlation = gene_correlations)
gene_results <- na.omit(gene_results)
gene_results <- gene_results[order(gene_results$Correlation, decreasing = TRUE),]
gene_results <- gene_results[gene_results$Gene != "SEMA4A", ]

top20 <- head(gene_results, 20)
# forces bars to appear sorted from highest to lowest
top20$Gene <- factor(top20$Gene, levels = rev(top20$Gene))

ggplot(top20, aes(x = Gene, y = Correlation)) + 
  geom_col(fill = "steelblue", color = "black", scale = "width") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Top Genes Correlated with SEMA4A",
    x = "Gene",
    y = "Correlation with SEMA4A"
  )

ggsave(
  "figures/Figure6_SEMA4A_top20_correlated_genes.png",
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/Figure6_SEMA4A_top20_correlated_genes.pdf",
  width = 8,
  height = 6
)

