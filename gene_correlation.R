install.packages("readr")
library(readr)
library(ggplot2)
library(tidyverse)

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

TARGET_seq_metadata <- read_tsv("data/TARGET-seq_metadata.tsv")

sema4a <- TARGET_seq_raw_counts %>% filter(Gene == "SEMA4A")
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

# %>% is like applying the filter to gene_results based on the condition inside
gene_results <- gene_results %>% filter(Gene != "SEMA4A")
gene_results <- na.omit(gene_results)
#take gene_results and sort the rows by the Correlation column in descending order
gene_results <- gene_results %>% arrange(desc(Correlation))
gene_results <- gene_results[gene_results$Gene != "SEMA4A", ]

top20 <- head(gene_results, 20)
# forces bars to appear sorted from highest to lowest
top20$Gene <- factor(top20$Gene, levels = rev(top20$Gene))

figure6 <- ggplot(top20, aes(x = Gene, y = Correlation)) + 
            geom_col(fill = "steelblue", color = "black") +
            coord_flip() +
            theme_bw(base_size = 12) +
            labs(
              title = "Top Genes Correlated with SEMA4A",
              x = "Gene",
              y = "Correlation with SEMA4A"
            )

save_plot(figure6, "Figure6_SEMA4A_top20_correlated_genes")

gck <- TARGET_seq_raw_counts %>% filter(Gene == "GCK")
gck_expression <- as.numeric(gck[1, -1])
TARGET_seq_metadata$GCK <- gck_expression

figure7 <- ggplot(TARGET_seq_metadata, aes(x = GCK, y = SEMA4A)) +
  geom_point(color = "black", alpha = 0.6, size = 0.2) +
  # transform the axes to spread out the bottom-left cluster
  scale_x_continuous(trans = "pseudo_log") +
  scale_y_continuous(trans = "pseudo_log") +
  # mini graphs for each cluster
  facet_wrap(~Cluster, scales = "free") +
  geom_smooth(method = "lm", color = "steelblue", se = FALSE) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6)
  ) +
  labs(
    title = "Correlation of GCK and SEMA4A Expression",
    x = "GCK Expression",
    y = "SEMA4A Expression"
  )
figure7
save_plot(figure7, "Figure7_SEMA4A_GCK_correlation_across_clusters")


cluster_means <- aggregate(SEMA4A ~ Cluster, data = TARGET_seq_metadata, FUN = mean)
figure8 <- ggplot(cluster_means, aes(x = reorder(Cluster, SEMA4A), y = SEMA4A)) +
  geom_segment(aes(xend = Cluster, y = 0, yend = SEMA4A), linewidth = 0.7, color = "grey") +
  geom_point(size = 3, color = "steelblue") +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(
    title = "Mean SEMA4A Expression by Cell Cluster",
    x = "Cell Cluster",
    y = "Mean SEMA4A Expression"
  )

figure8
save_plot(figure8, "Figure8_mean_SEMA4A_expression_across_clusters")

head(TARGET_seq_metadata[, c("Cluster", "SEMA4A", "GCK")])
# take targetseq metadata -> group rows by Cluster -> summarize each group by calculating the correlation between SEMA4A and GCK. 
# explicitly check whether either gene has any variation before calculating the correlation.
cluster_correlations <- TARGET_seq_metadata %>% group_by(Cluster) %>% summarise(Correlation =
      if (
        sd(SEMA4A) == 0 ||
        sd(GCK) == 0
      ) {
        NA
      } else {
        cor(SEMA4A, GCK)
      }
  )

#omitted the NA values 
cluster_correlations <- na.omit(cluster_correlations)
cluster_correlations <- cluster_correlations %>% arrange(desc(Correlation))
cluster_correlations

figure9 <- ggplot(cluster_correlations, aes(x = fct_reorder(Cluster, Correlation), y = Correlation)) + 
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Correlation between GCK and SEMA4A by Cell Cluster",
    x = "Cell Cluster",
    y = "Pearson Correlation"
  )

figure9
save_plot(figure9, "Figure9_correlation_GCK_SEMA4A_by_cluster")
