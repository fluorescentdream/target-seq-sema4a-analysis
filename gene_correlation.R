install.packages("readr")
install.packages("patchwork")
install.packages("ggplot2")
install.packages("tidyverse")
library(patchwork)
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

gene_results <- data.frame(
  Gene = TARGET_seq_raw_counts[[1]],
  Correlation = gene_correlations
)
# %>% is like applying the filter to gene_results based on the condition inside. take gene_results and sort the rows by the Correlation column in descending order
gene_results <- gene_results %>% na.omit() %>% filter(Gene != "SEMA4A") %>% arrange(desc(Correlation))

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

# which cell types express GCK?!? 
gck_means <- TARGET_seq_metadata %>%  filter(!is.na(Cluster)) %>% group_by(Cluster) %>% summarise(mean_gck = mean(GCK)) %>% arrange(desc(mean_gck))

figure10 <- ggplot(gck_means, aes(x = fct_reorder(Cluster, mean_gck), y = mean_gck)) + 
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Mean GCK Expression by Cell Cluster",
    x = "Cell Cluster",
    y = "Mean GCK Expression"
  )

figure10
save_plot(figure10, "Figure10_mean_gck_expression_by_cluster")

#hsc/mpp specific analysis - as hsc/gmp is where gck is most expressed
# part 1: % of cells expressing SEMA4A at all, by GCK detection status

# group_by splits data into cells where GCK was detected, and cells where is wasnt (TRUE/FALSE). every line after runs separately on each group
# n = n() counts how many cells in each group
# pct_SEMA4A_pos = mean(SEMA4A > 0) * 100 - SEMA4A > 0 creates TRUE/FALSE for each cell, based on if it expresses any/or is 0. 
# R treats T/F as 1/0 values, so taking the mean() of a bunch of TRUE/FALSE values gives the proportion that are TRUE. multiplying by 100 converts that proportion to a percentage
hsc_data <- TARGET_seq_metadata %>% filter(Cluster == "HSC/MPP")
hsc_data <- hsc_data %>% mutate(GCK_detected = factor(GCK > 0, labels = c("FALSE", "TRUE")))
prop_data <- hsc_data %>% group_by(GCK_detected) %>% summarise(n = n(), pct_SEMA4A_pos = mean(SEMA4A > 0) * 100)
prop_data

p1 <- ggplot(prop_data, aes(x = GCK_detected, y = pct_SEMA4A_pos, fill = GCK_detected)) +
  geom_col(width = 0.5) +
  # rounds the percentage to 1 decimal and glues a "%" symbol onto it as text. vjust = -0.5 nudges the text label slightly upward so it sits above the bar instead of overlapping it.
  geom_text(aes(label = paste0(round(pct_SEMA4A_pos, 1), "%")), vjust = -0.5) +
  theme_bw(base_size = 12) +
  labs(x = "GCK detection status", y = "% cells expressing SEMA4A",
       title = "A. Proportion of SEMA4A+ cells by GCK detection status") +
  theme(legend.position = "none", plot.title = element_text(size = 9))

# part 2: among SEMA4A-expressing cells only, what's the expression level?
hsc_data_pos <- hsc_data %>% filter(SEMA4A > 0)

p2 <- ggplot(hsc_data_pos, aes(x = GCK_detected, y = SEMA4A)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, fill = "steelblue", alpha = 0.3) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.3) +
  # plots every individual cell as a small dot, nudged randomly left/right (width = 0.15) so points don't stack in a perfectly straight vertical line and become unreadable.
  scale_y_continuous(trans = "pseudo_log", breaks = c(1, 10, 100, 1000, 3000)) +
  theme_bw(base_size = 12) +
  labs(x = "GCK detection status", y = "SEMA4A Expression (SEMA4A+ cells only, pseudo-log)",
       title = "B. SEMA4A Expression Level Among\nSEMA4A+ Cells, by GCK Status") +
  theme(plot.title = element_text(size = 9))

figure11 <- p1 + p2 + plot_layout(widths = c(1, 1.2))
figure11 <- figure11 +
  plot_annotation(
    title = "SEMA4A Expression by GCK Detection Status in HSC/MPP Cells",
    subtitle = "n = 4787 (GCK not detected), n = 164 (GCK detected)"
  )
figure11

save_plot(figure11, "Figure11_sema4a_gck_combined_hscmpp_cells")

# statistical backing for figure 11 - tests whether the two binary variables (GCK detection, SEMA4A detection) are statistically independent of each other, or whether there's a real association.
# result: not significant. GCK-detected cells are not meaningfully more or less likely to express SEMA4A at all.
table(hsc_data$GCK_detected, hsc_data$SEMA4A > 0) %>% fisher.test()

# statistical backing for figure 12 - runs on hsc_data_pos, the SEMA4A>0-only subset. so it's specifically testing "among cells that express SEMA4A, is the expression level different between GCK+ and GCK- cells?"
# also not significant at conventional thresholds, though it's noticeably lower than the Fisher's p-value and the original whole-dataset Wilcoxon
# so excluding zeros did sharpen the signal somewhat, consistent with what looked like a real shift in the boxplots. but it doesn't cross the line into statistical significance
wilcox.test(SEMA4A ~ GCK_detected, data = hsc_data_pos)