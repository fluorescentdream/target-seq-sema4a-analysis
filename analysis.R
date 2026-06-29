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
#View(TARGET_seq_raw_counts)

TARGET_seq_metadata <- read_tsv("data/TARGET-seq_metadata.tsv")
#View(TARGET_seq_metadata)

#dim(TARGET_seq_metadata)
#dim(TARGET_seq_raw_counts)
#colnames(TARGET_seq_raw_counts)[1:5]

sema4a <- TARGET_seq_raw_counts[which(TARGET_seq_raw_counts[[1]] == "SEMA4A"),]
#dim(sema4a)

sema4a_expression <- as.numeric(sema4a[1, -1])
#length(sema4a_expression)

TARGET_seq_metadata$SEMA4A <- sema4a_expression
#head(TARGET_seq_metadata)

summary(TARGET_seq_metadata$SEMA4A)

table(TARGET_seq_metadata$Genotype)
aggregate(SEMA4A ~ Genotype, data = TARGET_seq_metadata, FUN = mean)
# slightly elevated in CH-associated genotypes compared to WT controls
boxplot(SEMA4A ~ Genotype, data = TARGET_seq_metadata, las = 2, cex.axis = 0.7, main = "SEMA4A Expression by Genotype")

table(TARGET_seq_metadata$Cluster)
# clustered by cell type. plasma cells express sema4a most. then pDC. 
aggregate(SEMA4A ~ Cluster, data = TARGET_seq_metadata, FUN = mean)
boxplot(SEMA4A ~ Cluster, data = TARGET_seq_metadata, las = 2, cex.axis = 0.7, main = "SEMA4A Expression by Cell Cluster")

# CH vs control expression. CH seems more
table(TARGET_seq_metadata$Sample_type)
aggregate(SEMA4A ~ Sample_type, data = TARGET_seq_metadata, FUN = mean)
boxplot(SEMA4A ~ Sample_type, data = TARGET_seq_metadata, las = 2, cex.axis = 0.7, main = "SEMA4A Expression by Sample Type")

# top 20 SEMA4A-expressing cells
top <- TARGET_seq_metadata %>% arrange(desc(SEMA4A))
head(top$Cell, 20)

# tapply(values, groups, function). splits values into groups and applies a function
# which cell populations express SEMA4A the most? - plasma cells, endothelial, pDC top 3. 
sort(tapply(TARGET_seq_metadata$SEMA4A, TARGET_seq_metadata$Cluster, mean), decreasing = TRUE)

# largest genotype groups (WT-CH, WT-control, DNMT3A, and TET2) contain hundreds to thousands of cells, so their average SEMA4A expression values are more reliable
# HSC/MPP cells make up a large portion of the dataset across multiple genotypes
# because SEMA4A may naturally be more highly expressed in HSC/MPP cells, differences between genotypes could reflect differences in cell-type composition rather than a direct effect of the mutation itself
# although plasma cells showed the highest average SEMA4A expression, there were only 18 plasma cells in the dataset. 
# aditionally, the high SEMA4A expression observed in the DNMT3A-SF3B1 group cant be explained by plasma cells, since that genotype contains no plasma cells.
table(TARGET_seq_metadata$Genotype, TARGET_seq_metadata$Cluster)

aggregate(SEMA4A ~ Cluster, data = subset(TARGET_seq_metadata, Genotype == "DNMT3A"), mean)

# is sema4a expressed in most cells or only a subset ? most cells express it very little. 
figure1 <- ggplot(TARGET_seq_metadata, aes(x = SEMA4A)) +
            geom_histogram(bins = 50, fill = "steelblue", color = "black") +
            scale_x_continuous(trans = "log1p") +
            theme_bw() +
            labs(
              title = "Distribution of SEMA4A Expression",
              x = "SEMA4A Expression",
              y = "Number of Cells"
            )

save_plot(figure1, "Figure1_SEMA4A_distribution")


# highest expressing clusters on top because of fct_reorder
figure2 <- ggplot(TARGET_seq_metadata, aes(x = fct_reorder(Cluster, SEMA4A, .fun = mean), y=SEMA4A)) + 
            geom_violin(fill = "steelblue", color = "black", scale = "width") +
            coord_flip() +
            scale_y_continuous(trans = "log1p") +
            theme_bw(base_size = 12, base_family = "serif") +
            labs(
              title = "SEMA4A Expression by Cell Cluster",
              x = "Cell Cluster",
              y = "SEMA4A Expression"
            ) +
            geom_boxplot(width = 0.1, outlier.size = 0.3)

save_plot(figure2, "Figure2_SEMA4A_by_cluster")


figure3 <- ggplot(TARGET_seq_metadata, aes(x = fct_reorder(Genotype, SEMA4A, .fun = mean), y=SEMA4A)) + 
            geom_boxplot(fill = "steelblue", color = "black") +
            coord_flip() +
            scale_y_continuous(trans = "log1p") +
            theme_bw(base_size = 12, base_family = "serif") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(
              title = "SEMA4A Expression by Genotype",
              x = "Genotype",
              y = "SEMA4A Expression"
            )

save_plot(figure3, "Figure3_SEMA4A_by_genotype")


figure4 <- ggplot(TARGET_seq_metadata, aes(x = fct_reorder(Sample_type, SEMA4A, .fun = mean), y=SEMA4A)) + 
            geom_violin(fill = "steelblue", color = "black", scale = "width") +
            coord_flip() +
            scale_y_continuous(trans = "log1p") +
            theme_bw(base_size = 12, base_family = "serif") +
            labs(
              title = "SEMA4A Expression by Sample Type",
              x = "Sample Type",
              y = "Expression (log1p)"
            ) 

save_plot(figure4, "Figure4_SEMA4A_by_sample")


cluster_means <- aggregate(SEMA4A ~ Cluster, data = TARGET_seq_metadata, FUN = mean)
figure5 <- ggplot(cluster_means, aes(x = fct_reorder(Cluster, SEMA4A), y=SEMA4A)) + 
            geom_col(fill = "steelblue") +
            coord_flip() +
            theme_classic() +
            labs(
              title = "Mean SEMA4A Expression by Cell Cluster",
              x = "Cell Cluster",
              y = "Mean Expression"
            )

save_plot(figure5, "Figure5_mean_SEMA4A_by_cluster")

