# Microbiome Data Analysis and Visualization in R
# Author: [Arif Rahman]
# Description: This script performs microbiome diversity analysis, including relative abundance, heatmaps, dendrogram, PCoA, and statistical testing (PERMANOVA, Kruskal-Wallis).

# =======================
# Load Required Libraries
# =======================
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(openxlsx)
library(pheatmap)
library(ggpubr)
library(ggdendro)

# =======================
# Bar Plot: Microbiota Composition
# =======================

file_path <- ("your_path/mikrobiota.xlsx")
data <- read_excel(file_path, sheet = "Sheet1")

individual_order <- c("Ind1", "Ind2",...)

microbiota_columns <- c("Firmicutes", "Bacteroidetes",...... )

data_filtered <- data %>% 
  filter(SampleID %in% individual_order) %>% 
  mutate(SampleID = factor(SampleID, levels = individual_order))

data_filtered[is.na(data_filtered)] <- 0

data_long <- data_filtered %>%
  pivot_longer(cols = all_of(microbiota_columns), 
               names_to = "Microbiota", 
               values_to = "Abundance") %>%
  filter(Abundance > 0)

data_long <- data_long %>%
  group_by(SampleID) %>%
  mutate(Percentage = (Abundance / sum(Abundance)) * 100) %>%
  ungroup()

data_long <- data_long %>%
  group_by(SampleID) %>%
  mutate(Microbiota = factor(Microbiota, levels = unique(Microbiota[order(-Percentage)]))) %>%
  ungroup()

custom_colors <- c(
  "#2D6A4F", "#1B4332", "#A7C957", "#D9BF77", "#3A7CA5", 
  "#457B9D", "#A8DADC", "#F4A261", "#2A9D8F", "#E63946", 
  "#264653", "#6A0572", "#EEA47F", "#375A7F", "#A084CA", 
  "#B5838D", "#FFC107", "#D4A5A5", "#709775", "#A6A57A"
)

color_palette <- rep(custom_colors, length.out = length(unique(data_long$Microbiota)))

ggplot(data_long, aes(x = SampleID, y = Percentage, fill = Microbiota)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "",
    x = "",
    y = "Percentage Relative Abundance (%)",
    fill = "Microbiota phylum"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(ncol = 2))

getwd()
setwd("your_path")
ggsave("Composition of the gut microbiota.png", width = 12, height = 7.5, dpi = 500, units = "in")

# =======================
# Heatmap: Relative Abundance
# =======================
data_heatmap <- read_excel("your_path/mikrobiota.xlsx")
row.names(data_heatmap) <- data_heatmap$Sample_ID
microbiota_data <- data_heatmap[, -1]
microbiota_data[] <- lapply(microbiota_data, as.numeric)
rel_abund <- sweep(microbiota_data, 1, rowSums(microbiota_data), FUN = "/")
rel_abund_t <- t(rel_abund)

pheatmap(
  rel_abund_t,
  scale = "column",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of Relative Abundance of Orangutan Microbiota",
  fontsize_row = 8,
  fontsize_col = 10
)

# ========================
# Alpha Diversity (Shannon & ASV Richness)
# ========================
data_alpha <- read_excel("your_path/Shannon_ASV.xlsx")
row.names(data_alpha) <- data_alpha$Sampel_ID
data_alpha$Month <- as.factor(data_alpha$Month)
data_microbiota <- data_alpha[, -c(1, 2)]
data_microbiota[] <- lapply(data_microbiota, as.numeric)

shannon_index <- diversity(data_microbiota, index = "shannon")
asv_richness <- rowSums(data_microbiota > 0)

diversity_df <- data.frame(
  Individu = row.names(data_alpha),
  Shannon_Index = shannon_index,
  ASV_Richness = asv_richness,
  Month = data_alpha$Month
)

kruskal_shannon <- kruskal.test(Shannon_Index ~ Month, data = diversity_df)
kruskal_asv <- kruskal.test(ASV_Richness ~ Month, data = diversity_df)

print(kruskal_shannon)
print(kruskal_asv)

# Shannon Boxplot
boxplot_shannon <- ggplot(diversity_df, aes(x = Month, y = Shannon_Index, fill = Month)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(x = "Month", y = "Shannon Index") +
  theme_bw() +
  stat_compare_means(method = "kruskal.test")

ggsave("Shannon_Month.jpeg", plot = boxplot_shannon, width = 6, height = 5, dpi = 300)

# ASV Richness Boxplot
boxplot_asv <- ggplot(diversity_df, aes(x = Month, y = ASV_Richness, fill = Month)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(x = "Month", y = "ASV Richness") +
  theme_bw() +
  stat_compare_means(method = "kruskal.test")

ggsave("ASV_Month.jpeg", plot = boxplot_asv, width = 6, height = 5, dpi = 300)

# Combine Boxplots
combined_plot <- ggarrange(
  boxplot_shannon, boxplot_asv,
  labels = c("C", "D"),
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"
)

ggsave("Combined_Shannon_ASV_Month.jpeg", plot = combined_plot, width = 12, height = 7, dpi = 300)

# =======================
# Beta Diversity (PERMANOVA, Dendrogram & PCoA)
# =======================
data_meta <- read_excel("your_path/Dendogram_PCoA.xlsx")
metadata <- select(data_meta, SampleID, Sex, Month)
microbiota_matrix <- select(data_meta, where(is.numeric)) %>% as.matrix()

empty_rows <- rowSums(microbiota_matrix) == 0
microbiota_matrix <- microbiota_matrix[!empty_rows, ]
metadata <- metadata[!empty_rows, ]
rownames(microbiota_matrix) <- metadata$SampleID

microbiota_matrix <- microbiota_matrix + 1
microbiota_matrix <- microbiota_matrix / rowSums(microbiota_matrix)

bray_dist <- vegdist(microbiota_matrix, method = "bray")
permanova_result <- adonis2(bray_dist ~ Sex + Month, data = metadata, permutations = 999)
print(permanova_result)

# PCoA
pcoa_result <- cmdscale(bray_dist, eig = TRUE, k = 2)
pcoa_df <- data.frame(
  SampleID = rownames(microbiota_matrix),
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2],
  Sex = metadata$Sex,
  Month = metadata$Month
)

pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Sex, shape = Month)) +
  geom_point(size = 4) +
  geom_text(aes(label = SampleID), vjust = 2, size = 3, check_overlap = TRUE) +
  theme_bw() +
  labs(
    title = "PCoA Bray-Curtis Dissimilarity",
    x = paste0("PC1 (", round(pcoa_result$eig[1] / sum(pcoa_result$eig) * 100, 1), "%)"),
    y = paste0("PC2 (", round(pcoa_result$eig[2] / sum(pcoa_result$eig) * 100, 1), "%)")
  ) +
  scale_color_manual(values = c("salmon", "steelblue"))
print(pcoa_plot)
ggsave("PCoA_BrayCurtis_with_PERMANOVA.png", plot = pcoa_plot, width = 9, height = 6, dpi = 500)

# Dendrogram
hc <- hclust(bray_dist, method = "average")
dendro_data <- ggdendro::dendro_data(hc)
dendro_data$labels$label <- rownames(microbiota_matrix)[hc$order]

dendrogram_plot <- ggplot() +
  geom_segment(data = dendro_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dendro_data$labels, aes(x = x, y = 0, label = label), size = 3) +
  theme_minimal() +
  labs(title = "Dendrogram - Bray-Curtis", x = "Orangutan", y = "Distance")
print(dendrogram_plot)
ggsave("Dendrogram.png", plot = dendrogram_plot, width = 10, height = 7)
