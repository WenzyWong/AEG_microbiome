#######################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 03. Compare AEG microbiome with other GI cancers,
# including GC, AEG, ESCC
#
#######################################################
library(readxl)
library(ggvenn)
library(ggalluvial)
library(ggpubr)
library(ggplot2)
library(vegan) # diversity, vegist
library(ape) # pcoa
library(paletteer)
library(ComplexHeatmap)
library(RColorBrewer)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TOOL <- "/data/yzwang/git_project/AEG_microbiome/utils/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/F1/"

##############
# Data loading
cpm_escc <- readRDS(paste0(DIR_RDS, "gESCC_cpm_2cohorts.rds"))
cpm_stad <- readRDS(paste0(DIR_RDS, "gGC_CPM_RNA_WithoutCompFilt.rds"))
cpm_aeg <- readRDS(paste0(DIR_RDS, "gAEG_CPM_RNA_WithoutCompFilt.rds"))

count_escc <- readRDS(paste0(DIR_RDS, "gESCC_count_2cohorts.rds"))
count_stad <- readRDS(paste0(DIR_RDS, "gGC_count_RNA.rds"))
count_aeg <- readRDS(paste0(DIR_RDS, "gAEG_count.rds"))

#################
# Alpha diversity
shan_aeg <- apply(count_aeg, 2, diversity)
shan_stad <- apply(count_stad, 2, diversity)
shan_escc <- apply(count_escc, 2, diversity)

shan_violin <- data.frame(
  Shannon = c(shan_aeg, shan_escc, shan_stad),
  Cancer = factor(c(rep("AEG", length(shan_aeg)),
                    rep("ESCC", length(shan_escc)),
                    rep("STAD", length(shan_stad))),
                  levels = c("ESCC", "AEG", "STAD"))
)

shan_comp_list <- list(c("AEG", "ESCC"),
                       c("AEG", "STAD"),
                       c("ESCC", "STAD"))

p_shan_violin <-
  ggplot(shan_violin, aes(x = Cancer, y = Shannon, fill = Cancer)) +
  geom_violin() +
  scale_fill_manual(values = c("#60AB9EFF", "#485682FF", "#5C8447FF")) +
  geom_boxplot(width = .2, fill = "grey75",
               outlier.colour = NA) +
  geom_jitter(alpha = .8, color = "grey75",
              size = .2) +
  xlab("Cancer Type") +
  ylab("Shannon Index") + 
  stat_compare_means(comparisons = shan_comp_list) +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1)) 
pdf(file.path(DIR_RES, "A_alpha_diversity_among_cancer.pdf"), height = 3.3, width = 4)
print(p_shan_violin)
dev.off()

################
# Beta diversity
source(file.path(DIR_TOOL, "plot_beta_diversity.R"))

overlap <- intersect(rownames(cpm_aeg), rownames(cpm_stad))
overlap <- intersect(overlap, rownames(cpm_escc))

beta <- plot_beta_diversity(
  abundance_list = list(AEG = cpm_aeg, 
                        ESCC = cpm_escc, 
                        STAD = cpm_stad),
  colors = c("#485682", "#60AB9E", "#5C8447"),
  output_file = file.path(DIR_RES, "B_beta_diversity_all.pdf")
)

tumour_aeg <- cpm_aeg[overlap, grepl("C", colnames(cpm_aeg))]
tumour_escc <- cpm_escc[overlap, grepl("_T|_D", colnames(cpm_escc))]
tumour_stad <- cpm_stad[overlap, grepl("CT", colnames(cpm_stad))]
beta.tumour <- plot_beta_diversity(
  abundance_list = list(AEG = tumour_aeg, ESCC = tumour_escc, STAD = tumour_stad),
  colors = c("#B24745FF", "#80796BFF", "#DF8F44FF"),
  output_file = file.path(DIR_RES, "B_beta_diversity_tumours.pdf")
)

normal_aeg <- cpm_aeg[overlap, grepl("N", colnames(cpm_aeg))]
normal_escc <- cpm_escc[overlap, grepl("_N", colnames(cpm_escc))]
normal_stad <- cpm_stad[overlap, grepl("NT", colnames(cpm_stad))]
beta.normal <- plot_beta_diversity(
  abundance_list = list(AEG = normal_aeg, ESCC = normal_escc, STAD = normal_stad),
  colors = c("#485682", "#60AB9E", "#5C8447"),
  output_file = file.path(DIR_RES, "B_beta_diversity_normals.pdf")
)

##################
# Jaccard distance
# All samples
all_samples <- cbind(tumour_aeg, tumour_escc, tumour_stad, 
                     normal_aeg, normal_escc, normal_stad)
jaccard_dist <- as.matrix(dist(t(all_samples), method = "binary"))

groups_ord <- c(
  rep("AEG_T", ncol(tumour_aeg)),
  rep("ESCC_T", ncol(tumour_escc)),
  rep("STAD_T", ncol(tumour_stad)),
  rep("AEG_N", ncol(normal_aeg)),
  rep("ESCC_N", ncol(normal_escc)),
  rep("STAD_N", ncol(normal_stad))
)
names(groups_ord) <- colnames(all_samples)

group_anno <- HeatmapAnnotation(
  Group = groups_ord,
  col = list(Group = c(AEG_T = "#B24745FF",
                       ESCC_T = "#80796BFF",
                       STAD_T = "#DF8F44FF",
                       AEG_N = "#485682FF",
                       ESCC_N = "#60AB9EFF",
                       STAD_N = "#5C8447FF"))
)

pdf(file.path(DIR_RES, "C_jaccard_dist_all.pdf"), width = 7, height = 6)
Heatmap(jaccard_dist, 
        name = "Jaccard\nDistance",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        col = RColorBrewer::brewer.pal(name = "RdBu", n = 11),
        top_annotation = group_anno,
        show_row_names = F,
        show_column_names = F)
dev.off()

# Group means
groups <- cbind(
  AEG_T = rowMeans(tumour_aeg), AEG_N = rowMeans(normal_aeg),
  ESCC_T = rowMeans(tumour_escc), ESCC_N = rowMeans(normal_escc),
  STAD_T = rowMeans(tumour_stad), STAD_N = rowMeans(normal_stad)
)

jaccard_idx <- 1 - as.matrix(dist(t(groups), method = "binary"))

pdf(file.path(DIR_RES, "C_jaccard_index.pdf"), width = 5, height = 4.5)
Heatmap(jaccard_idx, 
        name = "Jaccard\nIndex",
        col = RColorBrewer::brewer.pal(name = "OrRd", n = 4),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE)
dev.off()

###########################
# Abundant genus comparison
abund_escc <- sort(apply(cpm_escc, MARGIN = 1, FUN = mean) / 1e+4, decreasing = T)
abund_aeg <- sort(apply(cpm_aeg, MARGIN = 1, FUN = mean) / 1e+4, decreasing = T)
abund_stad <- sort(apply(cpm_stad, MARGIN = 1, FUN = mean) / 1e+4, decreasing = T)

per_escc <- c(abund_escc[1:10], "Others" = 100 - sum(abund_escc[1:10]))
per_aeg <- c(abund_aeg[1:10], "Other" = 100 - sum(abund_aeg[1:10]))
per_stad <- c(abund_stad[1:10], "Other" = 100 - sum(abund_stad[1:10]))

per_escc <- data.frame(Genus = names(per_escc), Abundance = per_escc)
per_aeg <- data.frame(Genus = names(per_aeg), Abundance = per_aeg)
per_stad <- data.frame(Genus = names(per_stad), Abundance = per_stad)

abund_venn <- list(
  ESCC = per_escc$Genus,
  AEG = per_aeg$Genus,
  STAD = per_stad$Genus
)

pdf(file.path(DIR_RES, "D_venn_top10.pdf"), width = 3, height = 3)
ggvenn::ggvenn(abund_venn, c("ESCC", "AEG", "STAD"),
               show_percentage = F,
               fill_color = c("#60AB9EFF", "#485682FF", "#5C8447FF"))
dev.off()

col_escc <- c("Cutibacterium" = "#4269A099", 
              "Pseudomonas"= "#F89D59", 
              "Prevotella" = "#60AB9E88", 
              "Streptomyces" = "#FEEDAA",
              "Bacillus" = "#ADE0DC",
              "Escherichia" = "#F8CD5A",
              "Agrococcus" = "#8CBC6888", 
              "Staphylococcus" = "#9CCEE3", 
              "Stenotrophomonas" = "#EA9696",
              "Streptococcus" = "#9A342C",
              "Others" = "#DCDCDC")

col_aeg <- c("Staphylococcus" = "#64A4CC", 
             "Pasteurella" = "#9CCEE3", 
             "Bacillus" = "#ADE0DC", 
             "Klebsiella" = "#CAEAD2", 
             "Priestia" = "#D3F2D6", 
             "Providencia" = "#ECF6C8", 
             "Streptomyces" = "#FEEDAA", 
             "Sphingomonas" = "#FDC980", 
             "Pseudomonas" = "#F89D59", 
             "Burkholderia" = "#E75B3A", 
             "Other" = "#DCDCDC")

col_stad <- c("Pasteurella" = "#9CCEE3", 
              "Staphylococcus" = "#64A4CC", 
              "Bacillus" = "#ADE0DC", 
              "Klebsiella" = "#CAEAD2", 
              "Priestia" = "#D3F2D6",
              "Vibrio" = "#A6D384", 
              "Escherichia" = "#F8CD5A", 
              "Streptomyces" = "#FEEDAA", 
              "Stenotrophomonas" = "#EA9696", 
              "Streptococcus" = "#9A342C",
              "Other" = "#DCDCDC")

p_top_escc <-
  ggplot(per_escc) +
  geom_bar(aes(x = 1, y = Abundance, 
               fill = factor(Genus, levels = unique(Genus))), 
           stat = "identity", position = "stack") +
  scale_fill_manual(values = col_escc,
                    name = "Genus") +
  xlab("ESCC") +
  ylab("Relative abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing.x = element_blank())
pdf(file.path(DIR_RES, "D_abundance_top10_escc.pdf"), width = 3, height = 7)
print(p_top_escc)
dev.off()

p_top10_aeg <-
  ggplot(per_aeg) +
  geom_bar(aes(x = 1, y = Abundance, 
               fill = factor(Genus, levels = unique(Genus))), 
           stat = "identity", position = "stack") +
  scale_fill_manual(values = col_aeg,
                    name = "Genus") +
  xlab("AEG") +
  ylab("Relative abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing.x = element_blank())
pdf(file.path(DIR_RES, "D_abundance_top10_aeg.pdf"), width = 3, height = 7)
print(p_top_escc)
dev.off()

p_top10_stad <-
  ggplot(per_stad) +
  geom_bar(aes(x = 1, y = Abundance, 
               fill = factor(Genus, levels = unique(Genus))), 
           stat = "identity", position = "stack") +
  scale_fill_manual(values = col_stad,
                    name = "Genus") +
  xlab("STAD") +
  ylab("Relative abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing.x = element_blank())
pdf(file.path(DIR_RES, "D_abundance_top10_stad.pdf"), width = 3, height = 7)
print(p_top_escc)
dev.off()
