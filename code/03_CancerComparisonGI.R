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
DIR_TABLE <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_TOOL <- "/data/yzwang/git_project/AEG_microbiome/utils/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/F1/"

##############
# Data loading
cpm_escc <- readRDS(paste0(DIR_RDS, "gESCC_cpm.rds"))
cpm_escc <- cpm_escc[ , grep("GSE235537", colnames(cpm_escc))]
cpm_gc <- readRDS(paste0(DIR_RDS, "gGC_CPM_RNA_WithoutCompFilt.rds"))
cpm_aeg <- readRDS(paste0(DIR_RDS, "gAEG_CPM_RNA_WithoutCompFilt.rds"))

count_escc <- readRDS(paste0(DIR_RDS, "gESCC_count.rds"))
count_escc <- count_escc[ , grep("GSE235537", colnames(count_escc))]
count_gc <- readRDS(paste0(DIR_RDS, "gGC_count_RNA.rds"))
count_aeg <- readRDS(paste0(DIR_RDS, "gAEG_count.rds"))

ESCC_selected_data <- as.data.frame(read_excel(file.path(DIR_TABLE, "ESCC_selected_data.xlsx")))
rownames(ESCC_selected_data) <- ESCC_selected_data$Run
colnames(cpm_escc) <- gsub("GSE235537.", "", colnames(cpm_escc))
colnames(cpm_escc) <- paste0(ESCC_selected_data[colnames(cpm_escc), "Type"],
                             ESCC_selected_data[colnames(cpm_escc), "Patient"])

#################
# Alpha diversity
shan_aeg <- apply(count_aeg, 2, diversity)
shan_stad <- apply(count_gc, 2, diversity)
shan_escc <- apply(count_escc[ , grep("GSE235537", colnames(count_escc))], 2, diversity)

shan_violin <- data.frame(
  Shannon = c(shan_aeg, shan_escc, shan_stad),
  Cancer = factor(c(rep("AEG", length(shan_aeg)),
                    rep("ESCC", length(shan_escc)),
                    rep("STAD", length(shan_stad))),
                  levels = c("ESCC", "AEG", "STAD"))
)

p_shan_violin <-
  ggplot(shan_violin, aes(x = Cancer, y = Shannon, fill = Cancer)) +
  geom_violin() +
  scale_fill_manual(values = c("#60AB9E", "#485682", "#5C8447")) +
  geom_boxplot(width = .2, fill = "grey75",
               outlier.colour = NA) +
  geom_jitter(alpha = .8, color = "grey75",
              size = .2) +
  xlab("") +
  ylab("Shannon Index") + 
  stat_compare_means(label.x = 1.5, label.y = 4.5) +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        legend.position = "top",
        legend.direction = "horizontal") 
pdf(file.path(DIR_RES, "A_alpha_diversity_among_cancer.pdf"), height = 3.2, width = 3)
print(p_shan_violin)
dev.off()

################
# Beta diversity
source(file.path(DIR_TOOL, "plot_beta_diversity.R"))

overlap <- intersect(rownames(cpm_aeg), rownames(cpm_gc))
overlap <- intersect(overlap, rownames(cpm_escc))

tumour_aeg <- cpm_aeg[overlap, grepl("C", colnames(cpm_aeg))]
tumour_escc <- cpm_escc[overlap, grepl("T", colnames(cpm_escc))]
tumour_stad <- cpm_gc[overlap, grepl("CT", colnames(cpm_gc))]
beta.tumour <- plot_beta_diversity(
  abundance_list = list(AEG = tumour_aeg, ESCC = tumour_escc, STAD = tumour_stad),
  colors = c("#B24745FF", "#80796BFF", "#DF8F44FF"),
  output_file = file.path(DIR_RES, "B_Beta_diversity_tumours.pdf")
)

normal_aeg <- cpm_aeg[overlap, grepl("N", colnames(cpm_aeg))]
normal_escc <- cpm_escc[overlap, grepl("N", colnames(cpm_escc))]
normal_stad <- cpm_gc[overlap, grepl("NT", colnames(cpm_gc))]
beta.normal <- plot_beta_diversity(
  abundance_list = list(AEG = normal_aeg, ESCC = normal_escc, STAD = normal_stad),
  colors = c("#485682", "#60AB9E", "#5C8447"),
  output_file = file.path(DIR_RES, "B_Beta_diversity_normals.pdf")
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

pdf(file.path(DIR_RES, "C_Jaccard_dist_all.pdf"), width = 6, height = 7)
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

pdf(file.path(DIR_RES, "C_Jaccard_index.pdf"), width = 5, height = 4.5)
Heatmap(jaccard_idx, 
        name = "Jaccard\nIndex",
        col = RColorBrewer::brewer.pal(name = "OrRd", n = 4),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE)
dev.off()

#### 