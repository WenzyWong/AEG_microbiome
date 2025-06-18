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
DIR_TOOL <- "/data/yzwang/functions/"
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

###################
# Alpha diversity #
shan_aeg <- apply(count_aeg, 2, diversity)
shan_stad <- apply(count_gc, 2, diversity)
shan_escc <- apply(count_escc[ , grep("GSE235537", colnames(count_escc))], 2, diversity)

shan_violin <- data.frame(
  Shannon = c(shan_aeg, shan_escc, shan_stad),
  Cancer = c(rep("AEG", length(shan_aeg)),
             rep("ESCC", length(shan_escc)),
             rep("STAD", length(shan_stad)))
)

p_shan_violin <-
  ggplot(shan_violin, aes(x = factor(Cancer, levels = c("ESCC", "AEG", "STAD")), 
                        y = Shannon, fill = Cancer)) +
  geom_violin() +
  scale_fill_manual(values = c("#485682", "#60AB9E", "#5C8447")) +
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