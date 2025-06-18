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
overlap <- intersect(rownames(cpm_aeg), rownames(cpm_gc))
overlap <- intersect(overlap, rownames(cpm_escc))

over_aeg <- cpm_aeg[overlap, ]
over_escc <- cpm_escc[overlap, ]
over_stad <- cpm_gc[overlap, ]

distBray <- vegdist(t(cbind(over_aeg, over_escc, over_stad)), method = "bray")
distMtx <- as.matrix(distBray)

brayPCoA <- pcoa(distBray)

pcoaRes <- brayPCoA$values[ , "Relative_eig"]
pcoa1 <- as.numeric(sprintf("%.3f", pcoaRes[1])) * 100
pcoa2 <- as.numeric(sprintf("%.3f", pcoaRes[2])) * 100

pcoaDraw <- as.data.frame(brayPCoA$vectors)[1:2]
pcoaDraw$Samples <- colnames(cbind(over_aeg, over_escc, over_stad))
pcoaDraw$Group <- c(rep("AEG", ncol(over_aeg)),
                    rep("ESCC", ncol(over_escc)),
                    rep("STAD", ncol(over_stad)))

meanPcoa1 <- c(
  mean(pcoaDraw$Axis.1[pcoaDraw$Group == "AEG"]),
  mean(pcoaDraw$Axis.1[pcoaDraw$Group == "ESCC"]),
  mean(pcoaDraw$Axis.1[pcoaDraw$Group == "STAD"]))
meanPcoa2 <- c(
  mean(pcoaDraw$Axis.2[pcoaDraw$Group == "AEG"]),
  mean(pcoaDraw$Axis.2[pcoaDraw$Group == "ESCC"]),
  mean(pcoaDraw$Axis.2[pcoaDraw$Group == "STAD"]))

meanPoint <- data.frame(X = meanPcoa1, Y = meanPcoa2,
                        Group = c("AEG", "ESCC", "STAD"))

p_pcoa <-
  ggplot(pcoaDraw, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size = 0.7) +
  geom_segment(aes(x = Axis.1, xend = c(rep(meanPoint[1, 1], ncol(over_aeg)),
                                        rep(meanPoint[2, 1], ncol(over_escc)),
                                        rep(meanPoint[3, 1], ncol(over_stad))),
                   y = Axis.2, yend = c(rep(meanPoint[1, 2], ncol(over_aeg)),
                                        rep(meanPoint[2, 2], ncol(over_escc)),
                                        rep(meanPoint[3, 2], ncol(over_stad)))),
               linewidth = 0.2) +
  scale_color_manual(values = c("#485682", "#60AB9E", "#5C8447")) + 
  labs(x = paste("PCoA1(", pcoa1, "%)", sep = ""),
       y = paste("PCoA2(", pcoa2, "%)", sep = ""),
       title = "PCoA: Bray-Curtis Distance") +
  geom_hline(yintercept = 0, linetype = 4, color = "grey20", alpha = 0.6) + 
  geom_vline(xintercept = 0, linetype = 4, color = "grey20", alpha = 0.6) + 
  stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.05) +
  theme_classic() +
  theme(axis.text = element_text(colour = 1),
        legend.position = "top",
        legend.direction = "horizontal")
pdf(file.path(DIR_RES, "B_beta_diversity_among_cancer.pdf"), height = 3.4, width = 3)
print(p_pcoa)
dev.off()
