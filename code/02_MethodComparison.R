#####################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 02. Verify the feasibility of RNAseq-derived
# microbiome analysis by comparing 16S and RNAseq
#
#####################################################

library(paletteer)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TABLE <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_TOOL <- "/data/yzwang/functions/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/S1/"

#########################
# Dataset: gastric cancer
cpm_16s <- readRDS(paste0(DIR_RDS, "gGC_CPM_16S_comparison.rds"))
cpm_rna <- readRDS(paste0(DIR_RDS, "gGC_CPM_RNA_comparison.rds"))

# Overall correlation between 16S and RNA-seq
mean_16s <- log2(apply(cpm_16s[ , grep("C", colnames(cpm_16s))], 1, mean) + 1)
mean_rna <- log2(apply(cpm_rna[ , grep("C", colnames(cpm_rna))], 1, mean) + 1)

cor_mean_res <- psych::corr.test(mean_16s, mean_rna, 
                                 use = "pairwise", method = "spearman", adjust = "holm")

all(names(mean_16s) == names(mean_rna))
genus_abundant <- names(mean_16s)[mean_16s != 0 & mean_rna != 0]
length(genus_abundant)
cor_res <- psych::corr.test(t(cpm_16s[genus_abundant, ]), t(cpm_rna[genus_abundant, ]), 
                            use = "pairwise", method = "spearman", adjust = "holm")
cor_r_within <- sort(diag(as.matrix(cor_res$r)), decreasing = T)

draw_cor <- data.frame(
  Genus = names(mean_16s),
  Using_16S = mean_16s,
  Using_RNAseq = mean_rna
)

genus_highlight <- names(cor_r_within[1:5])
draw_highlight <- draw_cor[draw_cor$Genus %in% genus_highlight, ]

p_overall_cor <-
  ggplot(draw_cor, mapping = aes(x = Using_RNAseq, y = Using_16S)) + 
  geom_point(color = "grey80", size = 1.2) +
  geom_smooth(method = lm, color = "#4388BF", linewidth = 2, fill = "grey80") +
  ggtitle("GC tumour microbiome") +
  annotate("text", x = 1.6, y = 16, color = 1, 
           label = paste0("Rs=", round(cor_mean_res$r, 2),
                          "\np", format.pval(cor_mean_res$p.adj, 2))) +
  xlab("RNA-seq") +
  ylab("16S rRNA-seq") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        legend.title = element_blank(),
        legend.position = "bottom") +
  geom_point(draw_highlight, mapping = aes(x = Using_RNAseq, y = Using_16S),
             color = "#F48118", size = 3.5) +
  ggrepel::geom_label_repel(draw_highlight, mapping = aes(label = Genus),
                            color = "black", nudge_x = -1, nudge_y = -2,
                            force = 1, alpha = 0.7, fontface = "italic")

pdf(file.path(DIR_RES, "A_overall_correlation_16S_RNAseq.pdf"), width = 4, height = 4)
print(p_overall_cor)
dev.off()

# Case: Helicobacter
df_helicobacter <- data.frame(
  row.names = colnames(cpm_16s),
  method_16s = log2(as.numeric(cpm_16s["Helicobacter", ]) + 1), 
  method_rna = log2(as.numeric(cpm_rna["Helicobacter", ]) + 1))

cor_helicobacter <- psych::corr.test(df_helicobacter$method_16s, df_helicobacter$method_rna,
                                     use = "pairwise", method = "spearman")
p_density_cor_helicobacter <-
  ggplot(df_helicobacter, aes(x = method_16s, y = method_rna)) +
  geom_point(size = 1.2, colour = "grey80") +
  geom_smooth(method = lm, linewidth = 2, color = "#F48118", fill = "grey80") +
  stat_density2d(aes(colour = after_stat(level)), size = .5, alpha = 0.33) +
  scale_color_gradient(low = '#4F94CD85', high = '#A52A2A85') +
  ggtitle("Helicobacter - log2CPM") +
  annotate("text", x = 1.5, y = 18, colour = 1,
           label = paste0("Rs=",round(cor_helicobacter$r, 2),
                          "\np",format.pval(cor_helicobacter$p.adj, 2))) +
  xlab("RNA-seq") +
  ylab("16S rRNA") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        plot.title = element_text(face = "italic")) 
pdf(file.path(DIR_RES, "B_correlation_density_case_helicobacter.pdf"), width = 4, height = 3.5)
print(p_density_cor_helicobacter)
dev.off()

# Case: Haemophilus
df_haemophilus <- data.frame(
  row.names = colnames(cpm_16s),
  method_16s = log2(as.numeric(cpm_16s["Haemophilus", ]) + 1), 
  method_rna = log2(as.numeric(cpm_rna["Haemophilus", ]) + 1))

cor_haemophilus <- psych::corr.test(df_haemophilus$method_16s, df_haemophilus$method_rna,
                                    use = "pairwise", method = "spearman")
p_density_cor_haemophilus <-
  ggplot(df_haemophilus, aes(x = method_16s, y = method_rna)) +
  geom_point(size = 1.2, colour = "grey80") +
  geom_smooth(method = lm, linewidth = 2, color = "#F48118", fill = "grey80") +
  stat_density2d(aes(colour = after_stat(level)), size = .5, alpha = 0.33) +
  scale_color_gradient(low = '#4F94CD85', high = '#A52A2A85') +
  ggtitle("Haemophilus - log2CPM") +
  annotate("text", x = 1.5, y = 13.5, colour = 1,
           label = paste0("Rs=",round(cor_haemophilus$r, 2),
                          "\np",format.pval(cor_haemophilus$p.adj, 2))) +
  xlab("RNA-seq") +
  ylab("16S rRNA") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        plot.title = element_text(face = "italic")) 
pdf(file.path(DIR_RES, "B_correlation_density_case_haemophilus.pdf"), width = 4, height = 3.5)
print(p_density_cor_haemophilus)
dev.off()
