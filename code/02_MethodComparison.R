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
library(phyloseq)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TABLE <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_TOOL <- "/data/yzwang/functions/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/S1/"

source(paste0(DIR_TOOL, "wilcox_diff.R"))

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
  annotate("text", x = 2.2, y = 13, colour = 1,
           label = paste0("Rs=",round(cor_haemophilus$r, 2),
                          "\np=",format.pval(cor_haemophilus$p.adj, 2))) +
  xlab("RNA-seq") +
  ylab("16S rRNA") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        plot.title = element_text(face = "italic")) 
pdf(file.path(DIR_RES, "B_correlation_density_case_haemophilus.pdf"), 
    width = 4, height = 3.5)
print(p_density_cor_haemophilus)
dev.off()

##################
# Diff comparation
abund30_16s <- names(sort(rowSums(cpm_16s), decreasing = T))[1:30]
abund30_rna <- names(sort(rowSums(cpm_rna), decreasing = T))[1:30]

cutoffFC <- log2(1.3)
diff_16s <- wilcox_diff(cpm_16s, cpm_16s[ , grep("C", colnames(cpm_16s))],
                        cpm_16s[ , grep("N", colnames(cpm_16s))], 
                        T, 0.1, cutoffFC)
diff_rna <- wilcox_diff(cpm_rna, cpm_rna[ , grep("C", colnames(cpm_rna))],
                        cpm_rna[ , grep("N", colnames(cpm_rna))], 
                        T, 0.1, cutoffFC)

draw_diff <- merge(diff_16s, diff_rna, by = "ID")

draw_diff <- draw_diff %>%
  filter(ID %in% union(abund30_rna, abund30_16s)) %>%
  # Create helper variables for cleaner logic
  mutate(
    x_high = log2FC.x > cutoffFC,
    x_low = log2FC.x < -cutoffFC,
    x_mid_pos = log2FC.x > 0 & log2FC.x <= cutoffFC,
    x_mid_neg = log2FC.x < 0 & log2FC.x >= -cutoffFC,
    
    y_high = log2FC.y > cutoffFC,
    y_low = log2FC.y < -cutoffFC,
    y_mid_pos = log2FC.y > 0 & log2FC.y <= cutoffFC,
    y_mid_neg = log2FC.y < 0 & log2FC.y >= -cutoffFC
  ) %>%
  mutate(
    Colour = case_when(
      x_high & y_high ~ "brown4", # Both highly-up
      (x_high & y_mid_pos) | (y_high & x_mid_pos) ~ "brown2", # One highly-up, one moderately-up
      x_low & y_low ~ "steelblue4", # Both highly-dn
      (x_low & y_mid_neg) | (y_low & x_mid_neg) ~ "steelblue2", # One highly-dn, one moderately-dn
      abs(log2FC.x) < cutoffFC & abs(log2FC.y) < cutoffFC ~ "grey80", # Both NS
      TRUE ~ "grey80"
    )
  ) %>%
  select(-c(x_high, x_low, x_mid_pos, x_mid_neg, y_high, y_low, y_mid_pos, y_mid_neg)) %>%
  na.omit()

diff_highlight <- draw_diff %>%
  filter(ID %in% union(c("Staphylococcus", "Helicobacter", "Streptococcus", 
                         "Klebsiella", "Haemophilus"),
                       genus_highlight))
diff_highlight.up <- diff_highlight %>%
  filter(log2FC.x > cutoffFC | log2FC.y > cutoffFC)
diff_highlight.dn <- diff_highlight %>%
  filter(log2FC.x < -cutoffFC | log2FC.y < -cutoffFC)

p_diff_points <- 
  ggplot(draw_diff, aes(x = log2FC.y, y = log2FC.x)) +
  geom_point(colour = draw_diff$Colour,
             size = if_else(draw_diff$P.adj.x < 0.05 & draw_diff$P.adj.y < 0.05, 3, 2)) +
  geom_point(data = draw_diff[draw_diff$P.adj.x < 0.05 & draw_diff$P.adj.y < 0.05, ],
             aes(x = log2FC.y, y = log2FC.x),
             size = 3, shape = 1, colour = "black") +
  ggrepel::geom_label_repel(diff_highlight.up, mapping = aes(label = ID),
                            nudge_x = 1, nudge_y = 1, force = 1, 
                            alpha = 0.7, fontface = "italic") +
  ggrepel::geom_label_repel(diff_highlight.dn, mapping = aes(label = ID),
                            nudge_x = -1, nudge_y = -1, force = 1, 
                            alpha = 0.7, fontface = "italic") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -cutoffFC, linetype = "dashed") +
  geom_hline(yintercept = cutoffFC, linetype = "dashed") +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = -cutoffFC, linetype = "dashed") +
  geom_vline(xintercept = cutoffFC, linetype = "dashed") +
  xlab("RNA-seq") + 
  ylab("16S rRNA") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        axis.text.x = element_text(colour = "black"))
pdf(file.path(DIR_RES, "C_differential_points.pdf"), width = 3.5, height = 3.5)
print(p_diff_points)
dev.off()

#############################################
# 16S-detected differential genera in RNA-seq
rownames(draw_diff) <- draw_diff$ID
sig_names <- draw_diff$ID[draw_diff$Change.x != "NS"]
sig_diff <- data.frame(
  genus = rep(sig_names, 2),
  group = c(rep("16s", length(sig_names)),
            rep("rna", length(sig_names)))
) %>%
  mutate(log2FC = c(draw_diff[sig_names, "log2FC.x"],
                    draw_diff[sig_names, "log2FC.y"]),
         p.adj = c(draw_diff[sig_names, "P.adj.x"],
                   draw_diff[sig_names, "P.adj.y"]),
         change = c(draw_diff[sig_names, "Change.x"],
                    draw_diff[sig_names, "Change.y"])) %>%
  arrange(group, desc(log2FC))

sig_diff$fc_change <- ifelse(sig_diff$log2FC > 0, "up", "dn")

p_diff_lever <- 
  ggplot(sig_diff, aes(x = factor(genus, levels = unique(genus)), y = log2FC)) +
  geom_line(aes(group = genus), color = "grey25", linetype = "dashed", size = 0.75) +
  geom_point(aes(shape = group, color = fc_change), 
             size = 3.5) +
  geom_text(data = subset(sig_diff, change != "NS"),
            aes(x = factor(genus, levels = unique(genus)), 
                y = log2FC - 0.1),
            label = "*", size = 5, 
            color = "white", vjust = 0.5, hjust = 0.5) +
  theme(legend.position = "top") +
  scale_color_manual(values = c("up" = "#E4672D", "dn" = "#2C7094")) +
  geom_hline(yintercept = 0, color = "grey", linetype="dashed") +
  xlab("") +
  ylim(-3, 4.5) +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        axis.text.x = element_text(angle = 45, face = "italic",
                                   hjust = 1, vjust = 1, colour = "black"))
pdf(file.path(DIR_RES, "D_differential_lever.pdf"), width = 8, height = 4)
print(p_diff_lever)
dev.off()
