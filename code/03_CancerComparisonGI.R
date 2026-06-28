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
library(dplyr)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TOOL <- "/data/yzwang/git_project/AEG_microbiome/utils/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/F1_cancer_types/"

##############
# Data loading
cpm_escc <- readRDS(paste0(DIR_RDS, "gESCC_cpm_2cohorts.rds"))
cpm_stad <- readRDS(paste0(DIR_RDS, "gGC_CPM_RNA_WithoutCompFilt.rds"))
cpm_aeg <- readRDS(paste0(DIR_RDS, "gAEG_CPM_RNA_FiltMyco.rds"))

count_escc <- readRDS(paste0(DIR_RDS, "gESCC_count_2cohorts.rds"))
count_stad <- readRDS(paste0(DIR_RDS, "gGC_count_RNA.rds"))
count_aeg <- readRDS(paste0(DIR_RDS, "gAEG_Count_RNA_FiltMyco.rds"))

#################
# Alpha diversity
shan_aeg <- apply(count_aeg, 2, vegan::diversity)
shan_stad <- apply(count_stad, 2, vegan::diversity)
shan_escc <- apply(count_escc, 2, vegan::diversity)

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

pdf(file.path(DIR_RES, "A_alpha_diversity_among_cancer.pdf"), height = 3.3, width = 4)
ggplot(shan_violin, aes(x = Cancer, y = Shannon, fill = Cancer)) +
  geom_violin() +
  scale_fill_manual(values = c("#74ACCE", "#EEC549", "#8B87BA")) +
  geom_boxplot(width = .2, fill = "grey75",
               outlier.colour = NA) +
  geom_jitter(alpha = .8, color = "grey75",
              size = .2) +
  xlab("Cancer Type") +
  ylab("Shannon Index") +
  stat_compare_means(comparisons = shan_comp_list) +
  theme_test() +
  theme(aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))
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
  colors = c("#EEC549", "#74ACCE", "#8B87BA"),
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
        show_column_names = F,
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2")
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

per_escc <- c(abund_escc[1:20], "Others" = 100 - sum(abund_escc[1:20]))
per_aeg <- c(abund_aeg[1:20], "Other" = 100 - sum(abund_aeg[1:20]))
per_stad <- c(abund_stad[1:20], "Other" = 100 - sum(abund_stad[1:20]))

per_escc <- data.frame(Genus = names(per_escc), Abundance = per_escc)
per_aeg <- data.frame(Genus = names(per_aeg), Abundance = per_aeg)
per_stad <- data.frame(Genus = names(per_stad), Abundance = per_stad)

abund_venn <- list(
  ESCC = per_escc$Genus,
  AEG = per_aeg$Genus,
  STAD = per_stad$Genus
)

pdf(file.path(DIR_RES, "D_venn_top20.pdf"), width = 3, height = 3)
ggvenn::ggvenn(abund_venn, c("ESCC", "AEG", "STAD"),
               show_percentage = F,
               fill_color = c("#74ACCE", "#EEC549", "#8B87BA"))
dev.off()

# Palette shared with the AEG overall distribution (abund_aeg_distribution)
colAbund <- c(paletteer_d("khroma::discreterainbow")[c(10, 12:20, 23:27, 2, 4, 5, 7, 9)],
              "grey")

# AEG top 20: colours consistent with abund_aeg_distribution (by abundance rank)
aeg_top20 <- per_aeg$Genus[1:20]
col_aeg <- setNames(colAbund[1:20], aeg_top20)
col_aeg["Other"] <- "#DCDCDC"

# ESCC/STAD: genera shared with AEG top20 reuse AEG colours; the rest get
# Oranges/Greens (higher abundance -> darker, avoiding the lightest shades).
build_col <- function(per_df, brewer_name) {
  g <- per_df$Genus
  g <- g[!(g %in% c("Other", "Others"))]
  shared <- g[g %in% aeg_top20]
  novel <- g[!(g %in% aeg_top20)]
  cols <- character(0)
  cols[shared] <- col_aeg[shared]
  if (length(novel) > 0) {
    cols[novel] <- colorRampPalette(rev(brewer.pal(9, brewer_name)[2:9]))(length(novel))
  }
  cols["Other"] <- "#DCDCDC"
  cols["Others"] <- "#DCDCDC"
  cols
}

col_escc <- build_col(per_escc, "Oranges")
col_stad <- build_col(per_stad, "Greens")

######################################
# Combined top20 abundance (shared legend)
aeg_top20  <- per_aeg$Genus[1:20]
escc_g     <- setdiff(per_escc$Genus, c("Other", "Others"))
stad_g     <- setdiff(per_stad$Genus, c("Other", "Others"))
escc_novel <- setdiff(escc_g, aeg_top20)
stad_novel <- setdiff(stad_g, aeg_top20)
cross_g    <- intersect(escc_novel, stad_novel)   # non-AEG, in both ESCC & STAD
escc_only  <- setdiff(escc_novel, cross_g)
stad_only  <- setdiff(stad_novel, cross_g)

col_escc_only <- setNames(colorRampPalette(rev(brewer.pal(9, "Blues")[2:9]))(length(escc_only)), escc_only)
col_aeg_all   <- setNames(colAbund[1:20], aeg_top20)
col_stad_only <- setNames(colorRampPalette(rev(brewer.pal(9, "Purples")[2:9]))(length(stad_only)), stad_only)
col_cross     <- setNames(colorRampPalette(brewer.pal(9, "Greys")[5:8])(length(cross_g)), cross_g)

genus_levels <- c(escc_only, aeg_top20, stad_only, cross_g, "Other")
col_combined <- c(col_escc_only, col_aeg_all, col_stad_only, col_cross, "Other" = "#DCDCDC")
col_combined <- col_combined[genus_levels]

abund_combined <- rbind(
  data.frame(Cancer = "ESCC", per_escc),
  data.frame(Cancer = "AEG",  per_aeg),
  data.frame(Cancer = "STAD", per_stad)
)
rownames(abund_combined) <- NULL
abund_combined$Genus[abund_combined$Genus %in% c("Other", "Others")] <- "Other"
abund_combined$Cancer <- factor(abund_combined$Cancer, levels = c("ESCC", "AEG", "STAD"))

# Stack each bar by its OWN abundance (highest at top), Other forced to the bottom.
# Compute explicit y limits so the stacking order is decoupled from the grouped
# legend order (legend keeps: ESCC-only / AEG / STAD-only / cross / Other).
abund_combined <- abund_combined %>%
  mutate(ord_key = ifelse(Genus == "Other", -1, Abundance)) %>%
  arrange(Cancer, desc(ord_key)) %>%
  group_by(Cancer) %>%
  mutate(ymax = sum(Abundance) - cumsum(Abundance) + Abundance,
         ymin = sum(Abundance) - cumsum(Abundance)) %>%
  ungroup()
abund_combined$xc <- as.integer(abund_combined$Cancer)

pdf(file.path(DIR_RES, "D_abundance_top20_combined.pdf"), width = 6.5, height = 8)
print(
  ggplot(abund_combined) +
    geom_rect(aes(xmin = xc - 0.45, xmax = xc + 0.45,
                  ymin = ymin, ymax = ymax,
                  fill = factor(Genus, levels = genus_levels)),
             colour = NA) +
    scale_fill_manual(values = col_combined, name = "Genus") +
    scale_x_continuous(breaks = 1:3, labels = c("ESCC", "AEG", "STAD")) +
    xlab("") +
    ylab("Relative abundance (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(colour = 1),
          axis.text.y = element_text(colour = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.spacing.x = element_blank(),
          legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 7))
)
dev.off()


######################################
# AEG-enriched genera vs the other cancers
# Candidates = AEG top20 (by mean relative abundance). Keep those significantly
# HIGHER in AEG than BOTH ESCC and STAD: Wilcoxon on relative abundance,
# BH-adjusted padj < 0.05 in both comparisons, and AEG median above both.
aeg_cand <- names(sort(rowMeans(cpm_aeg), decreasing = TRUE))[1:20]

ra_vals <- function(mat, g) {
  if (g %in% rownames(mat)) as.numeric(mat[g, ]) / 1e4 else rep(0, ncol(mat))
}

aeg_test <- data.frame(Genus = aeg_cand)
for (i in seq_along(aeg_cand)) {
  g <- aeg_cand[i]
  a <- ra_vals(cpm_aeg, g); e <- ra_vals(cpm_escc, g); s <- ra_vals(cpm_stad, g)
  aeg_test$p_escc[i]   <- wilcox.test(a, e)$p.value
  aeg_test$p_stad[i]   <- wilcox.test(a, s)$p.value
  aeg_test$med_aeg[i]  <- median(a)
  aeg_test$med_escc[i] <- median(e)
  aeg_test$med_stad[i] <- median(s)
}
aeg_test$padj_escc <- p.adjust(aeg_test$p_escc, "BH")
aeg_test$padj_stad <- p.adjust(aeg_test$p_stad, "BH")
aeg_test$selected  <- with(aeg_test,
  padj_escc < 0.05 & padj_stad < 0.05 & med_aeg > med_escc & med_aeg > med_stad)
write.csv(aeg_test, file.path(DIR_RES, "E_aeg_enriched_stats.csv"), row.names = FALSE)

box_genus <- aeg_test$Genus[aeg_test$selected]          # already in AEG-abundance order
cat("AEG-enriched significant genera (n=", length(box_genus), "):\n", sep = "")
print(box_genus)

abund_box_comp <- do.call(rbind, lapply(box_genus, function(g) rbind(
  data.frame(Genus = g, Cancer = "ESCC", RA = ra_vals(cpm_escc, g)),
  data.frame(Genus = g, Cancer = "AEG",  RA = ra_vals(cpm_aeg,  g)),
  data.frame(Genus = g, Cancer = "STAD", RA = ra_vals(cpm_stad, g))
)))
abund_box_comp$Cancer <- factor(abund_box_comp$Cancer, levels = c("ESCC", "AEG", "STAD"))
abund_box_comp$Genus  <- factor(abund_box_comp$Genus, levels = box_genus)
abund_box_comp$`log2 CPM` <- log2(abund_box_comp$RA * 1e4 + 1)

ncol_facet <- min(length(box_genus), 4)
nrow_facet <- ceiling(length(box_genus) / ncol_facet)
pdf(file.path(DIR_RES, "E_box_aeg_specific_top_genera.pdf"),
    width = ncol_facet * 2 + 1, height = nrow_facet * 2.2 + 0.6)
p_box_genus <- ggboxplot(abund_box_comp, x = "Cancer", y = "log2 CPM",
          col = "Cancer", palette = c("#74ACCE", "#EEC549", "#8B87BA"),
          add = "jitter", facet.by = "Genus", ncol = ncol_facet) +
  stat_compare_means(comparisons = list(c("ESCC", "AEG"), c("AEG", "STAD")),
                     size = 2.4) +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        legend.position = "none")
print(p_box_genus)
dev.off()

##########################
# AEG overall distribution
abund_aeg_distribution <- cpm_aeg / 1e4
abund_aeg_distribution <- abund_aeg_distribution[order(rowMeans(abund_aeg_distribution), 
                                                       decreasing = TRUE), ]
abund_aeg_distribution <- abund_aeg_distribution[1:20, ]

abund_aeg_distribution <- rbind(abund_aeg_distribution, 
                                Others = 100 - apply(abund_aeg_distribution, MARGIN = 2, FUN = sum)) 
abund_aeg_distribution <- abund_aeg_distribution[, order(unlist(abund_aeg_distribution[1,]),
                                                         decreasing = T)]
abund_aeg_distribution$Genus <- rownames(abund_aeg_distribution)

long_distribution <- as.data.frame(reshape2::melt(abund_aeg_distribution, id.vars = c("Genus")))

pdf(file.path(DIR_RES, "G_aeg_distribution_genera_above_top20.pdf"),
    width = 6, height = 4.75)
ggplot(data = long_distribution, aes(x = variable, y = value, 
                                     alluvium = factor(Genus, levels = unique(Genus)))) +
  geom_alluvium(aes(fill = factor(Genus, levels = unique(Genus))),
                alpha = 1) +
  scale_fill_manual(values = colAbund,
                    name = "Genus") +
  xlab("Samples") +
  ylab("Relative abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing.x = element_blank())
dev.off()

###############################################################
# Supplementary: sequencing-depth QC & alpha diversity by tissue
###############################################################

## ---- QC: per-sample bacterial load (reads) & detected genera (by cancer) ----
qc_df <- rbind(
  data.frame(Sample = colnames(count_escc), Cancer = "ESCC",
             Reads = colSums(count_escc), Genera = colSums(count_escc > 0)),
  data.frame(Sample = colnames(count_aeg),  Cancer = "AEG",
             Reads = colSums(count_aeg),  Genera = colSums(count_aeg > 0)),
  data.frame(Sample = colnames(count_stad), Cancer = "STAD",
             Reads = colSums(count_stad), Genera = colSums(count_stad > 0))
)
qc_df$Cancer <- factor(qc_df$Cancer, levels = c("ESCC", "AEG", "STAD"))

cancer_cols <- c(ESCC = "#74ACCE", AEG = "#EEC549", STAD = "#8B87BA")
cancer_comp <- list(c("ESCC", "AEG"), c("AEG", "STAD"), c("ESCC", "STAD"))

# Bacterial load (per-sample bacterial reads). NOTE: raw reads are NOT normalised
# to host DNA input, so this largely reflects sequencing depth (ESCC ~13x shallower),
# not absolute bacterial load.
pdf(file.path(DIR_RES, "S_qc_bacterial_load_among_cancer.pdf"), width = 4, height = 4)
print(
  ggplot(qc_df, aes(x = Cancer, y = Reads, fill = Cancer)) +
    geom_violin(scale = "width", colour = NA, alpha = .8) +
    geom_boxplot(width = .18, fill = "grey90", outlier.colour = NA) +
    geom_jitter(width = .12, size = .2, colour = "grey40", alpha = .5) +
    scale_fill_manual(values = cancer_cols) +
    stat_compare_means(comparisons = cancer_comp, size = 2.4) +
    scale_y_log10() +
    xlab("") + ylab("Bacterial reads / sample") +
    labs(caption = "Raw bacterial reads (not host-normalised); reflects sequencing depth") +
    theme_test() +
    theme(aspect.ratio = 1,
          panel.border = element_rect(fill = NA, colour = 1),
          axis.text = element_text(colour = 1),
          plot.caption = element_text(size = 6, colour = "grey40"),
          legend.position = "none")
)
dev.off()

# Detected genera per sample (violin + box)
pdf(file.path(DIR_RES, "S_qc_detected_genera_among_cancer.pdf"), width = 3.6, height = 3.8)
print(
  ggplot(qc_df, aes(x = Cancer, y = Genera, fill = Cancer)) +
    geom_violin(scale = "width", colour = NA, alpha = .8) +
    geom_boxplot(width = .18, fill = "grey90", outlier.colour = NA) +
    geom_jitter(width = .12, size = .2, colour = "grey40", alpha = .5) +
    scale_fill_manual(values = cancer_cols) +
    stat_compare_means(comparisons = cancer_comp, size = 2.4) +
    xlab("") + ylab("Detected genera / sample") +
    theme_test() +
    theme(aspect.ratio = 1,
          panel.border = element_rect(fill = NA, colour = 1),
          axis.text = element_text(colour = 1),
          legend.position = "none")
)
dev.off()


## ---- Alpha diversity among cancers (pooled): Shannon / Simpson / Chao1 ----
alpha_index_df <- function(count_mat, cancer) {
  tm <- t(count_mat)
  est <- vegan::estimateR(tm)
  data.frame(Sample = colnames(count_mat), Cancer = cancer,
             Shannon = vegan::diversity(tm, "shannon"),
             Simpson = vegan::diversity(tm, "simpson"),
             Chao1   = est["S.chao1", ])
}
alpha_idx <- rbind(alpha_index_df(count_escc, "ESCC"),
                   alpha_index_df(count_aeg,  "AEG"),
                   alpha_index_df(count_stad, "STAD"))
alpha_idx$Cancer <- factor(alpha_idx$Cancer, levels = c("ESCC", "AEG", "STAD"))
alpha_idx_long <- reshape2::melt(alpha_idx, id.vars = c("Sample", "Cancer"),
                                 measure.vars = c("Shannon", "Simpson", "Chao1"),
                                 variable.name = "Metric", value.name = "Value")
alpha_idx_long$Metric <- factor(alpha_idx_long$Metric,
                                levels = c("Shannon", "Simpson", "Chao1"))

pdf(file.path(DIR_RES, "S_alpha_diversity_among_cancer_indices.pdf"), width = 8.5, height = 3.6)
print(
  ggplot(alpha_idx_long, aes(x = Cancer, y = Value, fill = Cancer)) +
    geom_violin(scale = "width", colour = NA, alpha = .8) +
    geom_boxplot(width = .18, fill = "grey90", outlier.colour = NA) +
    geom_jitter(width = .12, size = .2, colour = "grey40", alpha = .4) +
    scale_fill_manual(values = cancer_cols) +
    stat_compare_means(comparisons = cancer_comp, size = 2.2) +
    facet_wrap(~ Metric, scales = "free_y") +
    xlab("") + ylab("") +
    labs(caption = "Chao1 is depth-sensitive; ESCC is ~13x shallower (see QC).") +
    theme_test() +
    theme(aspect.ratio = 1,
          panel.border = element_rect(fill = NA, colour = 1),
          axis.text = element_text(colour = 1),
          plot.caption = element_text(size = 6, colour = "grey40"),
          legend.position = "none")
)
dev.off()
