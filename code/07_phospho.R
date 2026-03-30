####################################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 07. Using phosphoproteomics to understand potential mechanism
#
####################################################################
library(dplyr)
library(stringr)
library(tidyr)
library(broom)
library(purrr)
library(circlize)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(KSEAapp)

set.seed(42)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_FIG <- "/data/yzwang/project/AEG_seiri/results/F4/"
DIR_SUP <- "/data/yzwang/project/AEG_seiri/results/S4/"

##########################
# Preprocess original data
phospho <- read.delim(file.path(DIR_TAB, "Phosphoproteomics_iBAQ103_log2quantile_normlization_impute.txt"))

id_match <- read.delim(file.path(DIR_TAB, "Phospho_IDs_Protein_Gene.txt"))

for (i in 1:nrow(phospho)) {
  pro_id <- strsplit(phospho$Phospho_site[i], "_")[[1]][1]
  gene_name <- id_match[id_match$Protein.accession == pro_id, "Gene.name"]
  site <- strsplit(phospho$Phospho_site[i], "_")[[1]][2]
  phospho$Phospho_site[i] <- paste0(gene_name, "_", site)
}
phospho <- phospho[-grep("^_", phospho$Phospho_site), ]

colnames(phospho) <- sub("^([A-Z]+)(\\d{6})$", "\\10\\2", colnames(phospho))
pair_file <- read.csv(file.path(DIR_TAB, "ms_sample_name_pairing.csv"))
pair_file$InnerSample <- sub("^([A-Z]+)(\\d{6})$", "\\10\\2", pair_file$InnerSample)

phospho <- phospho[ , c(TRUE, colnames(phospho)[-1] %in% pair_file$InnerSample)]
sample_map <- setNames(pair_file$ZipSample, pair_file$InnerSample)
colnames(phospho)[-1] <- sample_map[colnames(phospho)[-1]]

rm_dup <- c(which(phospho$Phospho_site == "TMPO_S184")[1],
            which(phospho$Phospho_site == "TMPO_S306")[1])
phospho <- phospho[-rm_dup, ]
rownames(phospho) <- phospho$Phospho_site
phospho <- phospho[ , -1]

saveRDS(phospho, file.path(DIR_RDS, "Phosphoproteome_human_preprocessed.rds"))

#################
# Formal analysis
phospho <- readRDS(file.path(DIR_RDS, "Phosphoproteome_human_preprocessed.rds"))

mtx_cpm <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
common_sample <- intersect(colnames(phospho), colnames(mtx_cpm))
common_tumour <- common_sample[grep("C", common_sample)]

phos_tumour <- phospho[ , common_tumour] %>%
  filter(rowMeans(.) != min(phospho))
phos_tumour <- phos_tumour[!grepl("^-", rownames(phos_tumour)), ]

abund_ef <- data.frame(
  sample = common_tumour,
  t(mtx_cpm["Enterococcus_faecalis", common_tumour] / 1e+4)
)

phos_detect <- phos_tumour %>%
  as.data.frame() %>%
  mutate(feature = rownames(.)) %>%
  pivot_longer(-feature, names_to = "sample", values_to = "intensity") %>%
  mutate(present = as.integer(intensity > 3.891872))

#####################
# Logistic regression
logistic_ori <- phos_detect %>%
  left_join(abund_ef, by = "sample") %>%
  group_by(feature) %>%
  nest() %>%
  mutate(
    n_present = map_int(data, ~sum(.$present)),
    n_absent  = map_int(data, ~sum(.$present == 0)),
    prop_present = n_present / (n_present + n_absent)
  ) %>%
  filter(prop_present > 0.1, prop_present < 0.9) %>%
  mutate(
    fit = map(data, ~glm(present ~ Enterococcus_faecalis, data = ., family = binomial)),
    converged = map_lgl(fit, ~.$converged),
    tidied = map(fit, tidy)
  ) %>%
  filter(converged) %>%
  unnest(tidied) %>%
  filter(term == "Enterococcus_faecalis") %>%
  mutate(OR = exp(estimate)) %>%
  select(feature, estimate, std.error, OR, p.value, n_present, n_absent)

sig_discrete <- logistic_ori %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  filter(padj < 0.05)

write.csv(sig_discrete, file.path(DIR_TAB, "Phospho_discrete_significance.csv"))

# Landscape: all significant features
circos_data <- sig_discrete %>%
  mutate(
    log2OR = log2(OR),
    neg_log10_padj = -log10(padj),
    total = n_present + n_absent,
    presence_ratio = n_present / total,
    show_label = abs(estimate) > 4
  ) %>%
  arrange(-log2OR)

pdf(file.path(DIR_FIG, "A_Circos_phosphosites_presence_with_ef.pdf"), width = 12, height = 12)
circos.par(start.degree = 90, 
           gap.degree = 360 / nrow(circos_data) * 0.5,
           track.margin = c(0.01, 0.01),
           cell.padding = c(0.02, 0, 0.02, 0))

circos.initialize(factors = circos_data$feature, 
                  xlim = matrix(c(0, 1), ncol = 2, nrow = nrow(circos_data), byrow = TRUE))

# Track 1: log2OR
circos.track(factors = circos_data$feature, 
             y = circos_data$log2OR,
             panel.fun = function(x, y) {
               circos.barplot(value = y, pos = 0.5, 
                              col = ifelse(y > 0, "#E74C3C", "#3498DB"),
                              border = NA)
             },
             bg.border = NA, track.height = 0.3, ylim = range(circos_data$log2OR))

# Track 2: -log10(padj)
col_fun_p <- colorRamp2(c(min(circos_data$neg_log10_padj), 
                          max(circos_data$neg_log10_padj)), 
                        c("white", "darkred"))
circos.track(factors = circos_data$feature,
             y = circos_data$neg_log10_padj,
             panel.fun = function(x, y) {
               circos.rect(0, 0, 1, 1,
                           col = col_fun_p(y), border = NA)
             },
             bg.border = NA, track.height = 0.15, ylim = c(0, 1))

# Track 3: n_present
circos.track(factors = circos_data$feature,
             y = circos_data$n_present,
             panel.fun = function(x, y) {
               circos.barplot(value = y, pos = 0.5,
                              col = "#27AE60", border = NA)
             },
             bg.border = NA, track.height = 0.2, ylim = c(0, max(circos_data$n_present)))

# Track 4: n_absent
circos.track(factors = circos_data$feature,
             y = circos_data$n_absent,
             panel.fun = function(x, y) {
               circos.barplot(value = y, pos = 0.5,
                              col = "#F39C12", border = NA)
             },
             bg.border = NA, track.height = 0.2, ylim = c(0, max(circos_data$n_absent)))

legend("topright", 
       legend = c("log2OR (+)", "log2OR (-)", "-log10(padj)", "n_present", "n_absent"),
       fill = c("#E74C3C", "#3498DB", "darkred", "#27AE60", "#F39C12"),
       bty = "n", cex = 0.7)

circos.clear()
dev.off()

#########################
# Forest plot: top log2OR
n_forest <- 20
sig_discrete <- ungroup(sig_discrete)
top_log2or <- bind_rows(
  slice_max(sig_discrete, estimate, n = n_forest) %>% mutate(direction = "up"),
  slice_min(sig_discrete, estimate, n = n_forest) %>% mutate(direction = "down")
) %>%
  mutate(ci_low  = estimate - 1.96 * std.error,
         ci_high = estimate + 1.96 * std.error,
         feature = reorder(feature, estimate))
dim(top_log2or)

write.csv(top_log2or, file.path(DIR_TAB, "Phospho_top_sites_by_log2OR.csv"))

pdf(file.path(DIR_FIG, "B_Forest_log2OR.pdf"), width = 5, height = 8)
ggplot(top_log2or, aes(x = estimate, y = feature, color = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.3) +
  geom_point(size = 3) +
  scale_color_manual(values = c(up = "tomato2", down = "steelblue3")) +
  labs(x = "log2OR", y = NULL) +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

################
# Slope analysis
summary_sig <- phos_detect %>%
  filter(feature %in% sig_discrete$feature) %>%
  left_join(abund_ef, by = "sample") %>%
  mutate(ef_quartile = cut(Enterococcus_faecalis,
                           breaks = quantile(Enterococcus_faecalis, probs = 0:4/4, na.rm = TRUE),
                           include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4"))) %>%
  group_by(feature, ef_quartile) %>%
  summarise(presence_pct = mean(present, na.rm = TRUE), .groups = "drop")

mat_sig <- summary_sig %>%
  pivot_wider(names_from = ef_quartile, values_from = presence_pct) %>%
  { m <- as.matrix(.[, c("Q1","Q2","Q3","Q4")]); rownames(m) <- .$feature; m }

slope_df <- as.data.frame(t(apply(mat_sig, 1, function(y) {
  fit <- lm(y ~ c(1,2,3,4))
  c(slope = coef(fit)[2], pval = summary(fit)$coefficients[2, 4])
})))

quantile(abs(sig_discrete$estimate))
log2or_cutoff <- 2

slope_sig <- slope_df %>%
  filter(pval < 0.1) %>%
  rownames_to_column("feature") %>%
  rename(slope = 2) %>%
  inner_join(sig_discrete %>% select(feature, estimate, padj), by = "feature") %>%
  filter(abs(estimate) > log2or_cutoff)
dim(slope_sig)

###############################
# Case visualisation: top slope
top_features <- top_log2or %>%
  filter(feature %in% slope_sig$feature) %>%
  left_join(slope_sig %>% select(feature, slope, pval), by = "feature") %>%
  mutate(direction = ifelse(estimate > 0, "positive", "negative"))

feat_order <- top_features$feature[order(top_features$estimate)]

plot_percent <- phos_detect %>%
  filter(feature %in% top_features$feature) %>%
  left_join(abund_ef, by = "sample") %>%
  mutate(ef_quartile = cut(Enterococcus_faecalis,
                           breaks = quantile(Enterococcus_faecalis, probs = 0:4/4, na.rm = TRUE),
                           include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4")))

summary_data <- plot_percent %>%
  group_by(feature, ef_quartile) %>%
  summarise(presence_pct = mean(present, na.rm = TRUE), .groups = "drop") %>%
  left_join(top_features %>% select(feature, direction), by = "feature") %>%
  mutate(fill_var = paste0(direction, "_", ef_quartile))

annotation_data <- top_features %>%
  mutate(label = paste0("log2OR = ", round(estimate, 2),
                        "\nslope = ",  round(slope, 4),
                        "\np.adj = ",  format(padj, digits = 2, scientific = TRUE)),
         feature = factor(feature, levels = feat_order))

fill_colors <- c(
  setNames(brewer.pal(4, "Oranges"), paste0("positive_Q", 1:4)),
  setNames(brewer.pal(4, "Blues"),   paste0("negative_Q", 1:4))
)

plot_percent$feature <- factor(plot_percent$feature, levels = feat_order)
summary_data$feature <- factor(summary_data$feature, levels = feat_order)

facet_wrap(~ feature, scales = "free_y", ncol = 5)
n_features <- length(unique(plot_percent$feature))
n_col <- 5

n_row <- ceiling(n_features / n_col)

pdf(file.path(DIR_SUP, "A_Top_slope_features.pdf"), width = n_col * 2, height = n_row * 2)
ggplot() +
  geom_jitter(data = plot_percent,
              aes(x = ef_quartile, y = present),
              alpha = 0.3, color = "grey50", size = 1, height = 0.02, width = 0.1) +
  geom_line(data = summary_data,
            aes(x = ef_quartile, y = presence_pct, group = feature),
            color = "black", linewidth = 0.8) +
  geom_point(data = summary_data,
             aes(x = ef_quartile, y = presence_pct, fill = fill_var),
             color = "black", shape = 21, size = 4) +
  geom_text(data = annotation_data,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.2, size = 3, lineheight = 0.8) +
  scale_fill_manual(values = fill_colors) +
  facet_wrap(~ feature, scales = "free_y", ncol = n_col) +
  labs(x = "E.faecalis abundance (quartiles)", y = "Presence") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        strip.text = element_text(size = 9),
        legend.position = "none")
dev.off()

#################
# KSEA enrichment
ksea_input <- logistic_ori %>%
  ungroup() %>%
  mutate(
    Protein      = sub("_[STY]\\d+$", "", feature),
    Gene         = Protein,
    Peptide      = feature,
    Residue.Both = sub("^[^_]+_", "", feature),
    p            = p.value,
    FC           = abs(OR)
  ) %>%
  select(Protein, Gene, Peptide, Residue.Both, p, FC) %>%
  as.data.frame()

ksea_scores <- KSEA.Scores(
  KSData           = KSData,
  PX               = ksea_input,
  NetworKIN        = TRUE,
  NetworKIN.cutoff = 5
)

ksea_plot <- ksea_scores %>%
  filter(!is.na(p.value)) %>%
  mutate(
    sig = p.value < 0.05,
    direction = ifelse(z.score > 0, "activated", "inhibited"),
    Kinase = reorder(Kinase.Gene, z.score)
  ) %>%
  filter(sig)

pdf(file.path(DIR_FIG, "C_KSEA_kinase_activity.pdf"), width = 4, height = 0.2 * nrow(ksea_plot) + 1)
ggplot(ksea_plot, aes(x = z.score, y = Kinase, fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, color = "grey30") +
  scale_fill_manual(values = c(activated = "tomato2", inhibited = "steelblue3")) +
  labs(x = "KSEA z-score", y = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 9),
        panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))
dev.off()
