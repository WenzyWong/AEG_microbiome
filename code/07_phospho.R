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
library(ggridges)
library(ComplexHeatmap)
library(patchwork)

set.seed(42)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_FIG <- "/data/yzwang/project/AEG_seiri/results/F4/"
DIR_SUP <- "/data/yzwang/project/AEG_seiri/results/S4/"

##########################
# Preprocess original data
phospho <-
  read.delim(
    file.path(
      DIR_TAB,
      "Phosphoproteomics_iBAQ103_log2quantile_normlization_impute.txt"
    )
  )

id_match <-
  read.delim(file.path(DIR_TAB, "Phospho_IDs_Protein_Gene.txt"))

for (i in 1:nrow(phospho)) {
  pro_id <- strsplit(phospho$Phospho_site[i], "_")[[1]][1]
  gene_name <-
    id_match[id_match$Protein.accession == pro_id, "Gene.name"]
  site <- strsplit(phospho$Phospho_site[i], "_")[[1]][2]
  phospho$Phospho_site[i] <- paste0(gene_name, "_", site)
}
phospho <- phospho[-grep("^_", phospho$Phospho_site),]

colnames(phospho) <-
  sub("^([A-Z]+)(\\d{6})$", "\\10\\2", colnames(phospho))
pair_file <-
  read.csv(file.path(DIR_TAB, "ms_sample_name_pairing.csv"))
pair_file$InnerSample <-
  sub("^([A-Z]+)(\\d{6})$", "\\10\\2", pair_file$InnerSample)

phospho <-
  phospho[, c(TRUE, colnames(phospho)[-1] %in% pair_file$InnerSample)]
sample_map <- setNames(pair_file$ZipSample, pair_file$InnerSample)
colnames(phospho)[-1] <- sample_map[colnames(phospho)[-1]]

rm_dup <- c(
  which(phospho$Phospho_site == "TMPO_S184")[1],
  which(phospho$Phospho_site == "TMPO_S306")[1]
)
phospho <- phospho[-rm_dup,]
rownames(phospho) <- phospho$Phospho_site
phospho <- phospho[,-1]

saveRDS(phospho,
        file.path(DIR_RDS, "Phosphoproteome_human_preprocessed.rds"))

#################
# Formal analysis
phospho <- readRDS(file.path(DIR_RDS, "Phosphoproteome_human_preprocessed.rds"))

mtx_cpm <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
common_sample <- intersect(colnames(phospho), colnames(mtx_cpm))
common_tumour <- common_sample[grep("C", common_sample)]

phos_tumour <- phospho[, common_tumour] %>%
  filter(rowMeans(.) != min(phospho))
phos_tumour <- phos_tumour[!grepl("^-", rownames(phos_tumour)),]

abund_ef <- data.frame(sample = common_tumour,
                       t(mtx_cpm["Enterococcus_faecalis", common_tumour] / 1e+4))

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
    n_present = map_int(data, ~ sum(.$present)),
    n_absent  = map_int(data, ~ sum(.$present == 0)),
    prop_present = n_present / (n_present + n_absent)
  ) %>%
  filter(prop_present > 0.1, prop_present < 0.9) %>%
  mutate(
    fit = map(
      data,
      ~ glm(
        present ~ Enterococcus_faecalis,
        data = .,
        family = binomial
      )
    ),
    converged = map_lgl(fit, ~ .$converged),
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

pdf(
  file.path(DIR_FIG, "A_Circos_phosphosites_presence_with_ef.pdf"),
  width = 12,
  height = 12
)
circos.par(
  start.degree = 90,
  gap.degree = 360 / nrow(circos_data) * 0.5,
  track.margin = c(0.01, 0.01),
  cell.padding = c(0.02, 0, 0.02, 0)
)

circos.initialize(factors = circos_data$feature,
                  xlim = matrix(
                    c(0, 1),
                    ncol = 2,
                    nrow = nrow(circos_data),
                    byrow = TRUE
                  ))

# Track 1: log2OR
circos.track(
  factors = circos_data$feature,
  y = circos_data$log2OR,
  panel.fun = function(x, y) {
    circos.barplot(
      value = y,
      pos = 0.5,
      col = ifelse(y > 0, "#E74C3C", "#3498DB"),
      border = NA
    )
  },
  bg.border = NA,
  track.height = 0.3,
  ylim = range(circos_data$log2OR)
)

# Track 2: -log10(padj)
col_fun_p <- colorRamp2(c(
  min(circos_data$neg_log10_padj),
  max(circos_data$neg_log10_padj)
),
c("white", "darkred"))
circos.track(
  factors = circos_data$feature,
  y = circos_data$neg_log10_padj,
  panel.fun = function(x, y) {
    circos.rect(0, 0, 1, 1,
                col = col_fun_p(y), border = NA)
  },
  bg.border = NA,
  track.height = 0.15,
  ylim = c(0, 1)
)

# Track 3: n_present
circos.track(
  factors = circos_data$feature,
  y = circos_data$n_present,
  panel.fun = function(x, y) {
    circos.barplot(
      value = y,
      pos = 0.5,
      col = "#27AE60",
      border = NA
    )
  },
  bg.border = NA,
  track.height = 0.2,
  ylim = c(0, max(circos_data$n_present))
)

# Track 4: n_absent
circos.track(
  factors = circos_data$feature,
  y = circos_data$n_absent,
  panel.fun = function(x, y) {
    circos.barplot(
      value = y,
      pos = 0.5,
      col = "#F39C12",
      border = NA
    )
  },
  bg.border = NA,
  track.height = 0.2,
  ylim = c(0, max(circos_data$n_absent))
)

legend(
  "topright",
  legend = c(
    "log2OR (+)",
    "log2OR (-)",
    "-log10(padj)",
    "n_present",
    "n_absent"
  ),
  fill = c("#E74C3C", "#3498DB", "darkred", "#27AE60", "#F39C12"),
  bty = "n",
  cex = 0.7
)

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
  mutate(
    ci_low  = estimate - 1.96 * std.error,
    ci_high = estimate + 1.96 * std.error,
    feature = reorder(feature, estimate)
  )
dim(top_log2or)

write.csv(top_log2or, file.path(DIR_TAB, "Phospho_top_sites_by_log2OR.csv"))

pdf(file.path(DIR_FIG, "B_Forest_log2OR.pdf"), width = 5, height = 8)
ggplot(top_log2or, aes(x = estimate, y = feature, color = direction)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +
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
  mutate(
    ef_quartile = cut(
      Enterococcus_faecalis,
      breaks = quantile(
        Enterococcus_faecalis,
        probs = 0:4 / 4,
        na.rm = TRUE
      ),
      include.lowest = TRUE,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )
  ) %>%
  group_by(feature, ef_quartile) %>%
  summarise(presence_pct = mean(present, na.rm = TRUE),
            .groups = "drop")

mat_sig <- summary_sig %>%
  pivot_wider(names_from = ef_quartile, values_from = presence_pct) %>%
  {
    m <-
      as.matrix(.[, c("Q1", "Q2", "Q3", "Q4")])
    rownames(m) <- .$feature
    m
  }

slope_df <- as.data.frame(t(apply(mat_sig, 1, function(y) {
  fit <- lm(y ~ c(1, 2, 3, 4))
  c(slope = coef(fit)[2],
    pval = summary(fit)$coefficients[2, 4])
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
  mutate(
    ef_quartile = cut(
      Enterococcus_faecalis,
      breaks = quantile(
        Enterococcus_faecalis,
        probs = 0:4 / 4,
        na.rm = TRUE
      ),
      include.lowest = TRUE,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )
  )

summary_data <- plot_percent %>%
  group_by(feature, ef_quartile) %>%
  summarise(presence_pct = mean(present, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(top_features %>% select(feature, direction), by = "feature") %>%
  mutate(fill_var = paste0(direction, "_", ef_quartile))

annotation_data <- top_features %>%
  mutate(
    label = paste0("log2OR = ", round(estimate, 2),
      "\nslope = ", round(slope, 4),
      "\np.adj = ", format(padj, digits = 2, scientific = TRUE)
    ),
    feature = factor(feature, levels = feat_order)
  )

fill_colors <- c(setNames(brewer.pal(4, "Oranges"), paste0("positive_Q", 1:4)),
                 setNames(brewer.pal(4, "Blues"),   paste0("negative_Q", 1:4)))

plot_percent$feature <- factor(plot_percent$feature, levels = feat_order)
summary_data$feature <- factor(summary_data$feature, levels = feat_order)

facet_wrap( ~ feature, scales = "free_y", ncol = 5)
n_features <- length(unique(plot_percent$feature))
n_col <- 5
n_row <- ceiling(n_features / n_col)

pdf(file.path(DIR_SUP, "A_Top_slope_features.pdf"), 
    width = n_col * 2, height = n_row * 2)
ggplot() +
  geom_jitter(
    data = plot_percent,
    aes(x = ef_quartile, y = present),
    alpha = 0.3,
    color = "grey50",
    size = 1,
    height = 0.02,
    width = 0.1
  ) +
  geom_line(
    data = summary_data,
    aes(x = ef_quartile, y = presence_pct, group = feature),
    color = "black",
    linewidth = 0.8
  ) +
  geom_point(
    data = summary_data,
    aes(x = ef_quartile, y = presence_pct, fill = fill_var),
    color = "black",
    shape = 21,
    size = 4
  ) +
  geom_text(
    data = annotation_data,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.05,
    vjust = 1.2,
    size = 3,
    lineheight = 0.8
  ) +
  scale_fill_manual(values = fill_colors) +
  facet_wrap( ~ feature, scales = "free_y", ncol = n_col) +
  labs(x = "E.faecalis abundance (quartiles)", y = "Presence") +
  theme_minimal() +
  theme(
    panel.border    = element_rect(fill = NA, colour = 1),
    axis.text       = element_text(colour = 1),
    strip.text      = element_text(size = 9),
    legend.position = "none"
  )
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
    sig       = p.value < 0.05,
    direction = ifelse(z.score > 0, "activated", "inhibited"),
    Kinase    = reorder(Kinase.Gene, z.score)
  ) %>%
  filter(sig)

pdf(
  file.path(DIR_FIG, "C_KSEA_kinase_activity.pdf"),
  width = 4,
  height = 0.2 * nrow(ksea_plot) + 1
)
ggplot(ksea_plot, aes(x = z.score, y = Kinase, fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, color = "grey30") +
  scale_fill_manual(values = c(activated = "tomato2",inhibited = "steelblue3")) +
  labs(x = "KSEA z-score", y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y     = element_text(size = 9),
    panel.border    = element_rect(fill = NA, colour = 1),
    axis.text       = element_text(colour = 1)
  )
dev.off()

####################
# Clinical relavance
clinical <-
  readxl::read_excel(file.path(DIR_TAB, "AEG_clinical.xlsx"))

# Presence matrix
phos_mat <- phos_detect %>%
  filter(feature %in% sig_discrete$feature) %>%
  select(feature, sample, present) %>%
  pivot_wider(names_from = feature, values_from = present) %>%
  column_to_rownames("sample") %>%
  as.matrix()

clinical_aligned <- clinical %>%
  mutate(sample = paste0("C", No.)) %>%
  filter(sample %in% rownames(phos_mat)) %>%
  arrange(match(sample, rownames(phos_mat)))

phos_mat <- phos_mat[clinical_aligned$sample,]

# Encode clinical variables
fix_stage <-
  function(x) {
    x <- as.character(x)
    x[x == "\u2163"] <- "IV"
    x
  }

parse_tnm <- function(x) {
  list(T = as.integer(gsub(
    ".*T(\\d+).*", "\\1", gsub("[ab]", "", x)
  )),
  N = as.integer(gsub(
    ".*N(\\d+).*", "\\1", gsub("[ab]", "", x)
  )),
  M = as.integer(gsub(".*M(\\d+).*", "\\1", x)))
}

tnm <- parse_tnm(clinical_aligned$`TNM stage`)

anno_df <- clinical_aligned %>%
  transmute(
    Path_stage = as.integer(factor(
      fix_stage(`Pathological stage`),
      levels = c("IB", "IIA", "IIB", "IIIA", "IIIB", "IIIC", "IV")
    )),
    T_stage         = tnm$T,
    N_stage         = tnm$N,
    M_stage         = tnm$M,
    Survival_month  = `Survival time（Month）`,
    Status          = `Status（1=Dead, 0=Alive）`,
    Lauren          = as.integer(factor(`Lauren Classification`)),
    Differentiation = as.integer(factor(`Differentiated degree`)),
    Histology       = as.integer(factor(`Histological classification`)),
    Sex             = ifelse(Sex == "Male", 1, 0),
    Age             = Age,
    Smoking         = ifelse(Smoking == "Yes", 1, 0),
    Alcohol         = ifelse(Alcohol == "Yes", 1, 0)
  ) %>%
  as.data.frame()
rownames(anno_df) <- clinical_aligned$sample

# Spearman correlation
cor_mat <- matrix(
  NA,
  nrow = ncol(phos_mat),
  ncol = ncol(anno_df),
  dimnames = list(colnames(phos_mat), colnames(anno_df))
)
pval_mat <- cor_mat

for (f in colnames(phos_mat)) {
  for (v in colnames(anno_df)) {
    idx <- complete.cases(phos_mat[, f], anno_df[, v])
    if (sum(idx) >= 5) {
      res           <- cor.test(phos_mat[idx, f],
                                anno_df[idx, v],
                                method = "spearman",
                                exact = FALSE)
      cor_mat[f, v]  <- res$estimate
      pval_mat[f, v] <- res$p.value
    }
  }
}

# Filter features
p_thresh  <- 0.01
or_thresh <- 2

log2or_vec <- sig_discrete$estimate[match(rownames(pval_mat), sig_discrete$feature)]

keep_rows <- apply(pval_mat, 1, function(x) any(x < p_thresh, na.rm = TRUE)) &
  abs(log2or_vec) > or_thresh

cor_plot  <- cor_mat[keep_rows,]
pval_plot <- pval_mat[keep_rows,]
sig_mark  <- ifelse(pval_plot < 0.001, "***",
                    ifelse(pval_plot < 0.01,  "**",
                           ifelse(pval_plot < 0.05,  "*", "")))
nrow(cor_plot)

# Annotations
row_anno <- sig_discrete %>%
  filter(feature %in% rownames(cor_plot)) %>%
  select(feature, estimate) %>%
  column_to_rownames("feature")
row_anno <- row_anno[rownames(cor_plot), , drop = FALSE]

row_ha <- rowAnnotation(
  log2OR = row_anno$estimate,
  col    = list(log2OR = colorRamp2(
    c(min(row_anno$estimate), 0, max(row_anno$estimate)),
    c("steelblue3", "white", "tomato2")
  )),
  annotation_name_side = "top",
  width = unit(8, "mm")
)

col_type <- c(
  Path_stage      = "Stage",
  T_stage         = "Stage",
  N_stage         = "Stage",
  M_stage         = "Stage",
  Survival_month  = "Survival",
  Status          = "Survival",
  Lauren          = "Pathology",
  Differentiation = "Pathology",
  Histology       = "Pathology",
  Sex             = "Demographics",
  Age             = "Demographics",
  Smoking         = "Demographics",
  Alcohol         = "Demographics"
)

col_ha <- HeatmapAnnotation(
  Category = col_type[colnames(cor_plot)],
  col = list(
    Category = c(
      Stage        = "#E8A838",
      Survival     = "#3CB371",
      Pathology    = "#9370DB",
      Demographics = "#4682B4"
    )
  ),
  annotation_name_side = "left",
  annotation_height = unit(4, "mm")
)

# Draw heatmap
n_rows <- nrow(cor_plot)

pdf(file.path(DIR_FIG, "D_Clinical_correlation_heatmap.pdf"),
    width = 9, height = 10)
Heatmap(cor_plot, name = "Spearman r",
        col = colorRamp2(c(-0.5, 0, 0.5), rev(brewer.pal(11, "RdBu"))[c(2, 6, 10)]),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (!is.na(sig_mark[i, j]) && sig_mark[i, j] != "")
            grid.text(sig_mark[i, j], x, y, gp = gpar(fontsize = 8))
          },
  top_annotation       = col_ha,
  right_annotation     = row_ha,
  cluster_rows         = TRUE,
  cluster_columns      = FALSE,
  row_names_gp         = gpar(fontsize = max(5, min(9, 120 / n_rows))),
  column_names_gp      = gpar(fontsize = 9),
  column_names_rot     = 45,
  border               = TRUE,
  heatmap_legend_param = list(direction = "vertical")
)
dev.off()

# Revise heatmap to blocks
cor_t <- t(cor_plot)
pval_t <- t(pval_plot)

col_order <- sig_discrete %>%
  filter(feature %in% rownames(cor_plot)) %>%
  arrange(desc(estimate)) %>%
  pull(feature) %>%
  as.character()
row_order <- colnames(cor_plot)

# Clinical variable grouping
var_group <- c(
  Path_stage = "Stage", T_stage = "Stage", N_stage = "Stage", M_stage = "Stage",
  Survival_month = "Survival", Status = "Survival",
  Lauren = "Pathology", Differentiation = "Pathology", Histology = "Pathology",
  Sex = "Demographics", Age = "Demographics",
  Smoking = "Demographics", Alcohol = "Demographics"
)

group_colors <- c(Stage        = "#E8A838",
                  Survival     = "#3CB371",
                  Pathology    = "#9370DB",
                  Demographics = "#4682B4")

plot_df <- as.data.frame(cor_t) %>%
  rownames_to_column("clinical_var") %>%
  pivot_longer(-clinical_var, names_to = "feature", values_to = "r") %>%
  left_join(
    as.data.frame(pval_t) %>%
      rownames_to_column("clinical_var") %>%
      pivot_longer(-clinical_var, names_to = "feature", values_to = "pval"),
    by = c("clinical_var", "feature")
  ) %>%
  mutate(
    clinical_var = factor(clinical_var, levels = row_order),
    feature      = factor(feature,      levels = col_order),
    group        = var_group[clinical_var],
    group        = factor(group, levels = c("Stage", "Survival", "Pathology", "Demographics")),
    # point size scaled by -log10(p), capped for visual clarity
    neg_log_p    = pmin(-log10(pval), 4),
    sig          = pval < 0.05
  )

pgn_pal <- colorRampPalette(rev(brewer.pal(11, "PRGn")))(100)

p_main <- ggplot(plot_df, aes(x = feature, y = clinical_var)) +
  geom_point(aes(fill = r, size = neg_log_p),
             shape = 22, color = "grey30", stroke = 0.2) +
  scale_fill_gradientn(
    colors = pgn_pal,
    limits = c(-0.5, 0.5),
    name   = "Spearman r"
  ) +
  scale_size_continuous(
    range  = c(0.5, 6),
    limits = c(0, 4),
    breaks = c(-log10(0.05), -log10(0.01), -log10(0.001)),
    labels = c("0.05", "0.01", "0.001"),
    name   = "p value"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text.x       = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y       = element_text(size = 8),
    strip.text.y.left = element_text(angle = 0, size = 8, face = "bold"),
    strip.placement   = "outside",
    panel.spacing     = unit(2, "mm"),
    panel.grid        = element_blank(),
    legend.position   = "right"
  )

# log2OR annotation
log2or_df <- sig_discrete %>%
  filter(feature %in% col_order) %>%
  select(feature, estimate) %>%
  mutate(feature = factor(feature, levels = col_order))
p_anno <- ggplot(log2or_df, aes(x = feature, y = 1, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colors = rev(brewer.pal(11, "RdBu"))[c(2, 6, 10)],
    limits = c(-max(abs(log2or_df$estimate)), max(abs(log2or_df$estimate))),
    name   = "log2OR"
  ) +
  scale_x_discrete(limits = col_order) +
  labs(x = NULL, y = "log2OR") +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "right",
    plot.margin  = margin(0, 0, 0, 0)
  )

# calculate presence proportion in Q1 and Q4 for each feature
ef_breaks <- quantile(abund_ef$Enterococcus_faecalis, probs = 0:4/4, na.rm = TRUE)

presence_q <- phos_detect %>%
  filter(feature %in% col_order) %>%
  left_join(abund_ef, by = "sample") %>%
  mutate(ef_quartile = cut(Enterococcus_faecalis,
                           breaks = ef_breaks,
                           include.lowest = TRUE,
                           labels = c("Q1","Q2","Q3","Q4"))) %>%
  filter(ef_quartile %in% c("Q1", "Q4")) %>%
  group_by(feature, ef_quartile) %>%
  summarise(presence_pct = mean(present, na.rm = TRUE), .groups = "drop") %>%
  mutate(feature = factor(feature, levels = col_order))

p_presence <- ggplot(presence_q, aes(x = feature, y = presence_pct,
                                     color = ef_quartile)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c(Q1 = "steelblue3", Q4 = "tomato2"),
                     name = "EF quartile") +
  scale_x_discrete(limits = col_order) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0", "0.5", "1")) +
  labs(x = NULL, y = "Presence(%)") +
  theme_bw() +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.text.y        = element_text(size = 7),
    axis.title.y       = element_text(size = 11, angle = 90),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.key.size    = unit(3, "mm"),
    legend.text        = element_text(size = 7),
    plot.margin        = margin(0, 0, 0, 0)
  )

# Combine all three panels
p_combined <- p_anno / p_presence / p_main +
  plot_layout(heights = c(1, 3, 20), guides = "collect")
p_combined

n_features <- length(col_order)
n_vars     <- length(row_order)

pdf(file.path(DIR_FIG, "E_Clinical_correlation_heatmap.pdf"),
    width  = n_features * 0.18 + 2,
    height = n_vars     * 0.35 + 2)
print(p_combined)
dev.off()
