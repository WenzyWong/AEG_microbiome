#######################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 05. Analysis on microbe-host crosstalk,
# based on transcriptomics and proteomics
#
#######################################################
setwd("/data/yzwang/project/AEG_seiri/")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(compositions)
library(limma)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(paletteer)
library(msigdbr)
library(rEDM)
library(GENIE3)
library(foreach)
library(doParallel)
library(fgsea)

DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/F3_crosstalk/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"

#######################
# Import and preprocess
clinical <- readxl::read_excel(file.path(DIR_TAB, "AEG_clinical.xlsx"))

mtx_htpm_filt <- read.csv(file.path(DIR_TAB, "gene_tpm_matrix.csv"), row.names = 1) %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_name = gsub(".*\\|", "", gene_id)) %>%
  select(-gene_id) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum)) %>%
  filter(!grepl("^<class", gene_name)) %>%
  column_to_rownames("gene_name") %>%
  filter(rowMeans(.) > 1)

mtx_cpm    <- readRDS(file.path(DIR_RDS, "gAEG_CPM_RNA_FiltMyco.rds"))
mtx_cpm_sp <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))

mtx_hpro <- readxl::read_excel(file.path(DIR_TAB, "SupData13ProteinIntensity.xlsx"))
colnames(mtx_hpro)[3:ncol(mtx_hpro)] <- {
  raw <- colnames(mtx_hpro)[3:ncol(mtx_hpro)]
  ifelse(grepl("_T$", raw),
         paste0("C", gsub("AEG|_T$", "", raw)),
         paste0("N", gsub("AEG|_N$", "", raw)))
}
rownames(mtx_hpro) <- mtx_hpro$Protein_group

samples_common <- Reduce(intersect, list(
  grep("^C", colnames(mtx_htpm_filt), value = TRUE),
  grep("^C", colnames(mtx_hpro),      value = TRUE),
  grep("^C", colnames(mtx_cpm),       value = TRUE)
))

#######################
# Top 20 genera
top20_genera <- rownames(mtx_cpm)[order(rowMeans(mtx_cpm / 1e4), decreasing = TRUE)][1:20]
mean_abund   <- rowMeans(mtx_cpm[top20_genera, ] / 1e4)

mtx_clr <- t(as.matrix(clr(t(mtx_cpm[top20_genera, samples_common] + 0.5))))

#######################
# Shared helpers
spearman_pval <- function(r_mat, n) {
  2 * pt(-abs(r_mat * sqrt((n - 2) / (1 - r_mat^2))), df = n - 2)
}

# Protein correlation helper (reused for genus and species)
calc_prot_cor <- function(clr_mat, prot_mat, samples) {
  out <- matrix(NA_real_, nrow = nrow(clr_mat), ncol = nrow(prot_mat),
                dimnames = list(rownames(clr_mat), rownames(prot_mat)))
  for (g in seq_len(nrow(clr_mat))) {
    abund_vec <- clr_mat[g, samples]
    out[g, ] <- apply(prot_mat, 1, function(v) {
      idx <- !is.na(v)
      if (sum(idx) < 10) return(NA_real_)
      cor(abund_vec[idx], v[idx], method = "spearman")
    })
  }
  out
}

#######################
# RNA-seq correlation (genus)
mtx_tpm_cor <- log1p(as.matrix(mtx_htpm_filt[, samples_common]))
cv_gene     <- apply(mtx_tpm_cor, 1, function(x) sd(x) / (mean(x) + 1e-6))
mtx_tpm_cor <- mtx_tpm_cor[cv_gene > 0.2, ]

rna_cor_mat <- do.call(cbind, lapply(
  split(seq_len(nrow(mtx_tpm_cor)), ceiling(seq_len(nrow(mtx_tpm_cor)) / 2000)),
  function(idx) cor(t(mtx_clr), t(mtx_tpm_cor[idx, , drop = FALSE]),
                    method = "spearman", use = "pairwise.complete.obs")
))
saveRDS(rna_cor_mat, file.path(DIR_RDS, "rna_genus_gene_cor.rds"))

#######################
# Proteomics preprocessing (shared)
mtx_prot_num <- as.matrix(mutate(mtx_hpro[, samples_common], across(everything(), as.numeric)))
rownames(mtx_prot_num) <- rownames(mtx_hpro)
mtx_prot_num[mtx_prot_num == 0] <- NA
mtx_prot_log <- log2(mtx_prot_num)
mtx_prot_log <- mtx_prot_log[rowSums(!is.na(mtx_prot_log)) >= ceiling(0.5 * length(samples_common)), ]

#######################
# Proteomics correlation (genus)
prot_cor_mat <- calc_prot_cor(mtx_clr, mtx_prot_log, samples_common)
saveRDS(prot_cor_mat, file.path(DIR_RDS, "prot_genus_protein_cor.rds"))

#######################
# Circos lollipop
n_common <- length(samples_common)
rna_sig  <- abs(rna_cor_mat) >= 0.3 & spearman_pval(rna_cor_mat,  n_common) < 0.05
prot_sig <- abs(prot_cor_mat) >= 0.3 & spearman_pval(prot_cor_mat, n_common) < 0.05
rna_sig[is.na(rna_sig)]   <- FALSE
prot_sig[is.na(prot_sig)] <- FALSE

prot_gene_map <- setNames(mtx_hpro$Gene_names, mtx_hpro$Protein_group)

overlap_df <- data.frame(
  genus      = rownames(rna_cor_mat),
  count      = sapply(rownames(rna_cor_mat), function(g) {
    rna_genes  <- colnames(rna_cor_mat)[rna_sig[g, ]]
    prot_genes <- colnames(prot_cor_mat)[prot_sig[g, ]] %>%
      { prot_gene_map[.] } %>% na.omit() %>%
      paste(collapse = ";") %>% strsplit(";") %>%
      unlist() %>% trimws() %>% unique()
    length(intersect(rna_genes, prot_genes))
  }),
  prot_total = rowSums(prot_sig),
  stringsAsFactors = FALSE
) %>%
  mutate(pct = ifelse(prot_total > 0, count / prot_total * 100, 0)) %>%
  .[match(top20_genera, .$genus), ]

genera <- overlap_df$genus
abund_scaled <- setNames(
  (mean_abund - min(mean_abund)) / (max(mean_abund) - min(mean_abund)),
  top20_genera
)
genus_colors <- setNames(
  as.character(paletteer_d("khroma::discreterainbow")[c(10, 12:20, 23:27, 2, 4, 5, 7, 9)]),
  top20_genera
)
max_pct        <- ceiling(max(overlap_df$pct))
high_overlap_g <- overlap_df$genus[overlap_df$pct > 10]

pdf(file.path(DIR_RES, "Circos_genus_overlap_lollipop.pdf"), width = 10, height = 10)
circos.clear()
circos.par(start.degree = 90, gap.degree = 2,
           cell.padding = c(0, 0, 0, 0),
           points.overflow.warning = FALSE,
           track.margin = c(0.005, 0.005))

circos.initialize(factors = factor(genera, levels = genera),
                  xlim = matrix(rep(c(0, 1), length(genera)), ncol = 2, byrow = TRUE))

# Track 1: lollipop
circos.track(
  factors = factor(genera, levels = genera),
  ylim = c(0, max_pct), track.height = 0.38,
  bg.border = NA, bg.col = NA,
  panel.fun = function(x, y) {
    g   <- get.cell.meta.data("sector.index")
    row <- overlap_df[overlap_df$genus == g, ]
    is_high <- g %in% high_overlap_g
    circos.rect(
      get.cell.meta.data("xlim")[1], 0,
      get.cell.meta.data("xlim")[2], max_pct,
      col = NA, border = ifelse(is_high, "#c0392b", "grey85"),
      lwd = ifelse(is_high, 2.5, 1)
    )
    circos.segments(0.5, 0, 0.5, row$pct, lwd = 2.2, col = genus_colors[g])
    if (row$pct > 0)
      circos.points(0.5, row$pct, pch = 16,
                    cex = 1.5 + row$pct / max_pct * 1.1, col = genus_colors[g])
    circos.text(0.5, row$pct,
                labels = sprintf("%.1f%%\n(n=%d)", row$pct, row$count),
                adj = c(0.5, -0.3), cex = 0.45,
                niceFacing = TRUE, facing = "bending.inside", col = "black")
  }
)
circos.yaxis(side = "left", sector.index = genera[1],
             labels.cex = 0.42, at = pretty(c(0, max_pct), n = 4))

# Track 2: genus label
circos.track(
  factors = factor(genera, levels = genera),
  ylim = c(0, 1), track.height = 0.08, bg.border = NA,
  panel.fun = function(x, y) {
    g <- get.cell.meta.data("sector.index")
    circos.rect(0, 0, 1, 1, col = genus_colors[g], border = "white", lwd = 0.5)
    circos.text(0.5, 0.5, labels = gsub("_.*", "", g),
                facing = "bending.inside", niceFacing = TRUE,
                cex = 0.45, font = 2, col = "white")
  }
)

# Track 3: abundance
abund_col_fun <- colorRamp2(c(0, 1), c("#e8f4f8", "#1a5276"))
circos.track(
  factors = factor(genera, levels = genera),
  ylim = c(0, 1), track.height = 0.12,
  bg.border = "grey85", bg.col = "grey97",
  panel.fun = function(x, y) {
    g   <- get.cell.meta.data("sector.index")
    val <- abund_scaled[g]
    circos.rect(0, 0, 1, 1, col = abund_col_fun(val), border = NA)
    circos.text(0.5, 0.5,
                labels     = formatC(mean_abund[g], format = "f", digits = 2),
                adj        = c(0.5, 0.5), cex = 0.38,
                niceFacing = TRUE, facing = "bending.inside",
                col        = ifelse(val > 0.6, "white", "black"))
  }
)

circos.text(0, 0, labels = "Abund.\n(%)",
            sector.index = genera[1], track.index = 3,
            facing = "bending.inside", cex = 0.4, col = "black")

legend("center", legend = gsub("_", " ", genera), fill = genus_colors[genera],
       border = NA, cex = 0.55, ncol = 2, title = "Genus", bty = "n")
title(main = expression("Overlap genes (|" * italic(r)[s] * "| \u2265 0.3, p < 0.05)"),
      cex.main = 1.0, line = -1.5)
circos.clear()
dev.off()

#######################
# Species-level analysis
target_sp  <- read.csv(file.path(DIR_TAB, "Species_within_saving_module.csv"), row.names = 1)
sp_in_data <- intersect(target_sp$Species, rownames(mtx_cpm_sp))
mtx_cpm_sp_filt <- mtx_cpm_sp[sp_in_data, samples_common]
mtx_clr_sp <- t(as.matrix(clr(t(mtx_cpm_sp_filt + 0.5))))

# RNA correlation (species)
rna_cor_sp <- do.call(cbind, lapply(
  split(seq_len(nrow(mtx_tpm_cor)), ceiling(seq_len(nrow(mtx_tpm_cor)) / 2000)),
  function(idx) cor(t(mtx_clr_sp), t(mtx_tpm_cor[idx, , drop = FALSE]),
                    method = "spearman", use = "pairwise.complete.obs")
))
saveRDS(rna_cor_sp, file.path(DIR_RDS, "rna_species_gene_cor.rds"))

# Protein correlation (species)
prot_cor_sp <- calc_prot_cor(mtx_clr_sp, mtx_prot_log, samples_common)
saveRDS(prot_cor_sp, file.path(DIR_RDS, "prot_species_protein_cor.rds"))

################################################
# Heatmap: pathway-annotated species correlation
pathways_ordered <- c(
  # Metabolism
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_HYPOXIA",
  # Cell Cycle
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_P53_PATHWAY",
  # Metastasis
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_KRAS_SIGNALING_DN",
  # Immune
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)

pw_to_cat <- setNames(
  c(rep("Metabolism", 5), rep("Cell Cycle", 7),
    rep("Metastasis", 4), rep("Immune", 6)),
  pathways_ordered
)

# Hallmark gene sets
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  filter(gs_name %in% target_pathways) %>%
  select(gs_name, gene_symbol) %>%
  { split(.$gene_symbol, .$gs_name) }

gsea_results_sp <- lapply(setNames(rownames(rna_cor_sp), rownames(rna_cor_sp)), function(sp) {
  rank_vec <- sort(na.omit(rna_cor_sp[sp, ]), decreasing = TRUE)
  fgsea(pathways = hallmark_sets, stats = rank_vec,
        minSize = 10, maxSize = 500, nPermSimple = 1000, eps = 0) %>%
    filter(padj < 0.05) %>%
    mutate(species = sp)
})

core_genes_sp <- lapply(names(gsea_results_sp), function(sp) {
  res <- gsea_results_sp[[sp]]
  if (nrow(res) == 0) return(NULL)
  res %>%
    select(pathway, NES, leadingEdge) %>%
    mutate(gene = lapply(leadingEdge, identity)) %>%
    unnest(gene) %>%
    select(pathway, NES, gene) %>%
    mutate(species = sp)
}) %>% bind_rows()

gene_order_sp <- core_genes_sp %>%
  distinct(gene, pathway) %>%
  mutate(pathway = factor(pathway, levels = pathways_ordered)) %>%
  arrange(pathway, gene) %>%
  group_by(gene) %>% slice(1) %>% ungroup() %>%
  arrange(pathway, gene) %>% pull(gene) %>% unique()

gene_use_sp <- intersect(gene_order_sp, colnames(rna_cor_sp))
gene_use_sp <- gene_use_sp[
  colSums(!is.na(rna_cor_sp[, gene_use_sp, drop = FALSE])) >= ceiling(0.3 * nrow(rna_cor_sp))
]

rna_mat_sp <- rna_cor_sp[, gene_use_sp]

gene_primary_pw <- core_genes_sp %>%
  distinct(gene, pathway) %>%
  mutate(pathway = factor(pathway, levels = pathways_ordered)) %>%
  arrange(pathway) %>%
  group_by(gene) %>% slice(1) %>% ungroup() %>%
  filter(gene %in% gene_use_sp)

col_split_vec <- factor(
  pw_to_cat[as.character(gene_primary_pw$pathway[match(gene_use_sp, gene_primary_pw$gene)])],
  levels = c("Metabolism", "Cell Cycle", "Metastasis", "Immune")
)

prot_mat_sp <- matrix(NA_real_, nrow = nrow(rna_mat_sp), ncol = length(gene_use_sp),
                      dimnames = list(rownames(rna_mat_sp), gene_use_sp))
gene_name_map  <- setNames(mtx_hpro$Gene_names, rownames(mtx_hpro))
for (pg in colnames(prot_cor_sp)) {
  sym <- gene_name_map[pg]
  if (!is.na(sym) && sym %in% gene_use_sp) {
    for (sp in rownames(prot_mat_sp)) {
      if (is.na(prot_mat_sp[sp, sym]))
        prot_mat_sp[sp, sym] <- prot_cor_sp[sp, pg]
    }
  }
}

fc_cap <- 0.6
rna_mat_sp_capped  <- pmax(pmin(rna_mat_sp,  fc_cap), -fc_cap)
rna_mat_sp_capped[is.na(rna_mat_sp_capped)] <- 0

all_pathways_sp <- pathways_ordered[pathways_ordered %in% unique(core_genes_sp$pathway)]
pw_short_sp     <- gsub("HALLMARK_", "", all_pathways_sp) %>% str_to_sentence()
pw_colors_sp    <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(all_pathways_sp)),
  all_pathways_sp
)

col_anno_sp <- do.call(HeatmapAnnotation, c(
  lapply(setNames(all_pathways_sp, all_pathways_sp), function(pw) {
    vec <- ifelse(gene_use_sp %in% filter(core_genes_sp, pathway == pw)$gene, pw, NA_character_)
    anno_simple(vec, which = "column", col = setNames(pw_colors_sp[pw], pw),
                na_col = "white", height = unit(2, "mm"))
  }),
  list(show_annotation_name = FALSE, which = "column")
))

# Row annotations
sp_genus    <- gsub("_.*", "", rownames(rna_mat_sp))
sp_log2cpm  <- log2(rowMeans(mtx_cpm_sp_filt) + 1)
abund_col_sp <- colorRamp2(range(sp_log2cpm), c("#e8f4f8", "#1a5276"))

prot_anno_sp <- rowAnnotation(
  Proteome_overlap = anno_simple(
    ifelse(sp_genus %in% high_overlap_g, "yes", "no"),
    col = c("yes" = "tomato", "no" = "grey80"), na_col = "grey90", width = unit(4, "mm")
  ),
  log2CPM = anno_simple(
    sp_log2cpm[rownames(rna_mat_sp)],
    col = abund_col_sp, na_col = "grey90", width = unit(4, "mm")
  ),
  annotation_name_side = "bottom",
  annotation_name_gp   = gpar(fontsize = 7)
)

cor_col_fun <- colorRamp2(
  c(-fc_cap, -0.3, 0, 0.3, fc_cap),
  c("#313695", "#74add1", "white", "#f46d43", "#a50026")
)

ht_sp <- Heatmap(
  rna_mat_sp_capped, name = "RNA\nSpearman r", col = cor_col_fun,
  cluster_rows = TRUE, cluster_columns = FALSE,
  column_split = col_split_vec,
  column_gap = unit(1.5, "mm"),
  show_row_names = TRUE, show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8), na_col = "grey90",
  top_annotation = col_anno_sp, right_annotation = prot_anno_sp,
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  use_raster = FALSE
)

pdf(file.path(DIR_RES, "Heatmap_species_pathway_RNA_protein_cor.pdf"), width = 10, height = 8)
draw(ht_sp, annotation_legend_list = list(
  Legend(labels = pw_short_sp, legend_gp = gpar(fill = pw_colors_sp),
         title = "Pathway", ncol = 2,
         title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)),
  Legend(labels = c("> 10%", "<= 10%"), legend_gp = gpar(fill = c("tomato", "grey80")),
         title = "Proteome overlap",
         title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)),
  Legend(col_fun = abund_col_sp, title = "log2CPM",
         title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7),
         direction = "vertical")
), merge_legend = FALSE)
dev.off()

####################
# Species Rs scatter
rna_pval_sp <- spearman_pval(rna_cor_sp, n_common)

sig_df <- lapply(rownames(rna_cor_sp), function(sp) {
  data.frame(species = sp, gene = colnames(rna_cor_sp),
             Rs = rna_cor_sp[sp, ], pval = rna_pval_sp[sp, ],
             stringsAsFactors = FALSE)
}) %>%
  bind_rows() %>%
  filter(!is.na(Rs), !is.na(pval)) %>%
  mutate(sig = case_when(
    Rs >= 0.3 & pval < 0.05  ~ "Pos",
    Rs <= -0.3 & pval < 0.05 ~ "Neg",
    TRUE ~ "ns"
  ))

sp_order <- intersect(names(sort(sp_log2cpm, decreasing = TRUE)), rownames(rna_cor_sp))
sig_df$species <- factor(sig_df$species, levels = sp_order)

label_df <- sig_df %>%
  filter(sig != "ns") %>%
  group_by(species, sig) %>%
  slice_max(order_by = abs(Rs), n = 3) %>%
  ungroup()

bg_df <- data.frame(
  species = sp_order,
  xmin = seq_along(sp_order) - 0.5, xmax = seq_along(sp_order) + 0.5,
  fill = ifelse(seq_along(sp_order) %% 2 == 0, "grey92", "white")
)

sp_genus_order <- gsub("_.*", "", sp_order)
label_df_box <- data.frame(
  species    = factor(sp_order, levels = sp_order),
  label_fill = ifelse(sp_genus_order %in% high_overlap_g, "tomato", "grey"),
  stringsAsFactors = FALSE
)

y_cap      <- 0.8
jitter_pos <- position_jitter(width = 0.35, height = 0, seed = 42)

p_multi_volc <- ggplot() +
  geom_rect(data = bg_df,
            aes(xmin = xmin, xmax = xmax,
                ymin = -y_cap * 1.35, ymax = y_cap * 1.35, fill = fill),
            inherit.aes = FALSE) +
  scale_fill_identity() +
  geom_point(data = filter(sig_df, sig != "ns"),
             aes(x = species, y = pmax(pmin(Rs, y_cap), -y_cap), color = sig),
             position = jitter_pos, size = 1.0, alpha = 0.7) +
  scale_color_manual(values = c("Pos" = "#c0392b", "Neg" = "#2980b9"), name = NULL) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_hline(yintercept = c(0.3, -0.3), linewidth = 0.3, color = "grey60", linetype = "dashed") +
  geom_text_repel(
    data = label_df,
    aes(x = species, y = pmax(pmin(Rs, y_cap), -y_cap), label = gene, color = sig),
    position = jitter_pos, size = 2.2, fontface = "italic",
    max.overlaps = 20, segment.size = 0.3, segment.color = "grey50",
    show.legend = FALSE, box.padding = 0.2, force = 1.5, seed = 42
  ) +
  geom_tile(data = label_df_box,
            aes(x = species, y = 0, fill = label_fill),
            height = 0.08, width = 0.9, color = "white", linewidth = 0.3,
            inherit.aes = FALSE) +
  geom_text(data = label_df_box,
            aes(x = species, y = 0, label = gsub("_", " ", species)),
            size = 2.0, color = "white", fontface = "italic",
            inherit.aes = FALSE) +
  scale_y_continuous(limits = c(-y_cap * 1.35, y_cap * 1.35),
                     breaks = seq(-0.6, 0.6, 0.3)) +
  labs(x = "Species", y = "Spearman Rs") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_rect(color = 1),
    legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 9)
  )

ggsave(file.path(DIR_RES, "Lollipop_species_Rs_genes.pdf"),
       p_multi_volc, width = max(14, length(sp_order) * 0.9), height = 7)

########
# GENIE3
sig_genes <- unique(colnames(rna_cor_sp)[colSums(
  abs(rna_cor_sp) >= 0.3 & spearman_pval(rna_cor_sp, n_common) < 0.05,
  na.rm = TRUE
) > 0])

expr_combined <- rbind(
  mtx_clr_sp[, samples_common],
  mtx_tpm_cor[sig_genes, samples_common]
)

set.seed(42)
genie3_weight <- GENIE3(
  exprMatrix = as.matrix(expr_combined),
  regulators = rownames(mtx_clr_sp),
  targets    = sig_genes,
  treeMethod = "RF", K = "sqrt", nTrees = 1000, nCores = 4
)
saveRDS(genie3_weight, file.path(DIR_RDS, "genie3_species_gene_weight.rds"))

genie3_df <- as.data.frame(as.table(genie3_weight)) %>%
  setNames(c("species", "gene", "importance")) %>%
  mutate(importance = as.numeric(importance)) %>%
  filter(importance >= 0.001) %>%
  arrange(desc(importance))

write.csv(genie3_df, file.path(DIR_TAB, "GENIE3_species_gene_links.csv"), row.names = FALSE)
write.csv(
  genie3_df %>% group_by(species) %>% slice_max(importance, n = 20) %>% ungroup(),
  file.path(DIR_TAB, "GENIE3_species_gene_top20.csv"), row.names = FALSE
)

threshold_q99 <- quantile(as.vector(genie3_weight), 0.99)

sig_pairs <- which(
  abs(rna_cor_sp) >= 0.3 &
    spearman_pval(rna_cor_sp, length(samples_common)) < 0.05,
  arr.ind = TRUE
)
genie3_filtered <- genie3_df %>%
  filter(importance >= threshold_q99) %>%
  inner_join(
    data.frame(
      species    = rownames(rna_cor_sp)[sig_pairs[, 1]],
      gene       = colnames(rna_cor_sp)[sig_pairs[, 2]],
      spearman_r = rna_cor_sp[sig_pairs]
    ),
    by = c("species", "gene")
  )

#########################################################################
# Sankey bubble plot: genus -> species -> gene (GENIE3 top 5 per species)
library(ggnewscale)

rna_cor_sp_long <- as.data.frame(as.table(rna_cor_sp)) %>%
  setNames(c("species", "gene", "spearman_r")) %>%
  mutate(species     = as.character(species),
         gene        = as.character(gene),
         spearman_r  = as.numeric(spearman_r))

top5_base <- genie3_df %>%
  group_by(species) %>% slice_max(importance, n = 5) %>% ungroup() %>%
  left_join(rna_cor_sp_long, by = c("species", "gene")) %>%
  mutate(genus   = as.character(gsub("_.*", "", species)),
         gene    = as.character(gene),
         species = as.character(species))

genus_order_s <- top5_base %>% count(genus) %>% arrange(desc(n)) %>% pull(genus)
sp_order_s    <- top5_base %>% count(species) %>%
  arrange(match(gsub("_.*", "", species), genus_order_s)) %>% pull(species)
gene_order_s  <- top5_base %>% group_by(gene) %>%
  summarise(m = mean(importance)) %>% arrange(desc(m)) %>%
  pull(gene) %>% as.character()

# Sankey bubble plot
y_genus <- 3; y_species <- 2; y_gene <- 1

x_sp <- setNames(
  seq_along(sp_order_s) * (length(genus_order_s) * 3 / length(sp_order_s)),
  sp_order_s
)
x_g <- sapply(genus_order_s, function(gen) {
  members <- sp_order_s[gsub("_.*", "", sp_order_s) == gen]
  if (length(members) == 0) return(NA_real_)
  mean(x_sp[members])
})
x_gn <- setNames(
  seq_along(gene_order_s) * (max(x_sp, na.rm = TRUE) / length(gene_order_s)),
  gene_order_s
)

genus_order_s <- genus_order_s[!is.na(x_g)]
x_g           <- x_g[genus_order_s]

# Extend genus_colors for any missing genera
missing_genera <- setdiff(unique(gsub("_.*", "", sp_order_s)), names(genus_colors))
if (length(missing_genera) > 0) {
  genus_colors <- c(genus_colors, setNames(
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(missing_genera)),
    missing_genera
  ))
}

coord_sp   <- data.frame(species = sp_order_s,  x_sp_coord = x_sp[sp_order_s],   stringsAsFactors = FALSE)
coord_gene <- data.frame(gene    = gene_order_s, x_gn_coord = x_gn[gene_order_s], stringsAsFactors = FALSE)

top5_per_sp <- top5_base %>%
  left_join(coord_sp,   by = "species") %>%
  left_join(coord_gene, by = "gene")

df_genus <- data.frame(
  label = genus_order_s, x = x_g[genus_order_s], y = y_genus,
  stringsAsFactors = FALSE
)
df_species <- data.frame(
  label    = sp_order_s,
  x        = x_sp[sp_order_s],
  y        = y_species,
  genus    = gsub("_.*", "", sp_order_s),
  sp_label = gsub("_", " ", sp_order_s),
  stringsAsFactors = FALSE
)
df_gene <- data.frame(
  gene = gene_order_s, x = x_gn[gene_order_s], y = y_gene,
  stringsAsFactors = FALSE
) %>%
  left_join(top5_per_sp %>% group_by(gene) %>% summarise(importance = mean(importance)),
            by = "gene")

make_ribbon_v <- function(x0, y0, x1, y1, w0 = 0.3, w1 = 0.15, n = 60) {
  t  <- seq(0, 1, length.out = n)
  ym <- (y0 + y1) / 2
  by <- (1-t)^3*y0 + 3*(1-t)^2*t*ym + 3*(1-t)*t^2*ym + t^3*y1
  data.frame(
    x = c((1-t)*(x0+w0) + t*(x1+w1), rev((1-t)*(x0-w0) + t*(x1-w1))),
    y = c(by, rev(by))
  )
}

ribbons <- lapply(sp_order_s, function(sp) {
  gen  <- gsub("_.*", "", sp)
  n_sp <- sum(gsub("_.*", "", sp_order_s) == gen)
  make_ribbon_v(x_g[gen], y_genus, x_sp[sp], y_species,
                w0 = 0.8 / n_sp, w1 = 0.2) %>%
    mutate(species = sp, genus = gen, fill = genus_colors[gen])
}) %>% bind_rows()

x_range <- range(c(x_gn, x_sp, x_g))

p_sankey <- ggplot() +
  geom_polygon(data = ribbons,
               aes(x, y, group = species, fill = fill), alpha = 0.3, color = NA) +
  scale_fill_identity() +
  geom_segment(data = top5_per_sp,
               aes(x = x_sp_coord, y = y_species, xend = x_gn_coord, yend = y_gene,
                   color = spearman_r),
               linewidth = 0.4, alpha = 0.6) +
  scale_color_distiller(palette = "RdBu", direction = -1,
                        limits = c(-1, 1) * max(abs(top5_per_sp$spearman_r), na.rm = TRUE),
                        name = "Spearman r") +
  geom_tile(data = df_genus,
            aes(x, y, fill = genus_colors[label]),
            width = 1.5, height = 0.12, color = "white") +
  geom_text(data = df_genus,
            aes(x, y = y + 0.1, label = label), vjust = 0, size = 2.5, fontface = "bold") +
  geom_tile(data = df_species,
            aes(x, y, fill = genus_colors[genus]),
            width = 0.6, height = 0.1, alpha = 0.85, color = "white") +
  geom_text(data = df_species,
            aes(x, y = y - 0.08, label = sp_label),
            vjust = 1, size = 1.8, fontface = "italic", angle = 90, hjust = 1) +
  ggnewscale::new_scale_fill() +
  geom_point(data = df_gene,
             aes(x, y, size = importance, fill = importance),
             shape = 21, color = "grey30", stroke = 0.3) +
  scale_fill_distiller(palette = "Oranges", direction = 1, name = "Importance") +
  scale_size_continuous(range = c(1, 5), name = "Importance") +
  geom_text(data = df_gene,
            aes(x, y = y - 0.08, label = gene),
            vjust = .5, size = 1.8, fontface = "italic", angle = 90, hjust = 1) +
  annotate("text", x = min(x_range) - 1, y = c(y_genus, y_species, y_gene),
           label = c("Genus", "Species", "Gene"), fontface = "bold", size = 3, hjust = 1) +
  coord_cartesian(xlim = c(min(x_range) - 2, max(x_range) + 1),
                  ylim = c(y_gene - 0.8, y_genus + 0.5)) +
  theme_void(base_size = 11) +
  theme(legend.position = "right",
        legend.title    = element_text(size = 9, face = "bold"),
        plot.title      = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.margin     = margin(10, 10, 30, 10)) +
  labs(title = "Species-gene regulatory network (GENIE3)")

ggsave(file.path(DIR_RES, "Sankey_GENIE3_horizontal.pdf"), p_sankey,
       width  = 16, height = 6)
