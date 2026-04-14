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
library(compositions)
library(limma)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(paletteer)

DIR_RDS  <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_RES  <- "/data/yzwang/project/AEG_seiri/results/F3_crosstalk/"
DIR_TAB  <- "/data/yzwang/project/AEG_seiri/table_infos/"

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

mtx_cpm <- readRDS(file.path(DIR_RDS, "gAEG_CPM_RNA_FiltMyco.rds"))

mtx_hpro <- readxl::read_excel(file.path(DIR_TAB, "SupData13ProteinIntensity.xlsx"))
colnames(mtx_hpro)[3:ncol(mtx_hpro)] <- {
  raw <- colnames(mtx_hpro)[3:ncol(mtx_hpro)]
  ifelse(grepl("_T$", raw),
         paste0("C", gsub("AEG|_T$", "", raw)),
         paste0("N", gsub("AEG|_N$", "", raw)))
}
rownames(mtx_hpro) <- mtx_hpro$Protein_group

samples_common <- Reduce(intersect, list(
  colnames(mtx_htpm_filt)[grep("^C", colnames(mtx_htpm_filt))],
  colnames(mtx_hpro)[grep("^C", colnames(mtx_hpro))],
  colnames(mtx_cpm)[grep("^C", colnames(mtx_cpm))]
))

#######################
# Top 20 genera: global abundance across all samples
top20_genera <- rownames(mtx_cpm)[order(rowMeans(mtx_cpm / 1e4), decreasing = TRUE)][1:20]
mean_abund   <- rowMeans(mtx_cpm[top20_genera, ] / 1e4)

# CLR on common tumour samples only (for correlation)
mtx_clr <- t(as.matrix(clr(t(mtx_cpm[top20_genera, samples_common] + 0.5))))

#######################
# RNA-seq correlation
mtx_tpm_cor <- log1p(as.matrix(mtx_htpm_filt[, samples_common]))
cv_gene      <- apply(mtx_tpm_cor, 1, function(x) sd(x) / (mean(x) + 1e-6))
mtx_tpm_cor  <- mtx_tpm_cor[cv_gene > 0.2, ]

rna_cor_mat <- do.call(cbind, lapply(
  split(seq_len(nrow(mtx_tpm_cor)), ceiling(seq_len(nrow(mtx_tpm_cor)) / 2000)),
  function(idx) cor(t(mtx_clr), t(mtx_tpm_cor[idx, , drop = FALSE]),
                    method = "spearman", use = "pairwise.complete.obs")
))
saveRDS(rna_cor_mat, file.path(DIR_RDS, "rna_genus_gene_cor.rds"))

#######################
# Proteomics correlation
mtx_prot_num <- as.matrix(mutate(mtx_hpro[, samples_common], across(everything(), as.numeric)))
rownames(mtx_prot_num) <- rownames(mtx_hpro)
mtx_prot_num[mtx_prot_num == 0] <- NA
mtx_prot_log <- log2(mtx_prot_num)
mtx_prot_log <- mtx_prot_log[rowSums(!is.na(mtx_prot_log)) >= ceiling(0.5 * length(samples_common)), ]

prot_cor_mat <- matrix(NA_real_, nrow = nrow(mtx_clr), ncol = nrow(mtx_prot_log),
                       dimnames = list(rownames(mtx_clr), rownames(mtx_prot_log)))
for (g in seq_len(nrow(mtx_clr))) {
  abund_vec <- mtx_clr[g, samples_common]
  prot_cor_mat[g, ] <- apply(mtx_prot_log, 1, function(v) {
    idx <- !is.na(v)
    if (sum(idx) < 10) return(NA_real_)
    cor(abund_vec[idx], v[idx], method = "spearman")
  })
}
saveRDS(prot_cor_mat, file.path(DIR_RDS, "prot_genus_protein_cor.rds"))

rna_cor_mat  <- readRDS(file.path(DIR_RDS, "rna_genus_gene_cor.rds"))
prot_cor_mat <- readRDS(file.path(DIR_RDS, "prot_genus_protein_cor.rds"))

#######################
# Circos lollipop
spearman_pval <- function(r_mat, n) {
  2 * pt(-abs(r_mat * sqrt((n - 2) / (1 - r_mat^2))), df = n - 2)
}

rna_sig  <- abs(rna_cor_mat) >= 0.3 & spearman_pval(rna_cor_mat,  length(samples_common)) < 0.05
prot_sig <- abs(prot_cor_mat) >= 0.3 & spearman_pval(prot_cor_mat, length(samples_common)) < 0.05
rna_sig[is.na(rna_sig)]   <- FALSE
prot_sig[is.na(prot_sig)] <- FALSE

prot_gene_map  <- setNames(mtx_hpro$Gene_names, mtx_hpro$Protein_group)
prot_sig_counts <- rowSums(prot_sig)

overlap_counts <- sapply(rownames(rna_cor_mat), function(g) {
  rna_genes  <- colnames(rna_cor_mat)[rna_sig[g, ]]
  prot_genes <- colnames(prot_cor_mat)[prot_sig[g, ]] %>%
    { prot_gene_map[.] } %>% na.omit() %>%
    paste(collapse = ";") %>% strsplit(";") %>%
    unlist() %>% trimws() %>% unique()
  length(intersect(rna_genes, prot_genes))
})

overlap_df <- data.frame(
  genus      = names(overlap_counts),
  count      = overlap_counts,
  prot_total = prot_sig_counts[names(overlap_counts)],
  stringsAsFactors = FALSE
) %>%
  mutate(pct = ifelse(prot_total > 0, count / prot_total * 100, 0)) %>%
  .[match(top20_genera, .$genus), ]

genera       <- overlap_df$genus
abund_scaled <- setNames(
  (mean_abund - min(mean_abund)) / (max(mean_abund) - min(mean_abund)),
  top20_genera
)
genus_colors <- setNames(
  as.character(paletteer_d("khroma::discreterainbow")[c(10, 12:20, 23:27, 2, 4, 5, 7, 9)]),
  top20_genera
)
max_pct <- ceiling(max(overlap_df$pct))

# Abundance track: scale mean_abund to [0, 1] for bar height
abund_scaled <- setNames(
  (mean_abund - min(mean_abund)) / (max(mean_abund) - min(mean_abund)),
  top20_genera
)
high_overlap_g <- overlap_df$genus[overlap_df$pct > 10]

pdf(file.path(DIR_RES, "Circos_genus_overlap_lollipop.pdf"), width = 10, height = 10)
circos.clear()
circos.par(start.degree = 90, gap.degree = 2,
           cell.padding = c(0, 0, 0, 0),
           points.overflow.warning = FALSE,
           track.margin = c(0.005, 0.005))

circos.initialize(factors = factor(genera, levels = genera),
                  xlim = matrix(rep(c(0, 1), length(genera)), ncol = 2, byrow = TRUE))

# Track 1 (outermost): lollipop
circos.track(
  factors = factor(genera, levels = genera),
  ylim = c(0, max_pct), track.height = 0.38,
  bg.border = NA, bg.col = NA,
  panel.fun = function(x, y) {
    g        <- get.cell.meta.data("sector.index")
    row      <- overlap_df[overlap_df$genus == g, ]
    border_col <- ifelse(g %in% high_overlap_g, "#c0392b", "grey85")
    
    # Redraw background with per-sector border color
    circos.rect(
      xleft  = get.cell.meta.data("xlim")[1],
      ybottom = 0,
      xright = get.cell.meta.data("xlim")[2],
      ytop   = max_pct,
      col    = NA,
      border = border_col,
      lwd    = ifelse(g %in% high_overlap_g, 2.5, 1)
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

# Track 2: color block + genus label
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

# Track 3 (innermost): abundance as color intensity
circos.track(
  factors = factor(genera, levels = genera),
  ylim = c(0, 1), track.height = 0.12,
  bg.border = "grey85", bg.col = "grey97",
  panel.fun = function(x, y) {
    g   <- get.cell.meta.data("sector.index")
    val <- abund_scaled[g]
    circos.rect(0, 0, 1, 1,
                col    = abund_col_fun(val),
                border = NA)
    circos.text(0.5, 0.5,
                labels     = formatC(mean_abund[g], format = "f", digits = 2),
                adj        = c(0.5, 0.5),
                cex        = 0.38,
                niceFacing = TRUE,
                facing     = "bending.inside",
                col        = ifelse(val > 0.6, "white", "black"))
  }
)

# Abundance track y-axis label
circos.text(
  x                = 0,
  y                = 0,
  labels           = "Abund.\n(%)",
  sector.index     = genera[1],
  track.index      = 3,
  facing           = "bending.inside",
  cex              = 0.4,
  col              = "black"
)

legend("center",
       legend = gsub("_", " ", genera),
       fill   = genus_colors[genera],
       border = NA, cex = 0.55, ncol = 2,
       title  = "Genus", bty = "n")

title(main = expression("Overlap genes (|" * italic(r)[s] * "| \u2265 0.3, p < 0.05)"),
      cex.main = 1.0, line = -1.5)

circos.clear()
dev.off()

#########
# Heatmap
# Species-level correlation
mtx_cpm_sp <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))

target_sp <- read.csv(file.path(DIR_TAB, "Species_within_saving_module.csv"), row.names = 1)

# Filter to target species present in common samples
sp_in_data <- intersect(target_sp$Species, rownames(mtx_cpm_sp))
mtx_cpm_sp_filt <- mtx_cpm_sp[sp_in_data, samples_common]

mtx_clr_sp <- t(as.matrix(clr(t(mtx_cpm_sp_filt + 0.5))))

# RNA correlation: species x gene
rna_cor_sp <- do.call(cbind, lapply(
  split(seq_len(nrow(mtx_tpm_cor)), ceiling(seq_len(nrow(mtx_tpm_cor)) / 2000)),
  function(idx) cor(t(mtx_clr_sp), t(mtx_tpm_cor[idx, , drop = FALSE]),
                    method = "spearman", use = "pairwise.complete.obs")
))
saveRDS(rna_cor_sp, file.path(DIR_RDS, "rna_species_gene_cor.rds"))

# Protein correlation: species x protein
prot_cor_sp <- matrix(NA_real_, nrow = nrow(mtx_clr_sp), ncol = nrow(mtx_prot_log),
                      dimnames = list(rownames(mtx_clr_sp), rownames(mtx_prot_log)))
for (g in seq_len(nrow(mtx_clr_sp))) {
  abund_vec <- mtx_clr_sp[g, samples_common]
  prot_cor_sp[g, ] <- apply(mtx_prot_log, 1, function(v) {
    idx <- !is.na(v)
    if (sum(idx) < 10) return(NA_real_)
    cor(abund_vec[idx], v[idx], method = "spearman")
  })
}
saveRDS(prot_cor_sp, file.path(DIR_RDS, "prot_species_protein_cor.rds"))

rna_cor_sp  <- readRDS(file.path(DIR_RDS, "rna_species_gene_cor.rds"))
prot_cor_sp <- readRDS(file.path(DIR_RDS, "prot_species_protein_cor.rds"))

# Heatmap: pathway-annotated species correlation
gsea_results_sp <- lapply(rownames(rna_cor_sp), function(sp) {
  rank_vec <- rna_cor_sp[sp, ]
  rank_vec <- sort(rank_vec[!is.na(rank_vec)], decreasing = TRUE)
  
  fgsea(pathways    = hallmark_sets,
        stats       = rank_vec,
        minSize     = 10,
        maxSize     = 500,
        nPermSimple = 1000,
        eps         = 0) %>%
    filter(padj < 0.05) %>%
    mutate(species = sp)
})
names(gsea_results_sp) <- rownames(rna_cor_sp)

core_genes_sp <- lapply(names(gsea_results_sp), function(sp) {
  res <- gsea_results_sp[[sp]]
  if (nrow(res) == 0) return(NULL)
  res %>%
    select(pathway, NES, leadingEdge) %>%
    mutate(gene = lapply(leadingEdge, identity)) %>%
    unnest(gene) %>%
    select(pathway, NES, gene) %>%
    mutate(species = sp)
}) %>%
  bind_rows()

gene_order_sp <- core_genes_sp %>%
  distinct(gene, pathway) %>%
  mutate(pathway = factor(pathway, levels = sort(unique(pathway)))) %>%
  arrange(pathway, gene) %>%
  group_by(gene) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  arrange(pathway, gene) %>%
  pull(gene) %>%
  unique()

gene_use_sp <- gene_order_sp[gene_order_sp %in% colnames(rna_cor_sp)]
gene_use_sp <- gene_use_sp[
  rowSums(!is.na(rna_cor_sp[, gene_use_sp, drop = FALSE])) >=
    ceiling(0.3 * nrow(rna_cor_sp))
]

rna_mat_sp <- rna_cor_sp[rownames(rna_cor_sp), gene_use_sp]

prot_mat_sp <- matrix(NA_real_, nrow = nrow(rna_mat_sp), ncol = length(gene_use_sp),
                      dimnames = list(rownames(rna_mat_sp), gene_use_sp))
for (pg in colnames(prot_cor_sp)) {
  sym <- gene_name_map[pg]
  if (!is.na(sym) && sym %in% gene_use_sp) {
    for (sp in rownames(prot_mat_sp)) {
      if (is.na(prot_mat_sp[sp, sym]))
        prot_mat_sp[sp, sym] <- prot_cor_sp[sp, pg]
    }
  }
}

rna_mat_sp_capped  <- pmax(pmin(rna_mat_sp,  fc_cap), -fc_cap)
prot_mat_sp_capped <- pmax(pmin(prot_mat_sp, fc_cap), -fc_cap)

all_pathways_sp <- sort(unique(core_genes_sp$pathway))
pw_short_sp     <- gsub("HALLMARK_", "", all_pathways_sp) %>% str_to_sentence()
pw_colors_sp    <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(all_pathways_sp)),
  all_pathways_sp
)

pathway_anno_sp <- lapply(all_pathways_sp, function(pw) {
  vec <- ifelse(gene_use_sp %in% filter(core_genes_sp, pathway == pw)$gene, pw, NA_character_)
  anno_simple(vec, which = "column",
              col    = setNames(pw_colors_sp[pw], pw),
              na_col = "white",
              height = unit(2, "mm"))
})
names(pathway_anno_sp) <- all_pathways_sp

col_anno_sp <- do.call(HeatmapAnnotation, c(
  pathway_anno_sp,
  list(show_annotation_name = FALSE, which = "column")
))

sp_genus <- gsub("_.*", "", rownames(rna_mat_sp))

sp_log2cpm <- log2(rowMeans(mtx_cpm_sp_filt) + 1)

abund_col_sp <- colorRamp2(
  c(min(sp_log2cpm), max(sp_log2cpm)),
  c("#e8f4f8", "#1a5276")
)

prot_anno_sp <- rowAnnotation(
  Proteome_overlap = anno_simple(
    ifelse(sp_genus %in% high_overlap_g, "yes", "no"),
    col    = c("yes" = "tomato", "no" = "grey80"),
    na_col = "grey90",
    width  = unit(4, "mm")
  ),
  log2CPM = anno_simple(
    sp_log2cpm[rownames(rna_mat_sp)],
    col    = abund_col_sp,
    na_col = "grey90",
    width  = unit(4, "mm")
  ),
  annotation_name_side = "bottom",
  annotation_name_gp   = gpar(fontsize = 7)
)

pw_legend_sp <- Legend(
  labels    = pw_short_sp,
  legend_gp = gpar(fill = pw_colors_sp),
  title     = "Pathway",
  ncol      = 2,
  title_gp  = gpar(fontsize = 8, fontface = "bold"),
  labels_gp = gpar(fontsize = 7)
)

prot_legend <- Legend(
  labels    = c("> 10%", "<= 10%"),
  legend_gp = gpar(fill = c("tomato", "grey80")),
  title     = "Proteome overlap",
  title_gp  = gpar(fontsize = 8, fontface = "bold"),
  labels_gp = gpar(fontsize = 7)
)
abund_legend <- Legend(
  col_fun   = abund_col_sp,
  title     = "log2CPM",
  title_gp  = gpar(fontsize = 8, fontface = "bold"),
  labels_gp = gpar(fontsize = 7),
  direction = "vertical"
)

rna_mat_sp_capped[is.na(rna_mat_sp_capped)] <- 0

ht_sp <- Heatmap(
  rna_mat_sp_capped,
  name              = "RNA\nSpearman r",
  col               = cor_col_fun,
  cluster_rows      = TRUE,
  cluster_columns   = FALSE,
  show_row_names    = TRUE,
  show_column_names = FALSE,
  row_names_gp      = gpar(fontsize = 8),
  na_col            = "grey90",
  top_annotation    = col_anno_sp,
  right_annotation  = prot_anno_sp,
  column_title      = "Core enrichment genes (selected Hallmark pathways)",
  use_raster        = FALSE
)

pdf(file.path(DIR_RES, "Heatmap_species_pathway_RNA_protein_cor.pdf"),
    width = 10, height = 8)
draw(ht_sp,
     annotation_legend_list = list(pw_legend_sp, prot_legend, abund_legend),
     merge_legend = FALSE)
dev.off()

library(ggrepel)

# Significant RNA correlations per species
rna_pval_sp <- spearman_pval(rna_cor_sp, length(samples_common))

sig_df <- lapply(rownames(rna_cor_sp), function(sp) {
  r_vec <- rna_cor_sp[sp, ]
  p_vec <- rna_pval_sp[sp, ]
  data.frame(
    species  = sp,
    gene     = colnames(rna_cor_sp),
    Rs       = r_vec,
    pval     = p_vec,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  filter(!is.na(Rs), !is.na(pval)) %>%
  mutate(
    sig = case_when(
      Rs >= 0.3 & pval < 0.05  ~ "Pos",
      Rs <= -0.3 & pval < 0.05 ~ "Neg",
      TRUE                      ~ "ns"
    )
  )

# Species order: by abundance (high to low)
sp_order <- names(sort(sp_log2cpm, decreasing = TRUE))
sp_order <- sp_order[sp_order %in% rownames(rna_cor_sp)]
sig_df$species <- factor(sig_df$species, levels = sp_order)

# Top 3 up and down per species for labelling
label_df <- sig_df %>%
  filter(sig != "ns") %>%
  group_by(species, sig) %>%
  slice_max(order_by = abs(Rs), n = 3) %>%
  ungroup()

# Alternating background shading
bg_df <- data.frame(
  species = sp_order,
  xmin    = seq_along(sp_order) - 0.5,
  xmax    = seq_along(sp_order) + 0.5,
  fill    = ifelse(seq_along(sp_order) %% 2 == 0, "grey92", "white")
)

sp_genus_order <- gsub("_.*", "", sp_order)
label_fill <- ifelse(sp_genus_order %in% high_overlap_g, "tomato", "grey")
label_df_box <- data.frame(
  species    = factor(sp_order, levels = sp_order),
  label_fill = label_fill,
  stringsAsFactors = FALSE
)

# Y cap for display
y_cap <- 0.8
jitter_pos <- position_jitter(width = 0.35, height = 0, seed = 42)

p <- ggplot() +
  # Alternating background
  geom_rect(data = bg_df,
            aes(xmin = xmin, xmax = xmax, ymin = -y_cap * 1.35, ymax = y_cap * 1.35,
                fill = fill),
            inherit.aes = FALSE) +
  scale_fill_identity() +
  # Significant points
  geom_point(data = filter(sig_df, sig != "ns"),
             aes(x = species, y = pmax(pmin(Rs, y_cap), -y_cap), color = sig),
             position = jitter_pos,
             size = 1.0, alpha = 0.7) +
  scale_color_manual(
    values = c("Pos" = "#c0392b", "Neg" = "#2980b9"),
    name   = NULL
  ) +
  # Zero line
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  # Threshold lines
  geom_hline(yintercept =  0.3, linewidth = 0.3, color = "grey60", linetype = "dashed") +
  geom_hline(yintercept = -0.3, linewidth = 0.3, color = "grey60", linetype = "dashed") +
  # Gene labels
  geom_text_repel(
    data = label_df,
    aes(x = species, y = pmax(pmin(Rs, y_cap), -y_cap),
        label = gene, color = sig),
    position      = jitter_pos,
    size          = 2.2,
    fontface      = "italic",
    max.overlaps  = 20,
    segment.size  = 0.3,
    segment.color = "grey50",
    show.legend   = FALSE,
    box.padding   = 0.2,
    force         = 1.5,
    seed          = 42
  ) +
  # Species label boxes at y = 0
  geom_tile(data = label_df_box,
            aes(x = species, y = 0, fill = label_fill),
            height = 0.08, width = 0.9,
            color  = "white", linewidth = 0.3,
            inherit.aes = FALSE) +
  geom_text(data = label_df_box,
            aes(x = species, y = 0,
                label = gsub("_", " ", species)),
            size     = 2.0,
            color    = "white",
            fontface = "italic",
            angle    = 0,
            inherit.aes = FALSE) +
  # Axes
  scale_y_continuous(
    limits = c(-y_cap * 1.35, y_cap * 1.35),
    breaks = seq(-0.6, 0.6, 0.3)
  ) +
  labs(x = "Species", y = "Spearman Rs") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = 1),
    legend.position   = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "grey80"),
    legend.key.size   = unit(0.4, "cm"),
    legend.text       = element_text(size = 9)
  )

ggsave(file.path(DIR_RES, "Lollipop_species_Rs_genes.pdf"),
       p, width = max(14, length(sp_order) * 0.9), height = 7)
