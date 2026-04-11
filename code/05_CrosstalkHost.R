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
# Top 20 genera by mean relative abundance
mtx_cpm_filt  <- mtx_cpm[rowSums(mtx_cpm[, samples_common] > 1) >= ceiling(0.2 * length(samples_common)), samples_common]
top20_genera  <- rownames(mtx_cpm_filt)[order(rowMeans(mtx_cpm_filt / 1e4), decreasing = TRUE)][1:20]
mtx_clr       <- t(as.matrix(clr(t(mtx_cpm_filt[top20_genera, ] + 0.5))))

###########
# PERMANOVA
clinical_sub <- clinical %>%
  mutate(`No.` = paste0("C", `No.`)) %>%
  filter(No. %in% samples_common) %>%
  mutate(
    tumor_dims    = str_extract_all(`Primary tumor size`, "[0-9.]+"),
    tumor_max_dim = map_dbl(tumor_dims, ~ max(as.numeric(.))),
    T_stage       = str_extract(`TNM stage`, "T[0-9a-z]+"),
    N_stage       = str_extract(`TNM stage`, "N[0-9]+"),
    M_stage       = str_extract(`TNM stage`, "M[0-9]"),
    stage_group   = case_when(
      `Pathological stage` %in% c("I","IA","IB","II","IIA","IIB") ~ "Early",
      `Pathological stage` %in% c("III","IIIA","IIIB","IIIC","IV") ~ "Late",
      TRUE ~ NA_character_)
  )

meta_permanova <- clinical_sub %>%
  mutate(Age_group = if_else(Age >= median(Age, na.rm = TRUE), "High", "Low")) %>%
  select(No., Age_group, Sex, Smoking, Alcohol,
         `Siewert type`, `Borrmann classification`, `Lauren Classification`,
         stage_group, N_stage, M_stage, tumor_max_dim) %>%
  column_to_rownames("No.") %>%
  drop_na()

dist_sub <- as.dist(as.matrix(dist(t(mtx_clr)))[rownames(meta_permanova), rownames(meta_permanova)])

adonis2(
  dist_sub ~ stage_group + `Lauren Classification` + `Borrmann classification` +
    `Siewert type` + N_stage + M_stage + tumor_max_dim + Age_group + Sex + Smoking + Alcohol,
  data = meta_permanova, permutations = 999, by = "margin"
)

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
mtx_prot_log <- log2(replace(
  as.matrix(mutate(mtx_hpro[, samples_common], across(everything(), as.numeric))),
  as.matrix(mutate(mtx_hpro[, samples_common], across(everything(), as.numeric))) == 0, NA
))
rownames(mtx_prot_log) <- rownames(mtx_hpro)
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

#######################
# Dual heatmap
rna_cor_mat  <- readRDS(file.path(DIR_RDS, "rna_genus_gene_cor.rds"))
prot_cor_mat <- readRDS(file.path(DIR_RDS, "prot_genus_protein_cor.rds"))

gene_var  <- apply(rna_cor_mat, 2, var, na.rm = TRUE)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:200]

gene_name_map  <- setNames(mtx_hpro$Gene_names, rownames(mtx_hpro))
prot_vis       <- prot_cor_mat
colnames(prot_vis) <- ifelse(!is.na(gene_name_map[colnames(prot_vis)]),
                             gene_name_map[colnames(prot_vis)], colnames(prot_vis))

cor_col_fun <- colorRamp2(c(-0.6, -0.3, 0, 0.3, 0.6),
                          c("#313695","#74add1","white","#f46d43","#a50026"))

pdf(file.path(DIR_RES, "Heatmap_genus_RNA_protein_correlation.pdf"), width = 16, height = 6)
draw(
  Heatmap(pmax(pmin(rna_cor_mat[, top_genes], 0.6), -0.6),
          name = "RNA\nSpearman r", col = cor_col_fun,
          cluster_rows = TRUE, cluster_columns = TRUE,
          show_row_names = TRUE, show_column_names = FALSE,
          row_names_gp = gpar(fontsize = 8),
          column_title = "Genus - Gene TPM correlation",
          use_raster = TRUE, na_col = "grey90") +
    Heatmap(pmax(pmin(prot_vis, 0.6), -0.6),
            name = "Protein\nSpearman r", col = cor_col_fun,
            cluster_rows = TRUE, cluster_columns = TRUE,
            show_row_names = TRUE, show_column_names = FALSE,
            row_names_gp = gpar(fontsize = 8),
            column_title = "Genus - Protein Intensity correlation",
            use_raster = TRUE, na_col = "grey90")
)
dev.off()

#######################
# Circos lollipop
mean_abund <- rowMeans(mtx_cpm_filt[top20_genera, ] / 1e4)
genus_colors <- setNames(
  as.character(paletteer_d("khroma::discreterainbow")[c(10, 12:20, 23:27, 2, 4, 5, 7, 9)]),
  top20_genera
)

overlap_df <- data.frame(
  genus      = names(overlap_counts),
  count      = overlap_counts,
  prot_total = prot_sig_counts[names(overlap_counts)],
  stringsAsFactors = FALSE
) %>%
  mutate(pct = ifelse(prot_total > 0, count / prot_total * 100, 0)) %>%
  .[match(genera_by_abund, .$genus), ]

genera  <- overlap_df$genus
max_pct <- ceiling(max(overlap_df$pct))

# Abundance track: scale mean_abund to [0, 1] for bar height
abund_scaled <- setNames(
  (mean_abund - min(mean_abund)) / (max(mean_abund) - min(mean_abund)),
  top20_genera
)
top5_genera <- overlap_df$genus[order(overlap_df$pct, decreasing = TRUE)][1:5]

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
    border_col <- ifelse(g %in% top5_genera, "#c0392b", "grey85")
    
    # Redraw background with per-sector border color
    circos.rect(
      xleft  = get.cell.meta.data("xlim")[1],
      ybottom = 0,
      xright = get.cell.meta.data("xlim")[2],
      ytop   = max_pct,
      col    = NA,
      border = border_col,
      lwd    = ifelse(g %in% top5_genera, 2.5, 1)
    )
    
    circos.segments(0.5, 0, 0.5, row$pct, lwd = 2.2, col = genus_colors[g])
    if (row$pct > 0)
      circos.points(0.5, row$pct, pch = 16,
                    cex = 1.5 + row$pct / max_pct * 1.1, col = genus_colors[g])
    circos.text(0.5, row$pct,
                labels = sprintf("%.1f%%\n(n=%d)", row$pct, row$count),
                adj = c(0.5, -0.3), cex = 0.45,
                niceFacing = TRUE, facing = "clockwise", col = genus_colors[g])
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
                facing = "clockwise", niceFacing = TRUE,
                cex = 0.55, font = 2, col = "white")
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
                facing     = "clockwise",
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
