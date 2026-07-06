####################################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 06. Microbe abundance vs. predicted drug response.
#   Panels: IDWAS  +  CaDRReS  (CaDRReS replaces the former CARE panel)
#   Phosphoproteomics analysis has been removed (see 05/earlier steps).
#
# CaDRReS predictions are produced offline by
#   drugpred/prog/50_cadrres.py  ->  drugpred/out/cadrres_pred_all.csv
# (drugs x samples, ln IC50, higher = more resistant, i.e. same
#  convention as IDWAS; no sign flip needed.)
####################################################################
library(dplyr)
library(tidyr)
library(paletteer)
library(ComplexHeatmap)
library(circlize)
library(psych)
library(ggplot2)
library(RColorBrewer)
library(ggvenn)

set.seed(42)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS  <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TAB  <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_FIG  <- "/data/yzwang/project/AEG_seiri/results/F4_drug_response/"
DIR_PRED <- "/data/yzwang/project/AEG_seiri/drugpred/out/"   # CaDRReS predictions

mtx_cpm   <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
target_sp <- read.csv(file.path(DIR_TAB, "Species_within_saving_module.csv"))

###############
# Drug response
# IDWAS
idwas_mtx <- readRDS(file.path(DIR_RDS, "drugPred_Mtx.rds"))
abund_mtx <- mtx_cpm[target_sp$Species, grep("C", colnames(mtx_cpm))] / 1e+4
abund_mtx_all <- mtx_cpm[target_sp$Species, ] / 1e+4  # all samples (tumour+normal)

gdsc_target <- read.csv(file.path(DIR_TAB, "Drug_list2023.csv"))

apply_rename <- function(m, map) {
  rn <- rownames(m)
  rn[rn %in% names(map)] <- map[rn[rn %in% names(map)]]
  rownames(m) <- rn
  m
}

build_pathway_annotation <- function(drugs, gdsc_target) {
  pathway_all <- gdsc_target[gdsc_target$Name %in% drugs,
                             c("Name", "Target.pathway")]
  pathway_all <- pathway_all[!duplicated(pathway_all$Name), ]
  rename_map <- setNames(setdiff(drugs, pathway_all$Name),
                         setdiff(drugs, pathway_all$Name))

  for (drug in names(rename_map)) {
    hit <- grep(drug, gdsc_target$Synonyms)
    drug_name <- gdsc_target$Name[hit][1]
    drug_pathway <- gdsc_target$Target.pathway[hit][1]
    if (is.na(drug_name)) {
      pathway_all <- rbind(pathway_all, c(drug, "Unannotated"))
    } else {
      rename_map[drug] <- drug_name
      pathway_all <- rbind(pathway_all, c(drug_name, drug_pathway))
    }
  }

  pathway_all <- pathway_all[!duplicated(pathway_all$Name), ]
  rownames(pathway_all) <- pathway_all$Name
  list(pathway_all = pathway_all, rename_map = rename_map)
}

sig_in <- function(R, P) {
  apply(R, 1, function(x) any(abs(x) > 0.2, na.rm = TRUE)) &
    apply(P, 1, function(x) any(x < 0.05, na.rm = TRUE))
}

# Plot predicted drug responses: idwas
pdf(file.path(DIR_FIG, "A_heatmap_idwas_samples.pdf"), width = 6, height = 10)
Heatmap(scale(idwas_mtx, center = T), show_column_names = F,
        name = "Scaled response", column_title = "Samples",
        column_title_side = "bottom",
        col = circlize::colorRamp2(c(-10, 0, 10), c("#2166AC", "white", "#B2182B")),
        row_names_gp = grid::gpar(fontsize = 6))
dev.off()

rm_idwas <- c()
for (i in 1:nrow(idwas_mtx)) {
  if (length(na.omit(idwas_mtx[i, ])) < 0.5 * ncol(idwas_mtx)) {
    rm_idwas <- c(rm_idwas, i)
  }
}
idwas_mtx <- idwas_mtx[-rm_idwas, ]
dim(idwas_mtx)

###############
# CaDRReS  (replaces the former CARE panel)
# Load precomputed CaDRReS predictions (drugs x samples, ln IC50,
# higher = more resistant -- already aligned with IDWAS).
cadrres_pred_matrix <- as.matrix(read.csv(
  file.path(DIR_PRED, "cadrres_pred_all.csv"), row.names = 1, check.names = FALSE))
dim(cadrres_pred_matrix)

pdf(file.path(DIR_FIG, "Heatmap_cadrres_samples.pdf"), width = 8, height = 10)
Heatmap(cadrres_pred_matrix, name = "CaDRReS response", show_column_names = F)
dev.off()

# Target species
# in idwas
# per-species Spearman corr (each species adjusted over drugs, matching the EF call)
corr_by_species <- function(feat_mtx, ab_mtx, samp) {
  sps <- rownames(ab_mtx)
  cts <- lapply(sps, function(s)
    corr.test(t(as.matrix(feat_mtx[, samp])),
              t(as.matrix(ab_mtx[s, samp, drop = FALSE])),
              method = "spearman"))
  list(r     = `colnames<-`(sapply(cts, function(z) z$r[, 1]),     sps),
       p.adj = `colnames<-`(sapply(cts, function(z) z$p.adj[, 1]), sps))
}
common_sp_idwas <- intersect(colnames(idwas_mtx), colnames(abund_mtx_all))
cor_sp_idwas <- corr_by_species(idwas_mtx, abund_mtx_all, common_sp_idwas)
saveRDS(cor_sp_idwas, file.path(DIR_RDS, "Correlation_species_idwas.rds"))
cor_r_idwas <- cor_sp_idwas$r
cor_p_idwas <- cor_sp_idwas$p.adj

# Drug name normalisation
rownames(cor_r_idwas) <- gsub("\\.", "-", rownames(cor_r_idwas))
rownames(cor_p_idwas) <- gsub("\\.", "-", rownames(cor_p_idwas))

drug_anno <- build_pathway_annotation(rownames(cor_r_idwas), gdsc_target)
cor_r_idwas <- apply_rename(cor_r_idwas, drug_anno$rename_map)
cor_p_idwas <- apply_rename(cor_p_idwas, drug_anno$rename_map)
pathway_all <- drug_anno$pathway_all

# Filter
keep <- rownames(cor_r_idwas)[sig_in(cor_r_idwas, cor_p_idwas)]

# FDA-style filter from original code: drop rows containing "-" or digits, keep those with lowercase letters
fda_filter <- function(x) {
  x <- x[!grepl("\\-|[1-9]", x)]
  x <- x[grepl("[a-z]", x)]
  x
}
keep <- fda_filter(keep)
keep <- intersect(keep, rownames(cor_r_idwas))

cor_r_idwas <- cor_r_idwas[keep, , drop = FALSE]
cor_p_idwas <- cor_p_idwas[keep, , drop = FALSE]

# Right pathway annotation
rownames(pathway_all) <- pathway_all$Name
pathway_vec <- pathway_all[keep, "Target.pathway"]

upw_ <- unique(pathway_vec)
is_grey_ <- is.na(upw_) | upw_ %in% c("Other", "Unannotated")
col_pathway <- setNames(rep("grey85", length(upw_)), upw_)          # Unannotated / NA
col_pathway[!is_grey_] <- colorRampPalette(brewer.pal(8, "Set2")[1:7])(sum(!is_grey_))  # drop Set2's grey; reserve grey for "Other"
col_pathway[!is.na(upw_) & upw_ == "Other"] <- "grey70"            # "Other" -> grey

right_anno <- rowAnnotation(
  Pathway = pathway_vec,
  col = list(Pathway = col_pathway),
  show_annotation_name = TRUE
)

# Top annotation
build_top_anno <- function(R, P, sp_cols, show_legend = TRUE) {
  neg_v <- vapply(seq_len(ncol(R)), function(i)
    sum(R[, i] < -0.2 & P[, i] < 0.05, na.rm = TRUE), integer(1))
  pos_v <- vapply(seq_len(ncol(R)), function(i)
    sum(R[, i] >  0.2 & P[, i] < 0.05, na.rm = TRUE), integer(1))
  mean_abund <- rowMeans(log2(mtx_cpm[sp_cols, , drop = FALSE] + 1), na.rm = TRUE)  # true mean log2(CPM+1)

  HeatmapAnnotation(
    log2CPM = mean_abund,
    Positive  = anno_barplot(pos_v,
                             gp = gpar(border = NA, fill = "#701145FF", lty = "blank")),
    Negative  = anno_barplot(neg_v,
                             gp = gpar(border = NA, fill = "#008280FF", lty = "blank")),
    col = list(log2CPM = colorRamp2(range(mean_abund, na.rm = TRUE), c("white", "darkgreen"))),
    show_legend = show_legend,
    show_annotation_name = show_legend
  )
}
top_anno_idwas <- build_top_anno(cor_r_idwas, cor_p_idwas,
                                 colnames(cor_r_idwas), show_legend = TRUE)

# Significant star
make_cell_fun <- function(P) {
  function(j, i, x, y, w, h, fill) {
    p <- P[i, j]
    if (is.na(p)) return(invisible(NULL))
    sym <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else NULL
    if (!is.null(sym)) {
      gb   <- textGrob(sym)
      gb_h <- convertHeight(grobHeight(gb), "mm")
      grid.text(sym, x, y - gb_h * 0.2, gp = gpar(col = 1, cex = .5))
    }
  }
}

# Color scale & column-name bottom annotation
rng <- range(cor_r_idwas, na.rm = TRUE)
col_fun <- colorRamp2(seq(-max(abs(rng)), max(abs(rng)), length.out = 11),
                      rev(brewer.pal(11, "RdBu")))

cn_idwas <- colnames(cor_r_idwas)

bottom_anno_idwas <- HeatmapAnnotation(
  text = anno_text(cn_idwas, rot = 60, location = unit(1, "npc"), just = "right"),
  annotation_height = max_text_width(cn_idwas)
)

ht_idwas <- Heatmap(
  cor_r_idwas,
  name                = "Rs",
  col                 = col_fun,
  column_title        = "iDWAS",
  row_names_gp        = gpar(fontsize = 9),
  show_row_names      = TRUE,
  show_column_names   = FALSE,
  cluster_rows        = TRUE,
  cluster_columns     = TRUE,
  top_annotation      = top_anno_idwas,
  bottom_annotation   = bottom_anno_idwas,
  right_annotation    = right_anno,
  cell_fun            = make_cell_fun(cor_p_idwas),
  show_heatmap_legend = TRUE
)

pdf(file.path(DIR_FIG, "Heatmap_idwas.pdf"), width = 6, height = 8)
draw(ht_idwas,
     merge_legends          = TRUE,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right")
dev.off()

# Per-drug significance flags across all species
sig_pos_drug <- apply(cor_r_idwas > 0.2 & cor_p_idwas < 0.05, 1,
                      function(x) any(x, na.rm = TRUE))
sig_neg_drug <- apply(cor_r_idwas < -0.2 & cor_p_idwas < 0.05, 1,
                      function(x) any(x, na.rm = TRUE))

drug_sig_df <- data.frame(
  Drug    = rownames(cor_r_idwas),
  Pathway = pathway_vec,
  Pos     = sig_pos_drug,
  Neg     = sig_neg_drug,
  stringsAsFactors = FALSE
)

# Aggregate per pathway
pathway_count <- drug_sig_df %>%
  group_by(Pathway) %>%
  summarise(
    Positive = sum(Pos),
    Negative = sum(Neg),
    Total    = Positive + Negative,
    .groups  = "drop"
  ) %>%
  arrange(Total)

# Long format for stacked bar
pathway_long <- pathway_count %>%
  dplyr::select(Pathway, Positive, Negative, Total) %>%
  pivot_longer(cols = c("Negative", "Positive"),
               names_to = "Direction", values_to = "N") %>%
  mutate(
    Pathway   = factor(Pathway, levels = pathway_count$Pathway),
    Direction = factor(Direction, levels = c("Negative", "Positive"))
  )

cor_cols <- c(Negative = "#2E5C8A", Positive = "#E8A33D")

p_path <- ggplot(pathway_long, aes(x = N, y = Pathway, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_text(
    data = pathway_count,
    aes(x = Total, y = Pathway, label = Total),
    inherit.aes = FALSE,
    hjust = -0.3, size = 3
  ) +
  scale_fill_manual(values = cor_cols, name = "Correlation") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Number of drugs", y = NULL,
       title = "Drug count per target pathway") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    plot.title  = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(DIR_FIG, "Barplot_pathway_drug_count.pdf"),
  plot     = p_path,
  width    = 4,
  height   = 5
)

# in CaDRReS
common_sp_cadrres <- intersect(colnames(cadrres_pred_matrix), colnames(abund_mtx_all))
cor_sp_cadrres <- corr_by_species(cadrres_pred_matrix, abund_mtx_all, common_sp_cadrres)
saveRDS(cor_sp_cadrres, file.path(DIR_RDS, "Correlation_species_cadrres.rds"))
cor_r_cadrres <- cor_sp_cadrres$r
cor_p_cadrres <- cor_sp_cadrres$p.adj

# Drug name normalisation
rownames(cor_r_idwas)   <- gsub("\\.", "-", rownames(cor_r_idwas))
rownames(cor_p_idwas)   <- gsub("\\.", "-", rownames(cor_p_idwas))
rownames(cor_r_cadrres) <- gsub("\\.", "-", rownames(cor_r_cadrres))
rownames(cor_p_cadrres) <- gsub("\\.", "-", rownames(cor_p_cadrres))

drug_anno <- build_pathway_annotation(
  union(rownames(cor_r_idwas), rownames(cor_r_cadrres)),
  gdsc_target
)
cor_r_idwas   <- apply_rename(cor_r_idwas, drug_anno$rename_map)
cor_p_idwas   <- apply_rename(cor_p_idwas, drug_anno$rename_map)
cor_r_cadrres <- apply_rename(cor_r_cadrres, drug_anno$rename_map)
cor_p_cadrres <- apply_rename(cor_p_cadrres, drug_anno$rename_map)
pathway_all <- drug_anno$pathway_all

# Filter: significant in either panel
common_drugs <- intersect(rownames(cor_r_idwas), rownames(cor_r_cadrres))
keep <- common_drugs[
  sig_in(cor_r_idwas[common_drugs, , drop = FALSE],
         cor_p_idwas[common_drugs, , drop = FALSE]) |
    sig_in(cor_r_cadrres[common_drugs, , drop = FALSE],
           cor_p_cadrres[common_drugs, , drop = FALSE])
]

# Concordance filter: overall effect sign must agree across panels
common_sp <- intersect(colnames(cor_r_idwas), colnames(cor_r_cadrres))

concordant <- vapply(keep, function(d) {
  r1 <- cor_r_idwas[d, common_sp]
  r2 <- cor_r_cadrres[d, common_sp]
  ok <- is.finite(r1) & is.finite(r2)
  if (sum(ok) < 3) return(FALSE)

  r1 <- r1[ok]; r2 <- r2[ok]

  # 1) overall direction: mean signs must match (and be non-trivial)
  m1 <- mean(r1); m2 <- mean(r2)
  if (sign(m1) != sign(m2) || sign(m1) == 0) return(FALSE)

  # 2) sign agreement at the species level: majority of signs must match
  sign_agree <- mean(sign(r1) == sign(r2))
  if (sign_agree < 0.5) return(FALSE)

  # 3) pattern correlation
  rho <- suppressWarnings(cor(r1, r2, method = "pearson"))
  isTRUE(rho > 0)
}, logical(1))

keep <- keep[concordant]

keep <- intersect(keep, rownames(cor_r_idwas))
keep <- intersect(keep, rownames(cor_r_cadrres))

cor_r_idwas   <- cor_r_idwas[keep, , drop = FALSE]
cor_p_idwas   <- cor_p_idwas[keep, , drop = FALSE]
cor_r_cadrres <- cor_r_cadrres[keep,  , drop = FALSE]
cor_p_cadrres <- cor_p_cadrres[keep,  , drop = FALSE]

# Right pathway annotation
rownames(pathway_all) <- pathway_all$Name
pathway_vec <- pathway_all[keep, "Target.pathway"]

upw_ <- unique(pathway_vec)
is_grey_ <- is.na(upw_) | upw_ %in% c("Other", "Unannotated")
col_pathway <- setNames(rep("grey85", length(upw_)), upw_)          # Unannotated / NA
col_pathway[!is_grey_] <- colorRampPalette(brewer.pal(8, "Set2")[1:7])(sum(!is_grey_))  # drop Set2's grey; reserve grey for "Other"
col_pathway[!is.na(upw_) & upw_ == "Other"] <- "grey70"            # "Other" -> grey

right_anno <- rowAnnotation(
  Pathway = pathway_vec,
  col = list(Pathway = col_pathway),
  show_annotation_name = TRUE
)

# ---- Merged split-cell heatmap ------------------------------------
# One heatmap of drugs x species; every cell is split into two halves:
#   LEFT  half = IDWAS Rs,  RIGHT half = CaDRReS Rs (shared colour scale).
# Significance stars (*/**/***) drawn on each half.
sp_common <- intersect(colnames(cor_r_idwas), colnames(cor_r_cadrres))
R_i <- cor_r_idwas[keep, sp_common, drop = FALSE]
R_c <- cor_r_cadrres[keep, sp_common, drop = FALSE]
P_i <- cor_p_idwas[keep, sp_common, drop = FALSE]
P_c <- cor_p_cadrres[keep, sp_common, drop = FALSE]

rng <- range(c(R_i, R_c), na.rm = TRUE)
col_fun <- colorRamp2(seq(-max(abs(rng)), max(abs(rng)), length.out = 11),
                      rev(brewer.pal(11, "RdBu")))
fill_col <- function(v) { z <- col_fun(v); z[is.na(z)] <- "grey95"; z }

# row/column order from the mean pattern of the two methods
M_order <- (R_i + R_c) / 2
M_order[is.na(M_order)] <- 0

star_vec <- function(p) vapply(p, function(x)
  if (is.na(x)) "" else if (x < 0.001) "***" else if (x < 0.01) "**"
  else if (x < 0.05) "*" else "", character(1))

# top annotation: mean species abundance
mean_abund <- rowMeans(log2(mtx_cpm[sp_common, , drop = FALSE] + 1), na.rm = TRUE)  # true mean log2(CPM+1)
top_anno <- HeatmapAnnotation(
  log2CPM = mean_abund,
  col = list(log2CPM = circlize::colorRamp2(range(mean_abund, na.rm = TRUE),
                                            c("white", "darkgreen"))),
  show_annotation_name = TRUE
)

cell_w   <- 6.5            # column width (mm)
cell_h   <- 4.3            # row height (mm) ~ 2/3 of previous -> shorter figure
gap_frac <- 0.20          # inter-cell gap

# ---- drug -> pathway Sankey (merges Barplot_pathway_drug_count) ----
# Each drug row is linked by a ribbon to a pathway node on the right;
# node height is proportional to the number of drugs in that pathway.
pw_drug <- setNames(pathway_vec, keep)
pw_drug[is.na(pw_drug)] <- "Unannotated"
pw_levels <- names(sort(table(pw_drug), decreasing = TRUE))
pw_count  <- table(pw_drug)[pw_levels]
frac_pw   <- as.numeric(pw_count) / sum(pw_count);  names(frac_pw) <- pw_levels
top_pw    <- cumsum(c(0, frac_pw))[seq_along(frac_pw)]; names(top_pw) <- pw_levels
node_y    <- 1 - (top_pw + frac_pw / 2)             # node centre (npc, top = 1)
pw_col    <- col_pathway
pw_col[setdiff(pw_levels, names(pw_col))] <- "grey80"

anno_sankey <- AnnotationFunction(
  fun = function(index, k, n) {
    np <- length(index)
    drugs_disp <- rownames(M_order)[index]          # displayed top -> bottom
    y_drug <- 1 - (seq_len(np) - 0.5) / np           # match heatmap row positions
    p_vec  <- pw_drug[drugs_disp]
    # straight lines linking each drug row to its pathway node
    grid.segments(x0 = 0.05, y0 = y_drug, x1 = 0.70, y1 = node_y[p_vec],
                  gp = gpar(col = adjustcolor(pw_col[p_vec], alpha.f = 0.55), lwd = 0.8))
    for (p in pw_levels) {
      grid.rect(x = 0.80, y = node_y[p], width = 0.14, height = frac_pw[p] * 0.96,
                gp = gpar(fill = pw_col[p], col = "white", lwd = 0.5))
      grid.text(pw_count[p], x = 0.80, y = node_y[p],
                gp = gpar(fontsize = 6, col = "white", fontface = "bold"))
    }
  },
  which = "row", width = unit(3.4, "cm")
)

ht_both <- Heatmap(
  M_order,
  name              = "Rs",
  col               = col_fun,
  column_title      = "Microbe – drug-response correlation   (cell: left = IDWAS | right = CaDRReS)",
  column_title_gp   = gpar(fontsize = 10),
  rect_gp           = gpar(type = "none"),
  width             = ncol(M_order) * unit(cell_w, "mm"),
  height            = nrow(M_order) * unit(cell_h, "mm"),
  row_names_gp      = gpar(fontsize = 7),
  column_names_gp   = gpar(fontsize = 8),
  show_row_names    = TRUE,
  show_column_names = TRUE,
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  top_annotation    = top_anno,
  right_annotation  = rowAnnotation(pathway_link = anno_sankey,
                                    show_annotation_name = FALSE),
  layer_fun = function(j, i, x, y, w, h, fill) {
    ri <- ComplexHeatmap::pindex(R_i, i, j)
    rc <- ComplexHeatmap::pindex(R_c, i, j)
    pi <- ComplexHeatmap::pindex(P_i, i, j)
    pc <- ComplexHeatmap::pindex(P_c, i, j)
    bw <- w * (1 - gap_frac)
    bh <- h * (1 - gap_frac)
    grid.rect(x = x - bw/4, y = y, width = bw/2, height = bh,
              gp = gpar(fill = fill_col(ri), col = NA))   # left  = IDWAS
    grid.rect(x = x + bw/4, y = y, width = bw/2, height = bh,
              gp = gpar(fill = fill_col(rc), col = NA))   # right = CaDRReS
    grid.rect(x = x, y = y, width = bw, height = bh,
              gp = gpar(fill = NA, col = "grey70", lwd = 0.4))
    grid.text(star_vec(pi), x - bw/4, y, gp = gpar(fontsize = 4, col = "black"))
    grid.text(star_vec(pc), x + bw/4, y, gp = gpar(fontsize = 4, col = "black"))
  },
  show_heatmap_legend = TRUE
)

lgd_split <- Legend(title = "Cell half", labels = c("Left – IDWAS", "Right – CaDRReS"),
                    legend_gp = gpar(fill = c("grey75", "grey45")))
lgd_pw <- Legend(title = "Pathway (node = drug count)", labels = pw_levels,
                 legend_gp = gpar(fill = pw_col[pw_levels]))

pdf(file.path(DIR_FIG, "Heatmap_drug_both.pdf"), width = 10, height = 8)
draw(ht_both,
     merge_legends          = TRUE,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     annotation_legend_list = list(lgd_split, lgd_pw))
dev.off()

####
# EF
comman_sample_ef_all <- intersect(colnames(idwas_mtx), colnames(abund_mtx_all))
cor_ef_idwas <- corr.test(t(as.matrix(idwas_mtx[, comman_sample_ef_all])),
                          t(as.matrix(abund_mtx_all["Enterococcus_faecalis", comman_sample_ef_all])),
                          method = "spearman")
cor_ef_idwas_r <- cor_ef_idwas$r
cor_ef_idwas_p <- cor_ef_idwas$p.adj
rownames(cor_ef_idwas_r) <- gsub("\\.", "-", rownames(cor_ef_idwas_r))

cor_ef_idwas_res <- data.frame(
  drug = rownames(cor_ef_idwas_r),
  rs = cor_ef_idwas_r[ , 1],
  padj = cor_ef_idwas_p[ , 1]
) %>%
  mutate(sig = case_when(rs > 0.3 & padj < 0.05 ~ "Pos",
                         rs < -0.3 & padj < 0.05 ~ "Neg",
                         TRUE ~ "NS"))

cor_ef_idwas_highlight <- cor_ef_idwas_res %>%
  arrange(rs)
cor_ef_idwas_highlight <- cor_ef_idwas_highlight[c(1:3,
                                                   (nrow(cor_ef_idwas_highlight) - 2):
                                                     nrow(cor_ef_idwas_highlight)), ]

pdf(file.path(DIR_FIG, "A_volc_EF_cor_idwas.pdf"), width = 4, height = 3.5)
ggplot(cor_ef_idwas_res, aes(x = rs, y = -log10(padj), colour = sig)) +
  ggtitle("EF correlated drugs (IDWAS)") +
  geom_point(size = 2) +
  scale_color_manual(values = c("#126CAA", "grey", "#9A342C")) +
  geom_point(cor_ef_idwas_highlight, shape = 21, color = "#FFB900FF",
             mapping = aes(x = rs, y = -log10(padj)),
             size = 2.7) +
  xlim(-0.5, 0.5) +
  ggrepel::geom_label_repel(data = cor_ef_idwas_highlight,
                            aes(x = rs, y = -log10(padj),
                                label = drug),
                            color="grey27",
                            alpha = .8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
  geom_vline(xintercept = -0.3, linetype = "dashed") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))
dev.off()

# Correlation between abundance and CaDRReS-predicted responses
comman_sample_ef_cadrres <- intersect(colnames(cadrres_pred_matrix), colnames(abund_mtx_all))
cor_ef_cadrres <- corr.test(t(as.matrix(cadrres_pred_matrix[ , comman_sample_ef_cadrres])),
                            t(as.matrix(abund_mtx_all["Enterococcus_faecalis", comman_sample_ef_cadrres])),
                            method = "spearman")
cor_ef_cadrres_r <- cor_ef_cadrres$r
cor_ef_cadrres_p <- cor_ef_cadrres$p.adj

cor_ef_cadrres_res <- data.frame(
  drug = rownames(cor_ef_cadrres_r),
  rs = cor_ef_cadrres_r[ , 1],
  padj = cor_ef_cadrres_p[ , 1]
) %>%
  mutate(sig = case_when(rs > 0.3 & padj < 0.05 ~ "Pos",
                         rs < -0.3 & padj < 0.05 ~ "Neg",
                         TRUE ~ "NS"))

cor_ef_cadrres_highlight <- cor_ef_cadrres_res %>%
  arrange(rs)
cor_ef_cadrres_highlight <- cor_ef_cadrres_highlight[c(1:4,
                                                       (nrow(cor_ef_cadrres_highlight) - 3):
                                                         nrow(cor_ef_cadrres_highlight)), ]

pdf(file.path(DIR_FIG, "A_volc_EF_cor_cadrres.pdf"), width = 4, height = 3.5)
ggplot(cor_ef_cadrres_res, aes(x = rs, y = -log10(padj), colour = sig)) +
  ggtitle("EF correlated drugs (CaDRReS)") +
  geom_point(size = 2) +
  scale_color_manual(values = c("#126CAA", "grey", "#9A342C")) +
  geom_point(cor_ef_cadrres_highlight, shape = 21, color = "#FFB900FF",
             mapping = aes(x = rs, y = -log10(padj)),
             size = 2.7) +
  ggrepel::geom_label_repel(data = cor_ef_cadrres_highlight,
                            aes(x = rs, y = -log10(padj),
                                label = drug),
                            color="grey27",
                            alpha = .8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
  geom_vline(xintercept = -0.3, linetype = "dashed") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))
dev.off()

####################
# Calculate overlaps
process_drug <- function(res, sig_type) {
  drugs <- res$drug[res$sig == sig_type]
  drugs[!grepl("\\-|[1-9]", drugs)]
}

sensitive_idwas   <- process_drug(cor_ef_idwas_res, "Neg")
sensitive_cadrres <- process_drug(cor_ef_cadrres_res, "Neg")
resistant_idwas   <- process_drug(cor_ef_idwas_res, "Pos")
resistant_cadrres <- process_drug(cor_ef_cadrres_res, "Pos")

# Venn
venn_both <- list(
  "IDWAS Sensitive"   = sensitive_idwas,
  "CaDRReS Sensitive" = sensitive_cadrres,
  "CaDRReS Resistant" = resistant_cadrres,
  "IDWAS Resistant"   = resistant_idwas
)
pdf(file.path(DIR_FIG, "B_venn.pdf"), width = 6, height = 6)
ggvenn(venn_both,
       text_size = 4,
       text_color = "white",
       show_percentage = T,
       fill_color = paletteer_d("ggsci::default_jama")[c(3, 5, 2, 4)],
       fill_alpha = 0.7,
       stroke_color = "white")
dev.off()
