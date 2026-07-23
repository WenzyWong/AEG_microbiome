####################################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 06. Microbe abundance vs. predicted drug response.
#   Panels: pRRophetic  +  CaDRReS  (CaDRReS replaces the former CARE panel)
#   Phosphoproteomics analysis has been removed (see 05/earlier steps).
#
# CaDRReS predictions are produced offline by
#   drugpred/prog/50_cadrres.py  ->  drugpred/out/cadrres_pred_all.csv
# (drugs x samples, ln IC50, higher = more resistant, i.e. same
#  convention as pRRophetic; no sign flip needed.)
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
# pRRophetic
pRRophetic_mtx <- readRDS(file.path(DIR_RDS, "drugPred_Mtx.rds"))
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

# NOTE on method naming: the drug-response prediction matrix (drugPred_Mtx.rds)
# is imputed with pRRophetic (Geeleher et al., PLoS ONE 2014; ridge regression,
# GDSC expression -> IC50; newer package: oncoPredict). It was previously named
# after the downstream imputed-drug-wide-association framework, but the actual
# predictor is pRRophetic, so all variables and figure labels below use
# "pRRophetic".
#
# The per-sample pRRophetic and CaDRReS heatmaps are merged into a single
# split-cell heatmap (Heatmap_response_both_samples.pdf) below.

rm_pRRophetic <- c()
for (i in 1:nrow(pRRophetic_mtx)) {
  if (length(na.omit(pRRophetic_mtx[i, ])) < 0.5 * ncol(pRRophetic_mtx)) {
    rm_pRRophetic <- c(rm_pRRophetic, i)
  }
}
pRRophetic_mtx <- pRRophetic_mtx[-rm_pRRophetic, ]
dim(pRRophetic_mtx)

###############
# CaDRReS  (replaces the former CARE panel)
# Load precomputed CaDRReS predictions (drugs x samples, ln IC50,
# higher = more resistant -- already aligned with pRRophetic).
cadrres_pred_matrix <- as.matrix(read.csv(
  file.path(DIR_PRED, "cadrres_pred_all.csv"), row.names = 1, check.names = FALSE))
dim(cadrres_pred_matrix)

# ---- Merged per-sample prediction heatmap (split cell: pRRophetic | CaDRReS) ----
# One heatmap of drugs x samples; every cell is split into two halves:
#   LEFT half = pRRophetic (pRRophetic) z-score, RIGHT half = CaDRReS z-score.
# Each method is z-scored per drug (across samples) so the two methods -- which
# live on different absolute scales -- share a single colour scale.
zscore_rows <- function(m) {
  z <- t(apply(m, 1, function(v) {
    s <- sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) rep(NA_real_, length(v)) else (v - mean(v, na.rm = TRUE)) / s
  }))
  dimnames(z) <- dimnames(m); z
}

pRRophetic_s <- pRRophetic_mtx
rownames(pRRophetic_s) <- gsub("\\.", "-", rownames(pRRophetic_s))
drug_s <- intersect(rownames(pRRophetic_s), rownames(cadrres_pred_matrix))
samp_s <- intersect(colnames(pRRophetic_s), colnames(cadrres_pred_matrix))
Zi_s <- zscore_rows(pRRophetic_s[drug_s, samp_s, drop = FALSE])
Zc_s <- zscore_rows(cadrres_pred_matrix[drug_s, samp_s, drop = FALSE])

M_s <- (Zi_s + Zc_s) / 2                       # ordering/clustering matrix
M_s[is.na(M_s)] <- 0

z_lim  <- as.numeric(quantile(abs(c(Zi_s, Zc_s)), 0.98, na.rm = TRUE))  # robust symmetric limit
col_z  <- colorRamp2(seq(-z_lim, z_lim, length.out = 11), rev(brewer.pal(11, "RdBu")))
fill_z <- function(v) { z <- col_z(v); z[is.na(z)] <- "grey95"; z }
gap_s  <- 0.15

ht_samp <- Heatmap(
  M_s,
  name              = "Pred. response (z)",
  col               = col_z,
  column_title      = "Per-sample predicted drug response   (cell: left = pRRophetic | right = CaDRReS)",
  column_title_gp   = gpar(fontsize = 10),
  rect_gp           = gpar(type = "none"),
  row_names_gp      = gpar(fontsize = 6),
  show_row_names    = TRUE,
  show_column_names = FALSE,
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  layer_fun = function(j, i, x, y, w, h, fill) {
    zi <- ComplexHeatmap::pindex(Zi_s, i, j)
    zc <- ComplexHeatmap::pindex(Zc_s, i, j)
    bw <- w * (1 - gap_s)
    grid.rect(x = x - bw/4, y = y, width = bw/2, height = h,
              gp = gpar(fill = fill_z(zi), col = NA))   # left  = pRRophetic
    grid.rect(x = x + bw/4, y = y, width = bw/2, height = h,
              gp = gpar(fill = fill_z(zc), col = NA))   # right = CaDRReS
  },
  show_heatmap_legend = TRUE
)
lgd_split_s <- Legend(title = "Cell half", labels = c("Left - pRRophetic", "Right - CaDRReS"),
                      legend_gp = gpar(fill = c("grey75", "grey45")))

cairo_pdf(file.path(DIR_FIG, "Heatmap_response_both_samples.pdf"), width = 9, height = 12)
draw(ht_samp, merge_legends = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right", annotation_legend_list = list(lgd_split_s))
dev.off()

# Target species
# in pRRophetic
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
common_sp_pRRophetic <- intersect(colnames(pRRophetic_mtx), colnames(abund_mtx_all))
cor_sp_pRRophetic <- corr_by_species(pRRophetic_mtx, abund_mtx_all, common_sp_pRRophetic)
saveRDS(cor_sp_pRRophetic, file.path(DIR_RDS, "Correlation_species_pRRophetic.rds"))
cor_r_pRRophetic <- cor_sp_pRRophetic$r
cor_p_pRRophetic <- cor_sp_pRRophetic$p.adj

# Drug name normalisation
rownames(cor_r_pRRophetic) <- gsub("\\.", "-", rownames(cor_r_pRRophetic))
rownames(cor_p_pRRophetic) <- gsub("\\.", "-", rownames(cor_p_pRRophetic))

drug_anno <- build_pathway_annotation(rownames(cor_r_pRRophetic), gdsc_target)
cor_r_pRRophetic <- apply_rename(cor_r_pRRophetic, drug_anno$rename_map)
cor_p_pRRophetic <- apply_rename(cor_p_pRRophetic, drug_anno$rename_map)
pathway_all <- drug_anno$pathway_all

# Snapshot the FULL (unfiltered) pRRophetic correlation matrices BEFORE any drug
# selection, for the supplementary global heatmap (Heatmap_drug_both_unfiltered.pdf).
cor_r_pRRophetic_full <- cor_r_pRRophetic
cor_p_pRRophetic_full <- cor_p_pRRophetic

# Filter
keep <- rownames(cor_r_pRRophetic)[sig_in(cor_r_pRRophetic, cor_p_pRRophetic)]

# FDA-style filter from original code: drop rows containing "-" or digits, keep those with lowercase letters
fda_filter <- function(x) {
  x <- x[!grepl("\\-|[1-9]", x)]
  x <- x[grepl("[a-z]", x)]
  x
}
keep <- fda_filter(keep)
keep <- intersect(keep, rownames(cor_r_pRRophetic))

cor_r_pRRophetic <- cor_r_pRRophetic[keep, , drop = FALSE]
cor_p_pRRophetic <- cor_p_pRRophetic[keep, , drop = FALSE]

# NOTE: the significance + readable-name (FDA-style) filter above defines the
# pRRophetic drug set that seeds the merged split-cell heatmap below; the
# stand-alone pRRophetic-only heatmap and the per-pathway drug-count barplot
# were dropped because both are superseded by the merged Heatmap_drug_both.pdf
# (which carries the pRRophetic|CaDRReS split cells and the drug -> pathway Sankey).

# in CaDRReS
common_sp_cadrres <- intersect(colnames(cadrres_pred_matrix), colnames(abund_mtx_all))
cor_sp_cadrres <- corr_by_species(cadrres_pred_matrix, abund_mtx_all, common_sp_cadrres)
saveRDS(cor_sp_cadrres, file.path(DIR_RDS, "Correlation_species_cadrres.rds"))
cor_r_cadrres <- cor_sp_cadrres$r
cor_p_cadrres <- cor_sp_cadrres$p.adj

# Drug name normalisation
rownames(cor_r_pRRophetic)   <- gsub("\\.", "-", rownames(cor_r_pRRophetic))
rownames(cor_p_pRRophetic)   <- gsub("\\.", "-", rownames(cor_p_pRRophetic))
rownames(cor_r_cadrres) <- gsub("\\.", "-", rownames(cor_r_cadrres))
rownames(cor_p_cadrres) <- gsub("\\.", "-", rownames(cor_p_cadrres))

drug_anno <- build_pathway_annotation(
  union(rownames(cor_r_pRRophetic), rownames(cor_r_cadrres)),
  gdsc_target
)
cor_r_pRRophetic   <- apply_rename(cor_r_pRRophetic, drug_anno$rename_map)
cor_p_pRRophetic   <- apply_rename(cor_p_pRRophetic, drug_anno$rename_map)
cor_r_cadrres <- apply_rename(cor_r_cadrres, drug_anno$rename_map)
cor_p_cadrres <- apply_rename(cor_p_cadrres, drug_anno$rename_map)
pathway_all <- drug_anno$pathway_all

# Snapshot the FULL (unfiltered) CaDRReS matrices + the union-based pathway table
# for the supplementary global heatmap. Also re-harmonise the pRRophetic full copy
# with the union rename map so drug names match across the two full matrices.
cor_r_cadrres_full <- cor_r_cadrres
cor_p_cadrres_full <- cor_p_cadrres
cor_r_pRRophetic_full <- apply_rename(cor_r_pRRophetic_full, drug_anno$rename_map)
cor_p_pRRophetic_full <- apply_rename(cor_p_pRRophetic_full, drug_anno$rename_map)
pathway_all_full <- drug_anno$pathway_all

# Filter: significant in either panel
common_drugs <- intersect(rownames(cor_r_pRRophetic), rownames(cor_r_cadrres))
keep <- common_drugs[
  sig_in(cor_r_pRRophetic[common_drugs, , drop = FALSE],
         cor_p_pRRophetic[common_drugs, , drop = FALSE]) |
    sig_in(cor_r_cadrres[common_drugs, , drop = FALSE],
           cor_p_cadrres[common_drugs, , drop = FALSE])
]

# Concordance filter: overall effect sign must agree across panels
common_sp <- intersect(colnames(cor_r_pRRophetic), colnames(cor_r_cadrres))

concordant <- vapply(keep, function(d) {
  r1 <- cor_r_pRRophetic[d, common_sp]
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

keep <- intersect(keep, rownames(cor_r_pRRophetic))
keep <- intersect(keep, rownames(cor_r_cadrres))

cor_r_pRRophetic   <- cor_r_pRRophetic[keep, , drop = FALSE]
cor_p_pRRophetic   <- cor_p_pRRophetic[keep, , drop = FALSE]
cor_r_cadrres <- cor_r_cadrres[keep,  , drop = FALSE]
cor_p_cadrres <- cor_p_cadrres[keep,  , drop = FALSE]

# Right pathway annotation
rownames(pathway_all) <- pathway_all$Name
pathway_vec <- pathway_all[keep, "Target.pathway"]

upw_ <- unique(pathway_vec)
is_grey_ <- is.na(upw_) | upw_ %in% c("Other", "Unannotated")
col_pathway <- setNames(rep("grey85", length(upw_)), upw_)          # Unannotated / NA
pal_pw <- as.character(paletteer::paletteer_d("ggsci::nrc_npg"))    # 10 saturated ggsci colours (grey reserved for "Other")
n_col_ <- sum(!is_grey_)
col_pathway[!is_grey_] <- if (n_col_ <= length(pal_pw)) pal_pw[seq_len(n_col_)] else colorRampPalette(pal_pw)(n_col_)
col_pathway[!is.na(upw_) & upw_ == "Other"] <- "grey70"            # "Other" -> grey

right_anno <- rowAnnotation(
  Pathway = pathway_vec,
  col = list(Pathway = col_pathway),
  show_annotation_name = TRUE
)

# ---- Merged split-cell heatmap ------------------------------------
# One heatmap of drugs x species; every cell is split into two halves:
#   LEFT  half = pRRophetic Rs,  RIGHT half = CaDRReS Rs (shared colour scale).
# Significance stars (*/**/***) drawn on each half.
sp_common <- intersect(colnames(cor_r_pRRophetic), colnames(cor_r_cadrres))
R_i <- cor_r_pRRophetic[keep, sp_common, drop = FALSE]
R_c <- cor_r_cadrres[keep, sp_common, drop = FALSE]
P_i <- cor_p_pRRophetic[keep, sp_common, drop = FALSE]
P_c <- cor_p_cadrres[keep, sp_common, drop = FALSE]

rng <- range(c(R_i, R_c), na.rm = TRUE)
col_fun <- colorRamp2(seq(-max(abs(rng)), max(abs(rng)), length.out = 11),
                      rev(brewer.pal(11, "RdBu")))
fill_col <- function(v) { z <- col_fun(v); z[is.na(z)] <- "grey95"; z }

# row/column order from the mean pattern of the two methods
M_order <- (R_i + R_c) / 2
M_order[is.na(M_order)] <- 0

star_vec <- function(p) vapply(p, function(x)
  if (is.na(x)) "" else if (x < 0.05) "*" else "", character(1))  # single star for p < 0.05 (less crowded)

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
  column_title      = "Microbe – drug-response correlation   (cell: left = pRRophetic | right = CaDRReS)",
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
              gp = gpar(fill = fill_col(ri), col = NA))   # left  = pRRophetic
    grid.rect(x = x + bw/4, y = y, width = bw/2, height = bh,
              gp = gpar(fill = fill_col(rc), col = NA))   # right = CaDRReS
    grid.rect(x = x, y = y, width = bw, height = bh,
              gp = gpar(fill = NA, col = "grey70", lwd = 0.4))
    grid.text(star_vec(pi), x - bw/4, y, gp = gpar(fontsize = 4, col = "black"))
    grid.text(star_vec(pc), x + bw/4, y, gp = gpar(fontsize = 4, col = "black"))
  },
  show_heatmap_legend = TRUE
)

lgd_split <- Legend(title = "Cell half", labels = c("Left – pRRophetic", "Right – CaDRReS"),
                    legend_gp = gpar(fill = c("grey75", "grey45")))
lgd_pw <- Legend(title = "Pathway (node = drug count)", labels = pw_levels,
                 legend_gp = gpar(fill = pw_col[pw_levels]))

cairo_pdf(file.path(DIR_FIG, "Heatmap_drug_both.pdf"), width = 10, height = 8)
draw(ht_both,
     merge_legends          = TRUE,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     annotation_legend_list = list(lgd_split, lgd_pw))
dev.off()

# ---- Supplementary: UNFILTERED global split-cell heatmap -----------------------
# Same drugs x species split-cell heatmap as Heatmap_drug_both.pdf above, but with
# NO drug selection at all: no significance filter, no cross-method concordance
# filter, no readable-name (FDA-style) filter. Every drug shared by pRRophetic and
# CaDRReS is drawn, to show the global microbe - drug-response correlation landscape.
# Built entirely from the *_full snapshots taken before filtering, so it does not
# depend on / alter the filtered objects used for the main figure.
keep_all      <- intersect(rownames(cor_r_pRRophetic_full), rownames(cor_r_cadrres_full))
sp_common_all <- intersect(colnames(cor_r_pRRophetic_full), colnames(cor_r_cadrres_full))

R_i_all <- cor_r_pRRophetic_full[keep_all, sp_common_all, drop = FALSE]
R_c_all <- cor_r_cadrres_full[keep_all,   sp_common_all, drop = FALSE]
P_i_all <- cor_p_pRRophetic_full[keep_all, sp_common_all, drop = FALSE]
P_c_all <- cor_p_cadrres_full[keep_all,   sp_common_all, drop = FALSE]

rng_all      <- range(c(R_i_all, R_c_all), na.rm = TRUE)
col_fun_all  <- colorRamp2(seq(-max(abs(rng_all)), max(abs(rng_all)), length.out = 11),
                           rev(brewer.pal(11, "RdBu")))
fill_col_all <- function(v) { z <- col_fun_all(v); z[is.na(z)] <- "grey95"; z }

M_order_all <- (R_i_all + R_c_all) / 2
M_order_all[is.na(M_order_all)] <- 0

# pathway annotation for ALL drugs (colour bar; Sankey dropped -- too dense here)
pathway_vec_all <- pathway_all_full[keep_all, "Target.pathway"]
pathway_vec_all[is.na(pathway_vec_all)] <- "Unannotated"
upw_all     <- unique(pathway_vec_all)
is_grey_all <- upw_all %in% c("Other", "Unannotated")
col_pathway_all <- setNames(rep("grey85", length(upw_all)), upw_all)
n_col_all <- sum(!is_grey_all)
col_pathway_all[!is_grey_all] <- if (n_col_all <= length(pal_pw)) pal_pw[seq_len(n_col_all)] else colorRampPalette(pal_pw)(n_col_all)
col_pathway_all[upw_all == "Other"] <- "grey70"

# top annotation: mean species abundance (identical convention to the main heatmap)
mean_abund_all <- rowMeans(log2(mtx_cpm[sp_common_all, , drop = FALSE] + 1), na.rm = TRUE)
top_anno_all <- HeatmapAnnotation(
  log2CPM = mean_abund_all,
  col = list(log2CPM = circlize::colorRamp2(range(mean_abund_all, na.rm = TRUE),
                                            c("white", "darkgreen"))),
  show_annotation_name = TRUE)

right_anno_all <- rowAnnotation(
  Pathway = pathway_vec_all,
  col = list(Pathway = col_pathway_all),
  show_annotation_name = TRUE)

cell_w_all <- 6.5     # column width (mm)
cell_h_all <- 2.6     # row height (mm): shorter, since all drugs are shown

ht_both_all <- Heatmap(
  M_order_all,
  name              = "Rs",
  col               = col_fun_all,
  column_title      = "Global microbe - drug-response correlation (unfiltered)   (cell: left = pRRophetic | right = CaDRReS)",
  column_title_gp   = gpar(fontsize = 10),
  rect_gp           = gpar(type = "none"),
  width             = ncol(M_order_all) * unit(cell_w_all, "mm"),
  height            = nrow(M_order_all) * unit(cell_h_all, "mm"),
  row_names_gp      = gpar(fontsize = 4),
  column_names_gp   = gpar(fontsize = 8),
  show_row_names    = TRUE,
  show_column_names = TRUE,
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  top_annotation    = top_anno_all,
  right_annotation  = right_anno_all,
  layer_fun = function(j, i, x, y, w, h, fill) {
    ri <- ComplexHeatmap::pindex(R_i_all, i, j)
    rc <- ComplexHeatmap::pindex(R_c_all, i, j)
    pi <- ComplexHeatmap::pindex(P_i_all, i, j)
    pc <- ComplexHeatmap::pindex(P_c_all, i, j)
    bw <- w * (1 - gap_frac)
    bh <- h * (1 - gap_frac)
    grid.rect(x = x - bw/4, y = y, width = bw/2, height = bh,
              gp = gpar(fill = fill_col_all(ri), col = NA))   # left  = pRRophetic
    grid.rect(x = x + bw/4, y = y, width = bw/2, height = bh,
              gp = gpar(fill = fill_col_all(rc), col = NA))   # right = CaDRReS
    grid.rect(x = x, y = y, width = bw, height = bh,
              gp = gpar(fill = NA, col = "grey70", lwd = 0.3))
    grid.text(star_vec(pi), x - bw/4, y, gp = gpar(fontsize = 3, col = "black"))
    grid.text(star_vec(pc), x + bw/4, y, gp = gpar(fontsize = 3, col = "black"))
  },
  show_heatmap_legend = TRUE)

lgd_split_all <- Legend(title = "Cell half", labels = c("Left - pRRophetic", "Right - CaDRReS"),
                        legend_gp = gpar(fill = c("grey75", "grey45")))

n_drug_all <- nrow(M_order_all)
cairo_pdf(file.path(DIR_FIG, "Heatmap_drug_both_unfiltered.pdf"),
          width = 11, height = max(8, 0.12 * n_drug_all + 3))
draw(ht_both_all,
     merge_legends          = TRUE,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     annotation_legend_list = list(lgd_split_all))
dev.off()
cat("Supplementary unfiltered heatmap drugs:", n_drug_all,
    " species:", ncol(M_order_all), "\n")

# ---- Cross-cancer (ESCC + STAD) annotations on the AEG global heatmap ------------
# Combines the content of Fig_drug_direction_heatmap.pdf with the unfiltered AEG
# heatmap above: ESCC and STAD are attached to the SAME drug rows (keep_all) and the
# SAME 15 target species as drug-row annotations beside the main AEG heatmap --
#   * two mini heatmaps (per cancer, columns = pRRophetic | CaDRReS) of the per-drug
#     mean signed Spearman r across the target species (own diverging scale, since
#     these species-averaged Rs are smaller than the AEG per-cell values);
#   * two 100%-stacked bars = fraction of (target-species x method) cells whose
#     ESCC / STAD correlation sign matches AEG (replicated) vs discordant.
# ESCC/STAD drug-response predictions come from drugpred/out, microbe CPM matrices
# from drug_crosscancer/species; corr_by_species() is reused unchanged.
DIR_SP_CC <- "/data/yzwang/project/AEG_seiri/drug_crosscancer/species/"
rd_cc <- function(f) as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
cc_pred <- list(
  ESCC = list(pRRophetic = rd_cc(file.path(DIR_PRED, "oncopredict_pred_escc.csv")),
              CaDRReS    = rd_cc(file.path(DIR_PRED, "cadrres_pred_escc.csv"))),
  STAD = list(pRRophetic = rd_cc(file.path(DIR_PRED, "oncopredict_pred_gc.csv")),
              CaDRReS    = rd_cc(file.path(DIR_PRED, "cadrres_pred_gc.csv"))))
cc_ab <- list(
  ESCC = readRDS(file.path(DIR_SP_CC, "sESCC_CPM.rds")),
  STAD = rd_cc(file.path(DIR_SP_CC, "sGC_CPM.csv")))

aeg_R_cc <- list(pRRophetic = R_i_all, CaDRReS = R_c_all)   # keep_all x sp_common_all
cc_R <- list()
for (cx in names(cc_pred)) for (mt in names(cc_pred[[cx]])) {
  ab     <- cc_ab[[cx]]
  sp_use <- intersect(sp_common_all, rownames(ab))
  samp   <- intersect(colnames(cc_pred[[cx]][[mt]]), colnames(ab))
  cs <- corr_by_species(cc_pred[[cx]][[mt]], ab[sp_use, , drop = FALSE], samp)$r
  M  <- matrix(NA_real_, length(keep_all), length(sp_common_all),
               dimnames = list(keep_all, sp_common_all))
  dd <- intersect(keep_all, rownames(cs)); ss <- intersect(sp_common_all, colnames(cs))
  M[dd, ss] <- cs[dd, ss]
  cc_R[[paste(cx, mt, sep = "_")]] <- M
}

mean_r_cc <- function(M) rowMeans(M, na.rm = TRUE)
escc_mat <- cbind(pRRophetic = mean_r_cc(cc_R$ESCC_pRRophetic), CaDRReS = mean_r_cc(cc_R$ESCC_CaDRReS))
stad_mat <- cbind(pRRophetic = mean_r_cc(cc_R$STAD_pRRophetic), CaDRReS = mean_r_cc(cc_R$STAD_CaDRReS))

# per-drug replication %: sign agreement vs AEG across target species AND both methods
repl_pct_cc <- function(cx) vapply(keep_all, function(d) {
  ag <- logical(0)
  for (mt in c("pRRophetic", "CaDRReS")) {
    a <- aeg_R_cc[[mt]][d, ]; b <- cc_R[[paste(cx, mt, sep = "_")]][d, ]
    ok <- is.finite(a) & is.finite(b) & a != 0 & b != 0
    if (any(ok)) ag <- c(ag, sign(a[ok]) == sign(b[ok]))
  }
  if (length(ag) == 0) NA_real_ else 100 * mean(ag)
}, numeric(1))
escc_repl <- repl_pct_cc("ESCC"); stad_repl <- repl_pct_cc("STAD")
escc_bar  <- cbind(ifelse(is.na(escc_repl), 0, escc_repl), ifelse(is.na(escc_repl), 0, 100 - escc_repl))
stad_bar  <- cbind(ifelse(is.na(stad_repl), 0, stad_repl), ifelse(is.na(stad_repl), 0, 100 - stad_repl))

esc_col <- c("#74ACCE", "#E8830C"); sta_col <- c("#8B87BA", "#F2C300")   # replicated | discordant
mx_cc <- max(abs(c(escc_mat, stad_mat)), na.rm = TRUE)
col_fun_cc <- colorRamp2(c(-mx_cc, 0, mx_cc), c("#2166AC", "white", "#B2182B"))

# Pathway colours consistent with the filtered Heatmap_drug_both.pdf: reuse that
# figure's col_pathway for overlapping pathways; the unfiltered set has more
# pathways, so the non-overlapping ones keep col_pathway_all's colours.
col_pathway_cc <- col_pathway_all
ov_pw <- intersect(names(col_pathway_cc), names(col_pathway))
col_pathway_cc[ov_pw] <- col_pathway[ov_pw]

right_anno_cc <- rowAnnotation(
  Pathway = pathway_vec_all,
  col     = list(Pathway = col_pathway_cc),
  ESCC = anno_simple(escc_mat, col = col_fun_cc, gp = gpar(col = "grey90")),
  STAD = anno_simple(stad_mat, col = col_fun_cc, gp = gpar(col = "grey90")),
  `ESCC repl` = anno_barplot(escc_bar, gp = gpar(fill = esc_col, col = NA),
                             bar_width = 0.9, ylim = c(0, 100), width = unit(1.0, "cm")),
  `STAD repl` = anno_barplot(stad_bar, gp = gpar(fill = sta_col, col = NA),
                             bar_width = 0.9, ylim = c(0, 100), width = unit(1.0, "cm")),
  show_annotation_name = TRUE, annotation_name_gp = gpar(fontsize = 7),
  gap = unit(1, "mm"))

ht_both_all_cc <- Heatmap(
  M_order_all,
  name              = "Rs",
  col               = col_fun_all,
  column_title      = "AEG global heatmap with ESCC / STAD cross-cancer annotations   (cell: left = pRRophetic | right = CaDRReS)",
  column_title_gp   = gpar(fontsize = 10),
  rect_gp           = gpar(type = "none"),
  width             = ncol(M_order_all) * unit(cell_w_all * 0.75, "mm"),   # 75% of previous width
  height            = nrow(M_order_all) * unit(cell_h_all * 0.6, "mm"),   # 60% of previous height
  row_names_gp      = gpar(fontsize = 4),
  column_names_gp   = gpar(fontsize = 8),
  column_names_rot  = 45,
  show_row_names    = TRUE,
  show_column_names = TRUE,
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  top_annotation    = top_anno_all,
  right_annotation  = right_anno_cc,
  layer_fun = function(j, i, x, y, w, h, fill) {
    ri <- ComplexHeatmap::pindex(R_i_all, i, j)
    rc <- ComplexHeatmap::pindex(R_c_all, i, j)
    pi <- ComplexHeatmap::pindex(P_i_all, i, j)
    pc <- ComplexHeatmap::pindex(P_c_all, i, j)
    bw <- w * (1 - gap_frac); bh <- h * (1 - gap_frac)
    grid.rect(x - bw/4, y, bw/2, bh, gp = gpar(fill = fill_col_all(ri), col = NA))   # left  = pRRophetic
    grid.rect(x + bw/4, y, bw/2, bh, gp = gpar(fill = fill_col_all(rc), col = NA))   # right = CaDRReS
    grid.rect(x, y, bw, bh, gp = gpar(fill = NA, col = "grey70", lwd = 0.3))
    grid.text(star_vec(pi), x - bw/4, y, gp = gpar(fontsize = 3, col = "white"))
    grid.text(star_vec(pc), x + bw/4, y, gp = gpar(fontsize = 3, col = "white"))
  },
  show_heatmap_legend = TRUE)

lgd_cc_scale <- Legend(col_fun = col_fun_cc, title = "ESCC/STAD\nmean Rs")
lgd_repl     <- Legend(labels = c("ESCC replicated", "ESCC discordant",
                                  "STAD replicated", "STAD discordant"),
                       legend_gp = gpar(fill = c(esc_col, sta_col)),
                       title = "replication bars", ncol = 2)

cairo_pdf(file.path(DIR_FIG, "Heatmap_drug_both_unfiltered_crosscancer.pdf"),
          width = 9.75, height = max(5, 0.072 * n_drug_all + 1.8))   # 75% width, 60% height
draw(ht_both_all_cc,
     merge_legends          = TRUE,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     annotation_legend_list = list(lgd_split_all, lgd_cc_scale, lgd_repl))
dev.off()
cat("Cross-cancer combined heatmap: drugs", n_drug_all,
    " | ESCC repl", round(mean(escc_repl, na.rm = TRUE)),
    "%  STAD repl", round(mean(stad_repl, na.rm = TRUE)), "%\n")


# ---- Per-microbe drug-sensitivity volcano (median-split, both methods) ------
# For each microbe shown in the drug heatmap and each method, split the samples
# by that microbe's MEDIAN abundance into high/low groups; for every drug compare
# the predicted response between groups (Wilcoxon rank-sum, raw p).
#   Each drug's predicted values are first z-scored across samples (removes the
#   per-drug scale, so sensdiff is in SD units and comparable across drugs and
#   the two methods). The Wilcoxon p is unchanged by this monotonic transform.
#   sensdiff = median(high z) - median(low z).  Predicted value LOWER = more
#   sensitive, so sensdiff < 0  =>  high-abundance group is MORE SENSITIVE.
# Colours follow the heatmap convention (red = resistant [sensdiff > 0, right],
# blue = sensitive [sensdiff < 0, left]); no sensitivity-difference cutoff.
# The "N" (top-left of each panel) counts the drugs the high-abundance group
# is MORE SENSITIVE to (sensdiff < 0 & p < 0.05), mirroring the reference figure.
volc_methods <- list(
  pRRophetic = `rownames<-`(pRRophetic_mtx, gsub("\\.", "-", rownames(pRRophetic_mtx))),
  CaDRReS    = cadrres_pred_matrix
)

volc_list <- list()
for (meth in names(volc_methods)) {
  M    <- volc_methods[[meth]]
  samp <- intersect(colnames(M), colnames(abund_mtx_all))
  Mz   <- zscore_rows(M[, samp, drop = FALSE])   # per-drug z across samples: removes
                                                 # per-drug scale so sensdiff is in SD
                                                 # units and comparable across drugs/methods
  for (sp in rownames(abund_mtx_all)) {
    a   <- as.numeric(abund_mtx_all[sp, samp])
    med <- median(a, na.rm = TRUE)
    hi  <- samp[which(a >  med)]
    lo  <- samp[which(a <= med)]
    for (d in rownames(Mz)) {
      vh <- as.numeric(Mz[d, hi]); vh <- vh[is.finite(vh)]
      vl <- as.numeric(Mz[d, lo]); vl <- vl[is.finite(vl)]
      if (length(vh) < 3 || length(vl) < 3) next
      p <- tryCatch(suppressWarnings(wilcox.test(vh, vl)$p.value),
                    error = function(e) NA_real_)
      volc_list[[length(volc_list) + 1L]] <- data.frame(
        method = meth, microbe = sp, drug = d,
        sensdiff = median(vh) - median(vl), p = p,
        stringsAsFactors = FALSE)
    }
  }
}
volc <- do.call(rbind, volc_list)
volc <- volc[is.finite(volc$p) & is.finite(volc$sensdiff), ]
volc$microbe <- factor(volc$microbe, levels = rownames(abund_mtx_all))
volc$dir <- with(volc, ifelse(p < 0.05 & sensdiff > 0, "Resistant",
                       ifelse(p < 0.05 & sensdiff < 0, "Sensitive", "NS")))
volc$dir <- factor(volc$dir, levels = c("Resistant", "Sensitive", "NS"))

# Per-method, per-microbe N: number of drugs the high-abundance group is more
# SENSITIVE to (dir == "Sensitive") and more RESISTANT to (dir == "Resistant").
# Both directions are annotated in each facet.
volc_N <- aggregate(dir ~ method + microbe, data = volc,
                    FUN = function(x) sum(x == "Sensitive"))
names(volc_N)[3] <- "N"
volc_NR <- aggregate(dir ~ method + microbe, data = volc,
                     FUN = function(x) sum(x == "Resistant"))
names(volc_NR)[3] <- "N"

vol_red        <- brewer.pal(11, "RdBu")[2]    # pRRophetic resistant (saturated red)
vol_blue       <- brewer.pal(11, "RdBu")[10]   # pRRophetic sensitive (saturated blue)
vol_red_light  <- brewer.pal(11, "RdBu")[4]    # CaDRReS   resistant (light red)
vol_blue_light <- brewer.pal(11, "RdBu")[8]    # CaDRReS   sensitive (light blue)

# Merged layout: BOTH methods overlaid in ONE volcano per microbe, microbes laid
# out horizontally in 3 rows (facet_wrap(~ microbe, nrow = 3)). Colour encodes
# method x direction: pRRophetic keeps the original saturated red/blue, while
# CaDRReS's significant points use light red / light blue so the two methods stay
# distinguishable; non-significant points of either method share a single grey.
# (Both methods' sensdiff is in per-drug z SD units, so they share one x-axis.)
volc$method <- factor(volc$method, levels = c("pRRophetic", "CaDRReS"))
volc$grp <- "NS"
volc$grp[volc$method == "pRRophetic" & volc$dir == "Resistant"] <- "pRRophetic Resistant"
volc$grp[volc$method == "pRRophetic" & volc$dir == "Sensitive"] <- "pRRophetic Sensitive"
volc$grp[volc$method == "CaDRReS"    & volc$dir == "Resistant"] <- "CaDRReS Resistant"
volc$grp[volc$method == "CaDRReS"    & volc$dir == "Sensitive"] <- "CaDRReS Sensitive"
grp_levels <- c("pRRophetic Resistant", "pRRophetic Sensitive",
                "CaDRReS Resistant",    "CaDRReS Sensitive", "NS")
volc$grp <- factor(volc$grp, levels = grp_levels)
vol_cols <- c("pRRophetic Resistant" = vol_red,
              "pRRophetic Sensitive" = vol_blue,
              "CaDRReS Resistant"    = vol_red_light,
              "CaDRReS Sensitive"    = vol_blue_light,
              "NS"                   = "grey78")

# draw order: NS at the bottom, then CaDRReS, pRRophetic (saturated) on top
draw_rank <- ifelse(volc$grp == "NS", 0L,
              ifelse(volc$method == "CaDRReS", 1L, 2L))
volc <- volc[order(draw_rank), ]

# N labels per facet: SENSITIVE counts (blue, top-left) and RESISTANT counts
# (red, top-right); saturated colour = pRRophetic, light colour = CaDRReS.
mk_Nlab <- function(df) {
  df$method <- factor(df$method, levels = c("pRRophetic", "CaDRReS"))
  df$lab    <- paste0(df$method, ": N = ", df$N)
  df
}
volc_N  <- mk_Nlab(volc_N)    # sensitive
volc_NR <- mk_Nlab(volc_NR)   # resistant
volc_Ns_i <- volc_N [volc_N $method == "pRRophetic", ]  # sensitive, pRRophetic (blue)
volc_Ns_c <- volc_N [volc_N $method == "CaDRReS",    ]  # sensitive, CaDRReS   (light blue)
volc_Nr_i <- volc_NR[volc_NR$method == "pRRophetic", ]  # resistant, pRRophetic (red)
volc_Nr_c <- volc_NR[volc_NR$method == "CaDRReS",    ]  # resistant, CaDRReS   (light red)

p_volc <- ggplot(volc, aes(sensdiff, -log10(p))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_point(aes(colour = grp), size = 0.9, alpha = 0.8) +
  geom_text(data = volc_Ns_i, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            colour = vol_blue,       hjust = -0.1, vjust = 1.4, size = 2.3) +
  geom_text(data = volc_Ns_c, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            colour = vol_blue_light, hjust = -0.1, vjust = 2.9, size = 2.3) +
  geom_text(data = volc_Nr_i, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
            colour = vol_red,        hjust = 1.1,  vjust = 1.4, size = 2.3) +
  geom_text(data = volc_Nr_c, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
            colour = vol_red_light,  hjust = 1.1,  vjust = 2.9, size = 2.3) +
  scale_colour_manual(values = vol_cols, name = "Method & direction",
                      breaks = grp_levels[1:4]) +
  facet_wrap(~ microbe, nrow = 3, scales = "free") +
  labs(x = "Sensitivity difference (median high - low, per-drug z)",
       y = "-Log10(p.value)") +
  theme_bw(base_size = 8) +
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        strip.text       = element_text(size = 7.5, face = "italic"),
        strip.background  = element_rect(fill = "grey92", colour = NA),
        legend.position  = "right") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2.6, alpha = 1)))

n_micro <- nlevels(volc$microbe)
n_col   <- ceiling(n_micro / 3)
# Render with the base pdf() device (NOT cairo_pdf) so every axis number / text
# label is written as its own text block -> an independent, editable text object in
# Illustrator. cairo_pdf batches the glyphs so individual tick numbers cannot be
# selected separately. base pdf() alone leaves the fonts non-embedded ("Custom"
# encoding), so embed them afterwards with Ghostscript (embedFonts) -> independent
# text objects AND embedded fonts, avoiding the broken-text issue in Illustrator.
volc_pdf <- file.path(DIR_FIG, "Volcano_sensdiff_microbe_both.pdf")
pdf(volc_pdf, width = max(8, 2.1 * n_col + 1.4), height = 9, useDingbats = FALSE)
print(p_volc)
dev.off()
# -dPDFSETTINGS=/prepress is required for gs to actually embed the standard-14
# Helvetica (subsetted Type 1C); text stays real, editable, and per-label independent.
embedFonts(volc_pdf, options = "-dEmbedAllFonts=true -dPDFSETTINGS=/prepress")

####
# EF
comman_sample_ef_all <- intersect(colnames(pRRophetic_mtx), colnames(abund_mtx_all))
cor_ef_pRRophetic <- corr.test(t(as.matrix(pRRophetic_mtx[, comman_sample_ef_all])),
                          t(as.matrix(abund_mtx_all["Enterococcus_faecalis", comman_sample_ef_all])),
                          method = "spearman")
cor_ef_pRRophetic_r <- cor_ef_pRRophetic$r
cor_ef_pRRophetic_p <- cor_ef_pRRophetic$p.adj
rownames(cor_ef_pRRophetic_r) <- gsub("\\.", "-", rownames(cor_ef_pRRophetic_r))

cor_ef_pRRophetic_res <- data.frame(
  drug = rownames(cor_ef_pRRophetic_r),
  rs = cor_ef_pRRophetic_r[ , 1],
  padj = cor_ef_pRRophetic_p[ , 1]
) %>%
  mutate(sig = case_when(rs > 0.3 & padj < 0.05 ~ "Pos",
                         rs < -0.3 & padj < 0.05 ~ "Neg",
                         TRUE ~ "NS"))

# (EF volcano plots removed; replaced by the two-method dot-plot below.)

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

# ---- EF vs. predicted drug sensitivity dot-plot (pRRophetic | CaDRReS) ----
# Two method rows; each drug a column. Dot size = -log10(Padj); colour = Rs of
# EF abundance vs. predicted response (blue = negative Rs => EF associated with
# SENSITIVITY, red = positive Rs => resistance). Drugs whose direction agrees
# between the two methods (same Rs sign, significant in both) are boxed.
ef_anno <- build_pathway_annotation(
  union(cor_ef_pRRophetic_res$drug, cor_ef_cadrres_res$drug), gdsc_target)
ren_ef    <- ef_anno$rename_map
harmonise <- function(d) ifelse(d %in% names(ren_ef), ren_ef[d], d)

ef_i <- transform(cor_ef_pRRophetic_res,   drug = harmonise(drug))
ef_c <- transform(cor_ef_cadrres_res, drug = harmonise(drug))

# common drugs (present in both methods), readable names only
wide <- merge(ef_i[, c("drug", "rs", "padj", "sig")],
              ef_c[, c("drug", "rs", "padj", "sig")],
              by = "drug", suffixes = c("_i", "_c"))
wide <- wide[grepl("[a-z]", wide$drug), ]
wide <- wide[!duplicated(wide$drug), ]

# show only drugs significant (|rs|>0.3 & padj<0.05) in BOTH methods
show <- wide[wide$sig_i != "NS" & wide$sig_c != "NS", ]

# direction-concordant drugs: significant in BOTH, same sign of Rs
concord <- show$drug[show$sig_i != "NS" & show$sig_c != "NS" &
                     sign(show$rs_i) == sign(show$rs_c)]

ef_dot <- rbind(
  data.frame(drug = show$drug, rs = show$rs_i, padj = show$padj_i, method = "pRRophetic"),
  data.frame(drug = show$drug, rs = show$rs_c, padj = show$padj_c, method = "CaDRReS")
)
# order drugs by mean Rs: positive/resistant (red) left -> negative/sensitive (blue) right
lev <- names(sort(tapply(ef_dot$rs, ef_dot$drug, mean, na.rm = TRUE), decreasing = TRUE))
ef_dot$drug   <- factor(ef_dot$drug, levels = lev)
ef_dot$method <- factor(ef_dot$method, levels = c("CaDRReS", "pRRophetic"))  # pRRophetic drawn on top

rs_lim <- max(abs(ef_dot$rs), na.rm = TRUE)
p_dot <- ggplot(ef_dot, aes(x = drug, y = method)) +
  geom_point(aes(size = -log10(pmax(padj, 1e-10)), fill = rs),
             shape = 21, colour = "grey30") +
  scale_fill_gradientn(colours = rev(brewer.pal(11, "RdBu")),
                       limits = c(-rs_lim, rs_lim), name = "Rs") +
  scale_size_continuous(range = c(0.5, 6), name = "-log10(Padj)") +
  labs(x = NULL, y = NULL,
       title = "EF abundance vs. predicted drug sensitivity",
       subtitle = "drugs significant in both methods") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = 1),
        axis.text.y = element_text(colour = 1),
        panel.grid.minor = element_blank())
cairo_pdf(file.path(DIR_FIG, "A_dot_EF_cor_both.pdf"),
    width = max(6, 0.28 * length(lev) + 2), height = 3)
print(p_dot)
dev.off()

####################
# Calculate overlaps
# Reuse the SAME harmonised drug names (ef_i / ef_c) and readable-name filter as
# the dot-plot above, so the Venn's "significant in both" counts match the boxed
# direction-concordant drugs (previously they disagreed: the old process_drug
# used raw, un-harmonised names and dropped every name with a digit/hyphen).
readable <- function(d) unique(d[grepl("[a-z]", d)])
sensitive_pRRophetic   <- readable(ef_i$drug[ef_i$sig == "Neg"])
sensitive_cadrres <- readable(ef_c$drug[ef_c$sig == "Neg"])
resistant_pRRophetic   <- readable(ef_i$drug[ef_i$sig == "Pos"])
resistant_cadrres <- readable(ef_c$drug[ef_c$sig == "Pos"])

# Venn
venn_both <- list(
  "pRRophetic Sensitive"   = sensitive_pRRophetic,
  "CaDRReS Sensitive" = sensitive_cadrres,
  "CaDRReS Resistant" = resistant_cadrres,
  "pRRophetic Resistant"   = resistant_pRRophetic
)
cairo_pdf(file.path(DIR_FIG, "B_venn.pdf"), width = 6, height = 6)
ggvenn(venn_both,
       text_size = 4,
       text_color = "white",
       show_percentage = T,
       fill_color = paletteer_d("ggsci::default_jama")[c(3, 5, 2, 4)],
       fill_alpha = 0.7,
       stroke_color = "white")
dev.off()
