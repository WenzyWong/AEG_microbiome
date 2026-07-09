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

pdf(file.path(DIR_FIG, "Heatmap_response_both_samples.pdf"), width = 9, height = 12)
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

pdf(file.path(DIR_FIG, "Heatmap_drug_both.pdf"), width = 10, height = 8)
draw(ht_both,
     merge_legends          = TRUE,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     annotation_legend_list = list(lgd_split, lgd_pw))
dev.off()

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

# N = drugs the high-abundance group is more sensitive to (red label, per panel)
volc_N <- aggregate(dir ~ method + microbe, data = volc,
                    FUN = function(x) sum(x == "Sensitive"))
names(volc_N)[3] <- "N"
volc_N$lab <- paste0("N = ", volc_N$N)

vol_red  <- brewer.pal(11, "RdBu")[2]    # red  = resistant
vol_blue <- brewer.pal(11, "RdBu")[10]   # blue = sensitive
vol_cols <- c(Resistant = vol_red, Sensitive = vol_blue, NS = "grey78")

# Tall layout: one ROW per microbe, one COLUMN per method (2 methods).
# facet_grid(microbe ~ method, scales = "free") frees x per method-column (the
# two methods live on different SD-unit ranges) and y per microbe-row.
volc$method   <- factor(volc$method,   levels = c("pRRophetic", "CaDRReS"))
volc_N$method <- factor(volc_N$method, levels = c("pRRophetic", "CaDRReS"))

p_volc <- ggplot(volc, aes(sensdiff, -log10(p))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_point(aes(colour = dir), size = 0.9, alpha = 0.8) +
  geom_text(data = volc_N, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
            colour = "black", hjust = -0.15, vjust = 1.5, size = 2.6) +
  scale_colour_manual(values = vol_cols, name = NULL,
                      breaks = c("Resistant", "Sensitive")) +
  facet_grid(microbe ~ method, scales = "free") +
  labs(x = "Sensitivity difference (median high - low, per-drug z)",
       y = "-Log10(p.value)",
       title = "Per-microbe predicted drug sensitivity by abundance (median split)") +
  theme_bw(base_size = 8) +
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        strip.text.x     = element_text(size = 8, face = "bold"),
        strip.text.y     = element_text(size = 7, face = "italic", angle = 0),
        strip.background  = element_rect(fill = "grey92", colour = NA),
        legend.position  = "bottom")

n_micro <- nlevels(volc$microbe)
ggsave(file.path(DIR_FIG, "Volcano_sensdiff_microbe_both.pdf"),
       plot = p_volc,
       width  = 7,
       height = max(6, 1.7 * n_micro + 1),
       limitsize = FALSE)

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
pdf(file.path(DIR_FIG, "A_dot_EF_cor_both.pdf"),
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
pdf(file.path(DIR_FIG, "B_venn.pdf"), width = 6, height = 6)
ggvenn(venn_both,
       text_size = 4,
       text_color = "white",
       show_percentage = T,
       fill_color = paletteer_d("ggsci::default_jama")[c(3, 5, 2, 4)],
       fill_alpha = 0.7,
       stroke_color = "white")
dev.off()
