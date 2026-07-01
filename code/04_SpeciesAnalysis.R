#######################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 04. Analyse AEG-specific microbiome
#
#######################################################
set.seed(42)
library(vegan)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(survival)
library(survminer)
library(paletteer)
library(circlize)
library(ComplexHeatmap)
library(gridExtra)
library(tools)
library(stringr)
# Network analysis requirements
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(tidyfst)
library(pulsar)
library(ggClusterNet) # other requirements: ggraph, tidyfst
library(ggraph)
library(RColorBrewer)
# Regression
library(compositions)
library(glmnet)
# Differential
library(Maaslin2)
library(betapart)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/F2_species/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_TOOL <- "/data/yzwang/git_project/AEG_microbiome/utils/"

##############
# Data loading
mtx_gcpm <- readRDS(file.path(DIR_RDS, "gAEG_CPM_RNA_FiltMyco.rds"))
mtx_cpm <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
mtx_count <- readRDS(file.path(DIR_RDS, "sAEG_Count_RNA_FiltMyco.rds"))
# Clinical information (without filtering)
clinical <- readxl::read_excel(file.path(DIR_TAB, "AEG_clinical.xlsx"))

##############################
# Circle plot of genus-species
# Allocate positions for species
g_abund <- sort(apply(mtx_gcpm, MARGIN = 1, FUN = mean) / 1e+4, decreasing = T)
g_abund <- g_abund[1:20]
# Collect species rows that belong to any top-abundant genus
index_sp <- unlist(lapply(names(g_abund),
                          function(g) grep(g, rownames(mtx_cpm), fixed = TRUE)))

# Mean abundance (%) of each top-genus species; computed once and reused below
abund_sp_top <- apply(mtx_cpm[index_sp, ] / 1e+4, MARGIN = 1, FUN = mean) %>%
  data.frame(
    # Sanitise species names so they can be used as R formula terms
    Species = gsub("[[:punct:]]", "_", names(.)),
    Abundance = as.numeric(.),
    row.names = names(.)
  ) # Species within top 20 genera, regardless of their absolute abundance
# Notice: rownames(abund_sp_top) differs from abund_sp_top$Species.
# Pay attention to which one is used in downstream steps.

abund_sp_top <- abund_sp_top[abund_sp_top$Abundance > 0.01, ]

draw_hr <- data.frame(
  sample = paste0("C", clinical$No.),
  time = clinical$`Survival time（Month）`,
  state = clinical$`Status（1=Dead, 0=Alive）`
) %>%
  filter(sample %in% colnames(mtx_cpm))

draw_hr <- cbind(
  draw_hr,
  setNames(
    as.data.frame(t(log2(mtx_cpm[rownames(abund_sp_top), draw_hr$sample, drop = FALSE] + 1))),
    abund_sp_top$Species
  )
)

# Calculating HR
source(file.path(DIR_TOOL, "hr_calc.R"))
hr_res <- hr_calc(abund_sp_top$Species, draw_hr, 1.1, 0.9, 0.05)
# The warnings exist because all the tumour CPM values of these species equals to 0

abund_sp_top$Surv.HR <- hr_res$HR
abund_sp_top$Surv.P <- hr_res$Pvalue
abund_sp_top$Surv.Risk <- hr_res$Risk
abund_sp_top$Surv.HR[abund_sp_top$Surv.HR > 3] <- 3
# Removing the species with 0 count values in tumour samples, causing HR equals to NA
abund_sp_top <- na.omit(abund_sp_top)

# Number of identified species within each abundant genus
len_each_g <- vapply(names(g_abund),
                     function(g) length(grep(g, abund_sp_top$Species, fixed = TRUE)),
                     integer(1))

# Allocate genus tag and circular x-coordinate for each species
abund_sp_top$Genus <- rep(names(g_abund), times = len_each_g)
abund_sp_top$X     <- unlist(lapply(len_each_g,
                                    function(n) seq(0, 1, length.out = n)))

# --- Differential abundance via MaAsLin2 (paired tumour-normal) ----------------
# Run once here so the circos differential annotation and the DA figures below
# share the same model (replaces the earlier Wilcoxon test for consistency).
rel_ab_min <- 1e-4   # 0.01% mean relative abundance threshold

paired_ids <- intersect(
  gsub("^C", "", grep("^C", colnames(mtx_count), value = TRUE)),
  gsub("^N", "", grep("^N", colnames(mtx_count), value = TRUE))
)
paired_samples <- c(paste0("C", paired_ids), paste0("N", paired_ids))
meta_paired <- data.frame(
  sample = paired_samples,
  patient = gsub("^[CN]", "", paired_samples),
  group = factor(
    ifelse(grepl("^C", paired_samples), "Tumour", "Normal"),
    levels = c("Normal", "Tumour")
  ),
  row.names = paired_samples
)

maaslin_dir <- file.path(DIR_RES, "MaAsLin2_paired_TN")
dir.create(maaslin_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
maaslin_res <- Maaslin2(
  input_data = as.data.frame(mtx_count[, paired_samples, drop = FALSE]),
  input_metadata = meta_paired,
  output = maaslin_dir,
  fixed_effects = "group",
  random_effects = "patient",
  reference = "group,Normal",
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  min_prevalence = 0.10,
  min_abundance = 0.0,
  max_significance = 0.05,
  plot_heatmap = FALSE,
  plot_scatter = FALSE,
  cores = 4
)
saveRDS(maaslin_res, file.path(DIR_RDS, "MaAsLin2_paired_TN.rds"))

maaslin_tab <- maaslin_res$results %>%
  filter(metadata == "group") %>%
  transmute(
    taxon = feature,
    lfc = coef,
    se = stderr,
    p = pval,
    q = qval,
    diff = qval < 0.05,
    direction = case_when(
      qval < 0.05 & coef > 0 ~ "Up_in_Tumour",
      qval < 0.05 & coef < 0 ~ "Down_in_Tumour",
      TRUE ~ "NS"
    )
  )

# Per-taxon mean relative abundance (paired samples) for the > 0.01% figure filter
maaslin_tab$mean_rel <- rowMeans(
  sweep(mtx_count[, paired_samples, drop = FALSE], 2,
        colSums(mtx_count[, paired_samples, drop = FALSE]), "/")
)[maaslin_tab$taxon]

# Circos differential annotation from MaAsLin2 (match species by punctuation-normalised name)
maaslin_key <- gsub("[[:punct:]]", "_", maaslin_tab$taxon)
mi <- match(abund_sp_top$Species, maaslin_key)
abund_sp_top$Diff.Trend  <- maaslin_tab$direction[mi]
abund_sp_top$Diff.Padj   <- maaslin_tab$q[mi]
abund_sp_top$Diff.Log2FC <- maaslin_tab$lfc[mi]
abund_sp_top$Diff.Trend[is.na(abund_sp_top$Diff.Trend)]   <- "NS"
abund_sp_top$Diff.Log2FC[is.na(abund_sp_top$Diff.Log2FC)] <- 0

#saveRDS(abund_sp_top, file.path(DIR_RDS, "sAEG_CirclizeData_AbundSpTop_Genera20.rds"))
#abund_sp_top <- readRDS(file.path(DIR_RDS, "sAEG_CirclizeData_AbundSpTop_Genera20.rds"))

col_genera <- paste0(substr(paletteer_d("khroma::discreterainbow")[c(10, 12:20, 
                                                                     23:27, 2, 4, 5, 7, 9)], 
                            1, 7), "80") # Colours for genera. "80" reprensents alpha
col_sp <- substr(paletteer_d("khroma::discreterainbow")[c(10, 12:20, 
                                                          23:27, 2, 4, 5, 7, 9)], 
                 1, 7)

# Map a numeric metric onto a divergent colour gradient.
# Values above `hi` use the upper palette, below `lo` the lower one,
# and the in-between values share a neutral grey. Stronger absolute
# values (further from the neutral band) receive more saturated colours.
gradient_colours <- function(values, hi, lo, col_hi, col_lo, col_mid = "#EDEDED") {
  n_hi <- sum(values > hi)
  n_lo <- sum(values < lo)
  pal_hi <- colorRampPalette(c(col_hi, "#CDCDCD"))(n_hi)
  pal_lo <- colorRampPalette(c(col_lo, "#CDCDCD"))(n_lo)
  ord <- order(abs(values), decreasing = TRUE)
  out <- character(length(values))
  i_hi <- 0L; i_lo <- 0L
  for (k in seq_along(ord)) {
    v <- values[ord[k]]
    if (v > hi)       { i_hi <- i_hi + 1L; out[k] <- pal_hi[i_hi] }
    else if (v < lo)  { i_lo <- i_lo + 1L; out[k] <- pal_lo[i_lo] }
    else              { out[k] <- col_mid }
  }
  names(out) <- names(values)[ord]
  out[names(values)]
}

# Survival hazard ratio colours; keyed by rownames(abund_sp_top)
colour_hr <- gradient_colours(
  setNames(abund_sp_top$Surv.HR, rownames(abund_sp_top)),
  hi = 1.1, lo = 0.9, col_hi = "#CC6329", col_lo = "#409161"
)

# Differential log2 fold change colours; re-keyed by Species for the loop below
colour_log2fc <- gradient_colours(
  setNames(abund_sp_top$Diff.Log2FC, abund_sp_top$Species),
  hi = log2(1.5), lo = -log2(1.5), col_hi = "#9A342C", col_lo = "#2C6199"
)

# Drawing the most abundant genera section
# Allocating the canvas and sectors
pdf(file.path(DIR_RES, "Circlised_genus_species_filtered.pdf"), width = 8, height = 8)
par(mar = c(1, 1, 1, 1) * 11, cex = 0.6, xpd = NA)
sectors <- factor(names(g_abund), levels = names(g_abund))
circos.par(points.overflow.warning = FALSE,
           cell.padding = c(0, 0, 0, 0))
circos.initialize(factors = sectors, xlim = c(0, 1), sector.width = len_each_g)
circos.trackPlotRegion(factors = sectors, ylim = c(0, 12), 
                       track.height = 0.5, bg.border = "grey",
                       bg.col = col_genera)

# Drawing the species through loops
for (i in seq_along(g_abund)) {
  # Get species indices for current genus and sort by abundance
  sp_indices <- which(abund_sp_top$Genus == sectors[i])
  sp_order <- sp_indices[order(abund_sp_top$Abundance[sp_indices], decreasing = TRUE)]
  
  # Calculate uniform width for each species within this genus
  bar_width <- 1 / len_each_g[i] / 2  # Half width of each position
  
  for (j in seq_len(len_each_g[i])) {
    sp_idx <- sp_order[j]
    
    # Calculate the center position for this species
    x_center <- (j - 0.5) / len_each_g[i]
    
    # The length of each species-bar represents its relative abundance (log-transformed for better visualization)
    circos.rect(xleft = x_center - bar_width * 0.8,
                ybottom = 0,
                xright = x_center + bar_width * 0.8,
                ytop = log10(abund_sp_top$Abundance[sp_idx] + 1) * 
                  (12 / log10(max(abund_sp_top$Abundance) + 1)),
                sector.index = sectors[i],
                col = col_sp[i],
                border = NA)
    # The annotation circle representing survival risks of species
    circos.trackLines(sectors = sectors[i],
                      x = x_center,
                      y = -1,
                      col = colour_hr[sp_idx],
                      type = "h",
                      baseline = -2)
    # The annotation sub-circle representing survival significance (Surv.P)
    circos.trackLines(sectors = sectors[i],
                      x = x_center,
                      y = if_else(
                        abund_sp_top$Surv.Risk[sp_idx] != "NS", -2.4, -2.3
                      ),
                      col = if_else(
                        abund_sp_top$Surv.Risk[sp_idx] != "NS", "black", "white"
                      ),
                      type = "h",
                      baseline = -2.3)
    # The annotation circle representing differential log2foldchange between paired tumour & normal samples
    circos.trackLines(sectors = sectors[i],
                      x = x_center,
                      y = -4,
                      col = colour_log2fc[sp_idx],
                      type = "h",
                      baseline = -3)
    # The annotation sub-circle representing differential significance (Diff.Padj)
    circos.trackLines(sectors = sectors[i],
                      x = x_center,
                      y = if_else(
                        abund_sp_top$Diff.Trend[sp_idx] != "NS", -4.4, -4.3
                      ),
                      col = if_else(
                        abund_sp_top$Diff.Trend[sp_idx] != "NS", "black", "white"
                      ),
                      type = "h",
                      baseline = -4.3)
    
  }
  # The names of sectors (genera)
  circos.trackText(sectors = sectors[i],
                   x = 0.5, y = 7, cex = fontsize(14),
                   labels = names(g_abund)[i],
                   facing = "downward",
                   col = "black")
}

# Adding the exact scale of relative abundance
circos.yaxis(sector.index = sectors[11], side = "left", col = "grey30",
             tick.length = convert_x(.4, "mm", sectors[10]))

circos.clear()
dev.off()

# Drawing the legend of differential log2FC using ggplot2
# The only thing needed is the legend, so the main body of the plot can be ignored
pdf(file.path(DIR_RES, "ScaleSurv_Legend.pdf"), width = 6, height = 5)
ggplot(data = abund_sp_top, aes(colour = Surv.HR, 
                                x = Species, y = X)) +
  geom_point() + 
  scale_color_gradient2(aes(labs = c(min(Surv.HR), 1, max(Surv.HR))),
                        midpoint = 1,  
                        low = "#409161", 
                        mid = "#CDCDCD", 
                        high = "#CC6329")
dev.off()

pdf(file.path(DIR_RES, "ScaleDiff_Legend.pdf"), width = 6, height = 5)
ggplot(data = abund_sp_top, aes(colour = Diff.Log2FC, 
                                x = Species, y = X)) +
  geom_point() + 
  scale_color_gradient2(aes(labs = c(min(Diff.Log2FC), 0, max(Diff.Log2FC))),
                        midpoint = 0,  
                        low = "#2C6199", 
                        mid = "#CDCDCD", 
                        high = "#9A342C")
dev.off()

#######################################
# Alpha-diversity & ecological distance
shan_tumour <- apply(mtx_count[, grepl("C", colnames(mtx_count))],
                     MARGIN = 2, FUN = vegan::diversity)
dist_tumour <- clinical$`Distance from the tumor center to the esophagogastric junction()`
names(dist_tumour) <- paste0("C", clinical$No.)

overlap_samples <- intersect(names(shan_tumour), names(dist_tumour))
shan_tumour <- shan_tumour[overlap_samples]
dist_tumour <- dist_tumour[overlap_samples]

cor.test(shan_tumour, dist_tumour) # NS

summary(dist_tumour)

pdf(file.path(DIR_RES, "Test_alpha_distance.pdf"), width = 4, height = 3.2)
ggboxplot(data.frame(
  sample = overlap_samples,
  distance = dist_tumour,
  shannon = shan_tumour
) %>%
  mutate(dist_group = if_else(distance > median(distance), "Close", "Far")),
x = "dist_group", y = "shannon",
col = "dist_group", palette = "jco", add = "jitter") +
  stat_compare_means() + 
  xlab("Distance from tumour centre to junction") +
  ylab("Shannon index") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1)) # NS
dev.off()

# Alpha-diversity & siewart type
siewart_tumour <- clinical$`Siewert type`
names(siewart_tumour) <- paste0("C", clinical$No.)
siewart_tumour <- siewart_tumour[overlap_samples]

pdf(file.path(DIR_RES, "Test_alpha_siewart.pdf"), width = 4, height = 3.2)
ggboxplot(data.frame(
  sample = overlap_samples,
  siewart = siewart_tumour,
  shannon = shan_tumour
) %>%
  mutate(siewart = factor(siewart, levels = c("I", "II", "III"))),
x = "siewart", y = "shannon",
col = "siewart", palette = "jco", add = "jitter") +
  stat_compare_means(comparisons = list(c("I", "II"),
                                        c("II", "III"),
                                        c("I", "III"))) + 
  xlab("Siewart type") +
  ylab("Shannon index") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1)) # NS
dev.off()

#############################
# Alpha diversity & prognosis
# Note: hr_calc() was sourced earlier in this script
df_surv <- data.frame(
  sample = paste0("C", clinical$No.),
  time = clinical$`Survival time（Month）`,
  state = clinical$`Status（1=Dead, 0=Alive）`
) %>%
  filter(sample %in% overlap_samples) %>%
  arrange(., sample)
shan_tumour <- shan_tumour[df_surv$sample]

df_surv <- df_surv %>%
  mutate(shannon = shan_tumour,
         group = if_else(shan_tumour > median(shan_tumour), "High", "Low"))

sfit <- survfit(Surv(time, state) ~ df_surv$group, data = df_surv)
surv_res <- ggsurvplot(sfit, conf.int = T, pval = F, risk.table = F, 
                       legend.labs = c("High", "Low"), 
                       legend.title = "Diversity",  
                       surv.median.line = "hv",
                       palette = c("tomato3", "steelblue3"), 
                       title = "Alpha-diversity",
                       ggtheme = theme_test() +
                         theme(panel.border = element_rect(fill = NA, colour = 1),
                               axis.text = element_text(colour = 1)))
# Adding HR and its p-value to KM curve
df_hr <- merge(df_surv, clinical %>% mutate(No. = paste0("C", No.)),
               by.x = "sample", by.y = "No.", all.x = T, all.y = F) 
colnames(df_hr) <- gsub(" ", "_", colnames(df_hr))

df_hr <- df_hr %>%
  mutate(Pathological_stage = case_when(
    Pathological_stage == "IB" ~ 1.5,
    Pathological_stage == "IIA" ~ 2,
    Pathological_stage == "IIB" ~ 2.5,
    Pathological_stage == "IIIA" ~ 3,
    Pathological_stage == "IIIB" ~ 3.3,
    Pathological_stage == "IIIC" ~ 3.6,
    Pathological_stage %in% c("IV", "IV") ~ 4
  ),
  Differentiated_degree = case_when(
    Differentiated_degree == "Moderately differentiated" ~ 1,
    Differentiated_degree == "Low differentiated" ~ 2,
    Differentiated_degree == "Poorly differentiated" ~ 3
  ),
  Smoking = if_else(Smoking == "No", 0, 1),
  Alcohol = if_else(Alcohol == "No", 0, 1)
  )

hr_clinic <- hr_calc(c("shannon", "Age", "Smoking", "Alcohol",
                       "Pathological_stage", "Differentiated_degree"), 
                     df_hr, 1.1, 0.9, 0.05) %>%
  arrange(desc(HR)) %>%
  mutate(index = factor(rownames(.), levels = rownames(.)))

pdf(file.path(DIR_RES, "B_hr_alphadiv_clinic.pdf"), width = 4.5, height = 4)
ggplot(hr_clinic, aes(x = HR, y = index,
                      color = Risk)) +
  geom_point(shape = 15, size = 5) +
  scale_color_manual(values = c("Decrease" = "#126CAA", 
                                "NS" = "grey", 
                                "Increase" = "#9A342C")) +
  geom_errorbar(aes(xmin = HR.95L, xmax = HR.95H),
                width = .2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  xlab("Hazard Ratio (95% CI)") +
  ylab("Index") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        panel.grid = element_blank(),
        axis.text = element_text(colour = 1))
dev.off()

hr_shan <- hr_calc("shannon", df_hr, 1.1, 0.9, 0.05)

pdf(file.path(DIR_RES, "Alpha_km.pdf"), width = 3.8, height = 4)
surv_res$plot +
  ggplot2::annotate(
    "text",
    x = 140, y = 0.95,
    vjust = 1, hjust = 1,
    label = sprintf("HR = %.3f\np = %.3f", hr_shan$HR, hr_shan$Pvalue),
    size = 5
  )
dev.off()

##################
# Network analysis
ps.obj <- import_biom("./aeg.biom")
sample_data(ps.obj)$ID <- sample_data(ps.obj)$Id
sample_data(ps.obj)$group <- ifelse(grepl("C", sample_data(ps.obj)$Id), 
                                    "Tumour", 
                                    "Normal")
tab <- network.pip(ps = ps.obj, N = 200, # ra = 0.05,
                   big = FALSE, select_layout = FALSE,
                   layout_net = "model_maptree2",
                   r.threshold = 0.6, p.threshold = 0.05,
                   maxnode = 2, method = "kendall",
                   label = FALSE, lab = "elements",
                   group = "group", fill = "Species",
                   size = "igraph.degree", zipi = TRUE,
                   ram.net = TRUE, clu_method = "cluster_fast_greedy",
                   step = 100, R = 10, ncpus = 6
)
saveRDS(tab, file.path(DIR_RDS, "AEG_network.rds"))

cortab <- tab[[2]]$net.cor.matrix$cortab
saveRDS(cortab, file.path(DIR_RDS, "AEG_network_cor.rds"))

# Network plot
# Filter out viruses
net_dt <- tab[[2]]
node <- net_dt$net.cor.matrix$node %>%
  filter(Rank1 != "k__Viruses",
         igraph.degree != 0)
edge <- net_dt$net.cor.matrix$edge %>%
  filter(OTU_1 %in% node$ID & OTU_2 %in% node$ID)

# Replace blank ranks with "unclassified"
for (col in paste0("Rank", 1:7)) {
  node[[col]] <- ifelse(grepl("^[a-z]__$", node[[col]]), 
                        paste0(sub("__$", "", node[[col]]), "__unclassified"), 
                        node[[col]])
}

fmt_p <- function(p) ifelse(is.na(p), "= NA",
                            ifelse(p < 0.001, "< 0.001",
                                   paste0("= ", formatC(p, format = "f", digits = 3))))

# Zi-Pi (within-module connectivity vs participation coefficient), Normal vs Tumour.
# Significance: Wilcoxon rank-sum comparing the two networks' node distributions.
zp_dat <- as.data.frame(tab[[1]][[2]]$layers[[2]]$data)   # cols incl. z (Zi), p (Pi), group
zi_p <- wilcox.test(z ~ group, data = zp_dat)$p.value
pi_p <- wilcox.test(p ~ group, data = zp_dat)$p.value
pdf(file.path(DIR_RES, "Net_connectivity_tn.pdf"), width = 7.5, height = 5)
print(
  tab[[1]][[2]] +
    labs(x = "Participation coefficient (Pi, 0-1)",
         y = "Within-module connectivity z-score (Zi)",
         subtitle = paste0("Normal vs Tumour (Wilcoxon):  Pi p ", fmt_p(pi_p),
                           ",   Zi p ", fmt_p(zi_p))) +
    theme(plot.subtitle = element_text(size = 9))
)
dev.off()

# Degree distribution: observed network vs Erdos-Renyi (E-R) random, per network.
# x = node degree k, y = proportion of nodes P(k) (sums to 1 within each network).
# Significance: KS test of observed vs E-R degree distribution, within each facet.
er_dat <- as.data.frame(tab[[1]][[3]]$data)   # cols: network = P(k), group, ID = degree k, g = facet
recon_deg <- function(P, k) {                  # rebuild node-degree sample from P(k)
  P <- pmax(P, 0); pos <- P > 0
  if (!any(pos)) return(numeric(0))
  n <- round(1 / min(P[pos]))
  rep(k, round(P * n))
}
ks_df <- do.call(rbind, lapply(split(er_dat, er_dat$g), function(d) {
  obs <- d[d$group == "network", ]
  ern <- d[d$group != "network", ]
  pv <- tryCatch(suppressWarnings(
    ks.test(recon_deg(obs$network, obs$ID), recon_deg(ern$network, ern$ID))$p.value),
    error = function(e) NA_real_)
  data.frame(g = d$g[1], lab = paste0("KS (obs vs E-R) p ", fmt_p(pv)))
}))
pdf(file.path(DIR_RES, "Net_randomness_tn.pdf"), width = 7.5, height = 5)
print(
  tab[[1]][[3]] +
    scale_colour_manual(values = c("#0073C2FF", "#EFC000FF")) +
    scale_fill_manual(values = c("#0073C2FF", "#EFC000FF")) +
    labs(x = "Node degree (k)", y = "Proportion of nodes, P(k)") +
    geom_text(data = ks_df, aes(x = Inf, y = Inf, label = lab),
              inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 3) +
    theme(panel.border = element_rect(fill = NA, colour = 1),
          axis.text = element_text(colour = 1))
)
dev.off()

# Network destruction resistance
resis <- natural.con.microp(ps = ps.obj, corg = cortab,
                            norm = TRUE, end = 150, start = 0)
saveRDS(resis, file.path(DIR_RDS, "network_destruction_resistance.rds"))

pdf(file.path(DIR_RES, "Net_resistance.pdf"), width = 5, height = 4)
resis[[1]] +
  scale_colour_manual(values = c(Normal = "#126CAA", 
                                 Tumour = "#9A342C")) +
  labs(x = "Number of removed nodes",
       y = "Natural connectivity") +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))
dev.off()
write.csv(resis[[2]], file.path(DIR_TAB, "./Res2_network_resistance.csv"))

# Relative change
resis_data <- resis[[2]] %>%
  group_by(Group) %>%
  mutate(baseline = Natural.connectivity[1],
         changed_connectivity = Natural.connectivity / baseline) %>%
  ungroup()

pdf(file.path(DIR_RES, "Net_resistance_changed_0to100.pdf"), width = 5, height = 4)
ggplot(resis_data, aes(x = `Num.of.remove.nodes`, 
                       y = changed_connectivity, colour = Group)) +
  geom_point(alpha = 0.3) +
  xlim(0, 100) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  scale_colour_manual(values = c(Normal = "#126CAA", 
                                 Tumour = "#9A342C")) +
  labs(x = "Number of removed nodes",
       y = "Connectivity - relative change") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))
dev.off()

####################
# Network similarity
module <- module.compare.m(ps = ps.obj, corg = cortab, zipi = FALSE,
                           zoom = 0.2, padj = F, n = 3)
saveRDS(module, file.path(DIR_RDS, "AEG_network_module.rds"))

pdf(file.path(DIR_RES, "Net_module.pdf"), width = 4, height = 4)
module[[1]]
dev.off()

module_otu <- module[[2]]
module_otu$taxa_g <- ps.obj@tax_table[module_otu$ID, "Rank6"]
module_otu$taxa_s <- ps.obj@tax_table[module_otu$ID, "Rank7"]
module_otu$taxa_k <- ps.obj@tax_table[module_otu$ID, "Rank1"]
module_otu <- module_otu %>%
  filter(taxa_k != "k__Viruses",
         taxa_g != "g__",
         taxa_s != "s__") %>%
  mutate(standard_name = paste0(gsub("g__", "", taxa_g)  %>% 
                                  toTitleCase(.), "_", 
                                gsub("s__", "", taxa_s)),)
head(module_otu)
table(module_otu$group)

# Re-visualise the module relationships
# Expand nodes: each (ID, group) pair becomes a unique node
species_modules_expanded <- module_otu %>%
  select(ID, group) %>%
  distinct() %>%
  mutate(node_id = paste0(ID, "__", group))

# IDs that appear in both Normal and Tumour modules (conserved)
shared_ids <- species_modules_expanded %>%
  group_by(ID) %>%
  summarise(
    has_normal = any(grepl("Normal", group)),
    has_tumour = any(grepl("Tumour", group)),
    .groups = "drop"
  ) %>%
  filter(has_normal & has_tumour) %>%
  pull(ID)

# Node metadata with final node_type classification
nodes_expanded <- species_modules_expanded %>%
  left_join(
    module_otu %>% select(ID, standard_name) %>% distinct(),
    by = "ID"
  ) %>%
  mutate(node_type = case_when(
    ID %in% shared_ids      ~ "Conserved",
    grepl("Normal", group)  ~ "Normal-enriched",
    grepl("Tumour", group)  ~ "Tumour-enriched"
  ))

# Edges: within-module, using expanded node_ids
edges_within <- species_modules_expanded %>%
  inner_join(species_modules_expanded, by = "group", relationship = "many-to-many") %>%
  filter(node_id.x < node_id.y) %>%
  select(from = node_id.x, to = node_id.y)

# Edges: cross-module for shared nodes
edges_shared <- species_modules_expanded %>%
  filter(ID %in% shared_ids) %>%
  inner_join(
    species_modules_expanded %>% filter(ID %in% shared_ids),
    by = "ID", relationship = "many-to-many"
  ) %>%
  filter(group.x < group.y) %>%
  select(from = node_id.x, to = node_id.y) %>%
  distinct()

edges_all <- bind_rows(edges_within, edges_shared) %>%
  distinct() %>%
  mutate(
    edge_type = if_else(
      paste0(from, to) %in% paste0(edges_shared$from, edges_shared$to),
      "cross-model", "within-model"
    )
  )

g <- graph_from_data_frame(
  edges_all,
  directed = FALSE,
  vertices = nodes_expanded %>% select(node_id, standard_name, group, node_type) %>%
    rename(name = node_id)
)

node_colours <- c(
  "Normal-enriched" = "#4393C3",
  "Tumour-enriched" = "#DF8F44",
  "Conserved"       = "#79AF97"
)
layout <- create_layout(g, layout = 'kk')
pdf(file.path(DIR_RES, "C_net_module_revisualised.pdf"), width = 5, height = 4)
ggraph(layout) +
  geom_edge_link(aes(linetype = edge_type), alpha = 0.5, colour = "grey", width = 0.3) +
  scale_edge_linetype_manual(values = c("within-model" = "solid", "cross-model" = "dashed")) +
  geom_node_point(aes(colour = node_type), size = 1.5) +
  scale_colour_manual(values = node_colours, name = "Node type") +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank())
dev.off()

# Aggregate standard names by group in module_otu
module_species <- module_otu %>%
  group_by(group) %>%
  summarise(species = paste(standard_name, collapse = "|"),
            otu = paste(ID, collapse = "|"))

# Large Normal modules (>5 named-species nodes after the g__/s__ filter)
normal_top <- module_otu %>%
  filter(Group == "Normal") %>%
  count(group) %>%
  filter(n > 5) %>%
  pull(group)
tumour_module_map <- module_otu %>%
  filter(Group == "Tumour") %>%
  distinct(OTU = ID, tumour_module = group)
saving_module <- module_otu %>%
  filter(group %in% normal_top) %>%
  distinct(OTU = ID, Species = standard_name, normal_module = group) %>%
  left_join(tumour_module_map, by = "OTU") %>%
  mutate(tumour_module = replace_na(tumour_module, "none")) %>%
  merge(., node, by.x = "OTU", by.y = "ID")
saving_module$abundance <-
  (apply(mtx_cpm, MARGIN = 1, FUN = mean) / 1e+4)[saving_module$Species]
saving_module <- saving_module[!is.na(saving_module$abundance), ]

# Regression: Shannon index - candidate species
elnet_sp <- unique(saving_module$Species)
tumour_module_sp <- strsplit(module_species$species[module_species$group %in% 
                                                      paste0("Tumourmodel_", 1:5)], "\\|") %>%
  unlist(.) %>%
  unique(.)
elnet_sp <- elnet_sp[elnet_sp %in% tumour_module_sp]
x_raw <- as.matrix(log2(mtx_cpm[elnet_sp, names(shan_tumour)] + 1))
x_clr <- apply(x_raw, 2, function(v) clr(v)) %>% t(.)

# Ridge regression
set.seed(42)  # reproducible CV folds, independent of upstream RNG
elnet_fit <- cv.glmnet(
  x = x_clr,
  y = shan_tumour,
  alpha = 0,
  family = "gaussian",
  standardize = TRUE
)

coef_elnet <- coef(elnet_fit, s = "lambda.min") %>%
  as.matrix(.) %>%
  as.data.frame(.) %>%
  mutate(species = rownames(.)) %>%
  filter(species != "(Intercept)")

positive_vals <- coef_elnet$lambda.min[coef_elnet$lambda.min > 0]
coef_elnet <- coef_elnet %>%
  mutate(rank_contribute = ifelse(
    lambda.min > 0,
    rank(-positive_vals)[match(lambda.min, positive_vals)],
    NA_real_
  ))

# Lollipop: species contribution to Shannon diversity, coloured by normal module
contrib_df <- coef_elnet %>%
  left_join(saving_module %>% distinct(Species, normal_module),
            by = c("species" = "Species")) %>%
  mutate(coef = lambda.min)
mod_cols <- c("Normalmodel_1" = "#4393C3", "Normalmodel_4" = "#5AAE61")

lollipop_contribution <- function(df, fname, height, width = 5.5) {
  df <- df %>% mutate(species = factor(species, levels = species[order(coef)]))
  p <- ggplot(df, aes(x = coef, y = species, colour = normal_module)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_segment(aes(x = 0, xend = coef, yend = species), linewidth = 0.6) +
    scale_colour_manual(values = mod_cols, name = "Normal module") +
    labs(x = "Ridge coefficient (contribution to Shannon diversity)", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")
  if (any(df[["coef"]] == 0)) {
    p <- p + geom_point(aes(shape = coef == 0), size = 2.6) +
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1),
                         labels = c("selected", "shrunk to 0"), name = NULL)
  } else {
    p <- p + geom_point(size = 2.6)
  }
  ggsave(file.path(DIR_RES, fname), p, width = width, height = height)
}

# Main figure: positively contributing species
lollipop_contribution(contrib_df %>% filter(coef > 0),
                      "Elnet_shannon_contribution_positive.pdf", 3.2)
# Supplementary figure: zero- and negative-contribution species
lollipop_contribution(contrib_df %>% filter(coef <= 0),
                      "Elnet_shannon_contribution_zero_negative.pdf", 4.5, width = 8)

# Rank all parameters
library(patchwork)
# Multi-attribute ranking of positively-contributing species
# (degree/closeness/betweenness incl. normal-tumour stability; betweenness reversed: low = best)
cand_sp <- coef_elnet[["species"]][!is.na(coef_elnet[["rank_contribute"]])]
saving_normal <- saving_module %>% filter(Group == "Normal" & Species %in% cand_sp) %>%
  mutate(degree_rank = rank(-igraph.degree, ties.method = "min"),
         closeness_rank = rank(-igraph.closeness, ties.method = "min"),
         betweenness_rank = rank(igraph.betweenness, ties.method = "min"))
saving_tumour <- saving_module %>% filter(Group == "Tumour" & Species %in% cand_sp) %>%
  mutate(degree_rank = rank(-igraph.degree, ties.method = "min"),
         closeness_rank = rank(-igraph.closeness, ties.method = "min"),
         betweenness_rank = rank(igraph.betweenness, ties.method = "min"))
saving_info <- merge(
  saving_normal[, c("OTU", "Species", "abundance", "degree_rank", "closeness_rank", "betweenness_rank")],
  saving_tumour[, c("OTU", "Species", "degree_rank", "closeness_rank", "betweenness_rank")],
  by = c("OTU", "Species"), suffixes = c(".normal", ".tumour")) %>%
  mutate(degree_stability = rank(abs(degree_rank.normal - degree_rank.tumour)),
         degree_rank = (rank(degree_rank.normal) + rank(degree_rank.tumour)) / 2,
         closeness_stability = rank(abs(closeness_rank.normal - closeness_rank.tumour)),
         closeness_rank = (rank(closeness_rank.normal) + rank(closeness_rank.tumour)) / 2,
         betweenness_stability = rank(abs(betweenness_rank.normal - betweenness_rank.tumour)),
         betweenness_rank = (rank(betweenness_rank.normal) + rank(betweenness_rank.tumour)) / 2,
         diversity_contribute = coef_elnet[Species, "rank_contribute"],
         rank_abundance = rank(-abundance))

# Composite selection figure: lollipop | abundance | network bubbles
netvars <- c("degree_stability", "degree_rank", "closeness_stability",
             "closeness_rank", "betweenness_stability", "betweenness_rank")
si2 <- saving_info %>% mutate(across(all_of(c(netvars, "rank_abundance")), ~ rank(.x), .names = "{.col}_g"))
sp_levels <- saving_info %>% arrange(desc(diversity_contribute)) %>% pull(Species)
short_map <- setNames(paste0(toupper(substr(sp_levels, 1, 1)), ".", sub("^[^_]*_", "", sp_levels)), sp_levels)
mod_cols <- c("Normalmodel_1" = "#4393C3", "Normalmodel_4" = "#5AAE61")
tl <- "#0C8599"; th <- "#79C9C4"
loll <- saving_info %>% mutate(coef = coef_elnet[Species, "lambda.min"]) %>%
  left_join(saving_module %>% distinct(Species, normal_module), by = "Species") %>%
  mutate(Species = factor(Species, levels = sp_levels))
netg <- paste0(netvars, "_g")
bubN <- si2 %>% select(Species, all_of(netg)) %>%
  pivot_longer(all_of(netg), names_to = "attr", values_to = "rank") %>%
  mutate(attr = factor(attr, levels = netg), sz = max(rank) + 1 - rank, Species = factor(Species, levels = sp_levels))
bubA <- si2 %>% transmute(Species = factor(Species, levels = sp_levels), rank = rank_abundance_g,
                          sz = max(rank_abundance_g) + 1 - rank_abundance_g, attr = factor("Abundance"))
netlabs <- c("Degree stab.", "Degree", "Closeness stab.", "Closeness", "Betweenness stab.", "Betweenness")
pL <- ggplot(loll, aes(x = coef, y = Species, colour = normal_module)) +
  geom_segment(aes(x = 0, xend = coef, yend = Species), linewidth = 0.6) + geom_point(size = 2.4) +
  scale_colour_manual(values = mod_cols, name = "Normal module") + scale_x_reverse() + scale_y_discrete(labels = short_map) +
  labs(x = "Ridge coef", y = NULL, title = "Diversity contribution") +
  theme_bw(base_size = 8) + theme(panel.grid.major.y = element_blank(), legend.position = "bottom", plot.title = element_text(size = 8, hjust = 0.5))
pN <- ggplot(bubN, aes(x = attr, y = Species)) + geom_point(aes(size = sz, colour = rank)) +
  scale_colour_gradient(low = tl, high = th, name = "Rank") + scale_size(range = c(1.5, 5), guide = "none") +
  scale_x_discrete(labels = netlabs) + labs(x = NULL, y = NULL, title = "Network rank (stability + rank)") +
  theme_bw(base_size = 8) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.title = element_text(size = 8, hjust = 0.5))
pA <- ggplot(bubA, aes(x = attr, y = Species)) + geom_point(aes(size = sz, colour = rank)) +
  scale_colour_gradient(low = tl, high = th, guide = "none") + scale_size(range = c(1.5, 5), guide = "none") +
  labs(x = NULL, y = NULL, title = "Abundance") +
  theme_bw(base_size = 8) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 8, hjust = 0.5))
combined <- pL + pA + pN + plot_layout(widths = c(1.3, 0.45, 1.9))
ggsave(file.path(DIR_RES, "D_core_sp_selection_composite.pdf"), combined, width = 6.5, height = 4.6)

write.csv(saving_info[, c("Species", "degree_stability", "degree_rank", "closeness_stability",
                          "closeness_rank", "betweenness_stability", "betweenness_rank",
                          "rank_abundance", "diversity_contribute")],
          file.path(DIR_TAB, "Species_within_saving_module.csv"), row.names = FALSE)

# Candidate species (positively contributing) for the network plot
species_list <- as.character(saving_info[["Species"]])

pdf(file.path(DIR_RES, "Net_all_candidates.pdf"), width = 12, height = 5)
ggplot() +
  geom_segment(data = edge, alpha = 0.3,
               aes(x = X1, y = Y1, xend = X2, yend = Y2,
                   linewidth = weight)) +
  scale_linewidth_continuous(range = c(0.01, 0.05)) +
  geom_point(data = node, pch = 21, color = "gray40",
             aes(X1, X2, fill = Rank2, size = igraph.degree)) +
  scale_fill_manual(values = paletteer_d("ggsci::nrc_npg")) +
  facet_wrap(.~ label, scales = "free_y", nrow = 1) +
  ggrepel::geom_text_repel(
    data = node %>%
      mutate(standard_name = paste0(gsub("g__", "", Rank6), "_",
                                    gsub("s__", "", Rank7)),
             short_name = paste0(toupper(substr(gsub("g__", "", Rank6), 1, 1)),
                                 ".", gsub("s__", "", Rank7))) %>%
      filter(standard_name %in% species_list),
    aes(X1, X2, label = short_name),
    size = 2.6, fontface = "italic",
    max.overlaps = Inf, force = 8, force_pull = 0.4,
    box.padding = 0.8, point.padding = 0.3,
    min.segment.length = 0, segment.size = 0.25,
    segment.color = "grey45", segment.alpha = 0.85,
    seed = 42
  ) +
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL, expand = expansion(mult = 0.15)) +
  scale_y_continuous(breaks = NULL, expand = expansion(mult = 0.12)) +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_rect(colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
dev.off()

# Build correlation stats using samples shared by Shannon and species matrices
shared_samples <- intersect(names(shan_tumour), colnames(mtx_cpm))
cor_stats <- list()

for (sp in species_list) {
  abund <- mtx_cpm[sp, shared_samples] / 1e+4
  sub   <- data.frame(
    abund   = as.numeric(abund),
    shannon = shan_tumour[shared_samples]
  )
  if (sd(sub$abund, na.rm = TRUE) == 0) {
    cor_stats[[sp]] <- data.frame(species = sp, r = NA_real_, p = NA_real_)
  } else {
    ct <- cor.test(sub$abund, sub$shannon, method = "pearson")
    cor_stats[[sp]] <- data.frame(
      species = sp,
      r       = round(ct$estimate, 3),
      p       = ct$p.value
    )
  }
}

cor_stats_df <- do.call(rbind, cor_stats)
cor_stats_df$p.adj <- p.adjust(cor_stats_df$p, method = "BH")

# Plot one figure per species
for (sp in species_list) {
  sub <- data.frame(
    abund = as.numeric(mtx_cpm[sp, shared_samples] / 1e+4),
    shannon = shan_tumour[shared_samples]
  )
  stat <- cor_stats_df[cor_stats_df$species == sp, ]
  
  if (!is.na(stat$r)) {
    annot_label <- paste0("R = ", stat$r,
                          "\np = ", format.pval(stat$p.adj, 3))
    annot_x <- max(sub$abund)
    annot_y <- max(sub$shannon)
  }
  
  p <- ggplot(sub, aes(x = abund, y = shannon)) +
    geom_point(color = "tomato2", size = 1.5) +
    geom_smooth(method = lm, color = "#00A087",
                linewidth = 1.5, fill = "#4DBBD5") +
    { if (!is.na(stat$r))
      annotate("text", x = annot_x, y = annot_y,
               color = 1, hjust = 1, vjust = 1, size = 3,
               label = annot_label)
    } +
    xlab(paste0(gsub("_", " ", sp), " abundance (%)")) +
    ylab("Shannon index") +
    theme_test() +
    theme(panel.border = element_rect(fill = NA, colour = 1),
          axis.text    = element_text(colour = 1))
  pdf(file.path(DIR_RES, paste0("Cor_shannon_", sp, ".pdf")), width = 3.3, height = 3.3)
  print(p)
  dev.off()
}

################
# Build sample metadata for all C/N samples in mtx_cpm
meta_all <- data.frame(
  sample = colnames(mtx_cpm),
  patient = gsub("^[CN]", "", colnames(mtx_cpm)),
  group = ifelse(grepl("^C", colnames(mtx_cpm)), "Tumour", "Normal"),
  stringsAsFactors = FALSE
)

# Merge clinical variables (Tumour-specific clinical info propagated via patient ID)
clin_lite <- clinical %>%
  transmute(
    patient = sprintf("%03d", as.integer(No.)),
    Siewert = `Siewert type`,
    Distance = `Distance from the tumor center to the esophagogastric junction()`,
    Stage = `Pathological stage`,
    Differentiation = `Differentiated degree`,
    Age = Age,
    Smoking = Smoking,
    Alcohol = Alcohol
  )
meta_all <- meta_all %>% left_join(clin_lite, by = "patient")

# Bray-Curtis on relative-abundance-transformed CPM
mtx_rel <- sweep(mtx_cpm, 2, colSums(mtx_cpm), "/")
bc_dist <- vegdist(t(mtx_rel), method = "bray")

# PERMANOVA: group effect (with patient as strata for paired structure)
permanova_group <- adonis2(
  bc_dist ~ group,
  data = meta_all,
  permutations = 999,
  strata = meta_all$patient,
  by = "terms"
)

# PERMDISP: dispersion homogeneity
permdisp_obj <- betadisper(bc_dist, meta_all$group)
permdisp_res <- permutest(permdisp_obj, permutations = 999)
saveRDS(permdisp_res, file.path(DIR_RDS, "PERMDISP_group.rds"))

# PERMANOVA on tumour-only samples for clinical drivers
meta_t <- meta_all %>% filter(group == "Tumour") %>%
  mutate(
    Stage_num = case_when(
      Stage == "IB" ~ 1.5, Stage == "IIA" ~ 2, Stage == "IIB" ~ 2.5,
      Stage == "IIIA" ~ 3, Stage == "IIIB" ~ 3.3, Stage == "IIIC" ~ 3.6,
      Stage == "IV" ~ 4, TRUE ~ NA_real_
    ),
    Diff_num = case_when(
      Differentiation == "Moderately differentiated" ~ 1,
      Differentiation == "Low differentiated" ~ 2,
      Differentiation == "Poorly differentiated" ~ 3,
      TRUE ~ NA_real_
    )
  )

bc_dist_t <- vegdist(t(mtx_rel[, meta_t$sample]), method = "bray")
meta_t_complete <- meta_t %>%
  filter(!is.na(Stage_num), !is.na(Diff_num),
         !is.na(Siewert), !is.na(Distance), !is.na(Age))
keep_idx <- match(meta_t_complete$sample, meta_t$sample)
bc_dist_t_sub <- as.dist(as.matrix(bc_dist_t)[keep_idx, keep_idx])

adonis2(
  bc_dist_t_sub ~ Age + Stage_num + Diff_num + Siewert + Distance,
  data = meta_t_complete,
  permutations = 999,
  by = "margin"
)

# PCoA visualisation with envfit overlay
pcoa <- cmdscale(bc_dist, k = 2, eig = TRUE)
pcoa_df <- data.frame(
  PCo1 = pcoa$points[, 1],
  PCo2 = pcoa$points[, 2],
  sample = rownames(pcoa$points)
) %>% left_join(meta_all, by = "sample")

eig_pct <- round(pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)

# Nestedness vs turnover decomposition (Baselga framework) within tumour samples
mtx_pa <- (mtx_rel > 0) * 1
print(unlist(beta.multi(
  t(mtx_pa[, meta_all$sample[meta_all$group == "Tumour"]]),
  index.family = "sorensen"
))) # Not nested

# Pairwise nestedness for paired patients
paired_ids <- intersect(
  gsub("^C", "", grep("^C", colnames(mtx_count), value = TRUE)),
  gsub("^N", "", grep("^N", colnames(mtx_count), value = TRUE))
)

cross_nest <- lapply(paired_ids, function(p) {
  pair_mat <- t(mtx_pa[, c(paste0("N", p), paste0("C", p))])
  bp <- beta.pair(pair_mat, index.family = "sorensen")
  data.frame(
    patient = p,
    sorensen = as.numeric(bp$beta.sor),
    turnover = as.numeric(bp$beta.sim),
    nestedness = as.numeric(bp$beta.sne),
    nest_ratio = as.numeric(bp$beta.sne) / as.numeric(bp$beta.sor)
  )
}) %>% do.call(rbind, .)

# Direction check: which side has fewer species (the nested one)?
cross_nest$tumour_nested_in_normal <- sapply(paired_ids, function(p) {
  sum(mtx_pa[, paste0("C", p)]) < sum(mtx_pa[, paste0("N", p)])
})

summary(cross_nest$nest_ratio)
table(cross_nest$tumour_nested_in_normal)

# Compute species richness per sample (number of detected taxa)
richness <- colSums((mtx_rel > 0) * 1)
pcoa_df$richness <- richness[pcoa_df$sample]

# Build paired-line dataframe for connecting T-N pairs
paired_ids_pcoa <- intersect(
  gsub("^C", "", pcoa_df$sample[pcoa_df$group == "Tumour"]),
  gsub("^N", "", pcoa_df$sample[pcoa_df$group == "Normal"])
)
pair_lines <- lapply(paired_ids_pcoa, function(p) {
  t_pt <- pcoa_df[pcoa_df$sample == paste0("C", p), c("PCo1", "PCo2")]
  n_pt <- pcoa_df[pcoa_df$sample == paste0("N", p), c("PCo1", "PCo2")]
  data.frame(
    patient = p,
    x = n_pt$PCo1, y = n_pt$PCo2,
    xend = t_pt$PCo1, yend = t_pt$PCo2
  )
}) %>% do.call(rbind, .)

# Annotate direction: does tumour have fewer species than its paired normal?
pair_lines$tumour_smaller <- sapply(paired_ids_pcoa, function(p) {
  richness[paste0("C", p)] < richness[paste0("N", p)]
})

# Main PCoA panel: richness-encoded points + paired connectors
p_pcoa_enriched <- ggplot() +
  geom_segment(
    data = pair_lines,
    aes(x = x, y = y, xend = xend, yend = yend, colour = tumour_smaller),
    alpha = 0.35, linewidth = 0.3
  ) +
  scale_colour_manual(
    values = c(`TRUE` = "#9A342C", `FALSE` = "grey60"),
    name = "Tumour richness < Normal",
    labels = c(`TRUE` = "Yes (nested-like)", `FALSE` = "No")
  ) +
  ggnewscale::new_scale_colour() +
  stat_ellipse(
    data = pcoa_df,
    aes(PCo1, PCo2, fill = group),
    geom = "polygon", alpha = 0.12, colour = NA
  ) +
  geom_point(
    data = pcoa_df,
    aes(PCo1, PCo2, colour = group, size = richness),
    alpha = 0.75
  ) +
  scale_colour_manual(values = c(Normal = "#126CAA", Tumour = "#9A342C")) +
  scale_fill_manual(values = c(Normal = "#126CAA", Tumour = "#9A342C")) +
  scale_size_continuous(range = c(0.8, 3.5), name = "Species richness") +
  labs(
    x = paste0("PCo1 (", eig_pct[1], "%)"),
    y = paste0("PCo2 (", eig_pct[2], "%)"),
    title = sprintf(
      "PERMANOVA R2 = %.3f, p = %.3f",
      permanova_group$R2[1], permanova_group$`Pr(>F)`[1]
    )
  ) +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))

# Marginal density of nestedness ratio
p_nest <- ggplot(cross_nest, aes(x = nest_ratio)) +
  geom_density(fill = "#9A342C", alpha = 0.4, colour = "#9A342C") +
  geom_vline(xintercept = median(cross_nest$nest_ratio),
             linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey40") +
  annotate("text",
           x = median(cross_nest$nest_ratio),
           y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("median = %.2f",
                           median(cross_nest$nest_ratio)),
           size = 3) +
  xlim(0, 1) +
  labs(x = "Nestedness / Total dissimilarity\n(per paired T-N)",
       y = "Density") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))

# Direction summary as small inset
n_smaller <- sum(pair_lines$tumour_smaller)
n_total <- nrow(pair_lines)
p_dir <- ggplot(
  data.frame(
    cat = c("T < N\n(nested-like)", "T >= N"),
    n = c(n_smaller, n_total - n_smaller)
  ),
  aes(x = cat, y = n, fill = cat)
) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("#9A342C", "grey70")) +
  labs(x = NULL, y = "Patient pairs") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        legend.position = "none")

final_plot_tn <- p_pcoa_enriched + (p_nest / p_dir) +
  plot_layout(widths = c(2.2, 1))

pdf(file.path(DIR_RES, "Pcoa_nestedness_check.pdf"), width = 8, height = 4.5)
print(final_plot_tn)
dev.off()

paired_samples <- c(paste0("C", paired_ids), paste0("N", paired_ids))

meta_paired <- data.frame(
  sample = paired_samples,
  patient = gsub("^[CN]", "", paired_samples),
  group = factor(
    ifelse(grepl("^C", paired_samples), "Tumour", "Normal"),
    levels = c("Normal", "Tumour")
  ),
  row.names = paired_samples
)

# MaAsLin2 (paired tumour-normal) is run once near the top; maaslin_tab is reused here.

n_each <- 10

# Combine: up first (largest at top), then down (most negative at bottom)
top_da <- bind_rows(
  maaslin_tab %>%
    filter(diff, direction == "Up_in_Tumour", mean_rel > rel_ab_min) %>%
    arrange(desc(lfc)) %>%
    slice_head(n = n_each),
  maaslin_tab %>%
    filter(diff, direction == "Down_in_Tumour", mean_rel > rel_ab_min) %>%
    arrange(lfc) %>%
    slice_head(n = n_each)
) %>%
  mutate(
    taxon_short = gsub("^s__", "", taxon),
    taxon_short = factor(taxon_short, levels = rev(taxon_short))
  )

pdf(file.path(DIR_RES, "DA_top10each_barplot_split.pdf"),
    width = 6, height = 4)
print(
  ggplot(top_da, aes(lfc, taxon_short, fill = direction)) +
    geom_col(width = 0.7) +
    geom_errorbarh(aes(xmin = lfc - se, xmax = lfc + se),
                   width = 0.25, colour = "grey30") +
    geom_vline(xintercept = 0, colour = "black") +
    scale_fill_manual(values = c(Up_in_Tumour = "#9A342C",
                                 Down_in_Tumour = "#126CAA")) +
    labs(x = "MaAsLin2 coefficient (mean +/- SE)",
         y = NULL, fill = NULL) +
    theme_test() +
    theme(panel.border = element_rect(fill = NA, colour = 1),
          axis.text = element_text(colour = 1, size = 8),
          legend.position = "top")
)
dev.off()

volc_df <- maaslin_tab %>%
  mutate(neg_log10_q = -log10(q + 1e-300)) %>%
  arrange(lfc)

top_labels <- volc_df %>%
  filter(diff, direction %in% c("Up_in_Tumour", "Down_in_Tumour")) %>%
  group_by(direction) %>%
  slice_max(abs(lfc), n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  pull(taxon)

volc_df <- volc_df %>%
  mutate(label = ifelse(taxon %in% top_labels,
                        gsub("^s__", "", taxon),
                        NA_character_))

pdf(file.path(DIR_RES, "DA_volcano.pdf"), width = 4.8, height = 3.5)
print(
  ggplot(volc_df, aes(lfc, neg_log10_q, colour = direction)) +
    geom_point(alpha = 0.75, size = 1.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
               linetype = "dashed", colour = "grey50") +
    ggrepel::geom_text_repel(aes(label = label), size = 2.6,
                             max.overlaps = 20, segment.size = 0.2,
                             show.legend = FALSE) +
    scale_colour_manual(values = c(Up_in_Tumour = "#9A342C",
                                   Down_in_Tumour = "#126CAA",
                                   NS = "grey80")) +
    labs(x = "MaAsLin2 coefficient (Tumour vs Normal)",
         y = "-log10(q)",
         colour = NULL) +
    theme_test() +
    theme(panel.border = element_rect(fill = NA, colour = 1),
          axis.text = element_text(colour = 1),
          legend.position = "right")
)
dev.off()

# Reuse the paired tumour-normal samples already defined for MaAsLin2
mtx_int <- round(mtx_count[, c(paste0("N", paired_ids), paste0("C", paired_ids))])

# Rarefaction curves to confirm sequencing depth saturates richness
pdf(file.path(DIR_RES, "S_rarefaction.pdf"), width = 4.5, height = 4.5)
par(mar = c(4, 4, 1, 1), pty = "s")
rc <- rarecurve(t(mtx_int), step = 500, label = FALSE, col = "grey75",
                xaxt = "n", xlab = "Read counts", ylab = "Observed species")
axis(1, at = axTicks(1),
     labels = formatC(axTicks(1), format = "f", digits = 0, big.mark = ","))
# Fit a Michaelis-Menten saturation model to all pooled (reads, richness) points
fit_df <- do.call(rbind, lapply(rc, function(v)
  data.frame(reads = attr(v, "Subsample"), richness = as.numeric(v))))
mm <- nls(richness ~ SSmicmen(reads, Smax, K), data = fit_df)
Smax <- coef(mm)[["Smax"]]; K <- coef(mm)[["K"]]
med_depth <- median(colSums(mtx_int))           # typical sequencing depth
capt_pct  <- med_depth / (med_depth + K)        # fraction of Smax captured
fit_med   <- Smax * med_depth / (K + med_depth)
xx <- seq(0, max(fit_df$reads), length.out = 400)
lines(xx, Smax * xx / (K + xx), col = "#9A342C", lwd = 2.5)
abline(h = Smax, lty = 2, col = "grey40")
text(0.30 * max(xx), Smax, pos = 1, cex = 0.7, col = "grey20",
     labels = paste0("Estimated ceiling  Smax = ", formatC(Smax, format = "f", digits = 0, big.mark = ",")))
points(med_depth, fit_med, pch = 19, col = "#9A342C")
text(med_depth, fit_med, pos = 4, cex = 0.7, col = "grey20",
     labels = paste0("median depth ", formatC(med_depth, format = "f", digits = 0, big.mark = ","), " reads  (~", round(100 * capt_pct), "% captured)"))
dev.off()

# Network topology metrics per group
node_full <- net_dt[["net.cor.matrix"]][["node"]]
edge_full <- net_dt[["net.cor.matrix"]][["edge"]]

build_grp_graph <- function(grp) {
  nd <- node_full[node_full[["label"]] == grp &
                    node_full[["igraph.degree"]] > 0 &
                    node_full[["Rank1"]] != "k__Viruses" &   # drop viral nodes
                    node_full[["Rank6"]] != "g__" &          # keep genus-resolved
                    node_full[["Rank7"]] != "s__", ]         # keep species-resolved (downstream focus is species)
  ed <- edge_full[edge_full[["label"]] == grp &
                    edge_full[["OTU_1"]] %in% nd[["ID"]] &
                    edge_full[["OTU_2"]] %in% nd[["ID"]], ]
  igraph::graph_from_data_frame(
    ed[, c("OTU_1", "OTU_2")],
    vertices = data.frame(name = nd[["ID"]]),
    directed = FALSE
  ) %>% igraph::simplify()
}
groups_net <- as.character(unique(node_full[["label"]]))
graphs_net <- setNames(lapply(groups_net, build_grp_graph), groups_net)

topo_df <- lapply(groups_net, function(g) {
  ig <- graphs_net[[g]]
  data.frame(
    Group      = ifelse(grepl("^Tumour", g), "Tumour", "Normal"),  # map by actual label (avoids prior positional swap)
    Nodes      = igraph::vcount(ig),
    Edges      = igraph::ecount(ig),
    Density    = igraph::edge_density(ig),
    AvgDegree  = mean(igraph::degree(ig)),
    Clustering = igraph::transitivity(ig, type = "global"),
    Efficiency = igraph::global_efficiency(ig)
  )
}) %>% do.call(rbind, .)
write.csv(topo_df, file.path(DIR_TAB, "Res_network_topology.csv"), row.names = FALSE)

pdf(file.path(DIR_RES, "S_network_topology.pdf"), width = 6.5, height = 4)
ggplot(
  reshape2::melt(topo_df, id.vars = "Group",
                 variable.name = "metric", value.name = "value"),
  aes(Group, value, fill = Group)
) +
  geom_col(width = 0.6) +
  geom_text(aes(label = signif(value, 3)), vjust = -0.4, size = 2.5) +
  scale_fill_manual(values = c(Normal = "#126CAA", Tumour = "#9A342C")) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2, labeller = labeller(metric = c(Nodes = "Nodes (count)", Edges = "Edges (count)", Density = "Density (0-1)", AvgDegree = "Mean degree (per node)", Clustering = "Clustering (0-1)", Efficiency = "Global efficiency (0-1)"))) +
  xlab(NULL) + ylab(NULL) +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 7),
        legend.position = "none")
dev.off()

# Mean abundance vs |effect size| for all DA species
mtx_rel <- sweep(mtx_count[, paired_samples, drop = FALSE], 2,
                 colSums(mtx_count[, paired_samples, drop = FALSE]), "/")

abund_df <- data.frame(
  taxon = rownames(mtx_count),
  mean_abundance = rowMeans(mtx_rel),
  prevalence = rowMeans(mtx_count[, paired_samples, drop = FALSE] > 0)
) %>%
  inner_join(maaslin_tab, by = "taxon") %>%
  filter(mean_abundance > rel_ab_min)   # only species with relative abundance > 0.01%

pdf(file.path(DIR_RES, "S_DA_abundance_effect.pdf"), width = 5.5, height = 3.6)
ggplot(abund_df, aes(mean_abundance, abs(lfc),
                     colour = direction, size = prevalence)) +
  geom_point(alpha = 0.7) +
  ggrepel::geom_text_repel(
    data = subset(abund_df,
                  abs(lfc) > log2(1.5) &
                    mean_abundance > 1e-3 &
                    prevalence > 0.6),
    aes(label = taxon),
    size = 2.5, max.overlaps = 15, show.legend = FALSE
  ) +
  scale_x_log10(labels = scales::label_log()) +
  scale_size_continuous(range = c(0.6, 3)) +
  scale_colour_manual(values = c(Up_in_Tumour   = "#9A342C",
                                 Down_in_Tumour = "#126CAA",
                                 NS             = "grey80")) +
  labs(x = "Mean relative abundance (log10)",
       y = "|MaAsLin2 coefficient|") +
  theme_test()
dev.off()


############################################################################
# High- vs Low-Shannon-diversity TUMOUR network comparison
# ---------------------------------------------------------------------------
# Mirrors the tumour-normal network comparison above, but the grouping is the
# WITHIN-tumour median split of the Shannon index (identical to the split used
# for the prognostic KM curve, L410-412). Interpretation focuses on
# size-independent structural metrics (density, clustering, global efficiency,
# natural-connectivity robustness, module conservation), since node/edge counts
# are partly driven by the diversity metric used to define the groups.
#
# Colour key: High diversity = deep blue-purple, Low diversity = light purple
# (kept distinct from the red/blue tumour-normal palette above).
#
# Self-contained AND idempotent: the guard below reloads lightweight inputs;
# every heavy step (network inference, robustness, module compare, MaAsLin2) is
# cached to RDS and re-loaded if present, so re-running this block only redraws
# the figures and never re-runs the tumour-normal pipeline above.
############################################################################
if (!exists("mtx_count")) {
  set.seed(42)
  suppressMessages({
    library(vegan); library(tidyverse); library(tidyfst); library(pulsar)
    library(ggpubr); library(phyloseq); library(igraph); library(ggClusterNet)
    library(ggraph); library(RColorBrewer); library(Maaslin2); library(tools)
    library(reshape2); library(readxl)
  })
  setwd("/data/yzwang/project/AEG_seiri/")
  DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
  DIR_RES <- "/data/yzwang/project/AEG_seiri/results/F2_species/"
  DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
  mtx_cpm   <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
  mtx_count <- readRDS(file.path(DIR_RDS, "sAEG_Count_RNA_FiltMyco.rds"))
  clinical  <- readxl::read_excel(file.path(DIR_TAB, "AEG_clinical.xlsx"))
}

DIR_DIV <- file.path(DIR_RES, "diversity_network/")
dir.create(DIR_DIV, showWarnings = FALSE, recursive = TRUE)

# Colour palette: High = deep blue-purple, Low = light purple
div_pal  <- c(High = "#5E4FA2", Low = "#C3B1E1")
dir_pal  <- c(Up_in_High = "#5E4FA2", Down_in_High = "#C3B1E1", NS = "grey80")
rand_pal <- c("network" = "#5E4FA2", "E–R network" = "grey70")

short_sp <- function(s) {
  ifelse(is.na(s) | s == "", "",
         paste0(toupper(substr(sub("_.*", "", s), 1, 1)), ".", sub("^[^_]*_", "", s)))
}

# --- Reproduce the EXACT High/Low split used for the KM curve -----------------
shan_div <- apply(mtx_count[, grepl("C", colnames(mtx_count))],
                  MARGIN = 2, FUN = vegan::diversity)
dist_div <- clinical[["Distance from the tumor center to the esophagogastric junction()"]]
names(dist_div) <- paste0("C", clinical$No.)
ov_div   <- intersect(names(shan_div), names(dist_div))
shan_div <- shan_div[ov_div]
div_group <- ifelse(shan_div > median(shan_div), "High", "Low")
names(div_group) <- ov_div
cat("High/Low diversity tumour sample sizes:\n"); print(table(div_group))

# --- Tumour-only phyloseq carrying the diversity group -----------------------
ps.div <- import_biom("./aeg.biom")
sample_data(ps.div)$ID <- sample_data(ps.div)$Id
keep_ids <- sample_names(ps.div)[as.character(sample_data(ps.div)$Id) %in% names(div_group)]
ps.div   <- prune_samples(keep_ids, ps.div)
sample_data(ps.div)$group <- div_group[as.character(sample_data(ps.div)$Id)]

# species-name lookups (used for ZiPi relabel and volcano bacteria filter)
tt <- as.data.frame(as.matrix(ps.div@tax_table))
tt$sp     <- paste0(gsub("g__", "", tt$Rank6), "_", gsub("s__", "", tt$Rank7))
name2king <- setNames(tt$Rank1, tt$sp)

# --- Network construction (cached) -------------------------------------------
rds_net <- file.path(DIR_RDS, "AEG_network_diversity.rds")
if (file.exists(rds_net)) {
  tab_div <- readRDS(rds_net)
} else {
  set.seed(42)
  tab_div <- network.pip(ps = ps.div, N = 200,
                         big = FALSE, select_layout = FALSE,
                         layout_net = "model_maptree2",
                         r.threshold = 0.6, p.threshold = 0.05,
                         maxnode = 2, method = "kendall",
                         label = FALSE, lab = "elements",
                         group = "group", fill = "Species",
                         size = "igraph.degree", zipi = TRUE,
                         ram.net = TRUE, clu_method = "cluster_fast_greedy",
                         step = 100, R = 10, ncpus = 6)
  saveRDS(tab_div, rds_net)
}
cortab_div <- tab_div[[2]]$net.cor.matrix$cortab
saveRDS(cortab_div, file.path(DIR_RDS, "AEG_network_diversity_cor.rds"))
net_dt_div <- tab_div[[2]]

# id -> species / kingdom maps from the network node table
node_tab <- net_dt_div$net.cor.matrix$node
node_tab$sp <- paste0(gsub("g__", "", node_tab$Rank6), "_", gsub("s__", "", node_tab$Rank7))
id2sp   <- setNames(node_tab$sp,    as.character(node_tab$ID))
id2king <- setNames(node_tab$Rank1, as.character(node_tab$ID))

# --- ZiPi (within-module connectivity z vs participation) --------------------
# Relabel hub/connector nodes with species names; drop viral nodes.
p_zipi <- tab_div[[1]][[2]]
for (li in c(2, 3)) {
  d <- p_zipi$layers[[li]]$data
  d$.king <- id2king[as.character(d$id)]
  d$.sp   <- id2sp[as.character(d$id)]
  d <- d[is.na(d$.king) | d$.king != "k__Viruses", ]      # drop viruses
  if (li == 3) d$label <- ifelse(d$roles != "Peripherals", short_sp(d$.sp), "")
  p_zipi$layers[[li]]$data <- d
}
pdf(file.path(DIR_DIV, "Net_ZiPi_div.pdf"), width = 6, height = 5)   # was Net_connectivity_div.pdf
print(p_zipi +
        theme(panel.border = element_rect(fill = NA, colour = 1),
              axis.text = element_text(colour = 1)))
dev.off()
if (file.exists(file.path(DIR_DIV, "Net_connectivity_div.pdf")))
  file.remove(file.path(DIR_DIV, "Net_connectivity_div.pdf"))  # misnamed earlier (it is the ZiPi plot)

# --- Network randomness (degree distribution: observed vs Erdos-Renyi) --------
pdf(file.path(DIR_DIV, "Net_randomness_div.pdf"), width = 6, height = 5)
print(tryCatch(
  tab_div[[1]][[3]] +
    scale_colour_manual(values = rand_pal) +
    scale_fill_manual(values = rand_pal) +
    labs(x = "Degree", y = "Frequency") +
    theme(panel.border = element_rect(fill = NA, colour = 1),
          axis.text = element_text(colour = 1)),
  error = function(e) tab_div[[1]][[3]]))
dev.off()

# --- Network destruction resistance (cached) ---------------------------------
rds_resis <- file.path(DIR_RDS, "network_destruction_resistance_diversity.rds")
if (file.exists(rds_resis)) {
  resis_div <- readRDS(rds_resis)
} else {
  resis_div <- natural.con.microp(ps = ps.div, corg = cortab_div,
                                  norm = TRUE, end = 150, start = 0)
  saveRDS(resis_div, rds_resis)
}

pdf(file.path(DIR_DIV, "Net_resistance_div.pdf"), width = 5, height = 4)
print(resis_div[[1]] +
        scale_colour_manual(values = div_pal) +
        labs(x = "Number of removed nodes", y = "Natural connectivity") +
        theme(panel.border = element_rect(fill = NA, colour = 1),
              axis.text = element_text(colour = 1)))
dev.off()
write.csv(resis_div[[2]], file.path(DIR_TAB, "Res_network_resistance_diversity.csv"))

resis_div_data <- resis_div[[2]] %>%
  group_by(Group) %>%
  mutate(baseline = Natural.connectivity[1],
         changed_connectivity = Natural.connectivity / baseline) %>%
  ungroup()

pdf(file.path(DIR_DIV, "Net_resistance_changed_div.pdf"), width = 5, height = 4)
print(ggplot(resis_div_data, aes(x = `Num.of.remove.nodes`,
                                 y = changed_connectivity, colour = Group)) +
        geom_point(alpha = 0.3) +
        xlim(0, 100) +
        geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
        scale_colour_manual(values = div_pal) +
        labs(x = "Number of removed nodes",
             y = "Connectivity - relative change") +
        theme_minimal() +
        theme(panel.border = element_rect(fill = NA, colour = 1),
              axis.text = element_text(colour = 1)))
dev.off()

# --- Module similarity / conservation (cached) -------------------------------
rds_mod <- file.path(DIR_RDS, "AEG_network_module_diversity.rds")
if (file.exists(rds_mod)) {
  module_div <- readRDS(rds_mod)
} else {
  module_div <- module.compare.m(ps = ps.div, corg = cortab_div, zipi = FALSE,
                                 zoom = 0.2, padj = F, n = 3)
  saveRDS(module_div, rds_mod)
}

pdf(file.path(DIR_DIV, "Net_module_div.pdf"), width = 4, height = 4)
print(module_div[[1]])
dev.off()

module_otu_div <- module_div[[2]]
module_otu_div$taxa_g <- ps.div@tax_table[module_otu_div$ID, "Rank6"]
module_otu_div$taxa_s <- ps.div@tax_table[module_otu_div$ID, "Rank7"]
module_otu_div$taxa_k <- ps.div@tax_table[module_otu_div$ID, "Rank1"]
module_otu_div <- module_otu_div %>%
  filter(taxa_k != "k__Viruses", taxa_g != "g__", taxa_s != "s__") %>%
  mutate(standard_name = paste0(gsub("g__", "", taxa_g) %>% toTitleCase(.), "_",
                                gsub("s__", "", taxa_s)))

sme_div <- module_otu_div %>%
  select(ID, group) %>% distinct() %>%
  mutate(node_id = paste0(ID, "__", group))

shared_ids_div <- sme_div %>%
  group_by(ID) %>%
  summarise(has_high = any(grepl("High", group)),
            has_low  = any(grepl("Low", group)),
            .groups = "drop") %>%
  filter(has_high & has_low) %>%
  pull(ID)

nodes_div <- sme_div %>%
  left_join(module_otu_div %>% select(ID, standard_name) %>% distinct(), by = "ID") %>%
  mutate(node_type = case_when(
    ID %in% shared_ids_div ~ "Conserved",
    grepl("High", group)   ~ "High-enriched",
    grepl("Low", group)    ~ "Low-enriched"
  ))

edges_within_div <- sme_div %>%
  inner_join(sme_div, by = "group", relationship = "many-to-many") %>%
  filter(node_id.x < node_id.y) %>%
  select(from = node_id.x, to = node_id.y)

edges_shared_div <- sme_div %>%
  filter(ID %in% shared_ids_div) %>%
  inner_join(sme_div %>% filter(ID %in% shared_ids_div),
             by = "ID", relationship = "many-to-many") %>%
  filter(group.x < group.y) %>%
  select(from = node_id.x, to = node_id.y) %>%
  distinct()

edges_all_div <- bind_rows(edges_within_div, edges_shared_div) %>%
  distinct() %>%
  mutate(edge_type = if_else(
    paste0(from, to) %in% paste0(edges_shared_div$from, edges_shared_div$to),
    "cross-model", "within-model"))

g_div <- graph_from_data_frame(
  edges_all_div, directed = FALSE,
  vertices = nodes_div %>% select(node_id, standard_name, group, node_type) %>%
    rename(name = node_id))

node_colours_div <- c("High-enriched" = "#5E4FA2",
                      "Low-enriched"  = "#C3B1E1",
                      "Conserved"     = "#9E9AC8")
layout_div <- create_layout(g_div, layout = 'kk')
pdf(file.path(DIR_DIV, "C_net_module_revisualised_div.pdf"), width = 5, height = 4)
print(ggraph(layout_div) +
        geom_edge_link(aes(linetype = edge_type), alpha = 0.5,
                       colour = "grey", width = 0.3) +
        scale_edge_linetype_manual(values = c("within-model" = "solid",
                                              "cross-model" = "dashed")) +
        geom_node_point(aes(colour = node_type), size = 1.5) +
        scale_colour_manual(values = node_colours_div, name = "Node type") +
        theme_void() +
        theme(legend.position = "right", legend.title = element_blank()))
dev.off()

# --- Topology metrics per diversity group ------------------------------------
node_full_div <- net_dt_div[["net.cor.matrix"]][["node"]]
edge_full_div <- net_dt_div[["net.cor.matrix"]][["edge"]]

build_grp_graph_div <- function(grp) {
  nd <- node_full_div[node_full_div[["label"]] == grp &
                        node_full_div[["igraph.degree"]] > 0 &
                        node_full_div[["Rank1"]] != "k__Viruses" &
                        node_full_div[["Rank6"]] != "g__" &
                        node_full_div[["Rank7"]] != "s__", ]
  ed <- edge_full_div[edge_full_div[["label"]] == grp &
                        edge_full_div[["OTU_1"]] %in% nd[["ID"]] &
                        edge_full_div[["OTU_2"]] %in% nd[["ID"]], ]
  igraph::graph_from_data_frame(
    ed[, c("OTU_1", "OTU_2")],
    vertices = data.frame(name = nd[["ID"]]),
    directed = FALSE) %>% igraph::simplify()
}
groups_net_div <- as.character(unique(node_full_div[["label"]]))
graphs_net_div <- setNames(lapply(groups_net_div, build_grp_graph_div), groups_net_div)

topo_df_div <- lapply(groups_net_div, function(g) {
  ig <- graphs_net_div[[g]]
  data.frame(
    Group      = ifelse(grepl("^High", g), "High", "Low"),
    Nodes      = igraph::vcount(ig),
    Edges      = igraph::ecount(ig),
    Density    = igraph::edge_density(ig),
    AvgDegree  = mean(igraph::degree(ig)),
    Clustering = igraph::transitivity(ig, type = "global"),
    Efficiency = igraph::global_efficiency(ig))
}) %>% do.call(rbind, .)
write.csv(topo_df_div, file.path(DIR_TAB, "Res_network_topology_diversity.csv"),
          row.names = FALSE)

pdf(file.path(DIR_DIV, "S_network_topology_div.pdf"), width = 6.5, height = 4)
print(ggplot(
  reshape2::melt(topo_df_div, id.vars = "Group",
                 variable.name = "metric", value.name = "value"),
  aes(Group, value, fill = Group)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = signif(value, 3)), vjust = -0.4, size = 2.5) +
    scale_fill_manual(values = div_pal) +
    facet_wrap(~ metric, scales = "free_y", nrow = 2,
               labeller = labeller(metric = c(
                 Nodes = "Nodes (count)", Edges = "Edges (count)",
                 Density = "Density (0-1)", AvgDegree = "Mean degree (per node)",
                 Clustering = "Clustering (0-1)",
                 Efficiency = "Global efficiency (0-1)"))) +
    xlab(NULL) + ylab(NULL) +
    theme_test() +
    theme(panel.border = element_rect(fill = NA, colour = 1),
          axis.text = element_text(colour = 1),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.background = element_blank(),
          strip.text = element_text(size = 7),
          legend.position = "none"))
dev.off()

# --- Differential abundance: High vs Low (MaAsLin2, unpaired; cached) --------
rel_ab_min <- 1e-4
div_samples <- names(div_group)
meta_div <- data.frame(
  sample = div_samples,
  group  = factor(div_group[div_samples], levels = c("Low", "High")),
  row.names = div_samples)

rds_maas <- file.path(DIR_RDS, "MaAsLin2_HighLow.rds")
if (file.exists(rds_maas)) {
  maaslin_res_div <- readRDS(rds_maas)
} else {
  maaslin_dir_div <- file.path(DIR_DIV, "MaAsLin2_HighLow")
  dir.create(maaslin_dir_div, showWarnings = FALSE, recursive = TRUE)
  set.seed(42)
  maaslin_res_div <- Maaslin2(
    input_data = as.data.frame(mtx_count[, div_samples, drop = FALSE]),
    input_metadata = meta_div, output = maaslin_dir_div,
    fixed_effects = "group", reference = "group,Low",
    normalization = "TSS", transform = "LOG", analysis_method = "LM",
    min_prevalence = 0.10, min_abundance = 0.0, max_significance = 0.05,
    plot_heatmap = FALSE, plot_scatter = FALSE, cores = 4)
  saveRDS(maaslin_res_div, rds_maas)
}

maaslin_tab_div <- maaslin_res_div$results %>%
  filter(metadata == "group") %>%
  transmute(
    taxon = feature, lfc = coef, se = stderr, p = pval, q = qval,
    diff = qval < 0.05,
    direction = case_when(
      qval < 0.05 & coef > 0 ~ "Up_in_High",
      qval < 0.05 & coef < 0 ~ "Down_in_High",
      TRUE ~ "NS"))
maaslin_tab_div$mean_rel <- rowMeans(
  sweep(mtx_count[, div_samples, drop = FALSE], 2,
        colSums(mtx_count[, div_samples, drop = FALSE]), "/")
)[maaslin_tab_div$taxon]
write.csv(maaslin_tab_div,
          file.path(DIR_TAB, "Res_DA_HighLow_diversity.csv"), row.names = FALSE)

# Composite ranking: among significant, abundant taxa, balance effect size and
# abundance via composite = |lfc| * log10(mean_rel / floor) (larger = stronger
# AND more abundant). Top 10 per direction.
n_each <- 10
comp_pool <- maaslin_tab_div %>%
  filter(diff, mean_rel > rel_ab_min) %>%
  mutate(composite = abs(lfc) * log10(mean_rel / rel_ab_min))

top_da_div <- bind_rows(
  comp_pool %>% filter(direction == "Up_in_High") %>%
    arrange(desc(composite)) %>% slice_head(n = n_each),
  comp_pool %>% filter(direction == "Down_in_High") %>%
    arrange(desc(composite)) %>% slice_head(n = n_each)) %>%
  arrange(lfc) %>%
  mutate(taxon_short = gsub("^s__", "", taxon),
         taxon_short = factor(taxon_short, levels = taxon_short))

pdf(file.path(DIR_DIV, "DA_top10each_barplot_split_div.pdf"), width = 6, height = 4)
print(ggplot(top_da_div, aes(lfc, taxon_short, fill = direction)) +
        geom_col(width = 0.7) +
        geom_errorbarh(aes(xmin = lfc - se, xmax = lfc + se),
                       height = 0.25, colour = "grey30") +
        geom_vline(xintercept = 0, colour = "black") +
        scale_fill_manual(values = dir_pal) +
        labs(x = "MaAsLin2 coefficient (mean +/- SE)",
             y = NULL, fill = NULL) +
        theme_test() +
        theme(panel.border = element_rect(fill = NA, colour = 1),
              axis.text = element_text(colour = 1, size = 8),
              legend.position = "top"))
dev.off()

# --- Volcano plot ------------------------------------------------------------
# Colour a point only if it clears BOTH thresholds: q < 0.05 AND |coef| > log2(1.5)
# (i.e. it must also pass the vertical effect-size lines); everything else grey.
# Labels: the 5 most-significant (smallest q) BACTERIAL species among those that
# clear both thresholds.
coef_thr <- log2(1.5)
volc_df <- maaslin_tab_div %>%
  mutate(neg_log10_q = -log10(q + 1e-300),
         kingdom = name2king[taxon],
         is_bac  = !is.na(kingdom) & kingdom == "k__Bacteria",
         sig     = diff & abs(lfc) > coef_thr,
         vcol    = case_when(
           sig & lfc > 0 ~ "Up_in_High",
           sig & lfc < 0 ~ "Down_in_High",
           TRUE          ~ "NS"))
lab_taxa <- volc_df %>%
  filter(sig, is_bac) %>%
  slice_min(q, n = 5, with_ties = FALSE) %>%
  pull(taxon)
volc_df <- volc_df %>%
  mutate(label = ifelse(taxon %in% lab_taxa, gsub("_", " ", taxon), NA_character_))

pdf(file.path(DIR_DIV, "DA_volcano_div.pdf"), width = 4.8, height = 3.5)
print(ggplot(volc_df, aes(lfc, neg_log10_q, colour = vcol)) +
        geom_point(alpha = 0.75, size = 1.6) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
        geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
                   linetype = "dashed", colour = "grey50") +
        ggrepel::geom_text_repel(aes(label = label), size = 2.6,
                                 max.overlaps = 20, segment.size = 0.2,
                                 show.legend = FALSE) +
        scale_colour_manual(values = dir_pal) +
        labs(x = "MaAsLin2 coefficient (High vs Low diversity)",
             y = "-log10(q)", colour = NULL) +
        theme_test() +
        theme(panel.border = element_rect(fill = NA, colour = 1),
              axis.text = element_text(colour = 1),
              legend.position = "right"))
dev.off()

cat("\nLabelled volcano bacteria (5 smallest q):\n"); print(lab_taxa)
cat("\nComposite-ranked DA barplot taxa:\n"); print(as.character(top_da_div$taxon_short))
cat("\n=== Diversity-network comparison (restyled) finished. Outputs in:", DIR_DIV, "===\n")
