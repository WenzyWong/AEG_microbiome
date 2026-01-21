#######################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 04. Analyse AEG-specific microbiome
#
#######################################################
#remotes::install_github("taowenmicro/EasyStat")
#remotes::install_github("taowenmicro/ggClusterNet")
#remotes::install_github("zdk123/SpiecEasi")

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
# Regression
library(compositions)
library(glmnet)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/F2/"
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
gAbund <- sort(apply(mtx_gcpm, MARGIN = 1, FUN = mean) / 1e+4, decreasing = T)
gAbund <- gAbund[1:20]
indexSp <- c() # The species index in the cpm matrix, assigned to each genus
for (g in names(gAbund)) {
  pairSp <- grep(paste("*", g, "*", sep = ""), rownames(mtx_cpm))
  indexSp <- c(indexSp, pairSp)
}

abundSpTop <- apply(mtx_cpm[indexSp, ] / 1e+4, MARGIN = 1, FUN = mean)
abundSpTop <- data.frame(
  Species = gsub("[[:punct:]]", "_", names(abundSpTop)), # Adjusting the species names to successfully establish formulas
  Abundance = apply(mtx_cpm[indexSp, ] / 1e+4, MARGIN = 1, FUN = mean)
) # Species within top30 genera, no matter how abundant they really are
# Notice: the rownames of abundSpTop is different from abundSpTop$Species.
# Using species name in the latter steps, you must pay attention to the choice.

abundSpTop <- abundSpTop[abundSpTop$Abundance > 0.01, ]

# drawHR1 contains the clinical information of the samples
drawHR1 <- data.frame(
  sample = paste0("C", clinical$No.),
  time = clinical$`Survival time（Month）`,
  state = clinical$`Status（1=Dead, 0=Alive）`
) %>%
  filter(sample %in% colnames(mtx_cpm))

# drawHR2 contains the log2CPM values of species
drawHR2 <- as.data.frame(t(log2(mtx_cpm[rownames(abundSpTop), drawHR1$sample] + 1)))
colnames(drawHR2) <- abundSpTop$Species
drawHR2$sample <- rownames(drawHR2) # Used for merge only

drawHR <- merge(drawHR1, drawHR2, by = "sample")

# Calculating HR
source(file.path(DIR_TOOL, "hr_calc.R"))
hrRes <- hr_calc(abundSpTop$Species, drawHR, 1.1, 0.9, 0.05)
# The warnings exist because all the tumour CPM values of these species equals to 0

abundSpTop$Surv.HR <- hrRes$HR
abundSpTop$Surv.P <- hrRes$Pvalue
abundSpTop$Surv.Risk <- hrRes$Risk
abundSpTop$Surv.HR[abundSpTop$Surv.HR > 3] <- 3
# Removing the species with 0 count values in tumour samples, causing HR equals to NA
abundSpTop <- na.omit(abundSpTop)

lenEachG <- c() # The number of identified species within each abundant genus
for (g in names(gAbund)) {
  # This temporary variant "pairSp" is different from the "pairSp" in the first loop
  pairSp <- grep(paste("*", g, "*", sep = ""), abundSpTop$Species) 
  lenEachG <- c(lenEachG, length(pairSp))
}
names(lenEachG) <- names(gAbund)

# Allocate positions for species-pairing genera
tmpGenera <- c() # Species-pairing genera
tmpX <- c() # Allocate x-coordinates by the number of indentified species within each genus
for (i in c(1:length(gAbund))) {
  tmpGenera <- c(tmpGenera, rep(names(gAbund)[i], lenEachG[i]))
  tmpX <- c(tmpX, seq(0, 1, length.out = lenEachG[i]))
}
abundSpTop$Genus <- tmpGenera
abundSpTop$X <- tmpX

# Calculating differential
source(file.path(DIR_TOOL, "wilcox_diff.R"))
diffSp <- wilcox_diff(mtx_cpm[gsub("[[:punct:]]", "_", rownames(mtx_cpm)) %in% 
                                abundSpTop$Species, ],
                      mtx_cpm[gsub("[[:punct:]]", "_", rownames(mtx_cpm)) %in% 
                                abundSpTop$Species, grep("C", colnames(mtx_cpm))],
                      mtx_cpm[gsub("[[:punct:]]", "_", rownames(mtx_cpm)) %in% 
                                abundSpTop$Species, grep("N", colnames(mtx_cpm))],
                      T, 0.05, log2(1.5))
colnames(diffSp)[1] <- "Species"
diffSp <- diffSp[order(diffSp$Species, decreasing = T), ]

# Pairing differential trend toward species pending circlize plot
abundSpTop$Diff.Trend <- diffSp$Change
abundSpTop$Diff.Padj <- diffSp$P.adj
abundSpTop$Diff.Log2FC <- diffSp$log2FC

saveRDS(abundSpTop, file.path(DIR_RDS, "sAEG_CirclizeData_AbundSpTop_Genera20.rds"))

colGenera <- paste0(substr(paletteer_d("khroma::discreterainbow")[c(10, 12:20, 
                                                                    23:27, 2, 4, 5, 7, 9)], 
                           1, 7), "80") # Colours for genera. "80" reprensents alpha
colSp <- substr(paletteer_d("khroma::discreterainbow")[c(10, 12:20, 
                                                         23:27, 2, 4, 5, 7, 9)], 
                1, 7)

# Assigning the scaled colours for survival HR risks
# Calculating the assignment portion within colour gradiant
lenHRl <- length(rownames(abundSpTop)[abundSpTop$Surv.HR > 1.1])
lenHRs <- length(rownames(abundSpTop)[abundSpTop$Surv.HR < 0.9])
colourHRl <- colorRampPalette(c("#CC6329", "#CDCDCD"))(lenHRl)
colourHRs <- colorRampPalette(c("#409161", "#CDCDCD"))(lenHRs)
colourHR <- c()
cntl <- 0
cnts <- 0
orderHR <- abundSpTop[order(abs(abundSpTop$Surv.HR), decreasing = T), ]
for (i in 1:nrow(orderHR)) {
  if (orderHR$Surv.HR[i] > 1.1) {
    cntl <- cntl + 1
    colourHR <- c(colourHR, colourHRl[cntl])
  } else if (orderHR$Surv.HR[i] < 0.9) {
    cnts <- cnts + 1
    colourHR <- c(colourHR, colourHRs[cnts])
  } else {
    colourHR <- c(colourHR, "#EDEDED")
  }
}
names(colourHR) <- rownames(orderHR)
colourHR <- colourHR[rownames(abundSpTop)]

# Assigning the scaled colours for differential log2 fold changes
# Calculating the assignment portion within colour gradiant
lenLog2FCl <- length(rownames(abundSpTop)[abundSpTop$Diff.Log2FC > log2(1.5)])
lenLog2FCs <- length(rownames(abundSpTop)[abundSpTop$Diff.Log2FC < -log2(1.5)])
colourLog2FCl <- colorRampPalette(c("#9A342C", "#CDCDCD"))(lenLog2FCl)
colourLog2FCs <- colorRampPalette(c("#2C6199", "#CDCDCD"))(lenLog2FCs)

colourLog2FC <- c()
cntl <- 0
cnts <- 0
orderDiff <- abundSpTop[order(abs(abundSpTop$Diff.Log2FC), decreasing = T), ]
for (i in 1:nrow(orderDiff)) {
  if (orderDiff$Diff.Log2FC[i] > log2(1.5)) {
    cntl <- cntl + 1
    colourLog2FC <- c(colourLog2FC, colourLog2FCl[cntl])
  } else if (orderDiff$Diff.Log2FC[i] < -log2(1.5)) {
    cnts <- cnts + 1
    colourLog2FC <- c(colourLog2FC, colourLog2FCs[cnts])
  } else {
    colourLog2FC <- c(colourLog2FC, "#EDEDED")
  }
}
names(colourLog2FC) <- rownames(orderDiff)
colourLog2FC <- colourLog2FC[abundSpTop$Species]

# Drawing the most abundant genera section
# Allocating the canvas and sectors
pdf(file.path(DIR_RES, "A_circlised_genus_species_filtered.pdf"), width = 8, height = 8)
par(mar = c(1, 1, 1, 1) * 11, cex = 0.6, xpd = NA)
sectors <- factor(names(gAbund), levels = names(gAbund))
circos.par(points.overflow.warning = FALSE,
           cell.padding = c(0, 0, 0, 0))
circos.initialize(factors = sectors, xlim = c(0, 1), sector.width = lenEachG)
circos.trackPlotRegion(factors = sectors, ylim = c(0, 12), 
                       track.height = 0.5, bg.border = "grey",
                       bg.col = colGenera)

# Drawing the species through loops
cnt <- 0
for (i in 1:length(gAbund)) {
  # Get species indices for current genus and sort by abundance
  sp_indices <- which(abundSpTop$Genus == sectors[i])
  sp_order <- sp_indices[order(abundSpTop$Abundance[sp_indices], decreasing = TRUE)]
  
  # Calculate uniform width for each species within this genus
  bar_width <- 1 / lenEachG[i] / 2  # Half width of each position
  
  for (j in 1:lenEachG[i]) {
    cnt <- cnt + 1
    sp_idx <- sp_order[j]
    
    # Calculate the center position for this species
    x_center <- (j - 0.5) / lenEachG[i]
    
    # The length of each species-bar represents its relative abundance (log-transformed for better visualization)
    circos.rect(xleft = x_center - bar_width * 0.8,
                ybottom = 0,
                xright = x_center + bar_width * 0.8,
                ytop = log10(abundSpTop$Abundance[sp_idx] + 1) * 
                  (12 / log10(max(abundSpTop$Abundance) + 1)),
                sector.index = sectors[i],
                col = colSp[i],
                border = NA)
    # The annotation circle representing survival risks of species
    circos.trackLines(sectors = sectors[i],
                      x = x_center,
                      y = -1,
                      col = colourHR[sp_idx],
                      type = "h",
                      baseline = -2)
    # The annotation sub-circle representing survival significance (Surv.P)
    circos.trackLines(sectors = sectors[i],
                      x = x_center,
                      y = if_else(
                        abundSpTop$Surv.Risk[sp_idx] != "NS", -2.4, -2.3
                      ),
                      col = if_else(
                        abundSpTop$Surv.Risk[sp_idx] != "NS", "black", "white"
                      ),
                      type = "h",
                      baseline = -2.3)
    # The annotation circle representing differential log2foldchange between paired tumour & normal samples
    circos.trackLines(sectors = sectors[i],
                      x = x_center,
                      y = -4,
                      col = colourLog2FC[sp_idx],
                      type = "h",
                      baseline = -3)
    # The annotation sub-circle representing differential significance (Diff.Padj)
    circos.trackLines(sectors = sectors[i],
                      x = x_center,
                      y = if_else(
                        abundSpTop$Diff.Trend[sp_idx] != "NS", -4.4, -4.3
                      ),
                      col = if_else(
                        abundSpTop$Diff.Trend[sp_idx] != "NS", "black", "white"
                      ),
                      type = "h",
                      baseline = -4.3)
    
  }
  # The names of sectors (genera)
  circos.trackText(sectors = sectors[i],
                   x = 0.5, y = 7, cex = fontsize(14),
                   labels = names(gAbund)[i],
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
ggplot(data = abundSpTop, aes(colour = Surv.HR, 
                              x = Species, y = X)) +
  geom_point() + 
  scale_color_gradient2(aes(labs = c(min(Surv.HR), 1, max(Surv.HR))),
                        midpoint = 1,  
                        low = "#409161", 
                        mid = "#CDCDCD", 
                        high = "#CC6329")
dev.off()

pdf(file.path(DIR_RES, "ScaleDiff_Legend.pdf"), width = 6, height = 5)
ggplot(data = abundSpTop, aes(colour = Diff.Log2FC, 
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
mtx_count_t <- mtx_count[ , grepl("C", colnames(mtx_count))]
shan_tumour <- apply(mtx_count_t, 2, vegan::diversity)
dist_tumour <- clinical$`Distance from the tumor center to the esophagogastric junction()`
names(dist_tumour) <- paste0("C", clinical$No.)

overlap_samples <- intersect(names(shan_tumour), names(dist_tumour))
shan_tumour <- shan_tumour[overlap_samples]
dist_tumour <- dist_tumour[overlap_samples]

cor.test(shan_tumour, dist_tumour) # NS

summary(dist_tumour)

dist_groups <- data.frame(
  sample = overlap_samples,
  distance = dist_tumour,
  shannon = shan_tumour
) %>%
  mutate(dist_group = if_else(distance > median(distance), "Close", "Far"))

pdf(file.path(DIR_RES, "Test_alpha_distance.pdf"), width = 4, height = 3.2)
ggboxplot(dist_groups, x = "dist_group", y = "shannon",
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

siewart_groups <- data.frame(
  sample = overlap_samples,
  siewart = siewart_tumour,
  shannon = shan_tumour
) %>%
  mutate(siewart = factor(siewart, levels = c("I", "II", "III")))
pdf(file.path(DIR_RES, "Test_alpha_siewart.pdf"), width = 4, height = 3.2)
ggboxplot(siewart_groups, x = "siewart", y = "shannon",
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
source(file.path(DIR_TOOL, "hr_calc.R"))
df_surv <- data.frame(
  sample = paste0("C", clinical$No.),
  time = clinical$`Survival time（Month）`,
  state = clinical$`Status（1=Dead, 0=Alive）`
) %>%
  filter(sample %in% overlap_samples) %>%
  arrange(., sample)
shan_tumour <- shan_tumour[df_surv$sample]
df_surv$sample == names(shan_tumour)

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
hr_shan

pdf(file.path(DIR_RES, "B_alpha_km.pdf"), width = 3.8, height = 4)
surv_res$plot +
  ggplot2::annotate(
    "text",
    x = 140, y = 0.95,
    vjust = 1, hjust = 1,
    label = "HR = 0.439\np = 0.026",
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
rank_cols <- paste0("Rank", 1:7)
for (col in rank_cols) {
  node[[col]] <- ifelse(grepl("^[a-z]__$", node[[col]]), 
                        paste0(sub("__$", "", node[[col]]), "__unclassified"), 
                        node[[col]])
}

node_highlight <- node[node$ID %in% c("1351"), ] %>%
  mutate(standard_name = paste0(gsub("g__", "", Rank6) %>% 
                                  substr(., 1, 1) %>% 
                                  toupper(.), ".", 
                                gsub("s__", "", Rank7)),)

pdf(file.path(DIR_RES, "C_net.pdf"), width = 12, height = 5)
ggplot() + 
  geom_segment(data = edge, alpha = 0.3,
               aes(x = X1, y = Y1, xend = X2, yend = Y2, 
                   linewidth = weight)) +
  scale_linewidth_continuous(range = c(0.01, 0.05)) +
  geom_point(data = node, pch = 21, color = "gray40",
             aes(X1, X2, fill = Rank2, size = igraph.degree)) +
  scale_fill_manual(values = paletteer_d("ggsci::nrc_npg")) +
  facet_wrap(.~ label, scales = "free_y", nrow = 1) +
  geom_text(data = node_highlight, aes(X1, X2, label = standard_name)) +
  scale_size(range = c(0.8, 5)) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_rect(colour = NA),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())
dev.off()

# Other plots
plots <- tab[[1]]

pdf(file.path(DIR_RES, "C_net_connectivity_tn.pdf"), width = 6, height = 5)
plots[[2]]
dev.off()

pdf(file.path(DIR_RES, "C_net_randomness_tn.pdf"), width = 6, height =5)
plots[[3]] +
  scale_colour_manual(values = c("#0073C2FF", "#EFC000FF")) +
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF")) +
  labs(x = NULL,
       y = NULL) +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))
dev.off()

# Network destruction resistance
resis <- natural.con.microp(ps = ps.obj, corg = cortab,
                            norm = TRUE, end = 150, start = 0)
saveRDS(resis, file.path(DIR_RDS, "network_destruction_resistance.rds"))

pdf(file.path(DIR_RES, "C_net_resistance.pdf"), width = 5, height = 4)
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

pdf(file.path(DIR_RES, "C_net_resistance_changed_0to100.pdf"), width = 5, height = 4)
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
module <- module.compare.m(ps = NULL, corg = cortab, zipi = FALSE,
                           zoom = 0.2, padj = F, n = 3)
saveRDS(module, file.path(DIR_RDS, "AEG_network_module.rds"))

pdf(file.path(DIR_RES, "C_net_module.pdf"), width = 4, height = 4)
module[[1]]
dev.off()

module_otu <- module[[2]]
module_otu$taxa_g <- ps.obj@tax_table[module_otu$ID, "Rank6"]
module_otu$taxa_s <- ps.obj@tax_table[module_otu$ID, "Rank7"]
module_otu <- module_otu %>%
  filter(taxa_g != "g__",
         taxa_s != "s__") %>%
  mutate(standard_name = paste0(gsub("g__", "", taxa_g)  %>% 
                                  toTitleCase(.), "_", 
                                gsub("s__", "", taxa_s)),)

#mutate(standard_name = paste0(gsub("g__", "", taxa_g) %>% 
#                               substr(., 1, 1) %>% 
#                               toupper(.), ".", 
#                              gsub("s__", "", taxa_s)),)

head(module_otu)
table(module_otu$group)

# Re-visualise the module relationships
nodes <- module_otu %>%
  select(ID, standard_name) %>%
  distinct(ID, .keep_all = TRUE)
species_modules <- module_otu %>%
  select(ID, group) %>%
  distinct()

edges <- bind_rows(
  # Within-module connections
  species_modules %>%
    inner_join(species_modules, by = "group", relationship = "many-to-many") %>%
    filter(ID.x < ID.y) %>%
    select(from = ID.x, to = ID.y),
  
  # Cross-module connections (shared species)
  species_modules %>%
    inner_join(species_modules, by = "ID", relationship = "many-to-many") %>%
    filter(group.x != group.y, group.x < group.y) %>%
    select(from = ID, to = ID) %>%
    distinct()
) %>%
  distinct() %>%
  mutate(weight = 1)

nodes <- nodes %>%
  left_join(
    species_modules %>%
      group_by(ID) %>%
      summarise(
        # primary group: if normal exist, keep normal model number
        primary_group = {
          normal_groups <- group[grepl("Normal", group)]
          if (length(normal_groups) > 0) {
            paste(normal_groups, collapse = "; ")
          } else {
            paste(group, collapse = "; ")
          }
        },
        # secondary group: if tumour exist, keep tumour model number
        secondary_group = {
          tumour_groups <- group[grepl("Tumour", group)]
          if (length(tumour_groups) > 0) {
            paste(tumour_groups, collapse = "; ")
          } else {
            paste(group, collapse = "; ")
          }
        },
        group = paste(group, collapse = "; "),
        n_modules = n(),
        has_normal = any(grepl("Normal", group)),
        has_tumour = any(grepl("Tumour", group))
      ),
    by = "ID"
  )

g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

set.seed(42)
layout <- create_layout(g, layout = 'fr')

# Colours
blues_extended <- c(paletteer_d("RColorBrewer::Blues")[c(9)],
                    rev(paletteer_d("ggsci::blue_material")))
reds_extended <- rev(brewer.pal(9, "Reds")[2:6])

colour_map <- c(
  setNames(blues_extended, paste0("Normalmodel_", 1:11)),
  setNames(reds_extended, paste0("Tumourmodel_", 1:5))
)

pdf(file.path(DIR_RES, "C_net_module_revisulised.pdf"), width = 5, height = 4)
ggraph(layout) +
  geom_edge_link(alpha = 0.5, colour = "grey", width = 0.3) +
  # Main node color
  geom_node_point(aes(colour = primary_group), size = 3.5) +
  # Outer ring for shared nodes
  geom_node_point(data = function(x) filter(x, has_normal & has_tumour),
                  aes(x = x, y = y, fill = primary_group, colour = secondary_group), 
                  shape = 21, size = 3.8) +
  scale_colour_manual(values = colour_map) +
  scale_fill_manual(values = colour_map) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank())
dev.off()

# For legend saving
pdf(file.path(DIR_RES, "C_net_module_legend.pdf"), width = 5, height = 4)
ggraph(layout) +
  geom_edge_link(alpha = 0.5, colour = "grey", width = 0.3) +
  # Main node color
  geom_node_point(aes(colour = primary_group), size = 3.5) +
  scale_colour_manual(values = colour_map) +
  scale_fill_manual(values = colour_map) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank())
dev.off()

module_simi <- module[[3]]

module_simi$m1 <- module_simi$module1 %>% strsplit("model") %>%
  sapply(`[`, 1)
module_simi$m2 <- module_simi$module2 %>% strsplit("model") %>%
  sapply(`[`, 1)
module_simi$cross <- paste(module_simi$m1,module_simi$m2,sep = "_Vs_")

module_simi <- module_simi %>% 
  filter(module1 != "none")
head(module_simi)

# Aggregate standard names by group in module_otu
module_species <- module_otu %>%
  group_by(group) %>%
  summarise(species = paste(standard_name, collapse = "|"),
            otu = paste(ID, collapse = "|"))

# The biggest normal module that does not resemble any tumour module
saving_module <- data.frame(
  Species = strsplit(module_species$species[module_species$group == "Normalmodel_1"], "\\|")[[1]],
  OTU = strsplit(module_species$otu[module_species$group == "Normalmodel_1"], "\\|")[[1]]
) %>%
  merge(., node, by.x = "OTU", by.y = "ID")
abund_sp <- apply(mtx_cpm, MARGIN = 1, FUN = mean) / 1e+4
saving_module$abundance <- abund_sp[saving_module$Species]
saving_module <- na.omit(saving_module)

# Regression: Shannon index - candidate species
elnet_sp <- unique(saving_module$Species)
tumour_module_sp <- strsplit(module_species$species[module_species$group %in% 
                                                      paste0("Tumourmodel_", 1:5)], "\\|") %>%
  unlist(.) %>%
  unique(.)
elnet_sp <- elnet_sp[elnet_sp %in% tumour_module_sp]
x_raw <- as.matrix(log2(mtx_cpm[elnet_sp, names(shan_tumour)] + 1))
x_clr <- apply(x_raw, 2, function(v) clr(v)) %>% t(.)

# Elastic Net
elnet_fit <- cv.glmnet(
  x = x_clr,
  y = shan_tumour,
  alpha = 0.5,
  family = "gaussian",
  standardize = TRUE
)

coef_elnet <- coef(elnet_fit, s = "lambda.min") %>%
  as.matrix(.) %>%
  as.data.frame(.) %>%
  mutate(species = rownames(.)) %>%
  filter(species != "(Intercept)") %>%
  mutate(rank_contribute = rank(-lambda.min)) %>%
  filter(lambda.min > 1)

# Rank all parameters
# & Species %in% coef_elnet$species
saving_normal <- saving_module %>%
  filter(Group == "Normal" & Species %in% coef_elnet$species) %>%
  mutate(rank_degree = rank(-igraph.degree, ties.method = "min"),
         rank_closseness = rank(-igraph.closeness, ties.method = "min"))

saving_tumour <- saving_module %>%
  filter(Group == "Tumour"& Species %in% coef_elnet$species) %>%
  mutate(rank_degree = rank(-igraph.degree, ties.method = "min"),
         rank_closseness = rank(-igraph.closeness, ties.method = "min"))

saving_info <- merge(
  saving_normal[, c("OTU", "Species", "abundance",
                    "rank_degree", "rank_closseness")],
  saving_tumour[, c("OTU", "Species", 
                    "rank_degree", "rank_closseness")],
  by = c("OTU", "Species"),
  suffixes = c(".normal", ".tumour")
) %>%
  mutate(stability_degree = rank(abs(rank_degree.normal - rank_degree.tumour)),
         rank_degree = (rank(rank_degree.normal) + rank(rank_degree.tumour))/2,
         stability_closeness = rank(abs(rank_closseness.normal - rank_closseness.tumour)),
         rank_satbility = (rank(rank_closseness.normal) + rank(rank_closseness.tumour))/2)

saving_info <- saving_info %>%
  filter(Species %in% coef_elnet$species) %>%
  mutate(diversity_contribute = coef_elnet[Species, "rank_contribute"],
         rank_abundance = rank(-abundance)) %>%
  mutate(rank_score = sqrt(stability_degree * rank_degree) + 
           sqrt(stability_closeness * rank_satbility) + 
           diversity_contribute + rank_abundance,
         rank_total = rank(rank_score)) %>%
  arrange(rank_total)

rank_mtx <- saving_info[ , c("Species", "rank_abundance",
                             "stability_degree", "rank_degree",
                             "stability_closeness", "rank_satbility",
                             "diversity_contribute", "rank_total")] 
rownames(rank_mtx) <- saving_info$Species

long_ranks <- reshape2::melt(saving_info[ , c("Species", "rank_abundance",
                                             "stability_degree", "rank_degree",
                                             "stability_closeness", "rank_satbility",
                                             "diversity_contribute")])
long_total <- reshape2::melt(saving_info[ , c("Species", "rank_total")]) %>%
  arrange(desc(value)) %>%
  mutate(Species = factor(Species, level = Species))

long_ranks$Species <- factor(long_ranks$Species, level = long_total$Species)
colnames(long_ranks)[3] <- "Rank"
pdf(file.path(DIR_RES, "D_core_sp_selection_p1.pdf"), width = 4, height = 3)
ggplot(data = long_ranks, aes(x = variable, y = Species,col = Rank)) +
  geom_point(size = 2) +
  scale_color_distiller(palette = "Reds") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
colnames(long_total)[3] <- "Rank"
pdf(file.path(DIR_RES, "D_core_sp_selection_p2.pdf"), width = 5, height = 3)
ggplot(data = long_total, aes(y = rank(1/Rank), x = Species,
                              fill = Rank)) +
  geom_bar(stat = "identity") +
  scale_fill_distiller(palette = "Reds") +
  coord_flip() +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))
dev.off()

# Correlation between EF and shannon index
abund_ef <- mtx_cpm["Enterococcus_faecalis", names(shan_tumour)] / 1e+4
all(colnames(abund_ef) == names(shan_tumour))

cor_ef_shan <- psych::corr.test(t(abund_ef), shan_tumour)

df_ef_shan <- data.frame(
  row.names = names(shan_tumour),
  t(abund_ef),
  shannon = shan_tumour
)

pdf(file.path(DIR_RES, "E_cor_ef_shannon.pdf"), width = 3.3, height = 3.3)
ggplot(df_ef_shan, aes(x = Enterococcus_faecalis, y = shannon)) +
  geom_point(color = "tomato2", size = 1.5) +
  geom_smooth(method = lm, color = "#00A087",
              linewidth = 1.5, fill = "#4DBBD5") +
  annotate("text", x = 1.7, y = 2.8, color = 1,
           label = paste0("R = ", round(cor_ef_shan$r, 3),
                          "\np = ", format.pval(cor_ef_shan$p.adj, 3))) +
  xlab("E.faecalis abundance (%)") +
  ylab("Shannon index") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1),
        legend.title = element_blank(),
        legend.position = "bottom")
dev.off()
