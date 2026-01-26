#######################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 06. Explore drug-related effects of EF
#
#######################################################
library(ComplexHeatmap)
library(psych)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_FIG <- "/data/yzwang/project/AEG_seiri/results/F3/"
DIR_SUP <- "/data/yzwang/project/AEG_seiri/results/S3/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_TOOL <- "/data/yzwang/git_project/AEG_microbiome/utils/"

###########
# IDWAS
# Load data
idwas_drug <- readxl::read_excel(file.path(DIR_TAB, "IDWAS_Supplemental_Table_S2.xlsx"), 
                                 skip = 1) %>%
  select(Drug)

idwas_mtx <- readRDS(file.path(DIR_RDS, "drugPred_Mtx.rds"))
idwas_tumour_mtx <- idwas_mtx[ , grep("C", colnames(idwas_mtx))]

mtx_cpm <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
abund_ef <- data.frame(
  sample = colnames(mtx_cpm),
  t(mtx_cpm["Enterococcus_faecalis", ] / 1e+4)
)

abund <- sort(apply(mtx_cpm, MARGIN = 1, FUN = mean) / 1e+4, decreasing = T)
abundant_sp <- names(abund)[abund > 0.8]

abund_mtx <- mtx_cpm[abundant_sp, ] / 1e+4
dim(abund_mtx)

dim(idwas_mtx)

gdsc_target <- read.csv(file.path(DIR_TAB, "Drug_list2023.csv"))

# Plot predicted drug responses: idwas
draw_idwas_mtx <- scale(idwas_mtx, center = T)
col_zs <- circlize::colorRamp2(c(-10, 0, 10), c("#2166AC", "white", "#B2182B"))
pdf(file.path(DIR_SUP, "A_heatmap_idwas_samples.pdf"), width = 6, height = 10)
Heatmap(draw_idwas_mtx, show_column_names = F, 
        name = "Scaled response", column_title = "Samples",
        column_title_side = "bottom",
        col = col_zs,
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

# EF
cor_ef_idwas <- corr.test(t(as.matrix(idwas_mtx)), 
                          t(as.matrix(abund_mtx["Enterococcus_faecalis", ])),
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

#########
# CARE
library(rtracklayer) # For importing gtf file
library(tidyverse)
library(survival)
library(pROC)

anno_dt <- import("/data/yzwang/reference/gencode_ref/gencode_human_annotation.gtf") %>%
  as.data.frame()
coding_gene_list <- anno_dt$gene_name[anno_dt$gene_type == "protein_coding"] %>% unique(.)
length(coding_gene_list)
rm(anno_dt)

# Preprocessing human expression matrix
hexp_tpm <- readRDS(paste0(DIR_RDS, "AEG_humanTPM_Symbol.rds"))
hexp_tpm <- hexp_tpm[rownames(hexp_tpm) %in% coding_gene_list, ]
hexp_tpm <- hexp_tpm[rowMeans(hexp_tpm) > 1, ]
dim(hexp_tpm)

hexp_tpm <- hexp_tpm[common_genes, grep("C", colnames(hexp_tpm))]
hexp_norm <- t(apply(hexp_tpm, 1, function(gene_values) {
  (gene_values - mean(gene_values, na.rm = TRUE)) / sd(gene_values, na.rm = TRUE)
}))
gene_vars <- apply(hexp_norm, 1, var, na.rm = T)
top_var_genes <- names(gene_vars)[gene_vars > median(gene_vars, na.rm = TRUE)]
hexp_filt <- hexp_norm[top_var_genes, ]
dim(hexp_filt)

# Importing CARE scores
care_score <- read.table(file.path(DIR_TAB, "CARE_GDSC"), header = T, fill = T) %>%
  filter(Type == "t", )

meta_cols <- c("Response", "Type", "Target")
gene_cols <- setdiff(colnames(care_score), meta_cols)

care_score[, gene_cols] <- lapply(care_score[, gene_cols], function(x) {
  as.numeric(as.character(x))
})

care_aggregated <- care_score %>%
  group_by(Response) %>%
  summarise(
    Type = Type[1],
    Target = paste(Target, collapse = "/"),
    across(all_of(gene_cols), ~mean(.x, na.rm = TRUE)),
    .groups = 'drop'
  )

common_genes <- intersect(rownames(hexp_filt), colnames(care_score))
length(common_genes)

care_mtx <- care_aggregated[ , common_genes]
rownames(care_mtx) <- care_aggregated$Response

cor_care <- cor(hexp_filt, t(care_mtx))

Heatmap(cor_care)

######################
# All abundant species
all(colnames(idwas_mtx) == colnames(abund_mtx))
cor_sp_idwas <- corr.test(t(as.matrix(idwas_mtx)), t(as.matrix(abund_mtx)),
                         method = "spearman")
cor_r <- cor_sp_idwas$r
cor_p <- cor_sp_idwas$p.adj

dim(cor_r)

rownames(cor_r) <- gsub("\\.", "-", rownames(cor_r))
pathway_all <- gdsc_target[gdsc_target$Name %in% rownames(cor_r), 
                           c("Name", "Target.pathway")]
pathway_all <- pathway_all[!duplicated(pathway_all$Name), ]
dim(pathway_all)
potent_syn <- setdiff(rownames(cor_r), pathway_all$Name)

notfound <- c()
for (i in 1:length(potent_syn)) {
  tmp_name <- gdsc_target$Name[grep(potent_syn[i], gdsc_target$Synonyms)][1]
  tmp_path <- gdsc_target$Target.pathway[grep(potent_syn[i], gdsc_target$Synonyms)][1]
  if (is.na(tmp_name)) {
    notfound <- c(notfound, potent_syn[i])
    pathway_all <- rbind(pathway_all, c(potent_syn[i], "Unannotated"))
  } else {
    rownames(cor_r)[rownames(cor_r) == potent_syn[i]] <- tmp_name
    rownames(cor_p)[rownames(cor_p) == potent_syn[i]] <- tmp_name
    pathway_all <- rbind(pathway_all, c(tmp_name, tmp_path))
  }
}

rownames(pathway_all) <- pathway_all$Name
pathway_all <- pathway_all[rownames(cor_r), ]
pathway_all <- pathway_all[ , -1]
names(pathway_all) <- rownames(cor_r)
length(pathway_all)
length(unique(pathway_all))

col_pathway_all <- c(paletteer_d("pals::kelly"))
names(col_pathway_all) <- unique(pathway_all)

right_anno_all <- rowAnnotation(
  Pathway = pathway_all,
  col = list(Pathway = col_pathway_all)
)

dim(cor_r)

neg_all <- c()
pos_all <- c()
for (i in 1:ncol(cor_r)) {
  neg_tmp <- length(rownames(cor_r)[cor_r[, i] < -0.2 & cor_p[, i] < 0.05])
  neg_all <- c(neg_all, neg_tmp)
  pos_tmp <- length(rownames(cor_r)[cor_r[, i] > 0.2 & cor_p[, i] < 0.05])
  pos_all <- c(pos_all, pos_tmp)
}
names(neg_all) <- colnames(cor_r)
names(pos_all) <- colnames(cor_r)

top_anno_all <- HeatmapAnnotation(
  log2CPM = apply(log2(mtx_cpm[abundant_sp, ] + 1), MARGIN = 1, FUN = mean),
  Positive = anno_barplot(pos_all, gp = gpar(border = NA,fill = "#701145FF",lty = "blank")),
  Negative = anno_barplot(neg_all, gp = gpar(border = NA,fill = "#008280FF",lty = "blank")),
  col = list(log2CPM = circlize::colorRamp2(c(12, 17), 
                                            c("white", "darkgreen")))
)

cn <- colnames(cor_r)
Heatmap(cor_r, name = "Rs", col = rev(brewer.pal(n = 11, name = "RdBu")),
        row_names_gp = gpar(fontsize = 4), 
        show_column_names = T, right_annotation = right_anno_all,
        top_annotation = top_anno_all,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(cor_p[i, j] < 0.001) {
            gb = textGrob("***")
            gb_h = convertHeight(grobHeight(gb), "mm")
            grid.text("***", x, y - gb_h * 0.2,
                      gp = gpar(col = 1, cex = .5))
          } else if(cor_p[i, j] < 0.01) {
            gb = textGrob("**")
            gb_h = convertHeight(grobHeight(gb), "mm")
            grid.text("**", x, y - gb_h * 0.2,
                      gp = gpar(col = 1, cex = .5))
          } else if(cor_p[i, j] < 0.05) {
            gb = textGrob("*")
            gb_h = convertHeight(grobHeight(gb), "mm")
            grid.text("*", x, y - gb_h * 0.2,
                      gp = gpar(col = 1, cex = .5))
          }
        }
) 

cor_r[,"Enterococcus_faecalis"]
