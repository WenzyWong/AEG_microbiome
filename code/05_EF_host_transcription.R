#######################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 05. After EF was highlighted, explore its effect
# in host transcription
#
#######################################################
library(dplyr)
library(rtracklayer) # import gtf
library(DESeq2)
library(enrichplot)
library(paletteer)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/S2/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_TOOL <- "/data/yzwang/git_project/AEG_microbiome/utils/"

clinical <- readxl::read_excel(file.path(DIR_TAB, "AEG_clinical.xlsx"))

mtx_hcount <- readRDS(file.path(DIR_RDS, "hAEG_ExprCounts_Tumour.rds"))
rownames(mtx_hcount) <- mtx_hcount$Gene
mtx_hcount <- mtx_hcount[ , -1]
mtx_htpm <- read.delim(file.path(DIR_TAB, "hTumourTPM.txt"))

anno_coding <- import("/data/yzwang/reference/gtf/gencode_human_annotation.gtf") %>%
  as.data.frame(.) %>%
  filter(type == "gene" & gene_type == "protein_coding")
dim(anno_coding)

name_coding <- intersect(anno_coding$gene_name, rownames(mtx_hcount)) %>% 
  intersect(., rownames(mtx_htpm))
length(name_coding)

mtx_hcount <- mtx_hcount[name_coding, ]
mtx_htpm <- mtx_htpm[name_coding, colnames(mtx_hcount)]

mtx_cpm <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
abund_ef <- data.frame(
  sample = colnames(mtx_hcount),
  t(mtx_cpm["Enterococcus_faecalis", colnames(mtx_hcount)] / 1e+4)
)
qqnorm(abund_ef$Enterococcus_faecalis)
shapiro.test(abund_ef$Enterococcus_faecalis)
hist(abund_ef$Enterococcus_faecalis) # Normal distribution

summary(abund_ef$Enterococcus_faecalis)

# Differential analysis
# Using 1st Qu. and 3rd Qu. to define EF-high & EF-low groups
group_ef <- abund_ef %>%
  mutate(group = case_when(
    Enterococcus_faecalis < summary(Enterococcus_faecalis)[2] ~ "Low",
    Enterococcus_faecalis > summary(Enterococcus_faecalis)[5] ~ "High",
    TRUE ~ "Mid"
  ))

meta_ef <- group_ef$group
names(meta_ef) <- group_ef$sample

df_hcount <- cbind(gene.name = rownames(mtx_hcount), mtx_hcount)
dds <- DESeqDataSetFromMatrix(countData = df_hcount,
                              colData = group_ef,
                              design = ~ group,
                              tidy = TRUE)
dds <- DESeq(dds)

de_res <- results(dds, contrast = c("group", "High", "Low"))
de_res <- de_res[order(de_res$padj), ] %>% as.data.frame(.)

de_df <- data.frame(
  symbol = rownames(de_res),
  log2FC = de_res$log2FoldChange,
  p.adj = de_res$padj
) %>%
  na.omit(.) %>%
  mutate(change = case_when(
    log2FC > log2(1.5) & p.adj < 0.05 ~ "Up",
    log2FC < -log2(1.5) & p.adj < 0.05 ~ "Dn",
    TRUE ~ "NS"
  ))

draw_de <- de_df %>%
  mutate(log2FC = case_when(
    log2FC > 5 ~ 5,
    log2FC < -5 ~ 5,
    TRUE ~ log2FC
  ),
  p.adj = if_else(
    -log10(p.adj) > 5, 1e-5, p.adj
  ))

source(file.path(DIR_TOOL, "sort_gene.R"))
source(file.path(DIR_TOOL, "gsea_enrich.R"))
gsea_hallmark <- gsea_enrich("h.all", sort_gene(de_df), "Hs")
saveRDS(gsea_hallmark, file.path(DIR_RDS, "gsea_hallmark_ef_diff.rds"))
gsea_hallmark@result$ID[gsea_hallmark@result$NES < 0]
gsea_hallmark@result$ID[gsea_hallmark@result$NES > 0]

draw_hallmark <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                   "HALLMARK_INFLAMMATORY_RESPONSE",
                   "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                   "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
tolower(draw_hallmark)
pdf(file.path(DIR_RES, "E_gsea_curve.pdf"), width = 12, height = 7)
gseaplot2(gsea_hallmark, draw_hallmark, pvalue_table = T,
          color = paletteer_d("ggsci::default_jama")[1:4])
dev.off()

pdf(file.path(DIR_RES, "E_gsea_ridge.pdf"), width = 12, height = 7)
ridgeplot(gsea_hallmark)
dev.off()

# Matching genes in enriched terms
name_highlight <- strsplit(gsea_hallmark@result$core_enrichment
                           [gsea_hallmark@result$ID == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" |
                               gsea_hallmark@result$ID == "HALLMARK_IL6_JAK_STAT3_SIGNALING"|
                               gsea_hallmark@result$ID == "HALLMARK_IL2_STAT5_SIGNALING"],
                         "/") %>% unlist(.)
de_highlight <- de_df %>%
  filter(symbol %in% name_highlight & change != "NS") %>%
  arrange(p.adj)

ggplot(draw_de, aes(x = log2FC, y = -log10(p.adj), colour = change)) + 
  ggtitle("Tumour v.s. Normal") + 
  geom_point(size = 1, alpha = .8) + 
  xlim(-5, 5) +
  scale_color_manual(values = c("#126CAA", "grey", "#9A342C")) +
  geom_point(de_highlight, mapping = aes(x = log2FC, y = -log10(p.adj)),
            color = "#FFB900FF", size = 3) + 
  ggrepel::geom_label_repel(data = de_highlight[1:5, ],
                            aes(x = log2FC, y = -log10(p.adj), 
                                label = symbol),
                            color="grey27",
                            nudge_x = 1,
                            nudge_y = -1,
                            box.padding = 1.5,
                            alpha = .8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = log2(1.5), linetype = "dashed") +
  geom_vline(xintercept = -log2(1.5), linetype = "dashed") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))
