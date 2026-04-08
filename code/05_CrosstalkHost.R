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
library(DESeq2)
library(enrichplot)
library(paletteer)
library(patchwork)

DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/F3_crosstalk/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_TOOL <- "/data/yzwang/git_project/AEG_microbiome/utils/"

# Import and preprocess
clinical <- readxl::read_excel(file.path(DIR_TAB, "AEG_clinical.xlsx"))
dim(clinical)

mtx_hcount <- read.csv(file.path(DIR_TAB, "gene_count_matrix.csv"), row.names = 1) %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_name = gsub(".*\\|", "", gene_id)) %>%
  select(-gene_id) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum)) %>%
  filter(!grepl("^<class", gene_name)) %>%
  column_to_rownames("gene_name")
rownames(mtx_hcount)[1:3]
dim(mtx_hcount)
mtx_htpm <- read.csv(file.path(DIR_TAB, "gene_tpm_matrix.csv"), row.names = 1) %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_name = gsub(".*\\|", "", gene_id)) %>%
  select(-gene_id) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum)) %>%
  filter(!grepl("^<class", gene_name)) %>%
  column_to_rownames("gene_name")
rownames(mtx_htpm)[1:3]
dim(mtx_htpm)

mtx_htpm_filter <- mtx_htpm %>%
  filter(rowMeans(.) > 1)
dim(mtx_htpm_filter)

mtx_hcount_filter <- mtx_hcount[rownames(mtx_htpm_filter), ]
dim(mtx_hcount_filter)

saveRDS(mtx_hcount_filter, file.path(DIR_RDS, "hTumourCNT_filtered.rds"))
saveRDS(mtx_htpm_filter, file.path(DIR_RDS, "hTumourTPM_filtered.rds"))

mtx_cpm <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
dim(mtx_cpm)

mtx_hcount_filter <- readRDS(file.path(DIR_RDS, "hTumourCNT_filtered.rds"))
mtx_htpm_filter <- readRDS(file.path(DIR_RDS, "hTumourTPM_filtered.rds"))

#######################
# Differential analysis
tumour_samples <- colnames(mtx_hcount)[startsWith(colnames(mtx_hcount), "C")]

source(file.path(DIR_TOOL, "run_deseq_by_sp.R"))
target_sp <- read.csv(file.path(DIR_TAB, "Species_within_saving_module.csv"),
                      row.names = 1)


abund_target <- data.frame(
  sample = tumour_samples,
  t(mtx_cpm[target_sp$Species, tumour_samples] / 1e+4)
)

all_de_results <- lapply(
  target_sp$Species,
  run_deseq_by_sp,
  abund_df     = abund_target,
  count_mtx    = mtx_hcount[, tumour_samples],
  target_sp_df = target_sp
)

names(all_de_results) <- target_sp$Species

all_de_results <- Filter(Negate(is.null), all_de_results)

###############
# GSEA analysis
source(file.path(DIR_TOOL, "run_gsea_by_species.R"))
source(file.path(DIR_TOOL, "extract_ranking_df.R"))
source(file.path(DIR_TOOL, "draw_summary_ridgeplot.R"))

all_gsea_results <- lapply(
  names(all_de_results),
  run_gsea_by_species,
  de_results = all_de_results,
  dir_rds    = DIR_RDS
)
names(all_gsea_results) <- names(all_de_results)

for (gene_set in c("hallmark", "kegg")) {
  pdf(file.path(DIR_RES, paste0("summary_gsea_ridge_", gene_set, ".pdf")),
      width = 12, height = length(all_gsea_results) * 0.8)
  print(draw_summary_ridgeplot(all_gsea_results, gene_set))
  dev.off()
}

# Dot plot
source(file.path(DIR_TOOL, "draw_gsea_dotplot.R"))
for (gene_set in c("hallmark", "kegg")) {
  pdf(file.path(DIR_RES, paste0("summary_gsea_dot_", gene_set, ".pdf")),
      width = 8, height = length(all_gsea_results))
  print(draw_gsea_dotplot(all_gsea_results, gene_set))
  dev.off()
}

