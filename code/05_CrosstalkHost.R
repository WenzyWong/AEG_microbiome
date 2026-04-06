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
library(ggplot2)
library(DESeq2)
library(enrichplot)
library(paletteer)

DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/F3_crosstalk/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_TOOL <- "/data/yzwang/git_project/AEG_microbiome/utils/"

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
