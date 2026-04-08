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
source(file.path(DIR_TOOL, "run_deseq_by_sp.R"))
source(file.path(DIR_TOOL, "run_gsea_by_species.R"))
source(file.path(DIR_TOOL, "extract_ranking_df.R"))
source(file.path(DIR_TOOL, "draw_summary_ridgeplot.R"))
source(file.path(DIR_TOOL, "draw_gsea_dotplot.R"))

target_sp <- read.csv(file.path(DIR_TAB, "Species_within_saving_module.csv"),
                      row.names = 1)

sample_groups <- list(
  tumour = colnames(mtx_hcount)[startsWith(colnames(mtx_hcount), "C")],
  normal = colnames(mtx_hcount)[!startsWith(colnames(mtx_hcount), "C")]
)

for (grp in names(sample_groups)) {
  samples <- sample_groups[[grp]]
  dir_rds_grp <- file.path(DIR_RDS, grp)
  dir_res_grp <- file.path(DIR_RES, grp)
  dir.create(dir_rds_grp, showWarnings = FALSE)
  dir.create(dir_res_grp, showWarnings = FALSE)
  
  # DE analysis
  abund_target <- data.frame(
    sample = samples,
    t(mtx_cpm[target_sp$Species, samples] / 1e+4)
  )
  
  all_de_results <- lapply(
    target_sp$Species,
    run_deseq_by_sp,
    abund_df     = abund_target,
    count_mtx    = mtx_hcount[, samples],
    target_sp_df = target_sp
  )
  names(all_de_results) <- target_sp$Species
  all_de_results <- Filter(Negate(is.null), all_de_results)
  
  # GSEA analysis
  all_gsea_results <- lapply(
    names(all_de_results),
    run_gsea_by_species,
    de_results = all_de_results,
    dir_rds    = dir_rds_grp
  )
  names(all_gsea_results) <- names(all_de_results)
  
  # Plots
  for (gene_set in c("hallmark", "kegg")) {
    pdf(file.path(dir_res_grp, paste0("summary_gsea_ridge_", gene_set, ".pdf")),
        width = 12, height = length(all_gsea_results) * 0.8)
    print(draw_summary_ridgeplot(all_gsea_results, gene_set))
    dev.off()
    
    pdf(file.path(dir_res_grp, paste0("summary_gsea_dot_", gene_set, ".pdf")),
        width = 8, height = length(all_gsea_results))
    print(draw_gsea_dotplot(all_gsea_results, gene_set))
    dev.off()
  }
  
  # Save results
  saveRDS(all_de_results,   file.path(dir_rds_grp, "all_de_results.rds"))
  saveRDS(all_gsea_results, file.path(dir_rds_grp, "all_gsea_results.rds"))
}

gsea_tumour <- readRDS(file.path(DIR_RDS, "tumour", "all_gsea_results.rds"))
gsea_normal <- readRDS(file.path(DIR_RDS, "normal", "all_gsea_results.rds"))

n_sp    <- max(length(gsea_tumour), length(gsea_normal))
n_paths <- length(unique(collect_gsea_df(
  c(gsea_tumour, gsea_normal), "hallmark")$ID))

for (gene_set in c("hallmark", "kegg")) {
  n_paths <- length(unique(collect_gsea_df(
    c(gsea_tumour, gsea_normal), gene_set)$ID))
  
  pdf(file.path(DIR_RES, paste0("compare_gsea_dot_", gene_set, ".pdf")),
      width  = n_sp * 0.6 + 4,
      height = n_paths * 0.3)
  print(draw_gsea_dotplot(gsea_tumour, gsea_normal, gene_set = gene_set))
  dev.off()
}

##################
# Human proteomics
mtx_hpro <- readxl::read_excel(file.path(DIR_TAB, "SupData13ProteinIntensity.xlsx"))

sample_pre <- colnames(mtx_hpro)[3:ncol(mtx_hpro)]
sample_pre <- gsub("AEG", "", sample_pre) %>%
  gsub("_T$", "", .) %>%
  gsub("_N$", "", .) %>%
  { ifelse(grepl("_T$", colnames(mtx_hpro)[3:ncol(mtx_hpro)]), paste0("C", .), paste0("N", .)) }

colnames(mtx_hpro)[3:ncol(mtx_hpro)] <- sample_pre
