#######################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 03. Compare AEG microbiome with other GI cancers,
# including GC, AEG, ESCC
#
#######################################################
library(readxl)
library(ggvenn)
library(ggalluvial)
library(ggplot2)
library(vegan) # diversity, vegist
library(ape) # pcoa
library(paletteer)
library(ComplexHeatmap)
library(RColorBrewer)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TABLE <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_TOOL <- "/data/yzwang/functions/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/F1/"

##############
# Data loading
cpm_escc <- readRDS(paste0(DIR_RDS, "gESCC_cpm.rds"))
cpm_escc <- cpm_escc[ , grep("GSE235537", colnames(cpm_escc))]
cpm_gc <- readRDS(paste0(DIR_RDS, "gGC_CPM_RNA_WithoutCompFilt.rds"))
cpm_aeg <- readRDS(paste0(DIR_RDS, "gAEG_CPM_RNA_WithoutCompFilt.rds"))

count_escc <- readRDS(paste0(DIR_RDS, "gESCC_count.rds"))
count_escc <- count_escc[ , grep("GSE235537", colnames(count_escc))]
count_gc <- readRDS(paste0(DIR_RDS, "gGC_count_RNA.rds"))
count_aeg <- readRDS(paste0(DIR_RDS, "gAEG_count.rds"))

ESCC_selected_data <- as.data.frame(read_excel(file.path(DIR_TABLE, "ESCC_selected_data.xlsx")))
rownames(ESCC_selected_data) <- ESCC_selected_data$Run
colnames(cpm_escc) <- gsub("GSE235537.", "", colnames(cpm_escc))
colnames(cpm_escc) <- paste0(ESCC_selected_data[colnames(cpm_escc), "Type"],
                             ESCC_selected_data[colnames(cpm_escc), "Patient"])

