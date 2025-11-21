#######################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 05. After EF was highlighted, explore its effect
# in host transcription
#
#######################################################


setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/S2/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_TOOL <- "/data/yzwang/git_project/AEG_microbiome/utils/"

mtx_cpm <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
abund_ef <- mtx_cpm["Enterococcus_faecalis", names(shan_tumour)] / 1e+4
