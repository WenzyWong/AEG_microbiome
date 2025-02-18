#####################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 01. Preprocess sequencing data
#
#####################################################

library(readxl)
library(stringr)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TABLE <- "/data/yzwang/project/AEG_seiri/table_infos/"

contaminants <- colnames(read.table(paste0(DIR_TABLE, "contaminants.txt"), header = T))
### GC ###
# Read and process 16S abundance data
abund_16s <- read_excel(file.path(DIR_TABLE, "GC_16S_genus_abundance.xlsx")) |>
  as.data.frame()
rownames(abund_16s) <- gsub("g__", "", abund_16s$Genus)
abund_16s <- abund_16s[, -1]
abund_16s <- abund_16s[!rownames(abund_16s) %in% contaminants, ]

# Calculate 16S CPM
cpm_16s <- abund_16s * 1e4
saveRDS(cpm_16s, file.path(DIR_RDS, "gGC_CPM_16S_WithoutCompFilt.rds"))

# Read and process RNA-seq count data
count_rna <- read.csv(file.path(DIR_TABLE, "GC_RNA_genus_count.csv"))
rownames(count_rna) <- gsub("g__", "", unlist(strsplit(count_rna$X, "\\|"))[grep("g__", unlist(strsplit(count_rna$X, "\\|")))])
count_rna <- count_rna[, -1]
count_rna[is.na(count_rna)] <- 0
count_rna <- count_rna[!rownames(count_rna) %in% contaminants, ]
saveRDS(count_rna, file.path(DIR_RDS, "gGC_count_RNA.rds"))

# Calculate RNA-seq CPM
cpm_rna <- t(t(count_rna) * 1e6 / colSums(count_rna))
saveRDS(cpm_rna, file.path(DIR_RDS, "gGC_CPM_RNA_WithoutCompFilt.rds"))

# Create comparable matrices
common_rows <- sort(intersect(rownames(cpm_16s), rownames(cpm_rna)))
common_cols <- sort(intersect(colnames(cpm_16s), colnames(cpm_rna)))

comp_16s <- cpm_16s[common_rows, common_cols]
comp_rna <- cpm_rna[common_rows, common_cols]

dim(comp_16s)
dim(comp_rna)

saveRDS(comp_16s, file.path(DIR_RDS, "gGC_CPM_16S_comparison.rds"))
saveRDS(comp_rna, file.path(DIR_RDS, "gGC_CPM_RNA_comparison.rds"))
