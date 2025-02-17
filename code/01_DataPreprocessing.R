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
# Using 16S genus-level abundance matrix to generate 16S-based abundance matrix
abund_gc_16S <- as.data.frame(read_excel(paste0(DIR_TABLE, "GC_16S_genus_abundance.xlsx")))
rownames(abund_gc_16S) <- gsub("g__", "", abund_gc_16S$Genus)
abund_gc_16S <- abund_gc_16S[ , -1]
abund_gc_16S <- abund_gc_16S[!rownames(abund_gc_16S) %in% contaminants, ]

# Using RNA-seq genus-level count matrix to generate RNA-based abundance matrix
count_gc_RNA <- read.csv(paste0(DIR_TABLE, "GC_RNA_genus_count.csv"))

tmp_split <- unlist(strsplit(count_gc_RNA$X, "\\|"))
rownames(count_gc_RNA) <- gsub("g__", "", tmp_split[grep("g__", tmp_split)])

count_gc_RNA <- count_gc_RNA[ , -1]
count_gc_RNA[is.na(count_gc_RNA)] <- 0
count_gc_RNA <- count_gc_RNA[!rownames(count_gc_RNA) %in% contaminants, ]
saveRDS(count_gc_RNA, paste0(DIR_RDS, "GC_RNA_count.rds")) # Save for calculating Shannon indexes

abund_gc_RNA <- count_gc_RNA
for (i in 1:ncol(abund_gc_RNA)) {
  abund_gc_RNA[ , i] <- (abund_gc_RNA[ , i] * 1e2) / colSums(count_gc_RNA)[i]
}
saveRDS(abund_gc_RNA, paste0(DIR_RDS, "GC_RNA_abundance.rds"))
