####################################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 06. Using phosphoproteomics to understand potential mechanism
#
####################################################################
library(dplyr)
library(stringr)
library(tidyr)
library(broom)
library(purrr)
library(circlize)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(KSEAapp)
library(ggridges)
library(ComplexHeatmap)
library(patchwork)

set.seed(42)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_FIG <- "/data/yzwang/project/AEG_seiri/results/F4/"
DIR_SUP <- "/data/yzwang/project/AEG_seiri/results/S4/"

##########################
# Preprocess original data
phospho <-
  read.delim(
    file.path(
      DIR_TAB,
      "Phosphoproteomics_iBAQ103_log2quantile_normlization_impute.txt"
    )
  )

id_match <-
  read.delim(file.path(DIR_TAB, "Phospho_IDs_Protein_Gene.txt"))

for (i in 1:nrow(phospho)) {
  pro_id <- strsplit(phospho$Phospho_site[i], "_")[[1]][1]
  gene_name <-
    id_match[id_match$Protein.accession == pro_id, "Gene.name"]
  site <- strsplit(phospho$Phospho_site[i], "_")[[1]][2]
  phospho$Phospho_site[i] <- paste0(gene_name, "_", site)
}
phospho <- phospho[-grep("^_", phospho$Phospho_site),]

colnames(phospho) <- sub("^([A-Z]+)(\\d{6})$", "\\10\\2", colnames(phospho))
pair_file <- read.csv(file.path(DIR_TAB, "ms_sample_name_pairing.csv"))
pair_file$InnerSample <- sub("^([A-Z]+)(\\d{6})$", "\\10\\2", pair_file$InnerSample)

write.csv(pair_file[ , c("ZipSample", "InnerSample", "Match")],
          file.path(DIR_TAB, "Phosphoproteome_sample_name_match.csv"))

phospho <- phospho[, c(TRUE, colnames(phospho)[-1] %in% pair_file$InnerSample)]
sample_map <- setNames(pair_file$ZipSample, pair_file$InnerSample)
colnames(phospho)[-1] <- sample_map[colnames(phospho)[-1]]

rm_dup <- c(
  which(phospho$Phospho_site == "TMPO_S184")[1],
  which(phospho$Phospho_site == "TMPO_S306")[1]
)
phospho <- phospho[-rm_dup,]
rownames(phospho) <- phospho$Phospho_site
phospho <- phospho[,-1]

saveRDS(phospho, file.path(DIR_RDS, "Phosphoproteome_human_preprocessed.rds"))

#################
# Formal analysis
phospho <- readRDS(file.path(DIR_RDS, "Phosphoproteome_human_preprocessed.rds"))

mtx_cpm <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
common_sample <- intersect(colnames(phospho), colnames(mtx_cpm))
common_tumour <- common_sample[grep("C", common_sample)]

phos_tumour <- phospho[, common_tumour] %>%
  filter(rowMeans(.) != min(phospho))
phos_tumour <- phos_tumour[!grepl("^-", rownames(phos_tumour)),]

abund_ef <- data.frame(sample = common_tumour,
                       t(mtx_cpm["Enterococcus_faecalis", common_tumour] / 1e+4))

phos_detect <- phos_tumour %>%
  as.data.frame() %>%
  mutate(feature = rownames(.)) %>%
  pivot_longer(-feature, names_to = "sample", values_to = "intensity") %>%
  mutate(present = as.integer(intensity > 3.891872))

#####################
# Logistic regression
logistic_ori <- phos_detect %>%
  left_join(abund_ef, by = "sample") %>%
  group_by(feature) %>%
  nest() %>%
  mutate(
    n_present = map_int(data, ~ sum(.$present)),
    n_absent  = map_int(data, ~ sum(.$present == 0)),
    prop_present = n_present / (n_present + n_absent)
  ) %>%
  filter(prop_present > 0.1, prop_present < 0.9) %>%
  mutate(
    fit = map(
      data,
      ~ glm(
        present ~ Enterococcus_faecalis,
        data = .,
        family = binomial
      )
    ),
    converged = map_lgl(fit, ~ .$converged),
    tidied = map(fit, tidy)
  ) %>%
  filter(converged) %>%
  unnest(tidied) %>%
  filter(term == "Enterococcus_faecalis") %>%
  mutate(OR = exp(estimate)) %>%
  select(feature, estimate, std.error, OR, p.value, n_present, n_absent)

sig_discrete <- logistic_ori %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  filter(padj < 0.05)

write.csv(sig_discrete, file.path(DIR_TAB, "Phospho_discrete_significance.csv"))
