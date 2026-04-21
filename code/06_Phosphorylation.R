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
DIR_FIG <- "/data/yzwang/project/AEG_seiri/results/F5_phospho/"

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
phospho   <- readRDS(file.path(DIR_RDS, "Phosphoproteome_human_preprocessed.rds"))
mtx_cpm   <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
target_sp <- read.csv(file.path(DIR_TAB, "Species_within_saving_module.csv"),
                      row.names = 1)

sp_list <- target_sp$Species
sp_list <- intersect(sp_list, rownames(mtx_cpm))

common_sample <- intersect(colnames(phospho), colnames(mtx_cpm))
common_tumour <- common_sample[grep("C", common_sample)]

phos_tumour <- phospho[, common_tumour] %>%
  filter(rowMeans(.) != min(phospho))
phos_tumour <- phos_tumour[!grepl("^-", rownames(phos_tumour)), ]

phos_detect <- phos_tumour %>%
  as.data.frame() %>%
  mutate(feature = rownames(.)) %>%
  pivot_longer(-feature, names_to = "sample", values_to = "intensity") %>%
  mutate(present = as.integer(intensity > 3.891872))

# Pre-filter features by prevalence (shared across all species)
feat_prev <- phos_detect %>%
  group_by(feature) %>%
  summarise(
    n_present = sum(present),
    n_absent  = sum(present == 0),
    .groups   = "drop"
  ) %>%
  mutate(prop_present = n_present / (n_present + n_absent)) %>%
  filter(prop_present > 0.1, prop_present < 0.9)

phos_detect_f <- phos_detect %>%
  filter(feature %in% feat_prev$feature)

#####################################
# Logistic regression for each species
run_logit_one_sp <- function(sp) {
  abund <- data.frame(
    sample = common_tumour,
    abundance = as.numeric(mtx_cpm[sp, common_tumour]) / 1e+4
  )
  
  phos_detect_f %>%
    left_join(abund, by = "sample") %>%
    group_by(feature) %>%
    nest() %>%
    mutate(
      fit = map(
        data,
        ~ tryCatch(
          glm(present ~ abundance, data = .x, family = binomial),
          error = function(e) NULL
        )
      )
    ) %>%
    filter(!map_lgl(fit, is.null)) %>%
    mutate(
      converged = map_lgl(fit, ~ .$converged),
      tidied    = map(fit, broom::tidy)
    ) %>%
    filter(converged) %>%
    unnest(tidied) %>%
    filter(term == "abundance") %>%
    ungroup() %>%
    mutate(
      OR      = exp(estimate),
      padj    = p.adjust(p.value, method = "BH"),
      species = sp
    ) %>%
    select(species, feature, estimate, std.error, OR, p.value, padj)
}

logistic_all <- map_dfr(sp_list, run_logit_one_sp)

sig_all <- logistic_all %>% filter(padj < 0.05)

write.csv(logistic_all, file.path(DIR_TAB, "Phospho_discrete_logit_allSp.csv"),
          row.names = FALSE)
write.csv(sig_all, file.path(DIR_TAB, "Phospho_discrete_significance_allSp.csv"),
          row.names = FALSE)