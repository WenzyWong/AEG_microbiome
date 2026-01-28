####################################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 07. Using phosphoproteomics to understand potential mechanism
#
####################################################################
library(dplyr)
library(stringr)
library(tidyr)
library(broom)
library(purrr)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_FIG <- "/data/yzwang/project/AEG_seiri/results/F4/"

##########################
# Preprocess original data
phospho <- read.delim(file.path(DIR_TAB, "Phosphoproteomics_iBAQ103_log2quantile_normlization_impute.txt"))
dim(phospho)
id_match <- read.delim(file.path(DIR_TAB, "Phospho_IDs_Protein_Gene.txt"))
dim(id_match)

# Match IDs
for (i in 1:nrow(phospho)) {
  pro_id <- strsplit(phospho$Phospho_site[i], "_")[[1]][1]
  gene_name <- id_match[id_match$Protein.accession == pro_id, "Gene.name"]
  site <- strsplit(phospho$Phospho_site[i], "_")[[1]][2]
  phospho$Phospho_site[i] <- paste0(gene_name, "_", site)
}
phospho <- phospho[-grep("^_", phospho$Phospho_site), ]
dim(phospho)

# Pair file names
pair_file <- readRDS(paste0(DIR_RDS, "Proteome_file_name_pairing.rds"))

for (i in 2:ncol(phospho)) {
  tmp_sample <- colnames(phospho)[i]
  if (str_length(tmp_sample) < 8) {
    if (grepl("C", tmp_sample)) {
      tmp_sample <- paste0("C0", gsub("[C|N]", "", tmp_sample))
    } else {
      tmp_sample <- paste0("N0", gsub("[C|N]", "", tmp_sample))
    }
  }
  colnames(phospho)[i] <- tmp_sample
}

for (i in 1:nrow(pair_file)) {
  tmp_sample <- pair_file$OriName[i]
  if (str_length(tmp_sample) < 8) {
    if (grepl("C", tmp_sample)) {
      tmp_sample <- paste0("C0", gsub("[C|N]", "", tmp_sample))
    } else {
      tmp_sample <- paste0("N0", gsub("[C|N]", "", tmp_sample))
    }
  }
  pair_file$OriName[i] <- tmp_sample
}

phospho <- phospho[ , c(TRUE, colnames(phospho)[-1] %in% pair_file$OriName)]
dim(phospho)

rownames(pair_file) <- pair_file$OriName
colnames(phospho)[-1] <- pair_file[colnames(phospho)[-1], "TargetName"]
rm_dup <- c(which(phospho$Phospho_site == "TMPO_S184")[1],
            which(phospho$Phospho_site == "TMPO_S306")[1])
phospho <- phospho[-rm_dup, ]
dim(phospho)
rownames(phospho) <- phospho$Phospho_site
phospho <- phospho[ , -1]

saveRDS(phospho, file.path(DIR_RDS, "Phosphoproteome_human_preprocessed.rds"))


phospho <- file.path(DIR_RDS, "Phosphoproteome_human_preprocessed.rds")

mtx_cpm <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
common_sample <- intersect(colnames(phospho), colnames(mtx_cpm))
common_tumour <- common_sample[grep("C", common_sample)]
length(common_tumour)

phos_tumour <- phospho[ , common_tumour] %>%
  filter(rowMeans(.) != min(phospho))
phos_tumour <- phos_tumour[!grepl("^-", rownames(phos_tumour)), ]
dim(phos_tumour)

abund_ef <- data.frame(
  sample = common_tumour,
  t(mtx_cpm["Enterococcus_faecalis", common_tumour] / 1e+4)
)

cor_ef_phos <- psych::corr.test(abund_ef$Enterococcus_faecalis, t(as.matrix(phos_tumour)))

cor_res <- data.frame(
  point = colnames(cor_ef_phos$r),
  rs = as.numeric(cor_ef_phos$r),
  padj = as.numeric(cor_ef_phos$p.adj)
)

phos_detect <- phos_tumour %>%
  as.data.frame() %>%
  mutate(feature = rownames(.)) %>%
  pivot_longer(-feature, names_to = "sample", values_to = "intensity") %>%
  mutate(present = as.integer(intensity > 3.891872))

phos_cont <- phos_tumour %>%
  as.data.frame() %>%
  mutate(feature = rownames(.)) %>%
  pivot_longer(-feature, names_to = "sample", values_to = "intensity") %>%
  filter(intensity > 3.891872) %>%
  mutate(log_intensity = log2(intensity))

logistic_results <- phos_detect %>%
  left_join(abund_ef, by = "sample") %>%
  group_by(feature) %>%
  nest() %>%
  mutate(
    n_present = map_int(data, ~sum(.$present)),
    n_absent = map_int(data, ~sum(.$present == 0)),
    fit = map(data, ~{
      if(sum(.$present) < 2 || sum(.$present == 0) < 2) {
        return(NULL)
      }
      tryCatch(
        glm(present ~ Enterococcus_faecalis, data = ., family = binomial, 
            control = glm.control(maxit = 50)),
        warning = function(w) NULL,
        error = function(e) NULL
      )
    }),
    converged = map_lgl(fit, ~!is.null(.) && .$converged),
    tidied = map(fit, ~if(!is.null(.)) tidy(.) else tibble())
  ) %>%
  filter(converged) %>%
  unnest(tidied) %>%
  filter(term == "Enterococcus_faecalis") %>%
  mutate(
    OR = exp(estimate),
    p_adj = p.adjust(p.value, method = "BH")
  ) %>%
  select(feature, estimate, OR, p.value, p_adj)

cont_results <- phos_cont %>%
  left_join(abund_ef, by = "sample") %>%
  group_by(feature) %>%
  summarise(
    n = sum(!is.na(Enterococcus_faecalis) & !is.na(log_intensity)),
    rho = if(n > 2) cor(Enterococcus_faecalis, log_intensity, 
                        method = "spearman", use = "complete.obs") else NA_real_,
    p_val = if(n > 2) cor.test(Enterococcus_faecalis, log_intensity, 
                               method = "spearman", exact = FALSE)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_val, method = "BH"))

group_ef <- abund_ef %>%
  mutate(group = case_when(
    Enterococcus_faecalis < summary(Enterococcus_faecalis)[2] ~ "Low",
    Enterococcus_faecalis > summary(Enterococcus_faecalis)[5] ~ "High",
    TRUE ~ "Mid"
  ))


