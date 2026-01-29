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
library(ggplot2)

set.seed(42)

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
colnames(phospho) <- sub("^([A-Z]+)(\\d{6})$", "\\10\\2", colnames(phospho))
pair_file <- read.csv(file.path(DIR_TAB, "ms_sample_name_pairing.csv"))
pair_file$InnerSample <- sub("^([A-Z]+)(\\d{6})$", "\\10\\2", pair_file$InnerSample)

phospho <- phospho[ , c(TRUE, colnames(phospho)[-1] %in% pair_file$InnerSample)]
dim(phospho)

sample_map <- setNames(pair_file$ZipSample, pair_file$InnerSample)
colnames(phospho)[-1] <- sample_map[colnames(phospho)[-1]]

rm_dup <- c(which(phospho$Phospho_site == "TMPO_S184")[1],
            which(phospho$Phospho_site == "TMPO_S306")[1])
phospho <- phospho[-rm_dup, ]
dim(phospho)
rownames(phospho) <- phospho$Phospho_site
phospho <- phospho[ , -1]

saveRDS(phospho, file.path(DIR_RDS, "Phosphoproteome_human_preprocessed.rds"))

#################
# Formal analysis
phospho <- readRDS(file.path(DIR_RDS, "Phosphoproteome_human_preprocessed.rds"))

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
saveRDS(cor_ef_phos, file.path(DIR_RDS, "Correlation_ef_phosphosites.rds"))

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

# Logistic regression for detected/undetected signals
logistic_res <- phos_detect %>%
  left_join(abund_ef, by = "sample") %>%
  group_by(feature) %>%
  nest() %>%
  mutate(
    n_present = map_int(data, ~sum(.$present)),
    n_absent = map_int(data, ~sum(.$present == 0)),
    prop_present = n_present / (n_present + n_absent)
  ) %>%
  # Require at least 10 in each category AND between 10-90% detection
  filter(n_present >= 10, n_absent >= 10, 
         prop_present > 0.1, prop_present < 0.9) %>%
  mutate(
    fit = map(data, ~glm(present ~ Enterococcus_faecalis, data = ., 
                         family = binomial)),
    converged = map_lgl(fit, ~.$converged),
    tidied = map(fit, tidy)
  ) %>%
  filter(converged) %>%
  unnest(tidied) %>%
  filter(term == "Enterococcus_faecalis") %>%
  mutate(
    OR = exp(estimate),
    p_adj = p.adjust(p.value, method = "BH")
  ) %>%
  select(feature, estimate, OR, p.value, p_adj, n_present, n_absent)

write.csv(logistic_res, file.path(DIR_TAB, "Phos_discrete_sites.csv"))

# Spearman correaltion for detected intensities
cont_res <- phos_cont %>%
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

write.csv(cont_res, file.path(DIR_TAB, "Phos_continuous_sites.csv"))

merged_res <- logistic_res %>%
  rename(logistic_est = estimate,
         logistic_OR = OR,
         logistic_p = p.value,
         logistic_padj = p_adj) %>%
  full_join(
    cont_res %>%
      rename(spearman_rho = rho,
             cont_p = p_val,
             cont_padj = p_adj),
    by = "feature"
  )

write.csv(merged_res, file.path(DIR_TAB, "Phos_merged_res.csv"))

########################
# Significance selection
sig_discrete <- merged_res %>%
  filter(logistic_padj < 0.1)
dim(sig_discrete)

sig_continuous <- merged_res %>%
  filter(cont_padj < 0.1)
dim(sig_continuous)

sig_both <- merged_res %>%
  filter(logistic_padj < 0.1 & cont_padj < 0.1)
dim(sig_both)

plot_feature <- "CARMIL2_S1328"

ggplot(cont_res, aes(x = rho, y = -log10(p_adj))) +
  geom_point(aes(color = p_adj < 0.05 & abs(rho) > 0.3), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "blue") +
  scale_color_manual(values = c("grey60", "red3"), 
                     labels = c("Not significant", "Significant")) +
  labs(x = "Spearman correlation (ρ)", 
       y = "-log10(adjusted p-value)",
       color = NULL,
       title = "Phosphosite correlation with E. faecalis abundance") +
  theme_minimal()

ggplot(logistic_strict, aes(x = estimate, y = -log10(p_adj))) +
  geom_point(aes(color = p_adj < 0.05 & abs(estimate) > 0.5), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
  scale_color_manual(values = c("grey60", "blue3"), 
                     labels = c("Not significant", "Significant")) +
  labs(x = "Log odds ratio (coefficient)", 
       y = "-log10(adjusted p-value)",
       color = NULL,
       title = "Phosphosite detection vs E. faecalis abundance") +
  theme_minimal()
