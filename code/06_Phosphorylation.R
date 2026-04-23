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
library(uwot)       # UMAP
library(Rtsne)
library(paletteer)

set.seed(42)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_FIG <- "/data/yzwang/project/AEG_seiri/results/F5_phospho/"

##########################
# Preprocess original data
phospho <- read.delim(file.path(DIR_TAB,
                                "Phosphoproteomics_iBAQ103_log2quantile_normlization_impute.txt"))

id_match <- read.delim(file.path(DIR_TAB, "Phospho_IDs_Protein_Gene.txt"))

for (i in 1:nrow(phospho)) {
  pro_id <- strsplit(phospho$Phospho_site[i], "_")[[1]][1]
  gene_name <- id_match[id_match$Protein.accession == pro_id, "Gene.name"]
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

write.csv(logistic_all, file.path(DIR_TAB, "Phospho_discrete_logit_allSp.csv"),
          row.names = FALSE)

# Clustering based on sp-site association
mat_df <- logistic_all %>%
  mutate(
    log2OR = estimate / log(2),
    log2OR_shrunk = ifelse(p.value < 0.05, log2OR, 0)
  )

mat <- mat_df %>%
  select(feature, species, log2OR_shrunk) %>%
  pivot_wider(names_from = species, values_from = log2OR_shrunk) %>%
  tibble::column_to_rownames("feature") %>%
  as.matrix()

mat[is.na(mat)] <- 0

# Keep only features with at least one significant hit
keep_feat <- rowSums(mat != 0) >= 1
mat_f     <- mat[keep_feat, , drop = FALSE]
dim(mat_f)

# Cluster based on heatmap
ann_col <- data.frame(
  n_sig_feat = colSums(mat_f != 0),
  row.names  = colnames(mat_f)
)

lim <- quantile(abs(mat_f[mat_f != 0]), 0.98)
brk <- seq(-lim, lim, length.out = 101)

# Delegate clusters
hc_row <- hclust(as.dist(1 - cor(t(mat_f))), method = "ward.D2")
k_site <- 6
site_cluster <- cutree(hc_row, k = k_site)

hc_col <- hclust(as.dist(1 - cor(mat_f)), method = "ward.D2")
k_sp   <- 4
sp_cluster <- cutree(hc_col, k = k_sp)

table(site_cluster)
table(sp_cluster)

# UMAP: species as points
sp_umap <- umap(
  t(mat_f),
  n_neighbors = min(10, ncol(mat_f) - 1),
  min_dist    = 0.3,
  metric      = "cosine"
)
sp_umap_df <- data.frame(
  species = colnames(mat_f),
  UMAP1   = sp_umap[, 1],
  UMAP2   = sp_umap[, 2],
  cluster = factor(sp_cluster[colnames(mat_f)]),
  n_sig   = colSums(mat_f != 0)
)

pdf(file.path(DIR_FIG, "UMAP_species_by_phosphoProfile.pdf"),
    width = 7, height = 6)
ggplot(sp_umap_df, aes(UMAP1, UMAP2, color = cluster, size = n_sig)) +
  geom_point() +
  scale_colour_manual(values = paletteer_d("ggsci::default_nejm")[c(2:5)]) +
  ggrepel::geom_text_repel(aes(label = species), size = 3, max.overlaps = 30) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw(base_size = 11) +
  labs(title = "Species embedding by phosphosite-association profile")
dev.off()

# UMAP: phosphosites as points
site_umap <- umap(
  mat_f,
  n_neighbors = 15,
  min_dist    = 0.1,
  metric      = "cosine"
)
site_umap_df <- data.frame(
  feature = rownames(mat_f),
  UMAP1   = site_umap[, 1],
  UMAP2   = site_umap[, 2],
  cluster = factor(site_cluster),
  n_sig   = rowSums(mat_f != 0)
)

pdf(file.path(DIR_FIG, "UMAP_phosphosites_by_speciesProfile.pdf"),
    width = 7, height = 6)
ggplot(site_umap_df, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 2) +
  scale_colour_manual(values = paletteer_d("ggsci::signature_substitutions_cosmic")) +
  theme_bw(base_size = 11) +
  labs(title = "Phosphosite embedding by species-association profile")
dev.off()

# PCA as sanity check
pc <- prcomp(t(mat_f), scale. = FALSE)
var_exp <- pc$sdev^2 / sum(pc$sdev^2)

pc_df <- data.frame(
  species = rownames(pc$x),
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2],
  cluster = factor(sp_cluster[rownames(pc$x)])
)

pdf(file.path(DIR_FIG, "PCA_species_by_phosphoProfile.pdf"),
    width = 7, height = 6)
ggplot(pc_df, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    aes(label = species),
    size          = 3,
    max.overlaps  = Inf,
    box.padding   = 0.5,
    segment.color = "grey60",
    segment.size  = 0.3
  ) +
  theme_minimal(base_size = 11) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_exp[1] * 100),
    y = sprintf("PC2 (%.1f%%)", var_exp[2] * 100)
  )
dev.off()
