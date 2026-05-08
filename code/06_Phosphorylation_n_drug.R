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
library(OmnipathR)
library(forcats)
library(deldir)
library(sf)
library(ComplexHeatmap)
library(psych)
library(ggplot2)
library(RColorBrewer)
library(ggvenn)

set.seed(42)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_FIG <- "/data/yzwang/project/AEG_seiri/results/F4_phospho_n_drugs/"

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

sp_list <- intersect(target_sp$Species, rownames(mtx_cpm))

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

# Delegate clusters
hc_row <- hclust(as.dist(1 - cor(t(mat_f))), method = "ward.D2")
k_site <- 6
site_cluster <- cutree(hc_row, k = k_site)

hc_col <- hclust(as.dist(1 - cor(mat_f)), method = "ward.D2")
k_sp <- 4
sp_cluster <- cutree(hc_col, k = k_sp)

table(site_cluster)
table(sp_cluster)

write.csv(site_cluster, file.path(DIR_TAB, "Phospho_site_clusters.csv"))
write.csv(sp_cluster, file.path(DIR_TAB, "Phospho_species_clusters.csv"))

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
ggplot(sp_umap_df, aes(UMAP1, UMAP2, color = cluster)) +
  stat_ellipse(aes(group = cluster),
               type     = "norm",
               level    = 0.95,
               linetype = "dashed",
               linewidth = 0.5) +
  geom_point(aes(size = n_sig)) +
  scale_colour_manual(values = paletteer_d("ggsci::default_nejm")[c(2:5)]) +
  ggrepel::geom_text_repel(aes(label = species), size = 3, max.overlaps = Inf) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw(base_size = 11) +
  theme(axis.text = element_text(colour = 1))
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
  theme_bw(base_size = 11) +
  scale_colour_manual(values = paletteer_d("ggsci::signature_substitutions_cosmic")) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_exp[1] * 100),
    y = sprintf("PC2 (%.1f%%)", var_exp[2] * 100)
  )
dev.off()

############
# Annotation
# Parse phosphosite feature names
site_parsed <- data.frame(feature = rownames(mat_f)) %>%
  separate(feature, into = c("gene", "site"), sep = "_",
           remove = FALSE, extra = "merge") %>%
  mutate(
    residue  = str_extract(site, "^[STY]"),
    position = as.integer(str_extract(site, "(?<=^[STY])\\d+"))
  ) %>%
  dplyr::select(feature, gene, residue, position) %>%
  filter(!is.na(gene), !is.na(residue), !is.na(position))

# Use the gene that bears the phosphosite
gene_pathway <- import_omnipath_annotations(
  resources = c("SignaLink_pathway", "NetPath"),
  proteins  = unique(site_parsed$gene)
) %>%
  filter(label == "pathway") %>%
  transmute(gene = genesymbol, pathway = value) %>%
  distinct()

site_anno <- site_parsed %>%
  inner_join(gene_pathway, by = "gene", relationship = "many-to-many") %>%
  dplyr::select(feature, gene, pathway) %>%
  distinct()

write.csv(site_anno, file.path(DIR_TAB, "Phosphosite_pathway_annotation.csv"),
          row.names = FALSE)

sp_cluster_df <- data.frame(
  species    = names(sp_cluster),
  sp_cluster = as.integer(sp_cluster)
)

cluster_pathway_score <- mat_df %>%
  filter(p.value < 0.05) %>%
  inner_join(site_anno, by = "feature", relationship = "many-to-many") %>%
  inner_join(sp_cluster_df, by = "species") %>%
  group_by(sp_cluster, pathway) %>%
  summarise(
    score   = sum(abs(log2OR)),
    n_sites = n_distinct(feature),
    n_sp    = n_distinct(species),
    .groups = "drop"
  ) %>%
  filter(n_sites >= 3)

# Pick top pathway per cluster
cluster_label_df <- cluster_pathway_score %>%
  group_by(sp_cluster) %>%
  arrange(desc(score)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

assigned <- character(0)
top_per_cluster <- cluster_label_df %>%
  group_by(sp_cluster) %>%
  group_modify(~ {
    pick <- .x %>%
      filter(!pathway %in% assigned) %>%
      slice_head(n = 1)
    if (nrow(pick) > 0) assigned <<- c(assigned, pick$pathway)
    pick
  }) %>%
  ungroup() %>%
  dplyr::select(sp_cluster, pathway, score, n_sites)

cluster_label_vec <- setNames(top_per_cluster$pathway,
                              top_per_cluster$sp_cluster)
print(top_per_cluster)

sp_umap_df <- sp_umap_df %>%
  mutate(
    sp_cluster  = as.integer(as.character(cluster)),
    cluster_lab = paste0("Cluster ", sp_cluster, ": ",
                         cluster_label_vec[as.character(sp_cluster)])
  )

hub_colors <- paletteer_d("ggsci::default_nejm")[2:(1 + length(unique(sp_umap_df$cluster)))]

x_range <- range(sp_umap_df$UMAP1)
y_range <- range(sp_umap_df$UMAP2)
x_pad   <- diff(x_range) * 0.08
y_pad   <- diff(y_range) * 0.08

outline <- data.frame(
  x = c(x_range[1] - x_pad, x_range[2] + x_pad,
        x_range[2] + x_pad, x_range[1] - x_pad),
  y = c(y_range[1] - y_pad, y_range[1] - y_pad,
        y_range[2] + y_pad, y_range[2] + y_pad)
)

arrow_spec <- arrow(length = unit(0.25, "cm"), type = "closed")

centroids <- sp_umap_df %>%
  group_by(cluster) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2), .groups = "drop")

bbox <- c(
  xmin = x_range[1] - x_pad, xmax = x_range[2] + x_pad,
  ymin = y_range[1] - y_pad, ymax = y_range[2] + y_pad
)

dd <- deldir(centroids$UMAP1, centroids$UMAP2,
             rw = c(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]))
tiles <- tile.list(dd)

vor_df <- bind_rows(lapply(seq_along(tiles), function(i) {
  tile <- tiles[[i]]
  data.frame(
    UMAP1   = tile$x,
    UMAP2   = tile$y,
    cluster = centroids$cluster[i],
    group   = i
  )
}))

pdf(file.path(DIR_FIG, "UMAP_species_clustered.pdf"), width = 6.5, height = 6)
ggplot(sp_umap_df, aes(UMAP1, UMAP2)) +
  geom_polygon(data = vor_df, aes(UMAP1, UMAP2, fill = cluster, group = group),
               alpha = 0.2, color = NA) +
  geom_point(aes(color = cluster, size = n_sig)) +
  ggrepel::geom_text_repel(
    aes(label = species), size = 2.8,
    max.overlaps  = Inf, color = "black",
    segment.size  = 0.3, segment.color = "grey60"
  ) +
  annotate(
    "segment",
    x = x_range[1] - x_pad, xend = x_range[2] + x_pad,
    y = y_range[1] - y_pad, yend = y_range[1] - y_pad,
    arrow = arrow_spec, linewidth = 0.6, color = "black"
  ) +
  annotate(
    "segment",
    x = x_range[1] - x_pad, xend = x_range[1] - x_pad,
    y = y_range[1] - y_pad, yend = y_range[2] + y_pad,
    arrow = arrow_spec, linewidth = 0.6, color = "black"
  ) +
  scale_colour_manual(values = hub_colors) +
  scale_fill_manual  (values = hub_colors) +
  scale_size_continuous(range = c(2, 6)) +
  coord_cartesian(
    xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
    ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
    clip = "off"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    axis.text        = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank()
  ) +
  labs(x = "UMAP 1", y = "UMAP 2")
dev.off()

# Stack plot
site_grp_vec <- cutree(hc_row, k = 6)

assigned <- character(0)
grp_label <- data.frame(feature = names(site_grp_vec), grp = as.integer(site_grp_vec)) %>%
  inner_join(site_anno, by = "feature", relationship = "many-to-many") %>%
  count(grp, pathway, name = "n") %>%
  group_by(grp) %>%
  arrange(desc(n)) %>%
  group_modify(~ {
    pick <- .x %>% filter(!pathway %in% assigned) %>% slice_head(n = 1)
    if (nrow(pick) > 0) assigned <<- c(assigned, pick$pathway)
    pick
  }) %>%
  ungroup()

label_vec <- setNames(grp_label$pathway, as.character(grp_label$grp))
print(grp_label)

stacked_df <- mat_df %>%
  filter(p.value < 0.05) %>%
  inner_join(data.frame(feature = names(site_grp_vec), grp = as.integer(site_grp_vec)),
             by = "feature") %>%
  group_by(species, grp) %>%
  summarise(score = sum(abs(log2OR)), .groups = "drop_last") %>%
  mutate(prop = score / sum(score)) %>%
  ungroup() %>%
  inner_join(sp_cluster_df, by = "species") %>%
  mutate(
    grp = factor(grp, labels = paste0("G", sort(unique(grp)), ": ",
                                      label_vec[as.character(sort(unique(grp)))])),
    species = fct_reorder(species, sp_cluster)
  )

pdf(file.path(DIR_FIG, "Stacked_pathway_composition.pdf"),
    width = 8, height = max(5, length(unique(stacked_df$species)) * 0.25))
ggplot(stacked_df, aes(y = species, x = prop, fill = grp)) +
  geom_col(width = 0.8) +
  facet_grid(sp_cluster ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_brewer(palette = "Set2", name = "Pathway group") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold")) +
  labs(x = "Proportion", y = NULL)
dev.off()

###############
# Drug response
# IDWAS
# Load data
idwas_drug <- readxl::read_excel(file.path(DIR_TAB, "IDWAS_Supplemental_Table_S2.xlsx"), 
                                 skip = 1) %>%
  select(Drug)

idwas_mtx <- readRDS(file.path(DIR_RDS, "drugPred_Mtx.rds"))
idwas_tumour_mtx <- idwas_mtx[ , grep("C", colnames(idwas_mtx))]

abund <- sort(apply(mtx_cpm[, grep("C", colnames(mtx_cpm))], 
                    MARGIN = 1, FUN = mean) / 1e+4, decreasing = T)
abund_mtx_tumour <- mtx_cpm[target_sp$Species, grep("C", colnames(mtx_cpm))] / 1e+4

gdsc_target <- read.csv(file.path(DIR_TAB, "Drug_list2023.csv"))

# Plot predicted drug responses: idwas
draw_idwas_mtx <- scale(idwas_mtx, center = T)
col_zs <- circlize::colorRamp2(c(-10, 0, 10), c("#2166AC", "white", "#B2182B"))
pdf(file.path(DIR_SUP, "A_heatmap_idwas_samples.pdf"), width = 6, height = 10)
Heatmap(draw_idwas_mtx, show_column_names = F, 
        name = "Scaled response", column_title = "Samples",
        column_title_side = "bottom",
        col = col_zs,
        row_names_gp = grid::gpar(fontsize = 6))
dev.off()

rm_idwas <- c()
for (i in 1:nrow(idwas_mtx)) {
  if (length(na.omit(idwas_mtx[i, ])) < 0.5 * ncol(idwas_mtx)) {
    rm_idwas <- c(rm_idwas, i)
  }
}
idwas_mtx <- idwas_mtx[-rm_idwas, ]
dim(idwas_mtx)

# CARE
anno_dt <- rtracklayer::import("/data/yzwang/reference/gencode_ref/gencode_human_annotation.gtf") %>%
  as.data.frame()
coding_gene_list <- anno_dt$gene_name[anno_dt$gene_type == "protein_coding"] %>% unique(.)
length(coding_gene_list)
rm(anno_dt)

# Pre-processing human expression matrix
hexp_tpm <- readRDS(paste0(DIR_RDS, "AEG_humanTPM_Symbol.rds"))
hexp_tpm <- hexp_tpm[rownames(hexp_tpm) %in% coding_gene_list, ]
hexp_tpm <- hexp_tpm[rowMeans(hexp_tpm) > 1, ]
dim(hexp_tpm)

hexp_tpm <- hexp_tpm[, grep("C", colnames(hexp_tpm))]
hexp_norm <- t(apply(hexp_tpm, 1, function(gene_values) {
  (gene_values - mean(gene_values, na.rm = TRUE)) / sd(gene_values, na.rm = TRUE)
}))
gene_vars <- apply(hexp_norm, 1, var, na.rm = T)
top_var_genes <- names(gene_vars)[gene_vars > median(gene_vars, na.rm = TRUE)]
hexp_filt <- hexp_norm[top_var_genes, ]
dim(hexp_filt)

# Importing CARE scores
care_score <- read.table(file.path(DIR_TAB, "CARE_GDSC"), header = T, fill = T) %>%
  filter(Type == "t", )

meta_cols <- c("Response", "Type", "Target")
gene_cols <- setdiff(colnames(care_score), meta_cols)

care_score[, gene_cols] <- lapply(care_score[, gene_cols], function(x) {
  as.numeric(as.character(x))
})

care_aggregated <- care_score %>%
  group_by(Response) %>%
  summarise(
    Type = Type[1],
    Target = paste(Target, collapse = "/"),
    across(all_of(gene_cols), ~mean(.x, na.rm = TRUE)),
    .groups = 'drop'
  )

common_genes <- intersect(rownames(hexp_filt), colnames(care_score))
length(common_genes)

care_mtx <- care_aggregated[ , common_genes]
rownames(care_mtx) <- care_aggregated$Response

hexp_filt <- hexp_filt[common_genes, ]

# Using gene expression and CARE matrix to predict drug response per patient
n_top_genes <- 50

drug_prediction_models <- list()
care_feature_predictions <- list()

for(drug in rownames(care_mtx)) {
  # Select top genes by absolute CARE score
  drug_care <- care_mtx[drug, ]
  common_genes <- intersect(rownames(hexp_filt), names(drug_care))
  drug_care_subset <- drug_care[common_genes]
  
  top_gene_indices <- order(unlist(abs(drug_care_subset)), decreasing = TRUE)[1:n_top_genes]
  top_genes <- names(drug_care_subset)[top_gene_indices]
  
  # Feature matrix: patients × top genes
  feature_matrix <- t(hexp_filt[top_genes, , drop = FALSE])
  
  model_data <- as.data.frame(feature_matrix)
  model_data$Sample <- rownames(feature_matrix)
  
  pca_result <- prcomp(feature_matrix, scale. = TRUE, center = TRUE)
  n_pcs <- min(10, ncol(pca_result$x))
  pc_scores <- pca_result$x[, 1:n_pcs]
  
  # Weights: variance explained by each PC
  variance_explained <- pca_result$sdev[1:n_pcs]^2 / sum(pca_result$sdev^2)
  predicted_response <- pc_scores %*% variance_explained
  
  care_feature_predictions[[drug]] <- data.frame(
    Drug = drug,
    Sample = colnames(hexp_filt),
    Predicted_Response = as.numeric(predicted_response),
    N_Features = length(top_genes),
    stringsAsFactors = FALSE
  )
  
  drug_prediction_models[[drug]] <- list(
    top_genes = top_genes,
    pca_model = pca_result,
    variance_explained = variance_explained
  )
}
care_feature_pred_df <- bind_rows(care_feature_predictions)

care_feature_pred_df <- care_feature_pred_df %>%
  group_by(Drug) %>%
  mutate(Predicted_Response_Scaled = scale(Predicted_Response)[,1]) %>%
  ungroup()

# Convert long format to matrix: drugs × samples
care_pred_matrix <- care_feature_pred_df %>%
  dplyr::select(Drug, Sample, Predicted_Response_Scaled) %>%
  tidyr::pivot_wider(
    names_from = Sample,
    values_from = Predicted_Response_Scaled
  ) %>%
  tibble::column_to_rownames("Drug") %>%
  as.matrix()
dim(care_pred_matrix)
care_pred_matrix <- -care_pred_matrix # To align with IDWAS

write.csv(care_pred_matrix, file.path(DIR_TAB, "CARE_predicted_responses.csv"))

pdf(file.path(DIR_SUP, "Heatmap_care_samples.pdf"), width = 8, height = 10)
Heatmap(care_pred_matrix, name = "CARE response", show_column_names = F)
dev.off()

# Target species
# in idwas
cor_sp_idwas <- readRDS(file.path(DIR_RDS, "Correlation_species_idwas.rds"))
cor_r_idwas <- cor_sp_idwas$r
cor_p_idwas <- cor_sp_idwas$p.adj

# Drug name normalisation
rownames(cor_r_idwas) <- gsub("\\.", "-", rownames(cor_r_idwas))
rownames(cor_p_idwas) <- gsub("\\.", "-", rownames(cor_p_idwas))

drug_universe <- rownames(cor_r_idwas)

pathway_all <- gdsc_target[gdsc_target$Name %in% drug_universe,
                           c("Name", "Target.pathway")]
pathway_all <- pathway_all[!duplicated(pathway_all$Name), ]
potent_syn  <- base::setdiff(drug_universe, pathway_all$Name)

rename_map <- setNames(potent_syn, potent_syn)
notfound   <- c()
for (i in seq_along(potent_syn)) {
  hit <- grep(potent_syn[i], gdsc_target$Synonyms)
  tmp_name <- gdsc_target$Name[hit][1]
  tmp_path <- gdsc_target$Target.pathway[hit][1]
  if (is.na(tmp_name)) {
    notfound <- c(notfound, potent_syn[i])
    pathway_all <- rbind(pathway_all, c(potent_syn[i], "Unannotated"))
  } else {
    rename_map[potent_syn[i]] <- tmp_name
    pathway_all <- rbind(pathway_all, c(tmp_name, tmp_path))
  }
}

# apply rename
apply_rename <- function(m, map) {
  rn <- rownames(m)
  rn[rn %in% names(map)] <- map[rn[rn %in% names(map)]]
  rownames(m) <- rn
  m
}
cor_r_idwas <- apply_rename(cor_r_idwas, rename_map)
cor_p_idwas <- apply_rename(cor_p_idwas, rename_map)

rownames(pathway_all) <- pathway_all$Name
pathway_all <- pathway_all[!duplicated(pathway_all$Name), ]

# Filter
sig_in <- function(R, P) {
  apply(R, 1, function(x) any(abs(x) > 0.2, na.rm = TRUE)) &
    apply(P, 1, function(x) any(x < 0.05, na.rm = TRUE))
}
keep <- rownames(cor_r_idwas)[sig_in(cor_r_idwas, cor_p_idwas)]

# FDA-style filter from original code: drop rows containing "-" or digits, keep those with lowercase letters
fda_filter <- function(x) {
  x <- x[!grepl("\\-|[1-9]", x)]
  x <- x[grepl("[a-z]", x)]
  x
}
keep <- fda_filter(keep)
keep <- intersect(keep, rownames(cor_r_idwas))

cor_r_idwas <- cor_r_idwas[keep, , drop = FALSE]
cor_p_idwas <- cor_p_idwas[keep, , drop = FALSE]

# Right pathway annotation
rownames(pathway_all) <- pathway_all$Name
pathway_vec <- pathway_all[keep, "Target.pathway"]

col_pathway <- c(paletteer_d("ggsci::default_jco"),
                 paletteer_d("pals::kelly"))[seq_along(unique(pathway_vec))]
names(col_pathway) <- unique(pathway_vec)
col_pathway[is.na(names(col_pathway)) | names(col_pathway) == "Unannotated"] <- "grey80"

right_anno <- rowAnnotation(
  Pathway = pathway_vec,
  col = list(Pathway = col_pathway),
  show_annotation_name = TRUE
)

# Top annotation
build_top_anno <- function(R, P, sp_cols, show_legend = TRUE) {
  neg_v <- vapply(seq_len(ncol(R)), function(i)
    sum(R[, i] < -0.2 & P[, i] < 0.05, na.rm = TRUE), integer(1))
  pos_v <- vapply(seq_len(ncol(R)), function(i)
    sum(R[, i] >  0.2 & P[, i] < 0.05, na.rm = TRUE), integer(1))
  mean_abund <- rowMeans(abund_mtx_tumour[sp_cols, , drop = FALSE], na.rm = TRUE)
  
  HeatmapAnnotation(
    MeanAbund = mean_abund,
    Positive  = anno_barplot(pos_v,
                             gp = gpar(border = NA, fill = "#701145FF", lty = "blank")),
    Negative  = anno_barplot(neg_v,
                             gp = gpar(border = NA, fill = "#008280FF", lty = "blank")),
    col = list(MeanAbund = colorRamp2(c(0, 4), c("white", "darkgreen"))),
    show_legend = show_legend,
    show_annotation_name = show_legend
  )
}
top_anno_idwas <- build_top_anno(cor_r_idwas, cor_p_idwas,
                                 colnames(cor_r_idwas), show_legend = TRUE)

# Significant star
make_cell_fun <- function(P) {
  function(j, i, x, y, w, h, fill) {
    p <- P[i, j]
    if (is.na(p)) return(invisible(NULL))
    sym <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else NULL
    if (!is.null(sym)) {
      gb   <- textGrob(sym)
      gb_h <- convertHeight(grobHeight(gb), "mm")
      grid.text(sym, x, y - gb_h * 0.2, gp = gpar(col = 1, cex = .5))
    }
  }
}

# Color scale & column-name bottom annotation
rng <- range(cor_r_idwas, na.rm = TRUE)
col_fun <- colorRamp2(seq(-max(abs(rng)), max(abs(rng)), length.out = 11),
                      rev(brewer.pal(11, "RdBu")))

cn_idwas <- colnames(cor_r_idwas)

bottom_anno_idwas <- HeatmapAnnotation(
  text = anno_text(cn_idwas, rot = 60, location = unit(1, "npc"), just = "right"),
  annotation_height = max_text_width(cn_idwas)
)

ht_idwas <- Heatmap(
  cor_r_idwas,
  name                = "Rs",
  col                 = col_fun,
  column_title        = "iDWAS",
  row_names_gp        = gpar(fontsize = 9),
  show_row_names      = TRUE,
  show_column_names   = FALSE,
  cluster_rows        = TRUE,
  cluster_columns     = TRUE,
  top_annotation      = top_anno_idwas,
  bottom_annotation   = bottom_anno_idwas,
  right_annotation    = right_anno,
  cell_fun            = make_cell_fun(cor_p_idwas),
  show_heatmap_legend = TRUE
)

pdf(file.path(DIR_FIG, "Heatmap_idwas.pdf"), width = 6, height = 8)
draw(ht_idwas,
     merge_legends          = TRUE,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right")
dev.off()

library(ggplot2)
library(dplyr)
library(tidyr)

# Per-drug significance flags across all species
sig_pos_drug <- apply(cor_r_idwas > 0.2 & cor_p_idwas < 0.05, 1,
                      function(x) any(x, na.rm = TRUE))
sig_neg_drug <- apply(cor_r_idwas < -0.2 & cor_p_idwas < 0.05, 1,
                      function(x) any(x, na.rm = TRUE))

drug_sig_df <- data.frame(
  Drug    = rownames(cor_r_idwas),
  Pathway = pathway_vec,
  Pos     = sig_pos_drug,
  Neg     = sig_neg_drug,
  stringsAsFactors = FALSE
)

# Aggregate per pathway
pathway_count <- drug_sig_df %>%
  group_by(Pathway) %>%
  summarise(
    Positive = sum(Pos),
    Negative = sum(Neg),
    Total    = Positive + Negative,
    .groups  = "drop"
  ) %>%
  arrange(Total)

# Long format for stacked bar
pathway_long <- pathway_count %>%
  dplyr::select(Pathway, Positive, Negative, Total) %>%
  pivot_longer(cols = c("Negative", "Positive"),
               names_to = "Direction", values_to = "N") %>%
  mutate(
    Pathway   = factor(Pathway, levels = pathway_count$Pathway),
    Direction = factor(Direction, levels = c("Negative", "Positive"))
  )

cor_cols <- c(Negative = "#2E5C8A", Positive = "#E8A33D")

p_path <- ggplot(pathway_long, aes(x = N, y = Pathway, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_text(
    data = pathway_count,
    aes(x = Total, y = Pathway, label = Total),
    inherit.aes = FALSE,
    hjust = -0.3, size = 3
  ) +
  scale_fill_manual(values = cor_cols, name = "Correlation") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Number of drugs", y = NULL,
       title = "Drug count per target pathway") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    plot.title  = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(DIR_FIG, "Barplot_pathway_drug_count.pdf"),
  plot     = p_path,
  width    = 4,
  height   = 5
)

# in CARE
cor_sp_care <- readRDS(file.path(DIR_RDS, "Correlation_species_care.rds"))
cor_r_care <- cor_sp_care$r
cor_p_care <- cor_sp_care$p.adj

# Drug name normalisation
rownames(cor_r_idwas) <- gsub("\\.", "-", rownames(cor_r_idwas))
rownames(cor_p_idwas) <- gsub("\\.", "-", rownames(cor_p_idwas))
rownames(cor_r_care)  <- gsub("\\.", "-", rownames(cor_r_care))
rownames(cor_p_care)  <- gsub("\\.", "-", rownames(cor_p_care))

drug_universe <- union(rownames(cor_r_idwas), rownames(cor_r_care))

pathway_all <- gdsc_target[gdsc_target$Name %in% drug_universe,
                           c("Name", "Target.pathway")]
pathway_all <- pathway_all[!duplicated(pathway_all$Name), ]
potent_syn  <- setdiff(drug_universe, pathway_all$Name)

rename_map <- setNames(potent_syn, potent_syn)
notfound   <- c()
for (i in seq_along(potent_syn)) {
  hit <- grep(potent_syn[i], gdsc_target$Synonyms)
  tmp_name <- gdsc_target$Name[hit][1]
  tmp_path <- gdsc_target$Target.pathway[hit][1]
  if (is.na(tmp_name)) {
    notfound <- c(notfound, potent_syn[i])
    pathway_all <- rbind(pathway_all, c(potent_syn[i], "Unannotated"))
  } else {
    rename_map[potent_syn[i]] <- tmp_name
    pathway_all <- rbind(pathway_all, c(tmp_name, tmp_path))
  }
}

# apply rename to both matrices
apply_rename <- function(m, map) {
  rn <- rownames(m)
  rn[rn %in% names(map)] <- map[rn[rn %in% names(map)]]
  rownames(m) <- rn
  m
}
cor_r_idwas <- apply_rename(cor_r_idwas, rename_map)
cor_p_idwas <- apply_rename(cor_p_idwas, rename_map)
cor_r_care  <- apply_rename(cor_r_care,  rename_map)
cor_p_care  <- apply_rename(cor_p_care,  rename_map)

pathway_all <- pathway_all[!duplicated(pathway_all$Name), ]
rownames(pathway_all) <- pathway_all$Name

# Filter: significant in either panel
sig_in <- function(R, P) {
  apply(R, 1, function(x) any(abs(x) > 0.2, na.rm = TRUE)) &
    apply(P, 1, function(x) any(x < 0.05, na.rm = TRUE))
}
common_drugs <- intersect(rownames(cor_r_idwas), rownames(cor_r_care))
keep <- common_drugs[
  sig_in(cor_r_idwas[common_drugs, , drop = FALSE],
         cor_p_idwas[common_drugs, , drop = FALSE]) |
    sig_in(cor_r_care[common_drugs, , drop = FALSE],
           cor_p_care[common_drugs, , drop = FALSE])
]

# Concordance filter: overall effect sign must agree across panels
common_sp <- intersect(colnames(cor_r_idwas), colnames(cor_r_care))

concordant <- vapply(keep, function(d) {
  r1 <- cor_r_idwas[d, common_sp]
  r2 <- cor_r_care[d, common_sp]
  ok <- is.finite(r1) & is.finite(r2)
  if (sum(ok) < 3) return(FALSE)
  
  r1 <- r1[ok]; r2 <- r2[ok]
  
  # 1) overall direction: mean signs must match (and be non-trivial)
  m1 <- mean(r1); m2 <- mean(r2)
  if (sign(m1) != sign(m2) || sign(m1) == 0) return(FALSE)
  
  # 2) sign agreement at the species level: majority of signs must match
  sign_agree <- mean(sign(r1) == sign(r2))
  if (sign_agree < 0.5) return(FALSE)
  
  # 3) pattern correlation
  rho <- suppressWarnings(cor(r1, r2, method = "pearson"))
  isTRUE(rho > 0)
}, logical(1))

keep <- keep[concordant]

keep <- intersect(keep, rownames(cor_r_idwas))
keep <- intersect(keep, rownames(cor_r_care))

cor_r_idwas <- cor_r_idwas[keep, , drop = FALSE]
cor_p_idwas <- cor_p_idwas[keep, , drop = FALSE]
cor_r_care  <- cor_r_care[keep,  , drop = FALSE]
cor_p_care  <- cor_p_care[keep,  , drop = FALSE]

# Right pathway annotation
rownames(pathway_all) <- pathway_all$Name
pathway_vec <- pathway_all[keep, "Target.pathway"]

col_pathway <- c(paletteer_d("ggsci::default_jco"),
                 paletteer_d("pals::kelly"))[seq_along(unique(pathway_vec))]
names(col_pathway) <- unique(pathway_vec)
col_pathway[is.na(names(col_pathway)) | names(col_pathway) == "Unannotated"] <- "grey80"

right_anno <- rowAnnotation(
  Pathway = pathway_vec,
  col = list(Pathway = col_pathway),
  show_annotation_name = TRUE
)

# Panel
build_top_anno <- function(R, P, sp_cols, show_legend = TRUE) {
  mean_abund <- rowMeans(abund_mtx_tumour[sp_cols, , drop = FALSE], na.rm = TRUE)
  rng <- range(mean_abund, na.rm = TRUE)
  
  HeatmapAnnotation(
    log2CPM = mean_abund,
    col = list(log2CPM = colorRamp2(rng, c("white", "darkgreen"))),
    show_legend = show_legend,
    show_annotation_name = show_legend
  )
}
top_anno_idwas <- build_top_anno(cor_r_idwas, cor_p_idwas,
                                 colnames(cor_r_idwas), show_legend = TRUE)
top_anno_care  <- build_top_anno(cor_r_care,  cor_p_care,
                                 colnames(cor_r_care),  show_legend = FALSE)

# Significant star
make_cell_fun <- function(P) {
  function(j, i, x, y, w, h, fill) {
    p <- P[i, j]
    if (is.na(p)) return(invisible(NULL))
    sym <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else NULL
    if (!is.null(sym)) {
      gb   <- textGrob(sym)
      gb_h <- convertHeight(grobHeight(gb), "mm")
      grid.text(sym, x, y - gb_h * 0.2, gp = gpar(col = 1, cex = .5))
    }
  }
}

# Share color scale & column-name bottom annotation
rng <- range(c(cor_r_idwas, cor_r_care), na.rm = TRUE)
col_fun <- colorRamp2(seq(-max(abs(rng)), max(abs(rng)), length.out = 11),
                      rev(brewer.pal(11, "RdBu")))

cn_idwas <- colnames(cor_r_idwas)
cn_care  <- colnames(cor_r_care)

bottom_anno_idwas <- HeatmapAnnotation(
  text = anno_text(cn_idwas, location = unit(1, "npc"), just = "right"),
  annotation_height = max_text_width(cn_idwas)
)
bottom_anno_care <- HeatmapAnnotation(
  text = anno_text(cn_care, location = unit(1, "npc"), just = "right"),
  annotation_height = max_text_width(cn_care)
)

ht_idwas <- Heatmap(
  cor_r_idwas,
  name              = "Rs",
  col               = col_fun,
  column_title      = "iDWAS",
  row_names_gp      = gpar(fontsize = 9),
  show_column_names = FALSE,
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  top_annotation    = top_anno_idwas,
  bottom_annotation = bottom_anno_idwas,
  cell_fun          = make_cell_fun(cor_p_idwas),
  show_heatmap_legend = TRUE
)

ht_care <- Heatmap(
  cor_r_care,
  name                = "Rs",
  col                 = col_fun,
  column_title        = "CARE",
  row_names_gp        = gpar(fontsize = 9),
  show_row_names      = TRUE,
  show_column_names   = FALSE,
  cluster_rows        = FALSE,
  cluster_columns     = TRUE,
  top_annotation      = top_anno_care,
  bottom_annotation   = bottom_anno_care,
  right_annotation    = right_anno,
  cell_fun            = make_cell_fun(cor_p_care),
  show_heatmap_legend = FALSE
)

ht_list <- ht_idwas + ht_care

pdf(file.path(DIR_FIG, "Heatmap_drug_both.pdf"), width = 9, height = 7)
draw(ht_list,
     merge_legends          = TRUE,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     ht_gap = unit(4, "mm"))
dev.off()

####
# EF
cor_ef_idwas <- corr.test(t(as.matrix(idwas_mtx)), 
                          t(as.matrix(abund_mtx["Enterococcus_faecalis", ])),
                          method = "spearman")
cor_ef_idwas_r <- cor_ef_idwas$r
cor_ef_idwas_p <- cor_ef_idwas$p.adj
rownames(cor_ef_idwas_r) <- gsub("\\.", "-", rownames(cor_ef_idwas_r))

cor_ef_idwas_res <- data.frame(
  drug = rownames(cor_ef_idwas_r),
  rs = cor_ef_idwas_r[ , 1],
  padj = cor_ef_idwas_p[ , 1]
) %>%
  mutate(sig = case_when(rs > 0.3 & padj < 0.05 ~ "Pos",
                         rs < -0.3 & padj < 0.05 ~ "Neg",
                         TRUE ~ "NS"))

cor_ef_idwas_highlight <- cor_ef_idwas_res %>%
  arrange(rs)
cor_ef_idwas_highlight <- cor_ef_idwas_highlight[c(1:3, 
                                                   (nrow(cor_ef_idwas_highlight) - 2):
                                                     nrow(cor_ef_idwas_highlight)), ]

pdf(file.path(DIR_FIG, "A_volc_EF_cor_idwas.pdf"), width = 4, height = 3.5)
ggplot(cor_ef_idwas_res, aes(x = rs, y = -log10(padj), colour = sig)) + 
  ggtitle("EF correlated drugs (IDWAS)") + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#126CAA", "grey", "#9A342C")) +
  geom_point(cor_ef_idwas_highlight, shape = 21, color = "#FFB900FF",
             mapping = aes(x = rs, y = -log10(padj)),
             size = 2.7) + 
  xlim(-0.5, 0.5) +
  ggrepel::geom_label_repel(data = cor_ef_idwas_highlight,
                            aes(x = rs, y = -log10(padj), 
                                label = drug),
                            color="grey27",
                            alpha = .8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
  geom_vline(xintercept = -0.3, linetype = "dashed") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))
dev.off()

# Correlation between abundance and care-predicted responses
cor_ef_care <- corr.test(t(as.matrix(care_pred_matrix)), 
                         t(as.matrix(abund_mtx["Enterococcus_faecalis", colnames(care_pred_matrix)])),
                         method = "spearman")
cor_ef_care_r <- cor_ef_care$r
cor_ef_care_p <- cor_ef_care$p.adj

cor_ef_care_res <- data.frame(
  drug = rownames(cor_ef_care_r),
  rs = cor_ef_care_r[ , 1],
  padj = cor_ef_care_p[ , 1]
) %>%
  mutate(sig = case_when(rs > 0.3 & padj < 0.05 ~ "Pos",
                         rs < -0.3 & padj < 0.05 ~ "Neg",
                         TRUE ~ "NS"))

cor_ef_care_highlight <- cor_ef_care_res %>%
  arrange(rs)
cor_ef_care_highlight <- cor_ef_care_highlight[c(1:4, 
                                                 (nrow(cor_ef_care_highlight) - 3):
                                                   nrow(cor_ef_care_highlight)), ]

pdf(file.path(DIR_FIG, "A_volc_EF_cor_care.pdf"), width = 4, height = 3.5)
ggplot(cor_ef_care_res, aes(x = rs, y = -log10(padj), colour = sig)) + 
  ggtitle("EF correlated drugs (CARE)") + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#126CAA", "grey", "#9A342C")) +
  geom_point(cor_ef_care_highlight, shape = 21, color = "#FFB900FF",
             mapping = aes(x = rs, y = -log10(padj)),
             size = 2.7) + 
  ggrepel::geom_label_repel(data = cor_ef_care_highlight,
                            aes(x = rs, y = -log10(padj), 
                                label = drug),
                            color="grey27",
                            alpha = .8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
  geom_vline(xintercept = -0.3, linetype = "dashed") +
  theme_test() +
  theme(panel.border = element_rect(fill = NA, colour = 1),
        axis.text = element_text(colour = 1))
dev.off()

####################
# Calculate overlaps
process_drug <- function(res, sig_type) {
  drugs <- res$drug[res$sig == sig_type]
  drugs[!grepl("\\-|[1-9]", drugs)]
}

sensitive_idwas <- process_drug(cor_ef_idwas_res, "Neg")
sensitive_care <- process_drug(cor_ef_care_res, "Neg")
resistant_idwas <- process_drug(cor_ef_idwas_res, "Pos")
resistant_care <- process_drug(cor_ef_care_res, "Pos")

# Venn
venn_both <- list(
  "IDWAS Sensitive" = sensitive_idwas,
  "CARE Sensitive" = sensitive_care,
  "CARE Resistant" = resistant_care,
  "IDWAS Resistant" = resistant_idwas
)
pdf(file.path(DIR_FIG, "B_venn.pdf"), width = 6, height = 6)
ggvenn(venn_both,
       text_size = 4,
       text_color = "white",
       show_percentage = T,
       fill_color = paletteer_d("ggsci::default_jama")[c(3, 5, 2, 4)],
       fill_alpha = 0.7,
       stroke_color = "white")
dev.off()

# Dot plot showing clusters vs drug pathways
sp_cluster_df <- sp_cluster_df %>%
  dplyr::filter(species %in% colnames(cor_r_idwas)) %>%
  mutate(sp_cluster = as.integer(sp_cluster))

# Drug to pathway lookup (drop unannotated rows so the panel is interpretable)
drug_path_df <- pathway_all %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(Target.pathway), Target.pathway != "Unannotated") %>%
  dplyr::filter(Name %in% rownames(cor_r_idwas)) %>%
  distinct(Name, Target.pathway)

# Long-format correlation table: drug x species
r_long <- as.data.frame(cor_r_idwas) %>%
  tibble::rownames_to_column("Drug") %>%
  pivot_longer(-Drug, names_to = "species", values_to = "r")

p_long <- as.data.frame(cor_p_idwas) %>%
  tibble::rownames_to_column("Drug") %>%
  pivot_longer(-Drug, names_to = "species", values_to = "padj")

cor_long <- r_long %>%
  inner_join(p_long, by = c("Drug", "species")) %>%
  inner_join(drug_path_df, by = c("Drug" = "Name")) %>%
  inner_join(sp_cluster_df, by = "species") %>%
  dplyr::filter(is.finite(r), is.finite(padj))

# Aggregate at species-cluster x drug-pathway level
agg_df <- cor_long %>%
  group_by(sp_cluster, Target.pathway) %>%
  summarise(
    mean_r       = mean(r, na.rm = TRUE),
    n_total      = n(),
    n_sig_pos    = sum(r >  0.3 & padj < 0.05, na.rm = TRUE),
    n_sig_neg    = sum(r < -0.3 & padj < 0.05, na.rm = TRUE),
    n_sig        = n_sig_pos + n_sig_neg,
    n_sp         = n_distinct(species),
    n_drug       = n_distinct(Drug),
    sig_ratio    = n_sig / n_total,
    .groups = "drop"
  )

# Hypergeometric enrichment per cell
N_total <- sum(cor_long$padj < 0.05 & abs(cor_long$r) > 0.3, na.rm = TRUE)
N_pop   <- nrow(cor_long)

bg_path <- cor_long %>%
  group_by(Target.pathway) %>%
  summarise(K = n(), .groups = "drop")

bg_clu <- cor_long %>%
  group_by(sp_cluster) %>%
  summarise(M = sum(padj < 0.05 & abs(r) > 0.3, na.rm = TRUE),
            .groups = "drop")

agg_df <- agg_df %>%
  left_join(bg_path, by = "Target.pathway") %>%
  left_join(bg_clu,  by = "sp_cluster") %>%
  rowwise() %>%
  mutate(
    pval = phyper(n_sig - 1, K, N_pop - K, M, lower.tail = FALSE)
  ) %>%
  ungroup() %>%
  mutate(padj_enrich = p.adjust(pval, method = "BH"))

# Cluster labels (reuse top pathway per cluster from the upstream analysis)
if (exists("cluster_label_vec")) {
  agg_df <- agg_df %>%
    mutate(cluster_lab = paste0("C", sp_cluster, ": ",
                                cluster_label_vec[as.character(sp_cluster)]))
} else {
  agg_df <- agg_df %>% mutate(cluster_lab = paste0("Cluster ", sp_cluster))
}

# Order pathways by overall mean_r so the dotplot reads diagonally
pathway_order <- agg_df %>%
  group_by(Target.pathway) %>%
  summarise(score = sum(mean_r * n_sig, na.rm = TRUE), .groups = "drop") %>%
  arrange(score) %>%
  pull(Target.pathway)

agg_df <- agg_df %>%
  mutate(
    Target.pathway = factor(Target.pathway, levels = pathway_order),
    cluster_lab    = factor(cluster_lab,
                            levels = sort(unique(cluster_lab))),
    sig_flag       = ifelse(padj_enrich < 0.1, "*", "")
  )

write.csv(agg_df,
          file.path(DIR_TAB, "Phospho_spCluster_drugPathway_enrichment.csv"),
          row.names = FALSE)

# Symmetric color limit for diverging palette
lim <- max(abs(agg_df$mean_r), na.rm = TRUE)

p_dot <- ggplot(agg_df,
                aes(x = cluster_lab, y = Target.pathway,
                    size = n_sig, color = mean_r)) +
  geom_point() +
  geom_text(aes(label = sig_flag), color = "black",
            size = 4, vjust = 0.75, show.legend = FALSE) +
  scale_color_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits  = c(-lim, lim),
    name    = "Mean Rs"
  ) +
  scale_size_continuous(range = c(1, 8), name = "Sig. drug-species pairs") +
  labs(x = NULL, y = NULL,
       title = "Species cluster vs drug target pathway") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, color = "black"),
    axis.text.y      = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

ggsave(
  filename = file.path(DIR_FIG, "Dotplot_spCluster_vs_drugPathway.pdf"),
  plot     = p_dot,
  width    = 6,
  height   = 7
)

