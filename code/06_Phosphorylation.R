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

set.seed(42)

setwd("/data/yzwang/project/AEG_seiri/")
DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_FIG <- "/data/yzwang/project/AEG_seiri/results/F4_phospho/"

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
    legend.position  = "none",
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
