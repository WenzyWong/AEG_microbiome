#######################################################
# AEG microbiome
# Yunzhe Wang
#
# Step 05. Analysis on microbe-host crosstalk,
# based on transcriptomics and proteomics
#
#######################################################
setwd("/data/yzwang/project/AEG_seiri/")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(enrichplot)
library(paletteer)
library(patchwork)
library(compositions)
library(limma)
library(ComplexHeatmap)
library(circlize)

DIR_RDS <- "/data/yzwang/project/AEG_seiri/RDS/"
DIR_RES <- "/data/yzwang/project/AEG_seiri/results/F3_crosstalk/"
DIR_TAB <- "/data/yzwang/project/AEG_seiri/table_infos/"
DIR_TOOL <- "/data/yzwang/git_project/AEG_microbiome/utils/"

#######################
# Import and preprocess
clinical <- readxl::read_excel(file.path(DIR_TAB, "AEG_clinical.xlsx"))
dim(clinical)

mtx_hcount <- read.csv(file.path(DIR_TAB, "gene_count_matrix.csv"), row.names = 1) %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_name = gsub(".*\\|", "", gene_id)) %>%
  select(-gene_id) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum)) %>%
  filter(!grepl("^<class", gene_name)) %>%
  column_to_rownames("gene_name")
rownames(mtx_hcount)[1:3]
dim(mtx_hcount)
mtx_htpm <- read.csv(file.path(DIR_TAB, "gene_tpm_matrix.csv"), row.names = 1) %>%
  rownames_to_column("gene_id") %>%
  mutate(gene_name = gsub(".*\\|", "", gene_id)) %>%
  select(-gene_id) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum)) %>%
  filter(!grepl("^<class", gene_name)) %>%
  column_to_rownames("gene_name")
rownames(mtx_htpm)[1:3]
dim(mtx_htpm)

mtx_htpm_filt <- mtx_htpm %>%
  filter(rowMeans(.) > 1)
dim(mtx_htpm_filt)

mtx_hcount_filt <- mtx_hcount[rownames(mtx_htpm_filt), ]
dim(mtx_hcount_filt)

saveRDS(mtx_hcount_filt, file.path(DIR_RDS, "hTumourCNT_filtered.rds"))
saveRDS(mtx_htpm_filt, file.path(DIR_RDS, "hTumourTPM_filtered.rds"))

mtx_cpm <- readRDS(file.path(DIR_RDS, "sAEG_CPM_RNA_FiltMyco.rds"))
dim(mtx_cpm)

mtx_hcount_filt <- readRDS(file.path(DIR_RDS, "hTumourCNT_filtered.rds"))
mtx_htpm_filt <- readRDS(file.path(DIR_RDS, "hTumourTPM_filtered.rds"))

# Human proteomics
mtx_hpro <- readxl::read_excel(file.path(DIR_TAB, "SupData13ProteinIntensity.xlsx"))

sample_pre <- colnames(mtx_hpro)[3:ncol(mtx_hpro)]
sample_pre <- gsub("AEG", "", sample_pre) %>%
  gsub("_T$", "", .) %>%
  gsub("_N$", "", .) %>%
  { ifelse(grepl("_T$", colnames(mtx_hpro)[3:ncol(mtx_hpro)]), paste0("C", .), paste0("N", .)) }

colnames(mtx_hpro)[3:ncol(mtx_hpro)] <- sample_pre

# Overlapping (tumour)
samples_rna <- colnames(mtx_htpm_filt)[grep("C", colnames(mtx_htpm_filt))]
samples_prot <- colnames(mtx_hpro)
samples_bac <- colnames(mtx_cpm)

samples_common <- Reduce(intersect, list(samples_rna, samples_prot, samples_bac))
length(samples_common)

#######################
# Differential analysis
source(file.path(DIR_TOOL, "run_deseq_by_sp.R"))
source(file.path(DIR_TOOL, "run_gsea_by_species.R"))
source(file.path(DIR_TOOL, "extract_ranking_df.R"))
source(file.path(DIR_TOOL, "draw_summary_ridgeplot.R"))
source(file.path(DIR_TOOL, "draw_gsea_dotplot.R"))

target_sp <- read.csv(file.path(DIR_TAB, "Species_within_saving_module.csv"),
                      row.names = 1)

sample_groups <- list(
  tumour = colnames(mtx_hcount_filt)[startsWith(colnames(mtx_hcount_filt), "C")],
  normal = colnames(mtx_hcount_filt)[!startsWith(colnames(mtx_hcount_filt), "C")]
)

for (grp in names(sample_groups)) {
  samples <- sample_groups[[grp]]
  dir_rds_grp <- file.path(DIR_RDS, grp)
  dir_res_grp <- file.path(DIR_RES, grp)
  dir.create(dir_rds_grp, showWarnings = FALSE)
  dir.create(dir_res_grp, showWarnings = FALSE)
  
  # DE analysis
  abund_target <- data.frame(
    sample = samples,
    t(mtx_cpm[target_sp$Species, samples] / 1e+4)
  )
  
  all_de_results <- lapply(
    target_sp$Species,
    run_deseq_by_sp,
    abund_df     = abund_target,
    count_mtx    = mtx_hcount_filt[, samples],
    target_sp_df = target_sp
  )
  names(all_de_results) <- target_sp$Species
  all_de_results <- Filter(Negate(is.null), all_de_results)
  
  # GSEA analysis
  all_gsea_results <- lapply(
    names(all_de_results),
    run_gsea_by_species,
    de_results = all_de_results,
    dir_rds    = dir_rds_grp
  )
  names(all_gsea_results) <- names(all_de_results)
  
  # Plots
  for (gene_set in c("hallmark", "kegg")) {
    pdf(file.path(dir_res_grp, paste0("summary_gsea_ridge_", gene_set, ".pdf")),
        width = 12, height = length(all_gsea_results) * 0.8)
    print(draw_summary_ridgeplot(all_gsea_results, gene_set))
    dev.off()
    
    pdf(file.path(dir_res_grp, paste0("summary_gsea_dot_", gene_set, ".pdf")),
        width = 8, height = length(all_gsea_results))
    print(draw_gsea_dotplot(all_gsea_results, gene_set))
    dev.off()
  }
  
  # Save results
  saveRDS(all_de_results, file.path(dir_rds_grp, "all_de_results.rds"))
  saveRDS(all_gsea_results, file.path(dir_rds_grp, "all_gsea_results.rds"))
}

gsea_tumour <- readRDS(file.path(DIR_RDS, "tumour", "all_gsea_results.rds"))
gsea_normal <- readRDS(file.path(DIR_RDS, "normal", "all_gsea_results.rds"))

n_sp    <- max(length(gsea_tumour), length(gsea_normal))
n_paths <- length(unique(collect_gsea_df(
  c(gsea_tumour, gsea_normal), "hallmark")$ID))

for (gene_set in c("hallmark", "kegg")) {
  n_paths <- length(unique(collect_gsea_df(
    c(gsea_tumour, gsea_normal), gene_set)$ID))
  
  pdf(file.path(DIR_RES, paste0("compare_gsea_dot_", gene_set, ".pdf")),
      width  = n_sp * 0.6 + 4,
      height = n_paths * 0.3)
  print(draw_gsea_dotplot(gsea_tumour, gsea_normal, gene_set = gene_set))
  dev.off()
}

#############################################################
# Wether clinical traits contribute to microbial distribution
clinical_sub <- clinical %>%
  mutate(`No.` = paste0("C", `No.`)) %>%
  filter(No. %in% samples_common) %>%
  mutate(
    tumor_dims = str_extract_all(`Primary tumor size`, "[0-9.]+"),
    tumor_max_dim = map_dbl(tumor_dims, ~ max(as.numeric(.))),
    tumor_volume = map_dbl(tumor_dims, ~ prod(as.numeric(.))),
    T_stage = str_extract(`TNM stage`, "T[0-9a-z]+"),
    N_stage = str_extract(`TNM stage`, "N[0-9]+"),
    M_stage = str_extract(`TNM stage`, "M[0-9]"),
    stage_group = case_when(
      `Pathological stage` %in% c("I", "IA", "IB", "II", "IIA", "IIB") ~ "Early",
      `Pathological stage` %in% c("III", "IIIA", "IIIB", "IIIC", "IV")  ~ "Late",
      TRUE ~ NA_character_)
  )
dim(clinical_sub)

# Keep taxa present in at least 20% of samples with CPM > 1
min_samples <- ceiling(0.2 * length(samples_common))
keep <- rowSums(mtx_cpm[, samples_common] > 1) >= min_samples
mtx_cpm_filt <- mtx_cpm[keep, samples_common]

mtx_clr <- t(clr(t(mtx_cpm_filt + 0.5)))

pca_res <- prcomp(t(mtx_clr), scale. = FALSE)
pca_df <- as.data.frame(pca_res$x[, 1:3]) %>%
  tibble::rownames_to_column("No.") %>%
  left_join(clinical_sub, by = "No.")

# Example: color by TNM stage
ggplot(pca_df, aes(PC1, PC2, color = `TNM stage`)) +
  geom_point(size = 3) +
  theme_bw()

dist_aitchison <- dist(t(mtx_clr))

meta_permanova <- clinical_sub %>%
  mutate(Age_group = if_else(Age >= median(Age, na.rm = TRUE), "High", "Low")) %>%
  select(
    No.,
    Age_group,
    Sex,
    Smoking,
    Alcohol,
    `Siewert type`,
    `Borrmann classification`,
    `Lauren Classification`,
    stage_group,
    N_stage,
    M_stage,
    tumor_max_dim
  ) %>%
  column_to_rownames("No.") %>%
  drop_na()

samples_permanova <- rownames(meta_permanova)
dist_sub <- as.dist(as.matrix(dist_aitchison)[samples_permanova, samples_permanova])

adonis2(
  dist_sub ~ stage_group + `Lauren Classification` + `Borrmann classification` +
    `Siewert type` + N_stage + M_stage +
    tumor_max_dim + Age_group + Sex + Smoking + Alcohol,
  data = meta_permanova,
  permutations = 999,
  by = "margin"
) # Residue 78.5%, modest contribution

#############
# Proteomics
source(file.path(DIR_TOOL, "run_limma_by_sp.R"))

abund_target_common <- abund_target[samples_common, ]
rownames(mtx_hpro) <- mtx_hpro$Protein_group

all_de_prot <- lapply(target_sp$Species, function(sp) {
  run_limma_by_sp(
    species_name   = sp,
    abund_df       = abund_target_common,
    prot_mtx       = mtx_hpro[ , 3:ncol(mtx_hpro)],
    target_sp_df   = target_sp
  )
})
names(all_de_prot) <- target_sp$Species

all_de_prot <- lapply(all_de_prot, function(res) {
  if (is.null(res)) return(NULL)
  
  res$de_df <- res$de_df %>%
    mutate(
      idx           = as.integer(Protein_group),
      Protein_group = rownames(mtx_hpro)[idx],
      Gene_names    = mtx_hpro$Gene_names[idx]
    ) %>%
    dplyr::select(-idx)
  
  res
})

# Use genes in enriched pathway (RNA-seq) as the main source,
# with protein changes as supplements
extract_core_genes <- function(gsea_obj, gene_set = "hallmark") {
  gsea_df <- gsea_obj[[gene_set]]@result
  
  gsea_df %>%
    filter(p.adjust < 0.05) %>%
    select(ID, NES, core_enrichment) %>%
    mutate(gene = strsplit(core_enrichment, "/")) %>%
    unnest(gene) %>%
    select(pathway = ID, NES, gene)
}

core_genes_all <- lapply(names(gsea_tumour), function(sp) {
  extract_core_genes(gsea_tumour[[sp]], gene_set = "hallmark") %>%
    mutate(species = sp)
}) %>%
  bind_rows()

# Combining results
all_de_results <- readRDS(file.path(DIR_RDS, "tumour", "all_de_results.rds")) # tumour
rna_fc_all <- lapply(names(all_de_results), function(sp) {
  all_de_results[[sp]]$de_df %>%
    select(gene = symbol, rna_log2FC = log2FC, rna_change = change) %>%
    mutate(species = sp)
}) %>%
  bind_rows()
prot_fc_all <- lapply(names(all_de_prot), function(sp) {
  if (is.null(all_de_prot[[sp]])) return(NULL)
  all_de_prot[[sp]]$de_df %>%
    select(gene = Gene_names, prot_log2FC = log2FC, prot_change = change) %>%
    mutate(species = sp) %>%
    separate_rows(gene, sep = ";") %>%       # Split multi-gene protein groups
    mutate(gene = str_trim(gene)) %>%        # Remove whitespace
    distinct(gene, species, .keep_all = TRUE) # Keep one row per gene per species
}) %>%
  bind_rows()
vis_df <- core_genes_all %>%
  left_join(rna_fc_all,  by = c("gene", "species")) %>%
  left_join(prot_fc_all, by = c("gene", "species"),
            relationship = "many-to-many")

# Visualisation
# Species and gene order
species_order <- names(gsea_tumour)
rna_mat_raw <- vis_df %>%
  distinct(gene, species, .keep_all = TRUE) %>%
  select(gene, species, rna_log2FC) %>%
  pivot_wider(names_from = species, values_from = rna_log2FC) %>%
  column_to_rownames("gene")
species_use <- species_order[species_order %in% colnames(rna_mat_raw)]

# Pathway-based gene ordering: each gene sorted by its first appearing pathway
gene_pathway_df <- core_genes_all %>%
  distinct(gene, pathway) %>%
  mutate(pathway = factor(pathway, levels = sort(unique(pathway)))) %>%
  arrange(pathway, gene)

prot_genes <- unique(na.omit(unlist(
  lapply(all_de_prot, function(res) {
    if (is.null(res)) return(NULL)
    res$de_df$Gene_names
  })
)))

gene_order_primary <- gene_pathway_df %>%
  group_by(gene) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  arrange(pathway, gene) %>%
  pull(gene) %>%
  .[. %in% prot_genes]

# RNA log2FC matrix
gene_use <- gene_order_primary[gene_order_primary %in% rownames(rna_mat_raw)]
# Keep genes with data in at least 30% of species
gene_use <- gene_use[
  rowSums(!is.na(rna_mat_raw[gene_use, species_use])) >=
    ceiling(0.3 * length(species_use))
]

rna_mat  <- rna_mat_raw[gene_use, species_use] %>% as.matrix()

fc_cap <- 1.5

rna_mat_capped <- pmax(pmin(rna_mat, fc_cap), -fc_cap)

rna_col_fun <- colorRamp2(
  c(-fc_cap, -1, -0.5, 0, 0.5, 1, fc_cap),
  c("#313695", "#4575b4", "#abd9e9", "white", "#fdae61", "#f46d43", "#a50026")
)

# Proteomics log2FC matrix (same dimensions)
prot_mat_raw <- vis_df %>%
  distinct(gene, species, .keep_all = TRUE) %>%
  select(gene, species, prot_log2FC) %>%
  pivot_wider(names_from = species, values_from = prot_log2FC) %>%
  column_to_rownames("gene")

prot_mat <- matrix(NA_real_, nrow = nrow(rna_mat), ncol = ncol(rna_mat),
                   dimnames = dimnames(rna_mat))

common_genes <- intersect(rownames(prot_mat_raw), gene_use)
common_sp    <- intersect(colnames(prot_mat_raw), species_use)
prot_mat[common_genes, common_sp] <- as.matrix(
  prot_mat_raw[common_genes, common_sp]
)
prot_mat_capped <- pmax(pmin(prot_mat, fc_cap), -fc_cap)

all_pathways <- sort(unique(core_genes_all$pathway))

pathway_binary <- sapply(all_pathways, function(pw) {
  as.integer(gene_use %in% filter(core_genes_all, pathway == pw)$gene)
})
rownames(pathway_binary) <- gene_use

# Assign one color per pathway
n_pw     <- length(all_pathways)
pw_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_pw),
  all_pathways
)

# Build annotation: each pathway as a separate column, colored or white
pathway_anno_list <- lapply(all_pathways, function(pw) {
  vec <- ifelse(pathway_binary[, pw] == 1, pw, NA_character_)
  anno_simple(
    vec,
    which  = "row",
    col    = c(setNames(pw_colors[pw], pw)),
    na_col = "white",
    width  = unit(3, "mm")
  )
})
names(pathway_anno_list) <- all_pathways

row_anno <- do.call(rowAnnotation, c(
  pathway_anno_list,
  list(
    show_annotation_name = FALSE
  )
))

# Proteomics as right row annotation heatmap
prot_col_fun <- colorRamp2(
  c(-fc_cap, 0, fc_cap),
  c("#4575b4", "white", "#d73027")
)

prot_anno <- rowAnnotation(
  Proteomics = anno_simple(
    rowMeans(prot_mat_capped, na.rm = TRUE),
    col      = prot_col_fun,
    na_col   = "grey90",
    width    = unit(5, "mm")
  ),
  annotation_name_side = "top",
  annotation_name_gp   = gpar(fontsize = 8)
)

# Shorten pathway names for legend
pw_legend <- Legend(
  labels    = str_remove(all_pathways, "HALLMARK_") %>% tolower(.),
  legend_gp = gpar(fill = pw_colors),
  title     = "Pathway",
  ncol      = 2
)

rna_mat_capped[is.na(rna_mat_capped)] <- 0
ht <- Heatmap(
  as.matrix(rna_mat_capped),
  name              = "RNA\nlog2FC",
  col               = rna_col_fun,
  cluster_rows      = FALSE,
  cluster_columns   = TRUE,
  show_row_names    = FALSE,
  show_column_names = TRUE,
  row_names_gp      = gpar(fontsize = 6),
  column_names_gp   = gpar(fontsize = 8),
  column_names_rot  = 45,
  na_col            = "grey90",
  left_annotation   = row_anno,
  right_annotation  = prot_anno,
  column_title      = "RNA log2FC (core enrichment genes)",
  use_raster        = F
)

pdf(file.path(DIR_RES, "Heatmap_RNAseq_enriched_hallmark_with_protein_cor.pdf"))
draw(ht, annotation_legend_list = list(pw_legend),
     merge_legend = FALSE)
dev.off()
