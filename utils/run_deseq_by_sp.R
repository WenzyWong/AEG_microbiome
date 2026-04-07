run_deseq_by_sp <- function(species_name, abund_df, count_mtx, target_sp_df,
                                  log2fc_thresh = log2(1.5), padj_thresh = 0.05,
                                  log2fc_cap = 5, log10p_cap = 5) {
  # Check if species exists in target_sp
  if (!species_name %in% target_sp_df$Species) {
    warning(paste("Species not found in target_sp:", species_name))
    return(NULL)
  }
  
  # Check if species column exists in abund_df
  if (!species_name %in% colnames(abund_df)) {
    warning(paste("Species column not found in abund_df:", species_name))
    return(NULL)
  }
  
  # Group samples by species abundance quantile
  group_df <- abund_df %>%
    mutate(group = case_when(
      .data[[species_name]] < summary(.data[[species_name]])[2] ~ "Low",
      .data[[species_name]] > summary(.data[[species_name]])[5] ~ "High",
      TRUE ~ "Mid"
    ))
  
  # Check that High and Low both exist
  if (!all(c("High", "Low") %in% group_df$group)) {
    warning(paste("Insufficient group variation for species:", species_name))
    return(NULL)
  }
  
  # DESeq2
  df_hcount <- cbind(gene.name = rownames(count_mtx), count_mtx)
  
  dds <- DESeqDataSetFromMatrix(
    countData = df_hcount,
    colData   = group_df,
    design    = ~ group,
    tidy      = TRUE
  )
  dds <- DESeq(dds)
  
  de_res <- results(dds, contrast = c("group", "High", "Low"))
  de_res <- de_res[order(de_res$padj), ] %>% as.data.frame()
  
  # Build result data frame
  de_df <- data.frame(
    symbol = rownames(de_res),
    log2FC = de_res$log2FoldChange,
    p.adj  = de_res$padj
  ) %>%
    na.omit() %>%
    mutate(change = case_when(
      log2FC >  log2fc_thresh & p.adj < padj_thresh ~ "Up",
      log2FC < -log2fc_thresh & p.adj < padj_thresh ~ "Dn",
      TRUE ~ "NS"
    ))
  
  # Cap extreme values for visualization
  draw_de <- de_df %>%
    mutate(
      log2FC = case_when(
        log2FC >  log2fc_cap ~ log2fc_cap,
        log2FC < -log2fc_cap ~ -log2fc_cap,
        TRUE ~ log2FC
      ),
      p.adj = if_else(-log10(p.adj) > log10p_cap, 10^(-log10p_cap), p.adj)
    )
  
  return(list(
    species  = species_name,
    group_df = group_df,
    de_df    = de_df,
    draw_de  = draw_de,
    dds      = dds
  ))
}
