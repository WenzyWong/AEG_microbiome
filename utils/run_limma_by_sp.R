library(limma)

run_limma_by_sp <- function(species_name, abund_df, prot_mtx, target_sp_df,
                            log2fc_thresh = log2(1.5), padj_thresh = 0.05,
                            log2fc_cap = 5, log10p_cap = 5) {
  
  if (!species_name %in% target_sp_df$Species) {
    warning(paste("Species not found in target_sp:", species_name))
    return(NULL)
  }
  
  if (!species_name %in% colnames(abund_df)) {
    warning(paste("Species column not found in abund_df:", species_name))
    return(NULL)
  }
  
  # Group samples by quantile (same logic as RNA-seq version)
  group_df <- abund_df %>%
    mutate(group = case_when(
      .data[[species_name]] < summary(.data[[species_name]])[2] ~ "Low",
      .data[[species_name]] > summary(.data[[species_name]])[5] ~ "High",
      TRUE ~ "Mid"
    )) %>%
    filter(group %in% c("High", "Low"))
  
  if (!all(c("High", "Low") %in% group_df$group)) {
    warning(paste("Insufficient group variation for species:", species_name))
    return(NULL)
  }
  
  # Subset and log2-transform protein matrix
  samples_use <- rownames(group_df)
  mat <- prot_mtx[, samples_use]
  mat <- log2(mat + 1)
  
  # Remove proteins with excessive missing values (>50% in either group)
  high_samples <- rownames(group_df)[group_df$group == "High"]
  low_samples  <- rownames(group_df)[group_df$group == "Low"]
  
  keep <- rowMeans(!is.na(mat[, high_samples])) >= 0.5 &
    rowMeans(!is.na(mat[, low_samples]))  >= 0.5
  mat <- mat[keep, ]
  
  # Build design matrix
  group_factor <- factor(group_df$group, levels = c("Low", "High"))
  design <- model.matrix(~ group_factor)
  
  # Impute remaining missing values (minimum value imputation)
  mat <- apply(mat, 2, function(x) {
    x[is.na(x)] <- min(x, na.rm = TRUE)
    x
  })
  
  # limma
  fit  <- lmFit(mat, design)
  fit  <- eBayes(fit)
  
  de_res <- topTable(fit, coef = "group_factorHigh",
                     number = Inf, sort.by = "P") %>%
    tibble::rownames_to_column("Protein_group")
  
  de_df <- de_res %>%
    transmute(
      Protein_group,
      log2FC  = logFC,
      p.value = P.Value,
      p.adj   = adj.P.Val
    ) %>%
    na.omit() %>%
    mutate(change = case_when(
      log2FC >  log2fc_thresh & p.adj < padj_thresh ~ "Up",
      log2FC < -log2fc_thresh & p.adj < padj_thresh ~ "Dn",
      TRUE ~ "NS"
    ))
  
  draw_de <- de_df %>%
    mutate(
      log2FC = pmax(pmin(log2FC, log2fc_cap), -log2fc_cap),
      p.adj  = if_else(-log10(p.adj) > log10p_cap,
                       10^(-log10p_cap), p.adj)
    )
  
  return(list(
    species  = species_name,
    group_df = group_df,
    de_df    = de_df,
    draw_de  = draw_de,
    fit      = fit
  ))
}