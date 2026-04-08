collect_gsea_df <- function(gsea_list, gene_set = c("hallmark", "kegg")) {
  gene_set <- match.arg(gene_set)
  
  bind_rows(lapply(names(gsea_list), function(sp) {
    obj <- gsea_list[[sp]][[gene_set]]
    if (is.null(obj) || nrow(obj@result) == 0) return(NULL)
    
    obj@result %>%
      select(ID, NES, p.adjust, enrichmentScore, setSize) %>%
      mutate(species = sp)
  }))
}
draw_gsea_dotplot <- function(gsea_list, gene_set = c("hallmark", "kegg"),
                              padj_thresh = 0.05, n_top = NULL) {
  gene_set <- match.arg(gene_set)
  
  df <- collect_gsea_df(gsea_list, gene_set)
  if (is.null(df) || nrow(df) == 0) {
    message("No significant results to plot.")
    return(invisible(NULL))
  }
  
  if (!is.null(n_top)) {
    top_paths <- df %>%
      filter(p.adjust < padj_thresh) %>%
      count(ID, sort = TRUE) %>%
      slice_head(n = n_top) %>%
      pull(ID)
    df <- df %>% filter(ID %in% top_paths)
  }
  
  path_order <- df %>%
    filter(p.adjust < padj_thresh) %>%
    group_by(ID) %>%
    summarise(
      n_pos  = sum(NES > 0),
      n_neg  = sum(NES < 0),
      score  = n_pos - n_neg,
      .groups = "drop"
    ) %>%
    arrange(score) %>%
    pull(ID)
  
  # Pathways with no significant hits in any species go to bottom
  all_paths    <- unique(df$ID)
  unsig_paths  <- setdiff(all_paths, path_order)
  path_order   <- c(unsig_paths, path_order)
  
  # Fix species order to match input
  species_order <- names(gsea_list)
  
  df <- df %>%
    mutate(
      ID        = factor(ID, levels = path_order),
      species   = factor(species, levels = species_order),
      signif    = p.adjust < padj_thresh
    )
  
  nes_bound <- max(abs(df$NES), na.rm = TRUE)
  
  ggplot(df, aes(x = species, y = ID)) +
    geom_point(aes(size = -log10(p.adjust), color = NES, alpha = signif)) +
    scale_color_gradientn(
      colors = c("#3366CC", "white", "#CC3333"),
      limits = c(-nes_bound, nes_bound),
      name   = "NES"
    ) +
    scale_size_continuous(
      name  = expression(-log[10](p.adj)),
      range = c(1, 6)
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.2),
                       guide  = "none") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL,
         title = paste("GSEA Dotplot -", toupper(gene_set))) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 0, size = 8),
      axis.text.y      = element_text(size = 7),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position  = "right"
    )
}