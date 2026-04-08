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
draw_gsea_dotplot <- function(gsea_list_a, gsea_list_b,
                              label_a = "tumour", label_b = "normal",
                              gene_set = c("hallmark", "kegg"),
                              padj_thresh = 0.05, n_top = NULL) {
  gene_set <- match.arg(gene_set)
  
  df <- bind_rows(
    collect_gsea_df(gsea_list_a, gene_set) %>% mutate(group = label_a),
    collect_gsea_df(gsea_list_b, gene_set) %>% mutate(group = label_b)
  )
  
  if (nrow(df) == 0) {
    message("No results to plot.")
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
  
  # Y-axis order based on combined score across both groups
  path_order <- df %>%
    filter(p.adjust < padj_thresh) %>%
    group_by(ID) %>%
    summarise(score = sum(NES > 0) - sum(NES < 0), .groups = "drop") %>%
    arrange(score) %>%
    pull(ID)
  
  unsig_paths <- setdiff(unique(df$ID), path_order)
  path_order  <- c(unsig_paths, path_order)
  
  species_order <- union(names(gsea_list_a), names(gsea_list_b))
  
  df <- df %>%
    mutate(
      ID      = factor(ID, levels = path_order),
      species = factor(species, levels = species_order),
      group   = factor(group, levels = c(label_a, label_b)),
      signif  = p.adjust < padj_thresh
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
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.2), guide = "none") +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(labels = tolower) +
    facet_grid(. ~ group, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL,
         title = paste("GSEA Dotplot -", toupper(gene_set))) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x      = element_text(color = 1, angle = 90, hjust = 1, vjust = .5),
      axis.text.y      = element_text(color = 1),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold", size = 10),
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position  = "right"
    )
}