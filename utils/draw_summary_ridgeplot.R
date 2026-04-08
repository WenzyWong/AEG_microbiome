draw_summary_ridgeplot <- function(gsea_list, gene_set = c("hallmark", "kegg"),
                                   padj_thresh = 0.05, stat_cap = 5,
                                   fill_colors = NULL) {
  gene_set <- match.arg(gene_set)
  df <- extract_ranking_df(gsea_list, gene_set, padj_thresh)
  if (is.null(df) || nrow(df) == 0) {
    message("No significant results to plot.")
    return(invisible(NULL))
  }
  
  df <- df %>%
    group_by(species, pathway) %>%
    mutate(mean_stat = mean(stat, na.rm = TRUE)) %>%
    ungroup()
  
  species_vec <- unique(df$species)
  if (is.null(fill_colors)) {
    base_cols   <- as.character(paletteer_d("ggsci::default_jama"))
    fill_colors <- setNames(colorRampPalette(base_cols)(length(species_vec)), species_vec)
  }
  
  plots <- lapply(species_vec, function(sp) {
    df_sp <- df %>%
      filter(species == sp) %>%
      mutate(stat    = pmin(pmax(stat, -stat_cap), stat_cap),
             pathway = reorder(pathway, mean_stat))
    
    ggplot(df_sp, aes(x = stat, y = pathway, fill = species)) +
      ggridges::geom_density_ridges(alpha = 0.7, scale = 0.9,
                                    rel_min_height = 0.01,
                                    color = "white", linewidth = 0.3) +
      scale_fill_manual(values = fill_colors) +
      scale_x_continuous(limits = c(-stat_cap, stat_cap), name = NULL) +
      scale_y_discrete(labels = tolower) +
      labs(y = NULL, title = sp) +
      theme_bw(base_size = 9) +
      theme(plot.title         = element_text(face = "bold", size = 8, hjust = 0.5),
            axis.text.y        = element_text(size = 6),
            axis.text.x        = element_text(size = 7),
            panel.grid.major.y = element_blank(),
            legend.position    = "none")
  })
  
  wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title   = paste("GSEA Summary -", toupper(gene_set)),
      caption = "Ranked Statistic (log2FC)",
      theme   = theme(plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
                      plot.caption = element_text(size = 9, hjust = 0.5))
    )
}