#' Calculate and visualize beta diversity using PCoA with top density plot
#'
#' @param abundance_list A named list of abundance matrices (genes/features in rows, samples in columns)
#' @param colors A vector of colors for each group (optional, will use default palette if not provided)
#' @param method Distance method for vegdist (default: "bray")
#' @param ellipse_level Confidence level for ellipses (default: 0.95)
#' @param ellipse_alpha Alpha transparency for ellipses (default: 0.05)
#' @param point_size Size of points (default: 0.7)
#' @param line_width Width of connecting lines to centroids (default: 0.2)
#' @param show_density Logical, whether to show top density plot (default: TRUE)
#' @param density_alpha Alpha transparency for density fills (default: 0.05)
#' @param density_linewidth Line width for density curves (default: 1)
#' @param legend_position Position of legend: "right" (default), "top", "bottom", "left", or "none"
#' @param output_file Optional file path to save the plot (default: NULL)
#' @param plot_width Width of output plot in inches (default: 5)
#' @param plot_height Height of output plot in inches (default: 5)
#'
#' @return A list containing:
#'   - plot: Combined ggplot object
#'   - main_plot: Main PCoA plot
#'   - density_plot: Density plot
#'   - pcoa_result: PCoA results from ape::pcoa
#'   - distance_matrix: Distance matrix
#'   - plot_data: Data frame used for plotting
#'
#' @examples
#' result <- plot_beta_diversity(
#'   abundance_list = list(AEG = cpm_aeg, ESCC = cpm_escc, STAD = cpm_gc),
#'   colors = c("#485682", "#60AB9E", "#5C8447"),
#'   show_density = TRUE,
#'   output_file = "beta_diversity.pdf"
#' )
#' 
plot_beta_diversity <- function(abundance_list,
                                colors = NULL,
                                method = "bray",
                                ellipse_level = 0.95,
                                ellipse_alpha = 0.05,
                                point_size = 0.7,
                                line_width = 0.2,
                                show_density = TRUE,
                                density_alpha = 0.05,
                                density_linewidth = 1,
                                legend_position = "right",
                                output_file = NULL,
                                plot_width = 5,
                                plot_height = 5) {
  
  # Load required libraries
  require(vegan)
  require(ape)
  require(ggplot2)
  if (show_density) {
    require(cowplot)
  }
  
  # Input validation
  if (!is.list(abundance_list) || length(abundance_list) < 2) {
    stop("abundance_list must be a named list with at least 2 matrices")
  }
  
  if (is.null(names(abundance_list)) || any(names(abundance_list) == "")) {
    stop("All elements in abundance_list must be named")
  }
  
  n_groups <- length(abundance_list)
  group_names <- names(abundance_list)
  
  # Set default colors if not provided
  if (is.null(colors)) {
    # Use a color-blind friendly palette
    if (n_groups <= 8) {
      colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                  "#0072B2", "#D55E00", "#CC79A7", "#999999")[1:n_groups]
    } else {
      colors <- rainbow(n_groups)
    }
  } else if (length(colors) != n_groups) {
    stop(paste("Number of colors (", length(colors), 
               ") must match number of groups (", n_groups, ")", sep = ""))
  }
  
  # Find overlapping features across all matrices
  overlap <- rownames(abundance_list[[1]])
  for (i in 2:n_groups) {
    overlap <- intersect(overlap, rownames(abundance_list[[i]]))
  }
  
  if (length(overlap) == 0) {
    stop("No overlapping features found across all abundance matrices")
  }
  
  message(paste("Found", length(overlap), "overlapping features"))
  
  # Subset matrices to overlapping features
  subset_list <- lapply(abundance_list, function(x) x[overlap, , drop = FALSE])
  
  # Combine all matrices
  combined_matrix <- do.call(cbind, subset_list)
  
  # Calculate distance matrix
  dist_matrix <- vegdist(t(combined_matrix), method = method)
  dist_mtx <- as.matrix(dist_matrix)
  
  # Perform PCoA
  pcoa_result <- pcoa(dist_matrix)
  
  # Extract variance explained
  pcoa_var <- pcoa_result$values[, "Relative_eig"]
  pcoa1_var <- as.numeric(sprintf("%.3f", pcoa_var[1])) * 100
  pcoa2_var <- as.numeric(sprintf("%.3f", pcoa_var[2])) * 100
  
  # Prepare plotting data
  plot_data <- as.data.frame(pcoa_result$vectors[, 1:2])
  plot_data$Samples <- colnames(combined_matrix)
  
  # Assign groups
  plot_data$Group <- rep(group_names, sapply(subset_list, ncol))
  plot_data$Group <- factor(plot_data$Group, levels = group_names)
  
  # Calculate group centroids
  centroid_data <- data.frame(
    X = sapply(group_names, function(g) mean(plot_data$Axis.1[plot_data$Group == g])),
    Y = sapply(group_names, function(g) mean(plot_data$Axis.2[plot_data$Group == g])),
    Group = group_names
  )
  
  # Add centroid coordinates to plot_data for segments
  plot_data$centroid_x <- centroid_data$X[match(plot_data$Group, centroid_data$Group)]
  plot_data$centroid_y <- centroid_data$Y[match(plot_data$Group, centroid_data$Group)]
  
  # Determine legend direction based on position
  legend_direction <- ifelse(legend_position %in% c("top", "bottom"), "horizontal", "vertical")
  
  # Create main plot
  p_main <- ggplot(plot_data, aes(x = Axis.1, y = Axis.2, color = Group)) +
    geom_point(size = point_size) +
    geom_segment(aes(x = Axis.1, xend = centroid_x,
                     y = Axis.2, yend = centroid_y),
                 linewidth = line_width) +
    scale_color_manual(values = setNames(colors, group_names)) +
    labs(x = paste("PCoA1 (", pcoa1_var, "%)", sep = ""),
         y = paste("PCoA2 (", pcoa2_var, "%)", sep = ""),
         title = paste("PCoA:", tools::toTitleCase(method), "Distance")) +
    geom_hline(yintercept = 0, linetype = 4, color = "grey20", alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = 4, color = "grey20", alpha = 0.6) +
    stat_ellipse(geom = "polygon", level = ellipse_level, alpha = ellipse_alpha) +
    theme_classic() +
    theme(axis.text = element_text(colour = 1),
          legend.position = legend_position,
          legend.direction = legend_direction)
  
  # Create density plot if requested
  if (show_density) {
    # Get x-axis limits from main plot
    x_limits <- layer_scales(p_main)$x$range$range
    
    # Top density plot (PCoA1 only)
    p_density <- ggplot(plot_data, aes(x = Axis.1, fill = Group, color = Group)) +
      geom_density(alpha = density_alpha, linewidth = density_linewidth) +
      scale_fill_manual(values = setNames(colors, group_names)) +
      scale_color_manual(values = setNames(colors, group_names)) +
      scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
      labs(y = "Density") +
      theme_classic() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(colour = 1),
        legend.position = "none",
        plot.margin = margin(5, 5, 0, 5)
      )
    
    # Remove title from main plot for combined layout
    p_main_adjusted <- p_main + 
      labs(title = NULL) +
      theme(
        legend.position = "none",
        plot.margin = margin(0, 5, 5, 5)
      )
    
    # Extract legend with specified direction
    legend <- suppressWarnings(
      get_legend(
        p_main + 
          theme(
            legend.position = legend_position,
            legend.direction = legend_direction
          )
      )
    )
    # Combine density and main plots vertically
    p_combined <- plot_grid(
      p_density,
      p_main_adjusted,
      ncol = 1,
      align = "v",
      axis = "lr",
      rel_heights = c(0.8, 3)
    )
    
    # Add title
    title <- ggdraw() + 
      draw_label(
        paste("PCoA:", tools::toTitleCase(method), "Distance"),
        fontface = 'bold',
        x = 0,
        hjust = 0,
        size = 11
      ) +
      theme(
        plot.margin = margin(5, 5, 0, 5)
      )
    
    # Combine title with plots
    p_with_title <- plot_grid(
      title,
      p_combined,
      ncol = 1,
      rel_heights = c(0.1, 1)
    )
    
    # Add legend based on position
    if (legend_position == "right") {
      p_final <- plot_grid(
        p_with_title,
        legend,
        ncol = 2,
        rel_widths = c(1, 0.15)
      )
    } else if (legend_position == "left") {
      p_final <- plot_grid(
        legend,
        p_with_title,
        ncol = 2,
        rel_widths = c(0.15, 1)
      )
    } else if (legend_position == "top") {
      p_final <- plot_grid(
        legend,
        p_combined,
        ncol = 1,
        rel_heights = c(0.1, 1)
      )
      # Add title separately
      p_final <- plot_grid(
        title,
        p_final,
        ncol = 1,
        rel_heights = c(0.1, 1)
      )
    } else if (legend_position == "bottom") {
      p_final <- plot_grid(
        p_with_title,
        legend,
        ncol = 1,
        rel_heights = c(1, 0.1)
      )
    } else {
      # legend_position == "none"
      p_final <- p_with_title
    }
    
    final_plot <- p_final
    density_plot <- p_density
  } else {
    final_plot <- p_main
    density_plot <- NULL
  }
  
  # Save plot if output file is specified
  if (!is.null(output_file)) {
    pdf(output_file, height = plot_height, width = plot_width)
    print(final_plot)
    dev.off()
    message(paste("Plot saved to:", output_file))
  }
  
  # Return results
  return(list(
    plot = final_plot,
    main_plot = p_main,
    density_plot = density_plot,
    pcoa_result = pcoa_result,
    distance_matrix = dist_mtx,
    plot_data = plot_data,
    centroid_data = centroid_data
  ))
}

# Example usage with your original data (default: legend on right, vertical):
# result <- plot_beta_diversity(
#   abundance_list = list(AEG = cpm_aeg, ESCC = cpm_escc, STAD = cpm_gc),
#   colors = c("#485682", "#60AB9E", "#5C8447"),
#   show_density = TRUE,
#   output_file = file.path(DIR_RES, "B_beta_diversity_among_cancer.pdf")
# )
# 
# # With legend on top (horizontal):
# result2 <- plot_beta_diversity(
#   abundance_list = list(AEG = cpm_aeg, ESCC = cpm_escc, STAD = cpm_gc),
#   colors = c("#485682", "#60AB9E", "#5C8447"),
#   show_density = TRUE,
#   legend_position = "top",
#   output_file = file.path(DIR_RES, "B_beta_diversity_top_legend.pdf")
# )
# 
# # Without density plot (original style):
# result3 <- plot_beta_diversity(
#   abundance_list = list(AEG = cpm_aeg, ESCC = cpm_escc, STAD = cpm_gc),
#   colors = c("#485682", "#60AB9E", "#5C8447"),
#   show_density = FALSE,
#   legend_position = "top",
#   output_file = file.path(DIR_RES, "B_beta_diversity_simple.pdf"),
#   plot_width = 3,
#   plot_height = 3.4
# )