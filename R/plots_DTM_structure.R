#' Function to create plot of DTM and density structure
#'
#' @param dtm Dataframe with digital terrain model
#' @param grid Dataframe with voxelized density values
#' @param output_path The output_path to store the plots
#' @return Several radiation plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @importFrom scales rescale
#' @import ggplot2
#' @export
#'
plots_dtm_struct <- function(dtm, grid, output_path){

  #########
  # PLOTS #
  #########

  # Plot DTM

  dtm_plot = ggplot(dtm, aes(x = X - min(X), y = Y - min(Y), fill = Z - min(Z))) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma",
                                  guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) + # Mooie kleurenschaal
    coord_fixed() +
    theme_bw() +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 10))) +
    labs(title = "Digital terrain model (DTM), top view",
                  x = "west     ——     east (m)",
                  y = "north\n \n |\n \n south (m)",
                  fill = "Ground\nelevation (m)")
  print(dtm_plot)

  ggsave(paste0(output_path, '/TLS_dtm.png'), plot = dtm_plot, width = 9, height = 3, dpi = 300)


  # Plot structural grid

  # Ensure output directory exists
  if (!dir.exists('Output')) dir.create('Output', recursive = TRUE)

  # Calculate the average of the density over all y-slices for each combination of x and z
  average_density <- grid |>
    group_by(X, Z) |>
    summarise(mean_density = mean(density, na.rm = TRUE), .groups = 'keep')

  # background image
  img  <- png::readPNG("Output/horizontal_transect_2.png")
  img[,,4] <- img[,,4] *0.3 # make transparent
  bg_grob <- grid::rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))

  # Add light blue sky where mean_density == 0 and don't incorporate them in the legend
  average_density$color <- ifelse(average_density$mean_density == 0, "lightblue", NA)


  density_plot = ggplot(average_density, aes(x = X, y = Z, fill = mean_density)) +
    geom_tile(aes(fill = ifelse(mean_density == 0, NA, mean_density)), na.rm = TRUE) +
    annotation_custom(
      grob = bg_grob,
      xmin = 0, xmax = 140,
      ymin = 0, ymax = 40
    ) +
    annotate("segment", x = 5, xend = 145, y = 44, yend = 44,
             arrow = arrow(ends = "both", type = "closed", length = unit(0.3, "cm")),
             linewidth = 1.2) +
    annotate("text", x = 29,  y = 40,
             label = "towards forest core", size = 5) +
    annotate("text", x = 120, y = 40,
             label = "towards forest edge", size = 5) +
    scale_fill_gradient(
      low = "lightgreen",
      high = "darkgreen",
      na.value = "lightblue",
      guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")
    ) +
    labs(title = paste("Normalised density of forest structure"),
                  x = "Distance from forest core (m)", y = "Height (m)", fill = "Density \n(unitless)")+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 16)
    )

  print(density_plot)

  ggsave(paste0(output_path, '/TLS_grid_density.png'), plot = density_plot, width = 9, height = 3, dpi = 500)

}
