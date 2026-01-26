#' Function to create several longwave radiation plots
#'
#' @param lw_rad_2D Dataframe with longwave radiative values
#' @param output_path The output_path to store the plots
#' @return Several radiation plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @importFrom scales rescale
#' @import ggplot2
#' @export
#'
plots_lw <- function(lw_rad_2D, output_path){

  # 2D RTM plots
  ##############

  # background image
  img  <- png::readPNG("Output/horizontal_transect_2.png")
  img[,,4] <- img[,,4] *0.3 # make transparent
  bg_grob <- grid::rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))


  # Calculate averages across Y-slices
  final_avg_results_2D <- lw_rad_2D |>
    group_by(X, Z) |>
    summarise(avg_F_d_down = mean(F_d_down, na.rm = TRUE),
              avg_F_d_up = mean(F_d_up, na.rm = TRUE), .groups = 'keep')

  # Caption to be plotted below the plots
  caption = paste0(
    "Diffuse LW radiation = ", round(F_sky_lw,2), ' W/m²'
  )


  # F_d_down plot:
  F_d_down = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_down)) +
    geom_tile() +
    scale_x_continuous(
      breaks = seq(0, 150, by = 50),
      labels = rev(seq(0, 150, by = 50))
    ) +
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
    scale_fill_viridis_c(option = "inferno", limits = c(150, 450), breaks = c(200, 300, 400),
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("(a) Longwave downward radiation"),
         x = "\n", y = "Height (m)\n", fill = bquote("Flux"~(W~m^{-2})*"      " ) ) +#,
         #caption = caption) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20),
      plot.caption = element_text(hjust = 0, color = "black"),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 16)
    )

  print(F_d_down)

  ggsave(paste0(output_path, '/LW_2D_F_d_down.png'), plot = F_d_down, width = 9, height = 3, dpi = 500)


  # F_d_up plot:
  F_d_up = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_up)) +
    geom_tile() +
    scale_x_continuous(
      breaks = seq(0, 150, by = 50),
      labels = rev(seq(0, 150, by = 50))
    ) +
    annotation_custom(
      grob = bg_grob,
      xmin = 0, xmax = 140,
      ymin = 0, ymax = 40
    ) +
    scale_fill_viridis_c(option = "inferno", limits = c(150, 450), breaks = c(200, 300, 400),
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("(b) Longwave upward radiation"),
         x = "\nDistance from forest edge (m)", y = "\n", fill = bquote("Flux"~(W~m^{-2})*"      " ) )+ #,
         #caption = caption) +
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20),
      plot.caption = element_text(hjust = 0, color = "black"),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 16)
    )

  print(F_d_up)

  ggsave(paste0(output_path, '/LW_2D_F_d_up.png'), plot = F_d_up, width = 9, height = 3, dpi = 500)




}
