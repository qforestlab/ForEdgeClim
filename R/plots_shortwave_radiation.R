#' Function to create several shortwave radiation plots
#'
#' @param sw_rad_2D Dataframe with shortwave radiative values
#' @param output_path The output_path to store the plots
#' @return Several radiation plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @importFrom scales rescale
#' @import ggplot2
#' @export
#'
plots_sw <- function(sw_rad_2D, output_path){

  # 2D RTM plots
  ##############

  # background image
  img  <- png::readPNG("Output/horizontal_transect_2.png")
  img[,,4] <- img[,,4] *0.3 # make transparent
  bg_grob <- grid::rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))

  # Calculate averages across Y-slices
  final_avg_results_2D <- sw_rad_2D |>
    group_by(X, Z) |>
    summarise(avg_F_d_down = mean(F_d_down, na.rm = TRUE),
              avg_F_d_up = mean(F_d_up, na.rm = TRUE),
              avg_F_b_down = mean(F_b_down, na.rm = TRUE), .groups = 'keep')

  # Calculate values at ground surface, i.e., Z = 1
  final_ground_surface_results_2D <- sw_rad_2D |>
    filter(Z == 1, na.rm = TRUE)

  # Calculate values Y = 25 slice
  final_slice_results_2D <- sw_rad_2D |>
    filter(Y == 25, na.rm = TRUE)

  # Caption to be plotted below the plots
  caption = paste(
    "Direct radiation vertical =", round(F_sky_dir_v,2), 'W/m²',
    "| Direct radiation horizontal =", round(F_sky_dir_h,2), 'W/m²',
    "| Diffuse radiation =", round(F_sky_diff_init,2), 'W/m²'
  )

  # F_d_down plot:
  F_d_down = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_down)) +
    geom_tile()+
    annotation_custom(
      grob = bg_grob,
      xmin = 0, xmax = 140,
      ymin = 0, ymax = 40
    ) +
    scale_fill_viridis_c(option = "inferno",
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("(b) Shortwave diffuse downward radiation"),
         x = "\n", y = "Height (m)\n", fill = bquote("Flux"~(W~m^{-2})*"      ") )+#,
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

  ggsave(paste0(output_path, '/SW_2D_F_d_down.png'), plot = F_d_down, width = 9, height = 3, dpi = 500)


  # F_d_up plot:
  F_d_up = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_up)) +
    geom_tile()+
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
    scale_fill_viridis_c(option = "inferno",
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("(c) Shortwave diffuse upward radiation"),
         x = "\nDistance from forest core (m)", y = "\n", fill = bquote("Flux"~(W~m^{-2})*"      ") )+#,
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

  ggsave(paste0(output_path, '/SW_2D_F_d_up.png'), plot = F_d_up, width = 9, height = 3, dpi = 500)

  # F_b_down plot:
  F_b_down = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_b_down)) +
    geom_tile()+
    annotation_custom(
      grob = bg_grob,
      xmin = 0, xmax = 140,
      ymin = 0, ymax = 40
    ) +
    scale_fill_viridis_c(option = "inferno",
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("(a) Shortwave direct-beam radiation"),
         fill  = bquote("Flux"~(W~m^{-2})*"      "), x = "\n", y = "\n")+
    coord_fixed(ratio = 1) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 20),
      plot.caption = element_text(hjust = 0, color = "black"),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 16),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )

  ggsave(paste0(output_path, '/SW_2D_F_b_down.png'), plot = F_b_down, width = 9, height = 3, dpi = 500, bg = "transparent")


  # F_d_down + F_b_down plot:
  F_b_d_down = ggplot(final_avg_results_2D, aes(x = X, y = Z, fill = avg_F_d_down + avg_F_b_down)) +
    geom_tile() +
    annotation_custom(
      grob = bg_grob,
      xmin = 0, xmax = 140,
      ymin = 0, ymax = 40
    ) +
    scale_fill_viridis_c(option = "inferno",
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Direct beam & diffuse downward radiation"),
         x = "Distance from forest core (m)", y = "Height (m)", fill = bquote(Flux~(W~m^{-2})) )+ #,
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

  ggsave(paste0(output_path, '/SW_2D_F_b_d_down.png'), plot = F_b_d_down, width = 9, height = 3, dpi = 500)

  # 3D plots #
  ############

  # F_down plot for 3D plot
  F_down = ggplot(final_slice_results_2D, aes(x = X, y = Z, fill = F_b_down + F_d_down)) +
    geom_tile()+
    annotation_custom(
      grob = bg_grob,
      xmin = 0, xmax = 140,
      ymin = 0, ymax = 40
    ) +
    scale_fill_viridis_c(
      option = "inferno",
      limits = c(0, 900),
      oob = scales::squish,
      guide = guide_colorbar(
        barwidth = 1,
        barheight = 10,
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    labs(title = NULL,
         x = NULL, y = NULL, fill = NULL)+
    coord_fixed(ratio = 1) +
    theme_minimal() +
    scale_y_continuous(
      position = "right"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 20),
      legend.title = element_text(size = 18),
      legend.text  = element_text(size = 16),
      legend.position = "right",
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )

  ggsave(paste0(output_path, '/3D_SW_down.pdf'), plot = F_down, width = 9, height = 3, dpi = 500, bg = "transparent")

  # F_down ground surface plot:
  F_down_surface = ggplot(final_ground_surface_results_2D, aes(x = X, y = Y, fill = F_b_down + F_d_down)) +
    geom_tile()+
    scale_fill_viridis_c(
      option = "inferno",
      limits = c(0, 900),
      oob = scales::squish,
      guide = guide_colorbar(
        barwidth = 1,
        barheight = 10,
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    labs(title = NULL,
         x = NULL, y = NULL, fill = NULL)+
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.position = 'none',
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    )

  ggsave(paste0(output_path, '/3D_SW_down_ground.pdf'), plot = F_down_surface, width = 9, height = 3, dpi = 500, bg = "transparent")

}
