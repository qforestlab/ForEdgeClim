############################################################################################
# Script for Diagnosing Newton Iteration Convergence in ForEdgeClim
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################


library(ForEdgeClim)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

####################
# INPUT TIMESERIES #
####################

start_timeseries = Sys.time()
start_time <- as.POSIXct("2023-07-08 14:00:00", tz = "UTC")
end_time   <- as.POSIXct("2023-07-08 14:00:00", tz = "UTC")
datetime_series <- seq(start_time, end_time, by = "hour")

number_voxels = 10 # number of (randomly picked) voxels to investigate

TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_July2023.rds'

for (current_datetime in datetime_series) {

  current_datetime <- as.POSIXct(current_datetime, origin = "1970-01-01", tz = "UTC")
  message('Running model for ', current_datetime, ' ...')

  #################################################
  # INPUT DRIVERS, CONSTANTS AND MODEL PARAMETERS #
  #################################################

  create_input_drivers()
  create_physical_constants()
  create_model_parameters()

  # Create output path
  output_path <- paste0('Output/', format(current_datetime, "%Y-%m-%d_%H"))
  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

  ########################################################
  # EXTRACT THE OBSERVATIONS FROM THE INPUT DRIVER FILES #
  ########################################################

  import_RMI_observations(current_datetime)
  import_pyr_observations(current_datetime)
  import_soil_temperature(current_datetime)

  #######################
  # CREATION OF 3D GRID #
  #######################

  voxel_TLS <- readRDS(TLS_filtered_file)

  # Determine macro boundaries along X and Z
  macro_x <- max(voxel_TLS$grid$X, na.rm = TRUE)
  macro_z <- max(voxel_TLS$grid$Z, na.rm = TRUE)


  ##################
  # TRACE VOXELS   #
  ##################

  # pick random 10 voxels in the grid
  n_vox <- nrow(voxel_TLS$grid)
  set.seed(321)
  trace_idx <- round(seq(1, n_vox, length.out = number_voxels))

  #############
  # RUN MODEL #
  #############

  res <- run_foredgeclim(
    structure_grid = voxel_TLS$grid,
    datetime       = current_datetime,
    trace_idx      = trace_idx     # <--- switch on logging of iteration
  )

  saveRDS(res, file = file.path(output_path, 'model_results.rds'))

  ###############################
  # NEWTON CONVERGENCE PLOTS   #
  ###############################

  trace_df <- res$newton_trace

  if (!is.null(trace_df) && nrow(trace_df) > 0) {

    trace_df <- trace_df %>%
      dplyr::filter(!is.na(T_surface), !is.na(E_residual)) %>%
      dplyr::mutate(
        voxel_label = paste0("(", x, ",", y, ",", z, ")"),
        T_surface_C = T_surface - 273.15,
        abs_E       = abs(E_residual),
        # Distance to nearest macro boundary (max X or max Z)
        dist_macro  = pmin(macro_x - x, macro_z - z)
      )


    # Wide format is handier for separate plots
    trace_wide <- trace_df %>%
      dplyr::select(iter, voxel_label, T_surface_C, abs_E, dist_macro)


    # Left panel: surface temperature
    p_temp <- ggplot(
      trace_wide,
      aes(
        x      = iter,
        y      = T_surface_C,
        colour = dist_macro,
        group  = voxel_label
      )
    ) +
      geom_point() +
      geom_line(linewidth = 0.5) +
      scale_colour_gradientn(
        name  = "Distance to macroboundary (m)",
        colours = c("red", "orange", "lightblue", "blue"),
        values = scales::rescale(c(min(trace_wide$dist_macro),
                                   quantile(trace_wide$dist_macro, 0.25),
                                   quantile(trace_wide$dist_macro, 0.5),
                                   quantile(trace_wide$dist_macro, 0.75),
                                   max(trace_wide$dist_macro)))
      )+

      labs(
        x     = "iteration step",
        y     = expression("Forest surface temperature (" * degree * "C)"),
        title = "(a)"
      ) +
      scale_x_continuous(
        breaks = seq(1, max(trace_wide$iter), by = 2),
        limits = c(1, max(trace_wide$iter))
      ) +
      coord_cartesian(expand = TRUE) +
      theme_bw(base_size = 15)


    # Right panel: energy balance residual
    p_resid <- ggplot(
      trace_wide,
      aes(
        x      = iter,
        y      = abs_E,
        colour = dist_macro,
        group  = voxel_label
      )
    ) +
      geom_point() +
      geom_line(linewidth = 0.5) +
      scale_colour_gradientn(
        name  = "Distance to macroboundary (m)",
        colours = c("red", "orange", "lightblue", "blue"),
        values = scales::rescale(c(min(trace_wide$dist_macro),
                                   quantile(trace_wide$dist_macro, 0.25),
                                   quantile(trace_wide$dist_macro, 0.5),
                                   quantile(trace_wide$dist_macro, 0.75),
                                   max(trace_wide$dist_macro)))
      ) +

      labs(
        x = "iteration step",
        y = expression("Energy balance residual (W m"^{-2}*")"),
        title = "(b)"
      ) +
      scale_x_continuous(
        breaks = seq(1, max(trace_wide$iter), by = 2),
        limits = c(1, max(trace_wide$iter))
      ) +
      coord_cartesian(expand = TRUE) +
      theme_bw(base_size = 15)

    # Combine plots side by side
    p_conv <- p_temp + p_resid +
      plot_layout(nrow = 1, guides = "collect") &
      theme(
        plot.title = element_text(size = 25),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.key = element_rect(fill = NA),
        legend.title = element_text(size = 25, hjust = 0.5),
        legend.text  = element_text(size = 20),
        plot.margin = margin(10,10,10,10),
        axis.title       = element_text(size = 30),
        axis.text        = element_text(size = 25)
      )

    # Patchwork title
    p_conv <- p_conv +
      plot_annotation(
        title = "Newton convergence diagnostics",
        theme = theme(
          plot.title = element_text(size = 30, hjust = 0.5, margin = margin(b = 20))
        )
      )

    # legend
    p_conv <- p_conv &
      guides(
        colour = guide_colourbar(
          title          = "Distance to macroboundary (m)",
          title.position = "top",
          barwidth       = unit(6, "cm"),
          barheight      = unit(1, "cm"),
          frame.colour   = "black",
          frame.linewidth= 1,
          ticks.colour   = "black"
        )
      )

    print(p_conv)

    # Save
    ggsave(
      file.path(output_path, "newton_convergence.png"), p_conv, width  = 18, height = 11, dpi = 300
    )

  } else {
    message("No Newton trace available for plotting.")
  }

}

end_timeseries <- Sys.time()
message('Total running time timeseries = ',
        round(as.numeric(end_timeseries - start_timeseries, units = "secs"), 2), ' s')
