###############################################################################
#
# In this script multiple timepoints between model and observations are plot along the
# vertical transect line.
#
###############################################################################

library(dplyr)
library(ggplot2)
library(gganimate)
library(readxl)
library(lubridate)
library(ggimage)
library(png)
library(grid)
library(abind)

#########
# INPUT #
#########

model_times <- lubridate::ymd_h(c("2023-07-08 01", "2023-07-08 08", "2023-07-08 12"))

res_list   <- list(
  readRDS("Output/2023-07-08_01/model_results_calibrated.rds"),
  readRDS("Output/2023-07-08_08/model_results_calibrated.rds"),
  readRDS("Output/2023-07-08_12/model_results_calibrated.rds")
)

res_list_cal <- list(
  readRDS("Output/2023-07-08_01/model_results_calibrated_year.rds"),
  readRDS("Output/2023-07-08_08/model_results_calibrated_year.rds"),
  readRDS("Output/2023-07-08_12/model_results_calibrated_year.rds")
)


tomst_files <- c(
  "Data/model_TOMST_result_files/TOMST_filtered_height_temp_20230708_0100.csv",
  "Data/model_TOMST_result_files/TOMST_filtered_height_temp_20230708_0800.csv",
  "Data/model_TOMST_result_files/TOMST_filtered_height_temp_20230708_1200.csv"
)

tomst_times <- model_times

#########
# DATA  #
#########

# Model data (Seasonally calibrated fit) in long-form (x == 75)
model_df <- bind_rows(
  lapply(seq_along(res_list), function(i) {
    df   <- res_list[[i]]$micro_grid
    Tair <- res_list[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z <= 35, y == 18, x == 75) %>%
      transmute(
        z           = z,
        temperature = Tair - 273.15,
        model       = "Modelled",
        variant     = "Seasonally calibrated fit",
        time        = model_times[i]
      )
  })
)

# Model data (Annually calibrated fit) in long-form (x == 75)
model_df_cal <- bind_rows(
  lapply(seq_along(res_list_cal), function(i) {
    df   <- res_list_cal[[i]]$micro_grid
    Tair <- res_list_cal[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z <= 35, y == 18, x == 75) %>%
      transmute(
        z           = z,
        temperature = Tair - 273.15,
        model       = "Modelled",
        variant     = "Annually calibrated fit",
        time        = model_times[i]
      )
  })
)


# TOMST data in long-form
tomst_df <- bind_rows(
  lapply(seq_along(tomst_files), function(i) {
    read.csv(tomst_files[i]) %>%
      transmute(
        z           = height,
        temperature = Tair,
        model       = "Observed",
        time        = tomst_times[i]
      )
  })
)

anim_df <- bind_rows(model_df, model_df_cal, tomst_df) %>%
  mutate(
    variant  = ifelse(model == "Observed", "Observed", variant),
    time     = lubridate::floor_date(time, unit = "hour"),
    hour_lbl = format(time, "%Hh"),
    time_lbl = factor(hour_lbl,
                      levels = c("01h", "08h", "12h"),
                      labels = c("01h (~night)", "08h (~morning)", "12h (~noon)"))
  )


temp_min <- floor(min(anim_df$temperature, na.rm = TRUE))
temp_max <- ceiling(max(anim_df$temperature, na.rm = TRUE))

####################
# BACKGROUND IMAGE #
####################

bg_img <- png::readPNG("Output/vertical_transect.png")

if (dim(bg_img)[3] == 4) {
  bg_img[,,4] <- bg_img[,,4] * 0.4  # transparent
} else {
  alpha_channel <- array(0.4, dim = c(dim(bg_img)[1], dim(bg_img)[2]))
  bg_img <- abind(bg_img, alpha_channel, along = 3)
}

bg_grob <- rasterGrob(bg_img, width = unit(1, "npc"), height = unit(1, "npc"))

temp_mid <- (temp_min + temp_max) / 2
width_range <- 8
xmin_bg <- temp_mid - width_range / 2
xmax_bg <- temp_mid + width_range / 2

############
#   PLOT   #
############

cap_h <- 0.35

p_vert <- ggplot() +
  # background
  annotation_custom(
    grob = bg_grob,
    xmin = xmin_bg,
    xmax = xmax_bg,
    ymin = 0,
    ymax = 35
  ) +

  # Model as path (uncal + cal)
  geom_path(
    data = dplyr::filter(anim_df, model == "Modelled"),
    aes(
      x = temperature, y = z,
      colour = time_lbl,
      linetype = variant,
      group = interaction(model, variant, time_lbl)
    ),
    linewidth = 1
  ) +

  # TOMST errorbars (±0.5 °C)
  geom_segment(
    data = dplyr::filter(anim_df, model == "Observed"),
    inherit.aes = FALSE,
    aes(x = temperature - 0.5, xend = temperature + 0.5,
        y = z, yend = z,
        colour = time_lbl,
        group = interaction(model, time_lbl)),
    linewidth = 0.6, alpha = 0.9, show.legend = FALSE
  ) +
  geom_segment(
    data = dplyr::filter(anim_df, model == "Observed"),
    inherit.aes = FALSE,
    aes(x = temperature - 0.5, xend = temperature - 0.5,
        y = z - cap_h, yend = z + cap_h,
        colour = time_lbl,
        group = interaction(model, time_lbl)),
    linewidth = 0.6, alpha = 0.9, show.legend = FALSE
  ) +
  geom_segment(
    data = dplyr::filter(anim_df, model == "Observed"),
    inherit.aes = FALSE,
    aes(x = temperature + 0.5, xend = temperature + 0.5,
        y = z - cap_h, yend = z + cap_h,
        colour = time_lbl,
        group = interaction(model, time_lbl)),
    linewidth = 0.6, alpha = 0.9, show.legend = FALSE
  ) +

  # TOMST-points
  geom_point(
    data = dplyr::filter(anim_df, model == "Observed"),
    aes(
      x = temperature, y = z,
      colour = time_lbl,
      shape = model,
      group = interaction(model, time_lbl)
    ),
    size = 2.5
  ) +

  scale_colour_manual(
    name = "Time",
    values = c("01h (~night)" = "blue",
               "08h (~morning)" = "orange",
               "12h (~noon)"    = "red")
  ) +
  scale_linetype_manual(
    name = "Model output",
    values = c("Seasonally calibrated fit" = "solid",
               "Annually calibrated fit"   = "dashed")
  ) +
  scale_shape_manual(
    name = "Observations",
    values = c("Observed" = 16)
  ) +
  guides(
    colour   = guide_legend(order = 1, title = "Time", direction = "horizontal"),
    linetype = guide_legend(order = 2, title = "Model", direction = "horizontal",
                            keywidth = unit(3.5, "cm"),
                            override.aes = list(linewidth = 2)),
    shape    = guide_legend(order = 3, title = "Data", direction = "horizontal")
  ) +
  coord_cartesian(ylim = c(0, 35), xlim = c(temp_min, temp_max)) +
  labs(x = "Temperature (°C)", y = "Height (m)", title = "(b) Vertical air temperature") +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(size = 35),
    legend.position   = "top",
    legend.box        = "vertical",
    axis.title      = element_text(size = 35),
    axis.text       = element_text(size = 35),
    legend.title      = element_text(size = 35),
    legend.text       = element_text(size = 35),
    panel.background  = element_rect(fill = NA, colour = NA),
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA)
  )


p_vert <- p_vert +
  #coord_cartesian(ylim = c(0, 38), xlim = c(temp_min, temp_max), clip = "off") +
  annotate("segment",
           x = temp_min + 0.5, xend = temp_min + 2,
           y = 37, yend = 37,
           colour = "black", linewidth = 1) +
  annotate("segment",
           x = temp_min + 3, xend = temp_min + 4.5,
           y = 37, yend = 37,
           colour = "black", linewidth = 1) +
  theme(legend.position = "none",
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10))


ggsave("Output/Tair_vertical.png", p_vert, width = 9, height = 11, dpi = 300)
