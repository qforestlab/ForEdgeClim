###############################################################################
#
# In this script multiple timepoints between model and observations are plot along the
# horizontal transect line.
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
  "Data/model_TOMST_result_files/TOMST_filtered_distance_temp_20230708_0100.csv",
  "Data/model_TOMST_result_files/TOMST_filtered_distance_temp_20230708_0800.csv",
  "Data/model_TOMST_result_files/TOMST_filtered_distance_temp_20230708_1200.csv"
)

tomst_times <- model_times

req_height <- 1  # (m)
length_transect <- 135 # (m)

#########
# DATA  #
#########

# Model data (Seasonally calibrated fit) in long-form
model_df <- bind_rows(
  lapply(seq_along(res_list), function(i) {
    df   <- res_list[[i]]$micro_grid
    Tair <- res_list[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z == req_height, y == 18, x <= length_transect) %>%
      transmute(
        x           = x,
        temperature = Tair - 273.15,
        model       = "Modelled",
        variant     = "Seasonally calibrated fit",
        time        = model_times[i]
      )
  })
)

# Model data (Annually calibrated fit) in long-form
model_df_cal <- bind_rows(
  lapply(seq_along(res_list_cal), function(i) {
    df   <- res_list_cal[[i]]$micro_grid
    Tair <- res_list_cal[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z == req_height, y == 18, x <= length_transect) %>%
      transmute(
        x           = x,
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
        x           = 135 - D_edge,
        temperature = Tair,
        model       = "TOMST observations",
        time        = tomst_times[i]
      )
  })
)

# Combine
anim_df <- bind_rows(model_df, model_df_cal, tomst_df)

anim_df <- anim_df %>%
  mutate(
    variant = ifelse(model == "TOMST observations", "TOMST observations", variant)
  )


anim_df <- anim_df %>%
  mutate(time = floor_date(time, unit = "hour"),
         time_lbl = format(time, "%Y-%m-%d %Hh"))

label_df <- anim_df %>%
  distinct(time_lbl) %>%
  mutate(
    x = -Inf,
    y =  Inf
  )

anim_df$time_lbl <- factor(
  anim_df$time_lbl,
  levels = sort(unique(anim_df$time_lbl)),
  labels = c("01:00 UTC", "08:00 UTC", "12:00 UTC")
)

####################
# BACKGROUND IMAGE #
####################

bg_img  <- png::readPNG("Output/horizontal_transect.png")
bg_grob <- grid::rasterGrob(bg_img, width = unit(1, "npc"), height = unit(1, "npc"))

############
#   PLOT   #
############

p_all <- ggplot(anim_df, aes(
  x = x, y = temperature,
  colour = time_lbl,
  group = interaction(model, variant, time_lbl)
)) +
  # background
  annotation_custom(
    grob = bg_grob,
    xmin = 0, xmax = length_transect,
    ymin = 15, ymax = 32
  ) +

  # Model-lines (uncal + cal), with linetype legend for model variant
  geom_line(
    data = dplyr::filter(anim_df, model == "Modelled"),
    aes(linetype = variant),
    linewidth = 1
  ) +

  # TOMST errorbars (±0.5 °C)
  geom_errorbar(
    data = dplyr::filter(anim_df, model == "TOMST observations"),
    inherit.aes = FALSE,
    aes(x = x,
        ymin = temperature - 0.5,
        ymax = temperature + 0.5,
        colour = time_lbl,
        group = interaction(model, time_lbl)),
    width = 3,
    linewidth = 0.6,
    alpha = 0.9,
    show.legend = FALSE
  ) +

  # TOMST points (shape legend for TOMST observations)
  geom_point(
    data = dplyr::filter(anim_df, model == "TOMST observations"),
    aes(shape = model),
    size = 2.5
  ) +

  # scales and legends
  scale_colour_manual(
    name = "Time",
    values = c("01:00 UTC" = "blue",
               "08:00 UTC" = "orange",
               "12:00 UTC" = "red")
  ) +
  scale_linetype_manual(
    name = "Model output",
    values = c("Seasonally calibrated fit" = "solid",
               "Annually calibrated fit"   = "dashed")
  ) +
  scale_shape_manual(
    name = "Observations",
    values = c("TOMST observations" = 16)
  ) +
  guides(
    colour   = guide_legend(order = 1, title = "Time: ", direction = "horizontal", override.aes = list(size = 3)),
    linetype = guide_legend(order = 2, title = "Model: ", direction = "horizontal", keywidth = unit (2.5, "cm"), override.aes = list(linewidth = 2)),
    shape    = guide_legend(order = 3, title = "Data: ", direction = "horizontal", override.aes = list(size = 3))
  ) +
  coord_cartesian(xlim = c(0, length_transect), ylim = c(15, 32)) +
  labs(x = "Distance from forest core (m)\n", y = "Temperature (°C)", title = "(a) Horizontal air temperature") +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(size = 35),
    legend.position   = "bottom",
    legend.box        = "vertical",
    axis.title      = element_text(size = 35),
    axis.text       = element_text(size = 35),
    panel.background  = element_rect(fill = NA, colour = NA),
    legend.title      = element_text(size = 35),
    legend.text       = element_text(size = 35),
    legend.box.margin = margin(10, 10, 10, 10),
    legend.box.background = element_rect(
      fill      = NA,
      colour    = "black",
      linewidth = 1
    ),
    plot.margin       = margin(10, 20, 10, 10)
  )


ggsave("Output/Tair_horizontal.png", p_all, width = 16, height = 12, dpi = 300)
