###############################################################################
#
# In this script, a time series of model and observations is plotted for 6 time
# points on July 8 2023. For both model and observations, we compare edge vs
# core positions within the forest transect.
#
###############################################################################


library(dplyr)
library(lubridate)
library(ggplot2)
library(readxl)
library(purrr)
library(readr)
library(mgcv)   # voor s(..., bs = "cc")


#########
# INPUT #
#########

## Define hours of the day as character with leading zeros
hours     <- 0:23
hour_str  <- sprintf("%02d", hours)

model_times <- ymd_h(paste("2023-07-08", hour_str))
tomst_times <- model_times

# model_results files
res_list <- map(
  hour_str,
  ~ readRDS(
    file.path("Output", paste0("2023-07-08_", .x), "model_results.rds")
  )
)

# TOMST files
tomst_files <- file.path(
  "Data/model_TOMST_result_files",
  paste0("TOMST_filtered_distance_temp_20230708_", hour_str, "00.csv")
)

# model_times <- lubridate::ymd_h(c("2023-07-08 00", "2023-07-08 04", "2023-07-08 08",
#                                   "2023-07-08 12", "2023-07-08 16", "2023-07-08 20"))
#
# res_list   <- list(
#   readRDS("Data/model_result_files/model_results_0h.rds"),
#   readRDS("Data/model_result_files/model_results_4h.rds"),
#   readRDS("Data/model_result_files/model_results_8h.rds"),
#   readRDS("Data/model_result_files/model_results_12h.rds"),
#   readRDS("Data/model_result_files/model_results_16h.rds"),
#   readRDS("Data/model_result_files/model_results_20h.rds")
# )
#
# tomst_files <- c(
#   "Data/model_result_files/TOMST_filtered_distance_temp_0h.csv",
#   "Data/model_result_files/TOMST_filtered_distance_temp_4h.csv",
#   "Data/model_result_files/TOMST_filtered_distance_temp_8h.csv",
#   "Data/model_result_files/TOMST_filtered_distance_temp_12h.csv",
#   "Data/model_result_files/TOMST_filtered_distance_temp_16h.csv",
#   "Data/model_result_files/TOMST_filtered_distance_temp_20h.csv"
# )
#
# tomst_times <- model_times

req_height = 1 # (m)
length_transect = 135 # (m)


####################################
# DATA FROM MODEL AND OBSERVATIONS #
####################################

# Model core data in long-form
model_core_df <- bind_rows(
  lapply(seq_along(res_list), function(i) {
    df   <- res_list[[i]]$micro_grid
    Tair <- res_list[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z == req_height, y == 18, x <= min(x)) %>%
      transmute(
        x           = x,
        temperature = Tair - 273.15,
        id       = "Model-Core",
        Source       = "Modelled",
        Position = "Core",
        time        = model_times[i]
      )
  })
)

# Model edge data in long-form
model_edge_df <- bind_rows(
  lapply(seq_along(res_list), function(i) {
    df   <- res_list[[i]]$micro_grid
    Tair <- res_list[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z == req_height, y == 18, x == length_transect) %>%
      transmute(
        x           = x,
        temperature = Tair - 273.15,
        id       = "Model-Edge",
        Source       = "Modelled",
        Position = "Edge",
        time        = model_times[i]
      )
  })
)

# Macro data in long-form
macro_df <- bind_rows(
  lapply(seq_along(res_list), function(i) {
    df   <- res_list[[i]]$micro_grid
    Tair <- res_list[[i]]$air_temperature
    df %>%
      mutate(Tair = Tair) %>%
      filter(z == max(z), y == 18, x == max(x)) %>%
      transmute(
        x           = x,
        temperature = Tair - 273.15,
        id       = "Macro",
        Source       = "Modelled", # macro modelled = macro observed temperature (model input variable)
        Position = "Macro",
        time        = model_times[i]
      )
  })
)

# TOMST core data in long-form
tomst_core_df <- bind_rows(
  lapply(seq_along(tomst_files), function(i) {
    read.csv(tomst_files[i]) %>%
      transmute(
        x           = 135 - D_edge,
        temperature = Tair,
        id       = "TOMST-Core",
        Source       = "Observed",
        Position = "Core",
        time        = tomst_times[i]
      ) %>%
      filter(x == 0)
  })
)

# TOMST edge data in long-form
tomst_edge_df <- bind_rows(
  lapply(seq_along(tomst_files), function(i) {
    read.csv(tomst_files[i]) %>%
      transmute(
        x           = 135 - D_edge,
        temperature = Tair,
        id       = "TOMST-Edge",
        Source       = "Observed",
        Position = "Edge",
        time        = tomst_times[i]
      ) %>%
      filter(x == 135)
  })
)


# combine dataframes
combined_df <- bind_rows(model_core_df, model_edge_df, tomst_core_df, tomst_edge_df, macro_df)
combined_df$time <- as.POSIXct(combined_df$time, format = "%Y-%m-%d %H:%M:%S")
combined_df$Position <- factor(combined_df$Position, levels = c("Macro", "Edge", "Core"))

combined_df <- combined_df %>%
  mutate(hour = as.numeric(format(time, "%H")))



##############################################################
# PLOTTING TIMESERIES CORE VS EDGE AND MODEL VS OBSERVATIONS #
##############################################################



timeseries <- ggplot(combined_df, aes(x = hour, y = temperature)) +
  # points for TOMST observations
  geom_point(
    data = combined_df %>% filter(Source == "Observed"),
    aes(color = Position, shape = "Observed"),
    size = 3
  ) +
  # error bars for TOMST observations (± 0.5 °C)
  geom_errorbar(
    data = combined_df %>% filter(Source == "Observed"),
    aes(
      color = Position,
      ymin  = temperature - 0.5,
      ymax  = temperature + 0.5
    ),
    width = 0.2,
    linewidth = 0.6
  ) +
  # periodic smoothing line for modelled edge and core temperature
  geom_smooth(
    data = combined_df %>% filter(Source == "Modelled", Position != "Macro"),
    aes(color = Position, linetype = "Modelled"),
    method      = "gam",
    formula     = y ~ s(x, bs = "cc", k = 8),   # cyclic cubic spline
    method.args = list(knots = list(x = c(0, 24))),
    se          = FALSE,
    size = 1
  ) +
  # periodic smoothing line for macro temperature
  geom_smooth(
    data = combined_df %>% filter(Source == "Modelled", Position == "Macro"),
    aes(color = Position, linetype = "Modelled"),
    method      = "gam",
    formula     = y ~ s(x, bs = "cc", k = 8),   # cyclic cubic spline
    method.args = list(knots = list(x = c(0, 24))),
    se          = FALSE,
    size = 2.5
  ) +
  labs(
    x = "Hour (UTC)",
    y = "Temperature (°C)",
    color   = "Position",
    shape    = "Source",
    linetype = "Source"
  ) +
  scale_color_manual(
    values = c(
      "Macro" = "black",
      "Core"  = "blue",
      "Edge"  = "red"
    )
  ) +
  scale_shape_manual(
    values = c("Observed" = 16)
  ) +
  scale_linetype_manual(values = c(
    "Modelled" = "solid")) +
  scale_x_continuous(
    breaks = seq(0, 23, by = 4),
    limits = c(0, 23),
    expand = c(0.01, 0.01)
  ) +
  guides(
    color = guide_legend(order = 1, direction = "horizontal", title = "Position:", override.aes = list(
      size = 6,
      linewidth = 3
    )),
    linetype = guide_legend(title = "|  Source: ", order = 2, direction = "horizontal",
      override.aes = list(linewidth = 1,color = "black")
    ),
    shape = guide_legend(order = 3, direction = "horizontal", title = NULL, override.aes = list(size = 4))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 35, face = "bold"),
    plot.margin = margin(0, 25, 0, 25),
    axis.title = element_text(size = 35),
    axis.text  = element_text(size = 35),
    legend.position = "top",
    legend.box = "horizontal",
    legend.title = element_text(size = 35),
    legend.text  = element_text(size = 35),
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA)
  )

print(timeseries)

ggsave('Output/timeseries_core_vs_edge_AND_model_vs_observations.png',
       plot = timeseries, width = 16, height = 10, dpi = 300)
