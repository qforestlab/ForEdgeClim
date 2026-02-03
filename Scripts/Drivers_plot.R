############################################################################################
# This script visualizes the macroclimatic input drivers (macrotemperature and
# shortwave radiative values) used for the calibration process of ForEdgeClim.
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################


# ---- Packages ----
library(tidyverse)
library(lubridate)

# ---- Files (paths) ----
f_rmi <- "Data/RMI_Melle.csv"
f_pyr <- "Data/pyranometer_tower.dat"

season = "summer"


# ---- Choose 3 days to plot (UTC dates; format YYYY-MM-DD) ----
if (season == "winter"){
  days_to_plot <- as.Date(c("2025-01-13", "2025-01-20", "2025-01-23"))

} else if (season == "spring"){
  days_to_plot <- as.Date(c("2025-04-30", "2025-04-23", "2024-04-02"))

}else if (season == "summer"){
  days_to_plot <- as.Date(c("2023-07-07", "2023-07-31", "2023-07-20"))

}else { # autumn
  days_to_plot <- as.Date(c("2023-10-01", "2024-10-09", "2023-10-19"))

}

# =============================================================================
# 1) RMI_Melle (timestamps are UTC)
#    Robust parsing without parse_date_time()
# =============================================================================
rmi_raw <- readr::read_csv(f_rmi, show_col_types = FALSE)

rmi <- rmi_raw %>%
  mutate(
    timestamp_chr = as.character(timestamp),
    # try normal "YYYY-mm-dd HH:MM:SS"
    ts1 = ymd_hms(timestamp_chr, tz = "UTC", quiet = TRUE),
    # try ISO "YYYY-mm-ddTHH:MM:SS"
    ts2 = ymd_hms(str_replace(timestamp_chr, "T", " "), tz = "UTC", quiet = TRUE),
    timestamp = coalesce(ts1, ts2),
    date = as.Date(timestamp)
  ) %>%
  select(timestamp, date, temp)

# optional: see how many failed
# sum(is.na(rmi$timestamp))

rmi_hourly <- rmi %>%
  dplyr::filter(date %in% days_to_plot) %>%
  mutate(hour = floor_date(timestamp, unit = "hour")) %>%
  group_by(date, hour) %>%
  summarise(macro_temp = mean(temp, na.rm = TRUE), .groups = "drop")

# =============================================================================
# 2) Pyranometer tower (TIMESTAMP is summer time UTC+2 -> convert to UTC)
# =============================================================================
pyr <- readr::read_csv(
  f_pyr,
  skip = 4,
  col_names = c("TIMESTAMP", "RECORD", "BattV_Avg", "PTemp_C_Avg", "Total", "Diffuse"),
  show_col_types = FALSE,
  na = c("NAN", "NaN", "")
) %>%
  mutate(
    TIMESTAMP = ymd_hms(TIMESTAMP, tz = "Europe/Brussels", quiet = TRUE),
    TIMESTAMP = with_tz(TIMESTAMP, "UTC"),
    date = as.Date(TIMESTAMP),
    direct = Total - Diffuse
  ) %>%
  select(TIMESTAMP, date, Total, Diffuse, direct)

pyr_hourly <- pyr %>%
  dplyr::filter(date %in% days_to_plot) %>%
  mutate(hour = floor_date(TIMESTAMP, unit = "hour")) %>%
  group_by(date, hour) %>%
  summarise(
    diffuse_sw = mean(Diffuse, na.rm = TRUE)/2,
    direct_sw  = mean(direct,  na.rm = TRUE)/2,
    .groups = "drop"
  )
# =============================================================================
# 3) Combine + reshape for 2-panel plot (separate y-scales)
# =============================================================================
wide_df <- full_join(rmi_hourly, pyr_hourly, by = c("date", "hour")) %>%
  mutate(
    hour_of_day = hour(hour),
    day_type = factor(
      as.character(date),
      levels = as.character(days_to_plot),
      labels = c("Most sunny day", "Most cloudy day", "Most solar fluctuating day")
    )
  )

plot_long <- bind_rows(
  # Macrotemperature panel
  wide_df %>%
    transmute(
      hour_of_day, day_type,
      panel = "(a) Macrotemperature (°C)",
      series = "Macrotemperature",
      value = macro_temp
    ),
  # Radiation panel
  wide_df %>%
    transmute(
      hour_of_day, day_type,
      panel = "(b) Shortwave radiation (W m\u207B\u00B2)",  # W m⁻²
      series = "Direct-beam",
      value = direct_sw
    ),
  wide_df %>%
    transmute(
      hour_of_day, day_type,
      panel = "(b) Shortwave radiation (W m\u207B\u00B2)",
      series = "Diffuse",
      value = diffuse_sw
    )
) %>%
  mutate(
    panel  = factor(panel, levels = c("(a) Macrotemperature (°C)", "(b) Shortwave radiation (W m\u207B\u00B2)")),
    series = factor(series, levels = c("Macrotemperature", "Direct-beam", "Diffuse"))
  )

# =============================================================================
# 4) Plot: facets on top, clean legends, in-panel radiation legend with box
# =============================================================================

# ---- Mini-legend in radiation facet (rounded box via geom_label) ----
rad_panel_name <- "(b) Shortwave radiation (W m\u207B\u00B2)"
rad_max <- max(plot_long$value[plot_long$panel == rad_panel_name], na.rm = TRUE)
rad_min <- min(plot_long$value[plot_long$panel == rad_panel_name], na.rm = TRUE)
rad_rng <- rad_max - rad_min

# Anchor (top-left-ish of the legend box)
leg_x <- 16.4
leg_y <- rad_max - 0.030 * rad_rng

# Create a "blank" label that still has width/height:
# - use non-breaking spaces (\u00A0) so spaces aren't trimmed
# - add multiple lines to increase height
# Create a "blank" label that controls WIDTH only (no extra height)
blank_line <- paste(rep("\u00A0", 32), collapse = "")
blank_lab  <- paste(blank_line, blank_line, sep = "\n")

box_df <- tibble(
  panel = rad_panel_name,
  x = leg_x,
  y = leg_y,
  lab = blank_lab
)


# Row y positions (more spacing between text lines)
y_direct  <- leg_y - 0.040 * rad_rng
y_diffuse <- leg_y - 0.095 * rad_rng

# Line samples (aligned with text)
seg_df <- tibble(
  panel = rad_panel_name,
  x    = leg_x + 0.40,
  xend = leg_x + 1.90,
  y    = c(y_direct, y_diffuse),
  yend = c(y_direct, y_diffuse),
  lty  = c("solid", "dashed")
)

# Text (a bit further from the lines)
txt_df <- tibble(
  panel = rad_panel_name,
  x = leg_x + 2.20,
  y = c(y_direct, y_diffuse),
  lab = c("Direct-beam", "Diffuse")
)


# ---- Plot ----
p <- ggplot(
  plot_long,
  aes(
    x = hour_of_day,
    y = value,
    colour = day_type,
    linetype = series,
    group = interaction(day_type, series)
  )
) +
  geom_line(linewidth = 1.0, na.rm = TRUE) +

  scale_x_continuous(breaks = seq(0, 23, 3), limits = c(0, 23)) +

  facet_wrap(~ panel, nrow = 1, scales = "free_y", strip.position = "top") +

  scale_colour_manual(
    values = c(
      "Most sunny day"             = "red",
      "Most cloudy day"            = "blue",
      "Most solar fluctuating day" = "orange"
    )
  ) +

  scale_linetype_manual(
    values = c(
      "Macrotemperature" = "solid",
      "Direct-beam"      = "solid",
      "Diffuse"          = "dashed",
      "solid"            = "solid",
      "dashed"           = "dashed"
    )
  ) +

  labs(
    x = "Hour of day (UTC)",
    y = NULL,
    colour = "Day type",
    title = paste0(
      "Mean hourly drivers for selected calibration days during the ",
      season, " season"
    )
  ) +

  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_blank(), #element_text(size = 15),
    strip.text = element_text(size = 15),
    strip.text.x = element_text(size = 15, hjust = 0),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    strip.placement = "outside",
    legend.position = "bottom",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +

  guides(
    linetype = "none",
    colour = guide_legend(
      order = 1,
      override.aes = list(linetype = "solid")
    )
  ) +

  # Rounded mini-legend box (radiation facet only)
  geom_label(
    data = box_df,
    aes(x = x, y = y, label = lab),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    label.size = 0.4,
    label.r = unit(0.35, "lines"),
    label.padding = unit(0.35, "lines"),
    fill = "white", colour = "black",
    size = 3.2
  ) +


  # Line samples
  geom_segment(
    data = seg_df,
    aes(x = x, xend = xend, y = y, yend = yend, linetype = lty),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.9
  ) +

  # Text
  geom_text(
    data = txt_df,
    aes(x = x, y = y, label = lab),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0.5,
    size = 3.4,
    colour = "black"
  ) +


  coord_cartesian(clip = "off")

print(p)


ggsave(
  "Output/drivers.png",
  plot = p,
  width = 10,
  height = 4.5,
  dpi = 300
)
