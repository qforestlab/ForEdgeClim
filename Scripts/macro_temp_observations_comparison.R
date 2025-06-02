library(readr)
library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyverse)

#########
# INPUT #
#########

# Period with months for which to inspect the macrotemperature
start_month <- "2023-03"
end_month   <- "2025-02"

#########


# 3 input data sets with macrotemperature observations: RMI national weather station,
# ΔT BF5 sunshine detector and Plant Ecology temperature sensor
RMI_data <<- read.csv("Data/RMI_Melle.csv")
PE_data <<- read_delim("Data/Macro_temp_plant_eco.txt", delim = "\t", skip = 1,
  col_names = c("Date", "Time", "Temperature"), locale = locale(decimal_mark = ","),
  col_types = cols(
    Date        = col_character(),
    Time        = col_character(),
    Temperature = col_double()
  )
)
pyr_data <<-  read.csv("Data/pyranometer_tower.dat", skip = 4, header = FALSE, stringsAsFactors = FALSE)
# Add column names to pyr data
colnames(pyr_data) <- c("datetime", "record", "battery", "temp", "total", "diffuse")

# Only keep columns of interest
RMI_data = RMI_data[, c("timestamp", "temp")]
pyr_data = pyr_data[, c("datetime", "temp")]

# Set date and time to same coordinates (UTC), RMI is already in UTC
# pyr_data:
# Set datetime as a POSIXct (without forcing a time zone)
pyr_data$datetime <- ymd_hms(pyr_data$datetime)
# pyr_data is always in summer time (CEST = UTC + 2) => reduce by 2h to set into UTC
pyr_data$datetime_utc <- pyr_data$datetime - hours(2)
pyr_data = pyr_data[, c("datetime_utc", "temp")]
# PE_data:
# Merge Date and Time to POSIXct object in Brussel's time zone
PE_data <- PE_data %>%
  mutate(
    timestamp_local = dmy_hm(paste(Date, Time), tz = "Europe/Brussels"),
    timestamp_utc = with_tz(timestamp_local, tzone = "UTC")
  )
PE_data = PE_data[, c("timestamp_utc", "Temperature")]


# Make the timestamps similar and calculate hourly average

# RMI_data: has data every hour
RMI_data_clean <- RMI_data %>%
  mutate(time = ymd_hms(timestamp, tz = "UTC", quiet = TRUE),
         hour = floor_date(time, unit = "hour")) %>%
  select(hour, temp) %>%
  rename(Temperature = temp) %>%
  mutate(Source = "RMI") %>%
  arrange(hour)  # sort timestamps increasingly

# pyr_data: 15-min data -> hourly average
pyr_data_clean <- pyr_data %>%
  mutate(hour = floor_date(datetime_utc, unit = "hour")) %>%
  group_by(hour) %>%
  summarise(Temperature = mean(temp, na.rm = TRUE), .groups = "drop") %>%
  mutate(Source = "ΔT")

# PE_data: 5-min data -> hourly average
PE_data_clean <- PE_data %>%
  mutate(hour = floor_date(timestamp_utc, unit = "hour")) %>%
  group_by(hour) %>%
  summarise(Temperature = mean(Temperature, na.rm = TRUE), .groups = "drop") %>%
  mutate(Source = "PE")

# Combine datasets
combined_data <- bind_rows(RMI_data_clean, pyr_data_clean, PE_data_clean)



# Loop over months
###################

month_seq <- seq(
  ymd_hms(paste0(start_month, "-01 00:00:01"), tz = "UTC"),
  ymd_hms(paste0(end_month, "-01 0:00:01"), tz = "UTC"),
  by = "1 month"
)

for (month_start in month_seq) {
  # Get start and end of month
  month_start = as_datetime(month_start)
  month_end <- ceiling_date(as_datetime(month_start), "month") - hours(1)  # last hour of month

  # Get POSIXct series of every hour of this month
  hours_in_month <- seq(month_start, month_end, by = "1 hour")

  # Filter combined data for this month
  data_month <- combined_data %>%
    filter(hour >= month_start & hour <= month_end)

  data_month$Source <- factor(data_month$Source, levels = c("PE", "ΔT", "RMI")) # order in legend

  # Plot
  p <- ggplot(data_month, aes(x = hour, y = Temperature, color = Source)) +
    geom_line() +
    labs(
      title = paste0("Hourly averaged macrotemperature for ", format(month_start, "%B %Y")),
      x = "Date",
      y = "Temperature (°C)",
      color = "Source"
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white", color = "black")
    ) +
    scale_x_datetime(date_breaks = "5 days", date_labels = "%d-%b") +
    scale_color_manual(
      values = c(
        "RMI" = "cornflowerblue",
        "ΔT" = "orange",
        "PE"  = "darkgreen"
      )
    )

  # Save
  ggsave(
    filename = paste0("Output/macro_temps/Temperature_Comparison_", format(month_start, "%Y_%m"), ".png"),
    plot = p, width = 10, height = 6, dpi = 300
  )

  #print(p)
}


