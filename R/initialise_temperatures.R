#' Function to define the initial forest core temperature based on macro temperature, season and time of day
#'
#' @param T_macro Macrotemperature
#' @return core_temp Forest core temperature
#' @export
#'
core_temp <- function(T_macro) {
  if (season == "summer") {
    temp_diff <- if (N_or_D == 'day') -3 else 3  # Summer day: core is cooler; night: core is warmer # initial: 3
  } else if (season == "winter") {
    temp_diff <- if (N_or_D == 'day') -1 else 1  # Winter day: core is only a bit cooler; night: core is only a bit warmer # initial: 1
  } else if (season == "spring" || season == "autumn") {
    temp_diff <- if (N_or_D == 'day') -1.5 else 1.5  # Milder gradients for spring/autumn # initial: 1.5
  }

  temperature <- T_macro + temp_diff

  return(temperature)
}

#' Sinus for near-surface air temperature: time lagged wrt macro air temperature
#'
#' @param datetime Datetime object representing the current simulation time
#' @return TG_diff Lagging temperature correction to determine soil temperature
#' @export
#'
sin_lag <- function(datetime) {
  t_now <- format(datetime, "%H")
  t_now <- as.numeric(t_now)
  lag_hours <- 6  # Hours of lagging

  TG_mean <- (macro_temp_max + macro_temp_min) / 2 # mean near-surface temp
  A <- (macro_temp_max - macro_temp_min) / 2 # amplitude of sinus
  T <- 24  # Period of 1 day

  # Macro air temp on current time point and near-surface temperature via lagging of macro temperature
  TG_now <- TG_mean + A * cos((2 * pi / T) * (t_now - t_max))
  TG_lagged <- TG_mean + A * cos((2 * pi / T) * (t_now - lag_hours - t_max))
  TG_diff <- TG_now - TG_lagged
  return(TG_diff)
}

#' Function to determine night or day
#'
#' @param datetime Datetime object representing the current simulation time
#' @return 'night' or 'day'
#' @export
#'
is_night_or_day <- function(datetime) {

  hour <- as.integer(format(datetime, "%H"))

  if (hour >= 21 || hour < 6) {
    return("night")
  } else {
    return("day")
  }
}


#' Function to determine season
#'
#' @param datetime Datetime object representing the current simulation time
#' @return 'winter' or 'spring' or 'summer' or 'autumn'
#' @export
#'
get_season <- function(datetime) {

  month <- as.integer(format(datetime, "%m"))

  if (month %in% c(12, 1, 2)) {
    return("winter")
  } else if (month %in% c(3, 4, 5)) {
    return("spring")
  } else if (month %in% c(6, 7, 8)) {
    return("summer")
  } else if (month %in% c(9, 10, 11)) {
    return("autumn")
  }
}




