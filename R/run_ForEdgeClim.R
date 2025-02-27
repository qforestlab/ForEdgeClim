###############################################################################
#
#                               ForEdgeClim
#
# A microclimate model that solves the energy balance equation in fragmented
# forest areas: Rn = LE + H, where,
# Rn = netto radiation (shortwave & longwave)
# LE = latent heat flux (~ evapotranspiration)
# H = sensible heat flux
###############################################################################


#' Execute the microclimate model ForEdgeClim
#'
#' @param input_data Dataframe with input variables
#' @return Dataframe with simulated microclimate temperatures and fluxes
#' @export
run_foredgeclim <- function() {


  #########
  # INPUT #
  #########

  create_input_data(summer_day = TRUE)

  ###############################################################################

  print('🚩 𝙁𝙤𝙧𝙀𝙙𝙜𝙚𝘾𝙡𝙞𝙢')
  start = Sys.time()
  print('════════════════════════')

  ##################################
  # INITIALISATION OF TIME PERIODS #
  ##################################

  N_or_D <<- is_night_or_day(datetime)

  season <<- get_season(datetime)

  #######################
  # Creation of 3D grid #
  #######################

  print('3D grid 🌲🌳🌲🌳🌲')
  voxel_grid = generate_virtual_grid(season)

  #############################################
  # 2D SHORTWAVE RADIATIVE TRANSFER MODELLING #
  #############################################

  print('SW RTM ☀️')
  final_results_2D = shortwave_two_stream_RTM(datetime, lat, lon, voxel_grid,
                           F_sky_diff_init, F_sky_dir_init,
                           omega_g_v, omega_g_h,
                           Kd_v, Kb_v, Kd_h, Kb_h,
                           omega, betad, beta0) # A lot of these parameters are global parameters and actually do not need to be stated as arguments.


  den <<- final_results_2D$density
  net_SW <- final_results_2D$net_sw
  SW_down <- final_results_2D$F_d_down + final_results_2D$F_b_down
  SW_up <- final_results_2D$F_d_up
  x_coords <- final_results_2D$X
  y_coords <- final_results_2D$Y
  z_coords <- final_results_2D$Z

  # dimensions grid
  x_dim <<- max(x_coords)  # Length in x-direction
  y_dim <<- max(y_coords)  # Length in y-direction
  z_dim <<- max(z_coords)  # Length in z-direction

  #############################
  # INITIALISING TEMPERATURES #
  #############################

  print('Init temp 🌡️')

  core_temperature <<- core_temp(macro_temp)
  TG_diff = sin_lag(datetime)

  # First longwave RTM for initialising temperatures
  micro_grid <- data.frame(
    x = x_coords,
    y = y_coords,
    z = z_coords,
    temperature = macro_temp,
    T_soil = macro_temp
  )

  print('LW RTM ⛅ ')
  final_results_lw_2D = longwave_two_stream_RTM(voxel_grid, micro_grid, lw_two_stream,
                                                F_sky_lw, omega_g_lw_v, omega_g_lw_h, macro_temp,
                                                Kd_lw_v, Kd_lw_h, omega_lw, beta_lw, season, N_or_D)

  net_LW <- final_results_lw_2D$net_lw
  LW_down <- final_results_lw_2D$F_d_down
  LW_up <- final_results_lw_2D$F_d_up
  net_rad_ground <- ifelse(z_coords == 1, SW_down - SW_up + LW_down - LW_up, NA)

  # Net radiation forest surface
  net_radiation_init <- net_SW + net_LW

  # Define ground flux as a % of net rad at ground
  G = calculate_G(net_rad_ground)

  # Define T_soil from ground flux with ground conductance/convection g_soil
  # For this, first make a rudimentary guess on T air, weighted sigmoid weighting of the boundary temperatures
  beta = 0.01
  dis_macro_x = x_dim - x_coords+1
  dis_macro_z = z_dim - z_coords+1
  dis_core_x = x_coords
  dis_core_z = z_coords
  dis_mean = mean(c(dis_macro_z,dis_macro_x,dis_core_z,dis_core_x))
  w_macro_x = 1 / ( 1 + exp(beta*(dis_macro_x - dis_mean) ))
  w_macro_z = 1 / ( 1 + exp(beta*(dis_macro_z - dis_mean) ))
  w_core_x = 1 / ( 1 + exp(beta*(dis_core_x - dis_mean) ))
  w_core_z = 1 / ( 1 + exp(beta*(dis_core_z - dis_mean) ))
  T_air_vec = (macro_temp*(w_macro_x + w_macro_z) + core_temperature*(w_core_x+w_core_z))/(w_macro_x+w_macro_z+w_core_x+w_core_z)
  micro_grid$T_soil = G/g_soil + T_air_vec - TG_diff

  # Define micro_grid$temperature based on macro temp and initial net radiation + random noise
  micro_grid$temperature =  macro_temp + 0.01*net_radiation_init + rnorm(nrow(micro_grid), mean = 0, sd = 0.1)

  # Air temperature update based on linearisation as is done in microclimc (Maclean & Klinges, 2021)
  # This linearisation simplifies the complex nonlinear interactions between leaves and air, making the calculations more efficient.
  # This linearisation is mostly valid when there are small temperature differences between T air en T leaf, when there is
  # sufficient air streaming, when the net radiation or the humidity is not extreme, when the forest structure is quite homogeneous.
  # The linearisation calculates the air temperature as influenced by the macro, soil, core temperature boundaries and the surface temperature.
  T_air_vec = micro_grid$temperature + rnorm(nrow(micro_grid), mean = 0, sd = 0.1)
  # Define T_ground for every voxel in a vertical xy-column as T_soil
  micro_grid <- micro_grid |>
    dplyr::group_by(x, y) |>
    dplyr::mutate(T_ground = T_soil[z == 1][1]) |>
    dplyr::ungroup()
  # Empirical formulae for convective heat conductance coefficient (conductance) (W/m2/K):
  # gh_soil = 1.4* ( abs(micro_grid$T_ground - T_air_vec) )^(1/3)
  # gh_macro = 1.4 * ( abs(macro_temp - T_air_vec) )^(1/3)
  # gh_forest = 1.4 * ( abs(micro_grid$temperature - T_air_vec) )^(1/3)

  # Determine weights for each component in the linearisation via sigmoid weighting
  dis_core = x_coords
  dis_soil = z_coords
  dis_forest = 0.5 # size of a voxel edge = 1m, so distance from structure is on average 0.5m
  dis_mean = mean(c(dis_macro_x,dis_macro_z,dis_core,dis_soil,dis_forest))
  w_core = 1 / ( 1 + exp(beta*(dis_core - dis_mean) ))
  w_soil = 1 / ( 1 + exp(beta*(dis_soil - dis_mean) ))
  w_forest = 1 / ( 1 + exp(beta*(dis_forest - dis_mean) ))

  # Update air temperature, weighted by the conductance/convection coefficient and the distance to each boundary
  T_air_vec = ( (w_macro_x + w_macro_z)*g_macro * macro_temp + w_core*g_macro*core_temperature + w_soil*g_soil * micro_grid$T_ground + w_forest*g_forest * micro_grid$temperature ) / ( (w_macro_x + w_macro_z + w_core)*g_macro + w_soil*g_soil + w_forest*g_forest)


  #########################################
  # ITERATION TO CLOSE THE ENERGY BALANCE #
  #########################################

  print('E bal ⚖️')

  # Start iteration
  max_iter = 1000
  W = 1 # weighting step

  for (iter in 1:max_iter) {

    print(paste0('-> Iteration ', iter))


    ###################
    # NETTO RADIATION #
    ###################

    print('LW RTM ⛅ ')
    final_results_lw_2D = longwave_two_stream_RTM(voxel_grid, micro_grid, lw_two_stream,
                                                  F_sky_lw, omega_g_lw_v, omega_g_lw_h, macro_temp,
                                                  Kd_lw_v, Kd_lw_h, omega_lw, beta_lw, season, N_or_D) # A lot of these parameters are global parameters and actually do not need to be stated as arguments.


    net_LW = final_results_lw_2D$net_lw

    LW_down = final_results_lw_2D$F_d_down
    LW_up = final_results_lw_2D$F_d_up

    net_rad_ground = ifelse(z_coords == 1, SW_down - SW_up + LW_down - LW_up, NA)

    # Netto radiation for structure
    net_radiation = net_SW + net_LW


    ###############
    # Ground flux #
    ###############

    print('G 🟫')

    # G is positive if E is entering the soil
    G = calculate_G(net_rad_ground)
    micro_grid$T_soil = G/g_soil + T_air_vec - TG_diff
    micro_grid <- micro_grid |>
      dplyr::group_by(x, y) |>
      dplyr::mutate(T_ground = T_soil[z == 1][1]) |>
      dplyr::ungroup()


    ##############################
    # AIR TO AIR HEAT CONVECTION #
    ##############################

    print('C 💨')
    dt = 1 # convection step is 1, ie, the iteration step taken for temp convergence
    T_air_vec = T_air_vec - calculate_C(T_air_vec,micro_grid$T_ground)*dt/(Cp*V*rho)


    ######################
    # SENSIBLE HEAT FLUX #
    ######################

    print('H ♨️')
    sensible_flux = calculate_H(micro_grid$temperature,T_air_vec)


    ####################
    # LATENT HEAT FLUX #
    ####################

    print('LE 💧')
    latent_flux <- calculate_LE(micro_grid$temperature-273.15,net_radiation) # temp in calculate_LE must be in °C


    ##################
    # ENERGY BALANCE #
    ##################

    # Surface energy balance
    energy_balance_surf <- net_radiation - sensible_flux - latent_flux
    print(max(abs(energy_balance_surf)))

    # print(summary(net_radiation))
    # print(summary(sensible_flux))
    # print(summary(latent_flux))
    # print(summary(G))

    # Check convergence
    if (max(abs(energy_balance_surf)) < energy_balance_tolerance){
      print(paste0('Convergence is reached after ', iter, ' iterations'))
      break
    }

    # Use Newton's method to update surface temperature as is done in the SCOPE model
    micro_grid$temperature <- micro_grid$temperature - ifelse(energy_balance_surf == 0, 0, W*energy_balance_surf/de_dT(micro_grid$temperature,net_radiation) )

    micro_grid <- micro_grid |>
      dplyr::group_by(z) |>
      dplyr::mutate(mean_temp_z = mean(temperature[den>0], na.rm = TRUE)) |>
      dplyr::ungroup()
    micro_grid <- micro_grid |>
      dplyr::group_by(y) |>
      dplyr::mutate(mean_temp_y = mean(temperature[den>0], na.rm = TRUE)) |>
      dplyr::ungroup()
    micro_grid <- micro_grid |>
      dplyr::group_by(x) |>
      dplyr::mutate(mean_temp_x = mean(temperature[den>0], na.rm = TRUE)) |>
      dplyr::ungroup()
    micro_grid$mean_temp = (micro_grid$mean_temp_x + micro_grid$mean_temp_y + micro_grid$mean_temp_z)/3


    # Air temperature update based on linearisation as is done in microclimc (Maclean & Klinges, 2021)
    # As compared with the initial T_air_vec, the core temperature influence is not accounted for anymore.
    #T_air_vec = ( (w_macro_x + w_macro_z)*g_macro * macro_temp  + w_soil*g_soil * micro_grid$T_ground + ifelse(den == 0, w_forest*g_forest * micro_grid$mean_temp_z, w_forest*g_forest *micro_grid$temperature) ) / ( (w_macro_x + w_macro_z)*g_macro + w_soil*g_soil + w_forest*g_forest)
    T_air_vec <- ( (w_macro_x + w_macro_z)*g_macro * macro_temp + w_core*g_macro*core_temperature+ w_soil*g_soil * micro_grid$T_ground + ifelse(den == 0, w_forest*g_forest * micro_grid$mean_temp_z, w_forest*g_forest *micro_grid$temperature) )/ ( (w_macro_x + w_macro_z + w_core)*g_macro + w_soil*g_soil + w_forest*g_forest)

    # print('T_air:')
    # print(summary(T_air_vec-273.15))
    # print('T_forest:')
    # print(summary(micro_grid$temperature-273.15))

    #W <- min(1, max(0.1, 0.5 * energy_balance_tolerance / max(abs(energy_balance_surf))))
  }

  res <- list(
    micro_grid = micro_grid,
    air_temperature = T_air_vec,
    net_radiation = net_radiation,
    sensible_flux = sensible_flux,
    latent_flux = latent_flux,
    ground_flux = G
  )

  # Results
  cat("Final surface temperature distribution (°C):\n")
  print(summary(micro_grid$temperature - 273.15))
  cat("Final air temperature distribution (°C):\n")
  print(summary(T_air_vec- 273.15))
  cat("Final soil temperature distribution (°C):\n")
  print(summary(micro_grid$T_soil- 273.15))


  finish = Sys.time()
  print(paste0('Run time = ', round(as.numeric(finish - start, units = "secs"), 2), ' s'))
  print('════════════════════════')
  print('🏁 𝙁𝙤𝙧𝙀𝙙𝙜𝙚𝘾𝙡𝙞𝙢')

  return(res)

}



