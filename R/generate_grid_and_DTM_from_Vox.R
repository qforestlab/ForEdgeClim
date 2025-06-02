#' Generate a voxelised point density grid from a vox file
#'
#' This function reads a vox file and creates a 3D voxel grid where each voxel represents
#' a 1x1x1m^3 volume. The density of each voxel is computed using Plant Area Density (PAD) from AMAPVox.
#' The output variable in AMAPVox representing PAD is 'Attenuation unbiased maximum likelihood estimator based
#' on free path length (FPL MLE, F. Pimont)'. It takes into account beam attenuation and corrects for
#' the fact that not all beams can 'see' everything (occlustion).
#' The output is a data.table containing rescaled X, Y, and Z coordinates and a normalized density value.
#'
#' Additionally, a DTM is generated.
#'
#' @param vox_file The vox file
#' @return A data.table with columns: X, Y, Z, and density.
#' @importFrom AMAPVox readVoxelSpace
#' @importFrom data.table as.data.table setDT
#' @importFrom sp SpatialPointsDataFrame CRS
#' @importFrom lidR LAS grid_terrain normalize_height classify_ground csf
#' @importFrom raster rasterFromXYZ disaggregate rasterToPoints extract as.data.frame
#' @export
generate_DTM_grid_Vox <- function(vox_file) {

  #############
  # STRUCTURE #
  #############

  # Read voxel space
  vx <- readVoxelSpace(vox_file)
  vx_dt <- as.data.table(slot(vx, "data"))

  vx_dt$pad = vx_dt$attenuation_FPL_unbiasedMLE        # Plant Area Density
  # normalise pad with percentile clipping
  lower_bound <- quantile(vx_dt$pad, 0.01, na.rm = TRUE)
  upper_bound <- quantile(vx_dt$pad, 0.95, na.rm = TRUE)

  vx_dt$pad <- pmin(pmax((vx_dt$pad - lower_bound) / (upper_bound - lower_bound), 0), 1)

  vx_dt[is.na(pad), pad := 0]

  setnames(vx_dt, old = c("i", "j", "k"), new = c("X", "Y", "Z"))


  #######
  # DTM #
  #######

  point_df = vx_dt[vx_dt$pad > 0, c("X", "Y", "Z")]
  point_df$X <- as.numeric(as.character(point_df$X))
  point_df$Y <- as.numeric(as.character(point_df$Y))
  point_df$Z <- as.numeric(as.character(point_df$Z))


  # 4. Zet om naar een LAS-object
  las <- LAS(point_df)

  # 5. Exporteer als LAS-bestand
  writeLAS(las, "Data/vox_to_las.las")


  las <- readLAS("Data/vox_to_las.las")

  # Detecteer grondpunten (met bijv. PMF-algoritme)
  las <- classify_ground(las, pmf(ws = c(3, 5, 7, 9), th = c(1, 1.5, 2, 2.5)))

  # Interpoleer een DTM
  dtm <- grid_terrain(las, res = 1, algorithm = knnidw())

  # Omzetten naar dataframe
  dtm_df <- raster::as.data.frame(dtm, xy = TRUE)
  setnames(dtm_df, old = c("x", "y"), new = c("X", "Y"))


  #print(dtm_df)

  # # Get for every (x, y) combination the smallest Z-value where pad > 0
  # dtm <- vx_dt[pad > 0, .(ground_z = min(Z)), by = .(X, Y)]
  # dtm[, ground_z := ground_z - min(ground_z, na.rm = TRUE)]
  #
  # setnames(dtm, old = c("ground_z"), new = c("Z"))
  #
  #
  # ##############################
  # # CORRECT STRUCTURE WITH DTM #
  # ##############################



  # Add dtm to vx_dt based on X and Y
  vx_dt <- merge(vx_dt, dtm_df, by = c("X", "Y"), all.x = TRUE, suffixes = c("", "_ground"))

  # Correct Z so ground matches Z = 0
  vx_dt[, Z_corrected := Z - Z_ground]

  # # Do not keep Z lower than 0
  vx_dt <- vx_dt[Z_corrected >= 0]
  #
  # Only keep columns of interest
  vx_dt <- vx_dt[, .(X, Y, Z_corrected, pad)]
  setnames(vx_dt, old = c("Z_corrected", "pad"), new = c("Z", 'density'))

  # rescale vx_dt
  vx_dt <- vx_dt[X > 15 & X <= (170 - 5)]
  vx_dt[, X := X - min(X) + 1]
  vx_dt <- vx_dt[Y > 5 & Y <= (40 - 5)]
  vx_dt[, Y := Y - min(Y) + 1]
  vx_dt <- vx_dt[Z <= 41]
  vx_dt[, Z := Z - min(Z) + 1]
  #
  #
  # print(summary(vx_dt))

  ###################
  # FUNCTION OUTPUT #
  ###################

  print(summary(dtm_df))
  print(summary(vx_dt))

  # Return exact same structure as in TLS function
  return(list(
    dtm  = dtm_df,   # X, Y, Z ground elevation
    grid = vx_dt   # X, Y, Z height-above-ground, density
  ))

}



# # coordinate boundaries
# xmin <- -135.002; ymin <- 15.0005; zmin <- -15.001
# xmax <-  15.0015; ymax <- 59.9015; zmax <- 103.214

