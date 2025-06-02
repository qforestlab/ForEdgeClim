library(raster)
library(sp)

lines <- readLines("Data/2023-03-09_DTM_100cm.asc")

# Zoek alleen lijnen die beginnen met "Vertex"
vertex_lines <- grep("^Vertex", lines, value = TRUE)

# Gebruik regex om X, Y, Z eruit te halen
coords <- do.call(rbind, regmatches(vertex_lines,
                                    gregexpr("-?\\d+\\.\\d+", vertex_lines)))

# Zet om naar numerieke matrix
xyz <- as.data.frame(apply(coords, 2, as.numeric))
colnames(xyz) <- c("x", "y", "z")


r <- rasterFromXYZ(xyz)
writeRaster(r, "Data/2023-03-09_DTM_100cm_raster.asc", format = "ascii", overwrite=TRUE)
