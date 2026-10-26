library(terra)   # raster operations
library(sf)      # vector handling
library(spdep)   # LISA analysis
setwd("C:/LISA_ULEI")
r <- rast("ULEI_2004.tif")
r_ag <- aggregate(
  r,
  fact = 20,
  fun = mean,
  na.rm = TRUE
)
p <- as.polygons(r_ag, dissolve = FALSE)
p <- st_as_sf(p)
val_col <- names(p)[!names(p) %in% "geometry"]
names(p)[names(p) == val_col] <- "ULEI"
p <- p[!is.na(p$ULEI), ]
nb <- poly2nb(p, queen = TRUE)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
vals <- p$ULEI
lisa <- localmoran(vals, lw, zero.policy = TRUE)
p$Ii   <- lisa[, 1]  # Local Moran’s I
p$Z_Ii <- lisa[, 4]  # Z-score
p$P_Ii <- lisa[, 5]  # p-value
mean_val <- mean(vals, na.rm = TRUE)
lag_val  <- lag.listw(lw, vals, zero.policy = TRUE)
p$Cluster <- "Not significant"
p$Cluster[p$P_Ii <= 0.05 & vals >= mean_val & lag_val >= mean_val] <- "High-High"
p$Cluster[p$P_Ii <= 0.05 & vals <  mean_val & lag_val <  mean_val] <- "Low-Low"
p$Cluster[p$P_Ii <= 0.05 & vals >= mean_val & lag_val <  mean_val] <- "High-Low"
p$Cluster[p$P_Ii <= 0.05 & vals <  mean_val & lag_val >= mean_val] <- "Low-High"
st_write(
  p,
  "ULEI_LISA_2004.shp",
  delete_layer = TRUE
)
cat("✅ LISA completed successfully for NRSEI 2004\n")
p <- st_transform(p, 32646)
p$Area_km2 <- as.numeric(st_area(p)) / 1e6
head(p[, c("ULEI", "Cluster", "Area_km2")])
install.packages("dplyr")
library(dplyr)
df <- st_drop_geometry(p)
area_by_cluster <- aggregate(
  Area_km2 ~ Cluster,
  data = df,
  FUN = sum,
  na.rm = TRUE
)
area_by_cluster$Percent <-
  100 * area_by_cluster$Area_km2 / sum(area_by_cluster$Area_km2)

print(area_by_cluster)
cat("✅ LISA completed successfully for NULEI2 004\n")
