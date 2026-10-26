install.packages("raster")
install.packages("sf")      # For vector data (if needed)
install.packages("psych")   # For PCA and loadings
library(raster)
library(sf)
library(psych)
BEI <- raster("F:/ULEI_2014/BEI_14_PC1.tif")
API <- raster("F:/ULEI_2014/Dhaka_L8_API_MedianComposite_2014.tif")
UHIEII<- raster("F:/ULEI_2014/UHIEII_Dhaka_2014.tif")
NLI <- raster("F:/ULEI_2014/Dhaka_VIIRS_NLI_March2014.tif")
indices_stack <- stack(BEI, API, UHIEII, NLI)
names(indices_stack) <- c("BEI", "API", "UHIEII", "NLI")
normalize_raster <- function(r) {
  r_min <- cellStats(r, min)
  r_max <- cellStats(r, max)
  return((r - r_min) / (r_max - r_min))
}
indices_norm <- stack(lapply(1:nlayers(indices_stack),function(i) normalize_raster(indices_stack[[i]])))
names(indices_norm) <- names(indices_stack)
w_FAHP <-c(BEI = 0.30, API = 0.25, UHIEII = 0.20, NLI = 0.25)
mat <- as.data.frame(getValues(indices_norm))
mat <- mat[complete.cases(mat),]
pca <- principal(mat, nfactors = 1, rotate = "none")
pc1_loadings <- abs(pca$loadings[, 1])
w_PCA <- pc1_loadings / sum(pc1_loadings)
names(w_PCA) <- c("BEI", "API", "UHIEII", "NLI")
w_combined <- (w_FAHP + w_PCA) / 2
ULEI <- indices_norm[[1]] * w_combined["BEI"] +
  indices_norm[[2]] * w_combined["API"] +
  indices_norm[[3]] * w_combined["UHIEII"] +
  indices_norm[[4]] * w_combined["NLI"]
writeRaster(ULEI, filename = "F:/ULEI_2014", "ULEI_2014.tif", format = "GTiff", overwrite = TRUE)
