library(raster)

# Load rasters
ndvi <- raster("F:/BEI_2004/NDVI_O4_Landsat07.tif")
ibi <- raster("F:/BEI_2004/Dhaka_L5_IBI_JanJul2004_Full (1).tif")
wetness <- raster("F:/BEI_2004/Wetness_Index_Jan_Jul_2004.tif")

# Resample to common grid
template <- ndvi
ibi_rs <- resample(ibi, template, method="bilinear")
wetness_rs <- resample(wetness, template, method="bilinear")

# Normalize function
normalize_raster <- function(r) {
  (r - cellStats(r, min, na.rm=TRUE)) / (cellStats(r, max, na.rm=TRUE) - cellStats(r, min, na.rm=TRUE))
}

ndvi_norm <- normalize_raster(ndvi)
ibi_norm <- normalize_raster(ibi_rs)
wetness_norm <- normalize_raster(wetness_rs)

# Stack normalized rasters
stack_norm <- stack(ndvi_norm, ibi_norm, wetness_norm)
names(stack_norm) <- c("NDVI", "IBI", "Wetness")

# Convert to matrix (pixels x layers)
vals <- as.data.frame(getValues(stack_norm))

# Remove NA rows for PCA
vals_complete <- na.omit(vals)

# PCA
pca_result <- prcomp(vals_complete, center = TRUE, scale. = FALSE)
pc1_scores <- pca_result$x[, 1]   # BEI₀

# Normalize BEI₀ (PC1) according to your equation
pc1_min <- min(pc1_scores, na.rm = TRUE)
pc1_max <- max(pc1_scores, na.rm = TRUE)
bei_final_norm <- (pc1_scores - pc1_min) / (pc1_max - pc1_min)

# Assign normalized BEI to raster, preserving NA structure
bei_final <- raster(ndvi_norm)
vals_full_norm <- rep(NA, nrow(vals))
vals_full_norm[complete.cases(vals)] <- bei_final_norm
bei_final <- setValues(bei_final, vals_full_norm)

# Export final normalized BEI raster
writeRaster(bei_final, "F:/BEI_2004/BEI_final_normalized_2004.tif", overwrite=TRUE)

cat("Eigenvalues & Variance Explained\n")
print(summary(pca_result))
cat("\nPCA Loadings (Variable Contributions)\n")
print(pca_result$rotation)
cat("\nFirst 10 PC1 Scores\n")
print(head(pca_result$x[, 1], 10))
cat("\nFirst 10 rows of PCA score matrix\n")
print(head(pca_result$x, 10))
pca_var <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
cat("\nPercentage Variance Explained by Each PC\n")
print(round(pca_var, 2))
print(pca_result$rotation[,1])
