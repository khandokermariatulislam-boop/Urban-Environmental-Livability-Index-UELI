# ==============================================================================
# name: ulei_validation_sensitivity.R
# Purpose: ULEI validation (Confusion Matrix, AUC, Moran's I) 
#          & Sensitivity Analysis (OAT, Monte Carlo, PRCC)
# ==============================================================================

# 1. Define as a character vector (use c() not list())
libs <- c("terra", "sf", "gtools", "caret", "pROC", "spdep", 
          "ggplot2", "dplyr", "reshape2", "sensitivity")

# 2. Use a simple loop (cleaner and avoids the lapply argument mismatch)
for (lib in libs) {
  if (!require(lib, character.only = TRUE)) {
    install.packages(lib)
    library(lib, character.only = TRUE)
  }
}
# -------------------------
# 2. USER INPUTS
# -------------------------
setwd("C:/LISA_ULEI/validation")

# Criterion filenames (Order must match fahp_weights)
criteria_files <- c("F:/ULEI_2004/BEI_PC1_2004.tif", "F:/ULEI_2004/Dhaka_Median_API_Map_Mar-Jul_2004.tif", 
                    "F:/ULEI_2004/Dhaka_NLI_Mean_2004.tif", "F:/ULEI_2004/UHIEII_Dhaka_2004.tif")

# FAHP weights from your results
fahp_weights <- c(0.30, 0.25, 0.25, 0.20)

# Validation parameters
validation_points_file <- NA  # Replace with "points.shp" or "points.csv"
validation_obs_field   <- "observed"
ulei_breaks            <- c(0, 0.33, 0.66, 1.0)
high_threshold         <- 0.66

# Simulation settings
mc_runs   <- 500
oat_steps <- 25

# -------------------------
# 3. DATA PREPARATION
# -------------------------
# Load and align rasters
r_stack <- rast(criteria_files)
names(r_stack) <- c("BEI", "API", "NLI", "UHIEII")

# Min-Max Normalization (0-1)
for(i in 1:nlyr(r_stack)) {
  min_val <- min(values(r_stack[[i]]), na.rm = TRUE)
  max_val <- max(values(r_stack[[i]]), na.rm = TRUE)
  r_stack[[i]] <- (r_stack[[i]] - min_val) / (max_val - min_val)
}

# Function to compute ULEI
compute_ulei <- function(weights, rasters) {
  weighted.mean(rasters, w = weights)
}

# Generate baseline map
ulei_map <- compute_ulei(fahp_weights, r_stack)
writeRaster(ulei_map, "ULEI_Baseline.tif", overwrite=TRUE)

# -------------------------
# 4. VALIDATION
# -------------------------
if (!is.na(validation_points_file)) {
  vp_sf <- if(grepl(".csv$", validation_points_file)) {
    st_as_sf(read.csv(validation_points_file), coords = c("x","y"), crs = crs(r_stack))
  } else { st_read(validation_points_file) }
  
  # Extract and Classify
  vp_df <- as.data.frame(vp_sf)
  vp_df$pred_ulei <- extract(ulei_map, vect(vp_sf))[,2]
  vp_df$pred_class <- cut(vp_df$pred_ulei, breaks = ulei_breaks, include.lowest = TRUE, labels = c("Low", "Med", "High"))
  vp_df$obs_class  <- as.factor(vp_df[[validation_obs_field]])
  
  # Metrics
  print(confusionMatrix(vp_df$pred_class, vp_df$obs_class))
  
  # Moran's I (Spatial Autocorrelation)
  coords <- st_coordinates(vp_sf)
  lw <- nb2listw(knn2nb(knearneigh(coords, k = 8)))
  print(moran.test(vp_df$pred_ulei, lw))
}

# -------------------------
# 5. OAT SENSITIVITY (One-At-A-Time)
# -------------------------
oat_results <- list()
for(i in 1:nlyr(r_stack)) {
  # Vary weight from 20% to 500% of original, capped at 0.9
  w_seq <- seq(max(0.01, fahp_weights[i]*0.2), min(0.9, fahp_weights[i]*5), length.out = oat_steps)
  
  metrics <- data.frame(weight = w_seq, mean_ulei = NA, criterion = names(r_stack)[i])
  
  for(j in seq_along(w_seq)) {
    w_tmp <- fahp_weights
    w_tmp[i] <- w_seq[j]
    # Re-normalize others to sum to 1
    others <- setdiff(1:nlyr(r_stack), i)
    w_tmp[others] <- w_tmp[others] / sum(w_tmp[others]) * (1 - w_tmp[i])
    
    res <- compute_ulei(w_tmp, r_stack)
    metrics$mean_ulei[j] <- global(res, "mean", na.rm=TRUE)[1,1]
  }
  oat_results[[i]] <- metrics
}

ggplot(bind_rows(oat_results), aes(x = weight, y = mean_ulei, color = criterion)) +
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  labs(title = "OAT Sensitivity Analysis") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Centers and bolds the title
  )

# -------------------------
# 6. MONTE CARLO & PRCC
# -------------------------
set.seed(123)
# Sample weights from Dirichlet distribution
W <- rdirichlet(mc_runs, rep(1, nlyr(r_stack)))
colnames(W) <- names(r_stack)

# Extract values to matrix for speed
vals_mat <- as.matrix(values(r_stack))
valid_idx <- which(complete.cases(vals_mat))
vals_clean <- vals_mat[valid_idx, ]

# Run MC Simulations
mc_outputs <- numeric(mc_runs)
for(i in 1:mc_runs) {
  mc_outputs[i] <- mean(vals_clean %*% W[i, ], na.rm = TRUE)
}

# PRCC Calculation
pcc_res <- pcc(X = as.data.frame(W), y = mc_outputs, rank = TRUE, nboot = 500)

# Tidy PRCC Results for Plotting
prcc_df <- data.frame(
  criterion = names(r_stack),
  estimate  = pcc_res$PRCC$original,
  low       = pcc_res$PRCC$`min. c.i.`,
  high      = pcc_res$PRCC$`max. c.i.`
)

ggplot(prcc_df, aes(x = reorder(criterion, estimate), y = estimate)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.2) +
  coord_flip() + 
  theme_minimal() +
  labs(title = "Global Sensitivity: PRCC with 95% CI", 
       y = "Partial Rank Correlation", 
       x = "Criterion") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold") # Centers and bolds the title
  )

# -------------------------
# 7. EXPORT
# -------------------------
write.csv(prcc_df, "Sensitivity_Results_2004.csv", row.names = FALSE)


# 1. Extract the PRCC results table
# This includes the original estimate, bias, std. error, and 95% CIs
prcc_values <- as.data.frame(pcc_res$PRCC)

# 2. Add the criterion names as a column for better readability
prcc_values$Criterion <- rownames(prcc_values)

# 3. Reorder columns to put Criterion first
prcc_values <- prcc_values[, c("Criterion", "original", "min. c.i.", "max. c.i.")]

# 4. Rename columns for clarity
colnames(prcc_values) <- c("Criterion", "PRCC_Estimate", "Lower_CI", "Upper_CI")

# 5. Print the values to the console
print("--- PRCC Sensitivity Values ---")
print(prcc_values)

# 6. Optional: Save these values to a CSV file for your report
write.csv(prcc_values, "PRCC_Sensitivity_Values_2004.csv", row.names = FALSE)


# View the full results object
print(pcc_res)

# View just the summary table with the values, bias, and CIs
summary(pcc_res)

# To see the raw rank-transformed correlation matrix (internal)
# This shows how criteria interact with each other
print(pcc_res$PRCC)
