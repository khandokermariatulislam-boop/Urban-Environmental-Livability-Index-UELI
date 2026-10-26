# Global Moran's I for ULEI (robust, no dplyr required)
# - Auto-finds a ULEI raster in working directory if exact name not present
# - Builds polygons, weights, computes analytical Moran's I and permutation test
# - Exports CSV summary, Moran scatterplot PNG, and shapefile with results
#
# Usage: edit setwd() to your folder or run from that folder

library(terra)   # raster operations
library(sf)      # vector handling
library(spdep)   # Moran's I (global) and spatial weights

# ---------- Configuration ----------
wd <- "D:/Moran"    # change as needed (or run setwd() interactively)
expected_names <- c("ULEI_2024.tif")  # try these first
aggregate_factor <- 60
nsim <- 999              # permutation reps (increase for more precision)
seed <- 12345
output_prefix <- "ULEI_2024"  # used for filenames (adjust if using 2024)

# ---------- Prepare workspace and find raster ----------
setwd(wd)
message("Working directory: ", getwd())

tif_path <- NULL
for (nm in expected_names) {
  if (file.exists(nm)) {
    tif_path <- nm
    break
  }
}

if (is.null(tif_path)) {
  candidates <- list.files(pattern = "ULEI.*\\.tif$",
                           full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(candidates) == 0) {
    stop("No ULEI TIFF found in working directory. Check filename or set correct working directory.")
  } else {
    tif_path <- candidates[1]
    message("Using raster file: ", tif_path)
  }
} else {
  message("Found expected raster: ", tif_path)
}

# ---------- Read and aggregate raster ----------
r <- tryCatch(rast(tif_path), error = function(e) stop("Failed to open raster: ", e$message))
r_ag <- aggregate(r, fact = aggregate_factor, fun = mean, na.rm = TRUE)

# ---------- Convert aggregated raster to polygons ----------
p <- as.polygons(r_ag, dissolve = FALSE)
p <- st_as_sf(p)

# Ensure ULEI column name
val_col <- names(p)[!names(p) %in% "geometry"]
if (length(val_col) != 1) stop("Unexpected number of non-geometry columns in polygon object.")
names(p)[names(p) == val_col] <- "ULEI"

# Drop NA rows
p <- p[!is.na(p$ULEI), ]

# ---------- Spatial weights ----------
nb <- poly2nb(p, queen = TRUE)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

vals <- p$ULEI
n <- length(vals)
if (n < 3) stop("Too few observations to compute Moran's I (n < 3).")

# ---------- Analytical Moran's I (asymptotic test) ----------
global_analytical <- moran.test(vals, lw, zero.policy = TRUE)
print(global_analytical)

I_obs_analytical <- as.numeric(global_analytical$estimate["Moran I statistic"])
expected_I_analytical <- as.numeric(global_analytical$estimate["expected"])
z_analytical <- as.numeric(global_analytical$statistic)   # z-value from moran.test
p_analytical <- as.numeric(global_analytical$p.value)

# Derive variance from Z if possible: Var = ((I - E)/Z)^2 when Z != 0
var_I_analytical <- if (!is.na(z_analytical) && z_analytical != 0) {
  ((I_obs_analytical - expected_I_analytical) / z_analytical)^2
} else {
  NA_real_
}

# ---------- Permutation (Monte Carlo) Moran's I ----------
set.seed(seed)
global_mc <- moran.mc(vals, lw, nsim = nsim, zero.policy = TRUE, alternative = "two.sided")
print(global_mc)

I_obs_perm <- as.numeric(global_mc$statistic)
p_perm <- as.numeric(global_mc$p.value)

# Extract simulated values (spdep versions differ in slot name)
sim_vals <- NULL
if (!is.null(global_mc$res)) sim_vals <- global_mc$res
if (is.null(sim_vals) && !is.null(global_mc$simulate)) sim_vals <- global_mc$simulate
if (is.null(sim_vals) && !is.null(global_mc$sim)) sim_vals <- global_mc$sim

perm_mean <- if (!is.null(sim_vals)) mean(sim_vals) else NA_real_
perm_var  <- if (!is.null(sim_vals)) var(sim_vals) else NA_real_
perm_sd   <- if (!is.null(sim_vals)) sd(sim_vals) else NA_real_
z_perm    <- if (!is.null(sim_vals) && perm_sd > 0) (I_obs_perm - perm_mean) / perm_sd else NA_real_

# ---------- Moran scatterplot ----------
# ============================================================
# COLORFUL MORAN SCATTERPLOT (FITTED TO YOUR CODE)
# ============================================================
# ---------- Moran's I Calculation ----------
global_analytical <- moran.test(vals, lw, zero.policy = TRUE)
I_obs <- as.numeric(global_analytical$estimate["Moran I statistic"])

# ---------- Refined Moran Scatterplot (Blue Shade & Bold Dots) ----------

# 1. Data Preparation
z_vals <- as.numeric(scale(vals))
lag_vals <- lag.listw(lw, z_vals, zero.policy = TRUE)

# 1. Calculate symmetrical limits
# This finds the furthest point from zero in any direction
max_val <- max(abs(c(z_vals, lag_vals)), na.rm = TRUE)
plot_limit <- c(-max_val, max_val)

# 2. File Setup (600 DPI)
png("Centered_Moran_Scatterplot.png", width = 4800, height = 4800, res = 600) 

# Large margins to prevent huge fonts from clipping
par(mar = c(8, 9, 6, 2))

# Classify Moran quadrants
quad <- ifelse(z_vals >= 0 & lag_vals >= 0, "High-High",
               ifelse(z_vals <= 0 & lag_vals <= 0, "Low-Low",
                      ifelse(z_vals >= 0 & lag_vals <= 0, "High-Low",
                             "Low-High")))

# Define colors (journal standard)
quad_colors <- c(
  "High-High" = "#2C7BB6",  # red
  "Low-Low"   = "#2C7BB6",  # blue
  "High-Low"  = "#2C7BB6",  # orange
  "Low-High"  = "#2C7BB6"   # light blue
)

point_colors <- quad_colors[quad]

# Save colorful Moran scatterplot
png("Colorful_Moran_Scatterplot_ULEI_2024.png",
    width = 900, height = 700)

plot(
  z_vals,
  lag_vals,
  col = point_colors,
  pch = 19,
  cex = 0.8,
  cex.lab = 1.2,           # MASSIVE X and Y axis titles
  cex.axis = 1.2,          # MASSIVE numbers on the axes
  cex.main = 2,          # MASSIVE main title at the top
  xlab = "Standardized ULEI of 2024",
  ylab = "Spatial Lag of ULEI of 2024",
  main = paste("Global Moran's I:", round(I_obs, 3)),
  frame.plot = TRUE

)

# Reference lines
abline(h = 0, v = 0, lty = 2, col = "black", lwd = 2)

# Regression line (slope = Moran’s I)
abline(lm(lag_vals ~ z_vals), col = "red", lwd = 4)


dev.off()

cat("✅ Colorful Moran scatterplot saved as Colorful_Moran_Scatterplot_ULEI_2024.png\n")
getwd()

# ---------- Attach global results to polygon layer and save ----------
p$Global_I <- I_obs_analytical
p$Global_Z <- z_analytical
p$Global_p <- p_analytical
shp_name <- paste0("ULEI_polygons_with_GlobalI_2024", output_prefix, ".shp")
st_write(p, shp_name, delete_layer = TRUE)
message("Wrote shapefile with Global I: ", shp_name)

# ---------- Console summary ----------
cat("\n--- Global Moran's I Summary ---\n")
cat(sprintf("Analytical: I = %.6f, E[I] = %.6f, Var(I) = %s, Z = %s, p = %s\n",
            I_obs_analytical,
            expected_I_analytical,
            if (!is.na(var_I_analytical)) format(var_I_analytical, digits = 6) else "NA",
            if (!is.na(z_analytical)) format(z_analytical, digits = 6) else "NA",
            if (!is.na(p_analytical)) format(p_analytical, digits = 6) else "NA"))

cat(sprintf("Permutation (%d sims): I_obs = %.6f, perm_mean = %s, perm_sd = %s, Z_perm = %s, p_perm = %s\n",
            nsim,
            I_obs_perm,
            if (!is.na(perm_mean)) format(perm_mean, digits = 6) else "NA",
            if (!is.na(perm_sd)) format(perm_sd, digits = 6) else "NA",
            if (!is.na(z_perm)) format(z_perm, digits = 6) else "NA",
            if (!is.na(p_perm)) format(p_perm, digits = 6) else "NA"))

cat("\nCompleted successfully.\n")

