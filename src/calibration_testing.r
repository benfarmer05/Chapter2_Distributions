library(terra)
library(dplyr)
library(mgcv)
library(scam)
library(ggplot2)

# ===== STEP 1: Bin data to 100m grid =====

# Convert your data to a SpatVector with appropriate CRS
orbicella_vect <- vect(orbicella_model_data, 
                       geom = c("lon", "lat"), 
                       crs = "EPSG:4326")

# Project to a metric CRS (adjust UTM zone for your study area)
orbicella_vect <- project(orbicella_vect, "EPSG:32617")

# Create a 100m resolution raster template covering your data extent
ext_orbicella <- ext(orbicella_vect)
r_template <- rast(ext_orbicella, resolution = 100, crs = "EPSG:32617")

# Assign cell IDs to each point
cell_ids <- cells(r_template, orbicella_vect)[, "cell"]
orbicella_df <- as.data.frame(orbicella_vect, geom = "XY")
orbicella_df$cell_id <- cell_ids

# Aggregate by grid cell: ANY presence = 1, all absences = 0
# For abundance, take mean of cover_prop where present
orbicella_binned <- orbicella_df %>%
  group_by(cell_id) %>%
  summarise(
    present = as.integer(max(present, na.rm = TRUE)),
    cover_prop = mean(cover_prop[cover_prop > 0], na.rm = TRUE),  # Mean of presences only
    depth_bathy = mean(depth_bathy, na.rm = TRUE),
    aspect = mean(aspect, na.rm = TRUE),
    slope = mean(slope, na.rm = TRUE),
    complexity = mean(complexity, na.rm = TRUE),
    dir_at_max_hsig = mean(dir_at_max_hsig, na.rm = TRUE),
    range_SST = mean(range_SST, na.rm = TRUE),
    mean_Hsig = mean(mean_Hsig, na.rm = TRUE),
    mean_SST = mean(mean_SST, na.rm = TRUE),
    mean_PAR = mean(mean_PAR, na.rm = TRUE),
    lon = mean(x, na.rm = TRUE),
    lat = mean(y, na.rm = TRUE),
    .groups = "drop"
  )

# Check the aggregation
cat("Original data rows:", nrow(orbicella_model_data), "\n")
cat("Binned data rows:", nrow(orbicella_binned), "\n")

# ========================================
# PRESENCE/ABSENCE MODEL WITH CALIBRATION
# ========================================

cat("\n\n===== PRESENCE/ABSENCE MODEL =====\n\n")

# Fit GAM with binned data
orbicella_gam_binned <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                              s(slope) + s(complexity) +
                              s(dir_at_max_hsig, bs = 'cc') +
                              s(range_SST) + s(lon, lat),
                            data = orbicella_binned,
                            family = binomial())

summary(orbicella_gam_binned)

# Extract the data actually used in the model fit
model_data_used_pa <- orbicella_gam_binned$model

# Get fitted values
fitted.response.pa <- predict(orbicella_gam_binned, type = "response")
fitted.link.pa <- predict(orbicella_gam_binned, type = "link")

# Create calibration plot
pdata.pa <- data.frame(fits = fitted.response.pa, 
                       observed = model_data_used_pa$present)

calibration_plot_pa <- ggplot(pdata.pa, aes(fits, jitter(observed, 0.1))) +
  scale_x_continuous("GAM fitted values", limits = c(0, 1)) +
  scale_y_continuous("observed occurrences", limits = c(0, 1)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, col = "black") +
  geom_smooth(col = "grey60") +
  theme_bw() +
  ggtitle("P/A Calibration - Raw Predictions")

print(calibration_plot_pa)

# Calibration statistics
cali.gam.pa <- glm(model_data_used_pa$present ~ fitted.link.pa, 
                   family = binomial)

cat("\n=== P/A Calibration Statistics ===\n")
cat("Ideal: Intercept ≈ 0, Slope ≈ 1\n\n")
print(summary(cali.gam.pa))

# Fit calibration model
cal.scam.pa <- scam(model_data_used_pa$present ~ s(fitted.response.pa, bs = "mpi"), 
                    family = binomial)

# Get calibrated predictions
corr.fitted.response.pa <- predict(cal.scam.pa, type = "response")

# Plot calibrated predictions
pdata.corr.pa <- data.frame(corr = corr.fitted.response.pa,
                            observed = model_data_used_pa$present)

calibration_plot_corrected_pa <- ggplot(pdata.corr.pa, aes(corr, jitter(observed, 0.1))) +
  scale_x_continuous("GAM: calibrated values", limits = c(0, 1)) +
  scale_y_continuous("observed occurrences", limits = c(0, 1)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, col = "black") +
  geom_smooth(col = "grey60") +
  theme_bw() +
  ggtitle("P/A Calibration - Calibrated Predictions")

print(calibration_plot_corrected_pa)





