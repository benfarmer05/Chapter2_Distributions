    
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(terra)
  library(mgcv)
  library(cmocean)
  library(ggplot2)
  library(gridExtra)
  library(tidyverse)

  source(here("src/functions.R"))
  
  ################################## setup ##################################
  
  load(here("output", "all_combined_data.rda"))
  
  #load bathy_final
  # Load the spatial metadata
  spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
  #
  # Load the raster
  bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
  #
  # Apply the stored CRS
  terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
  
  #load derived bathy rasters
  load_spat_objects(directory = 'output/output_calculate_bathy_rasters/')
  source(here("src/functions.R"))
  load(here('output', 'output_create_habitat_grid/create_habitat_grid_workspace.RData'))
  
  #save information for exporting new objects later
  existing_objects <- ls(envir = .GlobalEnv)
  
  ################################## create GAM prediction grid ##################################
  
  # First, create your stack properly
  env_stack <- c(bathy_final, 
                 aspect_terra, 
                 slope_terra, 
                 slopeofslope_terra, 
                 TPI_terra, 
                 TRI, 
                 roughness, 
                 VRM,
                 maxcurv_multiscale,
                 meancurv_multiscale, 
                 planformcurv_multiscale,
                 profilecurv_multiscale,
                 VRM_multiscale)
  
  # Check geometry - compareGeom works on individual rasters, not the stack
  # Let's check that all layers match the first one (bathy_final)
  terra::compareGeom(bathy_final, aspect_terra, stopOnError = FALSE)
  terra::compareGeom(bathy_final, slope_terra, stopOnError = FALSE)
  # etc... or loop through them:
  
  for(i in 2:nlyr(env_stack)) {
    cat("Checking layer", i, ":", names(env_stack)[i], "\n")
    print(terra::compareGeom(env_stack[[1]], env_stack[[i]], stopOnError = FALSE))
  }
  
  # Alternative: just check basic properties
  terra::res(env_stack)  # Should all be the same
  terra::ext(env_stack)  # Should all be the same
  
  # Create prediction grid
  prediction_grid <- as.data.frame(env_stack, xy = TRUE)
  
  # # Remove rows with any NA values
  # prediction_grid <- prediction_grid[complete.cases(prediction_grid), ]
  
  # Check your grid
  dim(prediction_grid)
  head(prediction_grid)

  # Better names for your variables
  names(prediction_grid) <- c("x", "y", "bathymetry", "aspect", "slope", "slope_of_slope", 
                              "TPI", "TRI", "roughness", "VRM", "max_curv", "mean_curv", 
                              "planform_curv", "profile_curv", "VRM_multiscale")  
  
  ################################## test GAMs ##################################
  


  # Convert tibble to data.frame first
  species_df <- as.data.frame(combined_benthic_data_averaged_psu)
  
  # Now create SpatVector
  species_coords <- vect(species_df, 
                         geom = c("lon", "lat"),
                         crs = "EPSG:4326")  # WGS84 geographic
  
  # Transform to match your raster CRS
  species_coords_utm <- project(species_coords, crs(bathy_final))
  
  # Extract transformed coordinates
  utm_coords <- as.data.frame(geom(species_coords_utm)[, c("x", "y")])
  combined_benthic_data_averaged_psu$x_utm <- utm_coords$x
  combined_benthic_data_averaged_psu$y_utm <- utm_coords$y
  
  # Create simple environmental stack
  env_simple <- c(bathy_final, slope_terra)
  names(env_simple) <- c("bathymetry", "slope")
  
  # Extract environmental values at species locations
  species_env_values <- terra::extract(env_simple, 
                                       cbind(combined_benthic_data_averaged_psu$x_utm, 
                                             combined_benthic_data_averaged_psu$y_utm))
  
  # Add to dataframe
  combined_benthic_data_averaged_psu$bathymetry <- species_env_values$bathymetry
  combined_benthic_data_averaged_psu$slope <- species_env_values$slope
  
  # Add Y/N columns for bathymetry and slope presence
  combined_benthic_data_averaged_psu$bathymetry_present <- ifelse(is.na(combined_benthic_data_averaged_psu$bathymetry), "N", "Y")
  combined_benthic_data_averaged_psu$slope_present <- ifelse(is.na(combined_benthic_data_averaged_psu$slope), "N", "Y")
  
  # Check how many NAs you have
  # NOTE - this is actually kind of a lot. 
  sum(is.na(combined_benthic_data_averaged_psu$bathymetry))
  sum(is.na(combined_benthic_data_averaged_psu$slope))
  
  # PLOT NAs
  # NOTE - given the high amount of NAs that are occurring due to the landmask, could go upstream and test
  #         being less aggressive with applying a landmask over the Blondeau bathy. especially since I now
  #         have much better handling of inland water ... still wouldn't fix a totally messed up coastline, though
  # NOTE - you can see when looking at the north drop, just how much resolution mismatch may confound associations
  #         with bathymetric derivatives. may need to consider dropping most or all of the north drop observations
  # STOPPING POINT - 11 JULY 2025: make triple sure that I am applying the "right" landmask to PR, not just USVI side
  #       - because it appears that might be some funny areas in the PR landmask from Blondeau still
  #
  # Create a logical vector for NA locations
  na_mask <- is.na(combined_benthic_data_averaged_psu$bathymetry) | 
    is.na(combined_benthic_data_averaged_psu$slope)
  #
  # Get coordinates for plotting
  coords_na <- combined_benthic_data_averaged_psu[na_mask, c("x_utm", "y_utm")]
  coords_valid <- combined_benthic_data_averaged_psu[!na_mask, c("x_utm", "y_utm")]
  #
  # Create the plot
  
  # Define plot extent options
  # plot_extents = ext(260000, 310000, 2030000, 2062000) #for investigating north drop
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # plot_extents = ext(341000, 379000, 2057000, 2078000) # for investigating Anegada
  # plot_extents = ext(294000, 350000, 1950000, 1975000) #for investigating St Croix
  # plot_extents = ext(280000, 320000, 2000000, 2040000) #for investigating St Thomas
  # plot_extents = ext(240000, 280000, 2000000, 2040000) #for investigating Mona Island
  # plot_extents = ext(210000, 260000, 1995000, 2050000) #for investigating PR East
  # plot_extents = ext(55000, 230000, 2030000, 2058000) #for investigating PR North
  # plot_extents = ext(20000, 80000, 1970000, 2060000) #for investigating PR West
  plot_extents = ext(52000, 230000, 1970000, 2020000) #for investigating PR South
  
  # raster
  bathy_final_clamp <- clamp(bathy_final, lower = -50, upper = 0)
  plot(slope_terra,
  # plot(bathy_final_clamp,
       col = cmocean("deep")(100),
       ext = plot_extents,
       main = "Species locations over raster\n(Red = NA values, Black = Valid values)")
  #
  # Add points - valid locations in black
  points(coords_valid$x_utm, coords_valid$y_utm, 
         col = "black", pch = 16, cex = 0.5)
  #
  # Add points - NA locations in red
  points(coords_na$x_utm, coords_na$y_utm, 
         col = "red", pch = 16, cex = 0.5)
  #
  # Add legend
  legend("topright", 
         legend = c("Valid data", "NA values"), 
         col = c("black", "red"), 
         pch = 16, 
         cex = 0.8)
  #
  # Print summary
  cat("Total points:", nrow(combined_benthic_data_averaged_psu), "\n")
  cat("Valid points:", sum(!na_mask), "\n")
  cat("NA points:", sum(na_mask), "\n")
  cat("Percentage with NAs:", round(sum(na_mask)/nrow(combined_benthic_data_averaged_psu)*100, 1), "%\n")
  
  
  
  #TEST - taking out PR & NODICE entirely
  model_data_filtered <- combined_benthic_data_averaged_psu %>%
    filter(!dataset %in% c("PRCRMP", "NODICE")) %>%  # Remove PRCRMP and NODICE datasets
    filter(!grepl("_PR", PSU))  # Remove any PSU containing '_PR'
  
  cat("Remaining datasets:", paste(unique(model_data_filtered$dataset), collapse = ", "), "\n")
  
  psu_with_pr <- model_data_filtered$PSU[grepl("_PR", model_data_filtered$PSU)]
  if(length(psu_with_pr) > 0) {
    cat("PSUs with '_PR' still remaining:", paste(unique(psu_with_pr), collapse = ", "), "\n")
  } else {
    cat("No PSUs with '_PR' remaining\n")
  }
  
  model_data <- model_data_filtered[complete.cases(model_data_filtered[, c("bathymetry", "slope", "cover")]), ]
  
  
  # Check the new cover distribution
  summary(model_data$cover)
  prop_zeros_new <- sum(model_data$cover == 0) / nrow(model_data)
  cat("New proportion of zeros:", round(prop_zeros_new, 3), "\n")
  
  # Compare with original
  prop_zeros_orig <- sum(combined_benthic_data_averaged_psu$cover == 0, na.rm = TRUE) / nrow(combined_benthic_data_averaged_psu)
  cat("Original proportion of zeros:", round(prop_zeros_orig, 3), "\n")
  
  #TEST - taking out PR & NODICE entirely
  
  
  
  # # Remove rows with NA environmental values
  # model_data <- combined_benthic_data_averaged_psu[complete.cases(combined_benthic_data_averaged_psu[, c("bathymetry", "slope", "cover")]), ]
  
  
  
  
  
  
  # Option 1: Tweedie distribution (good for zero-inflated continuous data)
  gam_tweedie <- gam(cover ~ s(bathymetry) + s(slope), 
                     data = model_data,
                     family = tw())
  
  # Option 2: If cover is bounded (0-100%), use beta regression with zeros
  # First, transform cover to (0,1) range if it's percentage
  if(max(model_data$cover, na.rm = TRUE) > 1) {
    model_data$cover_prop <- model_data$cover / 100
  } else {
    model_data$cover_prop <- model_data$cover
  }
  
  # Beta regression (handles zeros with adjustment)
  gam_beta <- gam(cover_prop ~ s(bathymetry) + s(slope), 
                  data = model_data[model_data$cover_prop > 0, ],  # exclude zeros
                  family = betar())
  
  # Option 3: Two-part model (hurdle model)
  # Part 1: Presence/absence
  model_data$present <- ifelse(model_data$cover > 0, 1, 0)
  gam_presence <- gam(present ~ s(bathymetry) + s(slope), 
                      data = model_data,
                      family = binomial())
  
  # Part 2: Abundance given presence
  gam_abundance <- gam(cover ~ s(bathymetry) + s(slope), 
                       data = model_data[model_data$cover > 0, ],
                       family = Gamma(link = "log"))

  
  # Check model summaries
  summary(gam_tweedie)
  summary(gam_beta)
  summary(gam_presence)
  summary(gam_abundance)
  
  # Plot results
  par(mfrow = c(2, 2))
  plot(gam_tweedie, pages = 1, main = "Tweedie Model")
  plot(gam_presence, pages = 1, main = "Presence Model")  
  
  
  
  ################################## plot relationships ##################################
  
  # Plot both relationships on one page
  par(mfrow = c(1, 2))
  plot(gam_tweedie, select = 1, main = "Bathymetry effect on coral cover", 
       xlab = "Bathymetry (m)", ylab = "Smooth term")
  plot(gam_tweedie, select = 2, main = "Slope effect on coral cover", 
       xlab = "Slope", ylab = "Smooth term")
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  # Create prediction data for plotting
  bathy_range <- seq(min(model_data$bathymetry, na.rm = TRUE), 
                     max(model_data$bathymetry, na.rm = TRUE), 
                     length.out = 100)
  slope_range <- seq(min(model_data$slope, na.rm = TRUE), 
                     max(model_data$slope, na.rm = TRUE), 
                     length.out = 100)
  
  # Predictions for bathymetry (holding slope at median)
  pred_data_bathy <- data.frame(
    bathymetry = bathy_range,
    slope = median(model_data$slope, na.rm = TRUE)
  )
  
  # Predictions for slope (holding bathymetry at median)
  pred_data_slope <- data.frame(
    bathymetry = median(model_data$bathymetry, na.rm = TRUE),
    slope = slope_range
  )
  
  # Get predictions with standard errors
  pred_bathy <- predict(gam_tweedie, pred_data_bathy, se.fit = TRUE, type = "response")
  pred_slope <- predict(gam_tweedie, pred_data_slope, se.fit = TRUE, type = "response")
  
  # Create plots
  par(mfrow = c(1, 2))
  
  # Bathymetry plot
  plot(pred_data_bathy$bathymetry, pred_bathy$fit, type = "l", 
       xlab = "Bathymetry (m)", ylab = "Predicted coral cover (%)",
       main = "Coral cover vs Bathymetry", lwd = 2, col = "blue")
  lines(pred_data_bathy$bathymetry, pred_bathy$fit + 1.96*pred_bathy$se.fit, lty = 2, col = "blue")
  lines(pred_data_bathy$bathymetry, pred_bathy$fit - 1.96*pred_bathy$se.fit, lty = 2, col = "blue")
  
  # Add raw data points
  points(model_data$bathymetry, model_data$cover, pch = 16, cex = 0.3, col = "gray60")
  
  # Slope plot
  plot(pred_data_slope$slope, pred_slope$fit, type = "l", 
       xlab = "Slope", ylab = "Predicted coral cover (%)",
       main = "Coral cover vs Slope", lwd = 2, col = "red")
  lines(pred_data_slope$slope, pred_slope$fit + 1.96*pred_slope$se.fit, lty = 2, col = "red")
  lines(pred_data_slope$slope, pred_slope$fit - 1.96*pred_slope$se.fit, lty = 2, col = "red")
  
  # Add raw data points
  points(model_data$slope, model_data$cover, pch = 16, cex = 0.3, col = "gray60")
  
  par(mfrow = c(1, 1))

  
  
  
  
  
  # Create prediction dataframes as above, then:
  pred_bathy_df <- data.frame(
    bathymetry = pred_data_bathy$bathymetry,
    fitted = pred_bathy$fit,
    lower = pred_bathy$fit - 1.96*pred_bathy$se.fit,
    upper = pred_bathy$fit + 1.96*pred_bathy$se.fit
  )
  
  pred_slope_df <- data.frame(
    slope = pred_data_slope$slope,
    fitted = pred_slope$fit,
    lower = pred_slope$fit - 1.96*pred_slope$se.fit,
    upper = pred_slope$fit + 1.96*pred_slope$se.fit
  )
  
  # Bathymetry plot
  p1 <- ggplot(pred_bathy_df, aes(x = bathymetry)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "blue") +
    geom_line(aes(y = fitted), color = "blue", size = 1) +
    geom_point(data = model_data, aes(x = bathymetry, y = cover), 
               alpha = 0.3, size = 0.5) +
    labs(x = "Bathymetry (m)", y = "Predicted coral cover (%)", 
         title = "Coral cover vs Bathymetry") +
    theme_minimal()
  
  # Slope plot
  p2 <- ggplot(pred_slope_df, aes(x = slope)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "red") +
    geom_line(aes(y = fitted), color = "red", size = 1) +
    geom_point(data = model_data, aes(x = slope, y = cover), 
               alpha = 0.3, size = 0.5) +
    labs(x = "Slope", y = "Predicted coral cover (%)", 
         title = "Coral cover vs Slope") +
    theme_minimal()
  
  # Display both plots
  grid.arrange(p1, p2, ncol = 2)  
  
  
  
  
  
  ################################## troubleshoot ##################################
  
  # Basic data exploration
  table(model_data$cover == 0)  # How many exact zeros?
  summary(model_data$cover)
  length(unique(model_data$cover))  # How many unique values?
  
  # Check if you have spatial/temporal clustering
  plot(model_data$lon, model_data$lat, col = ifelse(model_data$cover > 0, "red", "black"), 
       pch = 16, cex = 0.5, main = "Spatial distribution of cover data")
  
  
  
  
  
  # Check if bathymetry values make ecological sense
  summary(model_data$bathymetry)
  hist(model_data$bathymetry, breaks = 30, main = "Bathymetry distribution")
  
  # Are the depths reasonable for coral? (usually 0-40m)
  range(model_data$bathymetry, na.rm = TRUE)
  
  # Check slope values
  summary(model_data$slope)
  hist(model_data$slope, breaks = 30, main = "Slope distribution")

  
  
  
  
  # Simple correlations - should show the expected patterns
  cor(model_data$cover, model_data$bathymetry, use = "complete.obs")
  cor(model_data$cover, model_data$slope, use = "complete.obs")
  
  # Scatterplots with trend lines
  plot(model_data$bathymetry, model_data$cover, main = "Raw data: Cover vs Depth")
  abline(lm(cover ~ bathymetry, data = model_data), col = "red")
  
  plot(model_data$slope, model_data$cover, main = "Raw data: Cover vs Slope")
  abline(lm(cover ~ slope, data = model_data), col = "red")  
  
  
  
  ################################## try more variables ##################################
  
  # Extract all the complexity variables at species locations
  env_complex <- c(bathy_final, slope_terra, TPI_terra, roughness, VRM,
                   maxcurv_multiscale, meancurv_multiscale, 
                   planformcurv_multiscale, profilecurv_multiscale)
  
  names(env_complex) <- c("bathymetry", "slope", "TPI", "roughness", "VRM",
                          "max_curv", "mean_curv", "planform_curv", "profile_curv")
  
  # Extract all environmental values at species locations
  species_env_complex <- terra::extract(env_complex, 
                                        cbind(model_data_filtered$x_utm, 
                                              model_data_filtered$y_utm))
  
  # Add all variables to your filtered dataset
  model_data_filtered$bathymetry <- species_env_complex$bathymetry
  model_data_filtered$slope <- species_env_complex$slope
  model_data_filtered$TPI <- species_env_complex$TPI
  model_data_filtered$roughness <- species_env_complex$roughness
  model_data_filtered$VRM <- species_env_complex$VRM
  model_data_filtered$max_curv <- species_env_complex$max_curv
  model_data_filtered$mean_curv <- species_env_complex$mean_curv
  model_data_filtered$planform_curv <- species_env_complex$planform_curv
  model_data_filtered$profile_curv <- species_env_complex$profile_curv
  
  # Create complete cases dataset with all variables
  complexity_vars <- c("bathymetry", "slope", "TPI", "roughness", "VRM", 
                       "max_curv", "mean_curv", "planform_curv", "profile_curv", "cover")
  
  model_data_complex <- model_data_filtered[complete.cases(model_data_filtered[, complexity_vars]), ]
  
  # Check how much data you have left
  cat("Observations with all complexity variables:", nrow(model_data_complex), "\n")
  
  # Fit expanded GAM models
  gam_complex <- gam(cover ~ s(bathymetry) + s(slope) + s(TPI) + s(roughness) + s(VRM) +
                       s(max_curv) + s(mean_curv) + s(planform_curv) + s(profile_curv), 
                     data = model_data_complex, 
                     family = tw())
  
  # Check the results
  summary(gam_complex)
  
  # Plot the relationships
  plot(gam_complex, pages = 2)  # Will create multiple pages
  
  # Compare model performance
  cat("Simple model AIC:", AIC(gam_tweedie), "\n")
  cat("Complex model AIC:", AIC(gam_complex), "\n")
  
  # Check correlations between variables first
  complexity_matrix <- model_data_complex[, c("bathymetry", "slope", "TPI", "roughness", "VRM", 
                                              "max_curv", "mean_curv", "planform_curv", "profile_curv")]
  cor_matrix <- cor(complexity_matrix, use = "complete.obs")
  print(round(cor_matrix, 2))
  
  # # If you want to see which variables are most important:
  # # Stepwise selection
  # gam_step <- step.Gam(gam_complex, scope = list("s(bathymetry)" = 1, "s(slope)" = 1, 
  #                                                "s(TPI)" = 1, "s(roughness)" = 1, "s(VRM)" = 1,
  #                                                "s(max_curv)" = 1, "s(mean_curv)" = 1, 
  #                                                "s(planform_curv)" = 1, "s(profile_curv)" = 1))

  
  
  
  
  # Two-part model with complexity
  model_data_complex$present <- ifelse(model_data_complex$cover > 0, 1, 0)
  
  gam_presence_complex <- gam(present ~ s(bathymetry) + s(slope) + s(TPI) + s(roughness) + s(VRM), 
                              data = model_data_complex,
                              family = binomial())
  
  gam_abundance_complex <- gam(cover ~ s(bathymetry) + s(slope) + s(TPI) + s(roughness) + s(VRM), 
                               data = model_data_complex[model_data_complex$cover > 0, ],
                               family = Gamma(link = "log"))
  
  summary(gam_presence_complex)
  summary(gam_abundance_complex)  
  
  
  
  
  
  
  # Create a reduced model with less correlated variables
  # Keep: bathymetry (depth), TPI (topographic position), VRM (rugosity)
  # These capture the main dimensions: depth, local position, and complexity
  
  gam_reduced <- gam(cover ~ s(bathymetry) + s(TPI) + s(VRM), 
                     data = model_data_complex, 
                     family = tw())
  
  summary(gam_reduced)
  
  # Compare AIC
  cat("Full complex model AIC:", AIC(gam_complex), "\n")
  cat("Reduced model AIC:", AIC(gam_reduced), "\n")
  
  # Plot the reduced model
  plot(gam_reduced, pages = 1)
  
  
  
  
  # Add spatial smooth to account for unmeasured spatial factors
  gam_spatial <- gam(cover ~ s(bathymetry) + s(TPI) + s(VRM) + s(lon, lat), 
                     data = model_data_complex, 
                     family = tw())
  
  summary(gam_spatial)
  
  
  # Maybe relationships vary by depth zone?
  gam_interaction <- gam(cover ~ s(bathymetry) + s(TPI) + s(VRM) + 
                           s(TPI, by = bathymetry) + s(VRM, by = bathymetry), 
                         data = model_data_complex, 
                         family = tw())
  summary(gam_interaction)
  
  
  
  ################################## make predictions ##################################
  
  # Convert prediction grid UTM coordinates back to lat/lon for the spatial term
  prediction_coords_utm <- vect(cbind(prediction_grid$x, prediction_grid$y), 
                                crs = crs(bathy_final))
  prediction_coords_latlon <- project(prediction_coords_utm, "EPSG:4326")
  coords_df <- as.data.frame(geom(prediction_coords_latlon)[, c("x", "y")])
  
  # Add lat/lon to prediction grid
  prediction_grid$lon <- coords_df$x
  prediction_grid$lat <- coords_df$y
  
  # Remove any rows with NAs in the required variables
  pred_vars <- c("bathymetry", "TPI", "VRM", "lon", "lat")
  prediction_grid_clean <- prediction_grid[complete.cases(prediction_grid[, pred_vars]), ]
  
  cat("Prediction grid size:", nrow(prediction_grid_clean), "cells\n")
  
  
  
  
  # # Make predictions with the spatial model
  # predictions <- predict(gam_spatial, prediction_grid_clean, type = "response")
  # 
  # # Add predictions to the grid
  # prediction_grid_clean$predicted_cover <- predictions
  # 
  # # Check prediction range
  # summary(prediction_grid_clean$predicted_cover)
  
  
  # # Simple approach with time estimates
  # start_time <- Sys.time()
  # 
  # # Test prediction on small subset to estimate time
  # test_subset <- prediction_grid_clean[1:1000, ]
  # test_start <- Sys.time()
  # test_pred <- predict(gam_spatial, test_subset, type = "response")
  # test_time <- as.numeric(difftime(Sys.time(), test_start, units = "secs"))
  # 
  # # Estimate total time
  # estimated_total <- (test_time / 1000) * nrow(prediction_grid_clean)
  # cat("Estimated total time:", round(estimated_total/60, 1), "minutes\n")
  # 
  # # Now do full prediction
  # predictions <- predict(gam_spatial, prediction_grid_clean, type = "response")
  # cat("Actual time:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes\n")
  
  
  # Resample your environmental stack to 200m resolution
  env_simple_200m <- aggregate(env_simple, fact = 4, fun = "mean")  # 4x aggregation = 200m
  
  # Check the new resolution
  res(env_simple_200m)  # Should show 200, 200
  
  # Create prediction grid from 200m rasters
  prediction_grid_200m <- as.data.frame(env_simple_200m, xy = TRUE)
  names(prediction_grid_200m) <- c("x", "y", "bathymetry", "slope")
  
  # Add TPI and VRM at 200m
  env_complex_200m <- aggregate(env_complex, fact = 4, fun = "mean")
  complex_grid_200m <- as.data.frame(env_complex_200m, xy = TRUE)
  names(complex_grid_200m) <- c("x", "y", "bathymetry", "slope", "TPI", "roughness", "VRM",
                                "max_curv", "mean_curv", "planform_curv", "profile_curv")
  
  # Keep just the variables you need
  prediction_grid_200m <- complex_grid_200m[, c("x", "y", "bathymetry", "TPI", "VRM")]
  
  # Remove NAs
  prediction_grid_200m <- prediction_grid_200m[complete.cases(prediction_grid_200m), ]
  
  cat("200m grid size:", nrow(prediction_grid_200m), "cells\n")
  cat("50m grid size was:", nrow(prediction_grid_clean), "cells\n")
  cat("Reduction factor:", round(nrow(prediction_grid_clean)/nrow(prediction_grid_200m), 1), "\n")
  
  # Convert UTM to lat/lon for the spatial term
  prediction_coords_utm_200m <- vect(cbind(prediction_grid_200m$x, prediction_grid_200m$y), 
                                     crs = crs(bathy_final))
  prediction_coords_latlon_200m <- project(prediction_coords_utm_200m, "EPSG:4326")
  coords_df_200m <- as.data.frame(geom(prediction_coords_latlon_200m)[, c("x", "y")])
  
  prediction_grid_200m$lon <- coords_df_200m$x
  prediction_grid_200m$lat <- coords_df_200m$y
  
  cat("Final 200m grid size:", nrow(prediction_grid_200m), "cells\n")
  
  # Step 4: Estimate time before full prediction
  start_time <- Sys.time()
  
  # Test prediction on small subset to estimate time
  test_subset_200m <- prediction_grid_200m[1:1000, ]
  test_start <- Sys.time()
  test_pred_200m <- predict(gam_spatial, test_subset_200m, type = "response")
  test_time <- as.numeric(difftime(Sys.time(), test_start, units = "secs"))
  
  # Estimate total time
  estimated_total_200m <- (test_time / 1000) * nrow(prediction_grid_200m)
  cat("200m grid size:", nrow(prediction_grid_200m), "cells\n")
  cat("Test time for 1000 cells:", round(test_time, 2), "seconds\n")
  cat("Estimated total time:", round(estimated_total_200m/60, 1), "minutes\n")
  
  # Compare to original estimate
  original_cells <- nrow(prediction_grid_clean)
  reduction_factor <- original_cells / nrow(prediction_grid_200m)
  cat("Reduction factor:", round(reduction_factor, 1), "x smaller\n")
  cat("Original estimate was 61 minutes\n")
  cat("New estimate should be ~", round(61/reduction_factor, 1), "minutes\n")
  
  # Ask user if they want to proceed
  cat("\nProceed with full prediction? (y/n)\n")
  user_input <- readline()
  
  if(tolower(user_input) == "y") {
    cat("Starting full prediction...\n")
    
    # Full prediction with timing
    full_start <- Sys.time()
    predictions_200m <- predict(gam_spatial, prediction_grid_200m, type = "response")
    prediction_grid_200m$predicted_cover <- predictions_200m
    full_end <- Sys.time()
    
    cat("Actual prediction time:", round(difftime(full_end, full_start, units = "mins"), 1), "minutes\n")
    cat("Prediction complete!\n")
    
  } else {
    cat("Prediction cancelled. Consider further reducing grid size or using sampling.\n")
  }
  
  # Convert to raster
  pred_raster_200m <- rast(prediction_grid_200m[, c("x", "y", "predicted_cover")], 
                           crs = crs(bathy_final))
  
  # Plot
  plot(pred_raster_200m, 
       main = "Predicted Coral Cover (%) - 200m resolution",
       col = viridis(100))
  
  # Add survey points
  points(model_data_complex$x_utm, model_data_complex$y_utm, 
         pch = 21, 
         bg = heat.colors(10)[cut(model_data_complex$cover, breaks = 10)],
         cex = 0.8)
  
  
  ################################## maps ##################################
  
  
  # library(viridis)
  library(sf)
  
  # Convert predictions back to raster
  pred_raster <- rast(prediction_grid_clean[, c("x", "y", "predicted_cover")], 
                      crs = crs(bathy_final))
  
  # Create a nice map
  # Option 1: Using terra plot
  plot(pred_raster, 
       main = "Predicted Coral Cover (%)",
       col = viridis(100),
       range = c(0, max(prediction_grid_clean$predicted_cover, na.rm = TRUE)))
  
  # Add your actual survey points
  points(model_data_complex$x_utm, model_data_complex$y_utm, 
         pch = 21, 
         bg = heat.colors(10)[cut(model_data_complex$cover, breaks = 10)],
         cex = 0.8)
  
  # Add legend for points
  legend("topright", 
         title = "Observed Cover",
         legend = c("0-10%", "10-20%", "20-30%", "30%+"),
         pch = 21,
         pt.bg = heat.colors(4),
         cex = 0.8)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  