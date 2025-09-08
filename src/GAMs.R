  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(terra)
  library(mgcv)
  library(cmocean)
  library(ggplot2)
  library(gridExtra)
  library(tidyverse)
  library(gratia)
  
  # library(gam.hp)
  # library(pROC)
  
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
  
  #load derived bathy rasters, and oceanographic rasters
  load_spat_objects(directory = 'output/output_calculate_bathy_rasters/')
  load_spat_objects(directory = 'output/output_calculate_ocean_rasters/')
  
  #pull habitat grid
  source(here("src/functions.R"))
  load(here('output', 'output_create_habitat_grid/create_habitat_grid_workspace.RData'))
  
  #save information for exporting new objects later
  existing_objects <- ls(envir = .GlobalEnv)
  
  
  ################################## SITE DATA SETUP ##################################
  
  # Convert tibble to data.frame first
  species_df <- as.data.frame(combined_benthic_data_averaged_psu)
  
  site_data <- combined_benthic_data_averaged_psu
  
  # Now create SpatVector
  species_coords <- vect(species_df,
                         geom = c("lon", "lat"),
                         crs = "EPSG:4326")  # WGS84 geographic
  
  # Transform to match your raster CRS
  species_coords_utm <- project(species_coords, crs(bathy_final))
  
  # Extract transformed coordinates
  utm_coords <- as.data.frame(geom(species_coords_utm)[, c("x", "y")])
  site_data$x_utm <- utm_coords$x
  site_data$y_utm <- utm_coords$y
  
  # Create environmental stacks
  env_simple <- c(bathy_final, slope_terra)
  names(env_simple) <- c("depth_bathy", "slope")
  
  env_complex <- c(bathy_final, aspect_terra, slope_terra, slopeofslope_terra, TPI_terra, VRM,
                   planformcurv_multiscale, SAPA, max_hsig_raster, dir_at_max_hsig_raster,
                   mean_hsig_raster, mean_sst_raster, range_sst_raster, range_par_raster,
                   mean_chlor_a_raster, mean_kd490_raster, mean_spm_raster, dist_to_land_raster,
                   distance_to_deep_raster, bov_full)
  
  names(env_complex) <- c("depth", "aspect", "slope", "complexity", "TPI", "VRM", "planform_curv",
                          "SAPA", "max_Hsig", "dir_at_max_hsig", "mean_Hsig", "mean_SST",
                          "range_SST", "range_PAR", "mean_chla", "mean_kd490", "mean_spm", "dist_to_land",
                          "dist_to_deep", "max_BOV")
  
  # Extract environmental values for simple variables
  site_env_values <- terra::extract(env_simple,
                                       cbind(site_data$x_utm,
                                             site_data$y_utm))
  
  # Add to dataframe
  site_data$depth_bathy <- site_env_values$depth_bathy
  site_data$slope <- site_env_values$slope
  
  # Add Y/N columns for bathymetry and slope presence
  site_data$bathymetry_present <- ifelse(is.na(site_data$depth_bathy), "N", "Y")
  site_data$slope_present <- ifelse(is.na(site_data$slope), "N", "Y")
  
  # Check how many NAs you have
  sum(is.na(site_data$depth_bathy))
  sum(is.na(site_data$slope))
  
  ################################## STATIC DEFINITIONS ##################################
  
  # Define variables that are consistent across all species
  complexity_vars <- c("depth", "aspect", "slope", "TPI", "planform_curv", 
                       "max_Hsig", "dir_at_max_hsig", "mean_Hsig",  "mean_SST",
                       "range_SST", "year", "date", "lat", "lon", "cover", "range_PAR",
                       "mean_chla", "mean_kd490", "mean_spm", "dist_to_land",
                       "dist_to_deep", "max_BOV")
  
  # Function to extract environmental data (reusable)
  extract_env_data <- function(species_data) {
    species_env_complex <- terra::extract(env_complex,
                                          cbind(species_data$x_utm,
                                                species_data$y_utm))
    return(species_env_complex)
  }
  
  # Function to add environmental variables to species data
  add_env_variables <- function(species_data, species_env_complex, model_data_filtered) {
    species_data$depth_bathy <- species_env_complex$depth
    species_data$aspect <- species_env_complex$aspect
    species_data$slope <- species_env_complex$slope
    species_data$complexity <- species_env_complex$complexity
    species_data$TPI <- species_env_complex$TPI
    species_data$VRM <- species_env_complex$VRM
    species_data$planform_curv <- species_env_complex$planform_curv
    species_data$SAPA <- species_env_complex$SAPA
    species_data$max_Hsig <- species_env_complex$max_Hsig
    species_data$dir_at_max_hsig <- species_env_complex$dir_at_max_hsig
    species_data$mean_Hsig <- species_env_complex$mean_Hsig
    species_data$mean_SST <- species_env_complex$mean_SST
    species_data$range_SST <- species_env_complex$range_SST
    species_data$range_PAR <- species_env_complex$range_PAR
    species_data$mean_chla <- species_env_complex$mean_chla
    species_data$mean_kd490 <- species_env_complex$mean_kd490
    species_data$mean_spm <- species_env_complex$mean_spm
    species_data$dist_to_land <- species_env_complex$dist_to_land
    species_data$dist_to_deep <- species_env_complex$dist_to_deep
    species_data$max_BOV <- species_env_complex$max_BOV
    species_data$year <- year(model_data_filtered$date)
    
    return(species_data)
  }
  
  # Function to create correlation matrix
  create_correlation_matrix <- function(model_data_complex) {
    complexity_matrix <- model_data_complex[, c("depth_bathy", "TPI", "slope", "planform_curv",
                                                "range_SST", "mean_SST", "dir_at_max_hsig", "max_Hsig",
                                                "mean_Hsig", "year", "range_PAR", "mean_chla", "mean_kd490",
                                                "mean_spm", "dist_to_land", "dist_to_deep", "max_BOV")]
    
    cor_matrix <- cor(complexity_matrix, use = "complete.obs")
    return(cor_matrix)
  }
  
  ################################## SIMPLE SITE GAMs ##################################
  
  # PLOT NAs - diagnostic section
  # Create a logical vector for NA locations
  na_mask <- is.na(site_data$depth_bathy) |
    is.na(site_data$slope)
  
  # Get coordinates for plotting
  coords_na <- site_data[na_mask, c("x_utm", "y_utm")]
  coords_valid <- site_data[!na_mask, c("x_utm", "y_utm")]
  
  # Define plot extent options (uncomment as needed)
  # plot_extents = ext(260000, 310000, 2030000, 2062000) #for investigating north drop
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # plot_extents = ext(341000, 379000, 2057000, 2078000) # for investigating Anegada
  # plot_extents = ext(294000, 350000, 1950000, 1975000) #for investigating St Croix
  # plot_extents = ext(280000, 320000, 2000000, 2040000) #for investigating St Thomas
  plot_extents = ext(240000, 280000, 2000000, 2040000) #for investigating Mona Island
  # plot_extents = ext(210000, 260000, 1995000, 2050000) #for investigating PR East
  # plot_extents = ext(55000, 230000, 2030000, 2058000) #for investigating PR North
  # plot_extents = ext(20000, 80000, 1970000, 2060000) #for investigating PR West
  # plot_extents = ext(52000, 230000, 1970000, 2020000) #for investigating PR South
  
  # Create plot
  bathy_final_clamp <- clamp(bathy_final, lower = -50, upper = 0)
  plot(bathy_final,
       col = cmocean("deep")(100),
       ext = plot_extents,
       main = "Species locations over raster\n(Red = NA values, Black = Valid values)")
  
  # Add points - valid locations in black
  points(coords_valid$x_utm, coords_valid$y_utm,
         col = "black", pch = 16, cex = 0.5)
  
  # Add points - NA locations in red
  points(coords_na$x_utm, coords_na$y_utm,
         col = "red", pch = 16, cex = 0.5)
  
  # Add legend
  legend("topright",
         legend = c("Valid data", "NA values"),
         col = c("black", "red"),
         pch = 16,
         cex = 0.8)
  
  # Print summary
  cat("Total points:", nrow(site_data), "\n")
  cat("Valid points:", sum(!na_mask), "\n")
  cat("NA points:", sum(na_mask), "\n")
  cat("Percentage with NAs:", round(sum(na_mask)/nrow(site_data)*100, 1), "%\n")
  
  # Filter data for simple models
  model_data_filtered <- site_data %>%
    filter(depth_bathy >= -60)
  
  cat("Remaining datasets:", paste(unique(model_data_filtered$dataset), collapse = ", "), "\n")
  
  psu_with_pr <- model_data_filtered$PSU[grepl("_PR", model_data_filtered$PSU)]
  if(length(psu_with_pr) > 0) {
    cat("PSUs with '_PR' still remaining:", paste(unique(psu_with_pr), collapse = ", "), "\n")
  } else {
    cat("No PSUs with '_PR' remaining\n")
  }
  
  model_data <- model_data_filtered[complete.cases(model_data_filtered[, c("depth_bathy", "slope", "cover")]), ]
  
  # Check the new cover distribution
  summary(model_data$cover)
  prop_zeros_new <- sum(model_data$cover == 0) / nrow(model_data)
  cat("New proportion of zeros:", round(prop_zeros_new, 3), "\n")
  
  # Compare with original
  prop_zeros_orig <- sum(site_data$cover == 0, na.rm = TRUE) / nrow(site_data)
  cat("Original proportion of zeros:", round(prop_zeros_orig, 3), "\n")
  
  # Option 1: Tweedie distribution (good for zero-inflated continuous data)
  gam_tweedie <- gam(cover ~ s(depth_bathy) + s(slope),
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
  gam_beta <- gam(cover_prop ~ s(depth_bathy) + s(slope),
                  data = model_data[model_data$cover_prop > 0, ],  # exclude zeros
                  family = betar())
  
  # Option 3: Two-part model (hurdle model)
  # Part 1: Presence/absence
  model_data$present <- ifelse(model_data$cover > 0, 1, 0)
  gam_presence <- gam(present ~ s(depth_bathy) + s(slope),
                      data = model_data,
                      family = binomial())
  
  # Part 2: Abundance given presence
  gam_abundance <- gam(cover ~ s(depth_bathy) + s(slope),
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
  
  # Plot both relationships on one page
  par(mfrow = c(1, 2))
  plot(gam_tweedie, select = 1, main = "Bathymetry effect on coral cover",
       xlab = "Bathymetry (m)", ylab = "Smooth term")
  plot(gam_tweedie, select = 2, main = "Slope effect on coral cover",
       xlab = "Slope", ylab = "Smooth term")
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
  
  # Create prediction data for plotting
  bathy_range <- seq(min(model_data$depth_bathy, na.rm = TRUE),
                     max(model_data$depth_bathy, na.rm = TRUE),
                     length.out = 100)
  slope_range <- seq(min(model_data$slope, na.rm = TRUE),
                     max(model_data$slope, na.rm = TRUE),
                     length.out = 100)
  
  # Predictions for bathymetry (holding slope at median)
  pred_data_bathy <- data.frame(
    depth_bathy = bathy_range,
    slope = median(model_data$slope, na.rm = TRUE)
  )
  
  # Predictions for slope (holding bathymetry at median)
  pred_data_slope <- data.frame(
    depth_bathy = median(model_data$depth_bathy, na.rm = TRUE),
    slope = slope_range
  )
  
  # Get predictions with standard errors
  pred_bathy <- predict(gam_tweedie, pred_data_bathy, se.fit = TRUE, type = "response")
  pred_slope <- predict(gam_tweedie, pred_data_slope, se.fit = TRUE, type = "response")
  
  # Create plots
  par(mfrow = c(1, 2))
  
  # Bathymetry plot
  plot(pred_data_bathy$depth_bathy, pred_bathy$fit, type = "l",
       xlab = "Bathymetry (m)", ylab = "Predicted coral cover (%)",
       main = "Coral cover vs Bathymetry", lwd = 2, col = "blue")
  lines(pred_data_bathy$depth_bathy, pred_bathy$fit + 1.96*pred_bathy$se.fit, lty = 2, col = "blue")
  lines(pred_data_bathy$depth_bathy, pred_bathy$fit - 1.96*pred_bathy$se.fit, lty = 2, col = "blue")
  
  # Add raw data points
  points(model_data$depth_bathy, model_data$cover, pch = 16, cex = 0.3, col = "gray60")
  
  # Slope plot
  plot(pred_data_slope$slope, pred_slope$fit, type = "l",
       xlab = "Slope", ylab = "Predicted coral cover (%)",
       main = "Coral cover vs Slope", lwd = 2, col = "red")
  lines(pred_data_slope$slope, pred_slope$fit + 1.96*pred_slope$se.fit, lty = 2, col = "red")
  lines(pred_data_slope$slope, pred_slope$fit - 1.96*pred_slope$se.fit, lty = 2, col = "red")
  
  # Add raw data points
  points(model_data$slope, model_data$cover, pch = 16, cex = 0.3, col = "gray60")
  
  par(mfrow = c(1, 1))
  
  # Alternative ggplot version
  # Create prediction dataframes
  pred_bathy_df <- data.frame(
    depth_bathy = pred_data_bathy$depth_bathy,
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
  
  # Bathymetry plot with ggplot
  p1 <- ggplot(pred_bathy_df, aes(x = depth_bathy)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "blue") +
    geom_line(aes(y = fitted), color = "blue", size = 1) +
    geom_point(data = model_data, aes(x = depth_bathy, y = cover),
               alpha = 0.3, size = 0.5) +
    labs(x = "Bathymetry (m)", y = "Predicted coral cover (%)",
         title = "Coral cover vs Bathymetry") +
    theme_minimal()
  
  # Slope plot with ggplot
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
  
  ################################## COMPLEX SITE GAMs ##################################
  
  # Prepare data for complex models using PSU-averaged data
  site_model_data_filtered <- site_data %>%
    filter(depth_bathy >= -60)
  
  # Extract all environmental values at species locations
  site_species_env_complex <- terra::extract(env_complex,
                                             cbind(site_model_data_filtered$x_utm,
                                                   site_model_data_filtered$y_utm))
  
  # Add all variables to your filtered dataset
  site_model_data_filtered$depth_bathy <- site_species_env_complex$depth
  site_model_data_filtered$aspect <- site_species_env_complex$aspect
  site_model_data_filtered$slope <- site_species_env_complex$slope
  site_model_data_filtered$complexity <- site_species_env_complex$complexity
  site_model_data_filtered$TPI <- site_species_env_complex$TPI
  site_model_data_filtered$VRM <- site_species_env_complex$VRM
  site_model_data_filtered$planform_curv <- site_species_env_complex$planform_curv
  site_model_data_filtered$SAPA <- site_species_env_complex$SAPA
  site_model_data_filtered$max_Hsig <- site_species_env_complex$max_Hsig
  site_model_data_filtered$dir_at_max_hsig <- site_species_env_complex$dir_at_max_hsig
  site_model_data_filtered$mean_Hsig <- site_species_env_complex$mean_Hsig
  site_model_data_filtered$mean_SST <- site_species_env_complex$mean_SST
  site_model_data_filtered$range_SST <- site_species_env_complex$range_SST
  site_model_data_filtered$range_PAR <- site_species_env_complex$range_PAR
  site_model_data_filtered$mean_chla <- site_species_env_complex$mean_chla
  site_model_data_filtered$mean_kd490 <- site_species_env_complex$mean_kd490
  site_model_data_filtered$mean_spm <- site_species_env_complex$mean_spm
  site_model_data_filtered$dist_to_land <- site_species_env_complex$dist_to_land
  site_model_data_filtered$dist_to_deep <- site_species_env_complex$dist_to_deep
  site_model_data_filtered$max_BOV <- site_species_env_complex$max_BOV
  site_model_data_filtered$year <- year(site_model_data_filtered$date)
  
  # Create complete cases dataset with all variables
  site_model_data_complex <- site_model_data_filtered[complete.cases(site_model_data_filtered[, complexity_vars]), ]
  
  # Check how much data you have left
  cat("Site observations with all complexity variables:", nrow(site_model_data_complex), "\n")
  
  # Fit expanded GAM models
  site_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
                            s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                            mean_Hsig + s(range_PAR) + s(mean_chla) + s(dist_to_land) +
                            s(dist_to_deep) + s(max_BOV),
                          data = site_model_data_complex,
                          family = tw())
  
  # Check the results
  summary(site_gam_all_tweedie)
  AIC(site_gam_all_tweedie)
  gam.check(site_gam_all_tweedie)
  draw(site_gam_all_tweedie)
  
  # Plot the relationships
  plot(site_gam_all_tweedie, pages = 2)
  
  # Compare model performance
  cat("Simple model AIC:", AIC(gam_tweedie), "\n")
  cat("Complex model AIC:", AIC(site_gam_all_tweedie), "\n")
  
  # Check correlations between variables
  site_cor_matrix <- create_correlation_matrix(site_model_data_complex)
  print(round(site_cor_matrix, 2))
  
  # Two-part model with complexity
  site_model_data_complex$present <- ifelse(site_model_data_complex$cover > 0, 1, 0)
  
  site_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                     s(complexity, k = 12) +
                                     max_Hsig +
                                     s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                     s(dist_to_deep, k = 12),
                                   data = site_model_data_complex,
                                   family = binomial())
  
  site_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                      s(TPI) +
                                      s(mean_SST, k = 12) + s(dir_at_max_hsig, bs = 'cc') +
                                      s(range_PAR, k = 12) + s(dist_to_land) +
                                      s(dist_to_deep) + s(max_BOV),
                                    data = site_model_data_complex[site_model_data_complex$cover > 0, ],
                                    family = Gamma(link = "log"))
  
  summary(site_gam_presence_binom)
  summary(site_gam_abundance_gamma)
  AIC(site_gam_presence_binom)
  AIC(site_gam_abundance_gamma)
  
  gam.check(site_gam_presence_binom)
  gam.check(site_gam_abundance_gamma)
  
  draw(site_gam_presence_binom)
  draw(site_gam_abundance_gamma)
  
  concurvity(site_gam_presence_binom)
  concurvity(site_gam_abundance_gamma)
  
  # Create a reduced model with less correlated variables
  site_gam_reduced <- gam(cover ~ s(depth_bathy) + TPI + s(planform_curv),
                          data = site_model_data_complex,
                          family = tw())
  
  summary(site_gam_reduced)
  plot(site_gam_reduced, pages = 1)
  
  
  
  ################################## SPP DATA SETUP ##################################
  
  # Convert tibble to data.frame first
  species_df <- as.data.frame(combined_benthic_data_averaged)
  
  spp_data <- combined_benthic_data_averaged
  
  # Now create SpatVector
  species_coords <- vect(species_df,
                         geom = c("lon", "lat"),
                         crs = "EPSG:4326")  # WGS84 geographic
  
  # Transform to match your raster CRS
  species_coords_utm <- project(species_coords, crs(bathy_final))
  
  # Extract transformed coordinates
  utm_coords <- as.data.frame(geom(species_coords_utm)[, c("x", "y")])
  spp_data$x_utm <- utm_coords$x
  spp_data$y_utm <- utm_coords$y
  
  
  
  
  # Extract environmental values for simple variables
  species_env_values <- terra::extract(env_simple,
                                       cbind(spp_data$x_utm,
                                             spp_data$y_utm))
  
  
  
  
  # Add to dataframe
  spp_data$depth_bathy <- species_env_values$depth_bathy
  spp_data$slope <- species_env_values$slope
  
  # Add Y/N columns for bathymetry and slope presence
  spp_data$bathymetry_present <- ifelse(is.na(spp_data$depth_bathy), "N", "Y")
  spp_data$slope_present <- ifelse(is.na(spp_data$slope), "N", "Y")
  
  ################################## AGARICIA ##################################
  
  agariciids <- spp_data %>%
    filter(depth_bathy >= -60) %>%
    filter(grepl("Agaricia", spp))
  
  agaricia_model_data_filtered <- agariciids
  
  # Extract environmental data
  agaricia_species_env_complex <- extract_env_data(agariciids)
  
  # Add environmental variables
  agariciids <- add_env_variables(agariciids, agaricia_species_env_complex, agaricia_model_data_filtered)
  
  agaricia_model_data_complex <- agariciids
  
  # Check how much data you have left
  cat("Agaricia observations with all complexity variables:", nrow(agaricia_model_data_complex), "\n")
  
  # Fit expanded GAM models
  agaricia_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
                                s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                s(mean_Hsig) + s(range_PAR) + s(mean_chla) + s(dist_to_land) +
                                s(dist_to_deep) + s(max_BOV),
                              data = agaricia_model_data_complex, family = tw())
  
  # Check the results
  summary(agaricia_gam_all_tweedie)
  AIC(agaricia_gam_all_tweedie)
  gam.check(agaricia_gam_all_tweedie)
  plot(agaricia_gam_all_tweedie, pages = 2)
  draw(agaricia_gam_all_tweedie)
  
  # Create correlation matrix
  agaricia_cor_matrix <- create_correlation_matrix(agaricia_model_data_complex)
  print(round(agaricia_cor_matrix, 2))
  
  # Two-part model with complexity
  agaricia_model_data_complex$present <- ifelse(agaricia_model_data_complex$cover > 0, 1, 0)
  
  # tic()
  # agaricia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(range_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = agaricia_model_data_complex,
  #                             select = TRUE,
  #                             family = binomial())
  # toc()
  
  tic()
  agaricia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                       s(complexity) +
                                       s(dir_at_max_hsig, bs = 'cc') + mean_Hsig +
                                       s(mean_SST) + s(mean_kd490) +
                                       s(max_BOV) +
                                       s(dist_to_land),
                                     data = agaricia_model_data_complex,
                                     select = TRUE,
                                     family = binomial())
  toc()
  
  agaricia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(VRM) + s(aspect, bs = 'cc') +
                                         s(planform_curv) + s(max_Hsig) + s(mean_spm) +
                                         s(TPI, k = 12) + s(complexity, k = 12) +
                                         s(mean_SST, k = 20) + s(dir_at_max_hsig, bs = 'cc') +
                                         s(dist_to_land) +
                                         s(dist_to_deep) + s(max_BOV),
                                       data = agaricia_model_data_complex,
                                     select = TRUE,
                                       family = binomial())
  
  # agaricia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(range_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = agaricia_model_data_complex[agaricia_model_data_complex$cover > 0, ],
  #                              family = Gamma(link = "log"))
  
  agaricia_gam_abundance_gamma <- gam(cover ~ depth_bathy + s(TPI) + s(VRM) +
                                          complexity +
                                          mean_SST +
                                          s(dir_at_max_hsig, bs = 'cc') +
                                          mean_Hsig + s(mean_chla) +
                                          s(dist_to_deep, k = 12),
                                        data = agaricia_model_data_complex[agaricia_model_data_complex$cover > 0, ],
                                        select = TRUE,
                                        family = Gamma(link = "log"))
  
  summary(agaricia_gam_presence_binom)
  summary(agaricia_gam_abundance_gamma)
  AIC(agaricia_gam_presence_binom)
  AIC(agaricia_gam_abundance_gamma)
  
  draw(agaricia_gam_presence_binom)
  draw(agaricia_gam_abundance_gamma)
  
  # Check if any smooths are hitting k limits
  gam.check(agaricia_gam_presence_binom)
  gam.check(agaricia_gam_abundance_gamma)
  
  # Look at concurvity
  concurvity(agaricia_gam_presence_binom, full = TRUE)
  concurvity(agaricia_gam_abundance_gamma, full = TRUE)
  
  #AUC / ROC
  agaricia_fitted_data <- agaricia_gam_presence_binom$model
  agaricia_roc_curve <- roc(agaricia_fitted_data$present, 
                             fitted(agaricia_gam_presence_binom))
  auc(agaricia_roc_curve)
  plot(agaricia_roc_curve, main = "ROC Curve for Madracis Presence Model")
  
  
  agaricia_gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                              data = agaricia_model_data_complex,
                              family = tw())
  
  summary(agaricia_gam_reduced)
  plot(agaricia_gam_reduced, pages = 1)
  
  # Model summaries
  summary(agaricia_gam_all_tweedie)
  summary(agaricia_gam_presence_binom)
  summary(agaricia_gam_abundance_gamma)
  summary(agaricia_gam_reduced)
  
  AIC(agaricia_gam_all_tweedie)
  AIC(agaricia_gam_presence_binom)
  AIC(agaricia_gam_abundance_gamma)
  AIC(agaricia_gam_reduced)
  
  ################################## MADRACIS ##################################
  
  madracis <- spp_data %>%
    filter(depth_bathy >= -60) %>%
    filter(grepl("Madracis", spp))
  
  madracis_model_data_filtered <- madracis
  
  # Extract environmental data
  madracis_species_env_complex <- extract_env_data(madracis)
  
  # Add environmental variables
  madracis <- add_env_variables(madracis, madracis_species_env_complex, madracis_model_data_filtered)
  
  madracis_model_data_complex <- madracis
  
  # Check how much data you have left
  cat("Madracis observations with all complexity variables:", nrow(madracis_model_data_complex), "\n")
  
  # Fit expanded GAM models
  madracis_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
                                    s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                    s(mean_Hsig) + s(range_PAR) + s(mean_chla) + s(dist_to_land) +
                                    s(dist_to_deep) + s(max_BOV),
                                  data = madracis_model_data_complex, family = tw())
  
  # Check the results
  summary(madracis_gam_all_tweedie)
  AIC(madracis_gam_all_tweedie)
  gam.check(madracis_gam_all_tweedie)
  plot(madracis_gam_all_tweedie, pages = 2)
  draw(madracis_gam_all_tweedie)
  
  # Create correlation matrix
  madracis_cor_matrix <- create_correlation_matrix(madracis_model_data_complex)
  print(round(madracis_cor_matrix, 2))
  
  # Two-part model with complexity
  madracis_model_data_complex$present <- ifelse(madracis_model_data_complex$cover > 0, 1, 0)
  
  # tic()
  # madracis_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(range_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = madracis_model_data_complex,
  #                             # select = TRUE,
  #                             family = binomial())
  # toc()
  madracis_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                       s(TPI) +
                                       s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
                                       s(mean_SST) + s(range_PAR, k = 15) + s(mean_kd490),
                                     data = madracis_model_data_complex,
                                     select = TRUE,
                                     family = binomial())
  
  # madracis_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(range_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = madracis_model_data_complex[madracis_model_data_complex$cover > 0, ],
  #                               # select = TRUE,
  #                              family = Gamma(link = "log"))
  madracis_gam_abundance_gamma <- gam(cover ~ s(complexity) +
                                        s(mean_Hsig) +
                                        s(mean_chla) +
                                        s(dist_to_land),
                                      data = madracis_model_data_complex[madracis_model_data_complex$cover > 0, ],
                                      select = TRUE,
                                      family = Gamma(link = "log"))
  
  summary(madracis_gam_presence_binom)
  summary(madracis_gam_abundance_gamma)
  AIC(madracis_gam_presence_binom)
  AIC(madracis_gam_abundance_gamma)
  
  draw(madracis_gam_presence_binom)
  draw(madracis_gam_abundance_gamma)
  
  # Check if any smooths are hitting k limits
  gam.check(madracis_gam_presence_binom)
  gam.check(madracis_gam_abundance_gamma)
  
  # Look at concurvity
  concurvity(madracis_gam_presence_binom, full = TRUE)
  concurvity(madracis_gam_abundance_gamma, full = TRUE)
  
  #AUC / ROC
  madracis_fitted_data <- madracis_gam_presence_binom$model
  madracis_roc_curve <- roc(madracis_fitted_data$present, 
                            fitted(madracis_gam_presence_binom))
  auc(madracis_roc_curve)
  plot(madracis_roc_curve, main = "ROC Curve for Madracis Presence Model")
  
  
  madracis_gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                              data = madracis_model_data_complex,
                              family = tw())
  
  summary(madracis_gam_reduced)
  plot(madracis_gam_reduced, pages = 1)
  
  # Model summaries
  summary(madracis_gam_all_tweedie)
  summary(madracis_gam_presence_binom)
  summary(madracis_gam_abundance_gamma)
  summary(madracis_gam_reduced)
  
  AIC(madracis_gam_all_tweedie)
  AIC(madracis_gam_presence_binom)
  AIC(madracis_gam_abundance_gamma)
  AIC(madracis_gam_reduced)
  
  ################################## PORITES ##################################
  
  porites <- spp_data %>%
    filter(depth_bathy >= -60) %>%
    filter(grepl("Porites", spp))
  
  porites_model_data_filtered <- porites
  
  # Extract environmental data
  porites_species_env_complex <- extract_env_data(porites)
  
  # Add environmental variables
  porites <- add_env_variables(porites, porites_species_env_complex, porites_model_data_filtered)
  
  porites_model_data_complex <- porites
  
  # Check how much data you have left
  cat("Porites observations with all complexity variables:", nrow(porites_model_data_complex), "\n")
  
  # Fit expanded GAM models
  porites_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
                                   s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                   s(mean_Hsig) + s(range_PAR) + s(mean_chla) + s(dist_to_land) +
                                   s(dist_to_deep) + s(max_BOV),
                                 data = porites_model_data_complex, family = tw())
  
  # Check the results
  summary(porites_gam_all_tweedie)
  AIC(porites_gam_all_tweedie)
  gam.check(porites_gam_all_tweedie)
  plot(porites_gam_all_tweedie, pages = 2)
  draw(porites_gam_all_tweedie)
  
  # Create correlation matrix
  porites_cor_matrix <- create_correlation_matrix(porites_model_data_complex)
  print(round(porites_cor_matrix, 2))
  
  # Two-part model with complexity
  porites_model_data_complex$present <- ifelse(porites_model_data_complex$cover > 0, 1, 0)
  
  porites_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                      s(TPI) +
                                      s(dir_at_max_hsig, bs = 'cc') +
                                      s(mean_kd490) +
                                      s(dist_to_deep) + s(max_BOV),
                                    data = porites_model_data_complex,
                                    family = binomial())
  
  porites_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') +
                                       s(mean_kd490) +
                                       s(dist_to_deep),
                                     data = porites_model_data_complex[porites_model_data_complex$cover > 0, ],
                                     family = Gamma(link = "log"))
  
  summary(porites_gam_presence_binom)
  summary(porites_gam_abundance_gamma)
  AIC(porites_gam_presence_binom)
  AIC(porites_gam_abundance_gamma)
  
  draw(porites_gam_presence_binom)
  draw(porites_gam_abundance_gamma)
  
  # Check if any smooths are hitting k limits
  gam.check(porites_gam_presence_binom)
  gam.check(porites_gam_abundance_gamma)
  
  # Look at concurvity
  concurvity(porites_gam_presence_binom, full = TRUE)
  concurvity(porites_gam_abundance_gamma, full = TRUE)
  
  porites_gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                             data = porites_model_data_complex,
                             family = tw())
  
  summary(porites_gam_reduced)
  plot(porites_gam_reduced, pages = 1)
  
  # Model summaries
  summary(porites_gam_all_tweedie)
  summary(porites_gam_presence_binom)
  summary(porites_gam_abundance_gamma)
  summary(porites_gam_reduced)
  
  AIC(porites_gam_all_tweedie)
  AIC(porites_gam_presence_binom)
  AIC(porites_gam_abundance_gamma)
  AIC(porites_gam_reduced)
  
  ################################## SIDERASTREA ##################################
  
  siderastrea <- spp_data %>%
    filter(depth_bathy >= -60) %>%
    filter(grepl("Siderastrea", spp))
  
  siderastrea_model_data_filtered <- siderastrea
  
  # Extract environmental data
  siderastrea_species_env_complex <- extract_env_data(siderastrea)
  
  # Add environmental variables
  siderastrea <- add_env_variables(siderastrea, siderastrea_species_env_complex, siderastrea_model_data_filtered)
  
  siderastrea_model_data_complex <- siderastrea
  
  # Check how much data you have left
  cat("Siderastrea observations with all complexity variables:", nrow(siderastrea_model_data_complex), "\n")
  
  # Fit expanded GAM models
  siderastrea_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
                                       s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                       s(mean_Hsig) + s(range_PAR) + s(mean_chla) + s(dist_to_land) +
                                       s(dist_to_deep) + s(max_BOV),
                                     data = siderastrea_model_data_complex, family = tw())
  
  # Check the results
  summary(siderastrea_gam_all_tweedie)
  AIC(siderastrea_gam_all_tweedie)
  gam.check(siderastrea_gam_all_tweedie)
  plot(siderastrea_gam_all_tweedie, pages = 2)
  draw(siderastrea_gam_all_tweedie)
  
  # Create correlation matrix
  siderastrea_cor_matrix <- create_correlation_matrix(siderastrea_model_data_complex)
  print(round(siderastrea_cor_matrix, 2))
  
  # Two-part model with complexity
  siderastrea_model_data_complex$present <- ifelse(siderastrea_model_data_complex$cover > 0, 1, 0)
  
  siderastrea_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                          s(complexity) +
                                          s(dir_at_max_hsig, bs = 'cc') +
                                          s(mean_kd490) +
                                          s(dist_to_deep, k = 12) + s(max_BOV),
                                        data = siderastrea_model_data_complex,
                                        family = binomial())
  
  siderastrea_gam_abundance_gamma <- gam(cover ~ s(max_Hsig) + s(dir_at_max_hsig, k = 12, bs = 'cc') +
                                           s(mean_kd490) +
                                           s(dist_to_deep, k = 12),
                                         data = siderastrea_model_data_complex[siderastrea_model_data_complex$cover > 0, ],
                                         family = Gamma(link = "log"))
  
  summary(siderastrea_gam_presence_binom)
  summary(siderastrea_gam_abundance_gamma)
  AIC(siderastrea_gam_presence_binom)
  AIC(siderastrea_gam_abundance_gamma)
  
  draw(siderastrea_gam_presence_binom)
  draw(siderastrea_gam_abundance_gamma)
  
  # Check if any smooths are hitting k limits
  gam.check(siderastrea_gam_presence_binom)
  gam.check(siderastrea_gam_abundance_gamma)
  
  # Look at concurvity
  concurvity(siderastrea_gam_presence_binom, full = TRUE)
  concurvity(siderastrea_gam_abundance_gamma, full = TRUE)
  
  siderastrea_gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                                 data = siderastrea_model_data_complex,
                                 family = tw())
  
  summary(siderastrea_gam_reduced)
  plot(siderastrea_gam_reduced, pages = 1)
  
  # Model summaries
  summary(siderastrea_gam_all_tweedie)
  summary(siderastrea_gam_presence_binom)
  summary(siderastrea_gam_abundance_gamma)
  summary(siderastrea_gam_reduced)
  
  AIC(siderastrea_gam_all_tweedie)
  AIC(siderastrea_gam_presence_binom)
  AIC(siderastrea_gam_abundance_gamma)
  AIC(siderastrea_gam_reduced)
  
  ################################## MONTASTRAEA ##################################
  
  montastraea <- spp_data %>%
    filter(depth_bathy >= -60) %>%
    filter(grepl("Montastraea", spp))
  
  montastraea_model_data_filtered <- montastraea
  
  # Extract environmental data
  montastraea_species_env_complex <- extract_env_data(montastraea)
  
  # Add environmental variables
  montastraea <- add_env_variables(montastraea, montastraea_species_env_complex, montastraea_model_data_filtered)
  
  montastraea_model_data_complex <- montastraea
  
  # Check how much data you have left
  cat("Montastraea observations with all complexity variables:", nrow(montastraea_model_data_complex), "\n")
  
  # Fit expanded GAM models
  montastraea_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
                                       s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                       s(mean_Hsig) + s(range_PAR) + s(mean_chla) + s(dist_to_land) +
                                       s(dist_to_deep) + s(max_BOV),
                                     data = montastraea_model_data_complex, family = tw())
  
  # Check the results
  summary(montastraea_gam_all_tweedie)
  AIC(montastraea_gam_all_tweedie)
  gam.check(montastraea_gam_all_tweedie)
  plot(montastraea_gam_all_tweedie, pages = 2)
  draw(montastraea_gam_all_tweedie)
  
  # Create correlation matrix
  montastraea_cor_matrix <- create_correlation_matrix(montastraea_model_data_complex)
  print(round(montastraea_cor_matrix, 2))
  
  # Two-part model with complexity
  montastraea_model_data_complex$present <- ifelse(montastraea_model_data_complex$cover > 0, 1, 0)
  
  montastraea_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                          complexity +
                                          s(dir_at_max_hsig, k = 12, bs = 'cc') +
                                          s(mean_SST) + s(mean_kd490) +
                                          s(dist_to_land, k = 15),
                                        data = montastraea_model_data_complex,
                                        family = binomial())
  
  montastraea_gam_abundance_gamma <- gam(cover ~ depth_bathy +
                                           s(VRM) +
                                           max_Hsig + s(mean_Hsig) +
                                           s(mean_SST, k = 12) + mean_kd490 +
                                           dist_to_deep + s(max_BOV, k = 12),
                                         data = montastraea_model_data_complex[montastraea_model_data_complex$cover > 0, ],
                                         family = Gamma(link = "log"))
  
  summary(montastraea_gam_presence_binom)
  summary(montastraea_gam_abundance_gamma)
  AIC(montastraea_gam_presence_binom)
  AIC(montastraea_gam_abundance_gamma)
  
  draw(montastraea_gam_presence_binom)
  draw(montastraea_gam_abundance_gamma)
  
  # Check if any smooths are hitting k limits
  gam.check(montastraea_gam_presence_binom)
  gam.check(montastraea_gam_abundance_gamma)
  
  # Look at concurvity
  concurvity(montastraea_gam_presence_binom, full = TRUE)
  concurvity(montastraea_gam_abundance_gamma, full = TRUE)
  
  cor(montastraea_model_data_complex[c("depth_bathy", "VRM", "max_Hsig", "mean_Hsig",
                                       "mean_SST", "mean_kd490", "dist_to_deep", "max_BOV")], use = "complete.obs")
  
  montastraea_gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                                 data = montastraea_model_data_complex,
                                 family = tw())
  
  summary(montastraea_gam_reduced)
  plot(montastraea_gam_reduced, pages = 1)
  
  # Model summaries
  summary(montastraea_gam_all_tweedie)
  summary(montastraea_gam_presence_binom)
  summary(montastraea_gam_abundance_gamma)
  summary(montastraea_gam_reduced)
  
  AIC(montastraea_gam_all_tweedie)
  AIC(montastraea_gam_presence_binom)
  AIC(montastraea_gam_abundance_gamma)
  AIC(montastraea_gam_reduced)
  
  ################################## ORBICELLA ##################################
  
  orbicellids <- spp_data %>%
    filter(depth_bathy >= -60) %>%
    filter(grepl("Orbicella", spp))
  
  orbicella_model_data_filtered <- orbicellids
  
  # Extract environmental data
  orbicella_species_env_complex <- extract_env_data(orbicellids)
  
  # Add environmental variables
  orbicellids <- add_env_variables(orbicellids, orbicella_species_env_complex, orbicella_model_data_filtered)
  
  orbicella_model_data_complex <- orbicellids
  
  # Check how much data you have left
  cat("Orbicella observations with all complexity variables:", nrow(orbicella_model_data_complex), "\n")
  
  # Fit expanded GAM models
  orbicella_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
                                     s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                     s(mean_Hsig) + s(range_PAR) + s(mean_chla) + s(dist_to_land) +
                                     s(dist_to_deep) + s(max_BOV),
                                   data = orbicella_model_data_complex, family = tw())
  
  # Check the results
  summary(orbicella_gam_all_tweedie)
  AIC(orbicella_gam_all_tweedie)
  gam.check(orbicella_gam_all_tweedie)
  plot(orbicella_gam_all_tweedie, pages = 2)
  draw(orbicella_gam_all_tweedie)
  
  # Create correlation matrix
  orbicella_cor_matrix <- create_correlation_matrix(orbicella_model_data_complex)
  print(round(orbicella_cor_matrix, 2))
  
  # Two-part model with complexity
  #
  #depth_bathy, aspect, slope, complexity, TPI, VRM, planform_curv, SAPA, max_Hsig,
  #   dir_at_max_Hsig, mean_Hsig, mean_SST, range_SST, range_PAR, mean_chla, mean_kd490, mean_spm,
  #   dist_to_land, dist_to_deep, max_BOV, year
  orbicella_model_data_complex$present <- ifelse(orbicella_model_data_complex$cover > 0, 1, 0)
  
  # tic()
  # orbicella_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(range_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = orbicella_model_data_complex,
  #                             select = TRUE,
  #                             family = binomial())
  # toc()
  orbicella_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(TPI) + s(complexity) +
                                        s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                        s(mean_Hsig) + s(range_PAR) + s(dist_to_land) +
                                        s(dist_to_deep) + s(max_BOV),
                                      data = orbicella_model_data_complex,
                                      # select = TRUE,
                                      family = binomial())
  
  # gam_abundance_complex <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(range_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                              data = model_data_complex[model_data_complex$cover > 0, ],
  #                              family = Gamma(link = "log"))
  orbicella_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(TPI) + complexity +
                                         s(range_SST, k = 12) + s(mean_SST, k = 12) +
                                         s(mean_Hsig, k = 12) + mean_chla +
                                         s(dist_to_deep, k = 12) +
                                         s(dist_to_land) + s(max_BOV),
                                       data = orbicella_model_data_complex[orbicella_model_data_complex$cover > 0, ],
                                       select = TRUE,
                                       family = Gamma(link = "log"))
  
  summary(orbicella_gam_presence_binom)
  summary(orbicella_gam_abundance_gamma)
  AIC(orbicella_gam_presence_binom)
  AIC(orbicella_gam_abundance_gamma)
  
  draw(orbicella_gam_presence_binom)
  draw(orbicella_gam_abundance_gamma)
  
  # Check if any smooths are hitting k limits
  gam.check(orbicella_gam_presence_binom)
  gam.check(orbicella_gam_abundance_gamma)
  
  # Look at concurvity
  concurvity(orbicella_gam_presence_binom, full = TRUE)
  concurvity(orbicella_gam_abundance_gamma, full = TRUE)
  
  #AUC / ROC
  orbicella_fitted_data <- orbicella_gam_presence_binom$model
  orbicella_roc_curve <- roc(orbicella_fitted_data$present, 
                            fitted(orbicella_gam_presence_binom))
  auc(orbicella_roc_curve)
  plot(orbicella_roc_curve, main = "ROC Curve for Madracis Presence Model")
  
  
  orbicella_gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                               data = orbicella_model_data_complex,
                               family = tw())
  
  summary(orbicella_gam_reduced)
  plot(orbicella_gam_reduced, pages = 1)
  
  # Model summaries
  summary(orbicella_gam_all_tweedie)
  summary(orbicella_gam_presence_binom)
  summary(orbicella_gam_abundance_gamma)
  summary(orbicella_gam_reduced)
  
  AIC(orbicella_gam_all_tweedie)
  AIC(orbicella_gam_presence_binom)
  AIC(orbicella_gam_abundance_gamma)
  AIC(orbicella_gam_reduced)
  
  ################################## SOLENASTREA ##################################
  
  # NOTE - may be too low sample size
  
  solenastrea <- spp_data %>%
    filter(depth_bathy >= -60) %>%
    filter(grepl("Solenastrea", spp))
  
  solenastrea_model_data_filtered <- solenastrea
  
  # Extract environmental data
  solenastrea_species_env_complex <- extract_env_data(solenastrea)
  
  # Add environmental variables
  solenastrea <- add_env_variables(solenastrea, solenastrea_species_env_complex, solenastrea_model_data_filtered)
  
  solenastrea_model_data_complex <- solenastrea
  
  # Check how much data you have left
  cat("Solenastrea observations with all complexity variables:", nrow(solenastrea_model_data_complex), "\n")
  
  # Fit expanded GAM models
  solenastrea_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
                                    s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                    s(mean_Hsig) + s(range_PAR) + s(mean_chla) + s(dist_to_land) +
                                    s(dist_to_deep) + s(max_BOV),
                                  data = solenastrea_model_data_complex, family = tw())
  
  # Check the results
  summary(solenastrea_gam_all_tweedie)
  AIC(solenastrea_gam_all_tweedie)
  gam.check(solenastrea_gam_all_tweedie)
  plot(solenastrea_gam_all_tweedie, pages = 2)
  draw(solenastrea_gam_all_tweedie)
  
  # Create correlation matrix
  solenastrea_cor_matrix <- create_correlation_matrix(solenastrea_model_data_complex)
  print(round(solenastrea_cor_matrix, 2))
  
  # Two-part model with complexity
  solenastrea_model_data_complex$present <- ifelse(solenastrea_model_data_complex$cover > 0, 1, 0)
  
  tic()
  solenastrea_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                s(slope) +
                                s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
                                s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
                                s(mean_SST) + s(range_PAR) + s(mean_chla) + s(mean_kd490) +
                                s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
                                s(range_SST) +
                                s(dist_to_land),
                              data = solenastrea_model_data_complex,
                              select = TRUE,
                              family = binomial())
  toc()
  solenastrea_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                          s(slope) +
                                          s(mean_SST, k = 15) +
                                          s(mean_spm, k = 15) +
                                          s(range_SST),
                                        data = solenastrea_model_data_complex,
                                        select = TRUE,
                                        family = binomial())
  
  solenastrea_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                s(slope) +
                                s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
                                s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
                                s(mean_SST) + s(range_PAR) + s(mean_chla) + s(mean_kd490) +
                                s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
                                s(range_SST) +
                                s(dist_to_land),
                                data = solenastrea_model_data_complex[solenastrea_model_data_complex$cover > 0, ],
                                select = TRUE,
                               family = Gamma(link = "log"))
  
  solenastrea_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                           s(slope) +
                                           s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
                                           s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
                                           s(mean_SST) + s(range_PAR) + s(mean_chla) + s(mean_kd490) +
                                           s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
                                           s(range_SST) +
                                           s(dist_to_land),
                                         data = solenastrea_model_data_complex[solenastrea_model_data_complex$cover > 0, ],
                                         # select = TRUE,
                                         family = Gamma(link = "log"))
  
  summary(solenastrea_gam_presence_binom)
  summary(solenastrea_gam_abundance_gamma)
  AIC(solenastrea_gam_presence_binom)
  AIC(solenastrea_gam_abundance_gamma)
  
  draw(solenastrea_gam_presence_binom)
  draw(solenastrea_gam_abundance_gamma)
  
  # Check if any smooths are hitting k limits
  gam.check(solenastrea_gam_presence_binom)
  gam.check(solenastrea_gam_abundance_gamma)
  
  # Look at concurvity
  concurvity(solenastrea_gam_presence_binom, full = TRUE)
  concurvity(solenastrea_gam_abundance_gamma, full = TRUE)
  
  #AUC / ROC
  solenastrea_fitted_data <- solenastrea_gam_presence_binom$model
  solenastrea_roc_curve <- roc(solenastrea_fitted_data$present, 
                            fitted(solenastrea_gam_presence_binom))
  auc(solenastrea_roc_curve)
  plot(solenastrea_roc_curve)
  
  
  solenastrea_gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                              data = solenastrea_model_data_complex,
                              family = tw())
  
  summary(solenastrea_gam_reduced)
  plot(solenastrea_gam_reduced, pages = 1)
  
  # Model summaries
  summary(solenastrea_gam_all_tweedie)
  summary(solenastrea_gam_presence_binom)
  summary(solenastrea_gam_abundance_gamma)
  summary(solenastrea_gam_reduced)
  
  AIC(solenastrea_gam_all_tweedie)
  AIC(solenastrea_gam_presence_binom)
  AIC(solenastrea_gam_abundance_gamma)
  AIC(solenastrea_gam_reduced)
  
  ################################## COLPOPHYLLIA ##################################
  
  colpophyllia <- spp_data %>%
    filter(depth_bathy >= -60) %>%
    filter(grepl("Colpophyllia", spp))
  
  colpophyllia_model_data_filtered <- colpophyllia
  
  # Extract environmental data
  colpophyllia_species_env_complex <- extract_env_data(colpophyllia)
  
  # Add environmental variables
  colpophyllia <- add_env_variables(colpophyllia, colpophyllia_species_env_complex, colpophyllia_model_data_filtered)
  
  colpophyllia_model_data_complex <- colpophyllia
  
  # Check how much data you have left
  cat("Colpophyllia observations with all complexity variables:", nrow(colpophyllia_model_data_complex), "\n")
  
  # Fit expanded GAM models
  colpophyllia_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
                                        s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                        s(mean_Hsig) + s(range_PAR) + s(mean_chla) + s(dist_to_land) +
                                        s(dist_to_deep) + s(max_BOV),
                                      data = colpophyllia_model_data_complex, family = tw())
  
  # Check the results
  summary(colpophyllia_gam_all_tweedie)
  AIC(colpophyllia_gam_all_tweedie)
  gam.check(colpophyllia_gam_all_tweedie)
  plot(colpophyllia_gam_all_tweedie, pages = 2)
  draw(colpophyllia_gam_all_tweedie)
  
  # Create correlation matrix
  colpophyllia_cor_matrix <- create_correlation_matrix(colpophyllia_model_data_complex)
  print(round(colpophyllia_cor_matrix, 2))
  
  # Two-part model with complexity
  colpophyllia_model_data_complex$present <- ifelse(colpophyllia_model_data_complex$cover > 0, 1, 0)
  
  colpophyllia_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                           SAPA +
                                           mean_SST +
                                           mean_spm + s(max_BOV) +
                                           s(range_SST),
                                         data = colpophyllia_model_data_complex,
                                         family = binomial())
  
  colpophyllia_gam_abundance_gamma <- gam(cover ~ s(slope) + TPI + VRM + planform_curv +
                                            mean_spm +
                                            s(dist_to_deep),
                                          data = colpophyllia_model_data_complex[colpophyllia_model_data_complex$cover > 0, ],
                                          family = Gamma(link = "log"))
  
  summary(colpophyllia_gam_presence_binom)
  summary(colpophyllia_gam_abundance_gamma)
  AIC(colpophyllia_gam_presence_binom)
  AIC(colpophyllia_gam_abundance_gamma)
  
  draw(colpophyllia_gam_presence_binom)
  draw(colpophyllia_gam_abundance_gamma)
  
  # Check if any smooths are hitting k limits
  gam.check(colpophyllia_gam_presence_binom)
  gam.check(colpophyllia_gam_abundance_gamma)
  
  # Look at concurvity
  concurvity(colpophyllia_gam_presence_binom, full = TRUE)
  concurvity(colpophyllia_gam_abundance_gamma, full = TRUE)
  
  cor(colpophyllia_model_data_complex[c("SAPA", "mean_SST", "mean_spm", "max_BOV", "depth")], use = "complete.obs")
  
  colpophyllia_gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                                  data = colpophyllia_model_data_complex,
                                  family = tw())
  
  summary(colpophyllia_gam_reduced)
  plot(colpophyllia_gam_reduced, pages = 1)
  
  # Model summaries
  summary(colpophyllia_gam_all_tweedie)
  summary(colpophyllia_gam_presence_binom)
  summary(colpophyllia_gam_abundance_gamma)
  summary(colpophyllia_gam_reduced)
  
  AIC(colpophyllia_gam_all_tweedie)
  AIC(colpophyllia_gam_presence_binom)
  AIC(colpophyllia_gam_abundance_gamma)
  AIC(colpophyllia_gam_reduced)
  
  
  
  ################################## DENDROGYRA ##################################
  
  dendrogyra <- spp_data %>%
    filter(depth_bathy >= -60) %>%
    filter(grepl("Dendrogyra", spp))
  
  dendrogyra_model_data_filtered <- dendrogyra
  
  # Extract environmental data
  dendrogyra_species_env_complex <- extract_env_data(dendrogyra)
  
  # Add environmental variables
  dendrogyra <- add_env_variables(dendrogyra, dendrogyra_species_env_complex, dendrogyra_model_data_filtered)
  
  dendrogyra_model_data_complex <- dendrogyra
  
  # Check how much data you have left
  cat("Dendrogyra observations with all complexity variables:", nrow(dendrogyra_model_data_complex), "\n")
  
  # Fit expanded GAM models
  dendrogyra_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
                                  s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                  s(mean_Hsig) + s(range_PAR) + s(mean_chla) + s(dist_to_land) +
                                  s(dist_to_deep) + s(max_BOV),
                                data = dendrogyra_model_data_complex, family = tw())
  
  # Check the results
  summary(dendrogyra_gam_all_tweedie)
  AIC(dendrogyra_gam_all_tweedie)
  gam.check(dendrogyra_gam_all_tweedie)
  plot(dendrogyra_gam_all_tweedie, pages = 2)
  draw(dendrogyra_gam_all_tweedie)
  
  # Create correlation matrix
  dendrogyra_cor_matrix <- create_correlation_matrix(dendrogyra_model_data_complex)
  print(round(dendrogyra_cor_matrix, 2))
  
  # Two-part model with complexity
  dendrogyra_model_data_complex$present <- ifelse(dendrogyra_model_data_complex$cover > 0, 1, 0)
  
  dendrogyra_gam_presence_binom <- gam(present ~ s(depth_bathy) + VRM +
                                           slope +
                                           SAPA + mean_spm +
                                           complexity +
                                           s(mean_SST, k = 12) +
                                           s(dist_to_land, k = 12) +
                                           dist_to_deep + max_BOV,
                                         data = dendrogyra_model_data_complex,
                                         family = binomial())
  
  dendrogyra_gam_abundance_gamma <- gam(cover ~ s(slope) + TPI + VRM + planform_curv +
                                            mean_spm +
                                            s(dist_to_deep),
                                          data = dendrogyra_model_data_complex[dendrogyra_model_data_complex$cover > 0, ],
                                          family = Gamma(link = "log"))
  
  summary(dendrogyra_gam_presence_binom)
  summary(dendrogyra_gam_abundance_gamma)
  AIC(dendrogyra_gam_presence_binom)
  AIC(dendrogyra_gam_abundance_gamma)
  
  draw(dendrogyra_gam_presence_binom)
  draw(dendrogyra_gam_abundance_gamma)
  
  # Check if any smooths are hitting k limits
  gam.check(dendrogyra_gam_presence_binom)
  gam.check(dendrogyra_gam_abundance_gamma)
  
  # Look at concurvity
  concurvity(dendrogyra_gam_presence_binom, full = TRUE)
  concurvity(dendrogyra_gam_abundance_gamma, full = TRUE)
  
  dendrogyra_gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                                data = dendrogyra_model_data_complex,
                                family = tw())
  
  summary(dendrogyra_gam_reduced)
  plot(dendrogyra_gam_reduced, pages = 1)
  
  # Model summaries
  summary(dendrogyra_gam_all_tweedie)
  summary(dendrogyra_gam_presence_binom)
  summary(dendrogyra_gam_abundance_gamma)
  summary(dendrogyra_gam_reduced)
  
  AIC(dendrogyra_gam_all_tweedie)
  AIC(dendrogyra_gam_presence_binom)
  AIC(dendrogyra_gam_abundance_gamma)
  AIC(dendrogyra_gam_reduced)
  

  ################################## DICHOCOENIA ##################################
  
  dichocoenia <- spp_data %>%
    filter(depth_bathy >= -60) %>%
    filter(grepl("Dichocoenia", spp))
  
  dichocoenia_model_data_filtered <- dichocoenia
  
  # Extract environmental data
  dichocoenia_species_env_complex <- extract_env_data(dichocoenia)
  
  # Add environmental variables
  dichocoenia <- add_env_variables(dichocoenia, dichocoenia_species_env_complex, dichocoenia_model_data_filtered)
  
  dichocoenia_model_data_complex <- dichocoenia
  
  # Check how much data you have left
  cat("Dichocoenia observations with all complexity variables:", nrow(dichocoenia_model_data_complex), "\n")
  
  # Fit expanded GAM models
  dichocoenia_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
                                    s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
                                    s(mean_Hsig) + s(range_PAR) + s(mean_chla) + s(dist_to_land) +
                                    s(dist_to_deep) + s(max_BOV),
                                  data = dichocoenia_model_data_complex, family = tw())
  
  # Check the results
  summary(dichocoenia_gam_all_tweedie)
  AIC(dichocoenia_gam_all_tweedie)
  gam.check(dichocoenia_gam_all_tweedie)
  plot(dichocoenia_gam_all_tweedie, pages = 2)
  draw(dichocoenia_gam_all_tweedie)
  
  # Create correlation matrix
  dichocoenia_cor_matrix <- create_correlation_matrix(dichocoenia_model_data_complex)
  print(round(dichocoenia_cor_matrix, 2))
  
  # Two-part model with complexity
  dichocoenia_model_data_complex$present <- ifelse(dichocoenia_model_data_complex$cover > 0, 1, 0)
  
  tic()
  dichocoenia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                s(slope) +
                                s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
                                s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
                                s(mean_SST) + s(range_PAR) + s(mean_chla) + s(mean_kd490) +
                                s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
                                s(range_SST) +
                                s(dist_to_land),
                              data = dichocoenia_model_data_complex,
                              select = TRUE,
                              family = binomial())
  toc()
  dichocoenia_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                          s(mean_SST) + s(range_PAR, k = 12) + s(dist_to_deep) +
                                          s(dist_to_land),
                                        data = dichocoenia_model_data_complex,
                                        select = TRUE,
                                        family = binomial())
  
  # dichocoenia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(range_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = dichocoenia_model_data_complex[dichocoenia_model_data_complex$cover > 0, ],
  #                               select = TRUE,
  #                              family = Gamma(link = "log"))
  dichocoenia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) +
                                           s(planform_curv) +
                                           s(range_PAR) + s(dist_to_land),
                                         data = dichocoenia_model_data_complex[dichocoenia_model_data_complex$cover > 0, ],
                                         select = TRUE,
                                         family = Gamma(link = "log"))
  
  summary(dichocoenia_gam_presence_binom)
  summary(dichocoenia_gam_abundance_gamma)
  AIC(dichocoenia_gam_presence_binom)
  AIC(dichocoenia_gam_abundance_gamma)
  
  draw(dichocoenia_gam_presence_binom)
  draw(dichocoenia_gam_abundance_gamma)
  
  # Check if any smooths are hitting k limits
  gam.check(dichocoenia_gam_presence_binom)
  gam.check(dichocoenia_gam_abundance_gamma)
  
  # Look at concurvity
  concurvity(dichocoenia_gam_presence_binom, full = TRUE)
  concurvity(dichocoenia_gam_abundance_gamma, full = TRUE)
  
  #AUC / ROC
  dichocoenia_fitted_data <- dichocoenia_gam_presence_binom$model
  dichocoenia_roc_curve <- roc(dichocoenia_fitted_data$present, 
                            fitted(dichocoenia_gam_presence_binom))
  auc(dichocoenia_roc_curve)
  plot(dichocoenia_roc_curve, main = "ROC Curve for Dichocoenia Presence Model")
  
  
  dichocoenia_gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                              data = dichocoenia_model_data_complex,
                              family = tw())
  
  summary(dichocoenia_gam_reduced)
  plot(dichocoenia_gam_reduced, pages = 1)
  
  # Model summaries
  summary(dichocoenia_gam_all_tweedie)
  summary(dichocoenia_gam_presence_binom)
  summary(dichocoenia_gam_abundance_gamma)
  summary(dichocoenia_gam_reduced)
  
  AIC(dichocoenia_gam_all_tweedie)
  AIC(dichocoenia_gam_presence_binom)
  AIC(dichocoenia_gam_abundance_gamma)
  AIC(dichocoenia_gam_reduced)
  
  
  ################################## DIPLORIA ##################################
  
  
  ################################## EUSMILIA ##################################
  
  
  ################################## MEANDRINA ##################################
  
  
  
  ################################## MYCETOPHYLLIA ##################################
  
  
  
  ################################## PSEUDODIPLORIA ##################################
  
  
  
  ################################## RARE HS ##################################
  
  
  
  # stopping point - 8 sep 2025
  #   - need to handle some predictor outliers
  #   - group rare HS, consider grouping others like solenastrea too
  #   - consider mean BOV & wave dir
  #   - look at some maps!!