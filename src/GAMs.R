    
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
  
  library(pROC)
  
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

  ################################## simple site GAMs ##################################
  
  # Convert tibble to data.frame first
  species_df <- as.data.frame(combined_benthic_data_averaged_psu)
  
  spp_data = combined_benthic_data_averaged_psu
  
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
  
  # Create simple environmental stack
  env_simple <- c(bathy_final, slope_terra)
  names(env_simple) <- c("depth_bathy", "slope")

  # Extract environmental values at species locations
  species_env_values <- terra::extract(env_simple,
                                       cbind(spp_data$x_utm,
                                             spp_data$y_utm))

  # Add to dataframe
  spp_data$depth_bathy <- species_env_values$depth_bathy
  spp_data$slope <- species_env_values$slope

  # Add Y/N columns for bathymetry and slope presence
  spp_data$bathymetry_present <- ifelse(is.na(spp_data$depth_bathy), "N", "Y")
  spp_data$slope_present <- ifelse(is.na(spp_data$slope), "N", "Y")

  # Check how many NAs you have
  # NOTE - this is actually kind of a lot.
  sum(is.na(spp_data$depth_bathy))
  sum(is.na(spp_data$slope))
  
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
  na_mask <- is.na(spp_data$depth_bathy) |
    is.na(spp_data$slope)
  #
  # Get coordinates for plotting
  coords_na <- spp_data[na_mask, c("x_utm", "y_utm")]
  coords_valid <- spp_data[!na_mask, c("x_utm", "y_utm")]
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
  plot_extents = ext(240000, 280000, 2000000, 2040000) #for investigating Mona Island
  # plot_extents = ext(210000, 260000, 1995000, 2050000) #for investigating PR East
  # plot_extents = ext(55000, 230000, 2030000, 2058000) #for investigating PR North
  # plot_extents = ext(20000, 80000, 1970000, 2060000) #for investigating PR West
  # plot_extents = ext(52000, 230000, 1970000, 2020000) #for investigating PR South

  # raster
  bathy_final_clamp <- clamp(bathy_final, lower = -50, upper = 0)
  plot(bathy_final,
  # plot(mean_hsig_raster, #for displaying loss if including 1 km wave data
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
  cat("Total points:", nrow(spp_data), "\n")
  cat("Valid points:", sum(!na_mask), "\n")
  cat("NA points:", sum(na_mask), "\n")
  cat("Percentage with NAs:", round(sum(na_mask)/nrow(spp_data)*100, 1), "%\n")

  #TEST - taking out PR & NODICE entirely
  # also, filter to just 2017-2018; remove anything below 50 m depth
  # - wait actually not doing that right now because there technically isn't date info for
  #     the collated TCRMP data as-is!
  model_data_filtered <- spp_data %>%
    # filter(!dataset %in% c("PRCRMP", "NODICE")) %>%  # Remove PRCRMP and NODICE datasets
    # filter(year(date) >= 2017) %>%
    # filter(depth_bathy >= -50) %>%
    filter(!grepl("_PR", PSU))  # Remove any PSU containing '_PR'

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
  prop_zeros_orig <- sum(spp_data$cover == 0, na.rm = TRUE) / nrow(spp_data)
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

  # Create prediction dataframes as above, then:
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

  # Bathymetry plot
  p1 <- ggplot(pred_bathy_df, aes(x = depth_bathy)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "blue") +
    geom_line(aes(y = fitted), color = "blue", size = 1) +
    geom_point(data = model_data, aes(x = depth_bathy, y = cover),
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


  ################################## complex site GAMs ##################################
  
  # Extract all the complexity variables at species locations
  env_complex <- c(bathy_final, aspect_terra, slope_terra, slopeofslope_terra, TPI_terra, VRM,
                   planformcurv_multiscale, SAPA, max_hsig_raster, mean_dir_raster,
                   mean_hsig_raster, mean_sst_raster, range_sst_raster)
  
  names(env_complex) <- c("depth", "aspect", "slope", "complexity", "TPI", "VRM", "planform_curv",
                          "SAPA", "max_Hsig", "mean_dir", "mean_Hsig", "mean_SST",
                          "range_SST")
  
  model_data_filtered <- spp_data %>%
    # filter(!dataset %in% c("PRCRMP", "NODICE")) %>%  # Remove PRCRMP and NODICE datasets
    # filter(!grepl("_PR", PSU))  # Remove any PSU containing '_PR'
    # filter(!dataset %in% c("NODICE")) %>%
    filter(year(date) >= 2017) %>%
    filter(depth_bathy >= -50)
  
  # Extract all environmental values at species locations
  species_env_complex <- terra::extract(env_complex,
                                        cbind(model_data_filtered$x_utm,
                                              model_data_filtered$y_utm))
  
  # Add all variables to your filtered dataset
  model_data_filtered$depth_bathy <- species_env_complex$depth
  model_data_filtered$aspect <- species_env_complex$aspect
  model_data_filtered$slope <- species_env_complex$slope
  model_data_filtered$complexity <- species_env_complex$complexity
  model_data_filtered$TPI <- species_env_complex$TPI
  model_data_filtered$VRM <- species_env_complex$VRM
  model_data_filtered$planform_curv <- species_env_complex$planform_curv
  model_data_filtered$SAPA <- species_env_complex$SAPA
  model_data_filtered$max_Hsig <- species_env_complex$max_Hsig
  model_data_filtered$mean_dir <- species_env_complex$mean_dir
  model_data_filtered$mean_Hsig <- species_env_complex$mean_Hsig
  model_data_filtered$mean_SST <- species_env_complex$mean_SST
  model_data_filtered$range_SST <- species_env_complex$range_SST
  
  # Create complete cases dataset with all variables
  complexity_vars <- c("depth", "aspect", "slope", "complexity", "TPI", "VRM", "planform_curv",
                       "SAPA", "max_Hsig", "mean_dir", "mean_Hsig", "mean_SST",
                       "range_SST", "cover")

  model_data_complex <- model_data_filtered[complete.cases(model_data_filtered[, complexity_vars]), ]

  # Check how much data you have left
  cat("Observations with all complexity variables:", nrow(model_data_complex), "\n")
  
  # Fit expanded GAM models
  # gam_complex <- gam(cover ~ s(bathymetry) + s(aspect, bs = 'cc') + s(slope) + s(TPI) + s(roughness) + s(VRM) +
  #                      s(max_curv, k = 15) + s(mean_curv) + s(planform_curv) + s(profile_curv),
  #                    data = model_data_complex,
  #                    family = tw())
  # gam_complex <- gam(cover ~ s(depth_bathy) + TPI + s(planform_curv),
  #                    data = model_data_complex,
  #                    family = tw())
  gam_complex <- gam(cover ~ s(depth_bathy) + TPI + s(slope) + s(complexity) + s(planform_curv) +
                      s(range_SST) + s(mean_SST) + s(mean_dir, bs = 'cc') + s(max_Hsig) +
                       s(mean_Hsig),
                     data = model_data_complex,
                     family = tw())
  
  # Check the results
  summary(gam_complex)
  gam.check(gam_complex)
  draw(gam_complex)

  # Plot the relationships
  plot(gam_complex, pages = 2)  # Will create multiple pages
  
  # Compare model performance
  cat("Simple model AIC:", AIC(gam_tweedie), "\n")
  cat("Complex model AIC:", AIC(gam_complex), "\n")

  # Check correlations between variables first
  # complexity_matrix <- model_data_complex[, c("bathymetry", "aspect", "slope", "TPI", "roughness", "VRM",
  #                                             "max_curv", "mean_curv", "planform_curv", "profile_curv")]
  complexity_matrix <- model_data_complex[, c("depth_bathy", "TPI", "slope", "planform_curv",
                                              "range_SST", "mean_SST", "mean_dir", "max_Hsig",
                                              "mean_Hsig")]
  cor_matrix <- cor(complexity_matrix, use = "complete.obs")
  print(round(cor_matrix, 2))
  
  # Two-part model with complexity
  model_data_complex$present <- ifelse(model_data_complex$cover > 0, 1, 0)

  # gam_presence_complex <- gam(present ~ s(bathymetry) + s(slope) + s(TPI) + s(roughness) + s(VRM),
  #                             data = model_data_complex,
  #                             family = binomial())
  gam_presence_complex <- gam(present ~ s(depth_bathy) + TPI + s(slope) + s(planform_curv) +
                                s(range_SST) + s(mean_SST) + s(mean_dir, bs = 'cc') + s(max_Hsig) +
                                s(mean_Hsig),
                              data = model_data_complex,
                              family = binomial())
  
  # gam_abundance_complex <- gam(cover ~ s(bathymetry) + s(slope) + s(TPI) + s(roughness) + s(VRM),
  #                              data = model_data_complex[model_data_complex$cover > 0, ],
  #                              family = Gamma(link = "log"))
  gam_abundance_complex <- gam(cover ~ s(depth_bathy) + TPI + s(slope) + s(planform_curv) +
                                 s(range_SST) + s(mean_SST) + s(mean_dir, bs = 'cc') + s(max_Hsig) +
                                 s(mean_Hsig),
                               data = model_data_complex[model_data_complex$cover > 0, ],
                               family = Gamma(link = "log"))
  
  summary(gam_presence_complex)
  summary(gam_abundance_complex)

  # Create a reduced model with less correlated variables
  # Keep: bathymetry (depth), TPI (topographic position), VRM (rugosity)
  # These capture the main dimensions: depth, local position, and complexity

  gam_reduced <- gam(cover ~ s(depth_bathy) + TPI + s(planform_curv),
                     data = model_data_complex,
                     family = tw())

  summary(gam_reduced)

    # Plot the reduced model
  plot(gam_reduced, pages = 1)

  ################################## orbicellid GAMs ##################################
  
  # Convert tibble to data.frame first
  species_df <- as.data.frame(combined_benthic_data_averaged)
  
  spp_data = combined_benthic_data_averaged
  
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
  
  # Extract all the complexity variables at species locations
  env_complex <- c(bathy_final, aspect_terra, slope_terra, slopeofslope_terra, TPI_terra, VRM,
                   planformcurv_multiscale, max_hsig_raster, mean_dir_raster,
                   mean_hsig_raster, mean_sst_raster, range_sst_raster)
  
  names(env_complex) <- c("depth", "aspect", "slope", "complexity", "TPI", "VRM", "planform_curv",
                          "max_Hsig", "mean_dir", "mean_Hsig", "mean_SST",
                          "range_SST")
  
  # Create simple environmental stack
  env_simple <- c(bathy_final, slope_terra)
  names(env_simple) <- c("depth_bathy", "slope")
  
  # Extract environmental values at species locations
  species_env_values <- terra::extract(env_simple,
                                       cbind(spp_data$x_utm,
                                             spp_data$y_utm))
  
  # Add to dataframe
  spp_data$depth_bathy <- species_env_values$depth_bathy
  spp_data$slope <- species_env_values$slope
  
  orbicellids <- spp_data %>%
    filter(!dataset %in% c("PRCRMP", "NODICE")) %>%  # Remove PRCRMP and NODICE datasets
    filter(year(date) >= 2017) %>%
    filter(depth_bathy >= -50) %>%
    filter(!grepl("_PR", PSU)) %>%  # Remove any PSU containing '_PR' %>%
    filter(grepl("Orbicella", spp))
  
  # Extract all environmental values at species locations
  species_env_complex <- terra::extract(env_complex,
                                        cbind(orbicellids$x_utm,
                                              orbicellids$y_utm))

  # Add all variables to your filtered dataset
  orbicellids$depth_bathy <- species_env_complex$depth
  orbicellids$aspect <- species_env_complex$aspect
  orbicellids$slope <- species_env_complex$slope
  orbicellids$complexity <- species_env_complex$complexity
  orbicellids$TPI <- species_env_complex$TPI
  orbicellids$VRM <- species_env_complex$VRM
  orbicellids$planform_curv <- species_env_complex$planform_curv
  orbicellids$max_Hsig <- species_env_complex$max_Hsig
  orbicellids$mean_dir <- species_env_complex$mean_dir
  orbicellids$mean_Hsig <- species_env_complex$mean_Hsig
  orbicellids$mean_SST <- species_env_complex$mean_SST
  orbicellids$range_SST <- species_env_complex$range_SST
  
  # Create complete cases dataset with all variables
  complexity_vars <- c("depth_bathy", "aspect", "slope", "complexity", "TPI", "VRM", "planform_curv",
                       "max_Hsig", "mean_dir", "mean_Hsig", "mean_SST", "range_SST",
                       "cover")
  
  model_data_complex <- orbicellids[complete.cases(orbicellids[, complexity_vars]), ]

  # Check how much data you have left
  cat("Observations with all complexity variables:", nrow(model_data_complex), "\n")
  
  # Fit expanded GAM models
  # gam_complex <- gam(cover ~ s(bathymetry) + s(aspect, bs = 'cc') + s(slope) + s(TPI) + s(roughness) + s(VRM) +
  #                      s(max_curv, k = 15) + s(mean_curv) + s(planform_curv) + s(profile_curv),
  #                    data = model_data_complex,
  #                    family = tw())
  # gam_complex <- gam(cover ~ s(depth_bathy) + TPI + s(planform_curv),
  #                    data = model_data_complex,
  #                    family = tw())
  gam_complex <- gam(cover ~ s(depth_bathy) + TPI + s(slope) + s(complexity) +
                       s(planform_curv) + s(range_SST) + s(mean_SST) + s(mean_dir, bs = 'cc') +
                       s(max_Hsig) + s(mean_Hsig),
                     data = model_data_complex,
                     family = tw())
  
  # STOPPING POINT - 18 July 2025
  #   - continue to consider waves and SST, especially the landmasking for waves
  
  # Check the results
  summary(gam_complex)
  gam.check(gam_complex)
  
  # Plot the relationships
  plot(gam_complex, pages = 2)  # Will create multiple pages
  draw(gam_complex)
  
  # Compare model performance
  cat("Simple model AIC:", AIC(gam_tweedie), "\n")
  cat("Complex model AIC:", AIC(gam_complex), "\n")
  
  # Check correlations between variables first
  # complexity_matrix <- model_data_complex[, c("bathymetry", "aspect", "slope", "TPI", "roughness", "VRM",
  #                                             "max_curv", "mean_curv", "planform_curv", "profile_curv")]
  complexity_matrix <- model_data_complex[, c("depth_bathy", "TPI", "slope", "complexity",
                                              "planform_curv", "range_SST", "mean_SST",
                                              "mean_dir", "max_Hsig", "mean_Hsig")]
  cor_matrix <- cor(complexity_matrix, use = "complete.obs")
  print(round(cor_matrix, 2))
  
  # Two-part model with complexity
  model_data_complex$present <- ifelse(model_data_complex$cover > 0, 1, 0)
  
  gam_presence_complex <- gam(present ~ s(depth_bathy) + TPI + s(slope) + s(complexity) +
                                s(planform_curv) + s(range_SST) + s(mean_SST) + s(mean_dir, bs = 'cc') +
                                s(max_Hsig) + s(mean_Hsig),
                              data = model_data_complex,
                              family = binomial())
  
  gam_abundance_complex <- gam(cover ~ s(depth_bathy) + TPI + s(slope) + s(complexity) +
                                 s(planform_curv) + s(range_SST) + s(mean_SST) + s(mean_dir, bs = 'cc') +
                                 s(max_Hsig) + s(mean_Hsig),
                               data = model_data_complex[model_data_complex$cover > 0, ],
                               family = Gamma(link = "log"))
  
  summary(gam_presence_complex)
  summary(gam_abundance_complex)
  
  gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                     data = model_data_complex,
                     family = tw())
  
  summary(gam_reduced)
  
  # Plot the reduced model
  plot(gam_reduced, pages = 1)
  
  ################################## spp GAMs ##################################
  
  # Define species list (starting with your priority species)
  species_list <- c("Agaricia", "Montastraea", "Pseudodiploria", "Porites",
                    "Orbicella", "Colpophyllia", "Dendrogyra",
                    "Siderastrea")
  
  # # Get all unique species from the dataset and add remaining ones
  # all_species <- unique(spp_data$spp)
  # remaining_species <- all_species[!grepl(paste(c("Orbicella", species_list), collapse = "|"), all_species)]
  # species_list <- c(species_list, remaining_species)
  
  # Initialize results storage
  gam_results <- list()
  
  # STOPPING POINT - 11 August 2025
  #   - need to collapse all Isophyllia
  #   - still in the middle of figuring out exactly how I want to pursue model selection
  #   - should probably test for spatial and maybe temporal autocorrelation at some point
  #   - Investigate increasing sample size (years / PR), why montastraea & pseudodip. are struggling, hurdle vs. beta vs. gamma vs. Tweedie
  
  # Loop through each species
  for(sp in species_list) {
    cat("\n==================== Processing", sp, "====================\n")
    
    # Filter data for current species
    species_data <- spp_data %>%
      # filter(!dataset %in% c("PRCRMP", "NODICE")) %>%
      # filter(!grepl("_PR", PSU)) %>%
      # filter(!dataset %in% c("NODICE")) %>%
      filter(year(date) >= 2017) %>%
      filter(depth_bathy >= -50) %>%
      filter(grepl(sp, spp))
    
    # Check if we have enough data
    if(nrow(species_data) < 10) {
      cat("Insufficient data for", sp, "- only", nrow(species_data), "observations\n")
      next
    }
    
    # Extract environmental values
    species_env_complex <- terra::extract(env_complex, cbind(species_data$x_utm, species_data$y_utm))
    
    # Add environmental variables
    species_data$depth_bathy <- species_env_complex$depth
    species_data$aspect <- species_env_complex$aspect
    species_data$slope <- species_env_complex$slope
    species_data$complexity <- species_env_complex$complexity
    species_data$TPI <- species_env_complex$TPI
    species_data$VRM <- species_env_complex$VRM
    species_data$planform_curv <- species_env_complex$planform_curv
    species_data$max_Hsig <- species_env_complex$max_Hsig
    species_data$mean_dir <- species_env_complex$mean_dir
    species_data$mean_Hsig <- species_env_complex$mean_Hsig
    species_data$mean_SST <- species_env_complex$mean_SST
    species_data$range_SST <- species_env_complex$range_SST
    
    # Create complete cases dataset
    model_data_complete <- species_data[complete.cases(species_data[, complexity_vars]), ]
    
    cat("Complete observations:", nrow(model_data_complete), "\n")
    
    if(nrow(model_data_complete) < 10) {
      cat("Insufficient complete data for", sp, "\n")
      next
    }
    
    # Fit complex GAM
    tryCatch({
      gam_complex <- gam(cover ~ s(depth_bathy) + TPI + s(slope) + s(complexity) +
                           s(planform_curv) + s(range_SST) + s(mean_SST) + s(mean_dir, bs = 'cc') +
                           s(max_Hsig) + s(mean_Hsig),
                         data = model_data_complete, family = tw())
      
      cat("Complex GAM AIC:", AIC(gam_complex), "\n")
      
      # Fit reduced GAM
      gam_reduced <- gam(cover ~ s(depth_bathy) + s(slope) + TPI + s(planform_curv),
                         data = model_data_complete, family = tw())
      
      cat("Reduced GAM AIC:", AIC(gam_reduced), "\n")
      
      # Two-part models
      model_data_complete$present <- ifelse(model_data_complete$cover > 0, 1, 0)
      
      gam_presence <- gam(present ~ s(depth_bathy) + TPI + s(slope) + s(complexity) +
                            s(planform_curv) + s(range_SST) + s(mean_SST) + s(mean_dir, bs = 'cc') +
                            s(max_Hsig) + s(mean_Hsig),
                          data = model_data_complete, family = binomial())
      
      gam_abundance <- gam(cover ~ s(depth_bathy) + TPI + s(slope) + s(complexity) +
                             s(planform_curv) + s(range_SST) + s(mean_SST) + s(mean_dir, bs = 'cc') +
                             s(max_Hsig) + s(mean_Hsig),
                           data = model_data_complete[model_data_complete$cover > 0, ],
                           family = Gamma(link = "log"))
      
      # Store results
      gam_results[[sp]] <- list(
        data = model_data_complete,
        complex = gam_complex,
        reduced = gam_reduced,
        presence = gam_presence,
        abundance = gam_abundance,
        n_obs = nrow(model_data_complete),
        n_presence = sum(model_data_complete$cover > 0)
      )
      
      # Print summaries for priority species
      if(sp %in% c("Agaricia", "Montastraea", "Pseudodiploria", "Porites")) {
        cat("\n--- Summary for", sp, "---\n")
        print(summary(gam_complex))
        cat("\n--- Plotting", sp, "---\n")
        plot(gam_complex, pages = 1, main = paste(sp, "GAM"))
      }
      
    }, error = function(e) {
      cat("Error fitting GAM for", sp, ":", e$message, "\n")
    })
  }
  
  # Print overall results summary
  cat("\n==================== SUMMARY ====================\n")
  for(sp in names(gam_results)) {
    result <- gam_results[[sp]]
    cat(sp, "- N obs:", result$n_obs, "- N presence:", result$n_presence, 
        "- Complex AIC:", round(AIC(result$complex), 1), 
        "- Reduced AIC:", round(AIC(result$reduced), 1), "\n")
  }
  
  
  
  
  # AGARICIA
  summary(gam_results[["Agaricia"]]$complex)
  gam.check(gam_results[["Agaricia"]]$complex)
  plot(gam_results[["Agaricia"]]$complex, pages = 2)
  draw(gam_results[["Agaricia"]]$complex)
  
  # Complex model
  summary(gam_results[["Agaricia"]]$complex)
  plot(gam_results[["Agaricia"]]$complex, pages = 2)
  
  # Reduced model  
  summary(gam_results[["Agaricia"]]$reduced)
  plot(gam_results[["Agaricia"]]$reduced, pages = 1)
  
  # Presence/absence model
  summary(gam_results[["Agaricia"]]$presence)
  
  # Abundance model
  summary(gam_results[["Agaricia"]]$abundance)
  
  
  
  
  # MONTASTRAEA
  summary(gam_results[["Montastraea"]]$complex)
  gam.check(gam_results[["Montastraea"]]$complex)
  plot(gam_results[["Montastraea"]]$complex, pages = 2)
  draw(gam_results[["Montastraea"]]$complex)
  
  
  # Complex model
  summary(gam_results[["Montastraea"]]$complex)
  plot(gam_results[["Montastraea"]]$complex, pages = 2)
  
  # Reduced model  
  summary(gam_results[["Montastraea"]]$reduced)
  plot(gam_results[["Montastraea"]]$reduced, pages = 1)
  
  # Presence/absence model
  summary(gam_results[["Montastraea"]]$presence)
  
  # Abundance model
  summary(gam_results[["Montastraea"]]$abundance)
  
  
  
  
  # PSEUDODIPLORIA
  summary(gam_results[["Pseudodiploria"]]$complex)
  gam.check(gam_results[["Pseudodiploria"]]$complex)
  plot(gam_results[["Pseudodiploria"]]$complex, pages = 2)
  draw(gam_results[["Pseudodiploria"]]$complex)
  
  
  # Complex model
  summary(gam_results[["Pseudodiploria"]]$complex)
  plot(gam_results[["Pseudodiploria"]]$complex, pages = 2)
  
  # Reduced model  
  summary(gam_results[["Pseudodiploria"]]$reduced)
  plot(gam_results[["Pseudodiploria"]]$reduced, pages = 1)
  
  # Presence/absence model
  summary(gam_results[["Pseudodiploria"]]$presence)
  
  # Abundance model
  summary(gam_results[["Pseudodiploria"]]$abundance)
  
  
  
  # PORITES
  summary(gam_results[["Porites"]]$complex)
  gam.check(gam_results[["Porites"]]$complex)
  plot(gam_results[["Porites"]]$complex, pages = 2)
  draw(gam_results[["Porites"]]$complex)
  
  
  # Complex model
  summary(gam_results[["Porites"]]$complex)
  plot(gam_results[["Porites"]]$complex, pages = 2)
  
  # Reduced model  
  summary(gam_results[["Porites"]]$reduced)
  plot(gam_results[["Porites"]]$reduced, pages = 1)
  
  # Presence/absence model
  summary(gam_results[["Porites"]]$presence)
  
  # Abundance model
  summary(gam_results[["Porites"]]$abundance)
  
  
  # ORBICELLA
  summary(gam_results[["Orbicella"]]$complex)
  gam.check(gam_results[["Orbicella"]]$complex)
  plot(gam_results[["Orbicella"]]$complex, pages = 2)
  draw(gam_results[["Orbicella"]]$complex)
  
  
  # Complex model
  summary(gam_results[["Orbicella"]]$complex)
  plot(gam_results[["Orbicella"]]$complex, pages = 2)
  
  # Reduced model  
  summary(gam_results[["Orbicella"]]$reduced)
  plot(gam_results[["Orbicella"]]$reduced, pages = 1)
  
  # Presence/absence model
  summary(gam_results[["Orbicella"]]$presence)
  
  # Abundance model
  summary(gam_results[["Orbicella"]]$abundance)
  
  
  # COLPOPHYLLIA
  summary(gam_results[["Colpophyllia"]]$complex)
  gam.check(gam_results[["Colpophyllia"]]$complex)
  plot(gam_results[["Colpophyllia"]]$complex, pages = 2)
  draw(gam_results[["Colpophyllia"]]$complex)
  
  
  # Complex model
  summary(gam_results[["Colpophyllia"]]$complex)
  plot(gam_results[["Colpophyllia"]]$complex, pages = 2)
  
  # Reduced model  
  summary(gam_results[["Colpophyllia"]]$reduced)
  plot(gam_results[["Colpophyllia"]]$reduced, pages = 1)
  
  # Presence/absence model
  summary(gam_results[["Colpophyllia"]]$presence)
  
  # Abundance model
  summary(gam_results[["Colpophyllia"]]$abundance)
  
  
  
  # DENDROGYRA
  summary(gam_results[["Dendrogyra"]]$complex)
  gam.check(gam_results[["Dendrogyra"]]$complex)
  plot(gam_results[["Dendrogyra"]]$complex, pages = 2)
  draw(gam_results[["Dendrogyra"]]$complex)
  
  
  # Complex model
  summary(gam_results[["Dendrogyra"]]$complex)
  plot(gam_results[["Dendrogyra"]]$complex, pages = 2)
  
  # Reduced model  
  summary(gam_results[["Dendrogyra"]]$reduced)
  plot(gam_results[["Dendrogyra"]]$reduced, pages = 1)
  
  # Presence/absence model
  summary(gam_results[["Dendrogyra"]]$presence)
  
  # Abundance model
  summary(gam_results[["Dendrogyra"]]$abundance)
  
  
  # SIDERASTREA
  summary(gam_results[["Siderastrea"]]$complex)
  gam.check(gam_results[["Siderastrea"]]$complex)
  plot(gam_results[["Siderastrea"]]$complex, pages = 2)
  draw(gam_results[["Siderastrea"]]$complex)
  
  
  # Complex model
  summary(gam_results[["Siderastrea"]]$complex)
  plot(gam_results[["Siderastrea"]]$complex, pages = 2)
  
  # Reduced model  
  summary(gam_results[["Siderastrea"]]$reduced)
  plot(gam_results[["Siderastrea"]]$reduced, pages = 1)
  
  # Presence/absence model
  summary(gam_results[["Siderastrea"]]$presence)
  
  # Abundance model
  summary(gam_results[["Siderastrea"]]$abundance)
  
  # ################################## make predictions ##################################
  # 
  # # NOTE - this is sparse right now compared to all available rasters, for computational
  # #         reasons. can return to this
  # # First, create your stack properly
  # env_stack <- env_complex
  # 
  # # Check geometry - compareGeom works on individual rasters, not the stack
  # # Let's check that all layers match the first one (bathy_final)
  # terra::compareGeom(bathy_final, aspect_terra, stopOnError = FALSE)
  # terra::compareGeom(bathy_final, slope_terra, stopOnError = FALSE)
  # # etc... or loop through them:
  # 
  # for(i in 2:nlyr(env_stack)) {
  #   cat("Checking layer", i, ":", names(env_stack)[i], "\n")
  #   print(terra::compareGeom(env_stack[[1]], env_stack[[i]], stopOnError = FALSE))
  # }
  # 
  # # Alternative: just check basic properties
  # terra::res(env_stack)  # Should all be the same
  # terra::ext(env_stack)  # Should all be the same
  # 
  # # Create prediction grid
  # prediction_grid <- as.data.frame(env_stack, xy = TRUE)
  # 
  # # # Remove rows with any NA values
  # # prediction_grid <- prediction_grid[complete.cases(prediction_grid), ]
  # 
  # # Check your grid
  # dim(prediction_grid)
  # head(prediction_grid)
  # 
  # # Better names for your variables
  # names(prediction_grid) <- c("x", "y", "depth", "aspect", "slope", "complexity",
  #                             "TPI", "VRM", "planform_curv", "max_Hsig", "mean_dir",
  #                             "mean_Hsig", "mean_SST", "range_SST")
  # 
  # # Convert prediction grid UTM coordinates back to lat/lon for the spatial term
  # prediction_coords_utm <- vect(cbind(prediction_grid$x, prediction_grid$y),
  #                               crs = crs(bathy_final))
  # prediction_coords_latlon <- project(prediction_coords_utm, "EPSG:4326")
  # coords_df <- as.data.frame(geom(prediction_coords_latlon)[, c("x", "y")])
  # 
  # # Add lat/lon to prediction grid
  # prediction_grid$lon <- coords_df$x
  # prediction_grid$lat <- coords_df$y
  # 
  # # Remove any rows with NAs in the required variables
  # # pred_vars <- c("bathymetry", "TPI", "VRM", "lon", "lat")
  # pred_vars <- c("depth", "TPI", "slope", "complexity", "planform_curv",
  #                "range_SST", "mean_SST", "mean_dir", "max_Hsig", "mean_Hsig", "lon", "lat")
  # prediction_grid_clean <- prediction_grid[complete.cases(prediction_grid[, pred_vars]), ]
  # 
  # cat("Prediction grid size:", nrow(prediction_grid_clean), "cells\n")
  # 
  # 
  # # # Make predictions with the spatial model
  # # predictions <- predict(gam_spatial, prediction_grid_clean, type = "response")
  # # 
  # # # Add predictions to the grid
  # # prediction_grid_clean$predicted_cover <- predictions
  # # 
  # # # Check prediction range
  # # summary(prediction_grid_clean$predicted_cover)
  # 
  # 
  # # # Simple approach with time estimates
  # # start_time <- Sys.time()
  # # 
  # # # Test prediction on small subset to estimate time
  # # test_subset <- prediction_grid_clean[1:1000, ]
  # # test_start <- Sys.time()
  # # test_pred <- predict(gam_complex, test_subset, type = "response")
  # # test_time <- as.numeric(difftime(Sys.time(), test_start, units = "secs"))
  # # 
  # # # Estimate total time
  # # estimated_total <- (test_time / 1000) * nrow(prediction_grid_clean)
  # # cat("Estimated total time:", round(estimated_total/60, 1), "minutes\n")
  # # 
  # # # Now do full prediction
  # # predictions <- predict(gam_spatial, prediction_grid_clean, type = "response")
  # # cat("Actual time:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes\n")
  # 
  # 
  # # Resample your environmental stack to 200m resolution
  # env_simple_200m <- aggregate(env_stack, fact = 4, fun = "mean")  # 4x aggregation = 200m
  # 
  # # Check the new resolution
  # res(env_simple_200m)  # Should show 200, 200
  # 
  # # Create prediction grid from 200m rasters
  # prediction_grid_200m <- as.data.frame(env_simple_200m, xy = TRUE)
  # names(prediction_grid_200m) <- c("x", "y", "bathymetry", "slope")
  # 
  # # Add TPI and VRM at 200m
  # env_complex_200m <- aggregate(env_complex, fact = 4, fun = "mean")
  # complex_grid_200m <- as.data.frame(env_complex_200m, xy = TRUE)
  # names(complex_grid_200m) <- c("x", "y", "bathymetry", "slope", "TPI", "roughness", "VRM",
  #                               "max_curv", "mean_curv", "planform_curv", "profile_curv")
  # 
  # # Keep just the variables you need
  # prediction_grid_200m <- complex_grid_200m[, c("x", "y", "bathymetry", "TPI", "VRM")]
  # 
  # # Remove NAs
  # prediction_grid_200m <- prediction_grid_200m[complete.cases(prediction_grid_200m), ]
  # 
  # cat("200m grid size:", nrow(prediction_grid_200m), "cells\n")
  # cat("50m grid size was:", nrow(prediction_grid_clean), "cells\n")
  # cat("Reduction factor:", round(nrow(prediction_grid_clean)/nrow(prediction_grid_200m), 1), "\n")
  # 
  # # Convert UTM to lat/lon for the spatial term
  # prediction_coords_utm_200m <- vect(cbind(prediction_grid_200m$x, prediction_grid_200m$y),
  #                                    crs = crs(bathy_final))
  # prediction_coords_latlon_200m <- project(prediction_coords_utm_200m, "EPSG:4326")
  # coords_df_200m <- as.data.frame(geom(prediction_coords_latlon_200m)[, c("x", "y")])
  # 
  # prediction_grid_200m$lon <- coords_df_200m$x
  # prediction_grid_200m$lat <- coords_df_200m$y
  # 
  # cat("Final 200m grid size:", nrow(prediction_grid_200m), "cells\n")
  # 
  # # Step 4: Estimate time before full prediction
  # start_time <- Sys.time()
  # 
  # # Test prediction on small subset to estimate time
  # test_subset_200m <- prediction_grid_200m[1:1000, ]
  # test_start <- Sys.time()
  # test_pred_200m <- predict(gam_spatial, test_subset_200m, type = "response")
  # test_time <- as.numeric(difftime(Sys.time(), test_start, units = "secs"))
  # 
  # # Estimate total time
  # estimated_total_200m <- (test_time / 1000) * nrow(prediction_grid_200m)
  # cat("200m grid size:", nrow(prediction_grid_200m), "cells\n")
  # cat("Test time for 1000 cells:", round(test_time, 2), "seconds\n")
  # cat("Estimated total time:", round(estimated_total_200m/60, 1), "minutes\n")
  # 
  # # Compare to original estimate
  # original_cells <- nrow(prediction_grid_clean)
  # reduction_factor <- original_cells / nrow(prediction_grid_200m)
  # cat("Reduction factor:", round(reduction_factor, 1), "x smaller\n")
  # cat("Original estimate was 61 minutes\n")
  # cat("New estimate should be ~", round(61/reduction_factor, 1), "minutes\n")
  # 
  # # Ask user if they want to proceed
  # cat("\nProceed with full prediction? (y/n)\n")
  # user_input <- readline()
  # 
  # if(tolower(user_input) == "y") {
  #   cat("Starting full prediction...\n")
  # 
  #   # Full prediction with timing
  #   full_start <- Sys.time()
  #   predictions_200m <- predict(gam_spatial, prediction_grid_200m, type = "response")
  #   prediction_grid_200m$predicted_cover <- predictions_200m
  #   full_end <- Sys.time()
  # 
  #   cat("Actual prediction time:", round(difftime(full_end, full_start, units = "mins"), 1), "minutes\n")
  #   cat("Prediction complete!\n")
  # 
  # } else {
  #   cat("Prediction cancelled. Consider further reducing grid size or using sampling.\n")
  # }
  # 
  # # Convert to raster
  # pred_raster_200m <- rast(prediction_grid_200m[, c("x", "y", "predicted_cover")],
  #                          crs = crs(bathy_final))
  # 
  # # Plot
  # plot(pred_raster_200m,
  #      main = "Predicted Coral Cover (%) - 200m resolution",
  #      col = viridis(100))
  # 
  # # Add survey points
  # points(model_data_complex$x_utm, model_data_complex$y_utm,
  #        pch = 21,
  #        bg = heat.colors(10)[cut(model_data_complex$cover, breaks = 10)],
  #        cex = 0.8)
  # 
  # 
  # ################################## maps ##################################
  # 
  # 
  # # library(viridis)
  # library(sf)
  # 
  # # Convert predictions back to raster
  # pred_raster <- rast(prediction_grid_clean[, c("x", "y", "predicted_cover")], 
  #                     crs = crs(bathy_final))
  # 
  # # Create a nice map
  # # Option 1: Using terra plot
  # plot(pred_raster, 
  #      main = "Predicted Coral Cover (%)",
  #      col = viridis(100),
  #      range = c(0, max(prediction_grid_clean$predicted_cover, na.rm = TRUE)))
  # 
  # # Add your actual survey points
  # points(model_data_complex$x_utm, model_data_complex$y_utm, 
  #        pch = 21, 
  #        bg = heat.colors(10)[cut(model_data_complex$cover, breaks = 10)],
  #        cex = 0.8)
  # 
  # # Add legend for points
  # legend("topright", 
  #        title = "Observed Cover",
  #        legend = c("0-10%", "10-20%", "20-30%", "30%+"),
  #        pch = 21,
  #        pt.bg = heat.colors(4),
  #        cex = 0.8)
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  
  
  ################################## 80/20 w/ Tweedie ##################################
  
  gam_train_test_split <- function(gam_model, gam_results_entry = NULL, train_percent = 80, 
                                   seed = 123, return_details = FALSE) {
    
    # Load required libraries
    library(mgcv)
    library(pROC)
    library(dplyr)
    library(ggplot2)
    library(viridis)
    library(gridExtra)
    
    set.seed(seed)
    
    # Extract model information and use model's own data
    model_formula <- formula(gam_model)
    model_family <- family(gam_model)
    model_data <- gam_model$model
    
    # If full gam_results entry provided, use that data for coordinates
    if(!is.null(gam_results_entry) && !is.null(gam_results_entry$data)) {
      full_data <- gam_results_entry$data
      # Match the model data with full data to get coordinates
      # Assuming the model data is a subset of the full data
      model_data_with_coords <- merge(model_data, full_data, by = intersect(names(model_data), names(full_data)))
    } else {
      model_data_with_coords <- model_data
    }
    
    # Get response variable name
    response_var <- all.vars(model_formula)[1]
    
    # Remove rows with missing values for the variables in the formula
    formula_vars <- all.vars(model_formula)
    complete_data <- model_data_with_coords[complete.cases(model_data_with_coords[, formula_vars]), ]
    
    cat("Using", nrow(complete_data), "complete observations from model's data\n")
    cat("Model family:", model_family$family, "\n")
    cat("Train/test split:", train_percent, "/", 100 - train_percent, "\n")
    cat("Response variable:", response_var, "\n")
    
    if(nrow(complete_data) < 20) {
      warning("Insufficient data for validation (< 20 observations)")
      return(NULL)
    }
    
    # Determine if this is predicting binary outcomes (0/1 only) vs continuous
    unique_values <- unique(complete_data[[response_var]])
    unique_values <- unique_values[!is.na(unique_values)]
    is_binary <- length(unique_values) == 2 && all(unique_values %in% c(0, 1))
    
    # Calculate appropriate metrics
    calculate_metrics <- function(actual, predicted) {
      if(is_binary) {
        # For binary outcomes - calculate AUC and classification metrics
        if(length(unique(actual)) > 1) {
          auc_val <- as.numeric(auc(actual, predicted))
          pred_class <- ifelse(predicted > 0.5, 1, 0)
          accuracy <- mean(actual == pred_class)
          sensitivity <- if(sum(actual == 1) > 0) sum(actual == 1 & pred_class == 1) / sum(actual == 1) else NA
          specificity <- if(sum(actual == 0) > 0) sum(actual == 0 & pred_class == 0) / sum(actual == 0) else NA
          
          return(list(
            auc = auc_val,
            accuracy = accuracy,
            sensitivity = sensitivity,
            specificity = specificity
          ))
        } else {
          return(list(auc = NA, accuracy = NA, sensitivity = NA, specificity = NA))
        }
      } else {
        # For continuous outcomes - calculate RMSE, MAE, R², and Pearson correlation
        rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
        mae <- mean(abs(actual - predicted), na.rm = TRUE)
        
        # Calculate R²
        if(var(actual, na.rm = TRUE) > 0) {
          r_squared <- cor(actual, predicted, use = "complete.obs")^2
        } else {
          r_squared <- NA
        }
        
        # Calculate Pearson correlation coefficient and test significance
        cor_test <- cor.test(actual, predicted, method = "pearson")
        pearson_r <- cor_test$estimate
        pearson_p <- cor_test$p.value
        
        # Determine significance level
        if(is.na(pearson_p)) {
          significance <- "Cannot determine"
        } else if(pearson_p < 0.001) {
          significance <- "***"
        } else if(pearson_p < 0.01) {
          significance <- "**"
        } else if(pearson_p < 0.05) {
          significance <- "*"
        } else {
          significance <- "ns"
        }
        
        return(list(
          rmse = rmse,
          mae = mae,
          r_squared = r_squared,
          pearson_r = pearson_r,
          pearson_p = pearson_p,
          significance = significance
        ))
      }
    }
    
    # Train-test split
    n <- nrow(complete_data)
    train_size <- floor((train_percent / 100) * n)
    train_idx <- sample(1:n, size = train_size)
    
    train_data <- complete_data[train_idx, ]
    test_data <- complete_data[-train_idx, ]
    
    cat("Training observations:", nrow(train_data), "\n")
    cat("Testing observations:", nrow(test_data), "\n")
    
    # Fit model on training data
    tryCatch({
      fit <- gam(model_formula, data = train_data, family = model_family)
      
      # Make predictions on test data
      predictions <- predict(fit, newdata = test_data, type = "response")
      actual <- test_data[[response_var]]
      
      # Calculate metrics
      metrics <- calculate_metrics(actual, predictions)
      
      # Print key results
      cat("\n=== RESULTS ===\n")
      if(is_binary) {
        cat("AUC:", round(metrics$auc, 3), "\n")
        cat("Accuracy:", round(metrics$accuracy, 3), "\n")
      } else {
        cat("R²:", round(metrics$r_squared, 3), "\n")
        cat("Pearson r:", round(metrics$pearson_r, 3), "\n")
        cat("P-value:", format(metrics$pearson_p, scientific = TRUE, digits = 3), "\n")
        cat("Significance:", metrics$significance, "\n")
        cat("RMSE:", round(metrics$rmse, 3), "\n")
        
        # Provide interpretation
        if(metrics$significance == "***") {
          cat("Interpretation: Correlation is highly significant (P < 0.001)\n")
        } else if(metrics$significance == "**") {
          cat("Interpretation: Correlation is significant (P < 0.01)\n")
        } else if(metrics$significance == "*") {
          cat("Interpretation: Correlation is significant (P < 0.05)\n")
        } else if(metrics$significance == "ns") {
          cat("Interpretation: Correlation is not significant (P ≥ 0.05)\n")
        }
      }
      
      # Create results dataframe with spatial coordinates for mapping
      test_results <- test_data
      test_results$predicted <- predictions
      test_results$residuals <- actual - predictions
      
      result <- list(
        n_total = nrow(complete_data),
        n_train = nrow(train_data),
        n_test = nrow(test_data),
        train_percent = train_percent,
        family = model_family$family,
        response_variable = response_var,
        metrics = metrics,
        test_results = test_results,  # Always include this for mapping
        is_binary = is_binary
      )
      
      if(return_details) {
        result$predictions <- predictions
        result$actual <- actual
        result$train_data <- train_data
        result$test_data <- test_data
        result$fitted_model <- fit
        result$formula <- model_formula
      }
      
      return(result)
      
    }, error = function(e) {
      cat("Error in model fitting:", e$message, "\n")
      return(NULL)
    })
  }
  
  ################################## 80/20 plot ##################################
  
  plot_validation_maps <- function(validation_result, coord_system = "latlon") {
    
    if(is.null(validation_result) || is.null(validation_result$test_results)) {
      stop("No test results found for mapping")
    }
    
    plot_data <- validation_result$test_results
    
    # Print available columns for debugging
    cat("Available columns in test results:\n")
    print(colnames(plot_data))
    
    # Determine coordinate columns with fallbacks
    if(coord_system == "latlon") {
      x_candidates <- c("lon", "longitude", "Longitude", "LON")
      y_candidates <- c("lat", "latitude", "Latitude", "LAT")
      x_lab <- "Longitude"
      y_lab <- "Latitude"
    } else {
      x_candidates <- c("x_utm", "X_utm", "utm_x", "UTM_X")
      y_candidates <- c("y_utm", "Y_utm", "utm_y", "UTM_Y")
      x_lab <- "UTM X"
      y_lab <- "UTM Y"
    }
    
    # Find the actual column names
    x_col <- x_candidates[x_candidates %in% colnames(plot_data)]
    y_col <- y_candidates[y_candidates %in% colnames(plot_data)]
    
    if(length(x_col) == 0 || length(y_col) == 0) {
      cat("Could not find coordinate columns. Trying the other coordinate system...\n")
      
      # Try the other coordinate system
      if(coord_system == "latlon") {
        x_candidates <- c("x_utm", "X_utm", "utm_x", "UTM_X")
        y_candidates <- c("y_utm", "Y_utm", "utm_y", "UTM_Y")
        x_lab <- "UTM X"
        y_lab <- "UTM Y"
      } else {
        x_candidates <- c("lon", "longitude", "Longitude", "LON")
        y_candidates <- c("lat", "latitude", "Latitude", "LAT")
        x_lab <- "Longitude"
        y_lab <- "Latitude"
      }
      
      x_col <- x_candidates[x_candidates %in% colnames(plot_data)]
      y_col <- y_candidates[y_candidates %in% colnames(plot_data)]
    }
    
    if(length(x_col) == 0 || length(y_col) == 0) {
      stop("Could not find any coordinate columns. Available columns: ", paste(colnames(plot_data), collapse = ", "))
    }
    
    # Use first match
    x_col <- x_col[1]
    y_col <- y_col[1]
    
    cat("Using coordinates:", x_col, "and", y_col, "\n")
    
    # Common theme for maps
    map_theme <- theme_minimal() +
      theme(
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)
      )
    
    # Get response variable name
    response_var <- validation_result$response_variable
    
    # Calculate shared color scale limits
    actual_values <- plot_data[[response_var]]
    predicted_values <- plot_data$predicted
    shared_min <- min(c(actual_values, predicted_values), na.rm = TRUE)
    shared_max <- max(c(actual_values, predicted_values), na.rm = TRUE)
    
    cat("Shared color scale range:", round(shared_min, 2), "to", round(shared_max, 2), "\n")
    
    # Plot 1: Actual cover values
    p1 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = response_var)) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_viridis_c(name = "Actual", option = "plasma", limits = c(shared_min, shared_max)) +
      labs(title = paste("Actual", response_var, "(Test Set)"), 
           x = x_lab, y = y_lab) +
      map_theme
    
    # Plot 2: Predicted cover values
    p2 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = "predicted")) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_viridis_c(name = "Predicted", option = "plasma", limits = c(shared_min, shared_max)) +
      labs(title = paste("Predicted", response_var), 
           x = x_lab, y = y_lab) +
      map_theme
    
    # Plot 3: Residuals (actual - predicted)
    p3 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = "residuals")) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_gradient2(name = "Residual", 
                            low = "red", mid = "white", high = "blue",
                            midpoint = 0) +
      labs(title = "Residuals (Actual - Predicted)", 
           x = x_lab, y = y_lab) +
      map_theme
    
    # Plot 4: Predicted vs Actual scatter plot
    p4 <- ggplot(plot_data, aes_string(x = "predicted", y = response_var)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      labs(title = "Predicted vs Actual", 
           x = paste("Predicted", response_var), 
           y = paste("Actual", response_var)) +
      theme_minimal()
    
    # Arrange plots
    grid.arrange(p1, p2, p3, p4, ncol = 2,
                 top = paste("Validation Results -", 
                             validation_result$family, "model"))
    
    # Return the individual plots as well
    invisible(list(actual = p1, predicted = p2, residuals = p3, scatter = p4))
  }
  
  ################################## 80/20 usage ##################################
  
  # Example 1: Validation with coordinates for mapping
  result <- gam_train_test_split(gam_results[["Orbicella"]]$complex,
                                 gam_results_entry = gam_results[["Orbicella"]],
                                 train_percent = 80)
  plot_validation_maps(result)
  
  # Example 2: Simple validation without coordinates (no mapping possible)
  # result <- gam_train_test_split(gam_results[["Agaricia"]]$complex, train_percent = 80)
  
  # Example 3: Validate all species models with mapping
  # for(sp in names(gam_results)) {
  #   cat("\n=== VALIDATING", sp, "===\n")
  #   result <- gam_train_test_split(gam_results[[sp]]$complex, 
  #                                 gam_results_entry = gam_results[[sp]], 
  #                                 train_percent = 80)
  #   if(!is.null(result)) {
  #     plot_validation_maps(result)
  #   }
  # }
  
  # Example 4: Generate summary table for multiple species
  # results_summary <- data.frame(
  #   Species = character(),
  #   Model_Type = character(),
  #   Pearson_r = numeric(),
  #   P_value = numeric(),
  #   Significance = character(),
  #   R_squared = numeric(),
  #   RMSE = numeric(),
  #   stringsAsFactors = FALSE
  # )
  # 
  # for(sp in names(gam_results)) {
  #   for(model_type in c("complex", "reduced")) {
  #     result <- gam_train_test_split(gam_results[[sp]][[model_type]], 
  #                                   gam_results_entry = gam_results[[sp]], 
  #                                   train_percent = 80)
  #     if(!is.null(result) && !result$is_binary) {
  #       results_summary <- rbind(results_summary, data.frame(
  #         Species = sp,
  #         Model_Type = model_type,
  #         Pearson_r = round(result$metrics$pearson_r, 3),
  #         P_value = result$metrics$pearson_p,
  #         Significance = result$metrics$significance,
  #         R_squared = round(result$metrics$r_squared, 3),
  #         RMSE = round(result$metrics$rmse, 3)
  #       ))
  #     }
  #   }
  # }
  # print(results_summary)
  
  
  ################################## GAM Train-Test Split Function ##################################
  
  gam_train_test_split <- function(gam_model, gam_results_entry = NULL, train_percent = 80, 
                                   seed = 123, return_details = FALSE, use_hurdle = FALSE) {
    
    # Load required libraries
    library(mgcv)
    library(pROC)
    library(dplyr)
    library(ggplot2)
    library(viridis)
    library(gridExtra)
    
    set.seed(seed)
    
    if(use_hurdle && !is.null(gam_results_entry)) {
      # Use two-part hurdle model approach
      return(gam_hurdle_validation(gam_results_entry, train_percent, seed, return_details))
    }
    
    # Original single-model approach
    # Extract model information and use model's own data
    model_formula <- formula(gam_model)
    model_family <- family(gam_model)
    model_data <- gam_model$model
    
    # If full gam_results entry provided, use that data for coordinates
    if(!is.null(gam_results_entry) && !is.null(gam_results_entry$data)) {
      full_data <- gam_results_entry$data
      # Match the model data with full data to get coordinates
      # Assuming the model data is a subset of the full data
      model_data_with_coords <- merge(model_data, full_data, by = intersect(names(model_data), names(full_data)))
    } else {
      model_data_with_coords <- model_data
    }
    
    # Get response variable name
    response_var <- all.vars(model_formula)[1]
    
    # Remove rows with missing values for the variables in the formula
    formula_vars <- all.vars(model_formula)
    complete_data <- model_data_with_coords[complete.cases(model_data_with_coords[, formula_vars]), ]
    
    cat("Using", nrow(complete_data), "complete observations from model's data\n")
    cat("Model family:", model_family$family, "\n")
    cat("Train/test split:", train_percent, "/", 100 - train_percent, "\n")
    cat("Response variable:", response_var, "\n")
    
    if(nrow(complete_data) < 20) {
      warning("Insufficient data for validation (< 20 observations)")
      return(NULL)
    }
    
    # Determine if this is predicting binary outcomes (0/1 only) vs continuous
    unique_values <- unique(complete_data[[response_var]])
    unique_values <- unique_values[!is.na(unique_values)]
    is_binary <- length(unique_values) == 2 && all(unique_values %in% c(0, 1))
    
    # Calculate appropriate metrics
    calculate_metrics <- function(actual, predicted) {
      if(is_binary) {
        # For binary outcomes - calculate AUC and classification metrics
        if(length(unique(actual)) > 1) {
          auc_val <- as.numeric(auc(actual, predicted))
          pred_class <- ifelse(predicted > 0.5, 1, 0)
          accuracy <- mean(actual == pred_class)
          sensitivity <- if(sum(actual == 1) > 0) sum(actual == 1 & pred_class == 1) / sum(actual == 1) else NA
          specificity <- if(sum(actual == 0) > 0) sum(actual == 0 & pred_class == 0) / sum(actual == 0) else NA
          
          return(list(
            auc = auc_val,
            accuracy = accuracy,
            sensitivity = sensitivity,
            specificity = specificity
          ))
        } else {
          return(list(auc = NA, accuracy = NA, sensitivity = NA, specificity = NA))
        }
      } else {
        # For continuous outcomes - calculate RMSE, MAE, R², and Pearson correlation
        rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
        mae <- mean(abs(actual - predicted), na.rm = TRUE)
        
        # Calculate R²
        if(var(actual, na.rm = TRUE) > 0) {
          r_squared <- cor(actual, predicted, use = "complete.obs")^2
        } else {
          r_squared <- NA
        }
        
        # Calculate Pearson correlation coefficient and test significance
        cor_test <- cor.test(actual, predicted, method = "pearson")
        pearson_r <- cor_test$estimate
        pearson_p <- cor_test$p.value
        
        # Determine significance level
        if(is.na(pearson_p)) {
          significance <- "Cannot determine"
        } else if(pearson_p < 0.001) {
          significance <- "***"
        } else if(pearson_p < 0.01) {
          significance <- "**"
        } else if(pearson_p < 0.05) {
          significance <- "*"
        } else {
          significance <- "ns"
        }
        
        return(list(
          rmse = rmse,
          mae = mae,
          r_squared = r_squared,
          pearson_r = pearson_r,
          pearson_p = pearson_p,
          significance = significance
        ))
      }
    }
    
    # Train-test split
    n <- nrow(complete_data)
    train_size <- floor((train_percent / 100) * n)
    train_idx <- sample(1:n, size = train_size)
    
    train_data <- complete_data[train_idx, ]
    test_data <- complete_data[-train_idx, ]
    
    cat("Training observations:", nrow(train_data), "\n")
    cat("Testing observations:", nrow(test_data), "\n")
    
    # Fit model on training data
    tryCatch({
      fit <- gam(model_formula, data = train_data, family = model_family)
      
      # Make predictions on test data
      predictions <- predict(fit, newdata = test_data, type = "response")
      actual <- test_data[[response_var]]
      
      # Calculate metrics
      metrics <- calculate_metrics(actual, predictions)
      
      # Print key results
      cat("\n=== RESULTS ===\n")
      if(is_binary) {
        cat("AUC:", round(metrics$auc, 3), "\n")
        cat("Accuracy:", round(metrics$accuracy, 3), "\n")
      } else {
        cat("R²:", round(metrics$r_squared, 3), "\n")
        cat("Pearson r:", round(metrics$pearson_r, 3), "\n")
        cat("P-value:", format(metrics$pearson_p, scientific = TRUE, digits = 3), "\n")
        cat("Significance:", metrics$significance, "\n")
        cat("RMSE:", round(metrics$rmse, 3), "\n")
        
        # Provide interpretation
        if(metrics$significance == "***") {
          cat("Interpretation: Correlation is highly significant (P < 0.001)\n")
        } else if(metrics$significance == "**") {
          cat("Interpretation: Correlation is significant (P < 0.01)\n")
        } else if(metrics$significance == "*") {
          cat("Interpretation: Correlation is significant (P < 0.05)\n")
        } else if(metrics$significance == "ns") {
          cat("Interpretation: Correlation is not significant (P ≥ 0.05)\n")
        }
      }
      
      # Create results dataframe with spatial coordinates for mapping
      test_results <- test_data
      test_results$predicted <- predictions
      test_results$residuals <- actual - predictions
      
      result <- list(
        n_total = nrow(complete_data),
        n_train = nrow(train_data),
        n_test = nrow(test_data),
        train_percent = train_percent,
        family = model_family$family,
        response_variable = response_var,
        metrics = metrics,
        test_results = test_results,  # Always include this for mapping
        is_binary = is_binary
      )
      
      if(return_details) {
        result$predictions <- predictions
        result$actual <- actual
        result$train_data <- train_data
        result$test_data <- test_data
        result$fitted_model <- fit
        result$formula <- model_formula
      }
      
      return(result)
      
    }, error = function(e) {
      cat("Error in model fitting:", e$message, "\n")
      return(NULL)
    })
  }
  
  ################################## Two-Part Hurdle Model Validation ##################################
  
  gam_hurdle_validation <- function(gam_results_entry, train_percent = 80, seed = 123, return_details = FALSE) {
    
    set.seed(seed)
    
    # Check if both presence and abundance models exist
    if(is.null(gam_results_entry$presence) || is.null(gam_results_entry$abundance)) {
      stop("Both presence and abundance models required for hurdle validation")
    }
    
    # Get the full dataset
    full_data <- gam_results_entry$data
    
    # Get formulas from both models
    presence_formula <- formula(gam_results_entry$presence)
    abundance_formula <- formula(gam_results_entry$abundance)
    
    # Get all variables needed
    presence_vars <- all.vars(presence_formula)
    abundance_vars <- all.vars(abundance_formula)
    all_vars <- unique(c(presence_vars, abundance_vars))
    
    # Filter to complete cases for all variables
    complete_data <- full_data[complete.cases(full_data[, all_vars]), ]
    
    cat("=== HURDLE MODEL VALIDATION ===\n")
    cat("Using", nrow(complete_data), "complete observations\n")
    cat("Train/test split:", train_percent, "/", 100 - train_percent, "\n")
    
    if(nrow(complete_data) < 20) {
      warning("Insufficient data for validation (< 20 observations)")
      return(NULL)
    }
    
    # Train-test split
    n <- nrow(complete_data)
    train_size <- floor((train_percent / 100) * n)
    train_idx <- sample(1:n, size = train_size)
    
    train_data <- complete_data[train_idx, ]
    test_data <- complete_data[-train_idx, ]
    
    cat("Training observations:", nrow(train_data), "\n")
    cat("Testing observations:", nrow(test_data), "\n")
    
    tryCatch({
      # Fit presence model on training data
      presence_fit <- gam(presence_formula, data = train_data, family = binomial())
      
      # Fit abundance model on training data (only where cover > 0)
      abundance_train_data <- train_data[train_data$cover > 0, ]
      abundance_fit <- gam(abundance_formula, data = abundance_train_data, family = Gamma(link = "log"))
      
      # Make predictions on test data
      # Step 1: Predict presence probability
      presence_prob <- predict(presence_fit, newdata = test_data, type = "response")
      
      # Step 2: Predict abundance (conditional on presence)
      abundance_pred <- predict(abundance_fit, newdata = test_data, type = "response")
      
      # Step 3: Combine using hurdle approach with threshold
      # Method 1: Expected value = P(presence) * E(abundance | presence)
      hurdle_predictions_raw <- presence_prob * abundance_pred
      
      # Apply threshold for more realistic zero predictions
      presence_threshold <- 0.5  # Can be adjusted
      probability_threshold <- 0.1  # Very low probability predictions become zero
      
      # Method A: Use presence classification (0.5 threshold)
      presence_binary <- ifelse(presence_prob > presence_threshold, 1, 0)
      classification_predictions <- presence_binary * abundance_pred
      
      # Method B: Use probability threshold (very low probabilities become zero)
      hurdle_predictions_threshold <- ifelse(presence_prob < probability_threshold, 0, hurdle_predictions_raw)
      
      # Use the classification approach as default (more zeros)
      hurdle_predictions <- classification_predictions
      
      actual <- test_data$cover
      
      # Calculate metrics for both approaches
      calculate_hurdle_metrics <- function(actual, predicted) {
        rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
        mae <- mean(abs(actual - predicted), na.rm = TRUE)
        
        if(var(actual, na.rm = TRUE) > 0) {
          r_squared <- cor(actual, predicted, use = "complete.obs")^2
        } else {
          r_squared <- NA
        }
        
        cor_test <- cor.test(actual, predicted, method = "pearson")
        pearson_r <- cor_test$estimate
        pearson_p <- cor_test$p.value
        
        if(is.na(pearson_p)) {
          significance <- "Cannot determine"
        } else if(pearson_p < 0.001) {
          significance <- "***"
        } else if(pearson_p < 0.01) {
          significance <- "**"
        } else if(pearson_p < 0.05) {
          significance <- "*"
        } else {
          significance <- "ns"
        }
        
        # Calculate how many zeros are correctly predicted
        actual_zeros <- sum(actual == 0)
        predicted_zeros <- sum(predicted == 0)
        correct_zeros <- sum(actual == 0 & predicted == 0)
        
        return(list(
          rmse = rmse,
          mae = mae,
          r_squared = r_squared,
          pearson_r = pearson_r,
          pearson_p = pearson_p,
          significance = significance,
          actual_zeros = actual_zeros,
          predicted_zeros = predicted_zeros,
          correct_zeros = correct_zeros
        ))
      }
      
      # Use the expected value approach as primary
      metrics <- calculate_hurdle_metrics(actual, hurdle_predictions)
      
      # Print results
      cat("\n=== HURDLE MODEL RESULTS ===\n")
      cat("R²:", round(metrics$r_squared, 3), "\n")
      cat("Pearson r:", round(metrics$pearson_r, 3), "\n")
      cat("P-value:", format(metrics$pearson_p, scientific = TRUE, digits = 3), "\n")
      cat("Significance:", metrics$significance, "\n")
      cat("RMSE:", round(metrics$rmse, 3), "\n")
      cat("Actual zeros:", metrics$actual_zeros, "/ Predicted zeros:", metrics$predicted_zeros, "\n")
      cat("Correctly predicted zeros:", metrics$correct_zeros, "\n")
      cat("Presence threshold used:", presence_threshold, "\n")
      
      # Provide interpretation
      if(metrics$significance == "***") {
        cat("Interpretation: Correlation is highly significant (P < 0.001)\n")
      } else if(metrics$significance == "**") {
        cat("Interpretation: Correlation is significant (P < 0.01)\n")
      } else if(metrics$significance == "*") {
        cat("Interpretation: Correlation is significant (P < 0.05)\n")
      } else if(metrics$significance == "ns") {
        cat("Interpretation: Correlation is not significant (P ≥ 0.05)\n")
      }
      
      # Create results dataframe
      test_results <- test_data
      test_results$predicted <- hurdle_predictions
      test_results$presence_prob <- presence_prob
      test_results$abundance_pred <- abundance_pred
      test_results$residuals <- actual - hurdle_predictions
      
      result <- list(
        n_total = nrow(complete_data),
        n_train = nrow(train_data),
        n_test = nrow(test_data),
        train_percent = train_percent,
        family = "Hurdle (Binomial + Gamma)",
        response_variable = "cover",
        metrics = metrics,
        test_results = test_results,
        is_binary = FALSE,
        model_type = "hurdle"
      )
      
      if(return_details) {
        result$presence_fit <- presence_fit
        result$abundance_fit <- abundance_fit
        result$presence_predictions <- presence_prob
        result$abundance_predictions <- abundance_pred
        result$hurdle_predictions <- hurdle_predictions
        result$hurdle_predictions_raw <- hurdle_predictions_raw
        result$hurdle_predictions_threshold <- hurdle_predictions_threshold
        result$presence_threshold <- presence_threshold
        result$probability_threshold <- probability_threshold
      }
      
      return(result)
      
    }, error = function(e) {
      cat("Error in hurdle model fitting:", e$message, "\n")
      return(NULL)
    })
  }
  
  ################################## Mapping Function ##################################
  
  plot_validation_maps <- function(validation_result, coord_system = "latlon") {
    
    if(is.null(validation_result) || is.null(validation_result$test_results)) {
      stop("No test results found for mapping")
    }
    
    plot_data <- validation_result$test_results
    
    # Print available columns for debugging
    cat("Available columns in test results:\n")
    print(colnames(plot_data))
    
    # Determine coordinate columns with fallbacks
    if(coord_system == "latlon") {
      x_candidates <- c("lon", "longitude", "Longitude", "LON")
      y_candidates <- c("lat", "latitude", "Latitude", "LAT")
      x_lab <- "Longitude"
      y_lab <- "Latitude"
    } else {
      x_candidates <- c("x_utm", "X_utm", "utm_x", "UTM_X")
      y_candidates <- c("y_utm", "Y_utm", "utm_y", "UTM_Y")
      x_lab <- "UTM X"
      y_lab <- "UTM Y"
    }
    
    # Find the actual column names
    x_col <- x_candidates[x_candidates %in% colnames(plot_data)]
    y_col <- y_candidates[y_candidates %in% colnames(plot_data)]
    
    if(length(x_col) == 0 || length(y_col) == 0) {
      cat("Could not find coordinate columns. Trying the other coordinate system...\n")
      
      # Try the other coordinate system
      if(coord_system == "latlon") {
        x_candidates <- c("x_utm", "X_utm", "utm_x", "UTM_X")
        y_candidates <- c("y_utm", "Y_utm", "utm_y", "UTM_Y")
        x_lab <- "UTM X"
        y_lab <- "UTM Y"
      } else {
        x_candidates <- c("lon", "longitude", "Longitude", "LON")
        y_candidates <- c("lat", "latitude", "Latitude", "LAT")
        x_lab <- "Longitude"
        y_lab <- "Latitude"
      }
      
      x_col <- x_candidates[x_candidates %in% colnames(plot_data)]
      y_col <- y_candidates[y_candidates %in% colnames(plot_data)]
    }
    
    if(length(x_col) == 0 || length(y_col) == 0) {
      stop("Could not find any coordinate columns. Available columns: ", paste(colnames(plot_data), collapse = ", "))
    }
    
    # Use first match
    x_col <- x_col[1]
    y_col <- y_col[1]
    
    cat("Using coordinates:", x_col, "and", y_col, "\n")
    
    # Common theme for maps
    map_theme <- theme_minimal() +
      theme(
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)
      )
    
    # Get response variable name
    response_var <- validation_result$response_variable
    
    # Calculate shared color scale limits
    actual_values <- plot_data[[response_var]]
    predicted_values <- plot_data$predicted
    shared_min <- min(c(actual_values, predicted_values), na.rm = TRUE)
    shared_max <- max(c(actual_values, predicted_values), na.rm = TRUE)
    
    cat("Shared color scale range:", round(shared_min, 2), "to", round(shared_max, 2), "\n")
    
    # Plot 1: Actual cover values
    p1 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = response_var)) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_viridis_c(name = "Actual", option = "plasma", limits = c(shared_min, shared_max)) +
      labs(title = paste("Actual", response_var, "(Test Set)"), 
           x = x_lab, y = y_lab) +
      map_theme
    
    # Plot 2: Predicted cover values
    p2 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = "predicted")) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_viridis_c(name = "Predicted", option = "plasma", limits = c(shared_min, shared_max)) +
      labs(title = paste("Predicted", response_var), 
           x = x_lab, y = y_lab) +
      map_theme
    
    # Plot 3: Residuals (actual - predicted)
    p3 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = "residuals")) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_gradient2(name = "Residual", 
                            low = "red", mid = "white", high = "blue",
                            midpoint = 0) +
      labs(title = "Residuals (Actual - Predicted)", 
           x = x_lab, y = y_lab) +
      map_theme
    
    # Plot 4: Predicted vs Actual scatter plot
    p4 <- ggplot(plot_data, aes_string(x = "predicted", y = response_var)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      labs(title = "Predicted vs Actual", 
           x = paste("Predicted", response_var), 
           y = paste("Actual", response_var)) +
      theme_minimal()
    
    # Arrange plots
    grid.arrange(p1, p2, p3, p4, ncol = 2,
                 top = paste("Validation Results -", 
                             validation_result$family, "model"))
    
    # Return the individual plots as well
    invisible(list(actual = p1, predicted = p2, residuals = p3, scatter = p4))
  }
  
  ################################## Usage Examples ##################################
  
  # Example 1: Standard single-model validation
  result <- gam_train_test_split(gam_results[["Montastraea"]]$complex,
                                 gam_results_entry = gam_results[["Montastraea"]],
                                 train_percent = 80)
  
  # Example 2: Hurdle model validation (uses presence + abundance models)
  result_hurdle <- gam_train_test_split(gam_model = NULL,  # Not used for hurdle
                                        gam_results_entry = gam_results[["Montastraea"]],
                                        train_percent = 80,
                                        use_hurdle = TRUE, seed = sample(1:10000, 1))
  plot_validation_maps(result_hurdle)
  
  # Example 3: Compare single vs hurdle approaches
  # single_result <- gam_train_test_split(gam_results[["Agaricia"]]$complex, 
  #                                      gam_results_entry = gam_results[["Agaricia"]])
  # hurdle_result <- gam_train_test_split(gam_model = NULL, 
  #                                      gam_results_entry = gam_results[["Agaricia"]], 
  #                                      use_hurdle = TRUE)
  # 
  # cat("Single model R²:", round(single_result$metrics$r_squared, 3), "\n")
  # cat("Hurdle model R²:", round(hurdle_result$metrics$r_squared, 3), "\n")
  # cat("Single model predicted zeros:", sum(single_result$test_results$predicted == 0), "\n")
  # cat("Hurdle model predicted zeros:", hurdle_result$metrics$predicted_zeros, "\n")
  
  # Example 4: Generate summary table for multiple species
  # results_summary <- data.frame(
  #   Species = character(),
  #   Model_Type = character(),
  #   Pearson_r = numeric(),
  #   P_value = numeric(),
  #   Significance = character(),
  #   R_squared = numeric(),
  #   RMSE = numeric(),
  #   stringsAsFactors = FALSE
  # )
  # 
  # for(sp in names(gam_results)) {
  #   for(model_type in c("complex", "reduced")) {
  #     result <- gam_train_test_split(gam_results[[sp]][[model_type]], 
  #                                   gam_results_entry = gam_results[[sp]], 
  #                                   train_percent = 80)
  #     if(!is.null(result) && !result$is_binary) {
  #       results_summary <- rbind(results_summary, data.frame(
  #         Species = sp,
  #         Model_Type = model_type,
  #         Pearson_r = round(result$metrics$pearson_r, 3),
  #         P_value = result$metrics$pearson_p,
  #         Significance = result$metrics$significance,
  #         R_squared = round(result$metrics$r_squared, 3),
  #         RMSE = round(result$metrics$rmse, 3)
  #       ))
  #     }
  #   }
  # }
  # print(results_summary)
  ################################## Save objects/workspace ##################################
  
  # #updated way to handle saving of new objects
  # save_new_objects("output/GAMs", existing_objects)
  