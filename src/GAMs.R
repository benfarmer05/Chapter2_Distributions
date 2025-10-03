  
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
  library(tictoc)
  library(pROC)
  library(leaflet)
  
  # library(gam.hp)
  
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
  
  
  ################################## functions ##################################
  
  # Define variables that are consistent across all species
  complexity_vars <- c("depth", "aspect", "slope", "TPI", "planform_curv", 
                       "max_Hsig", "dir_at_max_hsig", "mean_Hsig",  "mean_SST",
                       "range_SST", "year", "date", "lat", "lon", "cover", "mean_PAR",
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
  add_env_variables <- function(species_data, env_data) {
    species_data$depth_bathy <- env_data$depth
    species_data$aspect <- env_data$aspect
    species_data$slope <- env_data$slope
    species_data$complexity <- env_data$complexity
    species_data$TPI <- env_data$TPI
    species_data$VRM <- env_data$VRM
    species_data$planform_curv <- env_data$planform_curv
    species_data$SAPA <- env_data$SAPA
    species_data$max_Hsig <- env_data$max_Hsig
    species_data$dir_at_max_hsig <- env_data$dir_at_max_hsig
    species_data$mean_Hsig <- env_data$mean_Hsig
    species_data$mean_SST <- env_data$mean_SST
    species_data$range_SST <- env_data$range_SST
    species_data$mean_PAR <- env_data$mean_PAR
    species_data$mean_chla <- env_data$mean_chla
    species_data$mean_kd490 <- env_data$mean_kd490
    species_data$mean_spm <- env_data$mean_spm
    species_data$dist_to_land <- env_data$dist_to_land
    species_data$dist_to_deep <- env_data$dist_to_deep
    species_data$max_BOV <- env_data$max_BOV
    species_data$year <- year(species_data$date)
    
    return(species_data)
  }
  
  # Function to create correlation matrix
  create_correlation_matrix <- function(model_data_complex) {
    complexity_matrix <- model_data_complex[, c("depth_bathy", "TPI", "slope", "planform_curv",
                                                "range_SST", "mean_SST", "dir_at_max_hsig", "max_Hsig",
                                                "mean_Hsig", "year", "mean_PAR", "mean_chla", "mean_kd490",
                                                "mean_spm", "dist_to_land", "dist_to_deep", "max_BOV")]
    
    cor_matrix <- cor(complexity_matrix, use = "complete.obs")
    return(cor_matrix)
  }
  
  
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
                   mean_hsig_raster, mean_sst_raster, range_sst_raster, mean_par_raster,
                   mean_chlor_a_raster, mean_kd490_raster, mean_spm_raster, dist_to_land_raster,
                   distance_to_deep_raster, bov_full)
  
  names(env_complex) <- c("depth", "aspect", "slope", "complexity", "TPI", "VRM", "planform_curv",
                          "SAPA", "max_Hsig", "dir_at_max_hsig", "mean_Hsig", "mean_SST",
                          "range_SST", "mean_PAR", "mean_chla", "mean_kd490", "mean_spm", "dist_to_land",
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

  # ################################## SIMPLE SITE GAMs ##################################
  # 
  # # PLOT NAs - diagnostic section
  # # Create a logical vector for NA locations
  # na_mask <- is.na(site_data$depth_bathy) |
  #   is.na(site_data$slope)
  # 
  # # Get coordinates for plotting
  # coords_na <- site_data[na_mask, c("x_utm", "y_utm")]
  # coords_valid <- site_data[!na_mask, c("x_utm", "y_utm")]
  # 
  # # Define plot extent options (uncomment as needed)
  # # plot_extents = ext(260000, 310000, 2030000, 2062000) #for investigating north drop
  # # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # # plot_extents = ext(341000, 379000, 2057000, 2078000) # for investigating Anegada
  # # plot_extents = ext(294000, 350000, 1950000, 1975000) #for investigating St Croix
  # # plot_extents = ext(280000, 320000, 2000000, 2040000) #for investigating St Thomas
  # plot_extents = ext(240000, 280000, 2000000, 2040000) #for investigating Mona Island
  # # plot_extents = ext(210000, 260000, 1995000, 2050000) #for investigating PR East
  # # plot_extents = ext(55000, 230000, 2030000, 2058000) #for investigating PR North
  # # plot_extents = ext(20000, 80000, 1970000, 2060000) #for investigating PR West
  # # plot_extents = ext(52000, 230000, 1970000, 2020000) #for investigating PR South
  # 
  # # Create plot
  # bathy_final_clamp <- clamp(bathy_final, lower = -50, upper = 0)
  # plot(bathy_final,
  #      col = cmocean("deep")(100),
  #      ext = plot_extents,
  #      main = "Species locations over raster\n(Red = NA values, Black = Valid values)")
  # 
  # # Add points - valid locations in black
  # points(coords_valid$x_utm, coords_valid$y_utm,
  #        col = "black", pch = 16, cex = 0.5)
  # 
  # # Add points - NA locations in red
  # points(coords_na$x_utm, coords_na$y_utm,
  #        col = "red", pch = 16, cex = 0.5)
  # 
  # # Add legend
  # legend("topright",
  #        legend = c("Valid data", "NA values"),
  #        col = c("black", "red"),
  #        pch = 16,
  #        cex = 0.8)
  # 
  # # Print summary
  # cat("Total points:", nrow(site_data), "\n")
  # cat("Valid points:", sum(!na_mask), "\n")
  # cat("NA points:", sum(na_mask), "\n")
  # cat("Percentage with NAs:", round(sum(na_mask)/nrow(site_data)*100, 1), "%\n")
  # 
  # # Filter data for simple models
  # model_data_filtered <- site_data %>%
  #   filter(depth_bathy >= -60)
  # 
  # cat("Remaining datasets:", paste(unique(model_data_filtered$dataset), collapse = ", "), "\n")
  # 
  # psu_with_pr <- model_data_filtered$PSU[grepl("_PR", model_data_filtered$PSU)]
  # if(length(psu_with_pr) > 0) {
  #   cat("PSUs with '_PR' still remaining:", paste(unique(psu_with_pr), collapse = ", "), "\n")
  # } else {
  #   cat("No PSUs with '_PR' remaining\n")
  # }
  # 
  # model_data <- model_data_filtered[complete.cases(model_data_filtered[, c("depth_bathy", "slope", "cover")]), ]
  # 
  # # Check the new cover distribution
  # summary(model_data$cover)
  # prop_zeros_new <- sum(model_data$cover == 0) / nrow(model_data)
  # cat("New proportion of zeros:", round(prop_zeros_new, 3), "\n")
  # 
  # # Compare with original
  # prop_zeros_orig <- sum(site_data$cover == 0, na.rm = TRUE) / nrow(site_data)
  # cat("Original proportion of zeros:", round(prop_zeros_orig, 3), "\n")
  # 
  # # Option 1: Tweedie distribution (good for zero-inflated continuous data)
  # gam_tweedie <- gam(cover ~ s(depth_bathy) + s(slope),
  #                    data = model_data,
  #                    family = tw())
  # 
  # # Option 2: If cover is bounded (0-100%), use beta regression with zeros
  # # First, transform cover to (0,1) range if it's percentage
  # if(max(model_data$cover, na.rm = TRUE) > 1) {
  #   model_data$cover_prop <- model_data$cover / 100
  # } else {
  #   model_data$cover_prop <- model_data$cover
  # }
  # 
  # # Beta regression (handles zeros with adjustment)
  # gam_beta <- gam(cover_prop ~ s(depth_bathy) + s(slope),
  #                 data = model_data[model_data$cover_prop > 0, ],  # exclude zeros
  #                 family = betar())
  # 
  # # Option 3: Two-part model (hurdle model)
  # # Part 1: Presence/absence
  # model_data$present <- ifelse(model_data$cover > 0, 1, 0)
  # gam_presence <- gam(present ~ s(depth_bathy) + s(slope),
  #                     data = model_data,
  #                     family = binomial())
  # 
  # # Part 2: Abundance given presence
  # gam_abundance <- gam(cover ~ s(depth_bathy) + s(slope),
  #                      data = model_data[model_data$cover > 0, ],
  #                      family = Gamma(link = "log"))
  # 
  # # Check model summaries
  # summary(gam_tweedie)
  # summary(gam_beta)
  # summary(gam_presence)
  # summary(gam_abundance)
  # 
  # # Plot results
  # par(mfrow = c(2, 2))
  # plot(gam_tweedie, pages = 1, main = "Tweedie Model")
  # plot(gam_presence, pages = 1, main = "Presence Model")
  # 
  # # Plot both relationships on one page
  # par(mfrow = c(1, 2))
  # plot(gam_tweedie, select = 1, main = "Bathymetry effect on coral cover",
  #      xlab = "Bathymetry (m)", ylab = "Smooth term")
  # plot(gam_tweedie, select = 2, main = "Slope effect on coral cover",
  #      xlab = "Slope", ylab = "Smooth term")
  # 
  # # Reset plotting parameters
  # par(mfrow = c(1, 1))
  # 
  # # Create prediction data for plotting
  # bathy_range <- seq(min(model_data$depth_bathy, na.rm = TRUE),
  #                    max(model_data$depth_bathy, na.rm = TRUE),
  #                    length.out = 100)
  # slope_range <- seq(min(model_data$slope, na.rm = TRUE),
  #                    max(model_data$slope, na.rm = TRUE),
  #                    length.out = 100)
  # 
  # # Predictions for bathymetry (holding slope at median)
  # pred_data_bathy <- data.frame(
  #   depth_bathy = bathy_range,
  #   slope = median(model_data$slope, na.rm = TRUE)
  # )
  # 
  # # Predictions for slope (holding bathymetry at median)
  # pred_data_slope <- data.frame(
  #   depth_bathy = median(model_data$depth_bathy, na.rm = TRUE),
  #   slope = slope_range
  # )
  # 
  # # Get predictions with standard errors
  # pred_bathy <- predict(gam_tweedie, pred_data_bathy, se.fit = TRUE, type = "response")
  # pred_slope <- predict(gam_tweedie, pred_data_slope, se.fit = TRUE, type = "response")
  # 
  # # Create plots
  # par(mfrow = c(1, 2))
  # 
  # # Bathymetry plot
  # plot(pred_data_bathy$depth_bathy, pred_bathy$fit, type = "l",
  #      xlab = "Bathymetry (m)", ylab = "Predicted coral cover (%)",
  #      main = "Coral cover vs Bathymetry", lwd = 2, col = "blue")
  # lines(pred_data_bathy$depth_bathy, pred_bathy$fit + 1.96*pred_bathy$se.fit, lty = 2, col = "blue")
  # lines(pred_data_bathy$depth_bathy, pred_bathy$fit - 1.96*pred_bathy$se.fit, lty = 2, col = "blue")
  # 
  # # Add raw data points
  # points(model_data$depth_bathy, model_data$cover, pch = 16, cex = 0.3, col = "gray60")
  # 
  # # Slope plot
  # plot(pred_data_slope$slope, pred_slope$fit, type = "l",
  #      xlab = "Slope", ylab = "Predicted coral cover (%)",
  #      main = "Coral cover vs Slope", lwd = 2, col = "red")
  # lines(pred_data_slope$slope, pred_slope$fit + 1.96*pred_slope$se.fit, lty = 2, col = "red")
  # lines(pred_data_slope$slope, pred_slope$fit - 1.96*pred_slope$se.fit, lty = 2, col = "red")
  # 
  # # Add raw data points
  # points(model_data$slope, model_data$cover, pch = 16, cex = 0.3, col = "gray60")
  # 
  # par(mfrow = c(1, 1))
  # 
  # # Alternative ggplot version
  # # Create prediction dataframes
  # pred_bathy_df <- data.frame(
  #   depth_bathy = pred_data_bathy$depth_bathy,
  #   fitted = pred_bathy$fit,
  #   lower = pred_bathy$fit - 1.96*pred_bathy$se.fit,
  #   upper = pred_bathy$fit + 1.96*pred_bathy$se.fit
  # )
  # 
  # pred_slope_df <- data.frame(
  #   slope = pred_data_slope$slope,
  #   fitted = pred_slope$fit,
  #   lower = pred_slope$fit - 1.96*pred_slope$se.fit,
  #   upper = pred_slope$fit + 1.96*pred_slope$se.fit
  # )
  # 
  # # Bathymetry plot with ggplot
  # p1 <- ggplot(pred_bathy_df, aes(x = depth_bathy)) +
  #   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "blue") +
  #   geom_line(aes(y = fitted), color = "blue", size = 1) +
  #   geom_point(data = model_data, aes(x = depth_bathy, y = cover),
  #              alpha = 0.3, size = 0.5) +
  #   labs(x = "Bathymetry (m)", y = "Predicted coral cover (%)",
  #        title = "Coral cover vs Bathymetry") +
  #   theme_minimal()
  # 
  # # Slope plot with ggplot
  # p2 <- ggplot(pred_slope_df, aes(x = slope)) +
  #   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "red") +
  #   geom_line(aes(y = fitted), color = "red", size = 1) +
  #   geom_point(data = model_data, aes(x = slope, y = cover),
  #              alpha = 0.3, size = 0.5) +
  #   labs(x = "Slope", y = "Predicted coral cover (%)",
  #        title = "Coral cover vs Slope") +
  #   theme_minimal()
  # 
  # # Display both plots
  # grid.arrange(p1, p2, ncol = 2)
  # 
  # ################################## COMPLEX SITE GAMs ##################################
  # 
  # # Prepare data for complex models using PSU-averaged data
  # site_model_data_filtered <- site_data %>%
  #   filter(depth_bathy >= -60)
  # 
  # # Extract all environmental values at species locations
  # site_species_env_complex <- terra::extract(env_complex,
  #                                            cbind(site_model_data_filtered$x_utm,
  #                                                  site_model_data_filtered$y_utm))
  # 
  # # Add all variables to your filtered dataset
  # site_model_data_filtered$depth_bathy <- site_species_env_complex$depth
  # site_model_data_filtered$aspect <- site_species_env_complex$aspect
  # site_model_data_filtered$slope <- site_species_env_complex$slope
  # site_model_data_filtered$complexity <- site_species_env_complex$complexity
  # site_model_data_filtered$TPI <- site_species_env_complex$TPI
  # site_model_data_filtered$VRM <- site_species_env_complex$VRM
  # site_model_data_filtered$planform_curv <- site_species_env_complex$planform_curv
  # site_model_data_filtered$SAPA <- site_species_env_complex$SAPA
  # site_model_data_filtered$max_Hsig <- site_species_env_complex$max_Hsig
  # site_model_data_filtered$dir_at_max_hsig <- site_species_env_complex$dir_at_max_hsig
  # site_model_data_filtered$mean_Hsig <- site_species_env_complex$mean_Hsig
  # site_model_data_filtered$mean_SST <- site_species_env_complex$mean_SST
  # site_model_data_filtered$range_SST <- site_species_env_complex$range_SST
  # site_model_data_filtered$mean_PAR <- site_species_env_complex$mean_PAR
  # site_model_data_filtered$mean_chla <- site_species_env_complex$mean_chla
  # site_model_data_filtered$mean_kd490 <- site_species_env_complex$mean_kd490
  # site_model_data_filtered$mean_spm <- site_species_env_complex$mean_spm
  # site_model_data_filtered$dist_to_land <- site_species_env_complex$dist_to_land
  # site_model_data_filtered$dist_to_deep <- site_species_env_complex$dist_to_deep
  # site_model_data_filtered$max_BOV <- site_species_env_complex$max_BOV
  # site_model_data_filtered$year <- year(site_model_data_filtered$date)
  # 
  # # Create complete cases dataset with all variables
  # site_model_data_complex <- site_model_data_filtered[complete.cases(site_model_data_filtered[, complexity_vars]), ]
  # 
  # # Check how much data you have left
  # cat("Site observations with all complexity variables:", nrow(site_model_data_complex), "\n")
  # 
  # # Fit expanded GAM models
  # site_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
  #                           s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
  #                           mean_Hsig + s(mean_PAR) + s(mean_chla) + s(dist_to_land) +
  #                           s(dist_to_deep) + s(max_BOV),
  #                         data = site_model_data_complex,
  #                         family = tw())
  # 
  # # Check the results
  # summary(site_gam_all_tweedie)
  # AIC(site_gam_all_tweedie)
  # gam.check(site_gam_all_tweedie)
  # draw(site_gam_all_tweedie)
  # 
  # # Plot the relationships
  # plot(site_gam_all_tweedie, pages = 2)
  # 
  # # Compare model performance
  # cat("Simple model AIC:", AIC(gam_tweedie), "\n")
  # cat("Complex model AIC:", AIC(site_gam_all_tweedie), "\n")
  # 
  # # Check correlations between variables
  # site_cor_matrix <- create_correlation_matrix(site_model_data_complex)
  # print(round(site_cor_matrix, 2))
  # 
  # # Two-part model with complexity
  # site_model_data_complex$present <- ifelse(site_model_data_complex$cover > 0, 1, 0)
  # 
  # site_gam_presence_binom <- gam(present ~ s(depth_bathy) +
  #                                    s(complexity, k = 12) +
  #                                    max_Hsig +
  #                                    s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
  #                                    s(dist_to_deep, k = 12),
  #                                  data = site_model_data_complex,
  #                                  family = binomial())
  # 
  # site_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                     s(TPI) +
  #                                     s(mean_SST, k = 12) + s(dir_at_max_hsig, bs = 'cc') +
  #                                     s(mean_PAR, k = 12) + s(dist_to_land) +
  #                                     s(dist_to_deep) + s(max_BOV),
  #                                   data = site_model_data_complex[site_model_data_complex$cover > 0, ],
  #                                   family = Gamma(link = "log"))
  # 
  # summary(site_gam_presence_binom)
  # summary(site_gam_abundance_gamma)
  # AIC(site_gam_presence_binom)
  # AIC(site_gam_abundance_gamma)
  # 
  # gam.check(site_gam_presence_binom)
  # gam.check(site_gam_abundance_gamma)
  # 
  # draw(site_gam_presence_binom)
  # draw(site_gam_abundance_gamma)
  # 
  # concurvity(site_gam_presence_binom)
  # concurvity(site_gam_abundance_gamma)
  # 
  # # Create a reduced model with less correlated variables
  # site_gam_reduced <- gam(cover ~ s(depth_bathy) + TPI + s(planform_curv),
  #                         data = site_model_data_complex,
  #                         family = tw())
  # 
  # summary(site_gam_reduced)
  # plot(site_gam_reduced, pages = 1)
  # 
  # 
  # 
  ################################## SPP DATA SETUP ##################################
  
  #assign lat/lon information to occurrences, to facilitate extraction of environmental data at those
  #   lat/lons
  spp_data <- combined_benthic_data_averaged
  species_df <- as.data.frame(combined_benthic_data_averaged)
  species_coords <- vect(species_df,
                         geom = c("lon", "lat"),
                         crs = "EPSG:4326")  # WGS84 geographic
  species_coords_utm <- project(species_coords, crs(bathy_final))
  utm_coords <- as.data.frame(geom(species_coords_utm)[, c("x", "y")])
  spp_data$x_utm <- utm_coords$x
  spp_data$y_utm <- utm_coords$y
  
  #extract environmental data and apply it to 'spp_data'
  depth_filter <- terra::extract(bathy_final,
                                       cbind(spp_data$x_utm,
                                             spp_data$y_utm))
  names(depth_filter) <- "depth_bathy"
  spp_data$depth_bathy <- depth_filter$depth_bathy
  
  #create generic stack of environmental data at PSU's (using agaricia, but identical for any species)
  #   NOTE - here, I am filtering to <60 m depth. this is the place to adjust if wanting a deeper or
  #           shallower modeling limit
  spp_data = spp_data %>%
    filter(depth_bathy >= -60)
  filler <- spp_data %>%
    filter(grepl("Agaricia", spp))
  variables_at_PSUs <- extract_env_data(filler)
  variables_at_PSUs <- variables_at_PSUs %>% 
    filter(depth > -60)
  rm(filler)
  
  ################################## AGARICIA ##################################
  
  # NOTE - TPI vs. dist to land difficult to pin down as the one to keep for presence model
  #           -   dist to deep, meanSST, and meanchla seem to be tricky for abundance model
  
  agaricia_model_data = spp_data %>%
    filter(grepl("Agaricia", spp))
  agaricia_model_data = add_env_variables(agaricia_model_data, variables_at_PSUs)
  
  # # Fit expanded GAM models
  # agaricia_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
  #                               s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
  #                               s(mean_Hsig) + s(mean_PAR) + s(mean_chla) + s(dist_to_land) +
  #                               s(dist_to_deep) + s(max_BOV),
  #                             data = agaricia_model_data_complex, family = tw())
  # 
  # # Check the results
  # summary(agaricia_gam_all_tweedie)
  # AIC(agaricia_gam_all_tweedie)
  # gam.check(agaricia_gam_all_tweedie)
  # plot(agaricia_gam_all_tweedie, pages = 2)
  # draw(agaricia_gam_all_tweedie)
  
  # Two-part model with complexity
  agaricia_model_data$present <- ifelse(agaricia_model_data$cover > 0, 1, 0)
  
  
  n_absent <- sum(agaricia_model_data$present == 0)
  n_present <- sum(agaricia_model_data$present == 1)
  balance_ratio <- n_absent / n_present
  
  weights_vec <- ifelse(agaricia_model_data$present == 1, balance_ratio, 1)
  
  # tic()
  # agaricia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = agaricia_model_data,
  #                             select = TRUE,
  #                             family = binomial())
  # toc()
  # agaricia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                      s(complexity) +
  #                                      s(dir_at_max_hsig, bs = 'cc') + mean_Hsig +
  #                                      s(mean_SST) + s(mean_kd490) +
  #                                      s(max_BOV) +
  #                                      s(dist_to_land),
  #                                    data = agaricia_model_data,
  #                                    select = TRUE,
  #                                    family = binomial())
  # # non-weighted version; whittled down with worst concurvity
  # agaricia_gam_presence_binom <- gam(present ~ s(depth_bathy) + VRM + s(aspect, bs = 'cc') +
  #                                        s(TPI, k = 25) +
  #                                        s(mean_SST, k = 20) + s(dir_at_max_hsig, bs = 'cc') +
  #                                        s(max_BOV),
  #                                      data = agaricia_model_data,
  #                                    select = TRUE,
  #                                      family = binomial())
  # # weighted version; whittled down with observed/estimate concurvity (seems BEST so far)
  # # NOTE - dropped SAPA, TPI, spm, and kd490 b/c of extremes in partial effect curves
  # #           - also dropped distance to land and deep, because of weird predictions
  # agaricia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                      s(complexity) + s(planform_curv) + s(slope) +
  #                                      s(dir_at_max_hsig, bs = 'cc') + s(max_BOV) +
  #                                      s(mean_SST) + s(mean_PAR) + s(mean_chla) +
  #                                      s(range_SST),
  #                                    data = agaricia_model_data,
  #                                    weights = weights_vec,
  #                                    # select = TRUE,
  #                                    family = binomial())
  # non-weighted version; whittled down with worst observed/estimate concurvity (not done yet)
  agaricia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                       s(complexity) + s(planform_curv) +
                                       s(dir_at_max_hsig, bs = 'cc') +
                                       s(mean_SST) + s(mean_PAR) +
                                       s(max_BOV) +
                                       s(lon, lat),
                                     data = agaricia_model_data,
                                     # select = TRUE,
                                     family = binomial())
  
  # Use the model's data directly (avoids length mismatch)
  gam_pred_prob <- predict(agaricia_gam_presence_binom, type = "response")
  gam_pred_binary <- ifelse(gam_pred_prob > 0.5, 1, 0)
  
  # Get fitted values length
  observed_fitted <- agaricia_gam_presence_binom$y  # This matches the prediction length
  table(observed_fitted, gam_pred_binary)
  
  
  
  
  agaricia_model_data$cover_prop <- agaricia_model_data$cover / 100
  
  # agaricia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = agaricia_model_data[agaricia_model_data$cover > 0, ],
  #                               select = TRUE,
  #                              family = Gamma(link = "log"))
  # 
  #dist to deep, meanSST, meanchla
  agaricia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy, k = 3) +
                                          s(complexity) +
                                          s(dir_at_max_hsig, bs = 'cc') +
                                          s(mean_chla, k = 15),
                                        data = agaricia_model_data[agaricia_model_data$cover > 0, ],
                                        select = TRUE,
                                        family = Gamma(link = "log"))
  
  
  # agaricia_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                              data = agaricia_model_data[agaricia_model_data$cover_prop > 0, ],
  #                              # select = TRUE,
  #                              family = betar())
  # beta with cloglog using observed/estimate concurvity
  agaricia_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) +
                                       s(slope) +
                                       s(TPI) + s(VRM) + s(planform_curv) +
                                       s(dir_at_max_hsig, bs = 'cc') +
                                       s(mean_SST) + s(mean_PAR) + s(mean_chla, bs = 'ps'),
                                     data = agaricia_model_data[agaricia_model_data$cover_prop > 0, ],
                                     # select = TRUE,
                                     family = betar())
  
  pred <- predict(agaricia_gam_abundance_beta, type = "response")
  obs <- agaricia_gam_abundance_beta$model$cover_prop
  lim <- c(0, max(c(obs, pred)))
  plot(obs, pred, xlim = lim, ylim = lim)
  abline(0, 1, col = "red")
  r <- cor(obs, pred)
  text(0.05*max(lim), 0.95*max(lim), paste("R² =", round(r^2, 3), "\nR =", round(r, 3), ifelse(cor.test(obs, pred)$p.value < 0.05, "*", "ns")), adj = c(0, 1))
  
  
  summary(agaricia_gam_presence_binom)
  summary(agaricia_gam_abundance_gamma)
  summary(agaricia_gam_abundance_beta)
  AIC(agaricia_gam_presence_binom)
  AIC(agaricia_gam_abundance_gamma)
  AIC(agaricia_gam_abundance_beta)
  
  draw(agaricia_gam_presence_binom)
  draw(agaricia_gam_abundance_gamma)
  draw(agaricia_gam_abundance_beta)
  
  # Check if any smooths are hitting k limits
  gam.check(agaricia_gam_presence_binom)
  gam.check(agaricia_gam_abundance_gamma)
  gam.check(agaricia_gam_abundance_beta)
  
  # Look at concurvity
  concurvity(agaricia_gam_presence_binom, full = TRUE)
  concurvity(agaricia_gam_abundance_gamma, full = TRUE)
  concurvity(agaricia_gam_abundance_beta, full = TRUE)
  
  #AUC / ROC
  agaricia_roc_curve <- roc(agaricia_gam_presence_binom$model$present, 
                             fitted(agaricia_gam_presence_binom))
  auc(agaricia_roc_curve)
  plot(agaricia_roc_curve)
  
  # # Save models
  # saveRDS(agaricia_gam_presence_binom, 
  #         here("output", "output_GAMs", "agaricia_gam_presence_binom.rds"))
  # 
  # saveRDS(agaricia_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "agaricia_gam_abundance_gamma.rds"))
  
  ################################## MADRACIS ##################################
  
  madracis_model_data = spp_data %>%
    filter(grepl("Madracis", spp))
  madracis_model_data = add_env_variables(madracis_model_data, variables_at_PSUs)
  
  # Two-part model with complexity
  madracis_model_data$present <- ifelse(madracis_model_data$cover > 0, 1, 0)
  
  # tic()
  # madracis_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = madracis_model_data,
  #                             # select = TRUE,
  #                             family = binomial())
  # toc()
  madracis_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                       s(TPI) +
                                       s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
                                       s(mean_SST) + s(mean_kd490),
                                     data = madracis_model_data,
                                     select = TRUE,
                                     family = binomial())
  
  # madracis_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = madracis_model_data[madracis_model_data$cover > 0, ],
  #                               # select = TRUE,
  #                              family = Gamma(link = "log"))
  madracis_gam_abundance_gamma <- gam(cover ~ s(depth) + s(complexity) +
                                        s(mean_Hsig) +
                                        s(mean_chla) +
                                        s(dist_to_land),
                                      data = madracis_model_data[madracis_model_data$cover > 0, ],
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
  madracis_roc_curve <- roc(madracis_gam_presence_binom$model$present, 
                            fitted(madracis_gam_presence_binom))
  auc(madracis_roc_curve)
  plot(madracis_roc_curve)
  
  # # Save models
  # saveRDS(madracis_gam_presence_binom, 
  #         here("output", "output_GAMs", "madracis_gam_presence_binom.rds"))
  # 
  # saveRDS(madracis_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "madracis_gam_abundance_gamma.rds"))
  
  ################################## PORITES ##################################
  
  porites_model_data = spp_data %>%
    filter(grepl("Porites", spp))
  porites_model_data = add_env_variables(porites_model_data, variables_at_PSUs)
  
  # Two-part model with complexity
  porites_model_data$present <- ifelse(porites_model_data$cover > 0, 1, 0)
  
  n_absent <- sum(porites_model_data$present == 0)
  n_present <- sum(porites_model_data$present == 1)
  balance_ratio <- n_absent / n_present
  
  weights_vec <- ifelse(porites_model_data$present == 1, balance_ratio, 1)
  
  # tic()
  # porites_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = porites_model_data,
  #                             select = TRUE,
  #                             family = binomial())
  # toc()
  # # # non-weighted version; whittled down with worst concurvity
  # porites_gam_presence_binom <- gam(present ~ s(depth_bathy) +
  #                                     s(TPI, k = 12) +
  #                                     s(dir_at_max_hsig, bs = 'cc') +
  #                                     s(mean_kd490, k = 12) +
  #                                     s(max_BOV, k = 12),
  #                                   data = porites_model_data,
  #                                   select = TRUE,
  #                                   family = binomial())
  # non-weighted version; whittled down with observed/estimate concurvity (seems BEST so far)
  #   NOTE - dropped TPI, planform_curv and kd490 because of extreme partial effects
  porites_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                      s(complexity) + s(slope) +
                                      s(max_BOV) + s(dir_at_max_hsig, bs = 'cc') +
                                      s(mean_SST) + s(mean_spm) +
                                      s(range_SST) + s(dist_to_land),
                                    data = porites_model_data,
                                    # weights = weights_vec,
                                    # select = TRUE,
                                    family = binomial())
  
  gam_pred_prob <- predict(porites_gam_presence_binom, type = "response")
  gam_pred_binary <- ifelse(gam_pred_prob > 0.5, 1, 0)
  observed_fitted <- porites_gam_presence_binom$y  # This matches the prediction length
  table(observed_fitted, gam_pred_binary)
  
  
  porites_model_data$cover_prop <- porites_model_data$cover / 100
  
  # porites_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = porites_model_data[porites_model_data$cover > 0, ],
  #                               select = TRUE,
  #                               family = Gamma(link = "log"))
  # porites_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(dir_at_max_hsig, k = 8, bs = 'cc') +
  #                                      s(mean_kd490) +
  #                                      s(dist_to_deep),
  #                                    data = porites_model_data[porites_model_data$cover > 0, ],
  #                                    select = TRUE,
  #                                    family = Gamma(link = "log"))
  
  
  # porites_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = porites_model_data[porites_model_data$cover_prop > 0, ],
  #                               select = TRUE,
  #                               family = betar())
  # beta with cloglog using observed/estimate concurvity
  porites_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                      s(TPI) + s(VRM) + s(planform_curv) +
                                      s(mean_Hsig) +
                                      s(mean_SST) + s(mean_kd490) +
                                      s(max_BOV) +
                                      s(dist_to_land),
                                    data = porites_model_data[porites_model_data$cover_prop > 0, ],
                                    # select = TRUE,
                                    family = betar())
  
  pred <- predict(porites_gam_abundance_beta, type = "response")
  obs <- porites_gam_abundance_beta$model$cover_prop
  lim <- c(0, max(c(obs, pred)))
  plot(obs, pred, xlim = lim, ylim = lim)
  abline(0, 1, col = "red")
  r <- cor(obs, pred)
  text(0.05*max(lim), 0.95*max(lim), paste("R² =", round(r^2, 3), "\nR =", round(r, 3), ifelse(cor.test(obs, pred)$p.value < 0.05, "*", "ns")), adj = c(0, 1))
  
  
  summary(porites_gam_presence_binom)
  summary(porites_gam_abundance_gamma)
  summary(porites_gam_abundance_beta)
  AIC(porites_gam_presence_binom)
  AIC(porites_gam_abundance_gamma)
  AIC(porites_gam_abundance_beta)
  
  draw(porites_gam_presence_binom)
  draw(porites_gam_abundance_gamma)
  draw(porites_gam_abundance_beta)
  
  # Check if any smooths are hitting k limits
  gam.check(porites_gam_presence_binom)
  gam.check(porites_gam_abundance_gamma)
  gam.check(porites_gam_abundance_beta)
  
  # Look at concurvity
  concurvity(porites_gam_presence_binom, full = TRUE)
  concurvity(porites_gam_abundance_gamma, full = TRUE)
  concurvity(porites_gam_abundance_beta, full = TRUE)
  
  #AUC / ROC
  porites_roc_curve <- roc(porites_gam_presence_binom$model$present, 
                            fitted(porites_gam_presence_binom))
  auc(porites_roc_curve)
  plot(porites_roc_curve)
  
  # # Save models
  # saveRDS(porites_gam_presence_binom, 
  #         here("output", "output_GAMs", "porites_gam_presence_binom.rds"))
  # 
  # saveRDS(porites_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "porites_gam_abundance_gamma.rds"))
  
  ################################## SIDERASTREA ##################################
  
  siderastrea_model_data = spp_data %>%
    filter(grepl("Siderastrea", spp))
  siderastrea_model_data = add_env_variables(siderastrea_model_data, variables_at_PSUs)
    
  # Two-part model with complexity
  siderastrea_model_data$present <- ifelse(siderastrea_model_data$cover > 0, 1, 0)
  
  # tic()
  # siderastrea_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = siderastrea_model_data,
  #                             select = TRUE,
  #                             family = binomial())
  # toc()
  siderastrea_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                          s(complexity) +
                                          s(dir_at_max_hsig, bs = 'cc') +
                                          s(mean_kd490) +
                                          s(dist_to_deep, k = 12) + s(max_BOV),
                                        data = siderastrea_model_data,
                                        select = TRUE,
                                        family = binomial())
  
  # siderastrea_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = siderastrea_model_data[siderastrea_model_data$cover > 0, ],
  #                               select = TRUE,
  #                               family = Gamma(link = "log"))
  siderastrea_gam_abundance_gamma <- gam(cover ~ s(depth) +
                                           s(max_Hsig) + s(dir_at_max_hsig, k = 12, bs = 'cc') +
                                           s(mean_kd490) +
                                           s(dist_to_deep, k = 12),
                                         data = siderastrea_model_data[siderastrea_model_data$cover > 0, ],
                                         select = TRUE,
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
  
  #AUC / ROC
  siderastrea_roc_curve <- roc(siderastrea_gam_presence_binom$model$present, 
                           fitted(siderastrea_gam_presence_binom))
  auc(siderastrea_roc_curve)
  plot(siderastrea_roc_curve)
  
  # # Save models
  # saveRDS(siderastrea_gam_presence_binom, 
  #         here("output", "output_GAMs", "siderastrea_gam_presence_binom.rds"))
  # 
  # saveRDS(siderastrea_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "siderastrea_gam_abundance_gamma.rds"))
  
  ################################## MONTASTRAEA ##################################
  
  # NOTE - presence model is not the best. I think MCAV just struggles for some reason
  #       - abundance model could take a look as well
  
  montastraea_model_data = spp_data %>%
    filter(grepl("Montastraea", spp))
  montastraea_model_data = add_env_variables(montastraea_model_data, variables_at_PSUs)
    
  # Two-part model with complexity
  montastraea_model_data$present <- ifelse(montastraea_model_data$cover > 0, 1, 0)
  
  n_absent <- sum(montastraea_model_data$present == 0)
  n_present <- sum(montastraea_model_data$present == 1)
  balance_ratio <- n_absent / n_present
  
  weights_vec <- ifelse(montastraea_model_data$present == 1, balance_ratio, 1)
  
  # tic()
  # montastraea_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = montastraea_model_data,
  #                             select = TRUE,
  #                             family = binomial())
  # toc()
  # # non-weighted version; whittled down with worst concurvity
  # montastraea_gam_presence_binom <- gam(present ~ s(depth_bathy) +
  #                               s(mean_SST) +
  #                               s(mean_PAR) + s(mean_kd490, k = 12) +
  #                                 s(dist_to_deep, k = 20),
  #                             data = montastraea_model_data,
  #                             select = TRUE,
  #                             family = binomial())
  # weighted version; whittled down with observed/estimate concurvity (seems BEST so far)
  # NOTE - dropped kd490 for extreme partial effect curve. monitor BOV for this
  montastraea_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                          s(slope) +
                                          s(complexity) + s(planform_curv) +
                                          s(dir_at_max_hsig, bs = 'cc') +
                                          s(mean_SST) + s(mean_PAR) + s(mean_chla) +
                                          s(max_BOV) +
                                          s(range_SST) +
                                          s(dist_to_land),
                                        data = montastraea_model_data,
                                        weights = weights_vec,
                                        # select = TRUE,
                                        family = binomial())
  
  # Use the model's data directly (avoids length mismatch)
  gam_pred_prob <- predict(montastraea_gam_presence_binom, type = "response")
  gam_pred_binary <- ifelse(gam_pred_prob > 0.5, 1, 0)
  
  # Get fitted values length
  observed_fitted <- montastraea_gam_presence_binom$y  # This matches the prediction length
  table(observed_fitted, gam_pred_binary)
  
  
  montastraea_model_data$cover_prop <- montastraea_model_data$cover / 100
  
  
  # siderastrea_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = siderastrea_model_data[siderastrea_model_data$cover > 0, ],
  #                               select = TRUE,
  #                               family = Gamma(link = "log"))
  # # # gamma with worst concurvity
  # montastraea_gam_abundance_gamma <- gam(cover ~ s(depth) +
  #                                          s(VRM) +
  #                                           s(mean_Hsig) +
  #                                          s(mean_SST) + s(mean_kd490) +
  #                                          s(max_BOV, k = 12),
  #                                        data = montastraea_model_data[montastraea_model_data$cover > 0, ],
  #                                        select = TRUE,
  #                                        family = Gamma(link = "log"))
  
  
  # montastraea_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = montastraea_model_data[montastraea_model_data$cover_prop > 0, ],
  #                               select = TRUE,
  #                               family = betar())
  # # # beta with observed/estimate concurvity
  montastraea_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) +
                                          s(slope) +
                                          s(planform_curv) +
                                          s(mean_Hsig) +
                                          s(mean_SST) +
                                          s(mean_spm) + s(max_BOV),
                                        data = montastraea_model_data[montastraea_model_data$cover_prop > 0, ],
                                        # select = TRUE,
                                        family = betar())
  
  pred <- predict(montastraea_gam_abundance_beta, type = "response")
  obs <- montastraea_gam_abundance_beta$model$cover_prop
  lim <- c(0, max(c(obs, pred)))
  plot(obs, pred, xlim = lim, ylim = lim)
  abline(0, 1, col = "red")
  r <- cor(obs, pred)
  text(0.05*max(lim), 0.95*max(lim), paste("R² =", round(r^2, 3), "\nR =", round(r, 3), ifelse(cor.test(obs, pred)$p.value < 0.05, "*", "ns")), adj = c(0, 1))
  
  
  
  summary(montastraea_gam_presence_binom)
  summary(montastraea_gam_abundance_gamma)
  summary(montastraea_gam_abundance_beta)
  AIC(montastraea_gam_presence_binom)
  AIC(montastraea_gam_abundance_gamma)
  AIC(montastraea_gam_abundance_beta)
  
  draw(montastraea_gam_presence_binom)
  draw(montastraea_gam_abundance_gamma)
  draw(montastraea_gam_abundance_beta)
  
  # Check if any smooths are hitting k limits
  gam.check(montastraea_gam_presence_binom)
  gam.check(montastraea_gam_abundance_gamma)
  gam.check(montastraea_gam_abundance_beta)
  
  # Look at concurvity
  concurvity(montastraea_gam_presence_binom, full = TRUE)
  concurvity(montastraea_gam_abundance_gamma, full = TRUE)
  concurvity(montastraea_gam_abundance_beta, full = TRUE)
  
  #AUC / ROC
  montastraea_roc_curve <- roc(montastraea_gam_presence_binom$model$present, 
                               fitted(montastraea_gam_presence_binom))
  auc(montastraea_roc_curve)
  plot(montastraea_roc_curve, main = "ROC Curve for Madracis Presence Model")
  
  # # Save models
  # saveRDS(montastraea_gam_presence_binom, 
  #         here("output", "output_GAMs", "montastraea_gam_presence_binom.rds"))
  # 
  # saveRDS(montastraea_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "montastraea_gam_abundance_gamma.rds"))
  
  ################################## ORBICELLA ##################################
  
  # NOTE - dropped mean_SST from presence model b/c of severe k issues
  #         - also, dropped quite a bit from the abundance model...could maybe get higher deviance
  #             than current, but stripped it down to prevent concurvity and k issues
  
  orbicella_model_data = spp_data %>%
    filter(grepl("Orbicella", spp))
  orbicella_model_data = add_env_variables(orbicella_model_data, variables_at_PSUs)
  
  # Two-part model with complexity
  orbicella_model_data$present <- ifelse(orbicella_model_data$cover > 0, 1, 0)
  
  
  
  # n_absent <- sum(orbicella_model_data$present == 0)
  # n_present <- sum(orbicella_model_data$present == 1)
  # balance_ratio <- n_absent / n_present
  # 
  # weights_vec <- ifelse(orbicella_model_data$present == 1, balance_ratio, 1)
  
  # tic()
  # orbicella_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = orbicella_model_data,
  #                             weights = weights_vec,
  #                             # select = TRUE,
  #                             family = binomial())
  # toc()
  # # weighted version; whittled down with observed/estimate concurvity (seems BEST so far)
  # orbicella_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                       s(planform_curv) + s(SAPA) +
  #                                       s(dir_at_max_hsig, bs = 'cc') +
  #                                       s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                                       s(range_SST) +
  #                                       s(dist_to_land),
  #                                     data = orbicella_model_data,
  #                                     weights = weights_vec,
  #                                     select = TRUE,
  #                                     family = binomial())
  # # weighted version; WITH lon/lat, whittled down with observed/estimate concurvity
  # orbicella_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                       s(planform_curv) + s(SAPA) +
  #                                       s(dir_at_max_hsig, bs = 'cc') +
  #                                       s(mean_PAR) + s(mean_kd490) +
  #                                       s(range_SST) + s(lon, lat),
  #                                     data = orbicella_model_data,
  #                                     weights = weights_vec,
  #                                     # select = TRUE,
  #                                     family = binomial())
  # # non-weighted version; whittled down with worst concurvity
  # orbicella_gam_presence_binom <- gam(present ~ s(depth_bathy, k = 15) + s(complexity) +
  #                                       s(range_SST, k = 15) + s(dir_at_max_hsig, k = 12, bs = 'cc') +
  #                                       s(max_BOV),
  #                                     data = orbicella_model_data,
  #                                     # weights = weights_vec,
  #                                     select = TRUE,
  #                                     family = binomial())
  # non-weighted version; whittled down with observed/estimate concurvity
  orbicella_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                        s(planform_curv) + s(SAPA) +
                                        s(dir_at_max_hsig, bs = 'cc') +
                                        s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
                                        s(range_SST) +
                                        s(dist_to_land),
                                      data = orbicella_model_data,
                                      # weights = weights_vec,
                                      # select = TRUE,
                                      family = binomial())
  # non-weighted version WITH spatial smooth; whittled down with observed/estimate concurvity
  orbicella_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                        s(slope) + s(complexity) +
                                        s(dir_at_max_hsig, bs = 'cc') +
                                        s(range_SST) + s(lon, lat),
                                      data = orbicella_model_data,
                                      # select = TRUE,
                                      family = binomial())
  # # weighted version; whittled down with worst concurvity
  # orbicella_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                       s(planform_curv) + s(SAPA) +
  #                                       s(dir_at_max_hsig, bs = 'cc') +
  #                                       s(mean_SST) + s(mean_PAR) +
  #                                       s(range_SST) +
  #                                       s(dist_to_land),
  #                                     data = orbicella_model_data,
  #                                     weights = weights_vec,
  #                                     select = TRUE,
  #                                     family = binomial())
  
  # Use the model's data directly (avoids length mismatch)
  gam_pred_prob <- predict(orbicella_gam_presence_binom, type = "response")
  gam_pred_binary <- ifelse(gam_pred_prob > 0.5, 1, 0)
  
  # Get fitted values length
  observed_fitted <- orbicella_gam_presence_binom$y  # This matches the prediction length
  table(observed_fitted, gam_pred_binary)
  
  
  orbicella_model_data$cover_prop <- orbicella_model_data$cover / 100
  
  # orbicella_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                              data = orbicella_model_data[orbicella_model_data$cover > 0, ],
  #                              family = Gamma(link = "log"))
  # # gamma with log link; whittled down with worst concurvity
  # orbicella_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) +
  #                                        s(range_SST, k = 12) +
  #                                        s(mean_chla),
  #                                      data = orbicella_model_data[orbicella_model_data$cover > 0, ],
  #                                      select = TRUE,
  #                                      family = Gamma(link = "log"))
  # # gamma with log link; whittled down with observed/estimate concurvity (not done)
  # orbicella_gam_abundance_gamma <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                        s(VRM) +
  #                                        s(mean_Hsig) +
  #                                        s(mean_SST) +
  #                                        s(dist_to_deep) + s(max_BOV) +
  #                                        s(range_SST) +
  #                                        s(dist_to_land),
  #                                      data = orbicella_model_data[orbicella_model_data$cover_prop > 0, ],
  #                                      # select = TRUE,
  #                                      family = Gamma(link = "log"))
  # #lots of options for Gamma - inverse seems pretty good ??
  # # family = scat() ??
  
  

  # pred <- predict(orbicella_gam_abundance_gamma, type = "response")
  # obs <- orbicella_gam_abundance_gamma$model$cover
  # lim <- c(0, max(c(obs, pred)))
  # plot(obs, pred, xlim = lim, ylim = lim)
  # # plot(obs, pred, xlim = 40, ylim = 40)
  # abline(0, 1, col = "red")
  # r <- cor(obs, pred)
  # text(0.05*max(lim), 0.95*max(lim), paste("R² =", round(r^2, 3), "\nR =", round(r, 3), ifelse(cor.test(obs, pred)$p.value < 0.05, "*", "ns")), adj = c(0, 1))
  
  
  
  
  # orbicella_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                              data = orbicella_model_data[orbicella_model_data$cover_prop > 0, ],
  #                              # select = TRUE,
  #                              family = betar())
  # beta with cloglog using observed/estimate concurvity, no spatial smooth
  orbicella_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                        s(slope) +
                                        s(mean_Hsig) +
                                        s(mean_SST) + s(mean_PAR) +
                                        s(range_SST), #s(lon, lat)
                                      data = orbicella_model_data[orbicella_model_data$cover_prop > 0, ],
                                      # select = TRUE,
                                      family = betar(link = "cloglog"))
  # # beta with cloglog using observed/estimate concurvity; includes spatial smooth and drops other variables (seems unideal)
  # orbicella_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                       s(planform_curv) +
  #                                       s(max_BOV) +
  #                                       s(lon, lat, k = 120),
  #                                     data = orbicella_model_data[orbicella_model_data$cover_prop > 0, ],
  #                                     # select = TRUE,
  #                                     family = betar(link = "cloglog"))
  
  # #logit transform
  # orbicella_gam_abundance_beta <- gam(I(qlogis((cover_prop + 0.001)/(1 + 0.002))) ~ 
  #                                        s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                        s(slope) +
  #                                        s(complexity) + s(planform_curv) +
  #                                        s(mean_Hsig) +
  #                                        s(mean_SST) + s(mean_PAR) + s(mean_kd490) +
  #                                        s(range_SST) +
  #                                        s(dist_to_land),
  #                                      data = orbicella_model_data[orbicella_model_data$cover_prop > 0, ],
  #                                      family = gaussian())
  # #one-inflated beta
  # library(gamlss)
  # orbicella_gam_abundance_beta <- gamlss(cover_prop ~ pb(depth_bathy) + pb(aspect) +
  #                                          pb(slope) + pb(complexity) + pb(planform_curv) +
  #                                          pb(mean_Hsig) + pb(mean_SST) + pb(mean_PAR) + 
  #                                          pb(mean_kd490) + pb(range_SST) + pb(dist_to_land),
  #                                        data = na.omit(orbicella_model_data[orbicella_model_data$cover_prop > 0, 
  #                                                                            c("cover_prop", "depth_bathy", "aspect", "slope", 
  #                                                                              "complexity", "planform_curv", "mean_Hsig", "mean_SST", 
  #                                                                              "mean_PAR", "mean_kd490", "range_SST", "dist_to_land", 
  #                                                                              "lon", "lat")]),
  #                                        family = BEINF1())  
  # #zero-one-inflated beta
  # orbicella_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                       s(slope) + s(complexity) + s(planform_curv) +
  #                                       s(mean_Hsig) + s(mean_SST) + s(mean_PAR) + 
  #                                       s(mean_kd490) + s(range_SST) + s(dist_to_land) + 
  #                                       s(lon, lat),
  #                                     data = orbicella_model_data[orbicella_model_data$cover_prop > 0, ],
  #                                     family = zoinb())
  

  
  pred <- predict(orbicella_gam_abundance_beta, type = "response")
  obs <- orbicella_gam_abundance_beta$model$cover_prop  # Changed from 'cover' to 'cover_prop'
  lim <- c(0, max(c(obs, pred)))
  plot(obs, pred, xlim = lim, ylim = lim)
  abline(0, 1, col = "red")
  bestfitline_model = lm(pred ~ obs)
  lines(obs, predict(bestfitline_model))  # Also fixed 'bestfitline' to 'bestfitline_model'
  r <- cor(obs, pred)
  text(0.05*max(lim), 0.95*max(lim), paste("R² =", round(r^2, 3), "\nR =", round(r, 3), 
                                           ifelse(cor.test(obs, pred)$p.value < 0.05, "*", "ns")), adj = c(0, 1))  
  
  
  # #for logit version:
  # # Predictions (back-transform from logit scale)
  # pred_logit <- predict(orbicella_gam_abundance_beta, type = "response")
  # pred <- plogis(pred_logit) * (1 + 0.002) - 0.001
  # #
  # # Observed values (back-transform from logit scale)
  # obs_logit <- orbicella_gam_abundance_beta$model$`I(qlogis((cover_prop + 0.001)/(1 + 0.002)))`
  # obs <- plogis(obs_logit) * (1 + 0.002) - 0.001
  # #
  # lim <- c(0, max(c(obs, pred)))
  # plot(obs, pred, xlim = lim, ylim = lim)
  # abline(0, 1, col = "red")
  # r <- cor(obs, pred)
  # text(0.05*max(lim), 0.95*max(lim), paste("R² =", round(r^2, 3), "\nR =", round(r, 3), ifelse(cor.test(obs, pred)$p.value < 0.05, "*", "ns")), adj = c(0, 1))
  # 
  # 
  # #for BEINF1 version:
  # pred <- predict(orbicella_gam_abundance_beta, type = "response")
  # obs <- orbicella_gam_abundance_beta$y
  # lim <- c(0, max(c(obs, pred)))
  # plot(obs, pred, xlim = lim, ylim = lim)
  # abline(0, 1, col = "red")
  # r <- cor(obs, pred)
  # text(0.05*max(lim), 0.95*max(lim), paste("R² =", round(r^2, 3), "\nR =", round(r, 3), ifelse(cor.test(obs, pred)$p.value < 0.05, "*", "ns")), adj = c(0, 1))
  
  
  
  
  summary(orbicella_gam_presence_binom)
  summary(orbicella_gam_abundance_gamma)
  summary(orbicella_gam_abundance_beta)
  AIC(orbicella_gam_presence_binom)
  AIC(orbicella_gam_abundance_gamma)
  AIC(orbicella_gam_abundance_beta)
  
  draw(orbicella_gam_presence_binom)
  draw(orbicella_gam_abundance_gamma)
  draw(orbicella_gam_abundance_beta)
  
  # Check if any smooths are hitting k limits
  par(mfrow=c(2,2)); gam.check(orbicella_gam_presence_binom); par(mfrow=c(1,1))
  par(mfrow=c(2,2)); gam.check(orbicella_gam_abundance_gamma); par(mfrow=c(1,1))
  par(mfrow=c(2,2)); gam.check(orbicella_gam_abundance_beta); par(mfrow=c(1,1))
  
  # Look at concurvity
  concurvity(orbicella_gam_presence_binom, full = TRUE)
  concurvity(orbicella_gam_abundance_gamma, full = TRUE)
  concurvity(orbicella_gam_abundance_beta, full = TRUE)
  
  #AUC / ROC
  orbicella_roc_curve <- roc(orbicella_gam_presence_binom$model$present, 
                            fitted(orbicella_gam_presence_binom))
  auc(orbicella_roc_curve)
  plot(orbicella_roc_curve)
  
  
  #plot residual map
  pred <- predict(orbicella_gam_abundance_beta, type = "response")
  obs <- orbicella_gam_abundance_beta$model$cover_prop
  residuals <- obs - pred
  plot(orbicella_gam_abundance_beta$model$lon, orbicella_gam_abundance_beta$model$lat, 
       cex = abs(residuals)*5, 
       col = ifelse(residuals > 0, "red", "blue"),
       xlab = "Longitude", ylab = "Latitude")
  legend("topright", legend = c("Underprediction", "Overprediction"), 
         col = c("red", "blue"), pch = 1)  
  
  # # Save models
  # saveRDS(orbicella_gam_presence_binom, 
  #         here("output", "output_GAMs", "orbicella_gam_presence_binom.rds"))
  # 
  # saveRDS(orbicella_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "orbicella_gam_abundance_gamma.rds"))
  
  
  
  
  
  # # draft to test out ad-hoc adjustment of predictions
  # # Your current code
  # pred <- predict(orbicella_gam_abundance_beta, type = "response")
  # obs <- orbicella_gam_abundance_beta$model$cover
  # bestfitline_model <- lm(pred ~ obs)
  # 
  # # Find flexion point (where best-fit line crosses 1:1 line)
  # a <- coef(bestfitline_model)[1]  # intercept
  # b <- coef(bestfitline_model)[2]  # slope
  # x0 <- a / (1 - b)  # flexion point
  # 
  # # User-chosen parameters
  # L <- 2   # asymptote (2%)
  # k <- 500      # steepness
  # 
  # # Apply logistic transformation
  # correction_factor <- L / (1 + exp(-k * (obs - x0)))
  # pred_corrected <- pred * correction_factor
  # 
  # # Calculate R-squared
  # r_squared_original <- cor(obs, pred)^2
  # r_squared_corrected <- cor(obs, pred_corrected)^2
  # 
  # # Plot both together
  # lim <- c(0, max(c(obs, pred, pred_corrected)))
  # plot(obs, pred, xlim = lim, ylim = lim, col = "lightgray", pch = 16,
  #      main = "Original vs Corrected Predictions")
  # points(obs, pred_corrected, col = "darkblue", pch = 16, cex = 0.8)
  # abline(0, 1, col = "red", lwd = 2)
  # 
  # # Add legend
  # legend("bottomright", legend = c("Original", "Corrected"),
  #        col = c("lightgray", "darkblue"), pch = 16)
  # 
  # # Add R-squared to plot
  # text(x = lim[2] * 0.05, y = lim[2] * 0.95,
  #      labels = paste("R² (original) =", round(r_squared_original, 3)),
  #      adj = 0, cex = 0.9, col = "gray")
  # text(x = lim[2] * 0.05, y = lim[2] * 0.88,
  #      labels = paste("R² (corrected) =", round(r_squared_corrected, 3)),
  #      adj = 0, cex = 0.9, col = "darkblue")
  # 
  # cat("Flexion point (x0):", x0, "\n")
  # cat("Original slope:", b, "\n")
  # cat("Final slope:", coef(lm(pred_corrected ~ obs))[2], "\n")
  
  
  
  
  #possibly better ad hoc correction which properly works with predictions, not
  #   observations
  # Calibration with obs ~ pred regression
  # Get predictions and observations
  pred <- predict(orbicella_gam_abundance_beta, type = "response")
  obs <- orbicella_gam_abundance_beta$model$cover_prop
  
  # Regress observations on predictions (correct for calibration)
  bestfitline_model <- lm(obs ~ pred)
  
  # Get calibration parameters
  a <- coef(bestfitline_model)[1]  # intercept
  b <- coef(bestfitline_model)[2]  # slope
  
  # Find flexion point (where calibration line crosses 1:1 line)
  x0 <- a / (1 - b)
  
  # What the calibration line says predictions should be
  pred_target <- a + b * pred
  
  # User-chosen parameters
  shrinkage <- 0    # how much to pull toward calibration line (0-1)
  L_low <- 0.1        # additional multiplicative correction below x0
  L_high <- 2       # additional multiplicative correction above x0
  k <- 50             # steepness of transition
  
  # Step 1: Shrink toward calibration line (reduces variance)
  pred_shrunk <- pred + shrinkage * (pred_target - pred)
  
  # Step 2: Apply multiplicative correction (adjusts bias at extremes)
  sigmoid <- 1 / (1 + exp(-k * (pred_shrunk - x0)))
  correction_factor <- L_low + (L_high - L_low) * sigmoid
  
  pred_corrected <- pred_shrunk * correction_factor
  
  # Constrain to valid range
  pred_corrected <- pmax(0.001, pmin(0.999, pred_corrected))
  
  # Calculate R-squared
  r_squared_original <- cor(obs, pred)^2
  r_squared_corrected <- cor(obs, pred_corrected)^2
  
  # Plot both together
  lim <- c(0, max(c(obs, pred, pred_corrected)))
  plot(pred, obs, xlim = lim, ylim = lim, col = "lightgray", pch = 16,
       xlab = "Predicted", ylab = "Observed",
       main = "Original vs Corrected Predictions")
  points(pred_corrected, obs, col = "darkblue", pch = 16, cex = 0.8)
  abline(0, 1, col = "red", lwd = 2)
  abline(a, b, col = "green", lwd = 2)  # Calibration line (best-fit line)
  
  # Add legend
  legend("bottomright", legend = c("Original", "Corrected"),
         col = c("lightgray", "darkblue"), pch = 16)
  
  # Add R-squared to plot
  text(x = lim[2] * 0.05, y = lim[2] * 0.95,
       labels = paste("R² (original) =", round(r_squared_original, 3)),
       adj = 0, cex = 0.9, col = "gray")
  text(x = lim[2] * 0.05, y = lim[2] * 0.88,
       labels = paste("R² (corrected) =", round(r_squared_corrected, 3)),
       adj = 0, cex = 0.9, col = "darkblue")
  
  cat("Flexion point (x0):", x0, "\n")
  cat("Shrinkage toward calibration line:", shrinkage, "\n")
  cat("L_low:", L_low, "| L_high:", L_high, "\n")
  cat("Final slope:", coef(lm(obs ~ pred_corrected))[2], "\n")
  
  
  
  
  
  
  #possible way to optimize the above
  # ===== AGGRESSIVE BIAS-FOCUSED OPTIMIZATION =====
  library(mgcv)
  
  optimize_sigmoid_correction <- function(
    fitted_model,
    L_low_range = c(0.05, 1.5),
    L_high_range = c(0.8, 4.0),     # Allow even more aggressive high-end correction
    k_range = c(10, 300),
    shrinkage_range = c(0, 1),
    focus = c("extremes", "high_only", "balanced"),  # NEW: Where to focus correction
    plot = TRUE
  ) {
    
    pred <- predict(fitted_model, type = "response")
    obs <- fitted_model$model[[1]]
    
    bestfitline_model <- lm(obs ~ pred)
    a <- coef(bestfitline_model)[1]
    b <- coef(bestfitline_model)[2]
    x0 <- a / (1 - b)
    
    focus <- match.arg(focus)
    
    # ===== AGGRESSIVE OBJECTIVE - FOCUS ON EXTREMES =====
    sigmoid_objective <- function(params) {
      L_low <- params[1]
      L_high <- params[2]
      k <- params[3]
      shrinkage <- params[4]
      
      # Apply correction
      pred_target <- a + b * pred
      pred_shrunk <- pred + shrinkage * (pred_target - pred)
      sigmoid_val <- 1 / (1 + exp(-k * (pred_shrunk - x0)))
      correction_factor <- L_low + (L_high - L_low) * sigmoid_val
      pred_test <- pred_shrunk * correction_factor
      pred_test <- pmax(0.001, pmin(0.999, pred_test))
      
      # Calculate residuals
      residuals <- obs - pred_test
      
      # Bin predictions
      n_bins <- 20  # More bins for finer resolution
      pred_bins <- cut(pred_test, 
                       breaks = quantile(pred_test, probs = seq(0, 1, length.out = n_bins + 1)),
                       include.lowest = TRUE, labels = FALSE)
      
      bin_bias <- numeric(n_bins)
      bin_counts <- numeric(n_bins)
      
      for (i in 1:n_bins) {
        mask <- pred_bins == i & !is.na(pred_bins)
        if (sum(mask) > 3) {
          bin_bias[i] <- mean(residuals[mask])
          bin_counts[i] <- sum(mask)
        }
      }
      
      if (focus == "extremes") {
        # Heavily weight lowest 3 and highest 3 bins
        weights <- c(rep(10, 3), rep(1, n_bins - 6), rep(10, 3))
        loss <- sum(abs(bin_bias) * weights * sqrt(bin_counts + 1)) / sum(weights)
        
      } else if (focus == "high_only") {
        # Focus almost exclusively on high predictions (where you have underprediction)
        # Exponentially increasing weights toward high end
        weights <- exp(seq(0, 3, length.out = n_bins))  # Exponential increase
        loss <- sum(abs(bin_bias) * weights * sqrt(bin_counts + 1)) / sum(weights)
        
      } else {  # balanced
        # Standard weighted approach
        weights <- sqrt(bin_counts + 1)
        loss <- sum(abs(bin_bias) * weights) / sum(weights)
      }
      
      # Small penalty only if slope gets too extreme (outside 0.8-1.2)
      slope <- coef(lm(obs ~ pred_test))[2]
      if (slope < 0.8 || slope > 1.2) {
        loss <- loss + 10 * abs(slope - 1)
      }
      
      return(loss)
    }
    
    cat("Optimizing with focus on:", focus, "\n\n")
    
    # ===== AGGRESSIVE STARTING POINTS =====
    starting_points <- list(
      c(0.1, 2.0, 50, 0.0),      # Your manual values
      c(0.05, 2.5, 80, 0.1),     # More aggressive
      c(0.15, 1.8, 100, 0.2),    # Steeper transition
      c(0.08, 3.0, 60, 0.0),     # Very aggressive high
      c(0.2, 2.2, 40, 0.3),      # Moderate
      c(0.05, 3.5, 120, 0.0),    # Extremely aggressive
      c(0.1, 2.5, 150, 0.1)      # High k, aggressive
    )
    
    best_result <- NULL
    best_loss <- Inf
    
    for (i in 1:length(starting_points)) {
      result <- optim(
        par = starting_points[[i]],
        fn = sigmoid_objective,
        method = "L-BFGS-B",
        lower = c(L_low_range[1], L_high_range[1], k_range[1], shrinkage_range[1]),
        upper = c(L_low_range[2], L_high_range[2], k_range[2], shrinkage_range[2]),
        control = list(maxit = 1000, factr = 1e7)  # Allow more iterations, looser convergence
      )
      
      if (result$value < best_loss) {
        best_loss <- result$value
        best_result <- result
      }
    }
    
    L_low <- best_result$par[1]
    L_high <- best_result$par[2]
    k <- best_result$par[3]
    shrinkage <- best_result$par[4]
    
    cat("Optimized parameters:\n")
    cat("  L_low:     ", round(L_low, 4), "\n")
    cat("  L_high:    ", round(L_high, 4), "\n")
    cat("  k:         ", round(k, 2), "\n")
    cat("  Shrinkage: ", round(shrinkage, 4), "\n")
    cat("  x0:        ", round(x0, 4), "\n\n")
    
    # Create correction function
    correction_fn <- function(new_pred) {
      pred_target <- a + b * new_pred
      pred_shrunk <- new_pred + shrinkage * (pred_target - new_pred)
      sigmoid_val <- 1 / (1 + exp(-k * (pred_shrunk - x0)))
      correction_factor <- L_low + (L_high - L_low) * sigmoid_val
      corrected <- pred_shrunk * correction_factor
      pmax(0.001, pmin(0.999, corrected))
    }
    
    pred_corrected <- correction_fn(pred)
    
    # Metrics
    r2_original <- cor(obs, pred)^2
    r2_corrected <- cor(obs, pred_corrected)^2
    slope_original <- coef(lm(obs ~ pred))[2]
    slope_corrected <- coef(lm(obs ~ pred_corrected))[2]
    intercept_corrected <- coef(lm(obs ~ pred_corrected))[1]
    
    # Bias by range
    high_pred_mask <- pred_corrected > quantile(pred_corrected, 0.8)
    low_pred_mask <- pred_corrected < quantile(pred_corrected, 0.2)
    
    high_bias <- mean(obs[high_pred_mask] - pred_corrected[high_pred_mask])
    low_bias <- mean(obs[low_pred_mask] - pred_corrected[low_pred_mask])
    
    if (plot) {
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
      lim <- c(0, max(c(obs, pred, pred_corrected)))
      
      # Plot 1: Original
      plot(pred, obs, xlim = lim, ylim = lim, pch = 16, col = rgb(0, 0, 0, 0.3),
           xlab = "Predicted", ylab = "Observed", main = "Original")
      abline(0, 1, col = "red", lwd = 2)
      abline(lm(obs ~ pred), col = "black", lwd = 1, lty = 2)
      text(lim[1], lim[2] * 0.95, 
           paste0("R² = ", round(r2_original, 3), 
                  "\nSlope = ", round(slope_original, 3)),
           adj = 0, cex = 0.9)
      
      # Plot 2: Corrected
      plot(pred_corrected, obs, xlim = lim, ylim = lim, pch = 16, 
           col = rgb(0, 0.4, 0.8, 0.3),
           xlab = "Predicted (Corrected)", ylab = "Observed",
           main = paste0("Optimized (", focus, ")"))
      abline(0, 1, col = "red", lwd = 2)
      abline(lm(obs ~ pred_corrected), col = "darkblue", lwd = 1, lty = 2)
      text(lim[1], lim[2] * 0.95,
           paste0("R² = ", round(r2_corrected, 3), 
                  "\nSlope = ", round(slope_corrected, 3),
                  "\nIntercept = ", round(intercept_corrected, 3)),
           adj = 0, cex = 0.8)
      
      # Plot 3: Overlay
      plot(pred, obs, xlim = lim, ylim = lim, col = rgb(0.7, 0.7, 0.7, 0.4), pch = 16,
           xlab = "Predicted", ylab = "Observed", main = "Original vs Corrected")
      points(pred_corrected, obs, col = rgb(0, 0.3, 0.8, 0.4), pch = 16, cex = 0.8)
      abline(0, 1, col = "red", lwd = 2)
      legend("bottomright", 
             legend = c("Original", "Corrected"),
             col = c(rgb(0.7, 0.7, 0.7, 0.8), rgb(0, 0.3, 0.8, 0.8)), 
             pch = 16, cex = 0.9)
      
      # Plot 4: Correction curve - EXTENDED RANGE
      pred_seq <- seq(0, max(pred) * 1.2, length.out = 300)  # Extend beyond data
      corrected_seq <- correction_fn(pred_seq)
      
      plot(pred_seq, corrected_seq, type = "l", lwd = 3, col = "darkblue",
           xlab = "Original Prediction", ylab = "Corrected Prediction",
           main = "Correction Curve",
           xlim = c(0, max(pred_seq)), ylim = c(0, max(corrected_seq)))
      abline(0, 1, col = "red", lwd = 2, lty = 2)
      abline(v = x0, col = "gray", lty = 3, lwd = 1.5)
      
      # Add reference lines
      abline(h = seq(0, max(corrected_seq), by = 0.1), col = "gray90", lty = 1)
      abline(v = seq(0, max(pred_seq), by = 0.05), col = "gray90", lty = 1)
      
      # Show actual data range
      rect(0, 0, max(pred), max(pred_corrected), 
           border = "orange", lwd = 2, lty = 2)
      
      text(x0, max(corrected_seq) * 0.9, 
           paste0("x0 = ", round(x0, 3)), pos = 4, cex = 0.8)
      text(max(pred_seq) * 0.7, max(corrected_seq) * 0.1,
           paste0("L_low = ", round(L_low, 2), "\nL_high = ", round(L_high, 2)),
           cex = 0.8, col = "darkblue")
      
      grid()
      
      par(mfrow = c(1, 1))
    }
    
    cat("=== CALIBRATION SUMMARY ===\n")
    cat("Original  - R²:", round(r2_original, 3), 
        "| Slope:", round(slope_original, 3), "\n")
    cat("Corrected - R²:", round(r2_corrected, 3), 
        "| Slope:", round(slope_corrected, 3),
        "| Intercept:", round(intercept_corrected, 3), "\n")
    cat("\nBias in top 20%: ", round(high_bias, 4), "\n")
    cat("Bias in bottom 20%: ", round(low_bias, 4), "\n\n")
    
    list(
      correction_function = correction_fn,
      parameters = list(L_low = L_low, L_high = L_high, k = k, 
                        shrinkage = shrinkage, x0 = x0),
      r2_original = r2_original,
      r2_corrected = r2_corrected,
      slope_original = slope_original,
      slope_corrected = slope_corrected,
      high_bias = high_bias,
      low_bias = low_bias,
      pred_corrected = pred_corrected,
      observations = obs,
      predictions_original = pred
    )
  }
  
  # ===== USAGE - TRY THESE IN ORDER =====
  
  # 1. Focus on high-value areas (your main concern)
  result_high <- optimize_sigmoid_correction(
    orbicella_gam_abundance_beta,
    L_low_range = c(0.05, 0.5),   # Don't overcorrect low end
    L_high_range = c(1.5, 4.0),   # Allow very aggressive high-end correction
    k_range = c(30, 200),
    focus = "high_only",
    plot = TRUE
  )
  result_high <- optimize_sigmoid_correction(
    orbicella_gam_abundance_beta,
    L_low_range = c(0.05, 0.1),   # Don't overcorrect low end
    L_high_range = c(2, 2.5),   # Allow very aggressive high-end correction
    k_range = c(50, 200),
    focus = "high_only",
    plot = TRUE
  )
  
  # 2. Balance both extremes aggressively
  # result_extremes <- optimize_sigmoid_correction(
  #   orbicella_gam_abundance_beta,
  #   L_low_range = c(0.05, 0.8),
  #   L_high_range = c(1.5, 4.0),
  #   k_range = c(30, 200),
  #   focus = "extremes",
  #   plot = TRUE
  # )
  # result_extremes <- optimize_sigmoid_correction(
  #   orbicella_gam_abundance_beta,
  #   L_low_range = c(0.05, 0.1),
  #   L_high_range = c(2, 2.5),
  #   k_range = c(50, 200),
  #   focus = "extremes",
  #   plot = TRUE
  # )
  
  # 3. If those are too aggressive, try balanced
  result_balanced <- optimize_sigmoid_correction(
    orbicella_gam_abundance_beta,
    L_low_range = c(0.1, 1.0),
    L_high_range = c(1.2, 3.0),
    focus = "balanced",
    plot = TRUE
  )
  
  
  
  # or just run it basically deterministic
  result_high <- optimize_sigmoid_correction(
    orbicella_gam_abundance_beta,
    L_low_range = c(0.1, 0.101),   # Don't overcorrect low end
    L_high_range = c(2, 2.01),   # Allow very aggressive high-end correction
    k_range = c(50, 50.1),
    focus = "high_only",
    plot = TRUE
  )
  
  
  # Apply to new data
  # new_pred_corrected <- result_high$correction_function(new_pred)
  
  
  # #possible other models to zero in on better tails in the predictions
  # # 1. Negative Binomial (convert proportions to counts first)
  # orbicella_model_data$cover_count <- round(orbicella_model_data$cover_prop * 100)
  # orbicella_gam_negbin <- gam(cover_count ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) + s(complexity) + s(planform_curv) +
  #                               s(mean_Hsig) + s(mean_SST) + s(mean_PAR) + s(mean_kd490) +
  #                               s(range_SST) + s(dist_to_land),
  #                             data = orbicella_model_data[orbicella_model_data$cover_prop > 0, ],
  #                             select = TRUE, family = nb())
  # pred <- predict(orbicella_gam_negbin, type = "response") / 100
  # obs <- orbicella_gam_negbin$model$cover_count / 100
  # summary(orbicella_gam_negbin)
  # 
  # # 2. Quantile Regression
  # library(qgam)
  # orbicella_qgam_50 <- qgam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                             s(slope) + s(complexity) + s(planform_curv) +
  #                             s(mean_Hsig) + s(mean_SST) + s(mean_PAR) + s(mean_kd490) +
  #                             s(range_SST) + s(dist_to_land) + s(lon, lat, k = 75),
  #                           qu = 0.6, data = orbicella_model_data[orbicella_model_data$cover_prop > 0, ])
  # pred <- predict(orbicella_qgam_50, type = "response")
  # obs <- orbicella_qgam_50$model$cover_prop
  # summary(orbicella_qgam_50)
  # concurvity(orbicella_qgam_50)
  # 
  # # 3. Yeo-Johnson Transformation
  # library(car)
  # lambda <- powerTransform(cover_prop ~ 1, data = orbicella_model_data[orbicella_model_data$cover_prop > 0, ])
  # orbicella_model_data$cover_yj <- yjPower(orbicella_model_data$cover_prop, lambda$lambda)
  # orbicella_gam_yj <- gam(cover_yj ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                           s(slope) + s(complexity) + s(planform_curv) +
  #                           s(mean_Hsig) + s(mean_SST) + s(mean_PAR) + s(mean_kd490) +
  #                           s(range_SST) + s(dist_to_land),
  #                         data = orbicella_model_data[orbicella_model_data$cover_prop > 0, ],
  #                         select = TRUE, family = gaussian())
  # pred_yj <- predict(orbicella_gam_yj, type = "response")
  # pred <- yjPower(pred_yj, 1/lambda$lambda)
  # obs <- orbicella_gam_yj$model$cover_yj
  # obs <- yjPower(obs, 1/lambda$lambda)
  # summary(orbicella_gam_yj)
  # 
  # # 4. Truncated (fit to middle 80% of data)
  # middle_data <- orbicella_model_data[orbicella_model_data$cover_prop > 0.02 & orbicella_model_data$cover_prop < 0.8, ]
  # orbicella_gam_trunc <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                              s(slope) + s(complexity) + s(planform_curv) +
  #                              s(mean_Hsig) + s(mean_SST) + s(mean_PAR) + s(mean_kd490) +
  #                              s(range_SST) + s(dist_to_land),
  #                            data = middle_data, select = TRUE, family = betar(link = "cloglog"))
  # pred <- predict(orbicella_gam_trunc, type = "response")
  # obs <- orbicella_gam_trunc$model$cover_prop
  # summary(orbicella_gam_trunc)
  # 
  # # 5. Random Forest
  # library(randomForest)
  # complete_data <- orbicella_model_data[orbicella_model_data$cover_prop > 0, ]
  # complete_data <- complete_data[complete.cases(complete_data[c("cover_prop", "depth_bathy", "aspect", "slope", "complexity", 
  #                                                               "planform_curv", "mean_Hsig", "mean_SST", "mean_PAR", 
  #                                                               "mean_kd490", "range_SST", "dist_to_land")]), ]
  # 
  # orbicella_rf <- randomForest(cover_prop ~ depth_bathy + aspect + slope + complexity + 
  #                                planform_curv + mean_Hsig + mean_SST + mean_PAR + 
  #                                mean_kd490 + range_SST + dist_to_land + lon + lat,
  #                              data = complete_data,
  #                              na.action = na.omit)
  # 
  # pred <- predict(orbicella_rf)
  # obs <- complete_data$cover_prop
  # 
  # 
  # # Plot (same for all)
  # lim <- c(0, max(c(obs, pred)))
  # plot(obs, pred, xlim = lim, ylim = lim)
  # abline(0, 1, col = "red")
  # bestfitline_model = lm(pred ~ obs)
  # lines(obs, predict(bestfitline_model))
  # r <- cor(obs, pred)
  # text(0.05*max(lim), 0.95*max(lim), paste("R² =", round(r^2, 3), "\nR =", round(r, 3), ifelse(cor.test(obs, pred)$p.value < 0.05, "*", "ns")), adj = c(0, 1))
  
  
  
  
  
  
  
  # ############################## SOLENASTREA ##################################
  # 
  # # NOTE - likely too low sample size. may need to drop entirely b/c no easy group to pool with
  # 
  # solenastrea <- spp_data %>%
  #   filter(depth_bathy >= -60) %>%
  #   filter(grepl("Solenastrea", spp))
  # 
  # solenastrea_model_data_filtered <- solenastrea
  # 
  # # Extract environmental data
  # solenastrea_species_env_complex <- extract_env_data(solenastrea)
  # 
  # # Add environmental variables
  # solenastrea <- add_env_variables(solenastrea, solenastrea_species_env_complex, solenastrea_model_data_filtered)
  # 
  # solenastrea_model_data <- solenastrea
  # 
  # # Check how much data you have left
  # cat("Solenastrea observations with all complexity variables:", nrow(solenastrea_model_data), "\n")
  # 
  # # # Fit expanded GAM models
  # # solenastrea_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
  # #                                   s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
  # #                                   s(mean_Hsig) + s(mean_PAR) + s(mean_chla) + s(dist_to_land) +
  # #                                   s(dist_to_deep) + s(max_BOV),
  # #                                 data = solenastrea_model_data, family = tw())
  # # 
  # # # Check the results
  # # summary(solenastrea_gam_all_tweedie)
  # # AIC(solenastrea_gam_all_tweedie)
  # # gam.check(solenastrea_gam_all_tweedie)
  # # plot(solenastrea_gam_all_tweedie, pages = 2)
  # # draw(solenastrea_gam_all_tweedie)
  # 
  # # Create correlation matrix
  # solenastrea_cor_matrix <- create_correlation_matrix(solenastrea_model_data)
  # print(round(solenastrea_cor_matrix, 2))
  # 
  # # Two-part model with complexity
  # solenastrea_model_data$present <- ifelse(solenastrea_model_data$cover > 0, 1, 0)
  # 
  # # tic()
  # # solenastrea_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  # #                               s(slope) +
  # #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  # #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  # #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  # #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  # #                               s(range_SST) +
  # #                               s(dist_to_land),
  # #                             data = solenastrea_model_data,
  # #                             select = TRUE,
  # #                             family = binomial())
  # # toc()
  # solenastrea_gam_presence_binom <- gam(present ~ s(depth_bathy) +
  #                                         s(slope) +
  #                                         s(mean_SST, k = 15) +
  #                                         s(mean_spm, k = 15) +
  #                                         s(range_SST),
  #                                       data = solenastrea_model_data,
  #                                       select = TRUE,
  #                                       family = binomial())
  # 
  # # solenastrea_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  # #                               s(slope) +
  # #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  # #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  # #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  # #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  # #                               s(range_SST) +
  # #                               s(dist_to_land),
  # #                              data = solenastrea_model_data[solenastrea_model_data$cover > 0, ],
  # #                              family = Gamma(link = "log"))
  # solenastrea_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                          s(slope) +
  #                                          s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                                          s(range_SST) +
  #                                          s(dist_to_land),
  #                                        data = solenastrea_model_data[solenastrea_model_data$cover > 0, ],
  #                                        # select = TRUE,
  #                                        family = Gamma(link = "log"))
  # 
  # summary(solenastrea_gam_presence_binom)
  # summary(solenastrea_gam_abundance_gamma)
  # AIC(solenastrea_gam_presence_binom)
  # AIC(solenastrea_gam_abundance_gamma)
  # 
  # draw(solenastrea_gam_presence_binom)
  # draw(solenastrea_gam_abundance_gamma)
  # 
  # # Check if any smooths are hitting k limits
  # gam.check(solenastrea_gam_presence_binom)
  # gam.check(solenastrea_gam_abundance_gamma)
  # 
  # # Look at concurvity
  # concurvity(solenastrea_gam_presence_binom, full = TRUE)
  # concurvity(solenastrea_gam_abundance_gamma, full = TRUE)
  # 
  # #AUC / ROC
  # solenastrea_fitted_data <- solenastrea_gam_presence_binom$model
  # solenastrea_roc_curve <- roc(solenastrea_fitted_data$present, 
  #                           fitted(solenastrea_gam_presence_binom))
  # auc(solenastrea_roc_curve)
  # plot(solenastrea_roc_curve)
  # 
  ################################## COLPOPHYLLIA ##################################
  
  colpophyllia_model_data = spp_data %>%
    filter(grepl("Colpophyllia", spp))
  colpophyllia_model_data = add_env_variables(colpophyllia_model_data, variables_at_PSUs)
    
  # Two-part model with complexity
  colpophyllia_model_data$present <- ifelse(colpophyllia_model_data$cover > 0, 1, 0)
  
  n_absent <- sum(colpophyllia_model_data$present == 0)
  n_present <- sum(colpophyllia_model_data$present == 1)
  balance_ratio <- n_absent / n_present
  
  weights_vec <- ifelse(colpophyllia_model_data$present == 1, balance_ratio, 1)
  
  # tic()
  # colpophyllia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = colpophyllia_model_data,
  #                             select = TRUE,
  #                             family = binomial())
  # toc()
  # # # non-weighted version; whittled down with worst concurvity
  # colpophyllia_gam_presence_binom <- gam(present ~ s(depth_bathy, k = 12) +
  #                                          s(dir_at_max_hsig, bs = 'cc') +
  #                                          mean_SST + mean_kd490 +
  #                                          s(max_BOV) +
  #                                          s(range_SST),
  #                                        data = colpophyllia_model_data,
  #                                        select = TRUE,
  #                                        family = binomial())
  # weighted version; whittled down with observed/estimate concurvity (seems BEST so far)
  #   NOTE - very strange partial effects curves. attempting to drop slope, range PAR, kd490, planform, chla
  #           - possible that weighting is actually making things worse here?
  colpophyllia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                           s(complexity) +
                                           s(mean_Hsig) +
                                           s(mean_SST) +
                                           s(range_SST),
                                         data = colpophyllia_model_data,
                                         weights = weights_vec,
                                         # select = TRUE,
                                         family = binomial())
  
  gam_pred_prob <- predict(colpophyllia_gam_presence_binom, type = "response")
  gam_pred_binary <- ifelse(gam_pred_prob > 0.5, 1, 0)
  observed_fitted <- colpophyllia_gam_presence_binom$y  # This matches the prediction length
  table(observed_fitted, gam_pred_binary)
  
  
  colpophyllia_model_data$cover_prop <- colpophyllia_model_data$cover / 100
  
  # colpophyllia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                              data = colpophyllia_model_data[colpophyllia_model_data$cover > 0, ],
  #                              select = TRUE,
  #                              family = Gamma(link = "log"))
  # colpophyllia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) +
  #                                           s(slope) +
  #                                           s(range_SST),
  #                                         data = colpophyllia_model_data[colpophyllia_model_data$cover > 0, ],
  #                                         select = TRUE,
  #                                         family = Gamma(link = "log"))
  
  
  
  # colpophyllia_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = colpophyllia_model_data[colpophyllia_model_data$cover_prop > 0, ],
  #                               select = TRUE,
  #                               family = betar())
  # # # beta with observed/estimate concurvity
  colpophyllia_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) +
                                           s(TPI) +
                                           s(max_BOV) +
                                           s(mean_SST) + s(range_SST) + s(mean_PAR) +
                                           s(dist_to_land),
                                         data = colpophyllia_model_data[colpophyllia_model_data$cover_prop > 0, ],
                                         # select = TRUE,
                                         family = betar())
  
  
  pred <- predict(colpophyllia_gam_abundance_beta, type = "response")
  obs <- colpophyllia_gam_abundance_beta$model$cover_prop
  lim <- c(0, max(c(obs, pred)))
  plot(obs, pred, xlim = lim, ylim = lim)
  abline(0, 1, col = "red")
  r <- cor(obs, pred)
  text(0.05*max(lim), 0.95*max(lim), paste("R² =", round(r^2, 3), "\nR =", round(r, 3), ifelse(cor.test(obs, pred)$p.value < 0.05, "*", "ns")), adj = c(0, 1))
  
  
  
  summary(colpophyllia_gam_presence_binom)
  summary(colpophyllia_gam_abundance_gamma)
  summary(colpophyllia_gam_abundance_beta)
  AIC(colpophyllia_gam_presence_binom)
  AIC(colpophyllia_gam_abundance_gamma)
  AIC(colpophyllia_gam_abundance_beta)
  
  draw(colpophyllia_gam_presence_binom)
  draw(colpophyllia_gam_abundance_gamma)
  draw(colpophyllia_gam_abundance_beta)
  
  # Check if any smooths are hitting k limits
  gam.check(colpophyllia_gam_presence_binom)
  gam.check(colpophyllia_gam_abundance_gamma)
  gam.check(colpophyllia_gam_abundance_beta)
  
  # Look at concurvity
  concurvity(colpophyllia_gam_presence_binom, full = TRUE)
  concurvity(colpophyllia_gam_abundance_gamma, full = TRUE)
  concurvity(colpophyllia_gam_abundance_beta, full = TRUE)
  
  #AUC / ROC
  colpophyllia_roc_curve <- roc(colpophyllia_gam_presence_binom$model$present, 
                             fitted(colpophyllia_gam_presence_binom))
  auc(colpophyllia_roc_curve)
  plot(colpophyllia_roc_curve)
  
  # # Save models
  # saveRDS(colpophyllia_gam_presence_binom, 
  #         here("output", "output_GAMs", "colpophyllia_gam_presence_binom.rds"))
  # 
  # saveRDS(colpophyllia_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "colpophyllia_gam_abundance_gamma.rds"))
  
  # ################################## DENDROGYRA ##################################
  # 
  # # NOTE - presence & abundance models seems very unstable. likely need to lump with 'RARE HS'
  # 
  # dendrogyra <- spp_data %>%
  #   filter(depth_bathy >= -60) %>%
  #   filter(grepl("Dendrogyra", spp))
  # 
  # dendrogyra_model_data_filtered <- dendrogyra
  # 
  # # Extract environmental data
  # dendrogyra_species_env_complex <- extract_env_data(dendrogyra)
  # 
  # # Add environmental variables
  # dendrogyra <- add_env_variables(dendrogyra, dendrogyra_species_env_complex, dendrogyra_model_data_filtered)
  # 
  # dendrogyra_model_data <- dendrogyra
  # 
  # # Check how much data you have left
  # cat("Dendrogyra observations with all complexity variables:", nrow(dendrogyra_model_data), "\n")
  # 
  # # # Fit expanded GAM models
  # # dendrogyra_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
  # #                                 s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
  # #                                 s(mean_Hsig) + s(mean_PAR) + s(mean_chla) + s(dist_to_land) +
  # #                                 s(dist_to_deep) + s(max_BOV),
  # #                               data = dendrogyra_model_data, family = tw())
  # # 
  # # # Check the results
  # # summary(dendrogyra_gam_all_tweedie)
  # # AIC(dendrogyra_gam_all_tweedie)
  # # gam.check(dendrogyra_gam_all_tweedie)
  # # plot(dendrogyra_gam_all_tweedie, pages = 2)
  # # draw(dendrogyra_gam_all_tweedie)
  # 
  # # Create correlation matrix
  # dendrogyra_cor_matrix <- create_correlation_matrix(dendrogyra_model_data)
  # print(round(dendrogyra_cor_matrix, 2))
  # 
  # # Two-part model with complexity
  # dendrogyra_model_data$present <- ifelse(dendrogyra_model_data$cover > 0, 1, 0)
  # 
  # # tic()
  # # dendrogyra_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  # #                               s(slope) +
  # #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  # #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  # #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  # #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  # #                               s(range_SST) +
  # #                               s(dist_to_land),
  # #                             data = dendrogyra_model_data,
  # #                             select = TRUE,
  # #                             family = binomial())
  # # toc()
  # dendrogyra_gam_presence_binom <- gam(present ~ s(depth_bathy) +
  #                                        s(SAPA) +
  #                                        s(range_SST) +
  #                                        max_BOV,
  #                                      data = dendrogyra_model_data,
  #                                      # select = TRUE,
  #                                      family = binomial())
  # 
  # # dendrogyra_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  # #                               s(slope) +
  # #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  # #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  # #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  # #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  # #                               s(range_SST) +
  # #                               s(dist_to_land),
  # #                              data = dendrogyra_model_data[dendrogyra_model_data$cover > 0, ],
  # #                              select = TRUE,
  # #                              family = Gamma(link = "log"))
  # dendrogyra_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                         s(slope) +
  #                                         s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                                         s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                                         s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                                         s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                                         s(range_SST) +
  #                                         s(dist_to_land),
  #                                       data = dendrogyra_model_data[dendrogyra_model_data$cover > 0, ],
  #                                       # select = TRUE,
  #                                       family = Gamma(link = "log"))
  # dendrogyra_gam_abundance_gamma <- gam(cover ~ s(slope) + TPI + VRM + planform_curv +
  #                                           mean_spm +
  #                                           s(dist_to_deep),
  #                                         data = dendrogyra_model_data[dendrogyra_model_data$cover > 0, ],
  #                                         family = Gamma(link = "log"))
  # 
  # summary(dendrogyra_gam_presence_binom)
  # summary(dendrogyra_gam_abundance_gamma)
  # AIC(dendrogyra_gam_presence_binom)
  # AIC(dendrogyra_gam_abundance_gamma)
  # 
  # draw(dendrogyra_gam_presence_binom)
  # draw(dendrogyra_gam_abundance_gamma)
  # 
  # # Check if any smooths are hitting k limits
  # gam.check(dendrogyra_gam_presence_binom)
  # gam.check(dendrogyra_gam_abundance_gamma)
  # 
  # # Look at concurvity
  # concurvity(dendrogyra_gam_presence_binom, full = TRUE)
  # concurvity(dendrogyra_gam_abundance_gamma, full = TRUE)
  # 
  ################################## DICHOCOENIA ##################################
  
  dichocoenia_model_data = spp_data %>%
    filter(grepl("Dichocoenia", spp))
  dichocoenia_model_data = add_env_variables(dichocoenia_model_data, variables_at_PSUs)
    
  # Two-part model with complexity
  dichocoenia_model_data$present <- ifelse(dichocoenia_model_data$cover > 0, 1, 0)
  
  # tic()
  # dichocoenia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = dichocoenia_model_data,
  #                             select = TRUE,
  #                             family = binomial())
  # toc()
  dichocoenia_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                          s(mean_SST) +
                                          s(dist_to_land),
                                        data = dichocoenia_model_data,
                                        select = TRUE,
                                        family = binomial())
  
  # dichocoenia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = dichocoenia_model_data[dichocoenia_model_data$cover > 0, ],
  #                               select = TRUE,
  #                              family = Gamma(link = "log"))
  dichocoenia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy, k = 12) +
                                           s(mean_PAR) +
                                           s(dist_to_deep),
                                         data = dichocoenia_model_data[dichocoenia_model_data$cover > 0, ],
                                         # select = TRUE,
                                         family = Gamma(link = "log"))
  # dichocoenia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) +
  #                                          s(planform_curv) +
  #                                          s(mean_PAR) + s(dist_to_land),
  #                                        data = dichocoenia_model_data[dichocoenia_model_data$cover > 0, ],
  #                                        select = TRUE,
  #                                        family = Gamma(link = "log"))
  
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
  dichocoenia_roc_curve <- roc(dichocoenia_gam_presence_binom$model$present, 
                            fitted(dichocoenia_gam_presence_binom))
  auc(dichocoenia_roc_curve)
  plot(dichocoenia_roc_curve)
  
  # # Save models
  # saveRDS(dichocoenia_gam_presence_binom, 
  #         here("output", "output_GAMs", "dichocoenia_gam_presence_binom.rds"))
  # 
  # saveRDS(dichocoenia_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "dichocoenia_gam_abundance_gamma.rds"))
  
  ################################## DIPLORIA ##################################
  
  # NOTE - some k issues with the presence model
  
  diploria_model_data = spp_data %>%
    filter(grepl("Diploria", spp))
  diploria_model_data = add_env_variables(diploria_model_data, variables_at_PSUs)
    
  # Two-part model with complexity
  #
  #depth_bathy, aspect, slope, complexity, TPI, VRM, planform_curv, SAPA, max_Hsig,
  #   dir_at_max_Hsig, mean_Hsig, mean_SST, range_SST, mean_PAR, mean_chla, mean_kd490, mean_spm,
  #   dist_to_land, dist_to_deep, max_BOV, year
  diploria_model_data$present <- ifelse(diploria_model_data$cover > 0, 1, 0)
  
  # tic()
  # diploria_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = diploria_model_data,
  #                             select = TRUE,
  #                             family = binomial())
  # toc()
  diploria_gam_presence_binom <- gam(present ~ s(depth_bathy) +
                                       slope +
                                       TPI + s(VRM) +
                                       s(dir_at_max_hsig, k = 12, bs = 'cc') +
                                       s(mean_SST, k = 12) + s(dist_to_deep) + s(max_BOV),
                                     data = diploria_model_data,
                                     select = TRUE,
                                     family = binomial())
  
  # diploria_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = diploria_model_data[diploria_model_data$cover > 0, ],
  #                               select = TRUE,
  #                               family = Gamma(link = "log"))
  diploria_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                        s(SAPA) +
                                        s(mean_chla) +
                                        s(dist_to_deep),
                                      data = diploria_model_data[diploria_model_data$cover > 0, ],
                                      select = TRUE,
                                      family = Gamma(link = "log"))
  
  summary(diploria_gam_presence_binom)
  summary(diploria_gam_abundance_gamma)
  AIC(diploria_gam_presence_binom)
  AIC(diploria_gam_abundance_gamma)
  
  draw(diploria_gam_presence_binom)
  draw(diploria_gam_abundance_gamma)
  
  # Check if any smooths are hitting k limits
  gam.check(diploria_gam_presence_binom)
  gam.check(diploria_gam_abundance_gamma)
  
  # Look at concurvity
  concurvity(diploria_gam_presence_binom, full = TRUE)
  concurvity(diploria_gam_abundance_gamma, full = TRUE)
  
  #AUC / ROC
  diploria_roc_curve <- roc(diploria_gam_presence_binom$model$present, 
                                  fitted(diploria_gam_presence_binom))
  auc(diploria_roc_curve)
  plot(diploria_roc_curve)
  
  # # Save models
  # saveRDS(diploria_gam_presence_binom, 
  #         here("output", "output_GAMs", "diploria_gam_presence_binom.rds"))
  # 
  # saveRDS(diploria_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "diploria_gam_abundance_gamma.rds"))
  
  
  # ################################## EUSMILIA ##################################
  # 
  # # NOTE - abundance model seems very unstable due to sample size. likely need to lump with 'RARE HS'
  # 
  # eusmilia <- spp_data %>%
  #   filter(depth_bathy >= -60) %>%
  #   filter(grepl("Eusmilia", spp))
  # 
  # eusmilia_model_data_filtered <- eusmilia
  # 
  # # Extract environmental data
  # eusmilia_species_env_complex <- extract_env_data(eusmilia)
  # 
  # # Add environmental variables
  # eusmilia <- add_env_variables(eusmilia, eusmilia_species_env_complex, eusmilia_model_data_filtered)
  # 
  # eusmilia_model_data <- eusmilia
  # 
  # # Check how much data you have left
  # cat("eusmilia observations with all complexity variables:", nrow(eusmilia_model_data), "\n")
  # 
  # # # Fit expanded GAM models
  # # eusmilia_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
  # #                                    s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
  # #                                    s(mean_Hsig) + s(mean_PAR) + s(mean_chla) + s(dist_to_land) +
  # #                                    s(dist_to_deep) + s(max_BOV),
  # #                                  data = eusmilia_model_data, family = tw())
  # # 
  # # # Check the results
  # # summary(eusmilia_gam_all_tweedie)
  # # AIC(eusmilia_gam_all_tweedie)
  # # gam.check(eusmilia_gam_all_tweedie)
  # # plot(eusmilia_gam_all_tweedie, pages = 2)
  # # draw(eusmilia_gam_all_tweedie)
  # 
  # # Create correlation matrix
  # eusmilia_cor_matrix <- create_correlation_matrix(eusmilia_model_data)
  # print(round(eusmilia_cor_matrix, 2))
  # 
  # # Two-part model with complexity
  # #
  # #depth_bathy, aspect, slope, complexity, TPI, VRM, planform_curv, SAPA, max_Hsig,
  # #   dir_at_max_Hsig, mean_Hsig, mean_SST, range_SST, mean_PAR, mean_chla, mean_kd490, mean_spm,
  # #   dist_to_land, dist_to_deep, max_BOV, year
  # eusmilia_model_data$present <- ifelse(eusmilia_model_data$cover > 0, 1, 0)
  # 
  # # tic()
  # # eusmilia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  # #                               s(slope) +
  # #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  # #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  # #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  # #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  # #                               s(range_SST) +
  # #                               s(dist_to_land),
  # #                             data = eusmilia_model_data,
  # #                             select = TRUE,
  # #                             family = binomial())
  # # toc()
  # eusmilia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                      s(complexity) +
  #                                      s(mean_kd490) +
  #                                      s(dist_to_deep) +
  #                                      range_SST,
  #                                    data = eusmilia_model_data,
  #                                    select = TRUE,
  #                                    family = binomial())
  # 
  # # eusmilia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  # #                               s(slope) +
  # #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  # #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  # #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  # #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  # #                               s(range_SST) +
  # #                               s(dist_to_land),
  # #                               data = eusmilia_model_data[eusmilia_model_data$cover > 0, ],
  # #                               select = TRUE,
  # #                               family = Gamma(link = "log"))
  # eusmilia_gam_abundance_gamma <- gam(cover ~ depth_bathy +
  #                                       s(max_Hsig),
  #                                     data = eusmilia_model_data[eusmilia_model_data$cover > 0, ],
  #                                     # select = TRUE,
  #                                     family = Gamma(link = "log"))
  # 
  # summary(eusmilia_gam_presence_binom)
  # summary(eusmilia_gam_abundance_gamma)
  # AIC(eusmilia_gam_presence_binom)
  # AIC(eusmilia_gam_abundance_gamma)
  # 
  # draw(eusmilia_gam_presence_binom)
  # draw(eusmilia_gam_abundance_gamma)
  # 
  # # Check if any smooths are hitting k limits
  # gam.check(eusmilia_gam_presence_binom)
  # gam.check(eusmilia_gam_abundance_gamma)
  # 
  # # Look at concurvity
  # concurvity(eusmilia_gam_presence_binom, full = TRUE)
  # concurvity(eusmilia_gam_abundance_gamma, full = TRUE)
  # 
  # #AUC / ROC
  # eusmilia_fitted_data <- eusmilia_gam_presence_binom$model
  # eusmilia_roc_curve <- roc(eusmilia_fitted_data$present, 
  #                           fitted(eusmilia_gam_presence_binom))
  # auc(eusmilia_roc_curve)
  # plot(eusmilia_roc_curve, main = "ROC Curve for Madracis Presence Model")
  # 
  # # Save models
  # saveRDS(eusmilia_gam_presence_binom, 
  #         here("output", "output_GAMs", "eusmilia_gam_presence_binom.rds"))
  # 
  # saveRDS(eusmilia_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "eusmilia_gam_abundance_gamma.rds"))
  # 
  # 
  # ################################## MEANDRINA ##################################
  # 
  # # NOTE - abundance model seems very unstable due to sample size. likely need to lump with 'RARE HS'
  # 
  # meandrina <- spp_data %>%
  #   filter(depth_bathy >= -60) %>%
  #   filter(grepl("Meandrina", spp))
  # 
  # meandrina_model_data_filtered <- meandrina
  # 
  # # Extract environmental data
  # meandrina_species_env_complex <- extract_env_data(meandrina)
  # 
  # # Add environmental variables
  # meandrina <- add_env_variables(meandrina, meandrina_species_env_complex, meandrina_model_data_filtered)
  # 
  # meandrina_model_data <- meandrina
  # 
  # # Check how much data you have left
  # cat("meandrina observations with all complexity variables:", nrow(meandrina_model_data), "\n")
  # 
  # # # Fit expanded GAM models
  # # meandrina_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
  # #                                    s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
  # #                                    s(mean_Hsig) + s(mean_PAR) + s(mean_chla) + s(dist_to_land) +
  # #                                    s(dist_to_deep) + s(max_BOV),
  # #                                  data = meandrina_model_data, family = tw())
  # # 
  # # # Check the results
  # # summary(meandrina_gam_all_tweedie)
  # # AIC(meandrina_gam_all_tweedie)
  # # gam.check(meandrina_gam_all_tweedie)
  # # plot(meandrina_gam_all_tweedie, pages = 2)
  # # draw(meandrina_gam_all_tweedie)
  # 
  # # Create correlation matrix
  # meandrina_cor_matrix <- create_correlation_matrix(meandrina_model_data)
  # print(round(meandrina_cor_matrix, 2))
  # 
  # # Two-part model with complexity
  # #
  # #depth_bathy, aspect, slope, complexity, TPI, VRM, planform_curv, SAPA, max_Hsig,
  # #   dir_at_max_Hsig, mean_Hsig, mean_SST, range_SST, mean_PAR, mean_chla, mean_kd490, mean_spm,
  # #   dist_to_land, dist_to_deep, max_BOV, year
  # meandrina_model_data$present <- ifelse(meandrina_model_data$cover > 0, 1, 0)
  # 
  # # tic()
  # # meandrina_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  # #                               s(slope) +
  # #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  # #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  # #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  # #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  # #                               s(range_SST) +
  # #                               s(dist_to_land),
  # #                             data = meandrina_model_data,
  # #                             select = TRUE,
  # #                             family = binomial())
  # # toc()
  # 
  # #maybe drop Hsig
  # meandrina_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                       s(VRM) +
  #                                       s(dir_at_max_hsig, bs = 'cc') +
  #                                       s(mean_SST) + s(mean_PAR) + s(mean_kd490) +
  #                                       s(range_SST),
  #                                     data = meandrina_model_data,
  #                                     select = TRUE,
  #                                     family = binomial())
  # 
  # # meandrina_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  # #                               s(slope) +
  # #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  # #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  # #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  # #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  # #                               s(range_SST) +
  # #                               s(dist_to_land),
  # #                               data = meandrina_model_data[meandrina_model_data$cover > 0, ],
  # #                               select = TRUE,
  # #                               family = Gamma(link = "log"))
  # 
  # #meanSST, rangePAR, meanchla, diratmax, maxBOV, disttodeep?
  # meandrina_gam_abundance_gamma <- gam(cover ~ s(TPI) + s(dir_at_max_hsig, bs = 'cc') +
  #                                        mean_chla + s(mean_PAR),
  #                                      data = meandrina_model_data[meandrina_model_data$cover > 0, ],
  #                                      select = TRUE,
  #                                      family = Gamma(link = "log"))
  # 
  # summary(meandrina_gam_presence_binom)
  # summary(meandrina_gam_abundance_gamma)
  # AIC(meandrina_gam_presence_binom)
  # AIC(meandrina_gam_abundance_gamma)
  # 
  # draw(meandrina_gam_presence_binom)
  # draw(meandrina_gam_abundance_gamma)
  # 
  # # Check if any smooths are hitting k limits
  # gam.check(meandrina_gam_presence_binom)
  # gam.check(meandrina_gam_abundance_gamma)
  # 
  # # Look at concurvity
  # concurvity(meandrina_gam_presence_binom, full = TRUE)
  # concurvity(meandrina_gam_abundance_gamma, full = TRUE)
  # 
  # #AUC / ROC
  # meandrina_fitted_data <- meandrina_gam_presence_binom$model
  # meandrina_roc_curve <- roc(meandrina_fitted_data$present, 
  #                           fitted(meandrina_gam_presence_binom))
  # auc(meandrina_roc_curve)
  # plot(meandrina_roc_curve, main = "ROC Curve for Madracis Presence Model")
  # 
  # # Save models
  # saveRDS(meandrina_gam_presence_binom, 
  #         here("output", "output_GAMs", "meandrina_gam_presence_binom.rds"))
  # 
  # saveRDS(meandrina_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "meandrina_gam_abundance_gamma.rds"))
  # 
  # 
  # 
  # ################################## MYCETOPHYLLIA ##################################
  # 
  # # NOTE - presence model might be unstable due to sample size. had to drop meanHsig because of k issues
  # #           - yeah. the abundance model is definitely unstable too. crazy overfitted. should
  # #               pool with the Rare HS group
  # 
  # mycetophyllia <- spp_data %>%
  #   filter(depth_bathy >= -60) %>%
  #   filter(grepl("Mycetophyllia", spp))
  # 
  # mycetophyllia_model_data_filtered <- mycetophyllia
  # 
  # # Extract environmental data
  # mycetophyllia_species_env_complex <- extract_env_data(mycetophyllia)
  # 
  # # Add environmental variables
  # mycetophyllia <- add_env_variables(mycetophyllia, mycetophyllia_species_env_complex, mycetophyllia_model_data_filtered)
  # 
  # mycetophyllia_model_data <- mycetophyllia
  # 
  # # Check how much data you have left
  # cat("mycetophyllia observations with all complexity variables:", nrow(mycetophyllia_model_data), "\n")
  # 
  # # # Fit expanded GAM models
  # # mycetophyllia_gam_all_tweedie <- gam(cover ~ s(depth_bathy) + s(TPI) + s(complexity) +
  # #                                    s(range_SST) + s(mean_SST) + s(dir_at_max_hsig, bs = 'cc') +
  # #                                    s(mean_Hsig) + s(mean_PAR) + s(mean_chla) + s(dist_to_land) +
  # #                                    s(dist_to_deep) + s(max_BOV),
  # #                                  data = mycetophyllia_model_data, family = tw())
  # # 
  # # # Check the results
  # # summary(mycetophyllia_gam_all_tweedie)
  # # AIC(mycetophyllia_gam_all_tweedie)
  # # gam.check(mycetophyllia_gam_all_tweedie)
  # # plot(mycetophyllia_gam_all_tweedie, pages = 2)
  # # draw(mycetophyllia_gam_all_tweedie)
  # 
  # # Create correlation matrix
  # mycetophyllia_cor_matrix <- create_correlation_matrix(mycetophyllia_model_data)
  # print(round(mycetophyllia_cor_matrix, 2))
  # 
  # # Two-part model with complexity
  # #
  # #depth_bathy, aspect, slope, complexity, TPI, VRM, planform_curv, SAPA, max_Hsig,
  # #   dir_at_max_Hsig, mean_Hsig, mean_SST, range_SST, mean_PAR, mean_chla, mean_kd490, mean_spm,
  # #   dist_to_land, dist_to_deep, max_BOV, year
  # mycetophyllia_model_data$present <- ifelse(mycetophyllia_model_data$cover > 0, 1, 0)
  # 
  # # tic()
  # # mycetophyllia_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  # #                               s(slope) +
  # #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  # #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  # #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  # #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  # #                               s(range_SST) +
  # #                               s(dist_to_land),
  # #                             data = mycetophyllia_model_data,
  # #                             select = TRUE,
  # #                             family = binomial())
  # # toc()
  # mycetophyllia_gam_presence_binom <- gam(present ~ s(depth_bathy) +
  #                                           s(VRM) +
  #                                           s(dir_at_max_hsig, bs = 'cc'),
  #                                         data = mycetophyllia_model_data,
  #                                         select = TRUE,
  #                                         family = binomial())
  # 
  # # mycetophyllia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  # #                               s(slope) +
  # #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  # #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  # #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  # #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  # #                               s(range_SST) +
  # #                               s(dist_to_land),
  # #                               data = mycetophyllia_model_data[mycetophyllia_model_data$cover > 0, ],
  # #                               select = TRUE,
  # #                               family = Gamma(link = "log"))
  # mycetophyllia_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) +
  #                                            s(TPI),
  #                                          data = mycetophyllia_model_data[mycetophyllia_model_data$cover > 0, ],
  #                                          # select = TRUE,
  #                                          family = Gamma(link = "log"))
  # 
  # summary(mycetophyllia_gam_presence_binom)
  # summary(mycetophyllia_gam_abundance_gamma)
  # AIC(mycetophyllia_gam_presence_binom)
  # AIC(mycetophyllia_gam_abundance_gamma)
  # 
  # draw(mycetophyllia_gam_presence_binom)
  # draw(mycetophyllia_gam_abundance_gamma)
  # 
  # # Check if any smooths are hitting k limits
  # gam.check(mycetophyllia_gam_presence_binom)
  # gam.check(mycetophyllia_gam_abundance_gamma)
  # 
  # # Look at concurvity
  # concurvity(mycetophyllia_gam_presence_binom, full = TRUE)
  # concurvity(mycetophyllia_gam_abundance_gamma, full = TRUE)
  # 
  # #AUC / ROC
  # mycetophyllia_fitted_data <- mycetophyllia_gam_presence_binom$model
  # mycetophyllia_roc_curve <- roc(mycetophyllia_fitted_data$present, 
  #                           fitted(mycetophyllia_gam_presence_binom))
  # auc(mycetophyllia_roc_curve)
  # plot(mycetophyllia_roc_curve, main = "ROC Curve for Madracis Presence Model")
  # 
  # # Save models
  # saveRDS(mycetophyllia_gam_presence_binom, 
  #         here("output", "output_GAMs", "mycetophyllia_gam_presence_binom.rds"))
  # 
  # saveRDS(mycetophyllia_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "mycetophyllia_gam_abundance_gamma.rds"))
  # 
  # 
  ################################## PSEUDODIPLORIA ##################################
  
  
  
  # # random forest test stuff:
  # library(randomForest)
  # 
  # 
  # 
  # pseudodiploria_rf_presence <- randomForest(as.factor(present) ~ depth_bathy + aspect +
  #                                              slope + complexity + TPI + VRM + planform_curv + SAPA +
  #                                              max_Hsig + dir_at_max_hsig + mean_Hsig +
  #                                              mean_SST + mean_PAR + mean_chla + mean_kd490 +
  #                                              mean_spm + dist_to_deep + max_BOV +
  #                                              range_SST + dist_to_land,
  #                                            data = pseudodiploria_model_data,
  #                                            ntree = 500,
  #                                            importance = TRUE,
  #                                            na.action = na.omit,
  #                                            sampsize = c("0" = 300, "1" = 300))
  # 
  # pseudodiploria_rf_presence <- randomForest(as.factor(present) ~ mean_SST + max_BOV + 
  #                                              depth_bathy + mean_Hsig + max_Hsig + 
  #                                              range_SST + dist_to_deep,  # Keep top 7
  #                                            data = pseudodiploria_model_data,
  #                                            ntree = 500,
  #                                            importance = TRUE,
  #                                            # sampsize = c("0" = 200, "1" = 200),
  #                                            na.action = na.omit)
  # 
  # # Basic model assessment
  # print(pseudodiploria_rf_presence)
  # # Out-of-bag error rate and confusion matrix
  # pseudodiploria_rf_presence$confusion
  # # Variable importance
  # importance(pseudodiploria_rf_presence)
  # varImpPlot(pseudodiploria_rf_presence)
  # 
  # # Predictions and confusion matrix
  # rf_pred <- predict(pseudodiploria_rf_presence, type = "response")
  # 
  # # Dynamic variable extraction from the model
  # model_vars <- all.vars(pseudodiploria_rf_presence$terms)
  # present_complete <- pseudodiploria_model_data$present[complete.cases(pseudodiploria_model_data[model_vars])]
  # 
  # table(present_complete, rf_pred)
  # # Classification accuracy
  # mean(rf_pred == as.factor(present_complete))
  # # ROC curve and AUC
  # library(pROC)
  # rf_prob <- predict(pseudodiploria_rf_presence, type = "prob")[,2]
  # roc_rf <- roc(present_complete, rf_prob)
  # auc(roc_rf)
  # plot(roc_rf)
  # # Confusion matrix plot
  # library(ggplot2)
  # confusion_df <- data.frame(
  #   Observed = factor(rep(c("Absent", "Present"), each = 2)),
  #   Predicted = factor(rep(c("Absent", "Present"), 2)),
  #   Count = as.vector(table(present_complete, rf_pred))
  # )
  # ggplot(confusion_df, aes(x = Predicted, y = Observed, fill = Count)) +
  #   geom_tile() +
  #   geom_text(aes(label = Count), size = 5) +
  #   scale_fill_gradient(low = "white", high = "blue") +
  #   theme_minimal() +
  #   labs(title = "Confusion Matrix: RF Presence/Absence")
  # # Simple accuracy by site
  # plot(as.numeric(present_complete), as.numeric(rf_pred), 
  #      xlab = "Observed (1=Absent, 2=Present)", ylab = "Predicted (1=Absent, 2=Present)",
  #      main = "Observed vs Predicted Presence")
  # abline(0, 1, col = "red")
  
  
  
  
  
  # NOTE - should consider removing SST for presence because of k issue
  #         - dropped max_Hsig, though seemingly important, because of high concurvity in abundance model
  #         - and dropped mean_PAR from abundance model, because of k issue
  
  pseudodiploria_model_data = spp_data %>%
    filter(grepl("Pseudodiploria", spp))
  pseudodiploria_model_data = add_env_variables(pseudodiploria_model_data, variables_at_PSUs)
  
  # Two-part model with complexity
  pseudodiploria_model_data$present <- ifelse(pseudodiploria_model_data$cover > 0, 1, 0)
  
  
  
  n_absent <- sum(pseudodiploria_model_data$present == 0)
  n_present <- sum(pseudodiploria_model_data$present == 1)
  balance_ratio <- n_absent / n_present
  
  weights_vec <- ifelse(pseudodiploria_model_data$present == 1, balance_ratio, 1)
  
  # tic()
  # pseudodiploria_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                             data = pseudodiploria_model_data,
  #                             weights = weights_vec,
  #                             # select = TRUE,
  #                             family = binomial())
  # toc()
  # # non-weighted version; whittled down with worst concurvity
  # pseudodiploria_gam_presence_binom <- gam(present ~ s(depth_bathy) +
  #                                            slope +
  #                                            s(VRM) +
  #                                            s(mean_Hsig) +
  #                                            s(mean_chla),
  #                                          data = pseudodiploria_model_data,
  #                                          # weights = weights_vec,
  #                                          select = TRUE,
  #                                          family = binomial())
  # # non-weighted version; whittled down with observed/estimate concurvity
  # pseudodiploria_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                            s(complexity) + s(planform_curv) +
  #                                            s(dir_at_max_hsig, bs = 'cc') +
  #                                            s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                                            s(max_BOV) +
  #                                            s(range_SST) +
  #                                            s(dist_to_land),
  #                                          data = pseudodiploria_model_data,
  #                                          # weights = weights_vec,
  #                                          # select = TRUE,
  #                                          family = binomial())
  # # weighted version; whittled down with observed/estimate concurvity
  # #   NOTE - this seems slightly better. runs longer...but nice because it simplifies
  # #           selecting a threshold (or, seems to)
  # #           - actually, I might have just been accidentally running this on the orbicella
  # #               data the whole time, so I need to reassess lol. shit
  # pseudodiploria_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                                            s(planform_curv) + s(SAPA) +
  #                                            s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') +
  #                                            s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                                            s(range_SST) +
  #                                            s(dist_to_land),
  #                                          data = pseudodiploria_model_data,
  #                                          weights = weights_vec,
  #                                          # select = TRUE,
  #                                          family = binomial())
  # weighted version re-done; whittled down with observed/estimate concurvity
  #dirmaxhig, bov, dist to land? hard to know what to do with these. also, SAPA.
  #   NOTE - should really consider dropping distance to land and/or deep here - creates some very strange
  #             banding around places like PR where there is a short distance between shoreline and dropoff
  pseudodiploria_gam_presence_binom <- gam(present ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                             s(planform_curv) +
                                             s(dir_at_max_hsig, bs = 'cc') +
                                             s(mean_SST) + s(mean_PAR) + s(mean_chla) +
                                             s(dist_to_deep) + s(max_BOV) +
                                             s(dist_to_land),
                                           data = pseudodiploria_model_data,
                                           weights = weights_vec,
                                           # select = TRUE,
                                           family = binomial())
  
  
  # Use the model's data directly (avoids length mismatch)
  gam_pred_prob <- predict(pseudodiploria_gam_presence_binom, type = "response")
  gam_pred_binary <- ifelse(gam_pred_prob > 0.5, 1, 0)
  
  # Get fitted values length
  observed_fitted <- pseudodiploria_gam_presence_binom$y  # This matches the prediction length
  table(observed_fitted, gam_pred_binary)
  
  pseudodiploria_model_data$cover_prop <- pseudodiploria_model_data$cover / 100
  
  # pseudodiploria_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = pseudodiploria_model_data[pseudodiploria_model_data$cover > 0, ],
  #                               select = TRUE,
  #                               family = Gamma(link = "log"))
  # pseudodiploria_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) +
  #                                             s(VRM) +
  #                                             mean_Hsig +
  #                                             s(mean_SST),
  #                                           data = pseudodiploria_model_data[pseudodiploria_model_data$cover > 0, ],
  #                                           select = TRUE,
  #                                           family = Gamma(link = "log"))
  #version with gamma (just log link to start), that uses observed/estimate concurvity
  #   note - inverse does seem to do better...maybe start with it instead
  pseudodiploria_gam_abundance_gamma <- gam(cover ~ s(depth_bathy) + s(aspect, bs = 'cc') +
                                              s(VRM) + s(planform_curv) +
                                              s(mean_Hsig) +
                                              s(mean_SST) + s(mean_PAR) + s(mean_kd490) +
                                              s(dist_to_deep) +
                                              s(dist_to_land),
                                            data = pseudodiploria_model_data[pseudodiploria_model_data$cover > 0, ],
                                            # select = TRUE,
                                            family = Gamma(link = "inverse"))
  
  
  pred <- predict(pseudodiploria_gam_abundance_gamma, type = "response")
  obs <- pseudodiploria_gam_abundance_gamma$model$cover
  lim <- c(0, max(c(obs, pred)))
  plot(obs, pred, xlim = lim, ylim = lim)
  # plot(obs, pred, xlim = 40, ylim = 40)
  abline(0, 1, col = "red")
  r <- cor(obs, pred)
  text(0.05*max(lim), 0.95*max(lim), paste("R² =", round(r^2, 3), "\nR =", round(r, 3), ifelse(cor.test(obs, pred)$p.value < 0.05, "*", "ns")), adj = c(0, 1))
  
  # pseudodiploria_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) + s(aspect, bs = 'cc') +
  #                               s(slope) +
  #                               s(complexity) + s(TPI) + s(VRM) + s(planform_curv) + s(SAPA) +
  #                               s(max_Hsig) + s(dir_at_max_hsig, bs = 'cc') + s(mean_Hsig) +
  #                               s(mean_SST) + s(mean_PAR) + s(mean_chla) + s(mean_kd490) +
  #                               s(mean_spm) + s(dist_to_deep) + s(max_BOV) +
  #                               s(range_SST) +
  #                               s(dist_to_land),
  #                               data = pseudodiploria_model_data[pseudodiploria_model_data$cover_prop > 0, ],
  #                               select = TRUE,
  #                               family = betar())
  #beta using observed/estimate concurvity
  # NOTE - could consider adding spatial smooth and/or year to models like this that struggle a bit
  pseudodiploria_gam_abundance_beta <- gam(cover_prop ~ s(depth_bathy) +
                                             s(slope) +
                                             s(mean_Hsig) +
                                             s(mean_SST) + s(mean_PAR) + s(mean_kd490) +
                                             s(dist_to_deep) +
                                             s(range_SST),
                                           data = pseudodiploria_model_data[pseudodiploria_model_data$cover_prop > 0, ],
                                           # select = TRUE,
                                           family = betar())
  
  pred <- predict(pseudodiploria_gam_abundance_beta, type = "response")
  obs <- pseudodiploria_gam_abundance_beta$model$cover
  lim <- c(0, max(c(obs, pred)))
  plot(obs, pred, xlim = lim, ylim = lim)
  # plot(obs, pred, xlim = 40, ylim = 40)
  abline(0, 1, col = "red")
  r <- cor(obs, pred)
  text(0.05*max(lim), 0.95*max(lim), paste("R² =", round(r^2, 3), "\nR =", round(r, 3), ifelse(cor.test(obs, pred)$p.value < 0.05, "*", "ns")), adj = c(0, 1))
  
  
  summary(pseudodiploria_gam_presence_binom)
  summary(pseudodiploria_gam_abundance_gamma)
  summary(pseudodiploria_gam_abundance_beta)
  AIC(pseudodiploria_gam_presence_binom)
  AIC(pseudodiploria_gam_abundance_gamma)
  AIC(pseudodiploria_gam_abundance_beta)
  
  draw(pseudodiploria_gam_presence_binom)
  draw(pseudodiploria_gam_abundance_gamma)
  draw(pseudodiploria_gam_abundance_beta)
  
  # Check if any smooths are hitting k limits
  gam.check(pseudodiploria_gam_presence_binom)
  gam.check(pseudodiploria_gam_abundance_gamma)
  gam.check(pseudodiploria_gam_abundance_beta)
  
  # Look at concurvity
  concurvity(pseudodiploria_gam_presence_binom, full = TRUE)
  concurvity(pseudodiploria_gam_abundance_gamma, full = TRUE)
  concurvity(pseudodiploria_gam_abundance_beta, full = TRUE)
  
  #AUC / ROC
  pseudodiploria_roc_curve <- roc(pseudodiploria_gam_presence_binom$model$present, 
                             fitted(pseudodiploria_gam_presence_binom))
  auc(pseudodiploria_roc_curve)
  plot(pseudodiploria_roc_curve)
  
  # # Save models
  # saveRDS(pseudodiploria_gam_presence_binom, 
  #         here("output", "output_GAMs", "pseudodiploria_gam_presence_binom.rds"))
  # 
  # saveRDS(pseudodiploria_gam_abundance_gamma, 
  #         here("output", "output_GAMs", "pseudodiploria_gam_abundance_gamma.rds"))
  
  
  
  
  ################################## RARE HS ##################################
  
  
  #Dendrogyra cylindrus, Eusmilia fastigiata, Favia fragum, Isophyllia spp., Manicina areolata,
  #   Meandrina spp., Mussa angulosa, Mycetophyllia spp., and Scolymia spp. (9 species)
  
  
  
  
  # stopping point - 11 sep 2025
  #
  #   - main goal now is actually seeing how long a full SDM output will take to chug with the hurdle 
  #       model approach. then, running all the 80/20 validations to create R-squared, MAE, RMSE, etc.
  #
  #   - may need to handle some predictor outliers...but for now, table that until I see whether the
  #       predictions are turning out well. along the same lines, can later consider mean BOV,
  #       mean wavedir, and mean PAR if needed
  #   - also will go back upstream to group the rare HS corals....but can wait on that for now
  
  
  # ################################## test prediction of Orbicella ##################################
  # 
  # # Manually specify required variables (using raster layer names)
  # presence_vars <- c("depth", "complexity", "range_SST", "dir_at_max_hsig", "max_BOV")
  # abundance_vars <- c("depth", "range_SST", "mean_chla")
  # required_vars <- unique(c(presence_vars, abundance_vars))
  # 
  # # Subset raster stack to only required variables
  # env_subset <- env_complex[[required_vars]]
  # 
  # # Convert subset raster to dataframe
  # env_df <- as.data.frame(env_subset, xy = TRUE, na.rm = TRUE)
  # 
  # # Rename depth column to match model expectation
  # names(env_df)[names(env_df) == "depth"] <- "depth_bathy"
  # 
  # # Take only 1/100th of the data for testing
  # sample_size <- ceiling(nrow(env_df) / 100)
  # set.seed(123)  # For reproducible sampling
  # sample_indices <- sample(nrow(env_df), sample_size)
  # env_df_sample <- env_df[sample_indices, ]
  # 
  # cat("Full dataset size:", nrow(env_df), "cells\n")
  # cat("Sample size (1/100th):", nrow(env_df_sample), "cells\n")
  # 
  # # Test predictions on sample
  # cat("Testing presence prediction...\n")
  # start_time <- Sys.time()
  # presence_prob_sample <- predict(orbicella_gam_presence_binom, newdata = env_df_sample, type = "response")
  # presence_time <- difftime(Sys.time(), start_time, units = "secs")
  # cat("Sample presence prediction:", round(presence_time, 2), "seconds\n")
  # 
  # cat("Testing abundance prediction...\n")
  # start_time <- Sys.time()
  # abundance_pred_sample <- predict(orbicella_gam_abundance_gamma, newdata = env_df_sample, type = "response")
  # abundance_time <- difftime(Sys.time(), start_time, units = "secs")
  # cat("Sample abundance prediction:", round(abundance_time, 2), "seconds\n")
  # 
  # # Create hurdle predictions
  # hurdle_pred_sample <- presence_prob_sample * abundance_pred_sample
  # 
  # # Add predictions to sample dataframe
  # env_df_sample$presence_prob <- presence_prob_sample
  # env_df_sample$abundance_pred <- abundance_pred_sample
  # env_df_sample$hurdle_pred <- hurdle_pred_sample
  # 
  # # Estimate full dataset time
  # total_sample_time <- as.numeric(presence_time + abundance_time)
  # estimated_full_time <- total_sample_time * 100
  # 
  # cat("\nSample results summary:\n")
  # cat("Presence probability range:", round(range(presence_prob_sample), 4), "\n")
  # cat("Abundance prediction range:", round(range(abundance_pred_sample), 6), "\n")
  # cat("Hurdle prediction range:", round(range(hurdle_pred_sample), 6), "\n")
  # 
  # cat("\nTime estimates:\n")
  # cat("Sample time:", round(total_sample_time, 2), "seconds\n")
  # cat("Estimated full dataset time:", round(estimated_full_time, 1), "seconds")
  # if(estimated_full_time > 60) {
  #   cat(" (", round(estimated_full_time/60, 1), " minutes)", sep="")
  # }
  # cat("\n")
  # 
  # # Quick plot of sample results
  # library(ggplot2)
  # library(viridis)
  # ggplot(env_df_sample, aes(x = x, y = y, color = presence_prob)) +
  #   geom_point(size = 0.5) +
  #   scale_color_viridis_c(name = "Presence\nProbability", limits = c(0.1, 1)) +
  #   coord_equal() +
  #   theme_minimal() +
  #   ggtitle("Orbicella Presence Probability")
  # 
  # ggplot(env_df_sample, aes(x = x, y = y, color = hurdle_pred)) +
  #   geom_point(size = 0.5) +
  #   scale_color_viridis_c(name = "Predicted\nCover", limits = c(0, 15)) +
  #   coord_equal() +
  #   theme_minimal() +
  #   ggtitle(paste("Hurdle Model Sample Prediction (n =", nrow(env_df_sample), ")"))
  # 
  # ggplot(env_df_sample, aes(x = x, y = y, color = abundance_pred)) +
  #   geom_point(size = 0.5) +
  #   scale_color_viridis_c(name = "Predicted\nCover", limits = c(0, 15)) +
  #   coord_equal() +
  #   theme_minimal() +
  #   ggtitle(paste("Raw Abundance Model Sample Prediction (n =", nrow(env_df_sample), ")"))
  # 
  # ################################## test prediction of Agaricia ##################################
  # 
  # # Manually specify required variables for agaricia models
  # presence_vars <- c("depth", "VRM", "aspect", "TPI", "mean_SST", "dir_at_max_hsig", "max_BOV")
  # abundance_vars <- c("depth", "complexity", "dir_at_max_hsig", "mean_chla")
  # required_vars <- unique(c(presence_vars, abundance_vars))
  # 
  # # Subset raster stack to only required variables
  # env_subset <- env_complex[[required_vars]]
  # 
  # # Convert subset raster to dataframe
  # env_df <- as.data.frame(env_subset, xy = TRUE, na.rm = TRUE)
  # 
  # # Rename depth column to match model expectation
  # names(env_df)[names(env_df) == "depth"] <- "depth_bathy"
  # 
  # # Take only 1/100th of the data for testing
  # sample_size <- ceiling(nrow(env_df) / 100)
  # set.seed(123)
  # sample_indices <- sample(nrow(env_df), sample_size)
  # env_df_sample <- env_df[sample_indices, ]
  # 
  # cat("Full dataset size:", nrow(env_df), "cells\n")
  # cat("Sample size (1/100th):", nrow(env_df_sample), "cells\n")
  # 
  # # Test predictions on sample - using AGARICIA models
  # cat("Testing presence prediction...\n")
  # start_time <- Sys.time()
  # presence_prob_sample <- predict(agaricia_gam_presence_binom, newdata = env_df_sample, type = "response")
  # presence_time <- difftime(Sys.time(), start_time, units = "secs")
  # cat("Sample presence prediction:", round(presence_time, 2), "seconds\n")
  # 
  # cat("Testing abundance prediction...\n")
  # start_time <- Sys.time()
  # abundance_pred_sample <- predict(agaricia_gam_abundance_gamma, newdata = env_df_sample, type = "response")
  # abundance_time <- difftime(Sys.time(), start_time, units = "secs")
  # cat("Sample abundance prediction:", round(abundance_time, 2), "seconds\n")
  # 
  # # Create hurdle predictions
  # hurdle_pred_sample <- presence_prob_sample * abundance_pred_sample
  # 
  # # Add predictions to sample dataframe
  # env_df_sample$presence_prob <- presence_prob_sample
  # env_df_sample$abundance_pred <- abundance_pred_sample
  # env_df_sample$hurdle_pred <- hurdle_pred_sample
  # 
  # # Estimate full dataset time
  # total_sample_time <- as.numeric(presence_time + abundance_time)
  # estimated_full_time <- total_sample_time * 100
  # 
  # cat("\nSample results summary:\n")
  # cat("Presence probability range:", round(range(presence_prob_sample), 4), "\n")
  # cat("Abundance prediction range:", round(range(abundance_pred_sample), 6), "\n")
  # cat("Hurdle prediction range:", round(range(hurdle_pred_sample), 6), "\n")
  # 
  # cat("\nTime estimates:\n")
  # cat("Sample time:", round(total_sample_time, 2), "seconds\n")
  # cat("Estimated full dataset time:", round(estimated_full_time, 1), "seconds")
  # if(estimated_full_time > 60) {
  #   cat(" (", round(estimated_full_time/60, 1), " minutes)", sep="")
  # }
  # cat("\n")
  # 
  # # Quick plot of sample results
  # library(ggplot2)
  # library(viridis)
  # ggplot(env_df_sample, aes(x = x, y = y, color = presence_prob)) +
  #   geom_point(size = 0.5) +
  #   scale_color_viridis_c(name = "Presence\nProbability", limits = c(0.1, 1)) +
  #   coord_equal() +
  #   theme_minimal() +
  #   ggtitle("Agaricia Presence Probability")
  # 
  # ggplot(env_df_sample, aes(x = x, y = y, color = hurdle_pred)) +
  #   geom_point(size = 0.5) +
  #   scale_color_viridis_c(name = "Predicted\nCover", limits = c(0, 5)) +  # Adjusted for agaricia range
  #   # scale_color_viridis_c(name = "Predicted\nCover") +  # Adjusted for agaricia range
  #   coord_equal() +
  #   theme_minimal() +
  #   ggtitle(paste("Agaricia Hurdle Model Sample Prediction (n =", nrow(env_df_sample), ")"))
  # 
  # 
  # 
  # ################################## test prediction of Pseudodiploria ##################################
  # 
  # # Manually specify required variables for pseudodiploria models
  # presence_vars <- c("depth", "slope", "VRM", "mean_Hsig", "mean_chla")
  # abundance_vars <- c("depth", "VRM", "mean_Hsig", "mean_SST")
  # required_vars <- unique(c(presence_vars, abundance_vars))
  # 
  # # Subset raster stack to only required variables
  # env_subset <- env_complex[[required_vars]]
  # 
  # # Convert subset raster to dataframe
  # env_df <- as.data.frame(env_subset, xy = TRUE, na.rm = TRUE)
  # 
  # # Rename depth column to match model expectation
  # names(env_df)[names(env_df) == "depth"] <- "depth_bathy"
  # 
  # # Take only 1/100th of the data for testing
  # sample_size <- ceiling(nrow(env_df) / 100)
  # set.seed(123)
  # sample_indices <- sample(nrow(env_df), sample_size)
  # env_df_sample <- env_df[sample_indices, ]
  # 
  # cat("Full dataset size:", nrow(env_df), "cells\n")
  # cat("Sample size (1/100th):", nrow(env_df_sample), "cells\n")
  # 
  # # Test predictions on sample - using PSEUDODIPLORIA models
  # cat("Testing presence prediction...\n")
  # start_time <- Sys.time()
  # presence_prob_sample <- predict(pseudodiploria_gam_presence_binom, newdata = env_df_sample, type = "response")
  # presence_time <- difftime(Sys.time(), start_time, units = "secs")
  # cat("Sample presence prediction:", round(presence_time, 2), "seconds\n")
  # 
  # cat("Testing abundance prediction...\n")
  # start_time <- Sys.time()
  # abundance_pred_sample <- predict(pseudodiploria_gam_abundance_gamma, newdata = env_df_sample, type = "response")
  # abundance_time <- difftime(Sys.time(), start_time, units = "secs")
  # cat("Sample abundance prediction:", round(abundance_time, 2), "seconds\n")
  # 
  # # Create hurdle predictions
  # hurdle_pred_sample <- presence_prob_sample * abundance_pred_sample
  # 
  # # Add predictions to sample dataframe
  # env_df_sample$presence_prob <- presence_prob_sample
  # env_df_sample$abundance_pred <- abundance_pred_sample
  # env_df_sample$hurdle_pred <- hurdle_pred_sample
  # 
  # # Estimate full dataset time
  # total_sample_time <- as.numeric(presence_time + abundance_time)
  # estimated_full_time <- total_sample_time * 100
  # 
  # cat("\nSample results summary:\n")
  # cat("Presence probability range:", round(range(presence_prob_sample), 4), "\n")
  # cat("Abundance prediction range:", round(range(abundance_pred_sample), 6), "\n")
  # cat("Hurdle prediction range:", round(range(hurdle_pred_sample), 6), "\n")
  # 
  # cat("\nTime estimates:\n")
  # cat("Sample time:", round(total_sample_time, 2), "seconds\n")
  # cat("Estimated full dataset time:", round(estimated_full_time, 1), "seconds")
  # if(estimated_full_time > 60) {
  #   cat(" (", round(estimated_full_time/60, 1), " minutes)", sep="")
  # }
  # cat("\n")
  # 
  # # Quick plot of sample results
  # library(ggplot2)
  # library(viridis)
  # # Plot presence probability only
  # ggplot(env_df_sample, aes(x = x, y = y, color = presence_prob)) +
  #   geom_point(size = 0.5) +
  #   scale_color_viridis_c(name = "Presence\nProbability", limits = c(0.02, 1)) +
  #   coord_equal() +
  #   theme_minimal() +
  #   ggtitle("Pseudodiploria Presence Probability")
  # ggplot(env_df_sample, aes(x = x, y = y, color = hurdle_pred)) +
  #   geom_point(size = 0.5) +
  #   scale_color_viridis_c(name = "Predicted\nCover", limits = c(0, 5)) +  # Adjusted for pseudodiploria range
  #   coord_equal() +
  #   theme_minimal() +
  #   ggtitle(paste("Pseudodiploria Hurdle Model Sample Prediction (n =", nrow(env_df_sample), ")"))
  # ggplot(env_df_sample, aes(x = x, y = y, color = abundance_pred)) +
  #   geom_point(size = 0.5) +
  #   scale_color_viridis_c(name = "Predicted\nCover", limits = c(0, 10)) +  # Adjusted for pseudodiploria range
  #   coord_equal() +
  #   theme_minimal() +
  #   ggtitle(paste("Pseudodiploria Raw Abundance Model Sample Prediction (n =", nrow(env_df_sample), ")"))
  # 
  # 
  # 
  
  ################################## hard hurdle test & map ##################################
  
  # STOPPING POINT - 17 sep 2025
  #   - need to add lon and lat as spatrasters to env_complex to be able to make any
  #       predictions using them. still not convinced they'll work anyways...
  
  # Species toggle - change this to switch between species
  # species <- "orbicella" #0.5
  # species <- "pseudodiploria" #0.45
  species <- "orbicella" #0.38
  # species <- "montastraea" #0.25
  # species <- "colpophyllia" #0.1
  # species <- "porites" #0.54
  
  # Dynamically construct model variable names based on species
  presence_model_name <- paste0(species, "_gam_presence_binom")
  abundance_model_name <- paste0(species, "_gam_abundance_beta")
  # abundance_model_name <- paste0(species, "_gam_abundance_gamma")
  
  # Get the actual model objects (assuming they exist in your environment)
  presence_model <- get(presence_model_name)
  abundance_model <- get(abundance_model_name)
  
  # Manually specify required variables (using raster layer names)
  presence_vars <- all.vars(formula(presence_model))[-1]
  abundance_vars <- all.vars(formula(abundance_model))[-1]
  required_vars <- unique(c(presence_vars, abundance_vars))
  
  # Subset raster stack to only required variables
  required_vars <- gsub("^depth_bathy$", "depth", required_vars)
  env_subset <- env_complex[[required_vars]]
  
  # Convert subset raster to dataframe
  env_df <- as.data.frame(env_subset, xy = TRUE, na.rm = TRUE)
  
  # Rename depth column to match model expectation
  names(env_df)[names(env_df) == "depth"] <- "depth_bathy"
  
  # Take only 1/100th of the data for testing
  sample_size <- ceiling(nrow(env_df) / 100)
  set.seed(123)  # For reproducible sampling
  sample_indices <- sample(nrow(env_df), sample_size)
  env_df_sample <- env_df[sample_indices, ]
  
  cat("Species:", species, "\n")
  cat("Full dataset size:", nrow(env_df), "cells\n")
  cat("Sample size (1/100th):", nrow(env_df_sample), "cells\n")
  
  # Test predictions on sample
  cat("Testing presence prediction...\n")
  start_time <- Sys.time()
  presence_prob_sample <- predict(presence_model, newdata = env_df_sample, type = "response")
  presence_time <- difftime(Sys.time(), start_time, units = "secs")
  cat("Sample presence prediction:", round(presence_time, 2), "seconds\n")
  
  # Convert presence probability to binary (hard hurdle)
  presence_threshold <- 0.5 # 0.1 # 0.25 # 0.38 # 0.45 # 0.5
  presence_binary_sample <- ifelse(presence_prob_sample > presence_threshold, 1, 0)
  
  cat("Testing abundance prediction...\n")
  start_time <- Sys.time()
  abundance_pred_sample <- predict(abundance_model, newdata = env_df_sample, type = "response")
  abundance_time <- difftime(Sys.time(), start_time, units = "secs")
  cat("Sample abundance prediction:", round(abundance_time, 2), "seconds\n")
  
  # Create hurdle predictions (only abundance where presence predicted)
  hurdle_pred_sample <- ifelse(presence_binary_sample == 1, abundance_pred_sample, 0)
  
  # Add predictions to sample dataframe
  env_df_sample$presence_prob <- presence_prob_sample
  env_df_sample$presence_binary <- presence_binary_sample
  env_df_sample$abundance_pred <- abundance_pred_sample
  env_df_sample$hurdle_pred <- hurdle_pred_sample
  
  # Estimate full dataset time
  total_sample_time <- as.numeric(presence_time + abundance_time)
  estimated_full_time <- total_sample_time * 100
  
  cat("\nSample results summary:\n")
  cat("Presence probability range:", round(range(presence_prob_sample), 4), "\n")
  cat("Abundance prediction range:", round(range(abundance_pred_sample), 6), "\n")
  cat("Hurdle prediction range:", round(range(hurdle_pred_sample), 6), "\n")
  
  cat("\nTime estimates:\n")
  cat("Sample time:", round(total_sample_time, 2), "seconds\n")
  cat("Estimated full dataset time:", round(estimated_full_time, 1), "seconds")
  if(estimated_full_time > 60) {
    cat(" (", round(estimated_full_time/60, 1), " minutes)", sep="")
  }
  cat("\n")
  
  # Quick plot of sample results
  ggplot(env_df_sample, aes(x = x, y = y, color = presence_prob)) +
    geom_point(size = 0.5) +
    scale_color_viridis_c(name = "Presence\nProbability") +
    coord_equal() +
    theme_minimal() +
    ggtitle(paste(stringr::str_to_title(species), "Presence Probability"))
  
  ggplot(env_df_sample, aes(x = x, y = y, color = hurdle_pred)) +
    geom_point(size = 0.5) +
    scale_color_viridis_c(name = "Predicted\nCover", limits = c(0, 0.1)) +
    coord_equal() +
    theme_minimal() +
    ggtitle(paste("Hurdle Model Sample Prediction -", stringr::str_to_title(species), "(n =", nrow(env_df_sample), ")"))
  
  hist(env_df_sample$hurdle_pred[env_df_sample$hurdle_pred > 0 & env_df_sample$hurdle_pred < 1])
  
  ggplot(env_df_sample, aes(x = x, y = y, color = abundance_pred)) +
    geom_point(size = 0.5) +
    scale_color_viridis_c(name = "Predicted\nCover", limits = c(0, 0.1)) +
    coord_equal() +
    theme_minimal() +
    ggtitle(paste("Raw Abundance Model Sample Prediction -", stringr::str_to_title(species), "(n =", nrow(env_df_sample), ")"))  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Run full model predictions - THIS TAKES A LOT OF TIME / COMPUTING RESOURCES. 15-20+ MINUTES
  cat("\nRunning full model predictions...\n")
  start_time <- Sys.time()
  
  # Full presence prediction
  presence_threshold = 0.5 # 0.38 # 0.5
  presence_prob_full <- predict(presence_model, newdata = env_df, type = "response")
  presence_binary_full <- ifelse(presence_prob_full > presence_threshold, 1, 0)
  
  # Full abundance prediction
  abundance_pred_full <- predict(abundance_model, newdata = env_df, type = "response")
  
  # Full hurdle prediction
  hurdle_pred_full <- ifelse(presence_binary_full == 1, abundance_pred_full, 0)
  
  full_time <- difftime(Sys.time(), start_time, units = "secs")
  cat("Full prediction time:", round(full_time, 1), "seconds\n")
  
  # Add predictions to full dataframe
  env_df$hurdle_pred <- hurdle_pred_full
  
  # # Plot full results
  # ggplot(env_df, aes(x = x, y = y, color = hurdle_pred)) +
  #   geom_point(size = 0.1) +
  #   scale_color_viridis_c(name = "Predicted\nCover", limits = c(0, 0.2)) +
  #   coord_equal() +
  #   theme_minimal() +
  #   ggtitle(paste("Full Hurdle Model Prediction -", stringr::str_to_title(species)))
  
  
  
  # Convert predictions back to raster using the original env_df coordinates
  hurdle_raster <- rast(env_df[, c("x", "y", "hurdle_pred")], type="xyz")
  shallow_mask <- bathy_final > -60
  hurdle_raster <- resample(hurdle_raster, bathy_final)
  hurdle_raster <- mask(hurdle_raster, shallow_mask, maskvalues = FALSE)
  
  # Plot with terra
  hurdle_raster_clamped <- clamp(hurdle_raster, lower = 0, upper = 0.1, values = TRUE)
  # hurdle_raster_clamped <- clamp(hurdle_raster, lower = 0, upper = 1.0, values = TRUE)
  # plot(hurdle_raster_clamped, main = paste("Full Hurdle Model -", stringr::str_to_title(species)), 
  #      colNA = "gray90")  
  
  
  
  
  # Create leaflet map with hurdle predictions
  pal <- colorNumeric("viridis", values(hurdle_raster_clamped), na.color = "transparent")
  
  leaflet() %>%
    addTiles(group = "OpenStreetMap") %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addRasterImage(hurdle_raster_clamped, colors = pal, opacity = 0.7) %>%
    addLegend("bottomright", pal = pal, values = values(hurdle_raster_clamped),
              title = paste(stringr::str_to_title(species), "Cover (%)")) %>%
    addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                     options = layersControlOptions(collapsed = FALSE))
  
  
  
  
  
  
  
  
  #draft code to pull up another leaflet with spp-specific cover
  species <- "orbicella" #0.1
  
  # Get unique PSU locations with coral cover data for the selected species
  spp_leaflet_cover_data <- combined_benthic_data_averaged %>%
    filter(grepl(species, spp, ignore.case = TRUE)) %>%
    group_by(PSU) %>%
    summarise(
      lat = first(lat),
      lon = first(lon),
      dataset = first(dataset),
      avg_cover = mean(cover, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Create color palette for coral cover (yellow to red)
  color_pal <- colorNumeric(
    palette = c("#FFFF00", "#FF8C00", "#FF4500", "#FF0000"),  # Yellow to Red
    domain = spp_leaflet_cover_data$avg_cover,
    na.color = "gray"
  )
  
  # Create the interactive leaflet map
  coral_cover_map <- leaflet(spp_leaflet_cover_data) %>%
    addTiles(group = "OpenStreetMap") %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addCircleMarkers(
      lng = ~lon, 
      lat = ~lat,
      popup = ~paste0(
        "<b>PSU:</b> ", PSU, "<br>",
        "<b>Dataset:</b> ", dataset, "<br>",
        "<b>Average Coral Cover:</b> ", round(avg_cover, 2), "%<br>",
        "<b>Latitude:</b> ", round(lat, 5), "<br>",
        "<b>Longitude:</b> ", round(lon, 5)
      ),
      label = ~paste0("PSU ", PSU, ": ", round(avg_cover, 1), "% cover"),
      radius = 8,
      color = "white",
      fillColor = ~color_pal(avg_cover),
      fillOpacity = 0.8,
      weight = 2,
      stroke = TRUE
    ) %>%
    addLegend(
      "bottomright",
      pal = color_pal,
      values = ~avg_cover,
      title = "Coral Cover (%)",
      opacity = 1,
      labFormat = labelFormat(suffix = "%")
    ) %>%
    addLayersControl(
      baseGroups = c("OpenStreetMap", "Satellite"),
      options = layersControlOptions(collapsed = FALSE)
    ) %>%
    fitBounds(
      lng1 = min(spp_leaflet_cover_data$lon, na.rm = TRUE) - 0.1,
      lat1 = min(spp_leaflet_cover_data$lat, na.rm = TRUE) - 0.1,
      lng2 = max(spp_leaflet_cover_data$lon, na.rm = TRUE) + 0.1,
      lat2 = max(spp_leaflet_cover_data$lat, na.rm = TRUE) + 0.1
    )
  
  # Display the map
  coral_cover_map

  ################################## Save objects/workspace ##################################
  
  # # Get all objects matching the patterns
  # all_objects <- ls(envir = .GlobalEnv)
  # patterns <- c("_gam_presence_binom$", "_gam_abundance_gamma$", "_model_data$")
  # matching_objects <- grep(paste(patterns, collapse = "|"), all_objects, value = TRUE)
  # 
  # # Save each object
  # for (obj_name in matching_objects) {
  #   saveRDS(get(obj_name), here("output", "output_GAMs", paste0(obj_name, ".rds")))
  # }
  
  
  #updated way to handle saving of new objects
  save_new_objects("output/output_GAMs", existing_objects)
  