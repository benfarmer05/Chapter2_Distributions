  
  # .rs.restartR(clean = TRUE)
  
  library(sf)
  library(here)
  library(terra) 
  library(cmocean)

  source(here("src/functions.R"))
  
  ################################## setup ##################################
  
  #load next set
  load_spat_objects(directory = 'output/output_create_habitat_grid/')
  source(here("src/functions.R"))
  load(here('output', 'output_create_habitat_grid/create_habitat_grid_workspace.RData'))
  
  #load just bathy_final if desired (less memory than loading all of the rasters above)
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
  
  existing_objects <- ls(envir = .GlobalEnv)
  
  ################################## Import waves ##################################
  
  # VERSION where SWAN data was provided directly from Miguel Canals
  # 
  SWAN_target_res <- 50
  
  # Load the comprehensive wave summary data
  wave_summary <- read.csv(here("output", "swan_canals_comprehensive_summary_for_R.csv"))
  
  # Clean the data
  # wave_summary <- wave_summary %>%
  #   filter(!is.na(longitude) & !is.na(latitude))
  wave_summary <- wave_summary[!is.na(wave_summary$longitude) & !is.na(wave_summary$latitude), ]
  
  cat(sprintf("Loaded %d SWAN grid points\n", nrow(wave_summary)))
  
  # Create SpatVector from wave data
  wave_points <- vect(wave_summary, 
                      geom = c("longitude", "latitude"), 
                      crs = "EPSG:4326")
  
  # Project to match bathymetry
  wave_points_proj <- project(wave_points, crs(bathy_final))
  
  # Filter wave points to only water areas (non-NA bathymetry)
  cat("Filtering wave points to water areas only...\n")
  
  # Extract bathymetry values at each wave point location
  bathy_at_points <- extract(bathy_final, wave_points_proj)
  
  # Find points where bathymetry is NOT NA (i.e., water areas)
  water_points_idx <- !is.na(bathy_at_points[,2])
  
  cat(sprintf("Points over water: %d/%d (%.1f%%)\n", 
              sum(water_points_idx), length(water_points_idx), 
              sum(water_points_idx)/length(water_points_idx)*100))
  
  # Keep only points that are over water
  wave_points_water <- wave_points_proj[water_points_idx]
  
  cat(sprintf("After filtering to water: %d points\n", nrow(wave_points_water)))
  
  
  create_wave_raster <- function(wave_points, var_name, bathy_template, target_resolution) {
    cat(sprintf("Creating %dm raster for %s...\n", target_resolution, var_name))
    
    # Filter to valid values
    values <- values(wave_points)[[var_name]]
    valid_idx <- !is.na(values) & values >= 0
    
    if (grepl("dir", var_name, ignore.case = TRUE)) {
      valid_idx <- valid_idx & values <= 360
    }
    
    wave_points_valid <- wave_points[valid_idx]
    cat(sprintf("Using %d valid points\n", sum(valid_idx)))
    
    # Create template with target resolution
    bathy_ext <- ext(bathy_template)
    template_raster <- rast(bathy_ext, 
                            resolution = target_resolution,
                            crs = crs(bathy_template))
    
    n_cells <- ncell(template_raster)
    cat(sprintf("Grid: %d x %d cells (%.1f million)\n", 
                ncol(template_raster), nrow(template_raster), n_cells/1e6))
    
    # # Memory safety check
    # if (n_cells > 25e6) {
    #   stop(sprintf("Grid too large (%.1f million cells). Use resolution > 100m.", n_cells/1e6))
    # }
    # 
    # if (n_cells > 10e6) {
    #   warning(sprintf("Large grid (%.1f million cells) - may be slow", n_cells/1e6))
    # }
    
    # Rasterize wave points to template
    result <- rasterize(wave_points_valid, template_raster, 
                        field = var_name, 
                        fun = mean)
    
    # Gap filling with focal operations
    cat("Gap filling...\n")
    
    # For directional data, use circular-aware gap filling
    if (grepl("dir", var_name, ignore.case = TRUE)) {
      # Convert to radians for calculation
      result_rad <- result * pi / 180
      
      # Focal mean of sin and cos components with multiple passes
      sin_comp <- app(result_rad, sin)
      cos_comp <- app(result_rad, cos)
      
      # Iterative gap filling with increasing window sizes
      for (window_size in c(3, 5, 7, 9)) {
        sin_comp <- focal(sin_comp, w = window_size, fun = mean, na.policy = "only", na.rm = TRUE)
        cos_comp <- focal(cos_comp, w = window_size, fun = mean, na.policy = "only", na.rm = TRUE)
      }
      
      # Convert back to degrees
      result <- atan2(sin_comp, cos_comp) * 180 / pi
      result[result < 0] <- result[result < 0] + 360
      
    } else {
      # For non-directional data, use regular mean with multiple passes
      for (window_size in c(3, 5, 7, 9)) {
        result <- focal(result, w = window_size, fun = mean, na.policy = "only", na.rm = TRUE)
      }
    }
    
    # Resample bathymetry to match template and apply mask
    bathy_resampled <- resample(bathy_template, template_raster, method = "near")
    result <- mask(result, bathy_resampled)
    
    names(result) <- var_name
    return(result)
  }
  
  
  cat(sprintf("\n=== Creating all wave rasters at %dm resolution ===\n", SWAN_target_res))
  
  # Create all wave variable rasters
  mean_hsig_raster <- create_wave_raster(wave_points_water, "mean_hmean", bathy_final, SWAN_target_res)
  # max_mean_hsig_raster <- create_wave_raster(wave_points_water, "max_hmean", bathy_final, SWAN_target_res)
  # mean_max_hsig_raster <- create_wave_raster(wave_points_water, "mean_hmax", bathy_final, SWAN_target_res)
  max_hsig_raster <- create_wave_raster(wave_points_water, "max_hmax", bathy_final, SWAN_target_res)
  per_at_max_hsig_raster <- create_wave_raster(wave_points_water, "mean_tp", bathy_final, SWAN_target_res)
  dir_at_max_hsig_raster <- create_wave_raster(wave_points_water, "mean_dir", bathy_final, SWAN_target_res)
  
  # Get extent and properties from bathy_final
  bathy_ext <- ext(bathy_final)
  bathy_crs <- crs(bathy_final)
  bathy_res <- res(bathy_final)
  
  # Create template raster matching bathy_final exactly
  template_raster <- rast(bathy_ext, resolution = bathy_res, crs = bathy_crs)
  
  # Resample wave rasters to match bathy_final
  mean_hsig_raster <- resample(mean_hsig_raster, template_raster, method = "bilinear")
  # max_mean_hsig_raster <- resample(max_mean_hsig_raster, template_raster, method = "bilinear")
  # mean_max_hsig_raster <- resample(mean_max_hsig_raster, template_raster, method = "bilinear")
  max_hsig_raster <- resample(max_hsig_raster, template_raster, method = "bilinear")
  per_at_max_hsig_raster <- resample(per_at_max_hsig_raster, template_raster, method = "bilinear")
  dir_at_max_hsig_raster <- resample(dir_at_max_hsig_raster, template_raster, method = "bilinear")
  
  ################################## calculate BOV ##################################
  
  # # VERSION with no chunking; easier to diagnose issues
  # # Quick BOV Test on Small Subset (1000 cells south of St Thomas)
  # # Constants
  # g <- 9.81
  # 
  # # Make depth positive
  # depth_raster <- abs(bathy_final)
  # 
  # cat("Creating test subset south of St Thomas...\n")
  # 
  # # Define a bounding box south of St Thomas in UTM coordinates
  # # Based on your plot showing UTM coordinates (x: ~0-400000, y: ~1900000-2100000)
  # # Let's take a small area in the south-central part
  # # test_extent <- ext(280000, 305000, 2010000, 2030000)
  # # test_extent <- ext(260000, 310000, 2000000, 2050000)
  # # test_extent <- ext(200000, 320000, 2000000, 2060000)
  # # test_extent <- ext(20000, 340000, 2000000, 2060000)
  # test_extent <- ext(per_at_max_hsig_raster)  # 50km x 50km area
  # 
  # # Crop all rasters to test area
  # per_test <- crop(per_at_max_hsig_raster, test_extent)
  # 
  # plot(per_test)
  # 
  # depth_test <- crop(depth_raster, test_extent)
  # hsig_test <- crop(max_hsig_raster, test_extent)
  # 
  # cat(sprintf("Test area dimensions: %d x %d = %d cells\n", 
  #             nrow(per_test), ncol(per_test), ncell(per_test)))
  # 
  # # Optimized solve_L function (same as before)
  # solve_L_efficient <- function(T, h, g = 9.81) {
  #   if (is.na(T) || is.na(h) || T <= 0 || h <= 0) {
  #     return(NA)
  #   }
  #   
  #   # Define the dispersion relation function
  #   f <- function(L) {
  #     k <- 2 * pi / L
  #     return(((2 * pi) / T)^2 - g * k * tanh(k * h))
  #   }
  #   
  #   # Smart interval selection based on water depth regime
  #   L_shallow <- T * sqrt(g * h)  # Shallow water approximation
  #   L_deep <- g * T^2 / (2 * pi)  # Deep water approximation
  #   
  #   # Create intervals based on the physics
  #   L_min <- min(L_shallow, L_deep) * 0.1   
  #   L_max <- max(L_shallow, L_deep) * 5     
  #   
  #   # Primary interval
  #   interval_main <- c(L_min, L_max)
  #   
  #   # Try main interval first
  #   solution <- tryCatch({
  #     uniroot(f, interval = interval_main, tol = 1e-8)
  #   }, error = function(e) NULL)
  #   
  #   if (!is.null(solution)) {
  #     return(solution$root)
  #   }
  #   
  #   # Fallback intervals if main fails
  #   fallback_intervals <- list(
  #     c(0.1, 500),      
  #     c(0.01, 1000)     
  #   )
  #   
  #   for (interval in fallback_intervals) {
  #     solution <- tryCatch({
  #       uniroot(f, interval = interval, tol = 1e-8)
  #     }, error = function(e) NULL)
  #     
  #     if (!is.null(solution)) {
  #       return(solution$root)
  #     }
  #   }
  #   
  #   return(NA)
  # }
  # 
  # # Test the wavelength calculation on small subset
  # cat("Testing wavelength calculation...\n")
  # start_time <- Sys.time()
  # 
  # # Extract values for processing
  # per_vals <- values(per_test)
  # depth_vals <- values(depth_test)
  # 
  # # Pre-filter to valid cells only
  # valid_mask <- !is.na(per_vals) & !is.na(depth_vals) & per_vals > 0 & depth_vals > 0
  # valid_indices <- which(valid_mask)
  # valid_cells <- length(valid_indices)
  # cat(sprintf("Valid cells to process: %d\n", valid_cells))
  # 
  # # Process only valid cells
  # L_vals_valid <- mapply(solve_L_efficient, 
  #                        per_vals[valid_indices], 
  #                        depth_vals[valid_indices], 
  #                        MoreArgs = list(g = g), USE.NAMES = FALSE)
  # 
  # # Reconstruct full vector
  # L_vals <- rep(NA, length(per_vals))
  # L_vals[valid_indices] <- L_vals_valid
  # 
  #   # Create wavelength raster
  # L_test <- per_test
  # values(L_test) <- L_vals
  # names(L_test) <- "wavelength"
  # 
  # end_time <- Sys.time()
  # elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  # 
  # # Check results
  # na_count <- sum(is.na(L_vals) & !is.na(per_vals) & !is.na(depth_vals))
  # success_rate <- (valid_cells - na_count) / valid_cells
  # 
  # cat("\nTest Results:\n")
  # cat(sprintf("  Processing time: %.2f seconds\n", elapsed))
  # cat(sprintf("  Successful solutions: %d/%d (%.1f%%)\n", 
  #             valid_cells - na_count, valid_cells, success_rate * 100))
  # 
  # if (success_rate > 0) {
  #   wavelength_range <- range(L_vals, na.rm = TRUE)
  #   cat(sprintf("  Wavelength range: %.1f to %.1f m\n", 
  #               wavelength_range[1], wavelength_range[2]))
  #   
  #   # Quick sanity check
  #   sample_vals <- na.omit(L_vals)[1:min(10, sum(!is.na(L_vals)))]
  #   cat(sprintf("  Sample wavelengths: %s m\n", 
  #               paste(round(sample_vals, 1), collapse = ", ")))
  #   
  #   # Calculate BOV for test area
  #   cat("\nTesting BOV calculation...\n")
  #   k_h_test <- 2 * pi * depth_test / L_test
  #   cosh_term_test <- cosh(k_h_test)
  #   bov_test <- (hsig_test * g * per_test) / (2 * L_test * cosh_term_test)
  #   names(bov_test) <- "BOV_test"
  #   
  #   bov_range <- range(values(bov_test), na.rm = TRUE)
  #   cat(sprintf("  BOV range: %.3f to %.3f m/s\n", bov_range[1], bov_range[2]))
  #   
  #   # Detailed BOV analysis to find the outliers
  #   bov_vals <- values(bov_test)
  #   bov_valid <- bov_vals[!is.na(bov_vals)]
  #   
  #   cat("\nDetailed BOV Analysis:\n")
  #   cat(sprintf("  Mean BOV: %.3f m/s\n", mean(bov_valid)))
  #   cat(sprintf("  Median BOV: %.3f m/s\n", median(bov_valid)))
  #   cat(sprintf("  95th percentile: %.3f m/s\n", quantile(bov_valid, 0.95)))
  #   cat(sprintf("  99th percentile: %.3f m/s\n", quantile(bov_valid, 0.99)))
  #   cat(sprintf("  99.9th percentile: %.3f m/s\n", quantile(bov_valid, 0.999)))
  #   
  #   # Count extreme values
  #   extreme_count <- sum(bov_valid > 2.0)  # > 2 m/s is very high
  #   very_extreme_count <- sum(bov_valid > 10.0)  # > 10 m/s is unreasonable
  #   
  #   cat(sprintf("  Cells with BOV > 2.0 m/s: %d (%.2f%%)\n", 
  #               extreme_count, extreme_count/length(bov_valid)*100))
  #   cat(sprintf("  Cells with BOV > 10.0 m/s: %d (%.2f%%)\n", 
  #               very_extreme_count, very_extreme_count/length(bov_valid)*100))
  #   
  #   # Show the most extreme values and their conditions
  #   if (very_extreme_count > 0) {
  #     cat("\nExtreme value diagnosis:\n")
  #     extreme_indices <- which(bov_vals > 10.0)
  #     n_show <- min(5, length(extreme_indices))
  #     
  #     for (i in 1:n_show) {
  #       idx <- extreme_indices[i]
  #       bov_val <- bov_vals[idx]
  #       hsig_val <- values(hsig_test)[idx]
  #       per_val <- values(per_test)[idx]
  #       depth_val <- values(depth_test)[idx]
  #       L_val <- values(L_test)[idx]
  #       
  #       cat(sprintf("  Cell %d: BOV=%.1f, H=%.2f, T=%.1f, h=%.1f, L=%.1f\n",
  #                   i, bov_val, hsig_val, per_val, depth_val, L_val))
  #     }
  #   }
  #   
  #   # Estimate time for full domain
  #   cells_per_second <- valid_cells / elapsed
  #   total_cells_estimate <- sum(!is.na(values(per_at_max_hsig_raster)) & 
  #                                 !is.na(values(depth_raster)))
  #   estimated_time_minutes <- total_cells_estimate / cells_per_second / 60
  #   
  #   cat(sprintf("\nFull domain estimate:\n"))
  #   cat(sprintf("  Total cells to process: ~%d\n", total_cells_estimate))
  #   cat(sprintf("  Estimated processing time: %.1f minutes\n", estimated_time_minutes))
  #   
  #   # Plot test results with better visualization
  #   cat("\nPlotting test results...\n")
  #   par(mfrow = c(2, 3))
  #   plot(L_test, main = "Wavelength (m)")
  #   
  #   # # Plot BOV with reasonable range (clip extreme outliers)
  #   # bov_reasonable <- bov_test
  #   # bov_vals_clipped <- pmax(pmin(values(bov_test), 2.0), 0)  # Clip to 0-2 m/s
  #   # values(bov_reasonable) <- bov_vals_clipped
  #   # plot(bov_reasonable, main = "BOV (m/s) - Clipped to 0-2")
  #   
  #   # Dynamic upper limit based on 99th percentile
  #   bov_upper_limit <- quantile(bov_valid, 0.99, na.rm = TRUE)
  #   bov_reasonable <- bov_test
  #   bov_vals_clipped <- pmax(pmin(values(bov_test), bov_upper_limit), 0)
  #   values(bov_reasonable) <- bov_vals_clipped
  #   plot(bov_reasonable, main = sprintf("BOV (m/s) - 0 to %.2f", bov_upper_limit))
  #   
  #   # # Also show the original with full range
  #   # plot(bov_test, main = "BOV (m/s) - Full Range")
  #   
  #   # Create depth plot with outliers highlighted, limited color scale, and stars
  #   outlier_mask <- values(bov_test) > 5.0 & !is.na(values(bov_test))
  #   outlier_depths <- depth_test
  #   values(outlier_depths)[!outlier_mask] <- NA
  #   
  #   # Plot background in gray first
  #   plot(depth_test, main = "Depth (m) at BOV Outliers > 5 m/s", col = "gray", legend = FALSE)
  #   
  #   # Get the range of outlier depths and plot with limited color scale
  #   if(sum(outlier_mask, na.rm = TRUE) > 0) {
  #     outlier_range <- range(values(outlier_depths), na.rm = TRUE)
  #     plot(outlier_depths, add = TRUE, zlim = outlier_range)
  #     
  #     # Add red stars at outlier locations
  #     outlier_coords <- xyFromCell(bov_test, which(outlier_mask))
  #     points(outlier_coords, pch = 8, col = "red", cex = 0.5)
  #   }  
  #   
  #   plot(per_test, main = "Wave Period (s)")  
  #   # plot(depth_test, main = "Depth (m)")
  #   
  #   # Clip depth for better visualization
  #   depth_clipped <- depth_test
  #   values(depth_clipped) <- pmin(values(depth_test), 90)
  #   plot(depth_clipped, main = "Depth (m) - Clipped to 90m")
  #   
  #   plot(hsig_test, main = "Wave Height (m)")
  #   par(mfrow = c(1, 1))
  #   
  #   # Create histogram of BOV values to see distribution
  #   hist(bov_valid[bov_valid < 5], breaks = 50,
  #        main = "BOV Distribution (< 5 m/s)", xlab = "BOV (m/s)")
  #   
  #   cat("\nReasonable BOV range should be ~0.01-2.0 m/s for most conditions\n")
  #   cat("Values > 5 m/s suggest numerical issues (very shallow water, extreme waves, or wavelength errors)\n")
  #   
  # } else {
  #   cat("  ERROR: No successful solutions! Check input data.\n")
  #   
  #   # Debug info
  #   cat("\nDebug info:\n")
  #   cat(sprintf("  Period range: %.2f to %.2f s\n", 
  #               min(per_vals, na.rm = TRUE), max(per_vals, na.rm = TRUE)))
  #   cat(sprintf("  Depth range: %.2f to %.2f m\n", 
  #               min(depth_vals, na.rm = TRUE), max(depth_vals, na.rm = TRUE)))
  # }

  # VERSION with chunking
  # Chunked BOV Processing for Full Domain (8 chunks)
  # Constants
  g <- 9.81
  
  # Make depth positive
  depth_raster <- abs(bathy_final)
  
  cat("Starting chunked BOV processing for full domain...\n")
  
  # Get full extent and divide into 8 chunks (2x4 grid)
  full_extent <- ext(per_at_max_hsig_raster)
  x_range <- c(full_extent[1], full_extent[2])
  y_range <- c(full_extent[3], full_extent[4])
  
  # Create 8 chunks (2 rows x 4 columns)
  n_x_chunks <- 4
  n_y_chunks <- 2
  
  x_breaks <- seq(x_range[1], x_range[2], length.out = n_x_chunks + 1)
  y_breaks <- seq(y_range[1], y_range[2], length.out = n_y_chunks + 1)
  
  cat(sprintf("Full domain extent: %.0f to %.0f, %.0f to %.0f\n", 
              x_range[1], x_range[2], y_range[1], y_range[2]))
  
  # Optimized solve_L function
  solve_L_efficient <- function(T, h, g = 9.81) {
    if (is.na(T) || is.na(h) || T <= 0 || h <= 0) {
      return(NA)
    }
    
    # Define the dispersion relation function
    f <- function(L) {
      k <- 2 * pi / L
      return(((2 * pi) / T)^2 - g * k * tanh(k * h))
    }
    
    # Smart interval selection based on water depth regime
    L_shallow <- T * sqrt(g * h)  # Shallow water approximation
    L_deep <- g * T^2 / (2 * pi)  # Deep water approximation
    
    # Create intervals based on the physics
    L_min <- min(L_shallow, L_deep) * 0.1   
    L_max <- max(L_shallow, L_deep) * 5     
    
    # Primary interval
    interval_main <- c(L_min, L_max)
    
    # Try main interval first
    solution <- tryCatch({
      uniroot(f, interval = interval_main, tol = 1e-8)
    }, error = function(e) NULL)
    
    if (!is.null(solution)) {
      return(solution$root)
    }
    
    # Fallback intervals if main fails
    fallback_intervals <- list(
      c(0.1, 500),      
      c(0.01, 1000)     
    )
    
    for (interval in fallback_intervals) {
      solution <- tryCatch({
        uniroot(f, interval = interval, tol = 1e-8)
      }, error = function(e) NULL)
      
      if (!is.null(solution)) {
        return(solution$root)
      }
    }
    
    return(NA)
  }
  
  # Initialize lists to store results
  L_chunks <- list()
  bov_chunks <- list()
  total_start_time <- Sys.time()
  
  # Process each chunk
  chunk_num <- 0
  for (i in 1:n_y_chunks) {
    for (j in 1:n_x_chunks) {
      chunk_num <- chunk_num + 1
      
      # Define chunk extent
      chunk_extent <- ext(x_breaks[j], x_breaks[j+1], 
                          y_breaks[i], y_breaks[i+1])
      
      cat(sprintf("\n=== Processing Chunk %d of 8 ===\n", chunk_num))
      cat(sprintf("Extent: %.0f to %.0f, %.0f to %.0f\n", 
                  chunk_extent[1], chunk_extent[2], 
                  chunk_extent[3], chunk_extent[4]))
      
      chunk_start_time <- Sys.time()
      
      # Crop rasters to chunk
      per_chunk <- crop(per_at_max_hsig_raster, chunk_extent)
      depth_chunk <- crop(depth_raster, chunk_extent)
      hsig_chunk <- crop(max_hsig_raster, chunk_extent)
      
      cat(sprintf("Chunk dimensions: %d x %d = %d cells\n", 
                  nrow(per_chunk), ncol(per_chunk), ncell(per_chunk)))
      
      # Extract values for processing
      per_vals <- values(per_chunk)
      depth_vals <- values(depth_chunk)
      
      # Pre-filter to valid cells only
      valid_mask <- !is.na(per_vals) & !is.na(depth_vals) & per_vals > 0 & depth_vals > 0
      valid_indices <- which(valid_mask)
      valid_cells <- length(valid_indices)
      
      cat(sprintf("Valid cells to process: %d\n", valid_cells))
      
      if (valid_cells > 0) {
        # Process only valid cells
        L_vals_valid <- mapply(solve_L_efficient, 
                               per_vals[valid_indices], 
                               depth_vals[valid_indices], 
                               MoreArgs = list(g = g), USE.NAMES = FALSE)
        
        # Reconstruct full vector
        L_vals <- rep(NA, length(per_vals))
        L_vals[valid_indices] <- L_vals_valid
        
        # Create wavelength raster
        L_chunk <- per_chunk
        values(L_chunk) <- L_vals
        names(L_chunk) <- paste0("wavelength_chunk_", chunk_num)
        
        # Calculate BOV for chunk
        k_h_chunk <- 2 * pi * depth_chunk / L_chunk
        cosh_term_chunk <- cosh(k_h_chunk)
        bov_chunk <- (hsig_chunk * g * per_chunk) / (2 * L_chunk * cosh_term_chunk)
        names(bov_chunk) <- paste0("BOV_chunk_", chunk_num)
        
        # Store results
        L_chunks[[chunk_num]] <- L_chunk
        bov_chunks[[chunk_num]] <- bov_chunk
        
        # Check results
        na_count <- sum(is.na(L_vals) & !is.na(per_vals) & !is.na(depth_vals))
        success_rate <- (valid_cells - na_count) / valid_cells
        
        chunk_end_time <- Sys.time()
        chunk_elapsed <- as.numeric(difftime(chunk_end_time, chunk_start_time, units = "secs"))
        
        cat(sprintf("Chunk processing time: %.1f seconds\n", chunk_elapsed))
        cat(sprintf("Success rate: %.1f%%\n", success_rate * 100))
        
        if (success_rate > 0) {
          wavelength_range <- range(L_vals, na.rm = TRUE)
          bov_vals <- values(bov_chunk)
          bov_range <- range(bov_vals, na.rm = TRUE)
          
          cat(sprintf("Wavelength range: %.1f to %.1f m\n", 
                      wavelength_range[1], wavelength_range[2]))
          cat(sprintf("BOV range: %.3f to %.3f m/s\n", bov_range[1], bov_range[2]))
        }
      } else {
        cat("No valid cells in this chunk - skipping\n")
        L_chunks[[chunk_num]] <- NULL
        bov_chunks[[chunk_num]] <- NULL
      }
    }
  }
  
  total_end_time <- Sys.time()
  total_elapsed <- as.numeric(difftime(total_end_time, total_start_time, units = "mins"))
  
  cat(sprintf("\n=== PROCESSING COMPLETE ===\n"))
  cat(sprintf("Total processing time: %.1f minutes\n", total_elapsed))
  
  # Merge chunks back together
  cat("Merging chunks...\n")
  
  # Filter out NULL chunks
  valid_L_chunks <- L_chunks[!sapply(L_chunks, is.null)]
  valid_bov_chunks <- bov_chunks[!sapply(bov_chunks, is.null)]
  
  if (length(valid_L_chunks) > 0) {
    # Merge wavelength rasters
    L_full <- do.call(mosaic, valid_L_chunks)
    names(L_full) <- "wavelength_full"
    
    # Merge BOV rasters
    bov_full <- do.call(mosaic, valid_bov_chunks)
    names(bov_full) <- "BOV_full"
    
    cat("Merging complete!\n")
    
    # Final analysis
    bov_vals_all <- values(bov_full)
    bov_valid_all <- bov_vals_all[!is.na(bov_vals_all)]
    
    cat("\n=== FINAL RESULTS ===\n")
    cat(sprintf("Total valid BOV cells: %d\n", length(bov_valid_all)))
    cat(sprintf("BOV range: %.3f to %.3f m/s\n", 
                min(bov_valid_all), max(bov_valid_all)))
    cat(sprintf("BOV mean: %.3f m/s\n", mean(bov_valid_all)))
    cat(sprintf("BOV median: %.3f m/s\n", median(bov_valid_all)))
    
    # Count outliers
    outlier_count <- sum(bov_valid_all > 5.0)
    cat(sprintf("Cells with BOV > 5.0 m/s: %d (%.3f%%)\n", 
                outlier_count, outlier_count/length(bov_valid_all)*100))
    
    # # Save results
    # cat("\nSaving results...\n")
    # writeRaster(L_full, "wavelength_full_domain.tif", overwrite = TRUE)
    # writeRaster(bov_full, "BOV_full_domain.tif", overwrite = TRUE)
    # 
    # cat("Results saved as wavelength_full_domain.tif and BOV_full_domain.tif\n")
    # 
    # # Quick plot
    # cat("Creating summary plot...\n")
    # par(mfrow = c(1, 2))
    # plot(L_full, main = "Wavelength (m) - Full Domain")
    # 
    # # Plot BOV with reasonable range
    # bov_upper_limit <- quantile(bov_valid_all, 0.99, na.rm = TRUE)
    # bov_plot <- bov_full
    # values(bov_plot) <- pmax(pmin(values(bov_full), bov_upper_limit), 0)
    # plot(bov_plot, main = sprintf("BOV (m/s) - 0 to %.2f", bov_upper_limit))
    # par(mfrow = c(1, 1))
    
  } else {
    cat("ERROR: No valid chunks processed!\n")
  }
  
  # plot(bov_full)
  # 
  # BOV_clamp = clamp(bov_full, lower = 0, upper = 3)
  # plot(BOV_clamp)
  # 
  # bathy_clamper = clamp(bathy_final, lower = -50, upper = 0)
  # plot(bathy_clamper)
  
  ################################## Import ERDDAP waves ##################################
  
  #above wave section was data supplied directly from Dr. Miguel Canals, and is to be used for
  #   max BOV and mean/max Hsig. this section from ERDDAP is for mean wave direction only
  
  # VERSION where SWAN data is from ERDDAP (MATLAB output)
  
  SWAN_target_res <- 180
  
  DO_GAP_FILLING_180M <- TRUE   # Set to FALSE to disable gap filling for main 180m raster
  DO_GAP_FILLING_1KM <- FALSE   # Set to TRUE to enable gap filling for 1km low-res areas

  # Load the composite wave direction data from MATLAB
  wave_summary <- read.csv(here("output", "swan_composite_direction.csv"))
  
  # Clean the data
  wave_summary <- wave_summary[!is.na(wave_summary$lon) & !is.na(wave_summary$lat), ]
  
  cat(sprintf("Loaded %d SWAN grid points\n", nrow(wave_summary)))
  
  # Create SpatVector from wave data
  wave_points <- vect(wave_summary, 
                      geom = c("lon", "lat"), 
                      crs = "EPSG:4326")
  
  # Project to match bathymetry
  wave_points_proj <- project(wave_points, crs(bathy_final))
  
  cat(sprintf("Total points to process: %d\n", nrow(wave_points_proj)))
  
  
  create_wave_raster <- function(wave_points, var_name, bathy_template, target_resolution, 
                                 do_gap_filling = TRUE) {
    cat(sprintf("Creating %dm raster for %s...\n", target_resolution, var_name))
    
    # Filter to valid values
    values <- values(wave_points)[[var_name]]
    valid_idx <- !is.na(values) & values >= 0
    
    if (grepl("dir", var_name, ignore.case = TRUE)) {
      valid_idx <- valid_idx & values <= 360
    }
    
    wave_points_valid <- wave_points[valid_idx]
    cat(sprintf("Using %d valid points\n", sum(valid_idx)))
    
    # Create template with target resolution
    bathy_ext <- ext(bathy_template)
    template_raster <- rast(bathy_ext, 
                            resolution = target_resolution,
                            crs = crs(bathy_template))
    
    n_cells <- ncell(template_raster)
    cat(sprintf("Grid: %d x %d cells (%.1f million)\n", 
                ncol(template_raster), nrow(template_raster), n_cells/1e6))
    
    # Rasterize wave points to template
    # For directional data, use circular-aware rasterization
    if (grepl("dir", var_name, ignore.case = TRUE)) {
      cat("Using circular-aware rasterization for directional data...\n")
      
      # Convert angles to sin/cos components
      angles_rad <- values(wave_points_valid)[[var_name]] * pi / 180
      wave_points_valid$sin_comp <- sin(angles_rad)
      wave_points_valid$cos_comp <- cos(angles_rad)
      
      # Rasterize sin and cos components separately
      sin_rast <- rasterize(wave_points_valid, template_raster, 
                            field = "sin_comp", 
                            fun = mean)
      cos_rast <- rasterize(wave_points_valid, template_raster, 
                            field = "cos_comp", 
                            fun = mean)
      
      # Convert back to angles
      result <- atan2(sin_rast, cos_rast) * 180 / pi
      result[result < 0] <- result[result < 0] + 360
      
    } else {
      # For non-directional data, use regular mean
      result <- rasterize(wave_points_valid, template_raster, 
                          field = var_name, 
                          fun = mean)
    }
    
    # Gap filling with focal operations
    if (do_gap_filling) {
      cat("Gap filling...\n")
      
      # For directional data, use circular-aware gap filling
      if (grepl("dir", var_name, ignore.case = TRUE)) {
        # Convert to radians for calculation
        result_rad <- result * pi / 180
        
        # Focal mean of sin and cos components with multiple passes
        sin_comp <- app(result_rad, sin)
        cos_comp <- app(result_rad, cos)
        
        # Iterative gap filling with increasing window sizes
        for (window_size in c(3, 5, 7, 9)) {
          sin_comp <- focal(sin_comp, w = window_size, fun = mean, na.policy = "only", na.rm = TRUE)
          cos_comp <- focal(cos_comp, w = window_size, fun = mean, na.policy = "only", na.rm = TRUE)
        }
        
        # Convert back to degrees
        result <- atan2(sin_comp, cos_comp) * 180 / pi
        result[result < 0] <- result[result < 0] + 360
        
      } else {
        # For non-directional data, use regular mean with multiple passes
        for (window_size in c(3, 5, 7, 9)) {
          result <- focal(result, w = window_size, fun = mean, na.policy = "only", na.rm = TRUE)
        }
      }
    } else {
      cat("Skipping gap filling\n")
    }
    
    # Don't apply mask here - will be applied at the end at native bathymetry resolution
    
    names(result) <- var_name
    return(result)
  }
  
  
  cat(sprintf("\n=== Creating wave direction raster at %dm resolution ===\n", SWAN_target_res))
  
  # Create main wave direction raster at target resolution
  dir_erddap_raster <- create_wave_raster(wave_points_proj, "mean_dir", bathy_final, 
                                          SWAN_target_res, do_gap_filling = DO_GAP_FILLING_180M)
  
  # Create coarse 1km raster for patching low-resolution areas
  cat("\n=== Creating 1km raster for low-resolution areas ===\n")
  dir_1km_raster <- create_wave_raster(wave_points_proj, "mean_dir", bathy_final, 1000, 
                                       do_gap_filling = DO_GAP_FILLING_1KM)
  
  # Resample 1km raster to match the target resolution grid
  dir_1km_resampled <- resample(dir_1km_raster, dir_erddap_raster, method = "near")
  
  # Define low-resolution areas to patch in
  low_res_areas <- list(
    ext(-30000, 0, 2000000, 2025000),   # Mona Island
    ext(20000, 40000, 2030000, 2045000)  # Desecheo Island
  )
  
  # Patch in the 1km data for specified areas
  cat("\n=== Patching low-resolution areas ===\n")
  for (i in seq_along(low_res_areas)) {
    area_ext <- low_res_areas[[i]]
    cat(sprintf("Patching area %d: [%.0f, %.0f, %.0f, %.0f]\n", 
                i, area_ext[1], area_ext[2], area_ext[3], area_ext[4]))
    
    # Crop both rasters to the area
    patch_from_1km <- crop(dir_1km_resampled, area_ext)
    patch_to_target <- crop(dir_erddap_raster, area_ext)
    
    # Replace the target area with 1km data (where 1km has values)
    patch_merged <- cover(patch_from_1km, patch_to_target)
    
    # Put the patched area back into the full raster
    # Create a mask for the area
    mask_raster <- rast(dir_erddap_raster)
    values(mask_raster) <- 0
    mask_raster <- crop(mask_raster, area_ext)
    values(mask_raster) <- 1
    mask_raster <- extend(mask_raster, dir_erddap_raster, fill = 0)
    
    # Apply the patch
    dir_erddap_raster <- ifel(mask_raster == 1, 
                              extend(patch_merged, dir_erddap_raster, fill = NA),
                              dir_erddap_raster)
  }
  
  # Get extent and properties from bathy_final
  bathy_ext <- ext(bathy_final)
  bathy_crs <- crs(bathy_final)
  bathy_res <- res(bathy_final)
  
  # Create template raster matching bathy_final exactly
  template_raster <- rast(bathy_ext, resolution = bathy_res, crs = bathy_crs)
  
  # Resample wave raster to match bathy_final
  # Use nearest neighbor for directional data to preserve circular data integrity
  cat("\n=== Resampling to final resolution using nearest neighbor ===\n")
  dir_erddap_raster <- resample(dir_erddap_raster, template_raster, method = "near")
    
  ################################## distance from shore ##################################
  
  # Create landmask from bathy_crm_2024 which has positive values for land
  cat("Creating landmask from bathy_crm_2024...\n")
  # Load the raster with land elevations
  bathy_crm_2024 <- rast(here('data/exportImage.tiff'))
  # Create simple landmask: 1 for land (>0), 0 for water (<=0), NA stays NA
  landmask <- bathy_crm_2024 > 0
  landmask <- as.numeric(landmask)
  names(landmask) <- "landmask"
  cat("Landmask complete: 1=Land, 0=Water, NA=No Data\n")
  plot(landmask)
  
  cat("Aligning landmask with bathy_final grid structure...\n")
  # Assume bathy_final is already loaded - if not, load it first
  # bathy_final <- rast("path/to/your/bathy_final.tif")
  
  # Get the CRS of bathy_final
  target_crs <- crs(bathy_final)
  cat("Target CRS:", target_crs, "\n")
  
  # Create a template raster based on bathy_final structure
  bathy_template <- bathy_final
  # Set all valid bathymetry cells to 0 (water), keep NA where bathy is NA
  bathy_template[!is.na(bathy_template)] <- 0
  
  # Reproject and resample landmask to exactly match bathy_final grid
  landmask_aligned <- project(landmask, bathy_template, method = "near")
  
  # Now combine: where landmask shows land (1), set to 1; otherwise keep as 0 or NA
  # This creates a landmask that matches bathy_final's exact grid
  landmask_on_bathy <- bathy_template
  landmask_on_bathy[landmask_aligned == 1] <- 1
  
  names(landmask_on_bathy) <- "landmask_on_bathy_grid"
  
  cat("Landmask alignment complete\n")
  cat("Grid dimensions - bathy_final:", dim(bathy_final), "\n")
  cat("Grid dimensions - aligned landmask:", dim(landmask_on_bathy), "\n")
  
  # Plot comparison
  par(mfrow = c(1, 2))
  plot(bathy_final, main = "bathy_final")
  plot(landmask_on_bathy, main = "Landmask on bathy_final grid")
  par(mfrow = c(1, 1))
  
  cat("Calculating distance from each valid bathy_final cell to nearest land...\n")
  
  # Create binary raster for land cells only (1 = land, NA = everything else)
  land_cells_binary <- landmask_on_bathy
  land_cells_binary[land_cells_binary == 0] <- NA  # Set water to NA
  # land_cells_binary now has: 1 = land, NA = water/no-data
  
  # Calculate distance from all cells (including water) to nearest land cell
  distance_to_land <- distance(land_cells_binary)
  
  # Mask the distance raster to only show distances for valid bathy_final cells
  distance_for_bathy_cells <- mask(distance_to_land, bathy_final)
  names(distance_for_bathy_cells) <- "distance_to_land_m"
  
  cat("Distance calculation complete\n")
  cat("Distance units are in map units of the target CRS\n")
  
  # Plot the results
  par(mfrow = c(2, 2))
  plot(bathy_final, main = "bathy_final (original)")
  plot(landmask_on_bathy, main = "Landmask on bathy grid")
  plot(distance_to_land, main = "Distance to land (all cells)")
  plot(distance_for_bathy_cells, main = "Distance for valid bathy cells only")
  par(mfrow = c(1, 1))
  
  cat("Checking if distance calculation introduced additional NAs...\n")
  
  # Count NAs in original bathy_final
  na_count_bathy <- sum(is.na(values(bathy_final)))
  total_cells_bathy <- ncell(bathy_final)
  valid_cells_bathy <- total_cells_bathy - na_count_bathy
  
  # Count NAs in distance_for_bathy_cells  
  na_count_distance <- sum(is.na(values(distance_for_bathy_cells)))
  total_cells_distance <- ncell(distance_for_bathy_cells)
  valid_cells_distance <- total_cells_distance - na_count_distance
  
  cat("bathy_final: ", valid_cells_bathy, " valid cells,", na_count_bathy, "NAs\n")
  cat("distance_for_bathy_cells: ", valid_cells_distance, " valid cells,", na_count_distance, "NAs\n")
  
  if (na_count_distance > na_count_bathy) {
    additional_nas <- na_count_distance - na_count_bathy
    cat("WARNING: Distance calculation introduced", additional_nas, "additional NAs!\n")
  } else if (na_count_distance == na_count_bathy) {
    cat("GOOD: No additional NAs introduced - same NA pattern as bathy_final\n")
  } else {
    cat("UNEXPECTED: Distance raster has fewer NAs than bathy_final\n")
  }
  
  # Convert to km if units are in meters
  if (grepl("metre|meter|m", target_crs, ignore.case = TRUE)) {
    dist_to_land_raster <- distance_for_bathy_cells / 1000
    names(dist_to_land_raster) <- "distance_to_land_km"
    
    # Final plot in km
    plot(dist_to_land_raster, main = "Distance to Land for Bathy Cells (km)")
  }
  
  
  ################################## distance from deep water ##################################

  # Set depth threshold
  depth_threshold <- 70  # meters
  
  # Create binary raster: 1 for deep water (deeper than threshold), NA otherwise  
  # Assuming negative values for depth (e.g., -35m is deeper than -25m)
  deep_water <- bathy_final < -depth_threshold
  deep_water <- as.numeric(deep_water)
  deep_water[deep_water == 0] <- NA  # Set shallow areas to NA
  deep_water[is.na(bathy_final)] <- NA  # Preserve original NAs
  
  names(deep_water) <- paste0("deeper_than_", depth_threshold, "m")
  
  # Check how many deep water cells found
  n_deep_cells <- sum(!is.na(values(deep_water)))
  cat("Number of cells deeper than", depth_threshold, "m:", n_deep_cells, "\n")
  
  # Calculate distance to deep water
  if (n_deep_cells > 0) {
    distance_to_deep <- distance(deep_water)
    distance_to_deep_masked <- mask(distance_to_deep, bathy_final)
    
    # Convert to kilometers
    distance_to_deep_raster <- distance_to_deep_masked / 1000
    names(distance_to_deep_raster) <- paste0("distance_to_", depth_threshold, "m_depth_km")
    
    # Plot results
    par(mfrow = c(1, 3))
    plot(bathy_final, main = "Bathymetry")
    plot(deep_water, main = paste("Water >", depth_threshold, "m deep"))
    plot(distance_to_deep_raster, main = paste("Distance to", depth_threshold, "m depth (km)"))
    par(mfrow = c(1, 1))
    
    cat("Distance calculation complete!\n")
  } else {
    cat("No areas found deeper than", depth_threshold, "m\n")
  }
  
  plot(distance_to_deep_raster)
  
  ################################## lon & lat ##################################
  
  # Get coordinates from bathy_final
  coords <- xyFromCell(bathy_final, 1:ncell(bathy_final))
  
  # Create longitude raster
  lon_raster <- bathy_final
  values(lon_raster) <- coords[, 1]
  names(lon_raster) <- "longitude"
  
  # Create latitude raster
  lat_raster <- bathy_final
  values(lat_raster) <- coords[, 2]
  names(lat_raster) <- "latitude"
  
  # Mask to match bathy_final (only keep values where bathy exists)
  lon_raster <- mask(lon_raster, bathy_final)
  lat_raster <- mask(lat_raster, bathy_final)
  
  cat("Longitude and latitude rasters created\n")
  cat("Longitude range:", range(values(lon_raster), na.rm = TRUE), "\n")
  cat("Latitude range:", range(values(lat_raster), na.rm = TRUE), "\n")
  
  # Plot to verify
  par(mfrow = c(1, 2))
  plot(lon_raster, main = "Longitude")
  plot(lat_raster, main = "Latitude")
  par(mfrow = c(1, 1))
  
  ################################## Import SST ##################################
  
  sst_summary <- read.csv(here("output", "mur_sst_summary_for_R.csv"))
  
  # 2. Load coordinate information
  lon_vec <- read.csv(here("output", "mur_sst_longitude.csv"))$longitude
  lat_vec <- read.csv(here("output", "mur_sst_latitude.csv"))$latitude
  grid_info <- read.csv(here("output", "mur_sst_grid_info.csv"))
  
  # Get grid dimensions
  nlon <- length(lon_vec)
  nlat <- length(lat_vec)
  lon_range <- range(lon_vec)
  lat_range <- range(lat_vec)
  
  # Function to load matrix and create raster
  create_raster <- function(filename, var_name) {
    cat(sprintf("Loading %s...\n", filename))
    matrix_data <- as.matrix(read.csv(here("output", filename), header = FALSE))
    
    # Flip matrix vertically to match MATLAB orientation
    matrix_data <- matrix_data[nrow(matrix_data):1, ]
    
    # Create SpatRaster using terra
    r <- rast(matrix_data,
              extent = ext(min(lon_vec), max(lon_vec), min(lat_vec), max(lat_vec)),
              crs = "EPSG:4326")
    
    names(r) <- var_name
    return(r)
  }
  
  # Create rasters for SST variables
  mean_sst_raster <- create_raster("mean_sst_matrix.csv", "Mean_SST")
  max_sst_raster <- create_raster("max_sst_matrix.csv", "Max_SST")
  min_sst_raster <- create_raster("min_sst_matrix.csv", "Min_SST")
  std_sst_raster <- create_raster("std_sst_matrix.csv", "Std_SST")
  
  projected_crs = crs(bathy_final)
  
  #reproject
  mean_sst_raster = project(mean_sst_raster, projected_crs)
  max_sst_raster = project(max_sst_raster, projected_crs)
  min_sst_raster = project(min_sst_raster, projected_crs)
  std_sst_raster = project(std_sst_raster, projected_crs)
  
  # Resample sst rasters to match bathy_final
  mean_sst_raster <- resample(mean_sst_raster, template_raster, method = "bilinear")
  max_sst_raster <- resample(max_sst_raster, template_raster, method = "bilinear")
  min_sst_raster <- resample(min_sst_raster, template_raster, method = "bilinear")
  std_sst_raster <- resample(std_sst_raster, template_raster, method = "bilinear")
  
  ################################## Import PAR ##################################
  
  par_summary <- read.csv(here("output", "par_summary_for_R.csv"))
  
  # Load PAR coordinate information (should be same as SST, but loading separately for completeness)
  par_lon_vec <- read.csv(here("output", "par_longitude.csv"))$longitude
  par_lat_vec <- read.csv(here("output", "par_latitude.csv"))$latitude
  par_grid_info <- read.csv(here("output", "par_grid_info.csv"))
  
  # Function to create PAR rasters
  create_par_raster <- function(filename, var_name) {
    cat(sprintf("Loading PAR %s...\n", filename))
    matrix_data <- as.matrix(read.csv(here("output", filename), header = FALSE))
    
    # # Flip matrix vertically to match MATLAB orientation
    # matrix_data <- matrix_data[nrow(matrix_data):1, ]
    
    # Create SpatRaster using terra
    r <- rast(matrix_data,
              extent = ext(min(par_lon_vec), max(par_lon_vec), min(par_lat_vec), max(par_lat_vec)),
              crs = "EPSG:4326")
    # extent = ext(min(par_lon_vec), max(par_lon_vec), min(par_lat_vec), max(par_lat_vec)))
    
    names(r) <- var_name
    return(r)
  }
  
  # Create rasters for PAR variables
  mean_par_raster <- create_par_raster("mean_par_matrix.csv", "Mean_PAR")
  max_par_raster <- create_par_raster("max_par_matrix.csv", "Max_PAR")
  min_par_raster <- create_par_raster("min_par_matrix.csv", "Min_PAR")
  std_par_raster <- create_par_raster("std_par_matrix.csv", "Std_PAR")
  
  # Reproject PAR rasters
  mean_par_raster <- project(mean_par_raster, projected_crs)
  # mean_par_raster <- project(mean_par_raster, projected_crs, res = 25)
  max_par_raster <- project(max_par_raster, projected_crs)
  min_par_raster <- project(min_par_raster, projected_crs)
  std_par_raster <- project(std_par_raster, projected_crs)
  
  # Resample PAR rasters to match template
  mean_par_raster <- resample(mean_par_raster, template_raster, method = "bilinear")
  max_par_raster <- resample(max_par_raster, template_raster, method = "bilinear")
  min_par_raster <- resample(min_par_raster, template_raster, method = "bilinear")
  std_par_raster <- resample(std_par_raster, template_raster, method = "bilinear")
  
  ################################## Import Ocean Color Data ##################################
  
  ocean_color_summary <- read.csv(here("output", "ocean_color_summary_for_R.csv"))
  
  # Load ocean color coordinate information
  oc_lon_vec <- read.csv(here("output", "ocean_color_longitude.csv"))$longitude
  oc_lat_vec <- read.csv(here("output", "ocean_color_latitude.csv"))$latitude
  oc_grid_info <- read.csv(here("output", "ocean_color_grid_info.csv"))
  
  # Function to create ocean color rasters
  create_oc_raster <- function(filename, var_name) {
    cat(sprintf("Loading Ocean Color %s...\n", filename))
    matrix_data <- as.matrix(read.csv(here("output", filename), header = FALSE))
    
    # # Flip matrix vertically to match MATLAB orientation
    # matrix_data <- matrix_data[nrow(matrix_data):1, ]
    
    # Create SpatRaster using terra
    r <- rast(matrix_data,
              extent = ext(min(oc_lon_vec), max(oc_lon_vec), min(oc_lat_vec), max(oc_lat_vec)),
              crs = "EPSG:4326")
    
    names(r) <- var_name
    return(r)
  }
  
  ################################## Kd490 (Diffuse Attenuation Coefficient) ##################################
  
  # Create rasters for Kd490 variables  
  mean_kd490_raster <- create_oc_raster("mean_kd490_matrix.csv", "Mean_Kd490")
  max_kd490_raster <- create_oc_raster("max_kd490_matrix.csv", "Max_Kd490")
  min_kd490_raster <- create_oc_raster("min_kd490_matrix.csv", "Min_Kd490")
  std_kd490_raster <- create_oc_raster("std_kd490_matrix.csv", "Std_Kd490")
  
  # Reproject Kd490 rasters
  mean_kd490_raster <- project(mean_kd490_raster, projected_crs)
  max_kd490_raster <- project(max_kd490_raster, projected_crs)
  min_kd490_raster <- project(min_kd490_raster, projected_crs)
  std_kd490_raster <- project(std_kd490_raster, projected_crs)
  
  # Resample Kd490 rasters to match template
  mean_kd490_raster <- resample(mean_kd490_raster, template_raster, method = "bilinear")
  max_kd490_raster <- resample(max_kd490_raster, template_raster, method = "bilinear")
  min_kd490_raster <- resample(min_kd490_raster, template_raster, method = "bilinear")
  std_kd490_raster <- resample(std_kd490_raster, template_raster, method = "bilinear")
  
  ################################## Chlorophyll-a ##################################
  
  # Create rasters for Chlorophyll-a variables
  mean_chlor_a_raster <- create_oc_raster("mean_chlor_a_matrix.csv", "Mean_Chlor_a")
  max_chlor_a_raster <- create_oc_raster("max_chlor_a_matrix.csv", "Max_Chlor_a")
  min_chlor_a_raster <- create_oc_raster("min_chlor_a_matrix.csv", "Min_Chlor_a")
  std_chlor_a_raster <- create_oc_raster("std_chlor_a_matrix.csv", "Std_Chlor_a")
  
  # Reproject Chlorophyll-a rasters
  mean_chlor_a_raster <- project(mean_chlor_a_raster, projected_crs)
  max_chlor_a_raster <- project(max_chlor_a_raster, projected_crs)
  min_chlor_a_raster <- project(min_chlor_a_raster, projected_crs)
  std_chlor_a_raster <- project(std_chlor_a_raster, projected_crs)
  
  # Resample Chlorophyll-a rasters to match template
  mean_chlor_a_raster <- resample(mean_chlor_a_raster, template_raster, method = "bilinear")
  max_chlor_a_raster <- resample(max_chlor_a_raster, template_raster, method = "bilinear")
  min_chlor_a_raster <- resample(min_chlor_a_raster, template_raster, method = "bilinear")
  std_chlor_a_raster <- resample(std_chlor_a_raster, template_raster, method = "bilinear")
  
  ################################## Suspended Particulate Matter (SPM) ##################################
  
  # Create rasters for SPM variables
  mean_spm_raster <- create_oc_raster("mean_spm_matrix.csv", "Mean_SPM")
  max_spm_raster <- create_oc_raster("max_spm_matrix.csv", "Max_SPM")
  min_spm_raster <- create_oc_raster("min_spm_matrix.csv", "Min_SPM")
  std_spm_raster <- create_oc_raster("std_spm_matrix.csv", "Std_SPM")
  
  # Reproject SPM rasters
  mean_spm_raster <- project(mean_spm_raster, projected_crs)
  max_spm_raster <- project(max_spm_raster, projected_crs)
  min_spm_raster <- project(min_spm_raster, projected_crs)
  std_spm_raster <- project(std_spm_raster, projected_crs)
  
  # Resample SPM rasters to match template
  mean_spm_raster <- resample(mean_spm_raster, template_raster, method = "bilinear")
  max_spm_raster <- resample(max_spm_raster, template_raster, method = "bilinear")
  min_spm_raster <- resample(min_spm_raster, template_raster, method = "bilinear")
  std_spm_raster <- resample(std_spm_raster, template_raster, method = "bilinear")
  
  rm(template_raster)
  
  ################################## apply landmask ##################################
  
  # Create landmask from bathy_final 
  landmask <- !is.na(bathy_final)
  
  # Apply landmask to wave rasters
  mean_hsig_raster[landmask == 0] <- NA
  max_hsig_raster[landmask == 0] <- NA
  per_at_max_hsig_raster[landmask == 0] <- NA
  dir_at_max_hsig_raster[landmask == 0] <- NA
  dir_erddap_raster[landmask == 0] <- NA
  
  # Apply landmask to SST rasters
  mean_sst_raster[landmask == 0] <- NA
  max_sst_raster[landmask == 0] <- NA
  min_sst_raster[landmask == 0] <- NA
  std_sst_raster[landmask == 0] <- NA
  
  # Apply landmask to PAR rasters
  mean_par_raster[landmask == 0] <- NA
  max_par_raster[landmask == 0] <- NA
  min_par_raster[landmask == 0] <- NA
  std_par_raster[landmask == 0] <- NA
  
  # Apply landmask to ocean color rasters
  mean_kd490_raster[landmask == 0] <- NA
  max_kd490_raster[landmask == 0] <- NA
  min_kd490_raster[landmask == 0] <- NA
  std_kd490_raster[landmask == 0] <- NA
  
  mean_chlor_a_raster[landmask == 0] <- NA
  max_chlor_a_raster[landmask == 0] <- NA
  min_chlor_a_raster[landmask == 0] <- NA
  std_chlor_a_raster[landmask == 0] <- NA
  
  mean_spm_raster[landmask == 0] <- NA
  max_spm_raster[landmask == 0] <- NA
  min_spm_raster[landmask == 0] <- NA
  std_spm_raster[landmask == 0] <- NA
  
  #add SST range
  # NOTE - this should perhaps be done sooner upstream
  range_sst_raster = max_sst_raster - min_sst_raster
  
  # PAR range  
  range_par_raster <- max_par_raster - min_par_raster
  
  # Apply landmask to distance rasters
  dist_to_land_raster[landmask == 0] <- NA
  distance_to_deep_raster[landmask == 0] <- NA
  
  ################################## plots ##################################
  
  # Define plot extent options
  # plot_extents = ext(280000, 310000, 2010000, 2060000) #for investigating drops
  # plot_extents = ext(280000, 310000, 2000000, 2040000) #for investigating south of STT
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # plot_extents = ext(341000, 379000, 2057000, 2078000) # for investigating Anegada
  # plot_extents = ext(294000, 350000, 1950000, 1975000) #for investigating STX
  # plot_extents = ext(280000, 320000, 2000000, 2040000) #for investigating STT
  # plot_extents = ext(-30000, 0, 2000000, 2025000) #for investigating Mona Island
  # plot_extents = ext(20000, 40000, 2030000, 2045000) #for investigating Desecheo Island
  
  #plot SST
  plot(mean_sst_raster,
       col = cmocean("thermal")(100),
       ext = plot_extents)
  plot(max_sst_raster,
       col = cmocean("thermal")(100),
       ext = plot_extents)
  plot(min_sst_raster,
       col = cmocean("thermal")(100),
       ext = plot_extents)
  plot(std_sst_raster,
       col = cmocean("thermal")(100),
       ext = plot_extents)
  plot(range_sst_raster,
       col = cmocean("thermal")(100),
       ext = plot_extents)
  
  plot(mean_sst_raster,
       col = cmocean("thermal")(100))
  plot(max_sst_raster,
       col = cmocean("thermal")(100))
  plot(min_sst_raster,
       col = cmocean("thermal")(100))
  plot(std_sst_raster,
       col = cmocean("thermal")(100))
  plot(range_sst_raster,
       col = cmocean("thermal")(100))
  
  #plot waves
  plot(mean_hsig_raster,
       col = cmocean("amp")(100),
       ext = plot_extents)
  plot(max_hsig_raster,
       col = cmocean("amp")(100),
       ext = plot_extents)
  plot(dir_at_max_hsig_raster,
       col = cmocean("phase")(100),
       ext = plot_extents)
  plot(per_at_max_hsig_raster,
       col = cmocean("amp")(100),
       ext = plot_extents)
  plot(dir_erddap_raster,
       col = cmocean("phase")(100),
       ext = plot_extents)
  
  plot(mean_hsig_raster,
       col = cmocean("amp")(100))
  plot(max_hsig_raster,
       col = cmocean("amp")(100))
  plot(dir_at_max_hsig_raster,
       col = cmocean("phase")(100))
  plot(per_at_max_hsig_raster,
       col = cmocean("amp")(100))
  plot(dir_erddap_raster,
       col = cmocean("phase")(100))
  
  #plot PAR
  plot(mean_par_raster,
       col = cmocean("solar")(100),
       ext = plot_extents)
  plot(max_par_raster,
       col = cmocean("solar")(100),
       ext = plot_extents)
  plot(min_par_raster,
       col = cmocean("solar")(100),
       ext = plot_extents)
  plot(std_par_raster,
       col = cmocean("solar")(100),
       ext = plot_extents)
  plot(range_par_raster,
       col = cmocean("solar")(100),
       ext = plot_extents)
  
  plot(mean_par_raster,
       col = cmocean("solar")(100))
  plot(max_par_raster,
       col = cmocean("solar")(100))
  plot(min_par_raster,
       col = cmocean("solar")(100))
  plot(std_par_raster,
       col = cmocean("solar")(100))
  plot(range_par_raster,
       col = cmocean("solar")(100))
  
  #plot turbidity
  plot(mean_kd490_raster,
       col = cmocean("turbid")(100),
       ext = plot_extents)
  plot(max_kd490_raster,
       col = cmocean("turbid")(100),
       ext = plot_extents)
  plot(min_kd490_raster,
       col = cmocean("turbid")(100),
       ext = plot_extents)
  plot(std_kd490_raster,
       col = cmocean("turbid")(100),
       ext = plot_extents)
  # plot(range_kd490_raster,
  #      col = cmocean("turbid")(100),
  #      ext = plot_extents)
  
  plot(mean_kd490_raster,
       col = cmocean("turbid")(100))
  plot(max_kd490_raster,
       col = cmocean("turbid")(100))
  plot(min_kd490_raster,
       col = cmocean("turbid")(100))
  plot(std_kd490_raster,
       col = cmocean("turbid")(100))
  # plot(range_kd490_raster,
  #      col = cmocean("turbid")(100))
  
  mean_kd490_clamp <- clamp(mean_kd490_raster, lower = 0, upper = 0.5)
  plot(mean_kd490_clamp,
       col = cmocean("turbid")(100),
       ext = plot_extents)
  
  plot(mean_kd490_clamp,
       col = cmocean("turbid")(100))
  
  #plot chl-a
  plot(mean_chlor_a_raster,
       col = cmocean("algae")(100),
       ext = plot_extents)
  plot(max_chlor_a_raster,
       col = cmocean("algae")(100),
       ext = plot_extents)
  plot(min_chlor_a_raster,
       col = cmocean("algae")(100),
       ext = plot_extents)
  plot(std_chlor_a_raster,
       col = cmocean("algae")(100),
       ext = plot_extents)
  # plot(range_chlor_a_raster,
  #      col = cmocean("algae")(100),
  #      ext = plot_extents)
  
  plot(mean_chlor_a_raster,
       col = cmocean("algae")(100))
  plot(max_chlor_a_raster,
       col = cmocean("algae")(100))
  plot(min_chlor_a_raster,
       col = cmocean("algae")(100))
  plot(std_chlor_a_raster,
       col = cmocean("algae")(100))
  # plot(range_chlor_a_raster,
  #      col = cmocean("algae")(100))
  
  #plot particulates
  plot(mean_spm_raster,
       col = cmocean("matter")(100),
       ext = plot_extents)
  plot(max_spm_raster,
       col = cmocean("matter")(100),
       ext = plot_extents)
  plot(min_spm_raster,
       col = cmocean("matter")(100),
       ext = plot_extents)
  plot(std_spm_raster,
       col = cmocean("matter")(100),
       ext = plot_extents)
  # plot(range_spm_raster,
  #      col = cmocean("matter")(100),
  #      ext = plot_extents)
  
  plot(mean_spm_raster,
       col = cmocean("matter")(100))
  plot(max_spm_raster,
       col = cmocean("matter")(100))
  plot(min_spm_raster,
       col = cmocean("matter")(100))
  plot(std_spm_raster,
       col = cmocean("matter")(100))
  # plot(range_spm_raster,
  #      col = cmocean("matter")(100))
  
  mean_spm_clamp <- clamp(mean_spm_raster, lower = 0, upper = 1)
  plot(mean_spm_clamp,
       col = cmocean("matter")(100),
       ext = plot_extents)
  
  plot(mean_spm_clamp,
       col = cmocean("matter")(100))
  
  
  #plot distances
  plot(dist_to_land_raster,
       col = cmocean("solar")(100),
       ext = plot_extents)
  plot(distance_to_deep_raster,
       col = cmocean("solar")(100),
       ext = plot_extents)

  plot(dist_to_land_raster,
       col = cmocean("solar")(100))
  plot(distance_to_deep_raster,
       col = cmocean("solar")(100))
  
  ################################## Save objects/workspace ##################################
  
  # save_new_objects("output/output_calculate_ocean_rasters", existing_objects)
  