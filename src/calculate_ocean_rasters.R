  
  # .rs.restartR(clean = TRUE)
  
  library(sf)
  library(here)
  library(terra) 

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
  
  # # VERSION where SWAN data was pulled from ERDDAP (old)
  # #
  # # wave_summary <- read.csv(here("output", "swan_wave_summary_for_R.csv"))
  # 
  # # 2. Load coordinate information
  # lon_vec <- read.csv(here("output", "swan_longitude.csv"))$longitude
  # lat_vec <- read.csv(here("output", "swan_latitude.csv"))$latitude
  # grid_info <- read.csv(here("output", "swan_grid_info.csv"))
  # 
  # # Get grid dimensions
  # nlon <- length(lon_vec)
  # nlat <- length(lat_vec)
  # lon_range <- range(lon_vec)
  # lat_range <- range(lat_vec)
  # 
  # # Function to load matrix and create raster
  # create_raster <- function(filename, var_name) {
  #   cat(sprintf("Loading %s...\n", filename))
  #   matrix_data <- as.matrix(read.csv(here("output", filename), header = FALSE))
  #   
  #   # Flip matrix vertically to match MATLAB orientation
  #   matrix_data <- matrix_data[nrow(matrix_data):1, ]
  #   
  #   # Create SpatRaster using terra
  #   r <- rast(matrix_data,
  #             extent = ext(min(lon_vec), max(lon_vec), min(lat_vec), max(lat_vec)),
  #             crs = "EPSG:4326")
  #   
  #   names(r) <- var_name
  #   return(r)
  # }
  # 
  # # Create rasters for key variables
  # mean_hsig_raster <- create_raster("mean_hsig_matrix.csv", "Mean_Hsig")
  # max_hsig_raster <- create_raster("max_hsig_matrix.csv", "Max_Hsig")
  # mean_hswell_raster <- create_raster("mean_hswell_matrix.csv", "Mean_Hswell")
  # mean_dir_raster <- create_raster("mean_dir_matrix.csv", "Mean_Direction")
  # mean_per_raster <- create_raster("mean_per_matrix.csv", "Mean_Period")
  # 
  # projected_crs = crs(bathy_final)
  # 
  # # Reproject to match bathy_final CRS first
  # mean_hsig_raster <- project(mean_hsig_raster, projected_crs)
  # max_hsig_raster <- project(max_hsig_raster, projected_crs)
  # mean_hswell_raster <- project(mean_hswell_raster, projected_crs)
  # mean_dir_raster <- project(mean_dir_raster, projected_crs)
  # mean_per_raster <- project(mean_per_raster, projected_crs)
  
  # VERSION where SWAN data was provided directly from Miguel Canals
  # 
  SWAN_target_res <- 200
  
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
    
    # Memory safety check
    if (n_cells > 25e6) {
      stop(sprintf("Grid too large (%.1f million cells). Use resolution > 100m.", n_cells/1e6))
    }
    
    if (n_cells > 10e6) {
      warning(sprintf("Large grid (%.1f million cells) - may be slow", n_cells/1e6))
    }
    
    # Rasterize wave points to template
    result <- rasterize(wave_points_valid, template_raster, 
                        field = var_name, 
                        fun = mean)
    
    # Gap filling with focal operations
    cat("Gap filling...\n")
    result <- focal(result, w = 3, fun = mean, na.policy = "only", na.rm = TRUE)
    
    # Additional gap filling pass for smoother results
    result <- focal(result, w = 3, fun = mean, na.policy = "only", na.rm = TRUE)
    
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
  
  
  
  
  
  
  # STOPPING POINT - 24 August 2025
  #   - took a stab at "max" BOV - could maybe use work. not sure if it makes sense that the
  #       deeper shelf areas alone have high BOV
  
  # Constants
  g <- 9.81
  
  # Make depth positive
  depth_raster <- abs(bathy_final)
  
  # Initial wavelength guess using deep water approximation
  L <- g * per_at_max_hsig_raster^2 / (2 * pi)
  
  # Iteratively solve dispersion relation (Newton-Raphson, 5 iterations should be enough)
  cat("Solving dispersion relation...\n")
  pb <- txtProgressBar(min = 0, max = 5, style = 3)
  
  for (i in 1:5) {
    k <- 2 * pi / L
    omega_sq <- (2 * pi / per_at_max_hsig_raster)^2
    f <- omega_sq - g * k * tanh(k * depth_raster)
    df_dL <- g * (2 * pi / L^2) * (tanh(2 * pi * depth_raster / L) - 
                                     (2 * pi * depth_raster / L) * (1 / cosh(2 * pi * depth_raster / L))^2)
    L <- L + f / df_dL
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Calculate BOV
  k_h <- 2 * pi * depth_raster / L
  cosh_term <- cosh(k_h)
  max_bov_raster <- (max_hsig_raster * g * per_at_max_hsig_raster) / (2 * L * cosh_term)
  
  # Set name
  names(max_bov_raster) <- "max_BOV"
  
  # Summary
  cat("BOV calculation complete!\n")
  cat(sprintf("BOV range: %.3f to %.3f m/s\n", 
              minmax(max_bov_raster)[1], minmax(max_bov_raster)[2]))
  
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
  
  
  
  ##### plots #####
  
  
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
  
  plot(mean_hsig_raster,
       col = cmocean("amp")(100))
  plot(max_hsig_raster,
       col = cmocean("amp")(100))
  plot(dir_at_max_hsig_raster,
       col = cmocean("phase")(100))
  plot(per_at_max_hsig_raster,
       col = cmocean("amp")(100))
  
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
  plot(range_kd490_raster,
       col = cmocean("turbid")(100),
       ext = plot_extents)
  
  plot(mean_kd490_raster,
       col = cmocean("turbid")(100))
  plot(max_kd490_raster,
       col = cmocean("turbid")(100))
  plot(min_kd490_raster,
       col = cmocean("turbid")(100))
  plot(std_kd490_raster,
       col = cmocean("turbid")(100))
  plot(range_kd490_raster,
       col = cmocean("turbid")(100))
  
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
  plot(range_chlor_a_raster,
       col = cmocean("algae")(100),
       ext = plot_extents)
  
  plot(mean_chlor_a_raster,
       col = cmocean("algae")(100))
  plot(max_chlor_a_raster,
       col = cmocean("algae")(100))
  plot(min_chlor_a_raster,
       col = cmocean("algae")(100))
  plot(std_chlor_a_raster,
       col = cmocean("algae")(100))
  plot(range_chlor_a_raster,
       col = cmocean("algae")(100))
  
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
  plot(range_spm_raster,
       col = cmocean("matter")(100),
       ext = plot_extents)
  
  plot(mean_spm_raster,
       col = cmocean("matter")(100))
  plot(max_spm_raster,
       col = cmocean("matter")(100))
  plot(min_spm_raster,
       col = cmocean("matter")(100))
  plot(std_spm_raster,
       col = cmocean("matter")(100))
  plot(range_spm_raster,
       col = cmocean("matter")(100))
  
  # STOPPING POINT - 23 August 2025
  #   - Okay, still need to put together BOVs, distance from mesophotic isobath,
  #       and distance from shoreline. other than that, looking really good.
  #   - strongly consider adding year as a variable, given that Viehman 2025 / NCRMP
  #       clearly demonstrate a large decline in coral cover from 2013-2017.
  #       I would not like to make this a full-blown spatiotemporal SDM, as I do
  #       not think that is appropriate for the project scope - but should explore
  #       to what degree adding in "time" is possible
  #   - will ideally incorporate mean BOV and mean wave direction from Canals
  
  mean_spm_clamp <- clamp(mean_spm_raster, lower = 0, upper = 1)
  plot(mean_spm_clamp,
       col = cmocean("matter")(100),
       ext = plot_extents)
  
  plot(mean_spm_clamp,
       col = cmocean("matter")(100))
  
  
  ################################## Save objects/workspace ##################################
  
  # save_new_objects("output/output_calculate_ocean_rasters", existing_objects)
  