  
  # .rs.restartR(clean = TRUE)
  
  library(sf)
  library(here)
  library(terra) 
  library(tidyterra)
  library(ggplot2)
  library(tmap)
  library(rayshader) #this requires installation of XQuartz on MacOS, and possibly OpenGL if it isn't installed
  library(scico)
  library(RColorBrewer)
  library(progress)
  library(spatialEco)
  library(MultiscaleDTM)
  library(cmocean)
  
  source(here("src/functions.R"))
  
  ################################## setup ##################################
  
  # #load initial set of objects
  # load_spat_objects(directory = 'output/output_import_merge_rasters_higher-res/')
  # load(here('output', 'output_import_merge_rasters_higher-res/import_merge_rasters_workspace.RData'))
  # 
  # source(here("src/functions.R"))
  
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
  
  existing_objects <- ls(envir = .GlobalEnv)
  
  # ################################## Pretty bathy plot ##################################
  # 
  # # color_palette <- scico(100, palette = "hawaii", direction = -1) #scico palette
  # color_palette <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(100) #YlGnBu palette, as used by me in QGIS. color ramp reversed, and continuous colors interpolated
  # 
  # # Plot the raster
  # plot(bathy_final,
  #      col = color_palette,
  #      # zlim = c(-50, 0),
  #      main = "Bathymetry (50m Resolution)",
  #      legend = TRUE)
  # 
  # # # Reproject to NAD83 (geographic)
  # # # NOTE - 4 JULY 2025
  # # #   - will need to go back upstream and project to geographic earlier,
  # # #       to avoid loss of resolution
  # # #   - 9 July 2025: I think I am just keeping everything in NCEI's native UTM actually
  # # bathy_merged_geo <- project(bathy_merged3_crm_reefdepth, "EPSG:4269")
  # # 
  # # #note the slight shift in orientation
  # # # NOTE - after looking at this, I think maybe we SHOULD project earlier in the pipeline. mesophotic ridge loses
  # # #         significant detail
  # # plot(bathy_merged_geo,
  # #      col = color_palette,
  # #      # zlim = c(-50, 0), #redundant
  # #      main = "Bathymetry (50m Resolution)",
  # #      legend = TRUE)
  # 
  # #plot just north of STT to verify clamping function above worked correctly
  # #
  # # plot_extents = ext(280000, 310000, 2000000, 2040000) #for investigating south of STT
  # # plot_extents = ext(260000, 290000, 2000000, 2040000) #for investigating MCD
  # # plot_extents = ext(305000, 330000, 2020000, 2035000) #for investigating STJ
  # # plot_extents = ext(317000, 350000, 2030000, 2050000) #for investigating Tortola
  # # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # # plot_extents = ext(341000, 379000, 2057000, 2078000) # for investigating Anegada
  # # plot_extents = ext(300000, 340000, 1940000, 1980000) #for investigating St Croix
  # # plot_extents = ext(320000, 322000, 1962000, 1964000) #for investigating Altona Lagoon, STX
  # plot_extents = ext(279000, 310000, 2010000, 2050000) #for investigating St Thomas
  # # plot_extents = ext(240000, 275000, 2000000, 2040000) #for investigating Culebra
  # # plot_extents = ext(120000, 220000, 2020000, 2060000) #for investigating northern PR
  # 
  # extent_area = plot_extents
  # 
  # # Crop the raster to the defined extent
  # bathy_cropped <- crop(bathy_final, extent_area)
  # 
  # # Define a color palette for the plot
  # color_palette <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(100)
  # 
  # # Plot the cropped raster
  # plot(bathy_cropped,
  #      col = color_palette,
  #      main = "Bathymetry North of St. Thomas Up to Latitude 19°",
  #      legend = TRUE)
  # 
  # #verify that the clamp worked and introduced NAs
  # values <- values(bathy_cropped)
  # deeper_than_50m <- any(values < -50, na.rm = TRUE)
  # 
  # if (deeper_than_50m) {
  #   cat("There are values deeper than 50 meters in the raster.")
  # } else {
  #   cat("All values are within the specified range (i.e., shallower than 50 meters).")
  # }
  # 
  # na_values <- any(is.na(values))
  # if (na_values) {
  #   cat("There are NA values in the raster.")
  # } else {
  #   cat("There are no NA values in the raster.")
  # }
  # 
  # # Create a mask raster where NA values are set to 1 and all other values are set to 0
  # na_mask <- is.na(bathy_cropped)
  # 
  # # Define a color palette for the mask
  # mask_palette <- c("transparent", "red")  # Transparent for non-NA, red for NA
  # 
  # # Plot the NA mask raster
  # plot(na_mask,
  #      col = mask_palette,
  #      main = "NA Mask of Bathymetry Raster",
  #      legend = FALSE)  # No legend for binary mask
  # 
  # # Create a binary raster where values shallower than 50 meters are set to 1 and others are set to 0
  # deep_mask <- bathy_cropped
  # deep_mask[] <- ifelse(values(bathy_cropped) < -50, 1, 0)
  # 
  # # Define a color palette for the binary raster
  # deep_palette <- c("transparent", "blue")  # Transparent for shallow areas, blue for deep areas. NOTE - maybe reverse?
  # 
  # # Plot the binary raster for values deeper than 50 meters
  # plot(deep_mask,
  #      col = deep_palette,
  #      main = "Binary Mask for Depths Deeper Than 50 Meters",
  #      legend = FALSE)  # No legend for binary mask
  # 
  # bathy_value <- extract(bathy_cropped, cbind(-64.95, 18.45))
  # 
  # # #ggplot method, with tidyterra. this maybe didn't work because of the size of the rasters, so I went with tmap & rayshader below
  # #
  # # # # Convert raster to a data frame for ggplot2
  # # # bathy_raster <- rast(here("output", "bathy_50m.tif"))
  # # # bathy_df <- as_tibble(bathy_raster, xy = TRUE)
  # # bathy_df <- as_tibble(merged_bathy, xy = TRUE)
  # #
  # # # Plot raster data
  # # ggplot(bathy_df) +
  # #   geom_raster(aes(x = x, y = y, fill = layer)) +
  # #   scale_fill_viridis_c(option = "D", limits = c(-50, 0)) + # Adjust limits to desired depth range
  # #   labs(
  # #     title = "Bathymetry",
  # #     fill = "Depth (meters)"
  # #   ) +
  # #   coord_fixed() + # Fix aspect ratio to make sure x and y are scaled equally
  # #   theme_minimal() +
  # #   theme(
  # #     axis.title = element_blank(),
  # #     axis.text = element_text(size = 8),
  # #     legend.title = element_text(size = 10),
  # #     legend.text = element_text(size = 8)
  # #   )
  # #
  # # ggplot(bathy_df) +
  # #   geom_raster(aes(x = x, y = y))
  # 
  # #tmap method
  # tm_shape(bathy_final) +
  #   tm_raster(style = 'cont', palette = color_palette) #this is outdated as of tmap v4, but works for continuous gradient
  # # tm_shape(bathy_merged_50m) +
  # #   tm_raster() +
  # #   tm_scale_continuous(values = color_palette) #this doesn't seem right, but makes a discrete contoured depth gradient
  # 
  # # bathy_map = tm_shape(merged_bathy) +
  # #   tm_raster(style = 'cont')
  # 
  # # output_file = file.path(here('output'), "bathy_50m.svg")
  # # tmap_save(bathy_map, filename=output_file, height=8.5, width=11, units="in", dpi=300)
  # # output_file = file.path(here('output'), "bathy_50m.pdf")
  # # tmap_save(bathy_map, filename=output_file, height=8.5, width=11, units="in", dpi=300)
  # # output_file = file.path(here('output'), "bathy_50m.jpeg")
  # # tmap_save(bathy_map, filename=output_file, height=8.5, width=11, units="in", dpi=300)
  # # output_file = file.path(here('output'), "bathy_50m.tiff")
  # # tmap_save(bathy_map, filename=output_file, height=8.5, width=11, units="in", dpi=300)
  # 
  # # # Save the merged raster
  # # output_dir <- here("output")
  # # output_file <- file.path(output_dir, "bathy_50m.tif")
  # # writeRaster(merged_bathy, filename = output_file, overwrite = TRUE)
  # 
  # #rayshader method
  # # in 2D
  # # elmat = raster_to_matrix(merged_bathy)
  # # elmat %>%
  # #   sphere_shade(texture = "desert") %>%
  # #   plot_map()
  # raster_to_matrix(bathy_final) |> height_shade() |> plot_map()
  # 
  # # # in 3D
  # # elmat %>%
  # #   # sphere_shade(texture = "desert") %>%
  # #   # add_water(detect_water(elmat), color = "desert") %>%
  # #   # add_shadow(ray_shade(elmat, zscale = 3), 0.5) %>%
  # #   add_shadow(ambient_shade(elmat), 0) %>%
  # #   plot_3d(elmat, zscale = 10, fov = 0, theta = 135, zoom = 0.75, phi = 45, windowsize = c(1000, 800))
  # # Sys.sleep(0.2)
  # # render_snapshot()
  # 
  ################################## Test resolution/artifacts ##################################
  
  # Calculate terrain attributes. may consider Benthic Terrain Modeler (but need Arc & potentially arcgisbinding package in R) or MultiscaleDTM package (https://cran.r-project.org/web/packages/MultiscaleDTM/readme/README.html)
  #
  slope_terra <- terrain(bathy_final, v = "slope", unit = 'degrees')
  slopeofslope_terra = terrain(slope_terra, v = "slope", unit = 'degrees')
  aspect_terra <- terrain(bathy_final, v = "aspect", unit = 'degrees')
  TPI_terra = terrain(bathy_final, v = 'TPI') #Wilson 2007 bathymetric-friendly method
  # flowdir <- terrain(bathy_final, v = "flowdir")
  roughness = terrain(bathy_final, v = "roughness")
  TRI = terrain(bathy_final, v = 'TRI')
  
  # NOTE - running 'totalcurv' too here would be ideal since it might matter. but can crash R, so return later
  # TPI_spatialeco = tpi(bathy_final)
  VRM <- vrm(bathy_final)
  # # TRI_spatialeco <- tri(bathy_final, s = 3)
  # totalcurv = curvature(bathy_final, type = 'total')
  # # planformcurv = curvature(bathy_final, type = 'planform') #takes too long to run
  # # profilecurv = curvature(bathy_final, type = 'profile') #takes too long to run
  
  # consider geodiv! (can get fractal dimension)
  
  # slope_multiscale <- SlpAsp(r = bathy_final, w = 3, metrics = c("slope")) #na.rm = TRUE could be used
  # slopeofslope_multiscale <- SlpAsp(r = slope_multiscale, w = 3, metrics = c("slope"))
  # aspect_multiscale <- SlpAsp(r = bathy_final, w = 3, metrics = c("aspect"))
  # TPI_multiscale <- TPI(bathy_final, w = c(3, 3))
  # BPI_multiscale <- BPI(bathy_final, w = c(1, 3))  # BPI is similar to TPI for bathymetry
  # DMV <- DMV(bathy_final, w = c(3, 3))
  # RelPos <- RelPos(bathy_final, w = c(5, 5))
  SAPA <- SAPA(bathy_final, w = c(3, 3))
  VRM_multiscale <- VRM(bathy_final, w = c(3, 3))
  # surfarea = SurfaceArea(bathy_final)
  # maxcurv_multiscale = Qfit(r = bathy_final, w = 3, metrics = 'maxc')
  meancurv_multiscale = Qfit(r = bathy_final, w = 3, metrics = 'meanc')
  planformcurv_multiscale = Qfit(r = bathy_final, w = 3, metrics = 'planc')
  profilecurv_multiscale = Qfit(r = bathy_final, w = 3, metrics = 'profc')
  # depth = Qfit(r = bathy_final, metrics = 'elev')
  
  # NOTE - I think slope, slope of slope, aspect, and TPI, roughness, TRI, and
  #         VRM are the variables I am likely moving forward with
  
  ################################## Import waves ##################################
  
  wave_summary <- read.csv(here("output", "swan_wave_summary_for_R.csv"))
  
  # 2. Load coordinate information
  lon_vec <- read.csv(here("output", "swan_longitude.csv"))$longitude
  lat_vec <- read.csv(here("output", "swan_latitude.csv"))$latitude
  grid_info <- read.csv(here("output", "swan_grid_info.csv"))
  
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
  
  # Create rasters for key variables
  mean_hsig_raster <- create_raster("mean_hsig_matrix.csv", "Mean_Hsig")
  max_hsig_raster <- create_raster("max_hsig_matrix.csv", "Max_Hsig")
  mean_hswell_raster <- create_raster("mean_hswell_matrix.csv", "Mean_Hswell")
  mean_dir_raster <- create_raster("mean_dir_matrix.csv", "Mean_Direction")
  mean_per_raster <- create_raster("mean_per_matrix.csv", "Mean_Period")
  
  projected_crs = crs(bathy_final)
  
  # Reproject to match bathy_final CRS first
  mean_hsig_raster <- project(mean_hsig_raster, projected_crs)
  max_hsig_raster <- project(max_hsig_raster, projected_crs)
  mean_hswell_raster <- project(mean_hswell_raster, projected_crs)
  mean_dir_raster <- project(mean_dir_raster, projected_crs)
  mean_per_raster <- project(mean_per_raster, projected_crs)
  
  # Get extent and properties from bathy_final
  bathy_ext <- ext(bathy_final)
  bathy_crs <- crs(bathy_final)
  bathy_res <- res(bathy_final)
  
  # Create template raster matching bathy_final exactly
  template_raster <- rast(bathy_ext, resolution = bathy_res, crs = bathy_crs)
  
  # Resample wave rasters to match bathy_final
  mean_hsig_raster <- resample(mean_hsig_raster, template_raster, method = "bilinear")
  max_hsig_raster <- resample(max_hsig_raster, template_raster, method = "bilinear")
  mean_hswell_raster <- resample(mean_hswell_raster, template_raster, method = "bilinear")
  mean_dir_raster <- resample(mean_dir_raster, template_raster, method = "bilinear")
  mean_per_raster <- resample(mean_per_raster, template_raster, method = "bilinear")
  
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
    
    # Flip matrix vertically to match MATLAB orientation
    matrix_data <- matrix_data[nrow(matrix_data):1, ]
    
    # Create SpatRaster using terra
    r <- rast(matrix_data,
              extent = ext(min(par_lon_vec), max(par_lon_vec), min(par_lat_vec), max(par_lat_vec)),
              crs = "EPSG:4326")
    
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
    
    # Flip matrix vertically to match MATLAB orientation
    matrix_data <- matrix_data[nrow(matrix_data):1, ]
    
    # Create SpatRaster using terra
    r <- rast(matrix_data,
              extent = ext(min(oc_lon_vec), max(oc_lon_vec), min(oc_lat_vec), max(oc_lat_vec)),
              crs = "EPSG:4326")
    
    names(r) <- var_name
    return(r)
  }
  
  ################################## Kd490 (Diffuse Attenuation Coefficient) ##################################
  
  # Check if Kd490 data exists in the summary
  if("mean_kd490" %in% names(ocean_color_summary)) {
    
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
    
    cat("Kd490 rasters created successfully.\n")
  } else {
    cat("Kd490 data not found in ocean color summary.\n")
  }
  
  ################################## Chlorophyll-a ##################################
  
  # Check if Chlorophyll-a data exists in the summary
  if("mean_chlor_a" %in% names(ocean_color_summary)) {
    
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
    
    cat("Chlorophyll-a rasters created successfully.\n")
  } else {
    cat("Chlorophyll-a data not found in ocean color summary.\n")
  }
  
  ################################## Suspended Particulate Matter (SPM) ##################################
  
  # Check if SPM data exists in the summary
  if("mean_spm" %in% names(ocean_color_summary)) {
    
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
    
    cat("SPM rasters created successfully.\n")
  } else {
    cat("SPM data not found in ocean color summary.\n")
  }
  
  rm(template_raster)
  
  ################################## apply landmask ##################################
  
  # Create landmask from bathy_final 
  landmask <- !is.na(bathy_final)
  
  # Apply landmask to wave rasters
  mean_hsig_raster[landmask == 0] <- NA
  max_hsig_raster[landmask == 0] <- NA
  mean_hswell_raster[landmask == 0] <- NA
  mean_dir_raster[landmask == 0] <- NA
  mean_per_raster[landmask == 0] <- NA
  
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
  
  # Apply landmask to Kd490 rasters (if they exist)
  if(exists("mean_kd490_raster")) {
    mean_kd490_raster[landmask == 0] <- NA
    max_kd490_raster[landmask == 0] <- NA
    min_kd490_raster[landmask == 0] <- NA
    std_kd490_raster[landmask == 0] <- NA
  }
  
  # Apply landmask to Chlorophyll-a rasters (if they exist)
  if(exists("mean_chlor_a_raster")) {
    mean_chlor_a_raster[landmask == 0] <- NA
    max_chlor_a_raster[landmask == 0] <- NA
    min_chlor_a_raster[landmask == 0] <- NA
    std_chlor_a_raster[landmask == 0] <- NA
  }
  
  # Apply landmask to SPM rasters (if they exist)
  if(exists("mean_spm_raster")) {
    mean_spm_raster[landmask == 0] <- NA
    max_spm_raster[landmask == 0] <- NA
    min_spm_raster[landmask == 0] <- NA
    std_spm_raster[landmask == 0] <- NA
  }
  
  #add SST range
  # NOTE - this should perhaps be done sooner upstream
  range_sst_raster = max_sst_raster - min_sst_raster
  
  # PAR range  
  range_par_raster <- max_par_raster - min_par_raster
  
  # Kd490 range (if exists)
  if(exists("mean_kd490_raster")) {
    range_kd490_raster <- max_kd490_raster - min_kd490_raster
  }
  
  # Chlorophyll-a range (if exists)
  if(exists("mean_chlor_a_raster")) {
    range_chlor_a_raster <- max_chlor_a_raster - min_chlor_a_raster
  }
  
  # SPM range (if exists)
  if(exists("mean_spm_raster")) {
    range_spm_raster <- max_spm_raster - min_spm_raster
  }
  
  ################################## plots ##################################
  
  # Define plot extent options
  # plot_extents = ext(280000, 310000, 2010000, 2060000) #for investigating drops
  # plot_extents = ext(280000, 310000, 2000000, 2040000) #for investigating south of STT
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # plot_extents = ext(341000, 379000, 2057000, 2078000) # for investigating Anegada
  # plot_extents = ext(294000, 350000, 1950000, 1975000) #for investigating St Croix
  # plot_extents = ext(280000, 320000, 2000000, 2040000) #for investigating St Thomas
  # plot_extents = ext(240000, 280000, 2000000, 2040000) #for investigating Mona Island

  # Slope
  plot(slope_terra,
       col = cmocean("deep")(100),
       ext = plot_extents,
       main = "Slope - Terra Raster with cmocean deep")

  # Slope of Slope
  plot(slopeofslope_terra,
       col = cmocean("deep")(100),
       ext = plot_extents,
       main = "Slope of Slope - Terra Raster with cmocean deep")

  # Roughness
  plot(roughness,
       col = cmocean("deep")(100),
       ext = plot_extents,
       main = "Roughness - Terra Raster with cmocean deep")

  # VRM
  VRM_clamp <- clamp(VRM, lower = -0.01, upper = 0.01)
  plot(VRM_clamp,
  # plot(VRM,
       col = cmocean("balance")(100),
       ext = plot_extents,
       main = "VRM - Terra Raster with cmocean balance")

  # # TPI spatialeco
  # TPI_spatialeco_clamp <- clamp(TPI_spatialeco, lower = -5, upper = 5)
  # # plot(TPI_spatialeco_clamp,
  # plot(TPI_spatialeco,
  #      col = cmocean("balance")(100),
  #      ext = plot_extents,
  #      main = "TPI spatialeco - Terra Raster with cmocean balance")
  #
  # # TPI multiscale
  # TPI_multiscale <- clamp(TPI_multiscale, lower = -5, upper = 5)
  # plot(TPI_multiscale,
  #      col = cmocean("balance")(100),
  #      ext = plot_extents,
  #      main = "TPI multiscale - Terra Raster with cmocean balance")

  # TPI terra
  TPI_terra_clamp <- clamp(TPI_terra, lower = -40, upper = 40)
  # plot(TPI_terra,
  plot(TPI_terra_clamp,
       col = cmocean("balance")(100),
       # ext = plot_extents,
       main = "TPI terra - Terra Raster with cmocean balance")

  # # BPI multiscale
  # BPI_multiscale_clamp <- clamp(BPI_multiscale, lower = -10, upper = 10)
  # plot(BPI_multiscale,
  # # plot(BPI_multiscale_clamp,
  #      col = cmocean("balance")(100),
  #      ext = plot_extents,
  #      main = "BPI multiscale - Terra Raster with cmocean balance")
  #
  # # DMV
  # DMV <- clamp(DMV, lower = -5, upper = 5)
  # plot(DMV,
  #      col = cmocean("balance")(100),
  #      ext = plot_extents,
  #      main = "DMV - Terra Raster with cmocean balance")

  # Total curvature (commented out)
  # totalcurv <- clamp(totalcurv, lower = -0.006, upper = 0.006)
  # plot(totalcurv,
  #      col = cmocean("balance")(100),
  #      ext = plot_extents,
  #      main = "Total Curvature - Terra Raster with cmocean balance")

  # Planform curvature multiscale
  planformcurv_multiscale_clamp <- clamp(planformcurv_multiscale, lower = -0.005, upper = 0.005)
  # plot(planformcurv_multiscale,
  plot(planformcurv_multiscale_clamp,
       col = cmocean("balance")(100),
       ext = plot_extents,
       main = "Planform Curvature multiscale - Terra Raster with cmocean balance")

  # Profile curvature multiscale
  profilecurv_multiscale <- clamp(profilecurv_multiscale, lower = -0.01, upper = 0.01)
  plot(profilecurv_multiscale,
       col = cmocean("balance")(100),
       ext = plot_extents,
       main = "Profile Curvature multiscale - Terra Raster with cmocean balance")

  # Mean curvature multiscale
  meancurv_multiscale_clamp <- clamp(meancurv_multiscale, lower = -0.004, upper = 0.004)
  # plot(meancurv_multiscale,
  plot(meancurv_multiscale_clamp,
       col = cmocean("balance")(100),
       ext = plot_extents,
       main = "Mean Curvature multiscale - Terra Raster with cmocean balance")

  # # Max curvature multiscale
  # maxcurv_multiscale_clamp <- clamp(maxcurv_multiscale, lower = -0.002, upper = 0.002)
  # # plot(maxcurv_multiscale,
  # plot(maxcurv_multiscale_clamp,
  #      col = cmocean("balance")(100),
  #      ext = plot_extents,
  #      main = "Max Curvature multiscale - Terra Raster with cmocean balance")

  # Bathymetry
  depth_OG_clamp <- clamp(bathy_final, lower = -50, upper = 0)
  # plot(bathy_final,
  plot(depth_OG_clamp,
       col = rev(cmocean("deep")(100)),  # Reversed color scheme
       ext = plot_extents,
       main = "Bathymetry - Terra Raster with cmocean deep (reversed)")
  
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
  
  plot(mean_hsig_raster,
       col = cmocean("amp")(100),
       ext = plot_extents)
  plot(max_hsig_raster,
       col = cmocean("amp")(100),
       ext = plot_extents)
  plot(mean_hswell_raster,
       col = cmocean("amp")(100),
       ext = plot_extents)
  plot(mean_dir_raster,
       col = cmocean("phase")(100),
       ext = plot_extents)
  plot(mean_per_raster,
       col = cmocean("amp")(100),
       ext = plot_extents)
  
  plot(mean_hsig_raster,
       col = cmocean("amp")(100))
  plot(max_hsig_raster,
       col = cmocean("amp")(100))
  plot(mean_hswell_raster,
       col = cmocean("amp")(100))
  plot(mean_dir_raster,
       col = cmocean("phase")(100))
  plot(mean_per_raster,
       col = cmocean("amp")(100))
  

  # NOTE / STOPPING POINT (15 Aug 2025) - something very wrong with PAR in northern PR. hopefully fixable!
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
  
  ################################## Kd490 Plots ##################################
  
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

  ################################## Chlorophyll-a Plots ##################################
  
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

  ################################## SPM Plots ##################################
  
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

  # NOTES - 18 May 2025
  #   - could also take the current 650-m resolution grid and slap a full grid over it, but with grid squares
  #       extending out to depths <50 m. This could be useful if we care about lining things up with what we already
  #       know "works" in the CMS. I think the part about lining up with NCRMP data matters less, because it's
  #       always a mess trying to get PR & USVI grids to be happy with each other anyways. when producing the GAMs,
  #       will just need to take *average* of outputted coral cover & tissue SA, rather than *sum*, within each
  #       650 x 650 m grid square
  #   - SO, I think the game plan is to ideally work with existing 650 m habitat grid, and find a way to *add* to it
  #       without adding back squares which were previously deleted because of bad CMS exit codes
  #   - ACTUALLY, thinking about this more I think we should just use the existing grid for CMS tests to confirm
  #       CMS is "good", then move forward with making a new habitat grid that is completely independent of what was
  #       done in the past. will need to nail down a routine to properly handle and whittle down the grid to remove
  #       squares with bad CMS exit codes
  #   - SO, all I need to do right now is 1.) move on with getting the CMS in order, and 2.) preparing a new, coarse
  #       grid use once we know the CMS works at all



  ################################## Save objects/workspace ##################################
  
  # #updated way to handle saving of new objects
  # save_new_objects("output/output_calculate_bathy_rasters", existing_objects)
    