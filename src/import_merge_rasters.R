  
  # .rs.restartR(clean = TRUE)
  
  library(here)
  library(terra)
  library(sf)
  
  # NOTE - loading & use of 'terra' with reading of geodatabase rasters requires (at least for me on M1 Macbook
  #         in August 2024) homebrew installation of proj and gdal (https://github.com/OSGeo/gdal/pull/7389).
  #         could cause issues with macports installation of other things for QGIS but we'll see
  
  # #must run this after updating spatial packages in homebrew, for spatial packages to link properly
  # install.packages("terra", type = "source") 
  # install.packages("sf", type = "source") 
  
  # # this was temporarily needed after updating M1 Macbook to Sequoia OS. not currently needed
  # Sys.setenv(PROJ_LIB = "/opt/homebrew/share/proj")
  # Sys.setenv(GDAL_DATA = "/opt/homebrew/share/gdal")
  
  source(here("src/functions.R"))
  
  # NOTE - consider using glue package to insert dates for output files
  #   https://stackoverflow.com/questions/73584505/how-to-write-a-file-into-a-specific-folder-using-here-package
  
  # potentially useful:
  #https://www.benjaminbell.co.uk/2019/08/bathymetric-maps-in-r-colour-palettes.html
  #https://stackoverflow.com/questions/20581746/increase-resolution-of-color-scale-for-values-close-to-zero
  
  ################################## Set-up ##################################
  
  # Determine a common CRS for the entire dataset
  # common_crs <- st_crs(4326)  # WGS84. suitable for global use. geographic CRS
  common_crs <- "EPSG:4269" # NAD83. like WGS84 but best for contiguous United States. geographic CRS
  # common_crs <- st_crs(32620)  # WGS 84 / UTM zone 20N (suitable for Puerto Rico and the Virgin Islands). projected CRS
  
  # projected_crs <- st_crs(32620)  # WGS 84 / UTM zone 20N (suitable for Puerto Rico and the Virgin Islands). projected CRS
  # projected_crs <- st_crs(26920)  # NAD83 / UTM Zone 20N (suitable for Puerto Rico and the Virgin Islands). projected CRS
  projected_crs <- "EPSG:26920"  # NAD83 / UTM Zone 20N
  
  # Specify the path to ESRI geodatabase provided by Jeremiah Blondeau
  gdb_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Bathymetry/NOAA_LIDAR_Blondeau/US_Caribbean_Bathy_Mocaics.gdb"
  
  ################################## Load geodatabase rasters ##################################
  
  # List all raster layers within the geodatabase
  # Check sub-datasets
  describe(gdb_path, sds = TRUE)
  # raster_layers <- rast(gdb_path)
  # names(raster_layers)  # Print the names of the raster layers
  
  # Load specific sub-datasets from the geodatabase
  # bathy_STTSTJ = rast(gdb_path, subds = "STTSTJ_2m")
  # bathy_STX <- rast(gdb_path, subds = "STX_2m")
  bathy_PR_East <- rast(gdb_path, subds = "PuertoRico_East_2m")
  # bathy_PR_South <- rast(gdb_path, subds = "PuertoRico_South_2m")
  # bathy_PR_West <- rast(gdb_path, subds = "PuertoRico_West_2m")
  # bathy_PR_North <- rast(gdb_path, subds = "PuertoRico_North_2m")
  
  bathy_PR_East <- project(bathy_PR_East, projected_crs)
  
  ################################## Bathymetry resampling  ##################################
  
  # # REDACTED
  # 
  # #downscale the resolution of the bathymetry grid
  # # NOTE - may eventually simply snap the bathymetry to the NCRMP grid itself
  # #
  # # Define the aggregation factor
  # agg_factor <- 25
  # 
  # # Aggregate mean of the bathymetry to 50 m resolution (this may change; was doing this for ease of use early on)
  # bathy_STTSTJ_agg <- aggregate(bathy_STTSTJ, fact = agg_factor, fun = mean, na.rm = TRUE)
  # bathy_STX_agg <- aggregate(bathy_STX, fact = agg_factor, fun = mean, na.rm = TRUE)
  # bathy_PR_East_agg <- aggregate(bathy_PR_East, fact = agg_factor, fun = mean, na.rm = TRUE)
  # bathy_PR_South_agg <- aggregate(bathy_PR_South, fact = agg_factor, fun = mean, na.rm = TRUE)
  # bathy_PR_West_agg <- aggregate(bathy_PR_West, fact = agg_factor, fun = mean, na.rm = TRUE)
  # bathy_PR_North_agg <- aggregate(bathy_PR_North, fact = agg_factor, fun = mean, na.rm = TRUE)
  # # 
  # # # Merge the aggregated rasters. merge is quicker but maybe mosaic is better?
  # # agg_rasters <- c(bathy_STTSTJ_agg, bathy_STX_agg, bathy_PR_East_agg, bathy_PR_South_agg, bathy_PR_West_agg, bathy_PR_North_agg)
  # # bathy_merged <- mosaic(agg_rasters, fun = mean, na.rm = TRUE)
  # bathy_merged_50m = merge(bathy_STTSTJ_agg, bathy_STX_agg, bathy_PR_East_agg, bathy_PR_South_agg, bathy_PR_West_agg, bathy_PR_North_agg)
  # # bathy_merged_2m = merge(bathy_STTSTJ, bathy_STX, bathy_PR_East, bathy_PR_South, bathy_PR_West, bathy_PR_North)
  
  #pull NOAA crm, produced in ~2019 (maybe ? from Dan Holstein) from presumably NOAA bathymetry database (https://www.ncei.noaa.gov/maps/grid-extract/)
  # - CRS: 26920 (NAD83 / UTM zone 20N)
  # - max 30 m resolution
  # This will be used for USVI, since it has nice resolution generally, less artifacts than the crm below, and full extent across Anegada
  crm_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Holstein_VI_Shapes/VI_Shapes/Bathy/crm_usvi.tif"
  describe(crm_path)
  bathy_crm_2019 = rast(crm_path)
  
  bathy_crm_2019 <- project(bathy_crm_2019, projected_crs)
  
  res(bathy_crm_2019)
  crs(bathy_crm_2019)
  crs(bathy_PR_East)
  # crs(bathy_STTSTJ)
  
  #pull NOAA 1 arc-second crm, downloaded in 2024 from NOAA bathymetry database (https://www.ncei.noaa.gov/maps/grid-extract/)
  # - CRS: 4326 (WGS 84; requires reprojection)
  # - max 90 m resolution
  # This will be used for PR, because it has seamless projection with the above when crossing from STT to Culebra/Vieques.
  #   - Would simply use this for entire PR/USVI domain extent, but above crm happens to be a bit better for our use case in the USVI
  #       side. Not sure how to reproduce the crm above unfortunately, but data is available
  # - Extents:
  #   -	North: 19.032
  #   -	South: 16.988
  #   -	East: -64.246
  #   -	West: -68.022
  crm_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Bathymetry/NOAA_other_bathy/CRM_export/exportImage.tiff"
  describe(crm_path)
  bathy_crm_2024 = rast(crm_path)
  res(bathy_crm_2024)
  
  bathy_crm_2024 <- project(bathy_crm_2024, projected_crs)
  
  crs(bathy_crm_2024)
  crs(bathy_crm_2019)
  res(bathy_crm_2019)
  res(bathy_crm_2024)
  
  #read the polygon release grid (substrate for downstream hydrodynamic connectivity simulations) shapefile. CRS: 26920 (NAD83 / UTM zone 20N)
  # reef_polys <- st_read("/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/output/Habitat_Grid/Reef_Polygons/polys_apr2025_operational.shp")
  reef_polys = vect("/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/output/Habitat_Grid/Reef_Polygons/polys_apr2025_operational.shp")
  
  #read the CMS-formatted hydrodynamic domain extent
  hydro_extent = vect(here("output/", "hydro_domain_extent.shp"))
  # hydro_extent_proj <- project(hydro_extent, crs(reef_polys))
  hydro_extent_proj <- project(hydro_extent, projected_crs)
  
  # #clip crm bathymetry to the extent of release grid; reduces processing time of next step
  # # NOTE - will also need to clip out overlap with STTSTJ raster - otherwise, weird artifacts over MCD are brought over in the merge
  # bathy_PR_East_clipped <- crop(bathy_PR_East, reef_polys)
  # # bathy_PR_east_agg_clipped <- crop(bathy_PR_East_agg, reef_polys)
  # bathy_crm_2024_clipped <- crop(bathy_crm_2024, reef_polys)
  # bathy_crm_2019_clipped = crop(bathy_crm_2019, reef_polys)
  
  #clip crm bathymetry to the extent of release grid; reduces processing time of next step
  # NOTE - will also need to clip out overlap with STTSTJ raster - otherwise, weird artifacts over MCD are brought over in the merge
  
  # First crop to reduce processing time (gets the right dimensions)
  bathy_PR_East_clipped <- crop(bathy_PR_East, hydro_extent_proj)
  bathy_crm_2024_clipped <- crop(bathy_crm_2024, hydro_extent_proj)
  bathy_crm_2019_clipped = crop(bathy_crm_2019, hydro_extent_proj)
  
  # Then mask to remove values outside the exact boundary
  bathy_PR_East_clipped <- mask(bathy_PR_East_clipped, hydro_extent_proj)
  bathy_crm_2024_clipped <- mask(bathy_crm_2024_clipped, hydro_extent_proj)
  bathy_crm_2019_clipped <- mask(bathy_crm_2019_clipped, hydro_extent_proj)
  
  # # plot to verify proper crop
  # plot(hydro_extent_proj, main = "Sea Mask with Hydrological Extent Overlay",
  #      col=c("white", "lightblue"), legend = FALSE)
  # plot(bathy_PR_East_clipped, add = TRUE, border = 'darkblue', lwd = 2)
  # plot(bathy_crm_2024_clipped, add = TRUE, border = 'darkblue', lwd = 2)
  # plot(bathy_crm_2019_clipped, add = TRUE, border = 'darkblue', lwd = 2)
  
  #retrieve and apply crm raster extents to fresh raster template, then resample to 50 m resolution using template  # NOTE -
  # NOTE - should resampling be done before or after merging with PR East?
  #      - should I actually be aggregating to a masked / cookie-cutter'd grid which is flush with the 50-m NCRMP sampling grid?
  #         - answer: I think no, since the Puerto Rico grid is a different projection and the different grids would never all align anyways
  e_crm <- ext(bathy_crm_2019_clipped)
  e_pr <- ext(bathy_PR_East_clipped)
  xmin_combined <- min(e_crm$xmin, e_pr$xmin)
  xmax_combined <- max(e_crm$xmax, e_pr$xmax)
  ymin_combined <- min(e_crm$ymin, e_pr$ymin)
  ymax_combined <- max(e_crm$ymax, e_pr$ymax)
  xmin <- floor(xmin_combined / 50) * 50
  xmax <- ceiling(xmax_combined / 50) * 50
  ymin <- floor(ymin_combined / 50) * 50
  ymax <- ceiling(ymax_combined / 50) * 50
  # xmin = 233900 #manual edit to greatly cut down the size of the raster over deep ocean
  # xmax = 385000 #380000 #manual edit to greatly cut down the size of the raster over deep ocean
  # ymin = 1945000 #1940000 #manual edit to greatly cut down the size of the raster over deep ocean
  # ymax = 2087500 #2080000 #manual edit to greatly cut down the size of the raster over deep ocean
  new_ext <- ext(xmin, xmax, ymin, ymax)
  template_raster <- rast(new_ext, resolution = 50, crs = crs(projected_crs))
  bathy_crm_2024_agg <- resample(bathy_crm_2024_clipped, template_raster, method = "average") #NOTE - resampling was done because 'aggregate' could not produce exact discrete 50 x 50 m resolution
  bathy_PR_East_agg <- resample(bathy_PR_East_clipped, template_raster, method = "average")
  bathy_crm_2019_agg = resample(bathy_crm_2019_clipped, template_raster, method = 'average')
  rm(template_raster) #drop this raster after it is no longer needed, since it is empty and messed up object saving function
  
  #plot briefly
  bathy_crm_2024_agg_reefdepth <- clamp(bathy_crm_2024_agg, lower=-50, upper=0, values=TRUE) #limit depth to 0 m; eliminate land elevation
  bathy_PR_East_agg_reefdepth = clamp(bathy_PR_East_agg, lower = -50, upper = 0, values = TRUE)
  bathy_crm_2019_agg_reefdepth = clamp(bathy_crm_2019_agg, lower = -50, upper = 0, values = TRUE)
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(xmin, xmax, ymin, ymax)
  plot(bathy_crm_2024_agg_reefdepth, 
       main="CRM 2024",
       # col=hcl.colors(50, "Blues", rev=TRUE),
       # ext = e_pr, #e_crm,
       ext = plot_extents, #e_crm,
       legend=TRUE)
  plot(bathy_PR_East_agg_reefdepth,
       main = 'Blondeau',
       ext = e_pr,
       # ext = plot_extents, #e_crm,
       legend = TRUE)
  plot(bathy_crm_2019_agg_reefdepth, #this shows that the 2019 crm is the clear winner, for MCD
       main = 'CRM 2019',
       # ext = e_pr,
       ext = plot_extents, #e_crm,
       legend = TRUE)
  
  res(bathy_crm_2024_agg)  # Should be exactly [1] 50 50
  res(bathy_PR_East_agg)  # Should be exactly [1] 50 50
  res(bathy_crm_2019_agg)  # Should be exactly [1] 50 50
  
  ################################## Merge 1: PR to MCD ##################################
  
  #eliminate fill values at edge of MCD
  bathy_crm_2019_agg_nofill <- ifel(bathy_crm_2019_agg < -100000, NA, bathy_crm_2019_agg)
  
  bathy_merge1_crm = merge(bathy_crm_2019_agg_nofill, bathy_crm_2024_agg)
  
  #plot briefly
  bathy_merged_crm_reefdepth <- clamp(bathy_merge1_crm, lower = -50, upper = 0, values = TRUE) #limit depth to 0 m; eliminate land elevation
  plot_extents = ext(xmin, xmax, ymin, ymax)
  plot(bathy_merged_crm_reefdepth, 
       main="Merge #1",
       # col=hcl.colors(50, "Blues", rev=TRUE),
       # ext = e_pr, #e_crm,
       ext = plot_extents, #e_crm,
       legend=TRUE)
  
  ################################## Merge 2: MCD to Coral Bay ##################################
  
  #eliminate fill values at edge of MCD
  bathy_crm_2019_agg_nofill <- ifel(bathy_crm_2019_agg < -100000, NA, bathy_crm_2019_agg)
  
  bathy_crm_2024_agg_reefdepth_nofill <- ifel(bathy_crm_2024_agg_reefdepth == 0, NA, bathy_crm_2024_agg_reefdepth)
  
  kiddel_point_lon = -64.71987
  kiddel_point_lat = 18.30761
  point_wgs84 <- st_sfc(st_point(c(kiddel_point_lon, kiddel_point_lat)), crs = 4326)
  point_utm <- st_transform(point_wgs84, 26920)
  kiddel_point_east_x <- st_coordinates(point_utm)[1]
  cat("Converted longitude", kiddel_point_lon, "to UTM x-coordinate:", kiddel_point_east_x, "\n")
  
  lon_rast <- init(bathy_crm_2024_agg_reefdepth_nofill, "x")
  west_mask <- lon_rast >= kiddel_point_east_x
  STJ_clip <- bathy_crm_2024_agg_reefdepth_nofill
  STJ_clip[!west_mask] <- NA
  
  lon_rast <- init(STJ_clip, "x")
  east_non_na <- max(values(lon_rast)[!is.na(values(STJ_clip))], na.rm = TRUE)
  east_trim_limit <- east_non_na - 1000  # trim 1 km west of easternmost valid value
  east_mask <- lon_rast <= east_trim_limit
  STJ_clip[!east_mask] <- NA
  
  bathy_merge2_crm = merge(STJ_clip, bathy_merged_crm_reefdepth)
  
  bathy_merge2_crm_nofill = ifel(bathy_merge2_crm == 0, NA, bathy_merge2_crm)
  
  #plot briefly
  bathy_merged2_crm_reefdepth <- clamp(bathy_merge2_crm_nofill, lower = -50, upper = 0, values = TRUE) #limit depth to 0 m; eliminate land elevation
  plot_extents = ext(xmin, xmax, ymin, ymax)
  plot(bathy_merged2_crm_reefdepth, 
       main="Merge #2",
       # col=hcl.colors(50, "Blues", rev=TRUE),
       # ext = e_pr, #e_crm,
       ext = plot_extents, #e_crm,
       legend=TRUE)
  
  ################################## Merge 3: STX ##################################
  
  VI_trough_lon = -64.767128
  VI_trough_lat = 18.046230
  point_wgs84 <- st_sfc(st_point(c(VI_trough_lon, VI_trough_lat)), crs = 4326)
  point_utm <- st_transform(point_wgs84, 26920)
  VI_trough_north_y <- st_coordinates(point_utm)[2]
  cat("Converted latitude", VI_trough_lat, "to UTM y-coordinate:", VI_trough_north_y, "\n")
  
  lat_rast <- init(bathy_crm_2019_agg_reefdepth, "y")
  south_mask <- lat_rast <= VI_trough_north_y
  STX_clip <- bathy_crm_2019_agg_reefdepth
  STX_clip[!south_mask] <- NA
  
  bathy_merge3_crm = merge(STX_clip, bathy_merged2_crm_reefdepth)
  
  bathy_merge2_crm_nofill = ifel(bathy_merge3_crm == 0, NA, bathy_merge3_crm)
  
  #plot briefly
  bathy_merged3_crm_reefdepth <- clamp(bathy_merge2_crm_nofill, lower = -50, upper = 0, values = TRUE) #limit depth to 0 m; eliminate land elevation
  plot_extents = ext(xmin, xmax, ymin, ymax)
  plot(bathy_merged3_crm_reefdepth, 
       main="Merge #3",
       # col=hcl.colors(50, "Blues", rev=TRUE),
       # ext = e_pr, #e_crm,
       ext = plot_extents, #e_crm,
       legend=TRUE)
  
  # #compare with 2019 bathy
  # bathy_crm_2019_agg_reefdepth_nofill = ifel(bathy_crm_2019_agg_reefdepth == 0, NA, bathy_crm_2019_agg_reefdepth)
  # plot(bathy_crm_2019_agg_reefdepth_nofill, #this shows that the 2019 crm is the clear winner, for MCD
  #      main = 'CRM 2019',
  #      # ext = e_pr,
  #      ext = plot_extents, #e_crm,
  #      legend = TRUE)
  
  ################################## Save objects/workspace  ##################################
  
  #remove raster files with very large memory which don't work well with saving and re-loading downstream
  # NOTE - can return to this if direct access to PR East is required!
  rm(bathy_PR_East_clipped)
  rm(bathy_PR_East)
  
  #save terra objects #and then workspace for use in downstream scripts
  save_spat_objects(output_dir = 'output/output_import_merge_rasters/') #call from functions.R
  # save.image(file = here("output", 'output_import_merge_rasters/import_merge_rasters_workspace.RData'))
  
  # Get all non-spatial objects
  non_spatial <- ls()[!sapply(ls(), function(x) inherits(get(x), c("SpatRaster", "SpatVector", "SpatExtent")))]
  
  # Save only non-spatial objects
  # NOTE - this helps with avoiding 'pointer' warnings/errors when loading everything again downstream
  save(list = non_spatial, file = here('output', 'output_import_merge_rasters/import_merge_rasters_workspace.RData'))
  
  