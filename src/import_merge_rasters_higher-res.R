  
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
  
  # potentially useful:
  #https://www.benjaminbell.co.uk/2019/08/bathymetric-maps-in-r-colour-palettes.html
  #https://stackoverflow.com/questions/20581746/increase-resolution-of-color-scale-for-values-close-to-zero
  
  ################################## Load geodatabase rasters ##################################
  
  # Specify the path to ESRI geodatabase provided by Jeremiah Blondeau
  gdb_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Bathymetry/NOAA_LIDAR_Blondeau/US_Caribbean_Bathy_Mocaics.gdb"
  
  # List all raster layers within the geodatabase
  describe(gdb_path, sds = TRUE)
  
  # Load specific sub-datasets from the geodatabase
  bathy_STTSTJ = rast(gdb_path, subds = "STTSTJ_2m")
  bathy_STX <- rast(gdb_path, subds = "STX_2m")
  bathy_PR_East <- rast(gdb_path, subds = "PuertoRico_East_2m")
  bathy_PR_South = rast(gdb_path, subds = "PuertoRico_South_2m")
  bathy_PR_West = rast(gdb_path, subds = "PuertoRico_West_2m")
  bathy_PR_North = rast(gdb_path, subds = "PuertoRico_North_2m")
  
  projected_crs = crs(bathy_STTSTJ)
  
  ################################## Bathymetry resampling  ##################################
  
  #pull NOAA crm, produced in ~2019 (maybe ? from Dan Holstein) from presumably NOAA bathymetry database (https://www.ncei.noaa.gov/maps/grid-extract/)
  # CRS: 26920 (NAD83 / UTM zone 20N)
  # 30 m max resolution
  # using this for its good landmask (Blondeau bathymetry has some issues along coastlines)
  #   NOTE - also splicing in the "good" bathy south of STT from these data, to the Blondeau data. that area was higher
  #           resolution for some reason than the current product
  # Can also be used for USVI in the backup merge, since it has nice resolution generally, less artifacts than the crm below, and full extent across Anegada
  crm_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Holstein_VI_Shapes/VI_Shapes/Bathy/crm_usvi.tif"
  describe(crm_path)
  bathy_crm_2019 = rast(crm_path)
  bathy_crm_2019 <- project(bathy_crm_2019, projected_crs)
  
  #pull NOAA 1 arc-second crm, downloaded in 2024 from NOAA bathymetry database (https://www.ncei.noaa.gov/maps/grid-extract/)
  # - CRS: 4326 (WGS 84; requires reprojection)
  # - max 90 m resolution
  # This can be used for PR in backup merge, because it has seamless projection with the above when crossing from STT to Culebra/Vieques.
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
  
  # # VERSION WHERE MODELING IS LIMITED TO VIEQUES EASTWARD
  # #read the CMS-formatted hydrodynamic domain extent
  # hydro_extent = vect(here("output/", "hydro_domain_extent.shp"))
  # hydro_extent_proj <- project(hydro_extent, projected_crs)
  # 
  # # First crop to reduce processing time (gets the right dimensions)
  # bathy_STTSTJ_clipped <- crop(bathy_STTSTJ, hydro_extent_proj)
  # bathy_STX_clipped <- crop(bathy_STX, hydro_extent_proj)
  # bathy_PR_East_clipped <- crop(bathy_PR_East, hydro_extent_proj)
  # bathy_PR_South_clipped <- crop(bathy_PR_South, hydro_extent_proj)
  # bathy_PR_West_clipped <- crop(bathy_PR_West, hydro_extent_proj)
  # bathy_PR_North_clipped <- crop(bathy_PR_North, hydro_extent_proj)
  # bathy_crm_2019_clipped = crop(bathy_crm_2019, hydro_extent_proj)
  # bathy_crm_2024_clipped <- crop(bathy_crm_2024, hydro_extent_proj)
  # 
  # # Then mask to remove values outside the exact boundary
  # bathy_STTSTJ_clipped <- mask(bathy_STTSTJ_clipped, hydro_extent_proj)
  # bathy_STX_clipped <- mask(bathy_STX_clipped, hydro_extent_proj)
  # bathy_PR_East_clipped <- mask(bathy_PR_East_clipped, hydro_extent_proj)
  # bathy_crm_2019_clipped <- mask(bathy_crm_2019_clipped, hydro_extent_proj)
  # bathy_crm_2024_clipped <- mask(bathy_crm_2024_clipped, hydro_extent_proj)
  
  # VERSION WHERE MODELING COVERS ENTIRE US CARIBBEAN
  #
  # First crop to reduce processing time (gets the right dimensions)
  bathy_STTSTJ_clipped <- bathy_STTSTJ
  bathy_STX_clipped <- bathy_STX
  bathy_PR_East_clipped <- bathy_PR_East
  bathy_PR_South_clipped <- bathy_PR_South
  bathy_PR_West_clipped <- bathy_PR_West
  bathy_PR_North_clipped <- bathy_PR_North
  bathy_crm_2019_clipped = bathy_crm_2019
  bathy_crm_2024_clipped <- bathy_crm_2024
  
  #retrieve and apply crm raster extents to fresh raster template, then resample to 50 m resolution using template
  e_crm <- ext(bathy_crm_2019_clipped)
  e_pr_east <- ext(bathy_PR_East_clipped)
  e_pr_south <- ext(bathy_PR_South_clipped)
  e_pr_west <- ext(bathy_PR_West_clipped)
  e_pr_north <- ext(bathy_PR_North_clipped)

  # Find the maximum combined extent of all five rasters
  xmin_combined <- min(e_crm$xmin, e_pr_east$xmin, e_pr_south$xmin, e_pr_west$xmin, e_pr_north$xmin)
  xmax_combined <- max(e_crm$xmax, e_pr_east$xmax, e_pr_south$xmax, e_pr_west$xmax, e_pr_north$xmax)
  ymin_combined <- min(e_crm$ymin, e_pr_east$ymin, e_pr_south$ymin, e_pr_west$ymin, e_pr_north$ymin)
  ymax_combined <- max(e_crm$ymax, e_pr_east$ymax, e_pr_south$ymax, e_pr_west$ymax, e_pr_north$ymax)
  
  # Round to nearest 50m grid
  xmin <- floor(xmin_combined / 50) * 50
  xmax <- ceiling(xmax_combined / 50) * 50
  ymin <- floor(ymin_combined / 50) * 50
  ymax <- ceiling(ymax_combined / 50) * 50
  
  # xmin_combined <- min(e_crm$xmin, e_pr$xmin)
  # xmax_combined <- max(e_crm$xmax, e_pr$xmax)
  # ymin_combined <- min(e_crm$ymin, e_pr$ymin)
  # ymax_combined <- max(e_crm$ymax, e_pr$ymax)
  # xmin <- floor(xmin_combined / 50) * 50
  # xmax <- ceiling(xmax_combined / 50) * 50
  # ymin <- floor(ymin_combined / 50) * 50
  # ymax <- ceiling(ymax_combined / 50) * 50
  # xmin = 233900 #manual edit to greatly cut down the size of the raster over deep ocean
  # xmax = 385000 #380000 #manual edit to greatly cut down the size of the raster over deep ocean
  # ymin = 1945000 #1940000 #manual edit to greatly cut down the size of the raster over deep ocean
  # ymax = 2087500 #2080000 #manual edit to greatly cut down the size of the raster over deep ocean
  
  new_ext <- ext(xmin, xmax, ymin, ymax)
  template_raster <- rast(new_ext, resolution = 50, crs = crs(projected_crs))
  
  bathy_STTSTJ_agg <- resample(bathy_STTSTJ_clipped, template_raster, method = "average")
  bathy_STX_agg <- resample(bathy_STX_clipped, template_raster, method = "average")
  bathy_PR_East_agg <- resample(bathy_PR_East_clipped, template_raster, method = "average")
  bathy_PR_South_agg <- resample(bathy_PR_South_clipped, template_raster, method = "average")
  bathy_PR_West_agg <- resample(bathy_PR_West_clipped, template_raster, method = "average")
  bathy_PR_North_agg <- resample(bathy_PR_North_clipped, template_raster, method = "average")
  bathy_crm_2019_agg = resample(bathy_crm_2019_clipped, template_raster, method = 'average')
  bathy_crm_2024_agg <- resample(bathy_crm_2024_clipped, template_raster, method = "average") #NOTE - resampling was done because 'aggregate' could not produce exact discrete 50 x 50 m resolution
  rm(template_raster) #drop this raster after it is no longer needed, since it is empty and messes up object saving function
  
  #combine Blondeau raster chunks
  bathy_Blondeau_agg = merge(bathy_STTSTJ_agg, bathy_STX_agg, bathy_PR_East_agg, bathy_PR_South_agg, bathy_PR_West_agg,
                             bathy_PR_North_agg)
  
  
  # STOPPING POINT - 30 JUNE 2025
  #   1.) splice in Holstein bathy for BVI
  #   2.) consider Edmunds / CSUN coral data
  #   3.) splice in fix for the "pit" in the MCD, and the "tear" between PR & STT (if needed)
  
  #plot briefly
  bathy_Blondeau_agg_plot = clamp(bathy_Blondeau_agg, lower = -50, upper = 0, values = TRUE)
  bathy_crm_2019_agg_plot = clamp(bathy_crm_2019_agg, lower = -50, upper = 0, values = TRUE)
  bathy_crm_2024_agg_plot <- clamp(bathy_crm_2024_agg, lower=-50, upper=0, values=TRUE) #limit depth to 0 m; eliminate land elevation
  plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  plot(bathy_Blondeau_agg_plot,
       main = 'Blondeau',
       # ext = e_pr,
       ext = plot_extents, #e_crm,
       legend = TRUE)
  plot(bathy_crm_2019_agg_plot,
       main = 'CRM 2019',
       # ext = e_pr,
       ext = plot_extents, #e_crm,
       legend = TRUE)
  plot(bathy_crm_2024_agg_plot, 
       main="CRM 2024",
       # col=hcl.colors(50, "Blues", rev=TRUE),
       # ext = e_pr, #e_crm,
       ext = plot_extents, #e_crm,
       legend=TRUE)
  
  res(bathy_PR_East_agg)  # Should be exactly [1] 50 50
  res(bathy_crm_2019_agg)  # Should be exactly [1] 50 50
  res(bathy_crm_2024_agg)  # Should be exactly [1] 50 50
  
  
  
  
  ################################## Apply "good" landmask ##################################
  
  #read in version of 90m crm that spans entire US Caribbean. has a mostly normal landmask; some inland water needs to be filled
  crm_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Bathymetry/NOAA_other_bathy/CRM_export/PR_wide_90m.tiff"
  describe(crm_path)
  bathy_crm_USCarib = rast(crm_path)
  bathy_crm_USCarib <- project(bathy_crm_USCarib, projected_crs)
  
  # Create landmask: values > 0 become NA (land), values <= 0 remain as water
  landmask_USCarib <- bathy_crm_USCarib
  landmask_USCarib[landmask_USCarib > 0] <- NA
  

  # # Extend the landmask to match the full extent of bathy_Blondeau_agg
  # # This fills the extended areas with NA, which won't mask the Blondeau data
  # landmask_extended <- extend(landmask_USCarib, bathy_Blondeau_agg)
  # 
  # # Resample to match resolution
  # landmask_resampled <- resample(landmask_extended, bathy_Blondeau_agg, method = "near")
  # 
  # # Apply the mask - areas where landmask is NA (land) will be masked out,
  # # areas where landmask was extended (also NA) will remain unmasked
  # bathy_Blondeau_masked <- mask(bathy_Blondeau_agg, landmask_resampled, maskvalue = NA, updatevalue = NA)
  
  
  # Extend the landmask to match the full extent of bathy_Blondeau_agg
  # BUT fill the extended areas with a water value (e.g., -999) instead of NA
  landmask_extended <- extend(landmask_USCarib, bathy_Blondeau_agg, fill = -999)
  
  # Resample to match resolution
  landmask_resampled <- resample(landmask_extended, bathy_Blondeau_agg, method = "near")
  
  # Now mask only where landmask has actual NA values (land areas)
  # Areas with -999 (extended water areas) won't be masked
  bathy_Blondeau_masked <- mask(bathy_Blondeau_agg, landmask_resampled, maskvalue = NA)
  
  
  
  bathy_Blondeau_masked_plot = clamp(bathy_Blondeau_masked, lower = -50, upper = 0, values = TRUE)
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  plot_extents = ext(300000, 340000, 1940000, 1980000) #for investigating St Croix
  # plot_extents = ext(290000, 330000, 2040000, 2080000) #for investigating St Thomas
  # plot_extents = ext(240000, 280000, 2000000, 2040000) #for investigating Mona Island
  plot(bathy_Blondeau_masked_plot,
       main = 'Blondeau',
       # ext = e_pr,
       ext = plot_extents, #e_crm,
       legend = TRUE)
  
  
  ################################## seal off inland bathy ##################################
  
  # Create a binary land/water raster (1 = water, NA = land)
  water_binary <- bathy_Blondeau_masked
  water_binary[water_binary > 0] <- NA  # land = NA
  water_binary[water_binary <= 0] <- 1  # water = 1
  
  # Find connected water patches
  water_patches <- patches(water_binary, directions = 8, allowGaps = FALSE)
  
  # Get patch frequencies
  patch_freq <- freq(water_patches)
  
  # Keep the largest patch AND any patches above a size threshold
  min_patch_size <- 1000  # adjust this threshold as needed (number of cells)
  large_patches <- patch_freq[patch_freq$count >= min_patch_size, "value"]
  
  # Create mask: keep large water bodies, fill small ones
  ocean_mask <- water_patches
  ocean_mask[!water_patches %in% large_patches] <- NA  # small water becomes NA
  ocean_mask[water_patches %in% large_patches] <- 1    # large water remains
  
  # Apply mask
  bathy_sealed <- mask(bathy_Blondeau_masked, ocean_mask)
  
  bathy_sealed_plot = clamp(bathy_sealed, lower = -50, upper = 0, values = TRUE)
  plot_extents = ext(280000, 310000, 2000000, 2040000) #for investigating south of STT
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # plot_extents = ext(300000, 340000, 1940000, 1980000) #for investigating St Croix
  # plot_extents = ext(290000, 330000, 2040000, 2080000) #for investigating St Thomas
  # plot_extents = ext(240000, 280000, 2000000, 2040000) #for investigating Mona Island
  plot(bathy_sealed_plot,
       main = 'Blondeau',
       # ext = e_pr,
       ext = plot_extents, #e_crm,
       legend = TRUE)
  
  # # Sometimes a simple morphological operation works well
  # # Create land mask (areas > 0)
  # land_mask <- bathy_Blondeau_masked > 0
  # 
  # # Buffer land areas inward to close small water gaps
  # land_buffered <- buffer(land_mask, width = 100)  # adjust width as needed
  # 
  # # Apply buffered land mask
  # bathy_sealed <- mask(bathy_Blondeau_masked, land_buffered, inverse = TRUE)
  
  ################################## Splice "good" bathy to Blondeau ##################################
  
  south_of_STT_clipped = rast(here('data/south_of_STT_clipped.tif'))

  # Make sure it's in the same projection as bathy_sealed
  south_of_STT_clipped <- project(south_of_STT_clipped, crs(projected_crs))
  
  # Recreate the template raster using the same approach as upstream
  # (You can reuse the same extent calculations if they're still in memory)
  new_ext <- ext(xmin, xmax, ymin, ymax)
  template_raster <- rast(new_ext, resolution = 50, crs = crs(projected_crs))
  
  # Resample to 50m grid using the template
  south_of_STT_agg <- resample(south_of_STT_clipped, template_raster, method = "average")
  
  # Clean up template
  rm(template_raster)
  
  # bathy_Blondeau_final <- merge(south_of_STT_clipped, bathy_sealed)
  bathy_Blondeau_final <- mosaic(south_of_STT_agg, bathy_sealed, fun = "first")
  
  
  bathy_Blondeau_final_plot = clamp(bathy_Blondeau_final, lower = -50, upper = 0, values = TRUE)
  plot_extents = ext(280000, 310000, 2000000, 2040000) #for investigating south of STT
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # plot_extents = ext(300000, 340000, 1940000, 1980000) #for investigating St Croix
  # plot_extents = ext(240000, 280000, 2000000, 2040000) #for investigating Mona Island
  plot(bathy_Blondeau_final_plot,
       main = 'Blondeau',
       # ext = e_pr,
       ext = plot_extents, #e_crm,
       legend = TRUE)
  
  south_of_STT_clipped_plot = clamp(south_of_STT_clipped, lower = -50, upper = 0, values = TRUE)
  plot_extents = ext(280000, 310000, 2000000, 2040000) #for investigating south of STT
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # plot_extents = ext(300000, 340000, 1940000, 1980000) #for investigating St Croix
  # plot_extents = ext(240000, 280000, 2000000, 2040000) #for investigating Mona Island
  plot(south_of_STT_clipped_plot,
       main = 'Blondeau',
       # ext = e_pr,
       ext = plot_extents, #e_crm,
       legend = TRUE)
  
  
  ################################## Merge 1: PR to MCD ##################################
  
  # NOTE - the below 3 merges are for a "back-up" version of merged bathymetry that doesn't have the artifact
  #         issues the Blondeau set does
  
  
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
  
  bathy_crm_2024_agg_reefdepth_nofill <- ifel(bathy_crm_2024_agg_plot == 0, NA, bathy_crm_2024_agg_plot)
  
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
  
  lat_rast <- init(bathy_crm_2019_agg_plot, "y")
  south_mask <- lat_rast <= VI_trough_north_y
  STX_clip <- bathy_crm_2019_agg_plot
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
  # bathy_crm_2019_agg_reefdepth_nofill = ifel(bathy_crm_2019_agg_plot == 0, NA, bathy_crm_2019_agg_plot)
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
  save_spat_objects(output_dir = 'output/output_import_merge_rasters_higher-res/') #call from functions.R
  
  # Get all non-spatial objects
  non_spatial <- ls()[!sapply(ls(), function(x) inherits(get(x), c("SpatRaster", "SpatVector", "SpatExtent")))]
  
  # Save only non-spatial objects
  # NOTE - this helps with avoiding 'pointer' warnings/errors when loading everything again downstream
  save(list = non_spatial, file = here('output', 'output_import_merge_rasters_higher-res/import_merge_rasters_workspace.RData'))
  
  