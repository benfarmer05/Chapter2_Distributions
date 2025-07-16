  
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
  
  # # 2024 Blondeau
  # # Specify the path to ESRI geodatabase provided by Jeremiah Blondeau
  # gdb_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Bathymetry/NOAA_LIDAR_Blondeau/US_Caribbean_Bathy_Mocaics.gdb"
  # 
  # # List all raster layers within the geodatabase
  # describe(gdb_path, sds = TRUE)
  # 
  # # Load specific sub-datasets from the geodatabase
  # bathy_STTSTJ = rast(gdb_path, subds = "STTSTJ_2m")
  # bathy_STX <- rast(gdb_path, subds = "STX_2m")
  # bathy_PR_East <- rast(gdb_path, subds = "PuertoRico_East_2m")
  # bathy_PR_South = rast(gdb_path, subds = "PuertoRico_South_2m")
  # bathy_PR_West = rast(gdb_path, subds = "PuertoRico_West_2m")
  # bathy_PR_North = rast(gdb_path, subds = "PuertoRico_North_2m")
  # 
  # projected_crs = crs(bathy_STTSTJ)
  # geographic_crs = crs("EPSG:4269")
  
  # 2025 Blondeau
  bathy_STTSTJ <- rast(here("data", "NOAA_LIDAR_Blondeau_2025", "PR_bath2m_StTJ2.tif"))
  bathy_STX <- rast(here("data", "NOAA_LIDAR_Blondeau_2025", "PR_bath2m_stx2.tif"))
  bathy_PR_East <- rast(here("data", "NOAA_LIDAR_Blondeau_2025", "PR_bath2m_east2.tif"))
  bathy_PR_North <- rast(here("data", "NOAA_LIDAR_Blondeau_2025", "PR_bath2m_north.tif"))
  bathy_PR_South <- rast(here("data", "NOAA_LIDAR_Blondeau_2025", "PR_bath2m_south.tif"))
  bathy_PR_West <- rast(here("data", "NOAA_LIDAR_Blondeau_2025", "PR_bath2m_west.tif"))
  
  projected_crs <- crs(bathy_PR_East)
  geographic_crs <- crs("EPSG:4269")
  
  ################################## Bathymetry resampling  ##################################
  
  #pull NOAA crm, produced in ~2019 (maybe ? from Dan Holstein) from presumably NOAA bathymetry database (https://www.ncei.noaa.gov/maps/grid-extract/)
  # CRS: 26920 (NAD83 / UTM zone 20N)
  # 30 m max resolution
  # using this for its good landmask (Blondeau bathymetry has some issues along coastlines)
  #   NOTE - also splicing in the "good" bathy south of STT from these data, to the Blondeau data. that area was higher
  #           resolution for some reason than the current product
  # Can also be used for USVI in the backup merge, since it has nice resolution generally, less artifacts than the crm below, and full extent across Anegada
  bathy_crm_2019 = rast(here('data/crm_usvi.tif'))
  # crm_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Holstein_VI_Shapes/VI_Shapes/Bathy/crm_usvi.tif"
  # describe(crm_path)
  # bathy_crm_2019 = rast(crm_path)
  bathy_crm_2019 <- project(bathy_crm_2019, projected_crs)
  # bathy_crm_2019 <- project(bathy_crm_2019, geographic_crs)
  bathy_crm_2019 = clamp(bathy_crm_2019, lower = -300, upper = 0, values = TRUE) #make the bathy sensible - eliminates fill values
  
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
  bathy_crm_2024 = rast(here('data/exportImage.tiff'))
  # crm_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Bathymetry/NOAA_other_bathy/CRM_export/exportImage.tiff"
  # describe(crm_path)
  # bathy_crm_2024 = rast(crm_path)
  res(bathy_crm_2024)
  bathy_crm_2024 <- project(bathy_crm_2024, projected_crs)
  # bathy_crm_2024 <- project(bathy_crm_2024, geographic_crs)
  
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
  
  # Filter out bathymetry values deeper than 300m for 2025 Blondeau data only
  bathy_STTSTJ_clipped <- clamp(bathy_STTSTJ_clipped, lower = -300, upper = 0, values = FALSE)
  bathy_STX_clipped <- clamp(bathy_STX_clipped, lower = -300, upper = 0, values = FALSE)
  bathy_PR_East_clipped <- clamp(bathy_PR_East_clipped, lower = -300, upper = 0, values = FALSE)
  bathy_PR_South_clipped <- clamp(bathy_PR_South_clipped, lower = -300, upper = 0, values = FALSE)
  bathy_PR_West_clipped <- clamp(bathy_PR_West_clipped, lower = -300, upper = 0, values = FALSE)
  bathy_PR_North_clipped <- clamp(bathy_PR_North_clipped, lower = -300, upper = 0, values = FALSE)
  
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
  
  #combine Blondeau raster chunks
  bathy_Blondeau_agg = merge(bathy_STTSTJ_agg, bathy_STX_agg, bathy_PR_East_agg, bathy_PR_South_agg, bathy_PR_West_agg,
                             bathy_PR_North_agg)
  
  bathy_Blondeau_native_agg = merge(bathy_STTSTJ_clipped, bathy_STX_clipped, bathy_PR_East_clipped, bathy_PR_South_clipped,
                                    bathy_PR_West_clipped, bathy_PR_North_clipped)
  
  
  # NOTE - 3 JULY 2025
  #   2.) consider Edmunds / CSUN coral data [probably not necessary honestly]
  #   3.) splice in fix for the "pit" in the MCD, and the "tear" between PR & STT (if needed)
  
  #plot briefly
   
  
  
  ################################## Apply "good" landmask ##################################
  
  # NOTE - consider if the methods below needlessly introduce more NAs around land than are necessary
  
  #read in version of 90m crm that spans entire US Caribbean. has a mostly normal landmask; some inland water needs to be filled
  bathy_crm_USCarib = rast(here('data/PR_wide_90m.tiff'))
  # crm_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Bathymetry/NOAA_other_bathy/CRM_export/PR_wide_90m.tiff"
  # describe(crm_path)
  # bathy_crm_USCarib = rast(crm_path)
  bathy_crm_USCarib <- project(bathy_crm_USCarib, projected_crs)
  
  # Create landmask: values > 0 become NA (land), values <= 0 remain as water
  landmask_USCarib <- bathy_crm_USCarib
  landmask_USCarib[landmask_USCarib > 0] <- NA
  
  # # Extend the landmask to match the full extent of bathy_Blondeau_agg
  # # BUT fill the extended areas with a water value (e.g., -999) instead of NA
  # landmask_extended <- extend(landmask_USCarib, bathy_Blondeau_agg, fill = -999)
  # landmask_native_extended <- extend(landmask_USCarib, bathy_Blondeau_native_agg, fill = -999)
  
  # # Resample to match resolution
  # landmask_resampled <- resample(landmask_extended, bathy_Blondeau_agg, method = "near")
  # landmask_native_resampled <- resample(landmask_native_extended, bathy_Blondeau_native_agg, method = "near")
  # 
  # # apply the mask
  # real_landmask_area <- !is.na(resample(landmask_USCarib, bathy_Blondeau_agg, method = "near"))
  # bathy_Blondeau_masked <- bathy_Blondeau_agg
  # # bathy_Blondeau_masked[real_landmask_area & is.na(landmask_resampled)] <- NA
  # bathy_Blondeau_masked[landmask_resampled != -999 & is.na(landmask_resampled)] <- NA
  # 
  # # bathy_Blondeau_masked <- mask(bathy_Blondeau_agg, landmask_resampled, maskvalue = NA)
  # bathy_Blondeau_native_masked <- bathy_Blondeau_native_agg
  # bathy_Blondeau_native_masked[!is.na(landmask_native_resampled) & landmask_native_resampled > 0] <- NA
  # # bathy_Blondeau_native_masked <- mask(bathy_Blondeau_native_agg, landmask_native_resampled, maskvalue = NA)
  
  # Simple approach: just resample the landmask to match Blondeau resolution
  landmask_resampled <- resample(landmask_USCarib, bathy_Blondeau_agg, method = "near")
  # landmask_native_resampled <- resample(landmask_USCarib, bathy_Blondeau_native_agg, method = "near")
  
  # Apply the mask - only where landmask has data will masking occur
  bathy_Blondeau_masked <- bathy_Blondeau_agg
  # bathy_Blondeau_native_masked <- bathy_Blondeau_native_agg
  
  # Create a mask for where we actually have landmask coverage
  has_landmask_coverage <- !is.na(resample(bathy_crm_USCarib, bathy_Blondeau_agg, method = "near"))
  # has_landmask_coverage_native <- !is.na(resample(bathy_crm_USCarib, bathy_Blondeau_native_agg, method = "near"))
  
  # Apply landmask only where we have coverage AND the landmask indicates land
  bathy_Blondeau_masked[has_landmask_coverage & is.na(landmask_resampled)] <- NA
  # bathy_Blondeau_native_masked[has_landmask_coverage_native & is.na(landmask_native_resampled)] <- NA
  
  # Add Mona Island landmask
  # Create landmask for Mona Island from bathy_crm_2024
  landmask_mona <- bathy_crm_2024
  landmask_mona[landmask_mona > 0] <- NA  # Set land (>0) to NA
  
  # Clip to Mona Island area using UTM coordinates
  mona_extent <- ext(-25000, 0, 2000000, 2020000)  # Adjust coordinates as needed
  landmask_mona_clipped <- crop(landmask_mona, mona_extent)
  
  # Resample Mona landmask
  landmask_mona_resampled <- resample(landmask_mona_clipped, bathy_Blondeau_agg, method = "near")
  # landmask_mona_native_resampled <- resample(landmask_mona_clipped, bathy_Blondeau_native_agg, method = "near")
  
  # Create a mask for where the clipped Mona landmask has any data at all
  has_mona_data <- !is.na(crop(bathy_crm_2024, mona_extent))
  has_mona_data_resampled <- resample(has_mona_data, bathy_Blondeau_agg, method = "near")
  # has_mona_data_native_resampled <- resample(has_mona_data, bathy_Blondeau_native_agg, method = "near")
  
  # Apply Mona landmask only in areas where we have Mona data AND it's land
  bathy_Blondeau_masked[has_mona_data_resampled & is.na(landmask_mona_resampled)] <- NA
  # bathy_Blondeau_native_masked[has_mona_data_native_resampled & is.na(landmask_mona_native_resampled)] <- NA
  
  
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
       # ext = plot_extents, #e_crm,
       legend = TRUE)
  
  ################################## Splice "good" bathy to Blondeau ##################################
  
  south_of_STT_clipped = rast(here('data/south_of_STT_clipped.tif'))

  # Make sure it's in the right projection
  south_of_STT_clipped <- project(south_of_STT_clipped, crs(projected_crs))
  
  # Recreate the template raster using the same approach as upstream
  # (You can reuse the same extent calculations if they're still in memory)
  new_ext <- ext(xmin, xmax, ymin, ymax)
  template_raster <- rast(new_ext, resolution = 50, crs = crs(projected_crs))
  # template_native_raster <- rast(new_ext, resolution = 2, crs = crs(projected_crs))
  
  # Resample to 50m grid using the template
  south_of_STT_agg <- resample(south_of_STT_clipped, template_raster, method = "average")
  # south_of_STT_native_agg <- resample(south_of_STT_clipped, template_native_raster, method = "average")
  # south_of_STT_extended <- extend(south_of_STT_clipped, bathy_Blondeau_native_agg)
  south_of_STT_resampled <- resample(south_of_STT_clipped, bathy_Blondeau_native_agg, method = "bilinear")
  
  # bathy_Blondeau_spliced <- merge(south_of_STT_clipped, bathy_Blondeau_masked)
  bathy_Blondeau_spliced <- mosaic(south_of_STT_agg, bathy_Blondeau_masked, fun = "first")
  # bathy_Blondeau_native_spliced <- mosaic(south_of_STT_clipped, bathy_Blondeau_native_agg, fun = "first") # THIS TAKES FOREVER!!
  bathy_Blondeau_native_spliced <- cover(south_of_STT_resampled, bathy_Blondeau_native_agg)
  
  bathy_Blondeau_spliced_plot = clamp(bathy_Blondeau_spliced, lower = -50, upper = 0, values = TRUE)
  bathy_Blondeau_native_spliced_plot = clamp(bathy_Blondeau_native_spliced, lower = -50, upper = 0, values = TRUE)
  plot_extents = ext(280000, 310000, 2000000, 2040000) #for investigating south of STT
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # plot_extents = ext(300000, 340000, 1940000, 1980000) #for investigating St Croix
  # plot_extents = ext(240000, 280000, 2000000, 2040000) #for investigating Mona Island
  plot(bathy_Blondeau_native_spliced_plot,
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
  
  ################################## Tack BVI to Blondeau ##################################
  
  
  # STOPPING POINT - 15 JULY 2025
  #   - tedious work to get very very large rasters to chug and save/load progress...
  #   - currently figuring out this BVI step. I think there is an issue with extent mismatch for native
  
  
  # # METHOD WHICH FILLS IN - EVEN AT MESOPHOTIC RIDGE SOUTH OF STT
  # #
  # plot(bathy_Blondeau_final)
  # plot(bathy_crm_2019_agg)
  # blondeau_depth_range <- range(values(bathy_Blondeau_final), na.rm = TRUE)
  # bathy_crm_filtered <- bathy_crm_2019_agg
  # bathy_crm_filtered[bathy_crm_filtered < blondeau_depth_range[1] |
  #                      bathy_crm_filtered > blondeau_depth_range[2]] <- NA
  # bathy_Blondeau_final = mosaic(bathy_Blondeau_final, bathy_crm_filtered, fun = 'first')
  # plot(bathy_Blondeau_final)
  
  # METHOD WHICH ONLY ADDS NORTH/EAST (BVI) - NO FILLING OF EXISTING DATA
  #
  # Clean up the CRM data before clipping
  # Find the deepest value in Bl# METHOD WHICH ONLY ADDS NORTH/EAST (BVI) - NO FILLING OF EXISTING DATA
  #
  # Clean up the CRM data before clipping
  # Find the deepest value in Blondeau (most negative for bathymetry)
  deepest_blondeau <- min(values(bathy_Blondeau_spliced), na.rm = TRUE)
  print(paste("Deepest value in Blondeau:", deepest_blondeau))
  
  # Set CRM values deeper than Blondeau's deepest to NA
  # Also set values of 0 (land) to NA
  # ADDITIONAL: Set shallow values (> -2m) to NA as these are likely land/very shallow areas
  deepest_crm <- min(values(bathy_crm_2019_agg), na.rm = TRUE)
  print(paste("Deepest value in CRM:", deepest_crm))
  bathy_crm_cleaned <- bathy_crm_2019_agg
  bathy_crm_native_cleaned <- bathy_crm_2019
  # values(bathy_crm_cleaned)[values(bathy_crm_cleaned) < deepest_blondeau] <- NA
  # values(bathy_crm_cleaned)[values(bathy_crm_cleaned) >= -2] <- NA  # Remove land and very shallow areas
  values(bathy_crm_cleaned)[values(bathy_crm_cleaned) >= 0] <- NA  # Remove land and very shallow areas
  values(bathy_crm_cleaned)[values(bathy_crm_cleaned) <= deepest_crm] <- NA  # Remove land and very shallow areas
  
  # Create a binary mask where Blondeau has non-NA data
  blondeau_has_data <- !is.na(bathy_Blondeau_spliced)
  blondeau_native__has_data <- !is.na(bathy_Blondeau_native_spliced)
  
  # Apply this mask to remove CRM data where Blondeau has data
  # We want to KEEP CRM data where blondeau_has_data is FALSE (purple areas)
  # So we mask where blondeau_has_data is TRUE (yellow areas) - NO inverse needed
  bathy_crm_clipped <- mask(bathy_crm_cleaned, blondeau_has_data, inverse = FALSE, maskvalue = TRUE)
  bathy_crm_native_clipped <- mask(bathy_crm_cleaned, blondeau_native__has_data, inverse = FALSE, maskvalue = TRUE)
  plot(bathy_crm_clipped)
  
  # # SPATIAL CLEANUP: Remove CRM data south and west of intended extension area
  # # Based on your image, keep only CRM data that is:
  # # - East of approximately 330000 UTM AND/OR
  # # - North of approximately 2070000 UTM
  # # Get coordinates for each cell
  # coords <- xyFromCell(bathy_crm_clipped, 1:ncell(bathy_crm_clipped))
  # x_coords <- coords[, 1]  # Easting
  # y_coords <- coords[, 2]  # Northing
  # 
  # # Create a logical mask: keep only areas that are sufficiently east OR north
  # # Adjust these coordinates based on your specific area
  # # keep_mask <- (x_coords > 335500) | (y_coords > 2070000)
  # # keep_mask <- (x_coords > 330000) | (y_coords > 2070000)
  # keep_mask <- (x_coords > 330000) | (y_coords > 2065000)
  # 
  # # Apply the spatial mask
  # values(bathy_crm_clipped)[!keep_mask] <- NA
  # plot(bathy_crm_clipped)
  
  # Merge the datasets
  # The merge function will use Blondeau data where available, cleaned CRM data elsewhere
  bathy_Blondeau_BVI <- merge(bathy_Blondeau_spliced, bathy_crm_clipped)
  bathy_Blondeau_native_BVI <- merge(bathy_Blondeau_native_spliced, bathy_crm_native_clipped)
  plot(bathy_Blondeau_BVI)
  
  ################################## patch in USGS Anegada LIDAR ##################################
  
  #read in version of 90m crm that spans entire US Caribbean. has a mostly normal landmask; some inland water needs to be filled
  Anegada_LIDAR = rast(here('data/ANGD2014_EAARLB_z20_v09g12A_mosaic.tif'))
  # crm_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Bathymetry/USGS_Anegada/ANGD2014_EAARLB_z20_v09g12A_mosaic/ANGD2014_EAARLB_z20_v09g12A_mosaic.tif"
  # describe(crm_path)
  # Anegada_LIDAR = rast(crm_path)
  Anegada_LIDAR <- project(Anegada_LIDAR, projected_crs)
  
  #drop fill values and land
  values(Anegada_LIDAR)[values(Anegada_LIDAR) < deepest_blondeau] <- NA
  values(Anegada_LIDAR)[values(Anegada_LIDAR) > 0] <- NA
  
  #resample to 50 m res
  Anegada_LIDAR_agg <- resample(Anegada_LIDAR, template_raster, method = "average")
  
  # Clean up template
  rm(template_raster)
  
  plot_extents = ext(341000, 379000, 2057000, 2078000) # for investigating Anegada
  bathy_colors <- colorRampPalette(c("lightcyan", "cyan", "deepskyblue", "royalblue", "navy"))(100)
  
  # Plot the bathymetry first
  plot(Anegada_LIDAR_agg,  #bathy_merged3_crm_reefdepth
       main="Bathymetry with Sea Mask and Hydrological Extent", 
       col=bathy_colors,
       ext = plot_extents, #e_crm,
       legend=TRUE)
  
  plot(Anegada_LIDAR,  #bathy_merged3_crm_reefdepth
       main="Bathymetry with Sea Mask and Hydrological Extent", 
       col=bathy_colors,
       ext = plot_extents, #e_crm,
       legend=TRUE)
  
  #merge with main
  bathy_Blondeau_BVI_Anegada <- merge(Anegada_LIDAR_agg, bathy_Blondeau_BVI)
  
  plot(bathy_Blondeau_BVI,  #bathy_merged3_crm_reefdepth
       main="Bathymetry with Sea Mask and Hydrological Extent", 
       col=bathy_colors,
       ext = plot_extents, #e_crm,
       legend=TRUE)
  
  plot(bathy_Blondeau_BVI_Anegada,  #bathy_merged3_crm_reefdepth
       main="Bathymetry with Sea Mask and Hydrological Extent", 
       col=bathy_colors,
       ext = plot_extents, #e_crm,
       legend=TRUE)
  
  ################################## seal off inland bathy ##################################
  
  # # Sometimes a simple morphological operation works well
  # # Create land mask (areas > 0)
  # land_mask <- bathy_Blondeau_masked > 0
  #
  # # Buffer land areas inward to close small water gaps
  # land_buffered <- buffer(land_mask, width = 100)  # adjust width as needed
  #
  # # Apply buffered land mask
  # bathy_sealed <- mask(bathy_Blondeau_masked, land_buffered, inverse = TRUE)


  # NOTE - since the below is done POST-resampling to 50 m, some inland water bodies that
  #         truly are connected to the ocean might appear to be disconnected. I think the
  #         below actually doesn't do a good job of dropping those from the seamask anyways
  #         (which is good ??), but worth considering. running this on raw bathy might be very
  #         resource intensive and not worth it. because at the end of the day, resampling
  #         will introduce some artifacts no matter what - really more of a visual issue than
  #         anything. depends how "realistic" we want the final coral cover maps

  solution_1_no_chunking <- function(water_binary) {
    # Set terra options for memory management
    terraOptions(memfrac = 0.8, tempdir = tempdir())
    
    # Try processing without chunking - terra is more memory efficient than raster
    water_patches <- patches(water_binary, directions = 8, allowGaps = FALSE)
    
    return(water_patches)
  }
  
  # Create a binary land/water raster (1 = water, NA = land)
  water_binary <- bathy_Blondeau_BVI_Anegada 
  water_binary[water_binary > 0] <- NA  # land = NA
  water_binary[water_binary <= 0] <- 1  # water = 1
  
  # Find connected water patches
  water_patches <- solution_1_no_chunking(water_binary)
  
  # Get patch frequencies
  patch_freq <- freq(water_patches)
  
  # DEPTH-BASED FUNCTION (add this function definition somewhere above)
  preserve_ocean_patches_depth <- function(water_patches, bathy_Blondeau_BVI_Anegada, 
                                           min_patch_size = 10000, depth_threshold = -10) {
    
    patch_freq <- freq(water_patches)
    
    # Large patches are kept automatically (regardless of depth)
    large_patches <- patch_freq[patch_freq$count >= min_patch_size, "value"]
    
    # Smaller patches: keep only if they have deep water (real ocean)
    smaller_patches <- patch_freq[patch_freq$count < min_patch_size & 
                                    patch_freq$count > 100, "value"]
    
    patches_to_keep <- large_patches
    
    for (patch_id in smaller_patches) {
      patch_mask <- water_patches == patch_id
      
      # Check depth - if it has deep water, it's ocean (not inland lake)
      patch_depths <- mask(bathy_Blondeau_BVI_Anegada, patch_mask)
      min_depth <- global(patch_depths, "min", na.rm = TRUE)$min
      
      if (!is.na(min_depth) && min_depth <= depth_threshold) {
        patches_to_keep <- c(patches_to_keep, patch_id)
        print(paste("Keeping patch", patch_id, "- deep water (min depth:", round(min_depth, 1), "m)"))
      } else {
        print(paste("Removing patch", patch_id, "- shallow/inland (min depth:", round(min_depth, 1), "m)"))
      }
    }
    
    return(patches_to_keep)
  }
  
  # REPLACE THIS LINE:
  # large_patches <- patch_freq[patch_freq$count >= min_patch_size, "value"]
  
  # WITH THIS:
  min_patch_size <- 10000  # or whatever threshold you want
  large_patches <- preserve_ocean_patches_depth(water_patches, 
                                                bathy_Blondeau_BVI_Anegada,
                                                min_patch_size = min_patch_size, 
                                                depth_threshold = -10)
  
  # Rest of your workflow stays the same:
  # Create mask: keep large water bodies, fill small ones
  ocean_mask <- water_patches
  ocean_mask[!water_patches %in% large_patches] <- NA  # small water becomes NA
  ocean_mask[water_patches %in% large_patches] <- 1    # large water remains
  
  # Apply mask
  bathy_sealed <- mask(bathy_Blondeau_BVI_Anegada, ocean_mask)
  
  # NOTE / STOPPING POINT - 7 JULY 2025
  #   - ideally would be able to define where the sealed "coastline" contours are, and just "fill in"
  #       anything within those contours. this may be harder than expected, though. would also be nice
  #       if I end up looking at "distance from coastline", etc. but I guess could just use landmask for
  #       that
  
  bathy_final = bathy_sealed
  
  
  
  
  # # STOPPING POINT - 15 JULY 2025
  # #   - figuring out how to apply the landmask and filled-in land to the native resolution raster efficiently
  # # Apply the same landmask pattern from bathy_Blondeau_masked to the native resolution data
  # bathy_Blondeau_native_masked <- bathy_Blondeau_native_agg
  # 
  # # Create a mask from where bathy_Blondeau_masked is NA
  # landmask_pattern <- is.na(bathy_Blondeau_masked)
  # 
  # # Resample this mask to match the native resolution
  # landmask_pattern_native <- resample(landmask_pattern, bathy_Blondeau_native_agg, method = "near")
  # 
  # # Apply the mask to the native data
  # bathy_Blondeau_native_masked[landmask_pattern_native] <- NA
  
  
  
  
  # ################################## Merge 1: PR to MCD ##################################
  # 
  # # NOTE - the below 3 merges are for a "back-up" version of merged bathymetry that doesn't have the artifact
  # #         issues the Blondeau set does
  # 
  # 
  # #eliminate fill values at edge of MCD
  # bathy_crm_2019_agg_nofill <- ifel(bathy_crm_2019_agg < -100000, NA, bathy_crm_2019_agg)
  # 
  # bathy_merge1_crm = merge(bathy_crm_2019_agg_nofill, bathy_crm_2024_agg)
  # 
  # #plot briefly
  # bathy_merged_crm_reefdepth <- clamp(bathy_merge1_crm, lower = -50, upper = 0, values = TRUE) #limit depth to 0 m; eliminate land elevation
  # plot_extents = ext(xmin, xmax, ymin, ymax)
  # plot(bathy_merged_crm_reefdepth, 
  #      main="Merge #1",
  #      # col=hcl.colors(50, "Blues", rev=TRUE),
  #      # ext = e_pr, #e_crm,
  #      ext = plot_extents, #e_crm,
  #      legend=TRUE)
  # 
  # ################################## Merge 2: MCD to Coral Bay ##################################
  # 
  # #eliminate fill values at edge of MCD
  # bathy_crm_2019_agg_nofill <- ifel(bathy_crm_2019_agg < -100000, NA, bathy_crm_2019_agg)
  # 
  # bathy_crm_2024_agg_reefdepth_nofill <- ifel(bathy_crm_2024_agg_plot == 0, NA, bathy_crm_2024_agg_plot)
  # 
  # kiddel_point_lon = -64.71987
  # kiddel_point_lat = 18.30761
  # point_wgs84 <- st_sfc(st_point(c(kiddel_point_lon, kiddel_point_lat)), crs = 4326)
  # point_utm <- st_transform(point_wgs84, 26920)
  # kiddel_point_east_x <- st_coordinates(point_utm)[1]
  # cat("Converted longitude", kiddel_point_lon, "to UTM x-coordinate:", kiddel_point_east_x, "\n")
  # 
  # lon_rast <- init(bathy_crm_2024_agg_reefdepth_nofill, "x")
  # west_mask <- lon_rast >= kiddel_point_east_x
  # STJ_clip <- bathy_crm_2024_agg_reefdepth_nofill
  # STJ_clip[!west_mask] <- NA
  # 
  # lon_rast <- init(STJ_clip, "x")
  # east_non_na <- max(values(lon_rast)[!is.na(values(STJ_clip))], na.rm = TRUE)
  # east_trim_limit <- east_non_na - 1000  # trim 1 km west of easternmost valid value
  # east_mask <- lon_rast <= east_trim_limit
  # STJ_clip[!east_mask] <- NA
  # 
  # bathy_merge2_crm = merge(STJ_clip, bathy_merged_crm_reefdepth)
  # 
  # bathy_merge2_crm_nofill = ifel(bathy_merge2_crm == 0, NA, bathy_merge2_crm)
  # 
  # #plot briefly
  # bathy_merged2_crm_reefdepth <- clamp(bathy_merge2_crm_nofill, lower = -50, upper = 0, values = TRUE) #limit depth to 0 m; eliminate land elevation
  # plot_extents = ext(xmin, xmax, ymin, ymax)
  # plot(bathy_merged2_crm_reefdepth, 
  #      main="Merge #2",
  #      # col=hcl.colors(50, "Blues", rev=TRUE),
  #      # ext = e_pr, #e_crm,
  #      ext = plot_extents, #e_crm,
  #      legend=TRUE)
  # 
  # ################################## Merge 3: STX ##################################
  # 
  # VI_trough_lon = -64.767128
  # VI_trough_lat = 18.046230
  # point_wgs84 <- st_sfc(st_point(c(VI_trough_lon, VI_trough_lat)), crs = 4326)
  # point_utm <- st_transform(point_wgs84, 26920)
  # VI_trough_north_y <- st_coordinates(point_utm)[2]
  # cat("Converted latitude", VI_trough_lat, "to UTM y-coordinate:", VI_trough_north_y, "\n")
  # 
  # lat_rast <- init(bathy_crm_2019_agg_plot, "y")
  # south_mask <- lat_rast <= VI_trough_north_y
  # STX_clip <- bathy_crm_2019_agg_plot
  # STX_clip[!south_mask] <- NA
  # 
  # bathy_merge3_crm = merge(STX_clip, bathy_merged2_crm_reefdepth)
  # 
  # bathy_merge2_crm_nofill = ifel(bathy_merge3_crm == 0, NA, bathy_merge3_crm)
  # 
  # #plot briefly
  # bathy_merged3_crm_reefdepth <- clamp(bathy_merge2_crm_nofill, lower = -50, upper = 0, values = TRUE) #limit depth to 0 m; eliminate land elevation
  # plot_extents = ext(xmin, xmax, ymin, ymax)
  # plot(bathy_merged3_crm_reefdepth, 
  #      main="Merge #3",
  #      # col=hcl.colors(50, "Blues", rev=TRUE),
  #      # ext = e_pr, #e_crm,
  #      ext = plot_extents, #e_crm,
  #      legend=TRUE)
  # 
  # # #compare with 2019 bathy
  # # bathy_crm_2019_agg_reefdepth_nofill = ifel(bathy_crm_2019_agg_plot == 0, NA, bathy_crm_2019_agg_plot)
  # # plot(bathy_crm_2019_agg_reefdepth_nofill, #this shows that the 2019 crm is the clear winner, for MCD
  # #      main = 'CRM 2019',
  # #      # ext = e_pr,
  # #      ext = plot_extents, #e_crm,
  # #      legend = TRUE)
  # 
  ################################## Save objects/workspace  ##################################
  
  # #remove raster files with very large memory which don't work well with saving and re-loading downstream
  # # NOTE - can return to this if direct access to PR East is required!
  # rm(bathy_PR_East_clipped)
  # rm(bathy_PR_East)
  
  #save terra objects #and then workspace for use in downstream scripts
  save_spat_objects(output_dir = 'output/output_import_merge_rasters_higher-res/') #call from functions.R
  
  # Get all non-spatial objects
  non_spatial <- ls()[!sapply(ls(), function(x) inherits(get(x), c("SpatRaster", "SpatVector", "SpatExtent")))]
  
  # Save only non-spatial objects
  # NOTE - this helps with avoiding 'pointer' warnings/errors when loading everything again downstream
  save(list = non_spatial, file = here('output', 'output_import_merge_rasters_higher-res/import_merge_rasters_workspace.RData'))
  