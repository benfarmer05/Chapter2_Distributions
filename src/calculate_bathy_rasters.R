  
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
  
  # load_spat_objects(directory = here("output")) #call function
  # load(here("output/create_habitat_grid_workspace.RData")) #load workspace from upstream script
  
  #load initial set of objects
  load_spat_objects(directory = 'output/output_import_merge_rasters_higher-res/')
  load(here('output', 'output_import_merge_rasters_higher-res/import_merge_rasters_workspace.RData'))
  
  #load next set
  load_spat_objects(directory = 'output/output_create_habitat_grid/')
  load(here('output', 'output_create_habitat_grid/create_habitat_grid_workspace.RData'))
  
  ################################## Pretty bathy plot ##################################
  
  # color_palette <- scico(100, palette = "hawaii", direction = -1) #scico palette
  color_palette <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(100) #YlGnBu palette, as used by me in QGIS. color ramp reversed, and continuous colors interpolated
  
  # Plot the raster
  plot(bathy_merged3_crm_reefdepth,
       col = color_palette,
       # zlim = c(-50, 0),
       main = "Bathymetry (50m Resolution)",
       legend = TRUE)
  
  # # redacted - projected earlier on. can return to this if need to
  # # Reproject to NAD83 (geographic)
  # bathy_merged_geo <- project(bathy_merged3_crm_reefdepth, common_crs)
  # 
  # #note the slight shift in orientation
  # # NOTE - after looking at this, I think maybe we SHOULD project earlier in the pipeline. mesophotic ridge loses
  # #         significant detail
  # plot(bathy_merged_geo,
  #      col = color_palette,
  #      # zlim = c(-50, 0), #redundant
  #      main = "Bathymetry (50m Resolution)",
  #      legend = TRUE)
  
  #plot just north of STT to verify clamping function above worked correctly
  #
  # Define the extent for the region north of St. Thomas up to latitude 19°
  extent_area <- ext(c(-65.1, -64.75, 18.25, 18.6))
  
  # Crop the raster to the defined extent
  bathy_cropped <- crop(bathy_merged3_crm_reefdepth, extent_area)
  
  # Define a color palette for the plot
  color_palette <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(100)
  
  # Plot the cropped raster
  plot(bathy_cropped,
       col = color_palette,
       main = "Bathymetry North of St. Thomas Up to Latitude 19°",
       legend = TRUE)
  
  #verify that the clamp worked and introduced NAs
  values <- values(bathy_cropped)
  deeper_than_50m <- any(values < -50, na.rm = TRUE)
  
  if (deeper_than_50m) {
    cat("There are values deeper than 50 meters in the raster.")
  } else {
    cat("All values are within the specified range (i.e., shallower than 50 meters).")
  }

  na_values <- any(is.na(values))
  if (na_values) {
    cat("There are NA values in the raster.")
  } else {
    cat("There are no NA values in the raster.")
  }
  
  # Create a mask raster where NA values are set to 1 and all other values are set to 0
  na_mask <- is.na(bathy_cropped)
  
  # Define a color palette for the mask
  mask_palette <- c("transparent", "red")  # Transparent for non-NA, red for NA
  
  # Plot the NA mask raster
  plot(na_mask,
       col = mask_palette,
       main = "NA Mask of Bathymetry Raster",
       legend = FALSE)  # No legend for binary mask
  
  # Create a binary raster where values shallower than 50 meters are set to 1 and others are set to 0
  deep_mask <- bathy_cropped
  deep_mask[] <- ifelse(values(bathy_cropped) < -50, 1, 0)
  
  # Define a color palette for the binary raster
  deep_palette <- c("transparent", "blue")  # Transparent for shallow areas, blue for deep areas. NOTE - maybe reverse?
  
  # Plot the binary raster for values deeper than 50 meters
  plot(deep_mask,
       col = deep_palette,
       main = "Binary Mask for Depths Deeper Than 50 Meters",
       legend = FALSE)  # No legend for binary mask

  bathy_value <- extract(bathy_cropped, cbind(-64.95, 18.45))

  # #ggplot method, with tidyterra. this maybe didn't work because of the size of the rasters, so I went with tmap & rayshader below
  #
  # # # Convert raster to a data frame for ggplot2
  # # bathy_raster <- rast(here("output", "bathy_50m.tif"))
  # # bathy_df <- as_tibble(bathy_raster, xy = TRUE)
  # bathy_df <- as_tibble(merged_bathy, xy = TRUE)
  #
  # # Plot raster data
  # ggplot(bathy_df) +
  #   geom_raster(aes(x = x, y = y, fill = layer)) +
  #   scale_fill_viridis_c(option = "D", limits = c(-50, 0)) + # Adjust limits to desired depth range
  #   labs(
  #     title = "Bathymetry",
  #     fill = "Depth (meters)"
  #   ) +
  #   coord_fixed() + # Fix aspect ratio to make sure x and y are scaled equally
  #   theme_minimal() +
  #   theme(
  #     axis.title = element_blank(),
  #     axis.text = element_text(size = 8),
  #     legend.title = element_text(size = 10),
  #     legend.text = element_text(size = 8)
  #   )
  #
  # ggplot(bathy_df) +
  #   geom_raster(aes(x = x, y = y))
  
  #tmap method
  tm_shape(bathy_merged3_crm_reefdepth) +
    tm_raster(style = 'cont', palette = color_palette) #this is outdated as of tmap v4, but works for continuous gradient
  # tm_shape(bathy_merged_50m) +
  #   tm_raster() +
  #   tm_scale_continuous(values = color_palette) #this doesn't seem right, but makes a discrete contoured depth gradient

  # bathy_map = tm_shape(merged_bathy) +
  #   tm_raster(style = 'cont')

  # output_file = file.path(here('output'), "bathy_50m.svg")
  # tmap_save(bathy_map, filename=output_file, height=8.5, width=11, units="in", dpi=300)
  # output_file = file.path(here('output'), "bathy_50m.pdf")
  # tmap_save(bathy_map, filename=output_file, height=8.5, width=11, units="in", dpi=300)
  # output_file = file.path(here('output'), "bathy_50m.jpeg")
  # tmap_save(bathy_map, filename=output_file, height=8.5, width=11, units="in", dpi=300)
  # output_file = file.path(here('output'), "bathy_50m.tiff")
  # tmap_save(bathy_map, filename=output_file, height=8.5, width=11, units="in", dpi=300)

  # # Save the merged raster
  # output_dir <- here("output")
  # output_file <- file.path(output_dir, "bathy_50m.tif")
  # writeRaster(merged_bathy, filename = output_file, overwrite = TRUE)

  #rayshader method
  # in 2D
  # elmat = raster_to_matrix(merged_bathy)
  # elmat %>%
  #   sphere_shade(texture = "desert") %>%
  #   plot_map()
  raster_to_matrix(bathy_merged3_crm_reefdepth) |> height_shade() |> plot_map()

  # # in 3D
  # elmat %>%
  #   # sphere_shade(texture = "desert") %>%
  #   # add_water(detect_water(elmat), color = "desert") %>%
  #   # add_shadow(ray_shade(elmat, zscale = 3), 0.5) %>%
  #   add_shadow(ambient_shade(elmat), 0) %>%
  #   plot_3d(elmat, zscale = 10, fov = 0, theta = 135, zoom = 0.75, phi = 45, windowsize = c(1000, 800))
  # Sys.sleep(0.2)
  # render_snapshot()
  
  ################################## Test resolution/artifacts ##################################
  
  # STOPPING POINT - 17 Oct 2024
  #   - A few issues
  #   - 1.) No matter whether I use 2-m or downsampled 50-m resolution bathy raster, all derived products have issues. some of these are
  #           major, like at the MCD and along coastlines. Should figure out if it would be easier to simply use VI_Shapes bathy or reach
  #           out to Katharine Egan to simply re-use what she had. But she also had a far more limited domain...probably not helpful.
  #           I guess just try and get QGIS working again and take a look at the geodatabase to see what I can do with it. My gut is that
  #           a combination of VI_Shapes and other bathy I can mosaic together from NOAA sources will be the best route forward
  #   - 2.) Once the above is figured out, a lot is in place, but I still unfortunately need to consider how much that 50-m downsampling
  #           negatively impacts our capabilities for inference. Egan (2021) I'm pretty sure downsampled before deriving (Check), but
  #           studies I really like, like Young and Carr, did not do this and I think benefited from it. Will have to see
  #     3.) Final thing is assembling all the parts in terms of predictor layers - that will take some serious effort
  
  # Calculate terrain attributes. may consider Benthic Terrain Modeler (but need Arc & potentially arcgisbinding package in R) or MultiscaleDTM package (https://cran.r-project.org/web/packages/MultiscaleDTM/readme/README.html)
  #
  slope_terra <- terrain(bathy_merged3_crm_reefdepth, v = "slope", unit = 'degrees')
  slopeofslope_terra = terrain(slope_terra, v = "slope", unit = 'degrees')
  aspect_terra <- terrain(bathy_merged3_crm_reefdepth, v = "aspect", unit = 'degrees')
  TPI_terra = terrain(bathy_merged3_crm_reefdepth, v = 'TPI') #Wilson 2007 bathymetric-friendly method
  # flowdir <- terrain(bathy_merged3_crm_reefdepth, v = "flowdir")
  roughness = terrain(bathy_merged3_crm_reefdepth, v = "roughness")
  TRI = terrain(bathy_merged3_crm_reefdepth, v = 'TRI')
  
  # TPI_spatialeco = tpi(bathy_merged3_crm_reefdepth)
  VRM <- vrm(bathy_merged3_crm_reefdepth)
  # TRI_spatialeco <- tri(bathy_merged3_crm_reefdepth, s = 3)
  totalcurv = curvature(bathy_merged3_crm_reefdepth, type = 'total')
  # planformcurv = curvature(bathy_merged3_crm_reefdepth, type = 'planform') #takes too long to run
  # profilecurv = curvature(bathy_merged3_crm_reefdepth, type = 'profile') #takes too long to run
  
  # consider geodiv! (can get fractal dimension)
  
  # slope_multiscale <- SlpAsp(r = bathy_merged3_crm_reefdepth, w = 3, metrics = c("slope")) #na.rm = TRUE could be used
  # slopeofslope_multiscale <- SlpAsp(r = slope_multiscale, w = 3, metrics = c("slope"))
  # aspect_multiscale <- SlpAsp(r = bathy_merged3_crm_reefdepth, w = 3, metrics = c("aspect"))
  # TPI_multiscale <- TPI(bathy_merged3_crm_reefdepth, w = c(3, 3))
  # BPI_multiscale <- BPI(bathy_merged3_crm_reefdepth, w = c(1, 3))  # BPI is similar to TPI for bathymetry
  # DMV <- DMV(bathy_merged3_crm_reefdepth, w = c(3, 3))
  # RelPos <- RelPos(bathy_merged3_crm_reefdepth, w = c(5, 5))
  # SAPA <- SAPA(bathy_merged3_crm_reefdepth, w = c(5, 5))
  VRM_multiscale <- VRM(bathy_merged3_crm_reefdepth, w = c(5, 5))
  # surfarea = SurfaceArea(bathy_merged3_crm_reefdepth)
  maxcurv_multiscale = Qfit(r = bathy_merged3_crm_reefdepth, w = 3, metrics = 'maxc')
  meancurv_multiscale = Qfit(r = bathy_merged3_crm_reefdepth, w = 3, metrics = 'meanc')
  planformcurv_multiscale = Qfit(r = bathy_merged3_crm_reefdepth, w = 3, metrics = 'planc')
  profilecurv_multiscale = Qfit(r = bathy_merged3_crm_reefdepth, w = 3, metrics = 'profc')
  # depth = Qfit(r = bathy_merged3_crm_reefdepth, metrics = 'elev')
  
  #spatialeco for curvature?
  
  
  # NOTE - I think slope, slope of slope, aspect, and TPI, roughness, TRI, and
  #         VRM are the variables I am likely moving forward with
  
  # # TEST
  # 
  # bathy_crm_2019_clipped
  # test = clamp(bathy_crm_2019_clipped, lower = -50, upper = 0, values = TRUE)
  # plot(test)
  # 
  # slope_test <- terrain(test, v = "slope", unit = 'degrees')
  # slopeofslope_test = terrain(slope, v = "slope", unit = 'degrees')
  # aspect_test <- terrain(test, v = "aspect", unit = 'degrees')
  # flowdir_test <- terrain(test, v = "flowdir")
  # roughness_test = terrain(test, v = "roughness")
  # TPI_test = terrain(test, v = 'TPI') #Wilson 2007 bathymetric-friendly method
  # TRI_test = terrain(test, v = 'TRI')
  # VRM_test <- vrm(test)
  # VRM_largest_test <- vrm(test, s = 9)
  # TPI_multiscale_test <- TPI(test, w = c(5, 5))
  # DMV_test <- DMV(test, w = c(5, 5))
  # BPI_test <- BPI(test, w = c(5, 5))  # BPI is same as TPI for bathymetry
  # RelPos_test <- RelPos(test, w = c(5, 5))
  # SAPA_150m_test <- SAPA(test, w = c(3, 3))
  # SAPA_250m_test <- SAPA(test, w = c(5, 5))  # Surface Area to Planar Area ratio
  # SAPA_350m_test <- SAPA(test, w = c(7, 7))
  # VRM_multiscale_test <- VRM(test, w = c(5, 5))
  # 
  # 
  # 
  # plot(test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(slope_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(slopeofslope_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(aspect_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(flowdir_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(roughness_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(TPI_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(TRI_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(VRM_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(VRM_largest_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(TPI_multiscale_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(DMV_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(BPI_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(RelPos_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(SAPA_150m_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(SAPA_250m_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(SAPA_350m_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(VRM_multiscale_test, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # 
  # # TEST
  
  # plot(bathy_merged3_crm_reefdepth, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(slope, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(slopeofslope, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(aspect, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(flowdir, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(roughness, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(TPI, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(TRI, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(VRM, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(VRM_largest, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(TPI_spatialeco, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(TPI_multiscale, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(DMV, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(BPI, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(RelPos, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(SAPA_150m, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(SAPA_250m, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(SAPA_350m, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(VRM_multiscale, xlim = c(270000, 330000), ylim = c(2010000, 2040000))  
  
  # p1 <- ggplot() +
  #   geom_spatraster(data = TPI_spatialeco, maxcell = Inf) +
  #   scale_fill_cmocean(name = "balance", na.value = "transparent") +
  #   coord_sf(xlim = c(270000, 330000), ylim = c(2010000, 2040000)) +
  #   theme_minimal() +
  #   labs(title = "Terra Raster with cmocean 'balance'")
  
  
  plot(slope_terra, 
       col = cmocean("deep")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  plot(slopeofslope_terra, 
       col = cmocean("deep")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  plot(roughness, 
       col = cmocean("deep")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  VRM <- clamp(VRM, lower = -0.005, upper = 0.005)
  plot(VRM, 
       col = cmocean("balance")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  # Clamp values to -5 to 5 range for colorbar
  TPI_spatialeco <- clamp(TPI_spatialeco, lower = -5, upper = 5)
  plot(TPI_spatialeco, 
       col = cmocean("balance")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  TPI_multiscale <- clamp(TPI_multiscale, lower = -5, upper = 5)
  plot(TPI_multiscale, 
       col = cmocean("balance")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  TPI_terra <- clamp(TPI_terra, lower = -5, upper = 5)
  plot(TPI_terra, 
       col = cmocean("balance")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  BPI_multiscale <- clamp(BPI_multiscale, lower = -10, upper = 10)
  plot(BPI_multiscale, 
       col = cmocean("balance")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  DMV <- clamp(DMV, lower = -5, upper = 5)
  plot(DMV, 
       col = cmocean("balance")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  # totalcurv <- clamp(totalcurv, lower = -0.006, upper = 0.006)
  # plot(totalcurv, 
  #      col = cmocean("balance")(100),
  #      xlim = c(270000, 330000),
  #      ylim = c(2010000, 2040000),
  #      main = "Terra Raster with cmocean balance")
  
  planformcurv_multiscale <- clamp(planformcurv_multiscale, lower = -0.004, upper = 0.004)
  plot(planformcurv_multiscale, 
       col = cmocean("balance")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  profilecurv_multiscale <- clamp(profilecurv_multiscale, lower = -0.01, upper = 0.01)
  plot(profilecurv_multiscale, 
       col = cmocean("balance")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  meancurv_multiscale <- clamp(meancurv_multiscale, lower = -0.004, upper = 0.004)
  plot(meancurv_multiscale, 
       col = cmocean("balance")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  maxcurv_multiscale <- clamp(maxcurv_multiscale, lower = -0.002, upper = 0.002)
  plot(maxcurv_multiscale, 
       col = cmocean("balance")(100),
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  depth_OG <- clamp(bathy_merged3_crm_reefdepth, lower = -50, upper = 0)
  plot(depth_OG, 
       # col = cmocean("deep")(100),
       col = rev(cmocean("deep")(100)),  # Add rev() to reverse
       xlim = c(270000, 330000),
       ylim = c(2010000, 2040000),
       main = "Terra Raster with cmocean balance")
  
  # plot(slope, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(slopeofslope, xlim = c(280000, 330000), ylim = c(2010000, 2040000))
  # plot(aspect, xlim = c(280000, 330000), ylim = c(2010000, 2040000))
  # plot(flowdir, xlim = c(280000, 330000), ylim = c(2010000, 2040000))
  # plot(roughness, xlim = c(280000, 330000), ylim = c(2010000, 2040000))
  # plot(TPI, xlim = c(280000, 330000), ylim = c(2010000, 2040000))
  # plot(TRI, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(VRM, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(VRM_largest, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(TPI_multiscale, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(DMV, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(BPI, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(RelPos, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(SAPA_150m, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(SAPA_250m, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(SAPA_350m, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  # plot(VRM_multiscale, xlim = c(270000, 330000), ylim = c(2010000, 2040000))
  
  # NOTE - strikes me that BPI, slope of slope, roughness might be the most useful
  #         - remains to be seen whether these are worth calculating at native resolution. probably not, due to 
  #           our contiguous bathy product being a mosaic of resolutions. 50 m seems like a good compromise. but
  #           could consider 30 m (and possibly then resample to 50). move forward with this though!
  #         - SAPA in particular might become a lot more useful at higher res
  #         - BPI might be a good indicator of risk (or lack thereof) to sedimentation from upslope
  #         - roughness seems like a good metric, just unclear what it correlates to in the literature at the moment
  
  # STOPPING POINT - 20 June 2025
  #   - now just getting all the spatial layers figured out. need to have a grid, then snap the biological layers
  #       to it. then a matter of running the GAMs, and making predictions
  
  # #TEST
  # bathy_merged_crm_forslopes <- clamp(bathy_merge1_crm, lower = -50, upper = 0, values = TRUE)
  # # slope_merge1 = terrain(bathy_merge1_crm, v = "slope", unit = 'degrees') # THIS LOOKS CRAZY
  # slope_merge1 = terrain(bathy_merged_crm_forslopes, v = "slope", unit = 'degrees')
  # slopeofslope_merge1 = terrain(slope_merge1, v = "slope", unit = 'degrees')
  # plot(slope_merge1, col = color_palette, main = "Slope (degrees)", legend = TRUE)
  # plot(slopeofslope_merge1, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  # #TEST
  
  # #derive bathymetry products from original resolution (2 m) raster. clipped to just STTSTJ since even that takes a long time
  # slope_STTSTJ = terrain(bathy_STTSTJ, v = "slope", unit = 'degrees')
  # slopeofslope_STTSTJ = terrain(slope_STTSTJ, v = "slope", unit = 'degrees')
  # slope_PR_East = terrain(bathy_PR_East, v = "slope", unit = 'degrees')
  # 
  # color_palette <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(100) #YlGnBu palette, as used by me in QGIS. color ramp reversed, and continuous colors interpolated
  # plot(slope_STTSTJ, col = color_palette, main = "Slope (degrees)", legend = TRUE)
  # plot(slopeofslope_STTSTJ, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  # plot(slope_PR_East, col = color_palette, main = "Slope (degrees)", legend = TRUE)
  # 
  # 
  # slope_STTSTJ_50m = terrain(bathy_STTSTJ_agg, v = "slope", unit = 'degrees')
  # slopeofslope_STTSTJ_50m = terrain(slope_STTSTJ_50m, v = "slope", unit = 'degrees')
  # slope_PR_East_50m = terrain(bathy_PR_East_agg, v = "slope", unit = 'degrees')
  # slopeofslope_PR_East_50m = terrain(slope_PR_East_50m, v = "slope", unit = 'degrees')
  # 
  # plot(slope_STTSTJ_50m, col = color_palette, main = "Slope (degrees)", legend = TRUE)
  # plot(slopeofslope_STTSTJ_50m, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  # plot(slope_PR_East_50m, col = color_palette, main = "Slope (degrees)", legend = TRUE)
  # plot(slopeofslope_PR_East_50m, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  # plot(bathy_PR_East, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  # plot(bathy_PR_East_agg, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  # plot(bathy_STTSTJ, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  # plot(bathy_STTSTJ_agg, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  
  # STOPPING POINT - 17 OCT 2024
  #   - as sad as it is, I might need to seriously consider simply downscaling everything to a common resolution...which I was pretty much
  #       doing from the start already with the 50 m thing. ironically that might be perfect since it avoids totally losing slope and slope
  #       of slope data where resolution drops off. it's a bummer but I think I should probably consider this
  #         - NOTE: this is the route I am going for now. will do 50 m resolution as starting point for all downstream analyses and go from
  #                 there
  
  # # Roughness, TPI, TRI
  # roughness <- focal(bathy_merged, w = matrix(1, nrow = 3, ncol = 3), fun = function(x) max(x) - min(x))
  # tpi <- focal(bathy_merged, w = matrix(1, nrow = 3, ncol = 3), fun = function(x) x[5] - mean(x[-5]))
  # tri <- focal(bathy_merged, w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sum(abs(x[-5] - x[5])) / 8)
  
  # Reproject rasters to geographic CRS
  slope_geo <- project(slope, common_crs)
  slopeofslope_geo <- project(slopeofslope, common_crs)
  aspect_geo <- project(aspect, common_crs)
  flowdir_geo <- project(flowdir, common_crs)
  roughness_geo <- project(roughness, common_crs)
  tpi_geo <- project(TPI, common_crs)
  tri_geo <- project(TRI, common_crs)
  
  # Plot the results
  plot(bathy_merged3_crm_reefdepth, col = color_palette, main = "Depth (m)", legend = TRUE)
  plot(slope, col = color_palette, main = "Slope (degrees)", legend = TRUE)
  plot(slopeofslope, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE)
  plot(aspect, col = color_palette, main = "Aspect (degrees)", legend = TRUE)
  plot(flowdir, col = color_palette, main = "Flow Direction", legend = TRUE)
  plot(roughness, col = color_palette, main = "Roughness", legend = TRUE)
  plot(TPI, col = color_palette, main = "TPI", legend = TRUE)
  plot(TRI, col = color_palette, main = "TRI", legend = TRUE)
  
  # # Define the extent for the area between St. Thomas and St. Croix
  # # Adjust these coordinates as needed for your specific area of interest
  # extent_STT_STX <- ext(c(-65.1, -64.4, 17.6, 18.4))
  # 
  # # plot habitat grids
  # #50m NOAA sampleframe
  # ggplot() +
  #   geom_polygon(data = df.STTSTJ_grid, aes(x = long, y = lat, group = group,),
  #                fill = "green", color = "black", lwd = 1) +
  #   theme_bw()
  # 
  # #my custom 650m grid
  # grid_650m = readOGR(dsn = ".", layer = "STTSTJ_meso_merged_testforR")
  # #reprojection
  # grid_650m = spTransform(grid_650m, crs(STTSTJ_grid))
  # # df.grid_650 = fortify(grid_650m)
  # ggplot() +
  #   geom_polygon(data = df.grid_650, aes(x = long, y = lat, group = group,),
  #                fill = "pink", color = "black", lwd = 1) +
  #   theme_bw() #+
  # xlim(-65.2, -64.6) +
  #   ylim(18.15, 18.45)
  # 
  # 
  # ### testing local and global CRS ###
  # # NOTE: 8 NOV 2022: Need to repeat the below with all 2021 grids (that come from NCRMP package) in order to get all grids in the correct local CRS.
  # STTSTJ_grid_2019 = readOGR(dsn = ".", layer = "STTSTJ_2019_sampleframe_localCRS-forR")
  # crs(STTSTJ_grid_2019)
  # crs(STTSTJ_grid)
  # 
  # #reprojection
  # STTSTJ_grid_2021_localCRS = spTransform(STTSTJ_grid, crs(STTSTJ_grid_2019))
  # shapefile(x = STTSTJ_grid_2021_localCRS, file = "/Users/benja/Documents/Carib_Habitat/Carib_Habitat_QGIS/Sampleframes/STTSTJ_grid_2021_localCRS.shp")
  

  ################################## Export bathymetry ##################################
  
  # STOPPING POINT - 18 May 2025
  #   - here is where I need to be taking the merged bathy and using it to refine a rough habitat grid
  #   - could literally be as simple as taking the <50 m binary mask and producing a properly projected grid over it
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
  
  # #save terra objects and then workspace for use in downstream scripts
  # save_spat_objects() #call from functions.R
  # save.image(file = here("output", "calculate_bathy_rasters_workspace.RData"))
  