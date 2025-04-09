  
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
  source(here("src/functions.R"))
  
  # # Print GDAL, PROJ, and GEOS versions linked with sf
  # sf::sf_extSoftVersion()
  # 
  # # Print GDAL version linked with terra
  # terra::gdal()
  
  load(here("output/import_merge_rasters_workspace.RData")) #load workspace from upstream script
  load_spat_objects(directory = here("output")) #call function
  
  # # all various things I tried to get things working. not currently required on my M1 Macbook, newly updated
  # install.packages("sf", dependencies = TRUE) #this may be required after updating to Sequoia (I reinstalled R via homebrew due to issues with Sequoia and it broke some more things)
  # tools::package_dependencies(c("sf", "here", "terra", "tidyterra", "ggplot2", "tmap", "rayshader", "scico", "RColorBrewer"), recursive = TRUE) #check all dependencies for below packages
  # install.packages(c("sf", "here", "terra", "tidyterra", "ggplot2", "tmap", "rayshader", "scico", "RColorBrewer"), dependencies = TRUE) #more dependency installs, if required
  # Sys.getenv("LD_LIBRARY_PATH")
  # Sys.setenv(LD_LIBRARY_PATH = "/opt/homebrew/lib:/opt/homebrew/include:/opt/homebrew/share") #this may also be required, followed by possibly installing sf and terra from source (below)
  # Sys.setenv(PROJ_LIB = "/opt/homebrew/share/proj") #this was temporarily needed after updating M1 Macbook to Sequoia OS
  # Sys.setenv(GDAL_DATA = "/opt/homebrew/share/gdal") #this was temporarily needed after updating M1 Macbook to Sequoia OS
  
  # # What I had to do to get R spatial stuff working after updating to Sequoia was:
  # #   1.) Tried uninstalling everything macports, which was probably dumb because now QGIS has to be entirely reinstalled (and may
  # #       may not work now anyways through macports yet because of delays in dependency updates after Sequoia - we shall see)
  # #        - Update Oct 2024: looks like QGIS is working fine after installing the launcher and also macports version of it
  # #   2.) Tried installing R through homebrew, this was a cluster and most things installed poorly or not at all
  # #   3.) Tried installing R again through regular Mac ARM, and then installing spatial packages from source as below. this seems to work
  # # Install from source, forcing R to use the latest versions
  # install.packages("sf", type = "source") #this might not work yet for Sequoia? also maybe doesn't matter if the code runs I guess. for now
  # install.packages("terra", type = "source")
  # #check current R version
  # R.version.string
  
  ################################## Test simple crm-usvi ##################################
  
  
  
  # STOPPING POINT - 8 April 2025
  #   - okay new plan!! this might work!!
  #       - PR comes from crm_USVI which I pulled off NOAA servers in 2024
  #       - STT/STJ/STX come from crm_USVI which Dan pulled off NOAA servers in (?? 2019?)
  #       - Chose this because these are the best combination of resolution, no major artifacts, and extent
  #       - Will just have to figure out how to properly merge them together. this will involve clipping the 2024 crm at the PR/STT edge
  #
  #   - okay, end of day update and it's actually even more complicated than I thought
  #       - this dataviewer (https://www.ncei.noaa.gov/maps/bathymetry/?layers=nos_hydro&minx=-65.5684&maxx=-64.4215&miny=17.8471&maxy=18.906)
  #           shows CUDEM vs CRM, 3, 10, 30, and 90 m resolutions.
  #       - 30 m and 90 m res CRM are very solid, and also comprehensive (other than Anegada). I would simply use them for the entire
  #           domain, except for some slight problems:
  #         - Dan's crm (where/when did it come from ???) actually has much better resolution, and less artifacts, in the MCD & north drop
  #         - BUT, the recent 30 m crm has a bit better res south of STJ - sort of? toss up there
  #         - but DEFINITELY, the recent 30 m crm has far better resolution in Coral Bay & Haulover (STJ). also east end bay Tortola. but mosaicing artifacts introduced at east end
  #           - question is - is that resolution at Coral Bay / Haulover / East End Bay worth it?? would mean stitching crm-2024 for PR to MCD, with crm-2019 for MCD to Coral Bay, with crm-2024 for Coral Bay to Anegada
  #           - is that 3-part stitch even possible? maybe after merging and comprehensive resampling?
  #           - try it! if that seems impossible though, simply stitch crm-2024 for PR with crm-2019 for USVI
  #           - and worst case, even that doesn't work well, so I just straight up use crm-2024 for everything (artifacts and all). maybe just with Anegada tacked on with crm-2019
  #           - remember that, absolute worst case, I can always through out "training" points at NCRMP / DCRMP points that are over problematic sections of bathymetry
  
  
  #pull NOAA crm, produced in ~2019 (maybe ? from Dan Holstein) from presumably NOAA bathymetry database (https://www.ncei.noaa.gov/maps/grid-extract/)
  # - CRS: 26920 (NAD83 / UTM zone 20N)
  # - max 30 m resolution
  # This will be used for USVI, since it has nice resolution generally, less artifacts than the crm below, and full extent across Anegada
  crm_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Holstein_VI_Shapes/VI_Shapes/Bathy/crm_usvi.tif"
  describe(crm_path)
  bathy_crm_2019 = rast(crm_path)
  res(bathy_crm_2019)
  crs(bathy_crm_2019)
  crs(bathy_PR_East)
  crs(bathy_STTSTJ)
  
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
  # projected_crs <- st_crs(32620)  # WGS 84 / UTM zone 20N (suitable for Puerto Rico and the Virgin Islands). projected CRS
  # projected_crs <- st_crs(26920)  # NAD83 / UTM Zone 20N (suitable for Puerto Rico and the Virgin Islands). projected CRS
  projected_crs <- "EPSG:26920"  # NAD83 / UTM Zone 20N
  bathy_crm_2024 <- project(bathy_crm_2024, projected_crs)
  crs(bathy_crm_2024)
  crs(bathy_crm_2019)
  res(bathy_crm_2019)
  res(bathy_crm_2024)
  
  #read the polygon release grid (substrate for downstream hydrodynamic connectivity simulations) shapefile. CRS: 26920 (NAD83 / UTM zone 20N)
  # reef_polys <- st_read("/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/output/Habitat_Grid/Reef_Polygons/polys_apr2025_operational.shp")
  reef_polys = vect("/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/output/Habitat_Grid/Reef_Polygons/polys_apr2025_operational.shp")
  
  #clip crm bathymetry to the extent of release grid
  # NOTE - will also need to clip out overlap with STTSTJ raster - otherwise, weird artifacts over MCD are brought over in the merge
  # bathy_PR_East_clipped <- crop(bathy_PR_East, reef_polys)
  # bathy_PR_east_agg_clipped <- crop(bathy_PR_East_agg, reef_polys)
  bathy_crm_2024_clipped <- crop(bathy_crm_2024, reef_polys)
  
  #retrieve and apply crm raster extents to fresh raster template, then resample to 50 m resolution using template  # NOTE -
  # NOTE - should resampling be done before or after merging with PR East?
  #      - should I actually be aggregating to a masked / cookie-cutter'd grid which is flush with the 50-m NCRMP sampling grid?
  #         - answer: I think no, since the Puerto Rico grid is a different projection and the different grids would never all align anyways
  e_crm <- ext(bathy_crm)
  e_pr <- ext(bathy_PR_East_clipped)
  xmin_combined <- min(e_crm$xmin, e_pr$xmin)
  xmax_combined <- max(e_crm$xmax, e_pr$xmax)
  ymin_combined <- min(e_crm$ymin, e_pr$ymin)
  ymax_combined <- max(e_crm$ymax, e_pr$ymax)
  xmin <- floor(xmin_combined / 50) * 50
  # xmax <- ceiling(xmax_combined / 50) * 50
  xmax = 385000 #380000 #manual edit to greatly cut down the size of the raster over deep ocean
  # ymin <- floor(ymin_combined / 50) * 50
  # ymax <- ceiling(ymax_combined / 50) * 50
  ymin = 1945000 #1940000 #manual edit to greatly cut down the size of the raster over deep ocean
  ymax = 2087500 #2080000 #manual edit to greatly cut down the size of the raster over deep ocean
  new_ext <- ext(xmin, xmax, ymin, ymax)
  template_raster <- rast(new_ext, resolution = 50, crs = crs(bathy_crm))
  bathy_crm_agg <- resample(bathy_crm, template_raster, method = "average") #NOTE - resampling was done because 'aggregate' could not produce exact discrete 50 x 50 m resolution
  bathy_PR_East_agg <- resample(bathy_PR_East_clipped, template_raster, method = "average")
  
  res(bathy_crm_agg)  # Should be exactly [1] 50 50
  res(bathy_PR_East_agg)  # Should be exactly [1] 50 50
  
  #merge PR east with crm bathy - see how it looks
  # bathy_merged_crm = merge(bathy_PR_East, bathy_crm)
  bathy_merged_crm = merge(bathy_PR_East_agg, bathy_crm_agg)
  
  # NOTE - PR bathy instead has to come from Dan's 'VI_shapes' from 2019/2020 - it has less artifacts which would cause a lot of issues if introduced
  
  ################################## Set-up ##################################
  
  # After loading the workspace, re-write the SpatRasters from .tif (required because of the way terra works with R objects, I think)
  # bathy_merged = rast(here("output", "bathy_50m.tif"))
  # bathy_STTSTJ_agg <- rast(here("output", "bathy_STTSTJ_agg.tif"))
  # bathy_STX_agg <- rast(here("output", "bathy_STX_agg.tif"))
  # bathy_PR_East_agg <- rast(here("output", "bathy_PR_East_agg.tif"))
  # bathy_PR_South_agg <- rast(here("output", "bathy_PR_South_agg.tif"))
  # bathy_PR_West_agg <- rast(here("output", "bathy_PR_West_agg.tif"))
  # bathy_PR_North_agg <- rast(here("output", "bathy_PR_North_agg.tif"))
  # bathy_merged = readRDS(here("output", "bathy_50m.rds"))
  # bathy_STTSTJ_agg <- readRDS(here("output", "bathy_STTSTJ_agg.rds"))
  # bathy_STX_agg <- readRDS(here("output", "bathy_STX_agg.rds"))
  # bathy_PR_East_agg <- readRDS(here("output", "bathy_PR_East_agg.rds"))
  # bathy_PR_South_agg <- readRDS(here("output", "bathy_PR_South_agg.rds"))
  # bathy_PR_West_agg <- readRDS(here("output", "bathy_PR_West_agg.rds"))
  # bathy_PR_North_agg <- readRDS(here("output", "bathy_PR_North_agg.rds"))
  
  res(bathy_STTSTJ)
  res(bathy_merged_50m)
  res(bathy_crm)
  res(bathy_merged_crm)
  
  # Clip the raster to include only depths of 50 meters and shallower (this '50' is independent of the resolution of the bathy raster,
  #   which just happens to be 50 m horizontal)
  # bathy_merged <- clamp(bathy_merged, lower = -50, upper = 0, values = TRUE) #clamp does not introduce NAs, rather fills in -50's
  bathy_merged_50m <- ifel(bathy_merged_50m < -50, NA, bathy_merged_50m)
  bathy_merged_crm = ifel(bathy_merged_crm < -50 | bathy_merged_crm > 0, NA, bathy_merged_crm) #also set land as NA, since this isn't done by default (land is positive elevation DEM here)
  
  #verify that the above worked and introduced NAs
  values <- values(bathy_merged_50m)
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
  
  values <- values(bathy_merged_crm)
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
  na_mask <- is.na(bathy_merged_50m)
  
  # Define a color palette for the mask
  mask_palette <- c("transparent", "red")  # Transparent for non-NA, red for NA
  
  # # Plot the NA mask raster
  # plot(na_mask,
  #      col = mask_palette,
  #      main = "NA Mask of Bathymetry Raster",
  #      legend = FALSE)  # No legend for binary mask
  
  # Create a binary raster where values shallower than 50 meters are set to 1 and others are set to 0
  deep_mask <- bathy_merged_50m
  deep_mask[] <- ifelse(values(bathy_merged_50m) > -50, 1, 0)
  deep_mask_crm <- bathy_merged_crm
  deep_mask_crm[] <- ifelse(values(bathy_merged_crm) > -50, 1, 0)

  # Define a color palette for the binary raster
  deep_palette <- c("transparent", "blue")  # Transparent for shallow areas, blue for deep areas. NOTE - I think? wait maybe reverse
  
  # Plot the binary raster for values deeper than 50 meters
  plot(deep_mask,
       col = deep_palette,
       main = "Binary Mask for Depths Deeper Than 50 Meters",
       legend = FALSE)  # No legend for binary mask
  plot(deep_mask_crm,
       col = deep_palette,
       main = "Binary Mask for Depths Deeper Than 50 Meters",
       legend = FALSE)  # No legend for binary mask
  
  ################################## Plotting ##################################
  
  # color_palette <- scico(100, palette = "hawaii", direction = -1) #scico pallete
  color_palette <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(100) #YlGnBu palette, as used by me in QGIS. color ramp reversed, and continuous colors interpolated

  # Plot the raster
  plot(bathy_merged_50m,
       col = color_palette,
       # zlim = c(-50, 0),
       main = "Bathymetry (50m Resolution)",
       legend = TRUE)
  
  # Reproject to NAD83 (geographic)
  bathy_merged_geo <- project(bathy_merged_50m, common_crs)
  
  #note the slight shift in orientation
  plot(bathy_merged_geo,
       col = color_palette,
       # zlim = c(-50, 0), #redundant
       main = "Bathymetry (50m Resolution)",
       legend = TRUE)

  #plot just north of STT to verify clamping function above worked correctly
  #
  # Define the extent for the region north of St. Thomas up to latitude 19°
  extent_area <- ext(c(-65.1, -64.75, 18.25, 18.6))

  # Crop the raster to the defined extent
  bathy_cropped <- crop(bathy_merged_geo, extent_area)

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
  tm_shape(bathy_merged_50m) +
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
  raster_to_matrix(bathy_merged_50m) |> height_shade() |> plot_map()

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
  #derive bathymetry products from already-downsampled (2 m --> 5 m resolution) raster
  slope <- terrain(bathy_merged_50m, v = "slope", unit = 'degrees')
  slopeofslope = terrain(slope, v = "slope", unit = 'degrees')
  aspect <- terrain(bathy_merged_50m, v = "aspect", unit = 'degrees')
  flowdir <- terrain(bathy_merged_50m, v = "flowdir")
  roughness = terrain(bathy_merged_50m, v = "roughness")
  TPI = terrain(bathy_merged_50m, v = 'TRI') #Wilson 2007 bathymetric-friendly method
  TRI = terrain(bathy_merged_50m, v = 'TPI')
  
  #derive bathymetry products from original resolution (2 m) raster. clipped to just STTSTJ since even that takes a long time
  slope_STTSTJ = terrain(bathy_STTSTJ, v = "slope", unit = 'degrees')
  slopeofslope_STTSTJ = terrain(slope_STTSTJ, v = "slope", unit = 'degrees')
  slope_PR_East = terrain(bathy_PR_East, v = "slope", unit = 'degrees')
  
  color_palette <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(100) #YlGnBu palette, as used by me in QGIS. color ramp reversed, and continuous colors interpolated
  plot(slope_STTSTJ, col = color_palette, main = "Slope (degrees)", legend = TRUE)
  plot(slopeofslope_STTSTJ, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  plot(slope_PR_East, col = color_palette, main = "Slope (degrees)", legend = TRUE)
  
  
  slope_STTSTJ_50m = terrain(bathy_STTSTJ_agg, v = "slope", unit = 'degrees')
  slopeofslope_STTSTJ_50m = terrain(slope_STTSTJ_50m, v = "slope", unit = 'degrees')
  slope_PR_East_50m = terrain(bathy_PR_East_agg, v = "slope", unit = 'degrees')
  slopeofslope_PR_East_50m = terrain(slope_PR_East_50m, v = "slope", unit = 'degrees')
  
  plot(slope_STTSTJ_50m, col = color_palette, main = "Slope (degrees)", legend = TRUE)
  plot(slopeofslope_STTSTJ_50m, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  plot(slope_PR_East_50m, col = color_palette, main = "Slope (degrees)", legend = TRUE)
  plot(slopeofslope_PR_East_50m, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  plot(bathy_PR_East, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  plot(bathy_PR_East_agg, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  plot(bathy_STTSTJ, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  plot(bathy_STTSTJ_agg, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE) #STJ coastline is reef here ?? lord
  
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
  tpi_geo <- project(tpi, common_crs)
  tri_geo <- project(tri, common_crs)
  
  # Plot the results
  plot(slope, col = color_palette, main = "Slope (degrees)", legend = TRUE)
  plot(slopeofslope, col = color_palette, main = "Slope of slope (degrees)", legend = TRUE)
  plot(aspect, col = color_palette, main = "Aspect (degrees)", legend = TRUE)
  plot(flowdir, col = color_palette, main = "Flow Direction", legend = TRUE)
  plot(roughness, col = color_palette, main = "Roughness", legend = TRUE)
  plot(tpi, col = color_palette, main = "TPI", legend = TRUE)
  plot(tri, col = color_palette, main = "TRI", legend = TRUE)
  
  # Define the extent for the area between St. Thomas and St. Croix
  # Adjust these coordinates as needed for your specific area of interest
  extent_STT_STX <- ext(c(-65.1, -64.4, 17.6, 18.4))
  
  # plot habitat grids
  #50m NOAA sampleframe
  ggplot() +
    geom_polygon(data = df.STTSTJ_grid, aes(x = long, y = lat, group = group,),
                 fill = "green", color = "black", lwd = 1) +
    theme_bw()
  
  #my custom 650m grid
  grid_650m = readOGR(dsn = ".", layer = "STTSTJ_meso_merged_testforR")
  #reprojection
  grid_650m = spTransform(grid_650m, crs(STTSTJ_grid))
  # df.grid_650 = fortify(grid_650m)
  ggplot() +
    geom_polygon(data = df.grid_650, aes(x = long, y = lat, group = group,),
                 fill = "pink", color = "black", lwd = 1) +
    theme_bw() #+
  xlim(-65.2, -64.6) +
    ylim(18.15, 18.45)
  
  
  ### testing local and global CRS ###
  # NOTE: 8 NOV 2022: Need to repeat the below with all 2021 grids (that come from NCRMP package) in order to get all grids in the correct local CRS.
  STTSTJ_grid_2019 = readOGR(dsn = ".", layer = "STTSTJ_2019_sampleframe_localCRS-forR")
  crs(STTSTJ_grid_2019)
  crs(STTSTJ_grid)
  
  #reprojection
  STTSTJ_grid_2021_localCRS = spTransform(STTSTJ_grid, crs(STTSTJ_grid_2019))
  shapefile(x = STTSTJ_grid_2021_localCRS, file = "/Users/benja/Documents/Carib_Habitat/Carib_Habitat_QGIS/Sampleframes/STTSTJ_grid_2021_localCRS.shp")
  
  ### testing local and global CRS ###
  
  ################################## Aggregate raster data ##################################
  
  # aggregate raster data to polygon data
  #https://deepnote.com/@siew-sook-yan/R-Aggregate-raster-to-polygon-data-10a3150c-e88d-4776-bae1-4b6dcb9e3916
  #https://stackoverflow.com/questions/56110728/how-to-calculate-slope-and-aspect-ratio-values-by-giving-latitude-longitude-and
  # # mys = st_read("mys.gpkg")
  st.grid_650m = st_as_sf(grid_650m) #'st_geometry' requires 'sf' package
  crs(bathymetry_50m)
  crs(bathymetry_650m)
  crs(st.grid_650m)
  extent(bathymetry_50m)
  extent(bathymetry_650m)
  extent(st.grid_650m)
  
  slope = terrain(bathymetry_50m, opt = 'slope', unit = 'degrees')
  aspect = terrain(bathymetry_50m, opt = 'aspect', unit = 'degrees')
  flowdir = terrain(bathymetry_50m, opt = 'flowdir')
  rugosity = terrain(bathymetry_50m, opt = 'TRI') #using Terrain Ruggedness Index - should also investigate Benthic Terrain Modeler in Arc!
  
  #crop raster data to the extent of the habitat grid
  bathy_crop = crop(bathymetry_50m, extent(st.grid_650m))
  slope_crop = crop(slope, extent(st.grid_650m))
  aspect_crop = crop(aspect, extent(st.grid_650m))
  flowdir_crop = crop(flowdir, extent(st.grid_650m))
  rugosity_crop = crop(rugosity, extent(st.grid_650m))
  # df.bathy_crop = as.data.frame(bathy_crop, xy = T)
  
  #remove raster pixels outside of polygon - note I want to *keep* those pixels later when extrapolating to unknown space
  bathy_crop = mask(bathy_crop, mask = st.grid_650m)
  slope_crop = mask(slope_crop, mask = st.grid_650m)
  aspect_crop = mask(aspect_crop, mask = st.grid_650m)
  flowdir_crop = mask(flowdir_crop, mask = st.grid_650m)
  rugosity_crop = mask(rugosity_crop, mask = st.grid_650m)
  # ggplot() + #this takes 30ish seconds to run
  #   geom_tile(data = df.bathy_crop, aes(x = x, y = y, fill = final_merge_WGS1984)) +
  #   scale_fill_gradientn(colors = terrain.colors(50), limits = c(-80, 0)) +
  #   theme_bw() +
  #   xlim(-65.2, -64.6) +
  #   ylim(18.15, 18.45)
  
  plot(bathy_crop) #much quicker to run, but less control over plot
  plot(st_geometry(st.grid_650m), add = TRUE) #layers on top of previous map. may have issues when zooming in
  
  plot(slope_crop); plot(st_geometry(st.grid_650m), add = TRUE)
  plot(aspect_crop); plot(st_geometry(st.grid_650m), add = TRUE)
  plot(flowdir_crop); plot(st_geometry(st.grid_650m), add = TRUE)
  plot(rugosity_crop); plot(st_geometry(st.grid_650m), add = TRUE)
  
  #aggregate raster to spatial units. 'cellnumbers' requires 'sf' package, which speeds up processing time GREATLY
  # still may take ~1 min at 50 m resolution grid. 650 is basically instant
  bathy_cell = cellnumbers(bathy_crop, st.grid_650m)
  head(bathy_cell)
  
  slope_cell = cellnumbers(slope_crop, st.grid_650m)
  aspect_cell = cellnumbers(aspect_crop, st.grid_650m)
  flowdir_cell = cellnumbers(flowdir_crop, st.grid_650m)
  rugosity_cell = cellnumbers(rugosity_crop, st.grid_650m)
  
  #aggregate the values for all cell_ by object_
  bathy_agg = bathy_cell %>% mutate(bathy = raster::extract(bathy_crop, bathy_cell$cell_)) %>%
    group_by(object_) %>%
    summarise(bathys = mean(bathy, na.rm = TRUE)) #can choose to median or mean - not sure which is best
  nrow(bathy_agg)
  head(bathy_agg)
  
  slope_agg = slope_cell %>% mutate(slope = raster::extract(slope_crop, slope_cell$cell_)) %>%
    group_by(object_) %>%
    summarise(slopes = mean(slope, na.rm = TRUE))
  aspect_agg = aspect_cell %>% mutate(aspect = raster::extract(aspect_crop, aspect_cell$cell_)) %>%
    group_by(object_) %>%
    summarise(aspects = mean(aspect, na.rm = TRUE))
  flowdir_agg = flowdir_cell %>% mutate(flowdir = raster::extract(flowdir_crop, flowdir_cell$cell_)) %>%
    group_by(object_) %>%
    summarise(flowdirs = mean(flowdir, na.rm = TRUE))
  rugosity_agg = rugosity_cell %>% mutate(rugosity = raster::extract(rugosity_crop, rugosity_cell$cell_)) %>%
    group_by(object_) %>%
    summarise(rugositys = mean(rugosity, na.rm = TRUE))
  
  #add new column to data frame
  st.grid_650m$bathy = bathy_agg$bathys
  st.grid_650m$slope = slope_agg$slopes
  st.grid_650m$aspect = aspect_agg$aspects
  st.grid_650m$flowdir = flowdir_agg$flowdirs
  st.grid_650m$rugosity = rugosity_agg$rugositys
  
  #plot average bathymetry by grid square
  ggplot(st.grid_650m) +
    geom_sf(aes(fill = bathy)) +
    coord_sf() +
    scale_fill_gradientn(colors = terrain.colors(50), limits = c(-80, 0)) +
    theme_bw()
  ggplot(st.grid_650m) +
    geom_sf(aes(fill = slope)) +
    coord_sf() +
    scale_fill_gradientn(colors = terrain.colors(50), limits = c(0, 70)) +
    theme_bw()
  ggplot(st.grid_650m) +
    geom_sf(aes(fill = aspect)) +
    coord_sf() +
    scale_fill_gradientn(colors = terrain.colors(50), limits = c(0, 360)) +
    theme_bw()
  ggplot(st.grid_650m) +
    geom_sf(aes(fill = flowdir)) +
    coord_sf() +
    scale_fill_gradientn(colors = terrain.colors(50), limits = c(0, 120)) +
    theme_bw()
  ggplot(st.grid_650m) +
    geom_sf(aes(fill = rugosity)) +
    coord_sf() +
    scale_fill_gradientn(colors = terrain.colors(50), limits = c(0, 80)) +
    theme_bw()
  
  #https://www.rdocumentation.org/packages/raster/versions/3.0-2/topics/terrain
  #https://stackoverflow.com/questions/56110728/how-to-calculate-slope-and-aspect-ratio-values-by-giving-latitude-longitude-and
  #https://gis.stackexchange.com/questions/3310/seeking-spatial-r-tricks
  #https://medium.com/@dannymac_78271/exploring-the-relationship-between-slope-elevation-and-aspect-in-the-presidential-range-nh-fe3ef150f682
  #https://www.molecularecologist.com/2015/07/03/marmap/
  #https://www.researchgate.net/publication/255708640_Calculation_of_slope_angle_from_bathymetry_data_using_GIS_-_effects_of_computation_algorithms_data_resolution_and_analysis_scale
  #https://link.springer.com/referenceworkentry/10.1007/978-90-481-2639-2_141#:~:text=Rugosity%20is%20an%20estimate%20of,textural%20characteristics%20of%20a%20surface.
  #https://dusk.geo.orst.edu/djl/samoa/BTM_Exercise.pdf
  
  #test to see how R would perform doing the above with 50 m res (ongoing)
  st.grid_50m = st_as_sf(STTSTJ_grid) #'st_geometry' requires 'sf' package
  st.grid_650m = st.grid_50m
  #then run the above again
  
  ################################## Test response variable data ##################################
  
  STTSTJ_field_species$LAT_DEGREES = as.numeric(STTSTJ_field_species$LAT_DEGREES)
  STTSTJ_field_species$LON_DEGREES = as.numeric(STTSTJ_field_species$LON_DEGREES)
  
  #calling 'extract_grid_data'
  # 'grid' will be a regional NCRMP sampleframe, or custom one
  # FUNCTION: extract_grid_data <- function(field_sites, grid, region)
  STTSTJ_test2 = extract_grid_data(field_sites = STTSTJ_field_sites, grid = grid_650m, region = "STTSTJ")
  # STTSTJ_test = extract_grid_data(field_sites = STTSTJ_field_sites, grid = STTSTJ_grid, region = "STTSTJ")
  STTSTJ_test2$unique_ID = as.numeric(STTSTJ_test2$unique_ID)
  st.grid_650m$unique_ID = as.numeric(st.grid_650m$unique_ID)
  STTSTJ_test2 = dplyr::filter(STTSTJ_test2, cover_group %in% c("HARD CORALS"))
  
  #merge new bathymetry values to every row in 'STTSTJ_test2' by the Unique_ID
  coral = merge(STTSTJ_test2, st.grid_650m, by = "unique_ID", all.x = TRUE)
  saveRDS(coral, file = "coral.rds") #save clean data to an object file for BUGS
  coral = readRDS(file = "coral.rds") #used to restore the object file
  
  # # manual way to do 'extract_grid_data' (pulled from that script)
  # # subset coordinates from site df
  # xy  = STTSTJ_field_species[, c("LAT_DEGREES", "LON_DEGREES")]
  # # specify coordinates
  # sp::coordinates(xy) = ~LON_DEGREES+LAT_DEGREES
  #
  # # project coordinates to WGS84
  # sp::proj4string(xy) = sp::CRS("+init=epsg:4326")
  #
  # # overlay PR points on the reprojected grid
  # overlay_points = sp::over(xy, grid_650m, fn = NULL)
  #
  # # add rownames to the overlay df and the original df
  # overlay_points$rownames = rownames(overlay_points)
  # overlay_points$REGION = "STTSTJ"
  #
  # ## For ALL sampling geographies:
  # # add rownames to the original df for subsequent combining
  # STTSTJ_field_sites$rownames = rownames(STTSTJ_field_sites)
  #
  # # merge input sample file with grid overly df by rownames (same order going in)
  # samplesites_grid = dplyr::left_join(STTSTJ_field_sites, overlay_points,
  #                                      by = c("REGION", "rownames"))
  # samplesites_grid = dplyr::filter(samplesites_grid, cover_group %in% c("HARD CORALS"))
  
  ### playing with the response variable data ###
  
  ### OBJECTIVES ###
  # 1.) Retrieve NCRMP grids for PR and STT/STJ/STX (all regions)
  #      - Assess overlap (or lack thereof) between repeated measurements
  #      - How to move grid file over to GIS?
  # 2.) Retrieve benthic cover for all regions
  #      - How does this data "talk" to the grid data?
  #      - How to move cover raster (or vector) over to GIS?
  #      - Pool the grid at different resolutions for testing (e.g., 350 m; 650 m). Might require GIS
  #      - See how many "repeated measurements" occur after broadening the grid. May *want* these repeats
  
  ################################## Save objects/workspace ##################################
  
  # #save terra objects and then workspace for use in downstream scripts
  # save_spat_objects() #call from functions.R
  # save.image(file = here("output", "calculate_bathy_rasters_workspace.RData"))
  