
# NOTE - consider using glue package to insert dates for output files
#   https://stackoverflow.com/questions/73584505/how-to-write-a-file-into-a-specific-folder-using-here-package

# potentially useful:
#https://www.benjaminbell.co.uk/2019/08/bathymetric-maps-in-r-colour-palettes.html
#https://stackoverflow.com/questions/20581746/increase-resolution-of-color-scale-for-values-close-to-zero

# .rs.restartR(clean = TRUE)

library(here)
library(terra) #this requires (at least for me on M1 Macbook in August 2024) homebrew installation of proj and gdal (https://github.com/OSGeo/gdal/pull/7389). could cause issues with macports installation of other things for QGIS but we'll see
source(here("src/functions.R"))

# Determine a common CRS for the entire dataset
# common_crs <- st_crs(4326)  # WGS84. suitable for global use. geographic CRS
common_crs <- "EPSG:4269" # NAD83. like WGS84 but best for contiguous United States. geographic CRS
# common_crs <- st_crs(32620)  # WGS 84 / UTM zone 20N (suitable for Puerto Rico and the Virgin Islands). projected CRS

# Specify the path to ESRI geodatabase provided by Jeremiah Blondeau
gdb_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Bathymetry/NOAA_LIDAR_Blondeau/US_Caribbean_Bathy_Mocaics.gdb"

# List all raster layers within the geodatabase
# Check sub-datasets
describe(gdb_path, sds = TRUE)
# raster_layers <- rast(gdb_path)
# names(raster_layers)  # Print the names of the raster layers
# Load specific sub-datasets from the geodatabase
bathy_STTSTJ = rast(gdb_path, subds = "STTSTJ_2m")
bathy_STX <- rast(gdb_path, subds = "STX_2m")
bathy_PR_East <- rast(gdb_path, subds = "PuertoRico_East_2m")
bathy_PR_South <- rast(gdb_path, subds = "PuertoRico_South_2m")
bathy_PR_West <- rast(gdb_path, subds = "PuertoRico_West_2m")
bathy_PR_North <- rast(gdb_path, subds = "PuertoRico_North_2m")

#downscale the resolution of the bathymetry grid
# NOTE - may eventually simply snap the bathymetry to the NCRMP grid itself
#
# Define the aggregation factor
agg_factor <- 25

# Aggregate mean of the bathymetry to 50 m resolution (this may change; was doing this for ease of use early on)
bathy_STTSTJ_agg <- aggregate(bathy_STTSTJ, fact = agg_factor, fun = mean, na.rm = TRUE)
bathy_STX_agg <- aggregate(bathy_STX, fact = agg_factor, fun = mean, na.rm = TRUE)
bathy_PR_East_agg <- aggregate(bathy_PR_East, fact = agg_factor, fun = mean, na.rm = TRUE)
bathy_PR_South_agg <- aggregate(bathy_PR_South, fact = agg_factor, fun = mean, na.rm = TRUE)
bathy_PR_West_agg <- aggregate(bathy_PR_West, fact = agg_factor, fun = mean, na.rm = TRUE)
bathy_PR_North_agg <- aggregate(bathy_PR_North, fact = agg_factor, fun = mean, na.rm = TRUE)

# Save each input raster separately
# writeRaster(bathy_STTSTJ_agg, here("output", "bathy_STTSTJ_agg.tif"), overwrite = TRUE)
# writeRaster(bathy_STX_agg, here("output", "bathy_STX_agg.tif"), overwrite = TRUE)
# writeRaster(bathy_PR_East_agg, here("output", "bathy_PR_East_agg.tif"), overwrite = TRUE)
# writeRaster(bathy_PR_South_agg, here("output", "bathy_PR_South_agg.tif"), overwrite = TRUE)
# writeRaster(bathy_PR_West_agg, here("output", "bathy_PR_West_agg.tif"), overwrite = TRUE)
# writeRaster(bathy_PR_North_agg, here("output", "bathy_PR_North_agg.tif"), overwrite = TRUE)
# saveRDS(bathy_STTSTJ_agg, here("output", "bathy_STTSTJ_agg.rds"))
# saveRDS(bathy_STX_agg, here("output", "bathy_STX_agg.rds"))
# saveRDS(bathy_PR_East_agg, here("output", "bathy_PR_East_agg.rds"))
# saveRDS(bathy_PR_South_agg, here("output", "bathy_PR_South_agg.rds"))
# saveRDS(bathy_PR_West_agg, here("output", "bathy_PR_West_agg.rds"))
# saveRDS(bathy_PR_North_agg, here("output", "bathy_PR_North_agg.rds"))

# Call the function
save_spat_objects()

# # Merge the aggregated rasters. merge is quicker but maybe mosaic is better?
# agg_rasters <- c(bathy_STTSTJ_agg, bathy_STX_agg, bathy_PR_East_agg, bathy_PR_South_agg, bathy_PR_West_agg, bathy_PR_North_agg)
# bathy_merged <- mosaic(agg_rasters, fun = mean, na.rm = TRUE)
bathy_merged = merge(bathy_STTSTJ_agg, bathy_STX_agg, bathy_PR_East_agg, bathy_PR_South_agg, bathy_PR_West_agg, bathy_PR_North_agg)
# writeRaster(bathy_merged, here("output", "bathy_50m.tif"), overwrite = TRUE)
saveRDS(bathy_merged, here("output", "bathy_50m.rds"))

# # Save all relevant objects to an .RData file in the output directory
# save(
#   cover_site_STTSTJ, cover_spp_STTSTJ,
#   cover_site_PR, cover_spp_PR,
#   cover_site_STX, cover_spp_STX,
#   sampleframe_USVI_transformed, sampleframe_PR_transformed, sampleframe_DCRMP_transformed,
#   bathy_STTSTJ_agg, bathy_STX_agg,
#   bathy_PR_East_agg, bathy_PR_South_agg, bathy_PR_West_agg, bathy_PR_North_agg,
#   file = here("output", "aggregated_data.RData")
# )

# Save all objects in the current R session to an .RData file in the output directory
save.image(file = here("output", "workspace.RData"))
