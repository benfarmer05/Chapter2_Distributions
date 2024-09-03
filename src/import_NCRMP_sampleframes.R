
# NOTE - consider using glue package to insert dates for output files
#   https://stackoverflow.com/questions/73584505/how-to-write-a-file-into-a-specific-folder-using-here-package

# potentially useful:
#https://www.benjaminbell.co.uk/2019/08/bathymetric-maps-in-r-colour-palettes.html
#https://stackoverflow.com/questions/20581746/increase-resolution-of-color-scale-for-values-close-to-zero

# .rs.restartR(clean = TRUE)

library(sf)
library(here)
library(terra) #this requires (at least for me on M1 Macbook in August 2024) homebrew installation of proj and gdal (https://github.com/OSGeo/gdal/pull/7389). could cause issues with macports installation of other things for QGIS but we'll see
library(here)
library(tidyterra)
library(ggplot2)
library(tmap)
library(rayshader) #this requires installation of XQuartz on MacOS, and possibly OpenGL if it isn't installed

# # these were old file paths direct from either github or NCEI - maybe can delete. I think I was trying to find the source data
# #  (it's here somewhere though)
# cover_USVI_2013 = read.csv(here("data/NCRMP_USVI_2013_2021", "NCRMP_USVI2013_Benthic_Data01_BenthicCover.csv"))
# cover_USVI_2021 = read.csv(here("data/NCRMP_USVI_2013_2021", "NCRMP_USVI2021_Benthic_Data01_BenthicCover.csv"))

#retrieve data from cloned NCRMP benthics GitHub repo
load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_STTSTJ_2013_21_percent_cover_site.rda"))
load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_STTSTJ_2013_21_percent_cover_species.rda"))
load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_PRICO_2014_21_percent_cover_site.rda"))
load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_PRICO_2014_21_percent_cover_species.rda"))
load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_STX_2015_21_percent_cover_site.rda"))
load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_STX_2015_21_percent_cover_species.rda"))

# cover_site_STTSTJ = NCRMP_STTSTJ_2013_21_percent_cover_site
# cover_spp_STTSTJ = NCRMP_STTSTJ_2013_21_percent_cover_species
# cover_site_PR = NCRMP_PRICO_2014_21_percent_cover_site
# cover_spp_PR = NCRMP_PRICO_2014_21_percent_cover_species
# cover_site_STX = NCRMP_STX_2015_21_percent_cover_site
# cover_spp_STX = NCRMP_STX_2015_21_percent_cover_species
# 
# # Optionally, remove the original variable names if they are no longer needed
# rm(NCRMP_STTSTJ_2013_21_percent_cover_site, NCRMP_STTSTJ_2013_21_percent_cover_species, NCRMP_PRICO_2014_21_percent_cover_site,
#    NCRMP_PRICO_2014_21_percent_cover_species, NCRMP_STX_2015_21_percent_cover_site, NCRMP_STX_2015_21_percent_cover_species)
# 
# # Refactor YEAR column for each dataset
# cover_site_STTSTJ$YEAR <- as.factor(cover_site_STTSTJ$YEAR)
# cover_spp_STTSTJ$YEAR <- as.factor(cover_spp_STTSTJ$YEAR)
# 
# cover_site_PR$YEAR <- as.factor(cover_site_PR$YEAR)
# cover_spp_PR$YEAR <- as.factor(cover_spp_PR$YEAR)
# 
# cover_site_STX$YEAR <- as.factor(cover_site_STX$YEAR)
# cover_spp_STX$YEAR <- as.factor(cover_spp_STX$YEAR)
# 
# # Optionally, check levels of YEAR column in each dataset
# levels(cover_site_STTSTJ$YEAR)
# levels(cover_spp_STTSTJ$YEAR)
# levels(cover_site_PR$YEAR)
# levels(cover_spp_PR$YEAR)
# levels(cover_site_STX$YEAR)
# levels(cover_spp_STX$YEAR)

#load sampleframes. consider using sf/stars in the future for large datasets, especially data cubes, netCDF, etc. (https://github.com/r-spatial/stars/issues/633#issuecomment-1597356513)
sampleframe_USVI = vect(here("data/NCEI_inport/USVI/2023/1.1/data/0-data/Data_Sets/Sample_Frames/NCRMP_STTSTJ_2023_SampleFrame.shp"))
sampleframe_PR = vect(here("data/NCEI_inport/PR/2023/1.1/data/0-data/Data_Sets/Sample_Frames/NCRMP_PR2023_SAMPLE_FRAME.shp"))
sampleframe_DCRMP = vect(here("data/MesophoticSampGrid/MesophoticSampGrid.shp"))

# Determine a common CRS for the entire dataset
# common_crs <- st_crs(4326)  # WGS84. suitable for global use. geographic CRS
common_crs <- "EPSG:4269" # NAD83. like WGS84 but best for contiguous United States. geographic CRS
# common_crs <- st_crs(32620)  # WGS 84 / UTM zone 20N (suitable for Puerto Rico and the Virgin Islands). projected CRS

# Transform all datasets to the common CRS
sampleframe_USVI_transformed <- project(sampleframe_USVI, common_crs) #st_transform(sampleframe_USVI, common_crs)
sampleframe_PR_transformed <- project(sampleframe_PR, common_crs)
sampleframe_DCRMP_transformed <- project(sampleframe_DCRMP, common_crs)

# NOTE
#   - I need to consider what I want to accomplish. Merging these sampleframes (below is attempt) probably isn't going to cut it, from what
#     I am seeing. but do I even need to merge them? depends how I want to snap the bathy data to the frames, if at all. could just extract
#     what I need using the frames as a reference and deal with everything else downstream

# # Convert 'sf' objects to 'terra' SpatVector objects for better handling
# vect_USVI <- vect(sampleframe_USVI_transformed)
# vect_PR <- vect(sampleframe_PR_transformed)
# vect_DCRMP <- vect(sampleframe_DCRMP_transformed)
# 
# # Chunk-wise merge approach
# # Split the large PR dataset into chunks
# chunk_size <- 10000  # Define chunk size for processing
# num_chunks <- ceiling(nrow(vect_PR) / chunk_size)
# 
# # Initialize the progress bar
# pb <- txtProgressBar(min = 0, max = num_chunks, style = 3)
# 
# # Process in chunks
# for (i in 1:num_chunks) {
#   start_row <- ((i - 1) * chunk_size) + 1
#   end_row <- min(i * chunk_size, nrow(vect_PR))
# 
#   # Subset PR chunk
#   chunk_PR <- vect_PR[start_row:end_row, ]
# 
#   # Merge chunk with USVI and DCRMP using terra's union
#   merged_chunk <- union(vect_USVI, chunk_PR)
#   merged_chunk <- union(merged_chunk, vect_DCRMP)
# 
#   # Append merged chunk to the list
#   merged_chunks[[i]] <- merged_chunk
# 
#   # Update progress bar
#   setTxtProgressBar(pb, i)
# }
# 
# # Combine all chunks into a final merged SpatVector
# final_merged <- do.call(union, merged_chunks)
# 
# # Convert back to 'sf' object if needed
# final_merged_sf <- st_as_sf(final_merged)
# 
# # Save the merged output to a new shapefile
# st_write(final_merged_sf, here("output", "merged_sampleframes.shp"))

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

# Aggregate the rasters
bathy_STTSTJ_agg <- aggregate(bathy_STTSTJ, fact = agg_factor, fun = mean, na.rm = TRUE)
bathy_STX_agg <- aggregate(bathy_STX, fact = agg_factor, fun = mean, na.rm = TRUE)
bathy_PR_East_agg <- aggregate(bathy_PR_East, fact = agg_factor, fun = mean, na.rm = TRUE)
bathy_PR_South_agg <- aggregate(bathy_PR_South, fact = agg_factor, fun = mean, na.rm = TRUE)
bathy_PR_West_agg <- aggregate(bathy_PR_West, fact = agg_factor, fun = mean, na.rm = TRUE)
bathy_PR_North_agg <- aggregate(bathy_PR_North, fact = agg_factor, fun = mean, na.rm = TRUE)

# Save each input raster separately
writeRaster(bathy_STTSTJ_agg, here("output", "bathy_STTSTJ_agg.tif"), overwrite = TRUE)
writeRaster(bathy_STX_agg, here("output", "bathy_STX_agg.tif"), overwrite = TRUE)
writeRaster(bathy_PR_East_agg, here("output", "bathy_PR_East_agg.tif"), overwrite = TRUE)
writeRaster(bathy_PR_South_agg, here("output", "bathy_PR_South_agg.tif"), overwrite = TRUE)
writeRaster(bathy_PR_West_agg, here("output", "bathy_PR_West_agg.tif"), overwrite = TRUE)
writeRaster(bathy_PR_North_agg, here("output", "bathy_PR_North_agg.tif"), overwrite = TRUE)

# # Merge the aggregated rasters. merge is quicker but maybe mosaic is better?
# agg_rasters <- c(bathy_STTSTJ_agg, bathy_STX_agg, bathy_PR_East_agg, bathy_PR_South_agg, bathy_PR_West_agg, bathy_PR_North_agg)
# bathy_merged <- mosaic(agg_rasters, fun = mean, na.rm = TRUE)
bathy_merged = merge(bathy_STTSTJ_agg, bathy_STX_agg, bathy_PR_East_agg, bathy_PR_South_agg, bathy_PR_West_agg, bathy_PR_North_agg)
writeRaster(bathy_merged, here("output", "bathy_50m.tif"), overwrite = TRUE)

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
