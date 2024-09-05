

# .rs.restartR(clean = TRUE)

# library(sf)
library(here)
library(terra)
# library(tidyterra)
# library(ggplot2)
# library(tmap)
# library(rayshader) #this requires installation of XQuartz on MacOS, and possibly OpenGL if it isn't installed
# library(scico)
# library(RColorBrewer)

#import objects from from import_merge_rasters.R
load(here("output", "workspace.RData"))
load_spat_objects(directory = here("output")) #extra care is required for terra-produced objects with pointers

# # After loading the workspace, re-write the SpatRasters from .tif (required because of the way terra works with R objects, I think)
# # bathy_merged = rast(here("output", "bathy_50m.tif"))
# # bathy_STTSTJ_agg <- rast(here("output", "bathy_STTSTJ_agg.tif"))
# # bathy_STX_agg <- rast(here("output", "bathy_STX_agg.tif"))
# # bathy_PR_East_agg <- rast(here("output", "bathy_PR_East_agg.tif"))
# # bathy_PR_South_agg <- rast(here("output", "bathy_PR_South_agg.tif"))
# # bathy_PR_West_agg <- rast(here("output", "bathy_PR_West_agg.tif"))
# # bathy_PR_North_agg <- rast(here("output", "bathy_PR_North_agg.tif"))
# bathy_merged_2m = readRDS(here("output", "bathy_merged_2m.rds"))
# bathy_merged_50m = readRDS(here("output", "bathy_merged_50m.rds"))
# bathy_STTSTJ_agg <- readRDS(here("output", "bathy_STTSTJ_agg.rds"))
# bathy_STX_agg <- readRDS(here("output", "bathy_STX_agg.rds"))
# bathy_PR_East_agg <- readRDS(here("output", "bathy_PR_East_agg.rds"))
# bathy_PR_South_agg <- readRDS(here("output", "bathy_PR_South_agg.rds"))
# bathy_PR_West_agg <- readRDS(here("output", "bathy_PR_West_agg.rds"))
# bathy_PR_North_agg <- readRDS(here("output", "bathy_PR_North_agg.rds"))

# # these were old file paths direct from either github or NCEI - maybe can delete. I think I was trying to find the source data
# #  (it's here somewhere though)
# cover_USVI_2013 = read.csv(here("data/NCRMP_USVI_2013_2021", "NCRMP_USVI2013_Benthic_Data01_BenthicCover.csv"))
# cover_USVI_2021 = read.csv(here("data/NCRMP_USVI_2013_2021", "NCRMP_USVI2021_Benthic_Data01_BenthicCover.csv"))

# #retrieve data from cloned NCRMP benthics GitHub repo
# load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_STTSTJ_2013_21_percent_cover_site.rda"))
# load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_STTSTJ_2013_21_percent_cover_species.rda"))
# load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_PRICO_2014_21_percent_cover_site.rda"))
# load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_PRICO_2014_21_percent_cover_species.rda"))
# load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_STX_2015_21_percent_cover_site.rda"))
# load(here("data/NCRMP_benthics-master/ncrmp.benthics.analysis/data", "NCRMP_STX_2015_21_percent_cover_species.rda"))

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

# # Transform all datasets to the common CRS
# sampleframe_USVI_transformed <- project(sampleframe_USVI, common_crs) #st_transform(sampleframe_USVI, common_crs)
# sampleframe_PR_transformed <- project(sampleframe_PR, common_crs)
# sampleframe_DCRMP_transformed <- project(sampleframe_DCRMP, common_crs)


# TESTING # STOPPING POINT - 4 Sep 2024
# Example of adding a new ID field
sampleframe_USVI$ID <- 1:nrow(sampleframe_USVI)  # Assign unique IDs, for example
sample_frame_raster_USVI <- rasterize(sampleframe_USVI, bathy_merged_50m, field = "ID")

sampleframe_PR = project(sampleframe_PR, crs(sampleframe_USVI)) #reproject PR sample frame to that of USVI (NAD83 / UTM zone 20N)
sampleframe_PR$ID <- 1:nrow(sampleframe_PR)  # Assign unique IDs, for example
sample_frame_raster_PR <- rasterize(sampleframe_PR, bathy_merged, field = "ID")

# Compute zonal statistics
bathy_stats_USVI <- zonal(bathy_merged, sampleframe_USVI, fun = c("mean", "sd"), na.rm = TRUE)

# Print results
print(bathy_stats)




# STOPPING POINT - 4 Sept 2024

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

