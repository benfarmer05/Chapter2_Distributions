  
  # .rs.restartR(clean = TRUE)
  
  library(here)
  library(terra)
  
  #import objects from from import_merge_rasters.R
  load(here("output/import_merge_rasters_workspace.RData"))
  
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
  
  
  
  # # TESTING # STOPPING POINT - 4 SEP 2024
  # # Example of adding a new ID field
  # sampleframe_USVI$ID <- 1:nrow(sampleframe_USVI)  # Assign unique IDs, for example
  # frame_USVI_extent <- ext(sampleframe_USVI)
  # bathy_cropped_frame_USVI <- crop(bathy_merged_50m, frame_USVI_extent)
  # # bathy_hires_cropped_frame_USVI = crop(bathy_STTSTJ, frame_USVI_extent)
  # plot(bathy_cropped_frame_USVI)
  # # plot(bathy_hires_cropped_frame_USVI)
  # 
  # sample_frame_raster_USVI <- rasterize(sampleframe_USVI, bathy_cropped_frame_USVI, field = "ID")
  # plot(sampleframe_USVI)
  # plot(sample_frame_raster_USVI)
  # 
  # sampleframe_PR = project(sampleframe_PR, crs(sampleframe_USVI)) #reproject PR sample frame to that of USVI (NAD83 / UTM zone 20N)
  # sampleframe_PR$ID <- 1:nrow(sampleframe_PR)  # Assign unique IDs, for example
  # sample_frame_raster_PR <- rasterize(sampleframe_PR, bathy_merged_50m, field = "ID")
  # 
  # # Compute zonal statistics
  # bathy_stats_USVI <- zonal(bathy_merged_50m, sampleframe_USVI, fun = c("mean", "sd"), na.rm = TRUE)
  # 
  # # Print results
  # print(bathy_stats_USVI)
  
  
  # First, make sure your sample frame is properly prepared with unique IDs
  sampleframe_USVI$ID <- 1:nrow(sampleframe_USVI)
  
  # Create a raster representation of your sample frame
  sample_frame_raster_USVI <- rasterize(sampleframe_USVI, bathy_merged_50m, field = "ID")
  
  # Calculate mean bathymetry within each sample frame grid cell
  mean_bathy_by_cell <- zonal(bathy_merged_50m, sample_frame_raster_USVI, fun = "mean", na.rm = TRUE)
  
  # Check the structure of the zonal statistics result
  print(names(mean_bathy_by_cell))
  
  # Let's assume the correct column names are determined from the above step
  # The first column typically contains zone IDs, and the second contains the calculated statistic
  
  # Rename columns for clarity if needed
  names(mean_bathy_by_cell)[1] <- "zone_id"  # First column with zone identifiers
  names(mean_bathy_by_cell)[2] <- "mean_depth"  # Second column with mean values
  
  # Merge with correct column names
  sampleframe_with_bathy <- merge(sampleframe_USVI, mean_bathy_by_cell, 
                                  by.x = "ID", by.y = "zone_id")
  
  #set extent
  st_thomas_extent <- ext(280000, 295000, 2025000, 2040000)
  
  coastline = as.contour(bathy_merged_50m, level = -10)
  
  # Now plot the results
  plot(sampleframe_with_bathy, "mean_depth", 
       main = "Mean Bathymetry by Sample Frame Cell", 
       col = terrain.colors(100),
       ext = st_thomas_extent)
  
  plot(sampleframe_USVI,  
       main = "Mean Bathymetry by Sample Frame Cell", 
       ext = st_thomas_extent)
  
  
  # # Add coastlines on top
  # lines(coastline, col="black", lwd=1.5)
  
  # #interesting overlay of bathymetry and sampleframe
  # bathy_st_thomas <- crop(bathy_merged_50m, st_thomas_extent)
  # 
  # plot(bathy_st_thomas, 
  #      main = "St. Thomas - Bathymetry with Sample Frame Overlay",
  #      # col = blues9, # Or another color palette good for bathymetry
  #      ext = st_thomas_extent)
  # plot(sampleframe_with_bathy, "mean_depth", 
  #      col = terrain.colors(100, alpha=0.7), # Adding transparency
  #      add = TRUE)
  # legend("topright", title="Mean Depth (m)", 
  #        legend=round(seq(min(sampleframe_with_bathy$mean_depth, na.rm=TRUE), 
  #                         max(sampleframe_with_bathy$mean_depth, na.rm=TRUE), length.out=5), 1),
  #        fill=terrain.colors(5), bty="n")
  
  
  
  
  
  
  
  
  
  
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
  
