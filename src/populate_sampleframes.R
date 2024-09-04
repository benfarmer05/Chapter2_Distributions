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

