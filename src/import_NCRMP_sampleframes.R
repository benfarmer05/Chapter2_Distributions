
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

cover_site_STTSTJ = NCRMP_STTSTJ_2013_21_percent_cover_site
cover_spp_STTSTJ = NCRMP_STTSTJ_2013_21_percent_cover_species
cover_site_PR = NCRMP_PRICO_2014_21_percent_cover_site
cover_spp_PR = NCRMP_PRICO_2014_21_percent_cover_species
cover_site_STX = NCRMP_STX_2015_21_percent_cover_site
cover_spp_STX = NCRMP_STX_2015_21_percent_cover_species

# Optionally, remove the original variable names if they are no longer needed
rm(NCRMP_STTSTJ_2013_21_percent_cover_site, NCRMP_STTSTJ_2013_21_percent_cover_species, NCRMP_PRICO_2014_21_percent_cover_site,
   NCRMP_PRICO_2014_21_percent_cover_species, NCRMP_STX_2015_21_percent_cover_site, NCRMP_STX_2015_21_percent_cover_species)

# Refactor YEAR column for each dataset
cover_site_STTSTJ$YEAR <- as.factor(cover_site_STTSTJ$YEAR)
cover_spp_STTSTJ$YEAR <- as.factor(cover_spp_STTSTJ$YEAR)

cover_site_PR$YEAR <- as.factor(cover_site_PR$YEAR)
cover_spp_PR$YEAR <- as.factor(cover_spp_PR$YEAR)

cover_site_STX$YEAR <- as.factor(cover_site_STX$YEAR)
cover_spp_STX$YEAR <- as.factor(cover_spp_STX$YEAR)

# Optionally, check levels of YEAR column in each dataset
levels(cover_site_STTSTJ$YEAR)
levels(cover_spp_STTSTJ$YEAR)
levels(cover_site_PR$YEAR)
levels(cover_spp_PR$YEAR)
levels(cover_site_STX$YEAR)
levels(cover_spp_STX$YEAR)

#load sampleframes
sampleframe_USVI = st_read(here("data/NCEI_inport/USVI/2023/1.1/data/0-data/Data_Sets/Sample_Frames/NCRMP_STTSTJ_2023_SampleFrame.shp"))
sampleframe_PR = st_read(here("data/NCEI_inport/PR/2023/1.1/data/0-data/Data_Sets/Sample_Frames/NCRMP_PR2023_SAMPLE_FRAME.shp"))
sampleframe_DCRMP = st_read(here("data/MesophoticSampGrid/MesophoticSampGrid.shp"))

# Determine a common CRS for the entire dataset; consider using a geographic CRS like WGS84
common_crs <- st_crs(5070) # NAD83 / Conus Albers or use EPSG:4326 for WGS84 if desired

# Transform all datasets to the common CRS
sampleframe_USVI_transformed <- st_transform(sampleframe_USVI, common_crs)
sampleframe_PR_transformed <- st_transform(sampleframe_PR, common_crs)
sampleframe_DCRMP_transformed <- st_transform(sampleframe_DCRMP, common_crs)

# STOPPING POINT - 30 AUGUST 2024
#   - I need to consider what I want to accomplish. Merging these sampleframes probably isn't going to cut it, from what I am seeing.
#     but do I even need to merge them? depends how I want to snap the bathy data to the frames, if at all. could just extract what I need
#     using the frames as a reference and deal with everything else downstream
#
#   - yes okay the more I think about it, why not just use the native NCRMP grid(s) directly to calculate bathymetry, slope, aspect etc.
#       within their grid squares? this might be really computationally expensive, I don't know. but I'll have to figure out something
#   - and anyways, next steps will require
#       - test how bad the mosaicing / artifact patterns matter in Blondeau's NOAA product (produce slope/complexity and plot)
#       - compare with NOAA CRM exports. bring in mine and Holstein's to R, clip to 0-50 m, merge them (may need to clip PR against USVI
#         first). then produce slope/complexity and plot. compare plot with that of Blondeau's product
#       - lastly, and actually I'll do this first, is porting data from Allen Coral Atlas and seeing how that bathy looks! and turbidity

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






############# 29 August 2024

# file_temp = rast(dsn = '/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Bathymetry/NOAA_LIDAR_Blondeau/US_Caribbean_Bathy_Mocaics.gdb', subds = 'PuertoRico_West_2m')

# Specify the path to your ESRI geodatabase
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
# NOTE - may eventually simply snap the bathymetry to the NCRMP grid itself. probably makes more sense
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

# # Merge the aggregated rasters. merge would be quicker but maybe mosaic is better?
# agg_rasters <- c(bathy_STTSTJ_agg, bathy_STX_agg, bathy_PR_East_agg, bathy_PR_South_agg, bathy_PR_West_agg, bathy_PR_North_agg)
# merged_bathy <- mosaic(agg_rasters, fun = mean, na.rm = TRUE)
merged_bathy = merge(bathy_STTSTJ_agg, bathy_STX_agg, bathy_PR_East_agg, bathy_PR_South_agg, bathy_PR_West_agg, bathy_PR_North_agg)

# Save the merged raster
output_dir <- here("output")
output_file <- file.path(output_dir, "bathy_50m.tif")
writeRaster(merged_bathy, filename = output_file, overwrite = TRUE)

# PLOTTING
#
#plot method, with terra natively
raster_data <- rast(output_file)

# Define the color palette for the plot
color_palette <- rev(terrain.colors(100))

# Plot the raster
plot(raster_data,
     col = color_palette,
     zlim = c(-50, 0),
     main = "Bathymetry (50m Resolution)",
     legend = TRUE)

# #ggplot method, with tidyterra
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
tm_shape(merged_bathy) +
  tm_raster(style = 'cont')
bathy_map = tm_shape(merged_bathy) +
  tm_raster(style = 'cont')

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
raster_to_matrix(merged_bathy) |> height_shade() |> plot_map()

# # in 3D
# elmat %>%
#   # sphere_shade(texture = "desert") %>%
#   # add_water(detect_water(elmat), color = "desert") %>%
#   # add_shadow(ray_shade(elmat, zscale = 3), 0.5) %>%
#   add_shadow(ambient_shade(elmat), 0) %>%
#   plot_3d(elmat, zscale = 10, fov = 0, theta = 135, zoom = 0.75, phi = 45, windowsize = c(1000, 800))
# Sys.sleep(0.2)
# render_snapshot()

### Port in Allen Coral Atlas bathymetry & turbidity data. maybe try this using CLI in the future - https://openoceans.xyz/projects/pycoral/
# okay upon looking at the bathy - it's very grainy. I guess that makes sense for satellite imaging, but either way cannot use.
# note that it is in centimeters, and >depth = positive values.
#   - the turbidity raster also doesn't seem like it has great resolution either. very focused on super nearshore reefs. could be useful
#     for PR though ? not sure
bathy_ACA = rast(here("data/Allen_Coral_Atlas/PR_USVI-20240830214543/Bathymetry---composite-depth/bathymetry_0.tif"))
turbidity_ACA = rast(here("data/Allen_Coral_Atlas/PR_USVI-20240830214543/Turbidity-2019/turbidity-annual_0.tif"))
bathy_ACA <- clamp(bathy_ACA, lower = 0, upper = 50, values = TRUE)
tm_shape(bathy_ACA) +
  tm_raster(style = 'cont')
tm_shape(turbidity_ACA) +
  tm_raster(style = 'cont')

bathy_ACA_Culebra = rast(here("data/Allen_Coral_Atlas/Culebrita-20240830215241/Bathymetry---composite-depth/bathymetry_0.tif"))
turbidity_ACA_Culebra = rast(here("data/Allen_Coral_Atlas/Culebrita-20240830215241/Turbidity-2019/turbidity-annual_0.tif"))
bathy_ACA_Culebra <- clamp(bathy_ACA_Culebra, lower = 0, upper = 50, values = TRUE)
tm_shape(bathy_ACA_Culebra) +
  tm_raster(style = 'cont')
turbidity_ACA_Culebra <- clamp(turbidity_ACA_Culebra, lower = 0, upper = 100, values = TRUE)
tm_shape(turbidity_ACA_Culebra) +
  tm_raster(style = 'cont')


########## STOPPING POINT - 30 August 2024
# below was some old scripting from 2022 that I am adopting for the here and now

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


### playing with local and global CRS ###
  # NOTE: 8 NOV 2022: Need to repeat the below with all 2021 grids (that come from NCRMP package) in order to get all grids in the correct local CRS.
STTSTJ_grid_2019 = readOGR(dsn = ".", layer = "STTSTJ_2019_sampleframe_localCRS-forR")
crs(STTSTJ_grid_2019)
crs(STTSTJ_grid)

#reprojection
STTSTJ_grid_2021_localCRS = spTransform(STTSTJ_grid, crs(STTSTJ_grid_2019))
shapefile(x = STTSTJ_grid_2021_localCRS, file = "/Users/benja/Documents/Carib_Habitat/Carib_Habitat_QGIS/Sampleframes/STTSTJ_grid_2021_localCRS.shp")

### playing with local and global CRS ###

### AGGREGATING RASTER DATA ###
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

#a little test to see how R would perform doing the above with 50 m res (ongoing)
st.grid_50m = st_as_sf(STTSTJ_grid) #'st_geometry' requires 'sf' package
st.grid_650m = st.grid_50m
#then run the above again

### AGGREGATING RASTER DATA ###


### playing with the response variable data ###
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


