# Clear workspace
rm(list=ls())

# Load necessary libraries
# library(raster)
# library(ggplot2)
library(sf)
# library(tabularaster)
# library(dplyr)
library(here)


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


library(terra)

# Specify the path to your ESRI geodatabase
gdb_path <- "/Users/benja/Documents/Farmer_Ben_Dissertation/QGIS_Dissertation/data/Bathymetry/NOAA_LIDAR_Blondeau/US_Caribbean_Bathy_Mocaics.gdb"

# #the below did not work because GDAL nested in terra is not capable of this yet (I guess)
# # List all raster layers within the geodatabase
# raster_layers <- rast(gdb_path)
# names(raster_layers)  # Print the names of the raster layers
# 
# library(rgdal)
# subset(ogrDrivers(), grepl("GDB", name))
# ogrListLayers("/path/to/folder.gdb")

#with the above not working, I wanted to try and resample the rasters directly in QGIS, export as geoTIFFS, and import here
#   I tried starting a script in VSCode but was struggling, and didn't get far with the QGIS console either (and it seems janky).
#   so, I guess I'll just do the work in QGIS myself and return to an attempt at scripting later!


# STOPPING POINT
# - 20 August 2024


# # 'fortify' the data to get a dataframe format required by ggplot2 and other sundry packages
# grid_STTSTJ_df = fortify(grid_STTSTJ)
# df.PRICO_grid = fortify(PRICO_grid) #this may take forever
# df.STX_grid = fortify(STX_grid)
# object.size(df.STX_grid)/1000000 #this prints in megabytes, but will still say bytes
# object.size(PRICO_grid)/1000000
# rm(STTSTJ_grid) #clear up space; can also re-load these if needed
# rm(PRICO_grid)
# rm(STX_grid)


# # Sarah Groves note on grid: can resample to lo-res but then need to re-classify percent area of habitat (if of interest)
# #export shapefiles for Q
# shapefile(x = STTSTJ_grid, file = "/Users/benja/Documents/Carib_Habitat/Carib_Habitat_QGIS/Sampleframes/STTSTJ_2021_sample_frame.shp") #includes all metadata! very handy
# shapefile(x = PRICO_grid, file = "/Users/benja/Documents/Carib_Habitat/Carib_Habitat_QGIS/Sampleframes/PRICO_2021_sample_frame.shp") #includes all metadata! very handy
# shapefile(x = STX_grid, file = "/Users/benja/Documents/Carib_Habitat/Carib_Habitat_QGIS/Sampleframes/STX_2021_sample_frame.shp") #includes all metadata! very handy

# r = raster("final_merge_6Oct2022_NAs-for-land.tif")
# r = raster("final_merge_6Oct2022_zeros-for-land.tif")
bathymetry_50m = raster("final_merge_WGS1984.tif")

#resample to lower resolution using initial raster layer; defaults to mean value of cells (?)
bathymetry_650m = aggregate(bathymetry_50m, fact = 21.0964788648) #take ~30.18 res grid and resample to 650 m. napkin math-y
bathymetry_650m #not exact - 647 m res here (according to GIS). not sure how to address
# object.size(bathymetry_650m)/1000000 #this prints in megabytes, but will still say bytes

# plot(bathymetry_650m) #garb

df.bathymetry_50m = as.data.frame(bathymetry_50m, xy = T)
df.bathymetry_650m = as.data.frame(bathymetry_650m, xy = T)
# object.size(df.bathymetry_650m)/1000000 #this prints in megabytes, but will still say bytes

# potentially useful:
#https://www.benjaminbell.co.uk/2019/08/bathymetric-maps-in-r-colour-palettes.html
#https://stackoverflow.com/questions/20581746/increase-resolution-of-color-scale-for-values-close-to-zero

# #50 m - THIS EXCEEDS MEMORY ON A 16GB RAM LAPTOP
# ggplot() +
#   geom_tile(data = df.bathymetry_50m, aes(x = x, y = y, fill = final_merge_WGS1984)) +
#   scale_fill_gradientn(colors = terrain.colors(50), limits = c(-80, 0)) +
#   # scale_fill_gradientn(col=c(blue.col(50), terrain.colors(50)), breaks=s3,
#   #                        values = rescale(c(-10, -1, 0, 1, 10)),
#   #                        guide = "colorbar", limits=c(-10, 10))y
#   theme_bw() +
#   xlim(-65.2, -64.6) +
#   ylim(18.15, 18.45)

#650 m
#STTSTJ
ggplot() +
  geom_tile(data = df.bathymetry_650m, aes(x = x, y = y, fill = final_merge_WGS1984)) +
  scale_fill_gradientn(colors = terrain.colors(50), limits = c(-80, 0)) +
  # scale_fill_gradientn(col=c(blue.col(50), terrain.colors(50)), breaks=s3,
  #                        values = rescale(c(-10, -1, 0, 1, 10)),
  #                        guide = "colorbar", limits=c(-10, 10))y
  theme_bw() +
  xlim(-65.2, -64.6) +
  ylim(18.15, 18.45)

#PR - Vieques & Culebra
ggplot() +
  geom_tile(data = df.bathymetry_650m, aes(x = x, y = y, fill = final_merge_WGS1984)) +
  scale_fill_gradientn(colors = terrain.colors(50), limits = c(-80, 0)) +
  theme_bw() +
  xlim(-65.60, -64.5) +
  ylim(18.05, 18.60)

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


