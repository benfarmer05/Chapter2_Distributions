library(leaflet)

# Define your extents - comment out the ones you don't want to use
# plot_extents <- ext(280000, 340000, 2000000, 2050000)  # STTSTJ
# plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
# plot_extents = ext(200000, 400000, 2000000, 2040000) #for investigating MCD
plot_extents = ext(100000, 400000, 1950000, 2070000) #for looking at a lot of area
plot_extents = ext(266000, 345000, 1997000, 2060000) #for looking at a bit of area

# --- SLOPE ---
slope_terra_clamp <- clamp(slope_terra, lower = 0, upper = 20)
slope_cropped <- crop(slope_terra_clamp, plot_extents)
pal_slope <- colorNumeric(cmocean("deep")(100), values(slope_cropped), na.color = "transparent")

slope_map <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(slope_cropped, colors = pal_slope, opacity = 0.7) %>%
  addLegend("bottomright", pal = pal_slope, values = values(slope_cropped),
            title = "Slope (0-20)") %>%
  addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                   options = layersControlOptions(collapsed = FALSE))

# --- SLOPE OF SLOPE ---
slopeofslope_clamp <- clamp(slopeofslope_terra, lower = 0, upper = 15)
slopeofslope_cropped <- crop(slopeofslope_clamp, plot_extents)
pal_slopeofslope <- colorNumeric(cmocean("deep")(100), values(slopeofslope_cropped), na.color = "transparent")

slopeofslope_map <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(slopeofslope_cropped, colors = pal_slopeofslope, opacity = 0.7) %>%
  addLegend("bottomright", pal = pal_slopeofslope, values = values(slopeofslope_cropped),
            title = "Slope of Slope") %>%
  addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                   options = layersControlOptions(collapsed = FALSE))

# --- SAPA ---
SAPA_clamp <- clamp(SAPA, lower = 1, upper = 1.01)
SAPA_cropped <- crop(SAPA_clamp, plot_extents)
pal_SAPA <- colorNumeric(cmocean("deep")(100), values(SAPA_cropped), na.color = "transparent")

SAPA_map <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(SAPA_cropped, colors = pal_SAPA, opacity = 0.7) %>%
  addLegend("bottomright", pal = pal_SAPA, values = values(SAPA_cropped),
            title = "SAPA (1-1.01)") %>%
  addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                   options = layersControlOptions(collapsed = FALSE))

# --- ROUGHNESS ---
roughness_cropped <- crop(roughness, plot_extents)
pal_roughness <- colorNumeric(cmocean("deep")(100), values(roughness_cropped), na.color = "transparent")

roughness_map <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(roughness_cropped, colors = pal_roughness, opacity = 0.7) %>%
  addLegend("bottomright", pal = pal_roughness, values = values(roughness_cropped),
            title = "Roughness") %>%
  addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                   options = layersControlOptions(collapsed = FALSE))

# --- VRM ---
VRM_clamp <- clamp(VRM, lower = -0.01, upper = 0.01)
VRM_cropped <- crop(VRM_clamp, plot_extents)
pal_VRM <- colorNumeric(cmocean("balance")(100), values(VRM_cropped), na.color = "transparent")

VRM_map <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(VRM_cropped, colors = pal_VRM, opacity = 0.7) %>%
  addLegend("bottomright", pal = pal_VRM, values = values(VRM_cropped),
            title = "VRM (-0.01 to 0.01)") %>%
  addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                   options = layersControlOptions(collapsed = FALSE))

# --- TPI TERRA ---
TPI_terra_clamp <- clamp(TPI_terra, lower = -5, upper = 5)
TPI_cropped <- crop(TPI_terra_clamp, plot_extents)
pal_TPI <- colorNumeric(cmocean("balance")(100), values(TPI_cropped), na.color = "transparent")

TPI_map <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(TPI_cropped, colors = pal_TPI, opacity = 0.7) %>%
  addLegend("bottomright", pal = pal_TPI, values = values(TPI_cropped),
            title = "TPI (-40 to 40)") %>%
  addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                   options = layersControlOptions(collapsed = FALSE))

# --- PLANFORM CURVATURE ---
planformcurv_multiscale_clamp <- clamp(planformcurv_multiscale, lower = -0.005, upper = 0.005)
planformcurv_cropped <- crop(planformcurv_multiscale_clamp, plot_extents)
pal_planformcurv <- colorNumeric(cmocean("balance")(100), values(planformcurv_cropped), na.color = "transparent")

planformcurv_map <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(planformcurv_cropped, colors = pal_planformcurv, opacity = 0.7) %>%
  addLegend("bottomright", pal = pal_planformcurv, values = values(planformcurv_cropped),
            title = "Planform Curvature") %>%
  addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                   options = layersControlOptions(collapsed = FALSE))

# --- PROFILE CURVATURE ---
profilecurv_multiscale_clamp <- clamp(profilecurv_multiscale, lower = -0.01, upper = 0.01)
profilecurv_cropped <- crop(profilecurv_multiscale_clamp, plot_extents)
pal_profilecurv <- colorNumeric(cmocean("balance")(100), values(profilecurv_cropped), na.color = "transparent")

profilecurv_map <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(profilecurv_cropped, colors = pal_profilecurv, opacity = 0.7) %>%
  addLegend("bottomright", pal = pal_profilecurv, values = values(profilecurv_cropped),
            title = "Profile Curvature") %>%
  addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                   options = layersControlOptions(collapsed = FALSE))

# --- BATHYMETRY ---
depth_OG_clamp <- clamp(bathy_final, lower = -50, upper = 0)
bathy_cropped <- crop(depth_OG_clamp, plot_extents)
pal_bathy <- colorNumeric(rev(cmocean("deep")(100)), values(bathy_cropped), na.color = "transparent")

bathy_map <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(bathy_cropped, colors = pal_bathy, opacity = 0.7) %>%
  addLegend("bottomright", pal = pal_bathy, values = values(bathy_cropped),
            title = "Depth (m)") %>%
  addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                   options = layersControlOptions(collapsed = FALSE))

# Display maps (run each individually)
slope_map
slopeofslope_map
SAPA_map
roughness_map
VRM_map
TPI_map
planformcurv_map
profilecurv_map
bathy_map