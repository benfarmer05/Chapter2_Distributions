 
   # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(terra)
  library(leaflet)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(tidyterra)
  
  source(here("src/functions.R"))
  
  ################################## setup ##################################
  
  #load species observations
  load(here("output", "all_combined_data.rda"))
  
  #load bathy
  spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
  bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
  terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
  
  #pull all final species-specific raster maps; unwrap & restore CRS
  results_files <- list.files(here('output/output_maps'), 
                              pattern = '^results_light_.*\\.rds$', 
                              full.names = TRUE)
  all_results <- lapply(results_files, function(f) {
    result <- readRDS(f)
    result$raster <- terra::unwrap(result$raster)
    result$raster_raw <- terra::unwrap(result$raster_raw)
    return(result)
  })
  names(all_results) <- gsub('results_light_|\\.rds', '', basename(results_files))
  raster_list <- lapply(all_results, function(x) {
    r <- x$raster
    terra::crs(r) <- spatial_metadata$bathy_final$crs
    return(r)
  })
  
  #load the susceptibility-group summarized, and total, raster maps
  load_spat_objects(directory = 'output/output_summarize_maps/')
  source(here("src/functions.R"))
  load(here("output/output_summarize_maps", "output_summarize_maps_workspace.Rdata"))
  susc_rasters <- lapply(susc_rasters, function(x) {
    r <- terra::unwrap(x)
    terra::crs(r) <- spatial_metadata$bathy_final$crs
    return(r)
  })
  total_raster <- terra::unwrap(total_raster)
  terra::crs(total_raster) <- spatial_metadata$bathy_final$crs
  
  ################################## figure testing ##################################
  
  # Define plot extent options (now in lat/lon - WGS84)
  # Convert UTM extents to geographic coordinates for plotting
  
  # Helper function to convert UTM extent to lat/lon
  convert_extent_to_latlon <- function(xmin, xmax, ymin, ymax) {
    extent_utm <- ext(xmin, xmax, ymin, ymax)
    extent_vect <- as.polygons(extent_utm, crs = crs(bathy_final))
    extent_latlon <- project(extent_vect, "EPSG:4326")
    return(ext(extent_latlon))
  }
  
  # All your extent options converted to lat/lon
  plot_extents = convert_extent_to_latlon(280000, 310000, 2010000, 2060000) #for investigating drops
  plot_extents = convert_extent_to_latlon(280000, 310000, 2000000, 2040000) #for investigating south of STT
  plot_extents = convert_extent_to_latlon(270000, 290000, 2000000, 2040000) #for investigating MCD
  plot_extents = convert_extent_to_latlon(276000, 284500, 2025000, 2040000) #for investigating MCD issue
  plot_extents = convert_extent_to_latlon(300000, 340000, 2000000, 2050000) #for investigating STJ
  plot_extents = convert_extent_to_latlon(220000, 260000, 2000000, 2010000) #for investigating Vieques
  plot_extents = convert_extent_to_latlon(341000, 379000, 2057000, 2078000) # for investigating Anegada
  plot_extents = convert_extent_to_latlon(294000, 350000, 1950000, 1975000) #for investigating STX
  plot_extents = convert_extent_to_latlon(280000, 320000, 2000000, 2040000) #for investigating STT
  plot_extents = convert_extent_to_latlon(-30000, 0, 2000000, 2025000) #for investigating Mona Island
  plot_extents = convert_extent_to_latlon(20000, 40000, 2030000, 2045000) #for investigating Desecheo Island
  plot_extents = convert_extent_to_latlon(100000, 400000, 1950000, 2070000) #for looking at a lot of area
  
  
  plot_extents = ext(-66.7911391694081, -63.9426626207812, 17.60079552834, 18.7189860077195) #for investigating drops
  plot_extents = ext(-65.0830158006712, -64.7953619345193, 18.0775196566442, 18.4417215936095) #for investigating south of STT
  # plot_extents = ext(270000, 290000, 2000000, 2040000) #for investigating MCD
  # plot_extents = ext(276000, 284500, 2025000, 2040000) #for investigating MCD issue
  # plot_extents = ext(300000, 340000, 2000000, 2050000) #for investigating STJ
  # plot_extents = ext(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # plot_extents = ext(341000, 379000, 2057000, 2078000) # for investigating Anegada
  # plot_extents = ext(294000, 350000, 1950000, 1975000) #for investigating STX
  # plot_extents = ext(280000, 320000, 2000000, 2040000) #for investigating STT
  # plot_extents = ext(-30000, 0, 2000000, 2025000) #for investigating Mona Island
  # plot_extents = ext(20000, 40000, 2030000, 2045000) #for investigating Desecheo Island
  plot_extents = ext(100000, 400000, 1950000, 2070000) #for looking at a lot of area
  
  
  
  
  
  
  # Project total_raster to lat/lon for plotting
  total_raster_latlon <- project(total_raster, "EPSG:4326")
  
  # Clamp cover values to 0-50% for observation points
  psu_total_cover$total_cover_clamped <- pmin(psu_total_cover$total_cover, 50)
  
  # Create the plot
  p <- ggplot() +
    # Raster layer - directly from SpatRaster!
    # Multiply by 100 to convert 0-1 to 0-100%
    
    geom_spatraster(data = total_raster_latlon * 100,
                    aes(fill = hurdle_pred),
                    maxcell = 4e6) +
    
    # geom_spatraster(data = total_raster_latlon * 100,
    #                 aes(fill = hurdle_pred)) + #default downsampling
    
    # # Add observation points (colored by cover, fixed size)
    # geom_point(data = psu_total_cover, 
    #            aes(x = lon, y = lat, fill = total_cover_clamped),
    #            shape = 21, color = "white", size = 1.4,
    #            stroke = 0.3, alpha = 0.9) +
    
    # viridis turbo color scale for both raster and points
    scale_fill_viridis_c(option = "turbo",
                         name = "Coral Cover (%)",
                         limits = c(0, 50),
                         breaks = seq(0, 50, 10),
                         na.value = "transparent") +
    
    # Labels
    labs(title = "Total Coral Cover Across All Species",
         subtitle = "Model predictions (colored raster) with field observations (circles)",
         x = "Longitude (°W)",
         y = "Latitude (°N)") +
    
    # Coordinate system with geographic extent
    coord_sf(xlim = c(xmin(plot_extents), xmax(plot_extents)),
             ylim = c(ymin(plot_extents), ymax(plot_extents)),
             expand = FALSE) +
    
    # Clean theme for publication
    theme_bw(base_size = 12) +
    theme(
      # Panel
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      
      # Plot
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      plot.margin = margin(10, 10, 10, 10),
      
      # Axes
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(face = "bold", size = 11),
      axis.ticks = element_line(color = "black"),
      
      # Legend
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key = element_rect(fill = "white"),
      legend.spacing.y = unit(0.3, "cm")
    )
  
  print(p)
  
  
  
  
  
  
  
  
  
  
  
  library(leaflet)
  
  # Crop the raster to plot_extents
  total_raster_cropped <- crop(total_raster_latlon, plot_extents)
  
  # Clamp the raster values to 0-50% (converting from 0-1 to 0-0.5)
  total_raster_cropped_clamped <- clamp(total_raster_cropped, lower = 0, upper = 0.5, values = TRUE)
  
  # Clamp observation cover values
  psu_total_cover$total_cover_clamped <- pmin(psu_total_cover$total_cover, 50)
  
  # Filter observation points to only those within plot_extents
  psu_total_cover_filtered <- psu_total_cover %>%
    filter(lon >= xmin(plot_extents) & lon <= xmax(plot_extents) &
             lat >= ymin(plot_extents) & lat <= ymax(plot_extents))
  
  # Create color palette for 0-50% range
  pal <- colorNumeric("turbo", domain = c(0, 50), na.color = "transparent")
  
  # Create leaflet map
  leaflet_map <- leaflet() %>%
    addTiles(group = "OpenStreetMap") %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    
    # Add cropped raster (multiply by 100 to convert to percentage)
    addRasterImage(total_raster_cropped_clamped * 100, 
                   colors = pal, 
                   opacity = 0.7) %>%
    
    # Add observation points (filtered)
    addCircleMarkers(data = psu_total_cover_filtered, 
                     ~lon, ~lat, 
                     radius = 6,
                     fillColor = ~pal(total_cover_clamped), 
                     fillOpacity = 0.8,
                     color = "white", 
                     weight = 2, 
                     stroke = TRUE,
                     popup = ~paste0("<b>PSU:</b> ", PSU, 
                                     "<br><b>Cover:</b> ", round(total_cover, 2), "%")) %>%
    
    # Add legend
    addLegend("bottomright", 
              pal = pal, 
              values = seq(0, 50, by = 10),
              title = "Total Coral Cover (%)",
              opacity = 1) %>%
    
    # Set view to plot_extents
    fitBounds(lng1 = xmin(plot_extents), 
              lat1 = ymin(plot_extents),
              lng2 = xmax(plot_extents), 
              lat2 = ymax(plot_extents)) %>%
    
    # Layer control
    addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                     options = layersControlOptions(collapsed = FALSE))
  
  # Display the map
  leaflet_map
  
  
  
  
  
  
  
  
  
  
  library(tmap)
  
  # Crop to plot extents
  total_raster_cropped <- crop(total_raster_latlon, plot_extents)
  
  # Create map - just the raster
  tm <- tm_shape(total_raster_cropped * 100) +
    tm_raster(style = "cont",
              palette = "turbo",
              breaks = seq(0, 50, 10),
              title = "Coral Cover (%)") +
    tm_layout(title = "Total Coral Cover Across All Species",
              legend.position = c("right", "center"),
              frame = TRUE) +
    tm_graticules() +
    tm_compass(position = c("left", "top")) +
    tm_scale_bar(position = c("left", "bottom"))
  
  # Display
  tmap_mode("plot")
  tm
  
  # Save high quality
  tmap_save(tm, "coral_cover_map.png", width = 12, height = 10, dpi = 300)
  
  
  
  ################################## better test I think ##################################
  
  
  library(osmdata)
  library(sf)
  library(tmap)
  library(terra)
  
  # Define your bounding box
  bbox <- c(xmin(plot_extents), ymin(plot_extents), 
            xmax(plot_extents), ymax(plot_extents))
  
  # Get islands
  islands_osm <- opq(bbox = bbox) %>%
    add_osm_feature(key = "place", value = "island") %>%
    osmdata_sf()
  
  # Get coastline
  coast_osm <- opq(bbox = bbox) %>%
    add_osm_feature(key = "natural", value = "coastline") %>%
    osmdata_sf()
  
  # Extract just geometries for land shapes
  land_parts <- list()
  if(!is.null(islands_osm$osm_polygons) && nrow(islands_osm$osm_polygons) > 0) {
    land_parts[[length(land_parts) + 1]] <- st_geometry(islands_osm$osm_polygons)
  }
  if(!is.null(islands_osm$osm_multipolygons) && nrow(islands_osm$osm_multipolygons) > 0) {
    land_parts[[length(land_parts) + 1]] <- st_geometry(islands_osm$osm_multipolygons)
  }
  if(!is.null(coast_osm$osm_polygons) && nrow(coast_osm$osm_polygons) > 0) {
    land_parts[[length(land_parts) + 1]] <- st_geometry(coast_osm$osm_polygons)
  }
  if(!is.null(coast_osm$osm_multipolygons) && nrow(coast_osm$osm_multipolygons) > 0) {
    land_parts[[length(land_parts) + 1]] <- st_geometry(coast_osm$osm_multipolygons)
  }
  
  land_sf <- st_sf(geometry = do.call(c, land_parts))
  
  # Get island names with centroids for labels
  island_labels <- rbind(
    if(!is.null(islands_osm$osm_polygons) && nrow(islands_osm$osm_polygons) > 0) {
      islands_osm$osm_polygons[, c("name", "geometry")]
    },
    if(!is.null(islands_osm$osm_multipolygons) && nrow(islands_osm$osm_multipolygons) > 0) {
      islands_osm$osm_multipolygons[, c("name", "geometry")]
    }
  )
  
  # Keep only Saint Thomas
  island_labels <- island_labels[island_labels$name == "Saint Thomas", ]
  island_labels <- st_centroid(island_labels)
  
  # Crop raster to extent
  total_raster_cropped <- crop(total_raster_latlon, plot_extents)
  
  # Create tmap
  tmap_mode("plot")
  
  tm_shape(total_raster_cropped * 100, bbox = plot_extents) +
    tm_raster(style = "cont",
              palette = "viridis",
              breaks = seq(0, 50, 10),
              title = "Coral Cover (%)") +
    tm_shape(land_sf) +
    tm_fill(col = "gray90") +
    tm_borders(col = "black", lwd = 1) +
    tm_shape(island_labels) +
    tm_text("name", size = 0.8, col = "black", fontface = "bold") +
    tm_layout(title = "Total Coral Cover with Coastlines",
              legend.position = c("right", "center"),
              frame = TRUE) +
    tm_graticules() +
    tm_compass(position = c("left", "top")) +
    tm_scale_bar(position = c("left", "bottom"))  
  
  
  
  
  
  
  library(osmdata)
  library(sf)
  library(ggplot2)
  library(terra)
  library(tidyterra)
  
  # Define your bounding box
  bbox <- c(xmin(plot_extents), ymin(plot_extents), 
            xmax(plot_extents), ymax(plot_extents))
  
  # Get islands
  islands_osm <- opq(bbox = bbox) %>%
    add_osm_feature(key = "place", value = "island") %>%
    osmdata_sf()
  
  # Get coastline
  coast_osm <- opq(bbox = bbox) %>%
    add_osm_feature(key = "natural", value = "coastline") %>%
    osmdata_sf()
  
  # Extract just geometries for land shapes
  land_parts <- list()
  if(!is.null(islands_osm$osm_polygons) && nrow(islands_osm$osm_polygons) > 0) {
    land_parts[[length(land_parts) + 1]] <- st_geometry(islands_osm$osm_polygons)
  }
  if(!is.null(islands_osm$osm_multipolygons) && nrow(islands_osm$osm_multipolygons) > 0) {
    land_parts[[length(land_parts) + 1]] <- st_geometry(islands_osm$osm_multipolygons)
  }
  if(!is.null(coast_osm$osm_polygons) && nrow(coast_osm$osm_polygons) > 0) {
    land_parts[[length(land_parts) + 1]] <- st_geometry(coast_osm$osm_polygons)
  }
  if(!is.null(coast_osm$osm_multipolygons) && nrow(coast_osm$osm_multipolygons) > 0) {
    land_parts[[length(land_parts) + 1]] <- st_geometry(coast_osm$osm_multipolygons)
  }
  
  land_sf <- st_sf(geometry = do.call(c, land_parts))
  
  # Get island names with centroids for labels
  island_labels <- rbind(
    if(!is.null(islands_osm$osm_polygons) && nrow(islands_osm$osm_polygons) > 0) {
      islands_osm$osm_polygons[, c("name", "geometry")]
    },
    if(!is.null(islands_osm$osm_multipolygons) && nrow(islands_osm$osm_multipolygons) > 0) {
      islands_osm$osm_multipolygons[, c("name", "geometry")]
    }
  )
  
  # Keep only Saint Thomas
  island_labels <- island_labels[island_labels$name == "Saint Thomas", ]
  island_labels <- st_centroid(island_labels)
  
  # Crop raster to extent
  total_raster_cropped <- crop(total_raster_latlon, plot_extents)
  
  # Plot with raster underneath
  ggplot() +
    geom_spatraster(data = total_raster_cropped * 100) +
    scale_fill_viridis_c(name = "Coral Cover (%)", 
                         option = "turbo",
                         limits = c(0, 50),
                         na.value = "transparent") +
    geom_sf(data = land_sf, fill = "gray90", color = "black") +
    geom_sf_text(data = island_labels, aes(label = name), 
                 size = 3, fontface = "bold") +
    coord_sf(xlim = c(bbox[1], bbox[3]), 
             ylim = c(bbox[2], bbox[4])) +
    theme_minimal() +
    ggtitle("Total Coral Cover with Coastlines")