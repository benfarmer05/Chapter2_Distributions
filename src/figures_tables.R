  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(terra)
  library(leaflet)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(tidyterra)
  library(osmdata)
  library(sf)
  library(tmap)
  library(ggrepel)
  
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
  
  # Project total_raster to lat/lon for plotting
  total_raster_latlon <- project(total_raster, "EPSG:4326")
  
  # Crop to entire domain once
  domain_extent <- ext(-67.99, -63.99, 17.63, 18.83)
  total_raster_cropped <- crop(total_raster_latlon, domain_extent)
  
  # ################################## Download and save OSM data ##################################
  # 
  # # DOWNLOAD OSM DATA - Run this section once, then comment out
  # {
  #   cat("Downloading OSM data for full domain...\n")
  # 
  #   # Define bounding box for the ENTIRE domain
  #   bbox_full <- c(-67.99, 17.63, -63.99, 18.83)
  # 
  #   # Get OSM data for islands
  #   cat("Downloading island data...\n")
  #   islands_osm <- opq(bbox = bbox_full) %>%
  #     add_osm_feature(key = "place", value = "island") %>%
  #     osmdata_sf()
  # 
  #   # Get OSM data for coastlines
  #   cat("Downloading coastline data...\n")
  #   coast_osm <- opq(bbox = bbox_full) %>%
  #     add_osm_feature(key = "natural", value = "coastline") %>%
  #     osmdata_sf()
  # 
  #   # Convert to terra SpatVector objects and combine
  #   cat("Processing land vectors...\n")
  #   land_parts <- list()
  # 
  #   if(!is.null(islands_osm$osm_polygons) && nrow(islands_osm$osm_polygons) > 0) {
  #     land_parts[[length(land_parts) + 1]] <- vect(islands_osm$osm_polygons)
  #   }
  #   if(!is.null(islands_osm$osm_multipolygons) && nrow(islands_osm$osm_multipolygons) > 0) {
  #     land_parts[[length(land_parts) + 1]] <- vect(islands_osm$osm_multipolygons)
  #   }
  #   if(!is.null(coast_osm$osm_polygons) && nrow(coast_osm$osm_polygons) > 0) {
  #     land_parts[[length(land_parts) + 1]] <- vect(coast_osm$osm_polygons)
  #   }
  #   if(!is.null(coast_osm$osm_multipolygons) && nrow(coast_osm$osm_multipolygons) > 0) {
  #     land_parts[[length(land_parts) + 1]] <- vect(coast_osm$osm_multipolygons)
  #   }
  # 
  #   # Combine all land parts
  #   land_vect <- do.call(rbind, land_parts)
  # 
  #   # Get island labels from OSM polygons
  #   cat("Processing island labels...\n")
  #   island_features <- bind_rows(
  #     islands_osm$osm_polygons,
  #     islands_osm$osm_multipolygons
  #   ) %>%
  #     filter(!is.na(name)) %>%
  #     select(name, geometry)
  # 
  #   # Convert to terra and get centroids for labels
  #   island_vect <- vect(island_features)
  #   island_centroids <- centroids(island_vect)
  #   island_df <- as.data.frame(island_centroids, geom = "XY")
  # 
  #   # Save processed OSM data
  #   cat("Saving OSM data...\n")
  #   saveRDS(land_vect, here("output", "osm_land_vect.rds"))
  #   saveRDS(island_df, here("output", "osm_island_labels.rds"))
  # 
  #   cat("OSM data saved successfully!\n")
  # }
  # 
  ################################## Load saved OSM data ##################################
  
  cat("Loading saved OSM data...\n")
  land_vect <- readRDS(here("output", "osm_land_vect.rds"))
  island_df <- readRDS(here("output", "osm_island_labels.rds"))
  cat("OSM data loaded successfully!\n")
  
  ################################## Prepare labels ##################################
  
  # # Define plot extent options (now in lat/lon - WGS84)
  # # Convert UTM extents to geographic coordinates for plotting
  # 
  # # # Helper function to convert UTM extent to lat/lon
  # # convert_extent_to_latlon <- function(xmin, xmax, ymin, ymax) {
  # #   extent_utm <- ext(xmin, xmax, ymin, ymax)
  # #   extent_vect <- as.polygons(extent_utm, crs = crs(bathy_final))
  # #   extent_latlon <- project(extent_vect, "EPSG:4326")
  # #   return(ext(extent_latlon))
  # # }
  # # 
  # # plot_extents = convert_extent_to_latlon(280000, 310000, 2010000, 2060000) #for investigating drops
  # # plot_extents = convert_extent_to_latlon(280000, 310000, 2000000, 2040000) #for investigating south of STT
  # # plot_extents = convert_extent_to_latlon(270000, 290000, 2000000, 2040000) #for investigating MCD
  # # plot_extents = convert_extent_to_latlon(276000, 284500, 2025000, 2040000) #for investigating MCD issue
  # # plot_extents = convert_extent_to_latlon(300000, 340000, 2000000, 2050000) #for investigating STJ
  # # plot_extents = convert_extent_to_latlon(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # # plot_extents = convert_extent_to_latlon(341000, 379000, 2057000, 2078000) # for investigating Anegada
  # # plot_extents = convert_extent_to_latlon(294000, 350000, 1950000, 1975000) #for investigating STX
  # # plot_extents = convert_extent_to_latlon(280000, 320000, 2000000, 2040000) #for investigating STT
  # # plot_extents = convert_extent_to_latlon(-30000, 0, 2000000, 2025000) #for investigating Mona Island
  # # plot_extents = convert_extent_to_latlon(20000, 40000, 2030000, 2045000) #for investigating Desecheo Island
  # # plot_extents = convert_extent_to_latlon(100000, 400000, 1950000, 2070000) #for looking at a lot of area
  # 
  # 
  # # STT/STJ (mid-Vieques to mid-Tortola)
  # plot_extents = ext(-65.39305, -64.63481, 18.06478, 18.46411)
  # 
  # # #STX
  # # # plot_extents = ext(-64.96248, -64.39552, 17.60975, 17.87752)
  # # plot_extents = ext(-64.9474418065718, -64.4094822804147, 17.6145841629496, 17.8561058115635)
  # 
  # # #whole domain
  # # # plot_extents = ext(-67.9895422167158, -64.1456687288866, 17.6345841629496, 18.8254789170551)
  # # plot_extents = ext(-67.9895422167158, -63.9885660051412, 17.6345841629496, 18.8254789170551)
  
  # Filter for specific islands to label
  islands_to_label <- c("Saint Thomas", "Jost Van Dyke", "Saint John", "Tortola", 
                        "Virgin Gorda", "Anegada", "Saint Croix", "Isla de Culebra", 
                        "Isla de Vieques", "Puerto Rico", "Mona Island", "Isla de Mona")
  
  island_df_filtered <- island_df %>%
    filter(name %in% islands_to_label)
  
  ################################## Define regions ##################################
  
  # Define different plotting regions
  regions <- list(
    stt_stj = list(
      xlim = c(-65.39305, -64.63481),
      ylim = c(18.06478, 18.46411),
      name = "STT_STJ"
    ),
    stx = list(
      xlim = c(-64.9474418065718, -64.4094822804147),
      ylim = c(17.6145841629496, 17.8561058115635),
      name = "STX"
    ),
    whole = list(
      xlim = c(-67.9895422167158, -63.9885660051412),
      ylim = c(17.6345841629496, 18.8254789170551),
      name = "Full_Domain"
    )
  )
  
  # Select which region to plot
  current_region <- regions$whole  # Change this to switch regions
  
  ################################## Create plot ##################################
  
  textsize <- 9
  titlesize <- 10
  
  # Crop raster to current region with a small buffer
  buffer <- 0.1  # degrees (adjust as needed)
  region_extent <- ext(
    current_region$xlim[1] - buffer,
    current_region$xlim[2] + buffer,
    current_region$ylim[1] - buffer,
    current_region$ylim[2] + buffer
  )
  total_raster_region <- crop(total_raster_cropped, region_extent)
  
  # Plot with ggrepel
  fig <- ggplot() +
    geom_spatraster(data = total_raster_region * 100, maxcell = Inf) +  # maxcell = Inf prevents downsampling
    scale_fill_viridis_c(name = "Coral Cover (%)", 
                         option = "turbo",
                         limits = c(0, 50),
                         na.value = "transparent") +
    geom_spatvector(data = land_vect, fill = "gray90", color = "black") +
    geom_text_repel(data = island_df_filtered, aes(x = x, y = y, label = name),
                    size = 4, fontface = "bold", family = 'Georgia',
                    min.segment.length = 0,
                    box.padding = 0.5) +
    coord_sf(xlim = current_region$xlim,
             ylim = current_region$ylim,
             expand = FALSE) +
    theme_classic(base_family = "Georgia") + 
    theme(legend.position = "right",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text = element_text(size = 8),
          legend.text = element_text(size = textsize),
          legend.title = element_text(size = titlesize))
  
  ################################## Display and save ##################################
  
  # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  quartz(h = 5, w = 7.087)
  
  fig
  
  # Save the Quartz output directly as a PDF
  quartz.save(file = here("output", paste0("fig_", current_region$name, ".pdf")), type = "pdf")
  
  # ggplot-export to image
  # ggsave(filename = here("output", paste0("fig_", current_region$name, ".png")), 
  #        device = "png", width = 7.087, height = 5, dpi = 1200)
  ggsave(filename = here("output", paste0("fig_", current_region$name, ".pdf")), 
         plot = fig, device = "pdf", width = 7.087, height = 5)
  
  # Close the Quartz device
  dev.off()



  # # NOTE - can reference the below, but I think above is way to go for now. but will need
  # #           to fix labels showing from regions outside domain
  # 
  # #
  # # .rs.restartR(clean = TRUE)
  # rm(list=ls())
  # 
  # library(here)
  # library(terra)
  # library(leaflet)
  # library(dplyr)
  # library(ggplot2)
  # library(RColorBrewer)
  # library(tidyterra)
  # library(osmdata)
  # library(sf)
  # library(tmap)
  # library(ggrepel)
  # 
  # source(here("src/functions.R"))
  # 
  # ################################## setup ##################################
  # 
  # #load species observations
  # load(here("output", "all_combined_data.rda"))
  # 
  # #load bathy
  # spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
  # bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
  # terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
  # 
  # #pull all final species-specific raster maps; unwrap & restore CRS
  # results_files <- list.files(here('output/output_maps'), 
  #                             pattern = '^results_light_.*\\.rds$', 
  #                             full.names = TRUE)
  # all_results <- lapply(results_files, function(f) {
  #   result <- readRDS(f)
  #   result$raster <- terra::unwrap(result$raster)
  #   result$raster_raw <- terra::unwrap(result$raster_raw)
  #   return(result)
  # })
  # names(all_results) <- gsub('results_light_|\\.rds', '', basename(results_files))
  # raster_list <- lapply(all_results, function(x) {
  #   r <- x$raster
  #   terra::crs(r) <- spatial_metadata$bathy_final$crs
  #   return(r)
  # })
  # 
  # #load the susceptibility-group summarized, and total, raster maps
  # load_spat_objects(directory = 'output/output_summarize_maps/')
  # source(here("src/functions.R"))
  # load(here("output/output_summarize_maps", "output_summarize_maps_workspace.Rdata"))
  # susc_rasters <- lapply(susc_rasters, function(x) {
  #   r <- terra::unwrap(x)
  #   terra::crs(r) <- spatial_metadata$bathy_final$crs
  #   return(r)
  # })
  # total_raster <- terra::unwrap(total_raster)
  # terra::crs(total_raster) <- spatial_metadata$bathy_final$crs
  # 
  # ################################## figure testing ##################################
  # 
  # # Define plot extent options (now in lat/lon - WGS84)
  # # Convert UTM extents to geographic coordinates for plotting
  # 
  # # # Helper function to convert UTM extent to lat/lon
  # # convert_extent_to_latlon <- function(xmin, xmax, ymin, ymax) {
  # #   extent_utm <- ext(xmin, xmax, ymin, ymax)
  # #   extent_vect <- as.polygons(extent_utm, crs = crs(bathy_final))
  # #   extent_latlon <- project(extent_vect, "EPSG:4326")
  # #   return(ext(extent_latlon))
  # # }
  # # 
  # # plot_extents = convert_extent_to_latlon(280000, 310000, 2010000, 2060000) #for investigating drops
  # # plot_extents = convert_extent_to_latlon(280000, 310000, 2000000, 2040000) #for investigating south of STT
  # # plot_extents = convert_extent_to_latlon(270000, 290000, 2000000, 2040000) #for investigating MCD
  # # plot_extents = convert_extent_to_latlon(276000, 284500, 2025000, 2040000) #for investigating MCD issue
  # # plot_extents = convert_extent_to_latlon(300000, 340000, 2000000, 2050000) #for investigating STJ
  # # plot_extents = convert_extent_to_latlon(220000, 260000, 2000000, 2010000) #for investigating Vieques
  # # plot_extents = convert_extent_to_latlon(341000, 379000, 2057000, 2078000) # for investigating Anegada
  # # plot_extents = convert_extent_to_latlon(294000, 350000, 1950000, 1975000) #for investigating STX
  # # plot_extents = convert_extent_to_latlon(280000, 320000, 2000000, 2040000) #for investigating STT
  # # plot_extents = convert_extent_to_latlon(-30000, 0, 2000000, 2025000) #for investigating Mona Island
  # # plot_extents = convert_extent_to_latlon(20000, 40000, 2030000, 2045000) #for investigating Desecheo Island
  # # plot_extents = convert_extent_to_latlon(100000, 400000, 1950000, 2070000) #for looking at a lot of area
  # 
  # 
  # # STT/STJ (mid-Vieques to mid-Tortola)
  # plot_extents = ext(-65.39305, -64.63481, 18.06478, 18.46411)
  # 
  # # #STX
  # # # plot_extents = ext(-64.96248, -64.39552, 17.60975, 17.87752)
  # # plot_extents = ext(-64.9474418065718, -64.4094822804147, 17.6145841629496, 17.8561058115635)
  # 
  # # #whole domain
  # # # plot_extents = ext(-67.9895422167158, -64.1456687288866, 17.6345841629496, 18.8254789170551)
  # # plot_extents = ext(-67.9895422167158, -63.9885660051412, 17.6345841629496, 18.8254789170551)
  # 
  # # Define bounding box
  # bbox <- c(xmin(plot_extents), ymin(plot_extents), 
  #           xmax(plot_extents), ymax(plot_extents))
  # 
  # # Get OSM data (still need osmdata/sf for this part)
  # islands_osm <- opq(bbox = bbox) %>%
  #   add_osm_feature(key = "place", value = "island") %>%
  #   osmdata_sf()
  # 
  # coast_osm <- opq(bbox = bbox) %>%
  #   add_osm_feature(key = "natural", value = "coastline") %>%
  #   osmdata_sf()
  # 
  # # Convert to terra SpatVector objects and combine
  # land_parts <- list()
  # 
  # if(!is.null(islands_osm$osm_polygons) && nrow(islands_osm$osm_polygons) > 0) {
  #   land_parts[[length(land_parts) + 1]] <- vect(islands_osm$osm_polygons)
  # }
  # if(!is.null(islands_osm$osm_multipolygons) && nrow(islands_osm$osm_multipolygons) > 0) {
  #   land_parts[[length(land_parts) + 1]] <- vect(islands_osm$osm_multipolygons)
  # }
  # if(!is.null(coast_osm$osm_polygons) && nrow(coast_osm$osm_polygons) > 0) {
  #   land_parts[[length(land_parts) + 1]] <- vect(coast_osm$osm_polygons)
  # }
  # if(!is.null(coast_osm$osm_multipolygons) && nrow(coast_osm$osm_multipolygons) > 0) {
  #   land_parts[[length(land_parts) + 1]] <- vect(coast_osm$osm_multipolygons)
  # }
  # 
  # # Combine all land parts
  # land_vect <- do.call(rbind, land_parts)
  # 
  # # Get island labels
  # island_labels <- NULL
  # if(!is.null(islands_osm$osm_polygons) && nrow(islands_osm$osm_polygons) > 0) {
  #   island_labels <- vect(islands_osm$osm_polygons[, c("name", "geometry")])
  # }
  # if(!is.null(islands_osm$osm_multipolygons) && nrow(islands_osm$osm_multipolygons) > 0) {
  #   if(is.null(island_labels)) {
  #     island_labels <- vect(islands_osm$osm_multipolygons[, c("name", "geometry")])
  #   } else {
  #     island_labels <- rbind(island_labels, 
  #                            vect(islands_osm$osm_multipolygons[, c("name", "geometry")]))
  #   }
  # }
  # 
  # # get centroids
  # if(!is.null(island_labels)) {
  #   # island_labels <- island_labels[island_labels$name == "Saint Thomas", ]
  #   island_centroids <- centroids(island_labels)
  # }
  # 
  # # Project total_raster to lat/lon for plotting
  # total_raster_latlon <- project(total_raster, "EPSG:4326")
  # 
  # # Crop raster to extent
  # total_raster_cropped <- crop(total_raster_latlon, plot_extents)
  # 
  # 
  # # Get island labels from OSM polygons
  # island_features <- bind_rows(
  #   islands_osm$osm_polygons,
  #   islands_osm$osm_multipolygons
  # ) %>%
  #   filter(!is.na(name)) %>%
  #   select(name, geometry)
  # 
  # # Convert to terra and get centroids for labels
  # island_vect <- vect(island_features)
  # island_centroids <- centroids(island_vect)
  # island_df <- as.data.frame(island_centroids, geom = "XY")
  # 
  # # Filter for specific islands
  # islands_to_label <- c("Saint Thomas", "Jost Van Dyke", "Saint John", "Tortola", 
  #                       "Virgin Gorda", "Anegada", "Saint Croix", "Isla de Culebra", 
  #                       "Isla de Vieques", "Puerto Rico", "Mona Island", "Isla de Mona")
  # 
  # island_df_filtered <- island_df %>%
  #   filter(name %in% islands_to_label)
  # 
  # textsize = 9
  # titlesize = 10 #9  #text sizes in ggplot are actually in units of points, when specified using element_text
  # 
  # # Plot with ggrepel
  # 
  # # NOTE / STOPPING POINT - 24 Oct 2025: using 'coord_sf' directly in the plotting function
  # #         might be the way to go to easily adjust plot bounds between areas. but should probably
  # #         go ahead and download/save the OSM info for PR, then I can just load it quickly at the start
  # #         of this script. Takes a while to download. then can ideally just use that one OSM set for all figures
  # 
  # fig = ggplot() +
  #   geom_spatraster(data = total_raster_cropped * 100) +
  #   scale_fill_viridis_c(name = "Coral Cover (%)", 
  #                        option = "turbo",
  #                        limits = c(0, 50),
  #                        na.value = "transparent") +
  #   geom_spatvector(data = land_vect, fill = "gray90", color = "black") +
  #   geom_text_repel(data = island_df_filtered, aes(x = x, y = y, label = name),
  #                   size = 4, fontface = "bold", family = 'Georgia',
  #                   min.segment.length = 0,
  #                   box.padding = 0.5) +
  #   # coord_sf(xlim = c(xmin(total_raster_cropped), xmax(total_raster_cropped)), 
  #   #      ylim = c(ymin(total_raster_cropped), ymax(total_raster_cropped)),
  #   #      expand = FALSE) +
  #   coord_sf(xlim = c(bbox[1], bbox[3]),
  #            ylim = c(bbox[2], bbox[4])) +
  #   # theme_minimal()
  #   theme_classic(base_family = "Georgia") + 
  #   theme(legend.position = "right",
  #         axis.title.x = element_blank(),
  #         axis.title.y = element_blank(),
  #         strip.text = element_text(size = 8),
  #         legend.text = element_text(size = textsize),
  #         legend.title = element_text(size = titlesize)) #,
  #         # legend.key.height = unit(0, "cm"))
  # 
  # 
  # # Set a standard plot size. max is 7.087 inch wide by 9.45 inch tall
  # # NOTE - can try windows() or x11() instead of Quartz in Windows and Linux, respectively. with appropriate downstream modifications as needed
  # # quartz(h = 5, w = 3.35)
  # quartz(h = 5, w = 7.087)
  # # quartz(h = 6, w = 5)
  # 
  # fig
  # 
  # # Save the Quartz output directly as a PDF
  # quartz.save(file = here("output", "fig.pdf"), type = "pdf")
  # 
  # #ggplot-export to image
  # ggsave(filename = here("output", "fig.png"), device = "png", width = 7.087, height = 5, dpi = 1200)
  # 
  # # Close the Quartz device
  # dev.off()
  # 
  # 
  # # # Alternative: Manual label control
  # # island_labels_manual <- data.frame(
  # #   name = c("Saint Thomas", "Saint John", "Water Island"),
  # #   lon = c(-64.93, -64.73, -64.97),
  # #   lat = c(18.34, 18.33, 18.32)
  # # )
  # # 
  # # ggplot() +
  # #   geom_spatraster(data = total_raster_cropped * 100) +
  # #   scale_fill_viridis_c(name = "Coral Cover (%)",
  # #                        option = "viridis",
  # #                        limits = c(0, 50),
  # #                        na.value = "transparent") +
  # #   geom_spatvector(data = land_vect, fill = "gray90", color = "black") +
  # #   geom_text(data = island_labels_manual, aes(x = lon, y = lat, label = name),
  # #             size = 4, fontface = "bold") +
  # #   coord_sf(xlim = c(bbox[1], bbox[3]),
  # #            ylim = c(bbox[2], bbox[4])) +
  # #   theme_minimal() +
  # #   ggtitle("Total Coral Cover with Coastlines")
  # 
  # 
  # 
  # # plot(total_raster_cropped * 100, 
  # #      range = c(0, 50),
  # #      col = hcl.colors(100, "viridis"),
  # #      axes = TRUE,
  # #      plg = list(title = "Coral Cover (%)"))
  # # 
  # # # Add land
  # # plot(land_vect, col = "gray90", border = "black", add = TRUE)
  # # 
  # # # Add labels
  # # if(!is.null(island_labels)) {
  # #   text(island_centroids, labels = island_centroids$name, 
  # #        font = 2, cex = 0.8)
  # # }
  # 
  # # # Create tmap (tmap works with SpatVector objects!)
  # # tmap_mode("plot")
  # # 
  # # tm_shape(total_raster_cropped * 100, bbox = plot_extents) +
  # #   tm_raster(style = "cont",
  # #             palette = "viridis",
  # #             breaks = seq(0, 50, 10),
  # #             title = "Coral Cover (%)") +
  # #   tm_shape(land_vect) +
  # #   tm_fill(col = "gray90") +
  # #   tm_borders(col = "black", lwd = 1) +
  # #   tm_shape(island_labels) +
  # #   tm_text("name", size = 0.8, col = "black", fontface = "bold") +
  # #   tm_layout(title = "Total Coral Cover with Coastlines",
  # #             legend.position = c("right", "center"),
  # #             frame = TRUE) +
  # #   tm_graticules() +
  # #   tm_compass(position = c("left", "top")) +
  # #   tm_scale_bar(position = c("left", "bottom"))
  # 
  # # # Plot with raster underneath
  # # ggplot() +
  # #   geom_spatraster(data = total_raster_cropped * 100) +
  # #   scale_fill_viridis_c(name = "Coral Cover (%)", 
  # #                        option = "viridis",
  # #                        limits = c(0, 50),
  # #                        na.value = "transparent") +
  # #   geom_spatvector(data = land_vect, fill = "gray90", color = "black") +
  # #   geom_spatvector_text(data = island_labels, aes(label = name), 
  # #                        size = 3, fontface = "bold") +
  # #   coord_sf(xlim = c(bbox[1], bbox[3]), 
  # #            ylim = c(bbox[2], bbox[4])) +
  # #   theme_minimal() +
  # #   ggtitle("Total Coral Cover with Coastlines")
  # 
  # 