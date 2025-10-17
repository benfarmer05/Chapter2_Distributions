    
  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(terra)
  library(sf)
  library(ggplot2)
  library(leaflet)
  # library(dplyr)
  # library(RColorBrewer)
  
  source(here("src/functions.R"))
  
  
  ################################## setup ##################################
  
  #load 650 m habitat grid initially created from 0-60 m depth binary
  load_spat_objects(directory = 'output/output_create_habitat_grid/')
  source(here("src/functions.R"))
  load(here('output', 'output_create_habitat_grid/create_habitat_grid_workspace.RData'))
  
  #load bathy
  spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
  bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
  terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
  
  # Get all GAM map results files
  results_files <- list.files(here('output/output_maps'), 
                              pattern = '^results_light_.*\\.rds$', 
                              full.names = TRUE)
  
  # Read all results and unwrap rasters
  all_results <- lapply(results_files, function(f) {
    result <- readRDS(f)
    result$raster <- terra::unwrap(result$raster)
    result$raster_raw <- terra::unwrap(result$raster_raw)
    return(result)
  })
  names(all_results) <- gsub('results_light_|\\.rds', '', basename(results_files))
  
  ################################## restore CRS and sum rasters ##################################
  
  cat("Restoring CRS and calculating total cover...\n")
  
  # Extract all rasters and restore CRS to each
  raster_list <- lapply(all_results, function(x) {
    r <- x$raster
    terra::crs(r) <- spatial_metadata$bathy_final$crs
    return(r)
  })
  
  # Stack all rasters and sum them
  total_raster <- Reduce(`+`, raster_list)
  
  cat("Total raster created\n\n")
  
  
  ################################## check grid & centroid IDs match ##################################
  
  grid_SDMreef_no_land = grid_reefy_no_land
  
  # Extract and sort IDs
  grid_ids <- sort(grid_SDMreef_no_land$unique_ID)
  centroid_ids <- sort(centroids_all$unique_ID)
  
  # Check if identical
  identical(grid_ids, centroid_ids)
  
  # Quick visual check: random sample of 12
  set.seed(123)
  sample_ids <- sample(grid_ids, 12)
  
  grid_sample <- grid_SDMreef_no_land[grid_SDMreef_no_land$unique_ID %in% sample_ids, ]
  cent_sample <- centroids_all[centroids_all$unique_ID %in% sample_ids, ]
  
  # Plot
  plot(grid_sample, col = NA, border = "blue", lwd = 2)
  plot(cent_sample, col = "red", pch = 16, cex = 1.5, add = TRUE)
  text(crds(cent_sample), labels = cent_sample$unique_ID, pos = 4, cex = 0.7, col = "red")
  text(centroids(grid_sample), labels = grid_sample$unique_ID, pos = 2, cex = 0.7, col = "blue")
  legend("topright", legend = c("Grid squares", "Centroids"), 
         col = c("blue", "red"), lwd = c(2, NA), pch = c(NA, 16))
  
  
  ################################## filter 650m grid to predicted coral presence areas ##################################
  
  # Extract mean coral cover for each grid cell
  grid_coral_values <- extract(total_raster, grid_SDMreef_no_land, fun = mean, na.rm = TRUE)
  
  # Add mean coral cover to the grid
  grid_SDMreef_no_land$mean_coral_cover <- grid_coral_values[, 2]
  
  # Filter to cells with mean cover > X threshold (could do 0.005 or 0.01)
  has_coral <- grid_SDMreef_no_land$mean_coral_cover > 0.01 & !is.na(grid_SDMreef_no_land$mean_coral_cover)
  grid_SDMreef_no_land <- grid_SDMreef_no_land[has_coral]
  
  # # Filter centroids to match
  # centroids_with_coral <- centroids_all[centroids_all$unique_ID %in% grid_SDMreef_no_land$unique_ID]
  centroids_SDM_reef_no_land <- centroids_all[match(grid_SDMreef_no_land$unique_ID, centroids_all$unique_ID)]
  
  # Verify matching
  cat("  Original grid cells:", nrow(grid_SDMreef_no_land), "\n")
  cat("  Grid cells with predicted coral:", nrow(grid_SDMreef_no_land), "\n")
  cat("  Matching centroids:", nrow(centroids_SDM_reef_no_land), "\n")
  cat("  Reduction:", round((1 - nrow(grid_SDMreef_no_land)/nrow(grid_SDMreef_no_land)) * 100, 1), "%\n")
  cat("  Mean coral cover range:", round(range(grid_SDMreef_no_land$mean_coral_cover, na.rm=TRUE), 3), "\n")
  cat("  IDs match:", all(grid_SDMreef_no_land$unique_ID == centroids_SDM_reef_no_land$unique_ID), "\n\n")
  
  # Visualization with coral cover gradient
  grid_SDMreef_no_land_sf <- st_as_sf(grid_SDMreef_no_land)
  
  ggplot() +
    # geom_sf(data = st_as_sf(grid_SDMreef_no_land), fill = "lightgray", color = "gray80", size = 0.2) +
    geom_sf(data = grid_SDMreef_no_land_sf, aes(fill = mean_coral_cover * 100), color = "black", size = 0.1) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Mean Coral\nCover (%)") +
    geom_sf(data = st_as_sf(centroids_SDM_reef_no_land), color = "blue", size = 0.3, alpha = 0.5) +
    theme_minimal() +
    labs(title = "650m Grid: Mean Coral Cover per Cell",
         subtitle = paste0(nrow(grid_SDMreef_no_land), " cells with predicted coral"))
  
  
  
  
  # Optional: Save filtered outputs with coral cover attribute
  # writeVector(grid_SDMreef_no_land, here("output", "grid_SDMreef_no_land.shp"), overwrite = TRUE)
  # writeVector(centroids_SDM_reef_no_land, here("output", "centroids_SDM_reef_no_land"), overwrite = TRUE)
  
  
  # # Filter for coral cover > 0 but < a low threshold
  # low_coral <- grid_SDMreef_no_land[grid_SDMreef_no_land$mean_coral_cover > 0 & 
  #                                grid_SDMreef_no_land$mean_coral_cover < 0.01, ]
  # 
  # cat("Cells with 0% < coral cover < threshold %:", nrow(low_coral), "\n")
  # 
  # # Plot
  # low_coral_sf <- st_as_sf(low_coral)
  # 
  # ggplot() +
  #   geom_sf(data = st_as_sf(grid_SDMreef_no_land), fill = "lightgray", color = "gray80", size = 0.2) +
  #   geom_sf(data = low_coral_sf, aes(fill = mean_coral_cover * 100), color = "black", size = 0.3) +
  #   scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Mean Coral\nCover (%)") +
  #   theme_minimal() +
  #   labs(title = "Grid Cells with 0% < Coral Cover < 0.5%",
  #        subtitle = paste0(nrow(low_coral), " cells"))
  
  
  
  ################################## leaflet plot the grid ##################################
  
  # Transform to WGS84 for leaflet
  grid_SDMreef_no_land_wgs84 <- st_transform(grid_SDMreef_no_land_sf, 4326)
  centroids_with_coral_wgs84 <- st_transform(st_as_sf(centroids_SDM_reef_no_land), 4326)
  
  # Create color palette with viridis
  coral_pal <- colorNumeric(palette = "viridis", domain = grid_SDMreef_no_land_wgs84$mean_coral_cover * 100)
  
  # Create leaflet map
  leaflet() %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addTiles(group = "OpenStreetMap") %>%
    addPolygons(data = grid_SDMreef_no_land_wgs84,
                fillColor = ~coral_pal(mean_coral_cover * 100),
                fillOpacity = 0.6,
                color = "black",
                weight = 1,
                popup = ~paste0("<b>ID:</b> ", unique_ID, "<br><b>Coral Cover:</b> ", 
                                round(mean_coral_cover * 100, 2), "%")) %>%
    addCircleMarkers(data = centroids_with_coral_wgs84,
                     radius = 3,
                     color = "blue",
                     fillOpacity = 0.5,
                     stroke = FALSE) %>%
    addLegend("bottomright", pal = coral_pal, 
              values = grid_SDMreef_no_land_wgs84$mean_coral_cover * 100,
              title = "Mean Coral<br>Cover (%)",
              opacity = 1) %>%
    addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
                     options = layersControlOptions(collapsed = FALSE))  
  # ##### TEST to simply resample cover raster (we don't want this I don't think)
  # 
  # # Convert grid to raster to use as template
  # grid_rast <- rast(grid_reefy_no_land, res = 650)
  # 
  # # Resample total_raster to match the grid
  # total_resampled <- resample(total_raster, grid_rast, method = "bilinear")
  # 
  # # Plot for comparison
  # par(mfrow = c(1, 2))
  # plot(total_raster, main = "Original")
  # plot(total_resampled, main = "Resampled to Grid")
  
  
  
  ################################## write new 650_m csv for MATLAB / CMS ##################################
  
  geographic_crs <- crs("EPSG:4269")
  
  # Transform centroids to geographic coordinates
  centroids_FINALFORCMS_localCRS = centroids_SDM_reef_no_land
  grid_FINALFORCMS_localCRS = grid_SDMreef_no_land
  
  grid_FINALFORCMS_geographic = project(grid_FINALFORCMS_localCRS, geographic_crs)
  
  centroids_FINALFORCMS_geographic <- project(centroids_FINALFORCMS_localCRS, geographic_crs)
  centroids_FINALFORCMS_geographic <- centroids_FINALFORCMS_geographic[match(grid_FINALFORCMS_localCRS$unique_ID, centroids_FINALFORCMS_geographic$unique_ID), ]
  
  # Add vertices information
  # Function to extract vertices from each polygon
  get_vertices <- function(spatvec) {
    # Get coordinates for each polygon
    coords_list <- lapply(1:nrow(spatvec), function(i) {
      geom_coords <- geom(spatvec[i])[, c("x", "y")]
      # Get first 4 vertices (excluding the closing vertex which repeats the first)
      if(nrow(geom_coords) >= 5) {
        geom_coords[1:4, ]
      } else {
        geom_coords[1:4, ]
      }
    })
    
    # Create columns for each vertex
    v1_x <- sapply(coords_list, function(x) x[1, 1])
    v1_y <- sapply(coords_list, function(x) x[1, 2])
    v2_x <- sapply(coords_list, function(x) x[2, 1])
    v2_y <- sapply(coords_list, function(x) x[2, 2])
    v3_x <- sapply(coords_list, function(x) x[3, 1])
    v3_y <- sapply(coords_list, function(x) x[3, 2])
    v4_x <- sapply(coords_list, function(x) x[4, 1])
    v4_y <- sapply(coords_list, function(x) x[4, 2])
    
    return(data.frame(v1_x, v1_y, v2_x, v2_y, v3_x, v3_y, v4_x, v4_y))
  }
  
  # Add vertices to geographic version (we want geographic coords for output)
  vertices_geo <- get_vertices(grid_FINALFORCMS_geographic)
  grid_FINALFORCMS_geographic <- cbind(grid_FINALFORCMS_geographic, vertices_geo)
  
  # Verify columns
  cat("Geographic grid columns:", names(grid_FINALFORCMS_geographic), "\n\n")
  
  # Create CSV with centroids and their corresponding vertices
  centroid_coords <- geom(centroids_FINALFORCMS_geographic)[, c("x", "y")]
  
  centroids_vertices_data <- data.frame(
    unique_ID = centroids_FINALFORCMS_geographic$unique_ID,
    centroid_lon = centroid_coords[, 1],
    centroid_lat = centroid_coords[, 2],
    mean_coral_cover = grid_FINALFORCMS_geographic$mean_coral_cover,
    v1_lon = grid_FINALFORCMS_geographic$v1_x,
    v1_lat = grid_FINALFORCMS_geographic$v1_y,
    v2_lon = grid_FINALFORCMS_geographic$v2_x,
    v2_lat = grid_FINALFORCMS_geographic$v2_y,
    v3_lon = grid_FINALFORCMS_geographic$v3_x,
    v3_lat = grid_FINALFORCMS_geographic$v3_y,
    v4_lon = grid_FINALFORCMS_geographic$v4_x,
    v4_lat = grid_FINALFORCMS_geographic$v4_y
  )
  
  # Write CSV
  write.csv(centroids_vertices_data, here("output", "centroids_vertices_FINALFORCMS.csv"), row.names = FALSE)
  cat("Created CSV with", nrow(centroids_vertices_data), "centroids and vertices\n\n")
  
  # Read back in for verification
  verification_data <- read.csv(here("output", "centroids_vertices_FINALFORCMS.csv"))
  
  cat("=== VERIFICATION CHECKS ===\n\n")
  
  # Check 1: Row counts match
  cat("Original centroids:", nrow(centroids_FINALFORCMS_geographic), "\n")
  cat("CSV rows:", nrow(verification_data), "\n")
  cat("Counts match:", nrow(centroids_FINALFORCMS_geographic) == nrow(verification_data), "\n\n")
  
  # Check 2: IDs are identical and in same order
  cat("IDs identical:", identical(centroids_FINALFORCMS_geographic$unique_ID, verification_data$unique_ID), "\n")
  cat("IDs in same order:", all(centroids_FINALFORCMS_geographic$unique_ID == verification_data$unique_ID), "\n\n")
  
  # Check 3: Sample of values match
  cat("Sample comparison (first 5 rows):\n")
  print(head(verification_data, 5))
  cat("\n")
  
  # Visual verification with random sample
  set.seed(42)
  sample_n <- 8
  sample_idx <- sample(1:nrow(grid_FINALFORCMS_geographic), sample_n)
  
  # Get sample data from objects
  grid_sample <- grid_FINALFORCMS_geographic[sample_idx, ]
  cent_sample <- centroids_FINALFORCMS_geographic[sample_idx, ]
  
  # Get corresponding rows from CSV
  csv_sample <- verification_data[verification_data$unique_ID %in% grid_sample$unique_ID, ]
  
  cat("=== VISUAL VERIFICATION ===\n")
  cat("Plotting", sample_n, "random cells with vertices from both sources\n\n")
  
  # Plot with better visibility
  plot(grid_sample, col = NA, border = "blue", lwd = 3, main = "Verification: Grid, Centroids, and CSV Vertices")
  
  # Add vertices from CSV first (red dots at corners)
  for(i in 1:nrow(csv_sample)) {
    points(c(csv_sample$v1_lon[i], csv_sample$v2_lon[i], csv_sample$v3_lon[i], csv_sample$v4_lon[i]),
           c(csv_sample$v1_lat[i], csv_sample$v2_lat[i], csv_sample$v3_lat[i], csv_sample$v4_lat[i]),
           col = "red", pch = 16, cex = 2)
    # Label vertices with position number
    text(csv_sample$v1_lon[i], csv_sample$v1_lat[i], "1", pos = 2, cex = 0.6, col = "red")
    text(csv_sample$v2_lon[i], csv_sample$v2_lat[i], "2", pos = 4, cex = 0.6, col = "red")
    text(csv_sample$v3_lon[i], csv_sample$v3_lat[i], "3", pos = 4, cex = 0.6, col = "red")
    text(csv_sample$v4_lon[i], csv_sample$v4_lat[i], "4", pos = 2, cex = 0.6, col = "red")
  }
  
  # Add centroids on top (larger, different color)
  plot(cent_sample, col = "cyan", pch = 17, cex = 2.5, add = TRUE)
  plot(cent_sample, col = "darkgreen", pch = 17, cex = 2, add = TRUE)
  
  # Add labels showing both centroid and grid IDs
  cent_coords <- geom(cent_sample)[, c("x", "y")]
  grid_coords <- crds(centroids(grid_sample))
  
  for(i in 1:nrow(cent_sample)) {
    # Centroid ID above (green)
    text(cent_coords[i, "x"], cent_coords[i, "y"], 
         labels = paste0("C:", cent_sample$unique_ID[i]), 
         pos = 3, offset = 1.2, cex = 0.9, col = "darkgreen", font = 2)
    
    # Grid ID below (blue)
    text(grid_coords[i, "x"], grid_coords[i, "y"], 
         labels = paste0("G:", grid_sample$unique_ID[i]), 
         pos = 1, offset = 1.2, cex = 0.9, col = "blue", font = 2)
  }
  
  legend("topright", 
         legend = c("Grid borders", "Centroids (ID labeled)", "Vertices 1-4 (from CSV)"), 
         col = c("blue", "darkgreen", "red"), 
         lwd = c(3, NA, NA), 
         pch = c(NA, 17, 16),
         pt.cex = c(NA, 2, 2),
         cex = 0.9, bg = "white")
  
  cat("\n=== FINAL CHECK ===\n")
  cat("If vertices (red dots) align with grid corners (blue squares),\n")
  cat("and centroids (green) are inside their respective grids,\n")
  cat("then the CSV is correctly structured!\n")
  
  ################################## Save objects/workspace ##################################
  
  # #updated way to handle saving of new objects
  # save_new_objects("output/output_create_habitat_grid", existing_objects)
  