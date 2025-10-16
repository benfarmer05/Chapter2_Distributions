    
  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(terra)
  # library(leaflet)
  # library(dplyr)
  # library(ggplot2)
  # library(RColorBrewer)
  
  source(here("src/functions.R"))
  
  
  ################################## load spatial metadata ##################################
  
  #load 650 m habitat grid initially created from 0-60 m depth binary
  load_spat_objects(directory = 'output/output_create_habitat_grid/')
  source(here("src/functions.R"))
  load(here('output', 'output_create_habitat_grid/create_habitat_grid_workspace.RData'))
  
  #load bathy
  spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
  bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
  terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
  
  ################################## read in all model outputs ##################################
  
  cat("Reading model outputs...\n")
  
  # Get all results files
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
  
  cat("Loaded", length(all_results), "species models\n\n")  
  
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
  
  
  ################################## sum observed cover by PSU ##################################
  
  cat("Calculating observed total cover...\n")
  
  # Get all species names from your results
  all_species <- names(all_results)
  
  # Load benthic data (assuming it's in the workspace or needs to be loaded)
  load(here("output", "all_combined_data.rda"))
  
  # Sum observed cover across all species
  # Extract genus from species name and match with model species
  psu_total_cover <- combined_benthic_data_averaged %>%
    mutate(genus = tolower(sub(" .*", "", spp))) %>%
    filter(genus %in% all_species) %>%
    group_by(PSU) %>%
    summarise(lat = first(lat), 
              lon = first(lon), 
              total_cover = sum(cover, na.rm = TRUE), 
              .groups = 'drop') %>%
    mutate(total_cover_clamped = pmin(total_cover / 100, 1))
  
  cat("PSU total cover calculated\n")
  cat("Species matched:", length(unique(combined_benthic_data_averaged %>% 
                                          mutate(genus = tolower(sub(" .*", "", spp))) %>% 
                                          filter(genus %in% all_species) %>% 
                                          pull(genus))), "\n\n")
  
  
  
  
  ################################## filter 650m grid to areas with predicted coral cover ##################################
  
  cat("Filtering 650m grid to predicted coral presence areas...\n")
  
  # Load the 650m grid
  grid_reefy_no_land <- vect(here("output", "grid_reefy_no_land.shp"))  # Adjust path as needed
  
  # Create binary mask where total cover > 0
  coral_presence_mask <- total_raster_clamped > 0
  
  # Convert mask to polygons
  coral_presence_poly <- as.polygons(coral_presence_mask, dissolve = TRUE)
  
  # Filter grid to only cells that intersect with coral presence areas
  has_coral <- relate(grid_reefy_no_land, coral_presence_poly, "intersects")
  grid_with_coral <- grid_reefy_no_land[has_coral]
  
  cat("  Original grid cells:", nrow(grid_reefy_no_land), "\n")
  cat("  Grid cells with predicted coral:", nrow(grid_with_coral), "\n")
  cat("  Reduction:", round((1 - nrow(grid_with_coral)/nrow(grid_reefy_no_land)) * 100, 1), "%\n\n")
  
  # Quick visualization
  plot(grid_reefy_no_land, col = "lightgray", border = "gray", main = "650m Grid Filtered to Coral Areas")
  plot(grid_with_coral, col = "lightgreen", border = "darkgreen", add = TRUE)
  
  # Save filtered grid if desired
  # writeVector(grid_with_coral, here("output", "grid_with_coral.shp"), overwrite = TRUE)
  
  
  
  
  ################################## filter 650m grid to predicted coral presence areas ##################################
  
  cat("Filtering 650m grid to areas with predicted coral cover...\n")
  
  # Extract mean coral cover for each grid cell
  cat("  Calculating mean coral cover per grid cell...\n")
  grid_coral_values <- extract(total_raster_clamped, grid_reefy_no_land, fun = mean, na.rm = TRUE)
  
  # Add mean coral cover to the grid
  grid_reefy_no_land$mean_coral_cover <- grid_coral_values[, 2]
  
  # Filter to cells with mean cover > 0
  has_coral <- grid_reefy_no_land$mean_coral_cover > 0 & !is.na(grid_reefy_no_land$mean_coral_cover)
  grid_with_coral <- grid_reefy_no_land[has_coral]
  
  # Filter centroids to match
  centroids_with_coral <- centroids_all[centroids_all$unique_ID %in% grid_with_coral$unique_ID]
  
  # Verify matching
  cat("  Original grid cells:", nrow(grid_reefy_no_land), "\n")
  cat("  Grid cells with predicted coral:", nrow(grid_with_coral), "\n")
  cat("  Matching centroids:", nrow(centroids_with_coral), "\n")
  cat("  Reduction:", round((1 - nrow(grid_with_coral)/nrow(grid_reefy_no_land)) * 100, 1), "%\n")
  cat("  Mean coral cover range:", round(range(grid_with_coral$mean_coral_cover, na.rm=TRUE), 3), "\n")
  cat("  IDs match:", all(grid_with_coral$unique_ID == centroids_with_coral$unique_ID), "\n\n")
  cat("  IDs match:", all(grid_reefy_no_land$unique_ID == centroids_all$unique_ID), "\n\n")
  
  # Visualization with coral cover gradient
  library(sf)
  grid_with_coral_sf <- st_as_sf(grid_with_coral)
  
  ggplot() +
    geom_sf(data = st_as_sf(grid_reefy_no_land), fill = "lightgray", color = "gray80", size = 0.2) +
    geom_sf(data = grid_with_coral_sf, aes(fill = mean_coral_cover * 100), color = "black", size = 0.1) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Mean Coral\nCover (%)") +
    geom_sf(data = st_as_sf(centroids_with_coral), color = "blue", size = 0.3, alpha = 0.5) +
    theme_minimal() +
    labs(title = "650m Grid: Mean Coral Cover per Cell",
         subtitle = paste0(nrow(grid_with_coral), " cells with predicted coral"))
  
  
  
  
  # Optional: Save filtered outputs with coral cover attribute
  # writeVector(grid_with_coral, here("output", "grid_with_coral.shp"), overwrite = TRUE)
  # writeVector(centroids_with_coral, here("output", "centroids_with_coral.shp"), overwrite = TRUE)