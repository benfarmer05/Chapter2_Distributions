# Species Leaflet Viewer (Cached)
# Caches processed data for fast loading, builds maps fresh per session

rm(list = ls())

library(here)
library(leaflet)
library(terra)
library(dplyr)

################################## TOGGLES ##################################

MODE <- "total"       # "single", "all", or "total"
SPECIES <- "orbicella" # Species for single mode
SHOW_OBS <- FALSE       # Show observation points?
REBUILD_CACHE <- FALSE # Set TRUE to force rebuild cache

################################## CACHE SETUP ##################################

cache_dir <- here("temp")
cache_file <- here("temp", "species_leaflet_cache.rds")

if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

cat("Cache location:", cache_dir, "\n\n")

################################## BUILD OR LOAD CACHE ##################################

if (file.exists(cache_file) && !REBUILD_CACHE) {
  
  cat("Loading from cache...\n")
  cache <- readRDS(cache_file)
  available_species <- names(cache)[names(cache) != "total"]
  cat("Cached species:", paste(available_species, collapse = ", "), "\n")
  cat("Total raster:", ifelse("total" %in% names(cache), "yes", "no"), "\n\n")
  
} else {
  
  cat("Building cache (one-time per session)...\n\n")
  
  # Load source data
  load(here("output", "all_combined_data.rda"))
  
  # Load spatial metadata for CRS
  spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
  
  results_files <- list.files(
    here('output/output_maps'), 
    pattern = '^results_light_.*\\.rds$', 
    full.names = TRUE
  )
  available_species <- gsub('results_light_|\\.rds', '', basename(results_files))
  
  cache <- list()
  
  # Process individual species
  for (sp in available_species) {
    cat("  Processing", sp, "...")
    
    results <- readRDS(here('output/output_maps', paste0('results_light_', sp, '.rds')))
    
    obs <- combined_benthic_data_averaged %>%
      filter(grepl(sp, spp, ignore.case = TRUE)) %>%
      group_by(PSU) %>%
      summarise(lat = first(lat), lon = first(lon), 
                avg_cover = mean(cover, na.rm = TRUE), .groups = 'drop')
    
    cache[[sp]] <- list(
      raster = results$raster,
      raster_raw = results$raster_raw,
      observations = obs,
      threshold = results$threshold,
      max_cover = results$max_cover
    )
    
    cat(" done (", nrow(obs), " obs)\n", sep = "")
  }
  
  # Process total (stacked) raster
  cat("  Processing total (stacked species)...")
  
  # Load and unwrap first raster to initialize
  first_rast <- terra::unwrap(cache[[available_species[1]]]$raster_raw)
  total_rast <- first_rast
  
  # Stack all species
  for (sp in available_species[-1]) {
    sp_rast <- terra::unwrap(cache[[sp]]$raster_raw)
    sp_rast <- resample(sp_rast, total_rast)  # Ensure alignment
    total_rast <- total_rast + sp_rast
  }
  
  # Get total observations (all coral)
  total_obs <- combined_benthic_data_averaged %>%
    group_by(PSU) %>%
    summarise(lat = first(lat), lon = first(lon),
              avg_cover = sum(cover, na.rm = TRUE), .groups = 'drop')
  
  max_total <- max(values(total_rast), na.rm = TRUE)
  
  cache[["total"]] <- list(
    raster = terra::wrap(total_rast),
    raster_raw = terra::wrap(total_rast),
    observations = total_obs,
    threshold = NA,
    max_cover = max_total
  )
  
  cat(" done (max cover:", round(max_total * 100, 1), "%)\n", sep = "")
  
  saveRDS(cache, cache_file)
  cat("\nCache saved to:", cache_file, "\n\n")
}

################################## MAP BUILDER ##################################

build_map <- function(species_name, show_observations = TRUE) {
  
  if (!species_name %in% names(cache)) {
    cat("ERROR: Species not in cache:", species_name, "\n")
    return(NULL)
  }
  
  sp_data <- cache[[species_name]]
  
  cat("Building map for", species_name, "...\n")
  if (!is.na(sp_data$threshold)) cat("  Threshold:", sp_data$threshold, "\n")
  cat("  Max cover:", round(sp_data$max_cover * 100, 2), "%\n")
  cat("  Observations:", nrow(sp_data$observations), ifelse(show_observations, "(showing)", "(hidden)"), "\n")
  
  # Unwrap raster
  raster_data <- terra::unwrap(sp_data$raster)
  
  # Prep observations
  obs <- sp_data$observations %>%
    mutate(avg_cover_clamped = pmin(avg_cover / 100, sp_data$max_cover))
  
  # Color palette
  paper_colors <- colorRampPalette(c("#0000FF", "#00BFFF", "#00FFFF", 
                                     "#FFFF00", "#FF8C00", "#FF4500"))(256)
  pal <- colorNumeric(paper_colors, domain = c(0, sp_data$max_cover), na.color = "transparent")
  
  # Title
  title <- ifelse(species_name == "total", 
                  "Total Coral (%)", 
                  paste0(tools::toTitleCase(species_name), " (%)"))
  
  # Build map
  map <- leaflet() %>%
    addTiles(group = "OpenStreetMap") %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addRasterImage(raster_data, colors = pal, opacity = 0.7) %>%
    addLegend("bottomright", pal = pal, values = c(0, sp_data$max_cover),
              title = title,
              labFormat = labelFormat(transform = function(x) round(x * 100, 1), suffix = "%")) %>%
    addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                     options = layersControlOptions(collapsed = FALSE))
  
  # Add observations if toggled on
  if (show_observations && nrow(obs) > 0) {
    map <- map %>%
      addCircleMarkers(data = obs, ~lon, ~lat, radius = 8,
                       fillColor = ~pal(avg_cover_clamped), fillOpacity = 0.8,
                       color = "white", weight = 2, stroke = TRUE,
                       popup = ~paste0("<b>PSU:</b> ", PSU, 
                                       "<br><b>Cover:</b> ", round(avg_cover, 2), "%"))
  }
  
  return(map)
}

################################## RUN ##################################

if (MODE == "single") {
  
  if (!SPECIES %in% names(cache)) {
    stop(paste("Species not found:", SPECIES, 
               "\nAvailable:", paste(names(cache), collapse = ", ")))
  }
  
  cat("=== SINGLE SPECIES MODE ===\n")
  print(build_map(SPECIES, show_observations = SHOW_OBS))
  
} else if (MODE == "total") {
  
  cat("=== TOTAL CORAL MODE ===\n")
  print(build_map("total", show_observations = SHOW_OBS))
  
} else if (MODE == "all") {
  
  cat("=== ALL SPECIES MODE ===\n")
  cat("Press ENTER to advance, 'q' to quit\n\n")
  
  for (sp in names(cache)) {
    cat("\n--- ", toupper(sp), " ---\n")
    print(build_map(sp, show_observations = SHOW_OBS))
    
    user_input <- readline(prompt = paste0("[", sp, "] ENTER=next, q=quit: "))
    if (tolower(user_input) == "q") break
  }
  
  cat("\n=== FINISHED ===\n")
}

################################## HELPER INFO ##################################

cat("\n---\nCache info:\n")
cat("  Location:", cache_file, "\n")
cat("  Species:", length(cache), "\n")
cat("  To rebuild: set REBUILD_CACHE <- TRUE\n")