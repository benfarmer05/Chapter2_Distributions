  
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
  library(extrafont)
  library(cowplot)
  
  source(here("src/functions.R"))
  
  ################################## CONTROL TOGGLES ##################################
  
  # Toggle 1: Run mode - "single" for one species, "all" for all species
  RUN_MODE <- "all"  # Options: "single" or "all"
  
  # Toggle 2: Species selection (only used if RUN_MODE = "single")
  # SINGLE_SPECIES <- "dendrogyra"  # Options: see species list below
  # SINGLE_SPECIES <- "porites"  # Options: see species list below
  SINGLE_SPECIES <- "orbicella"  # Options: see species list below
  
  # Toggle 3: Raster constraint - TRUE for maxcell = Inf, FALSE for default downsampling
  UNCONSTRAINED_RASTER <- TRUE  # Options: TRUE or FALSE
  
  # Available species (for reference):
  # "total", "agaricia", "madracis", "porites", "siderastrea", "orbicella", 
  # "meandrina", "mycetophyllia", "pseudodiploria", "stephanocoenia", 
  # "montastraea", "dendrogyra", "dichocoenia", "diploria", "eusmilia", "colpophyllia"
  
  ################################## Setup ##################################
  
  load(here("output", "all_combined_data.rda"))
  spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
  bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
  terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
  
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
  
  results_files <- list.files(here('output/output_maps'), pattern = '^results_light_.*\\.rds$', full.names = TRUE)
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
  
  total_raster_latlon <- project(total_raster, "EPSG:4326")
  domain_extent <- ext(-67.99, -63.99, 17.63, 18.83)
  total_raster_cropped <- crop(total_raster_latlon, domain_extent)
  
  ################################## Load saved OSM data ##################################
  
  land_vect <- readRDS(here("output", "osm_land_vect.rds"))
  island_df <- readRDS(here("output", "osm_island_labels.rds"))
  
  ################################## Labels ##################################
  
  labels_stt_stj <- data.frame(
    name = c("Saint Thomas", "Saint John", "Culebra", "Vieques"),
    x = c(-64.93106, -64.75019, -65.28551, -65.34985),
    y = c(18.3485, 18.33596, 18.318, 18.13872)
  )
  labels_stx <- data.frame(name = c("Saint Croix"), x = c(-64.77806), y = c(17.73595))
  labels_whole <- data.frame(name = c("Puerto Rico"), x = c(-66.41416), y = c(18.23930))
  
  get_region_labels <- function(region_name) {
    switch(region_name,
           "STT_STJ" = labels_stt_stj,
           "STX" = labels_stx,
           "Full_Domain" = labels_whole,
           data.frame(name = character(), x = numeric(), y = numeric()))
  }
  
  ################################## Regions ##################################
  
  regions <- list(
    stt_stj = list(xlim = c(-65.39305, -64.65467), ylim = c(18.06478, 18.46411), name = "STT_STJ", title = "St. Thomas/St. John"),
    stx = list(xlim = c(-64.9474, -64.4094), ylim = c(17.6145, 17.8561), name = "STX", title = "St. Croix"),
    whole = list(xlim = c(-67.9895, -63.9886), ylim = c(17.6145, 18.8255), name = "Full_Domain", title = "Full Region")
  )
  
  ################################## Species ##################################
  
  # Define species name formatting
  get_species_display_name <- function(choice) {
    species_names <- list(
      "total" = "All coral",
      "agaricia" = "Agaricia spp.",
      "madracis" = "Madracis spp.",
      "porites" = "Porites spp.",
      "siderastrea" = "Siderastrea spp.",
      "orbicella" = "Orbicella spp.",
      "meandrina" = "Meandrina spp.",
      "mycetophyllia" = "Mycetophyllia spp.",
      "pseudodiploria" = "Pseudodiploria spp.",
      "stephanocoenia" = "Stephanocoenia spp.",
      "montastraea" = "Montastraea cavernosa",
      "dendrogyra" = "Dendrogyra cylindrus",
      "dichocoenia" = "Dichocoenia stokesii",
      "diploria" = "Diploria labyrinthiformis",
      "eusmilia" = "Eusmilia fastigiata",
      "colpophyllia" = "Colpophyllia natans"
    )
    
    return(species_names[[choice]])
  }
  
  # Helper function to determine if species name needs special formatting
  get_species_formatting <- function(species_name) {
    # Check if it ends with "spp."
    has_spp <- grepl("spp\\.$", species_name)
    
    if (has_spp) {
      # Split into genus and "spp."
      genus <- sub(" spp\\.$", "", species_name)
      return(list(
        type = "genus_spp",
        genus = genus,
        full_name = species_name
      ))
    } else if (species_name == "All coral") {
      # Not italic
      return(list(
        type = "non_italic",
        full_name = species_name
      ))
    } else {
      # Full species name (genus + species), all italic
      return(list(
        type = "full_italic",
        full_name = species_name
      ))
    }
  }
  
  # Get list of all available species
  all_species <- c("total", names(raster_list))
  
  ################################## Margin Settings ##################################
  
  margin_stt_t <- -20
  margin_stt_r <- -3
  margin_stt_b <- -30
  margin_stt_l <- -10
  
  margin_stx_t <- -20
  margin_stx_r <- 3
  margin_stx_b <- -30
  margin_stx_l <- 0
  
  margin_whole_t <- -30
  margin_whole_r <- 3
  margin_whole_b <- 0
  margin_whole_l <- 3
  
  bottom_height <- 1
  topleft_width <- 1.2
  
  legend_x <- 0.0
  legend_y <- -0.15
  
  ################################## Region Plot Function ##################################
  
  create_region_plot <- function(region, raster_data, land_data, max_value, show_legend = TRUE, 
                                 textsize = 2.5, titlesize = 10, use_maxcell_inf = TRUE) {
    region_labels <- get_region_labels(region$name)
    buffer <- 0.1
    region_extent <- ext(region$xlim[1]-buffer, region$xlim[2]+buffer, region$ylim[1]-buffer, region$ylim[2]+buffer)
    raster_region <- crop(raster_data, region_extent)
    
    # Set maxcell parameter based on toggle
    maxcell_param <- if(use_maxcell_inf) Inf else 5e5
    
    ggplot() +
      geom_spatraster(data = raster_region * 100, maxcell = maxcell_param) +
      scale_fill_viridis_c(name = "Coral cover (%)", option = "turbo", limits = c(0, max_value), na.value = "transparent") +
      geom_spatvector(data = land_data, fill = "gray90", color = "black") +
      geom_text(data = region_labels, aes(x = x, y = y, label = name), size = textsize, family = 'Georgia', color = 'black') +
      coord_sf(xlim = region$xlim, ylim = region$ylim, expand = FALSE) +
      theme_classic(base_family = "Georgia") +
      theme(
        legend.position = if(show_legend) "bottom" else "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 6),
        plot.margin = margin(-10, -10, -10, -10, 'pt'),
        panel.spacing = unit(0, 'pt')
      )
  }
  
  ################################## Main Processing Function ##################################
  
  process_species <- function(species_choice, use_maxcell_inf = TRUE) {
    
    cat(paste0("\n========================================\n"))
    cat(paste0("Processing: ", species_choice, "\n"))
    cat(paste0("Raster mode: ", if(use_maxcell_inf) "Unconstrained (maxcell = Inf)" else "Constrained (default)", "\n"))
    cat(paste0("========================================\n"))
    
    # Select raster and prepare data
    if (species_choice == "total") {
      selected_raster <- total_raster_latlon
      selected_raster_cropped <- total_raster_cropped
      species_name <- get_species_display_name(species_choice)
      species_formatting <- get_species_formatting(species_name)
    } else {
      selected_raster <- project(raster_list[[species_choice]], "EPSG:4326")
      selected_raster_cropped <- crop(selected_raster, domain_extent)
      species_name <- get_species_display_name(species_choice)
      species_formatting <- get_species_formatting(species_name)
    }
    
    max_coral_cover <- max(values(selected_raster_cropped), na.rm = TRUE) * 100
    cat(paste0("Species: ", species_name, " (max cover: ", round(max_coral_cover, 1), "%)\n"))
    
    # Create plots with individual margins
    plot_stt_stj <- create_region_plot(regions$stt_stj, selected_raster_cropped, land_vect, 
                                       max_coral_cover, FALSE, textsize = 2.0, use_maxcell_inf = use_maxcell_inf) + 
      theme(plot.margin = margin(margin_stt_t, margin_stt_r, margin_stt_b, margin_stt_l, 'pt'),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA))
    
    plot_stx <- create_region_plot(regions$stx, selected_raster_cropped, land_vect, 
                                   max_coral_cover, TRUE, textsize = 2.5, use_maxcell_inf = use_maxcell_inf) +
      scale_fill_viridis_c(name = "Cover (%)", 
                           option = "turbo", 
                           limits = c(0, max_coral_cover), 
                           breaks = c(0, max_coral_cover),
                           labels = c("0", round(max_coral_cover, 0)),
                           na.value = "transparent") +
      theme(plot.margin = margin(margin_stx_t, margin_stx_r, margin_stx_b, margin_stx_l, 'pt'),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = c(legend_x, legend_y),
            legend.justification = c("right", "bottom"),
            legend.direction = "horizontal",
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_rect(fill = "transparent", color = NA),
            legend.key = element_rect(fill = "transparent", color = NA),
            legend.key.size = unit(0.35, "cm"),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA))
    
    # Add species title with appropriate formatting
    if (species_formatting$type == "genus_spp") {
      # Italicized genus + regular "spp."
      # Calculate offset based on genus name length (approximate)
      genus_length <- nchar(species_formatting$genus)
      offset <- 0.014 * genus_length  # Adjusted multiplier for Georgia font
      
      plot_stx <- plot_stx +
        annotate("text", x = regions$stx$xlim[1] + 0.05, y = regions$stx$ylim[2] - 0.01,
                 label = species_formatting$genus, size = 4, fontface = "italic", family = "Georgia", hjust = 0) +
        annotate("text", x = regions$stx$xlim[1] + 0.05 + offset, y = regions$stx$ylim[2] - 0.01,
                 label = "spp.", size = 4, fontface = "plain", family = "Georgia", hjust = 0)
    } else if (species_formatting$type == "full_italic") {
      # Full species name in italics
      plot_stx <- plot_stx +
        annotate("text", x = regions$stx$xlim[1] + 0.05, y = regions$stx$ylim[2] - 0.01,
                 label = species_formatting$full_name, size = 4, fontface = "italic", family = "Georgia", hjust = 0)
    } else {
      # Non-italic (e.g., "All coral")
      plot_stx <- plot_stx +
        annotate("text", x = regions$stx$xlim[1] + 0.05, y = regions$stx$ylim[2] - 0.01,
                 label = species_formatting$full_name, size = 4, fontface = "bold", family = "Georgia", hjust = 0)
    }
    
    plot_whole <- create_region_plot(regions$whole, selected_raster_cropped, land_vect, 
                                     max_coral_cover, FALSE, textsize = 3.0, use_maxcell_inf = use_maxcell_inf) +
      theme(plot.margin = margin(margin_whole_t, margin_whole_r, margin_whole_b, margin_whole_l, 'pt'),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA))
    
    # Combine plots
    top_row <- plot_grid(plot_stt_stj, plot_stx, nrow = 1, align = 'h', rel_widths = c(topleft_width, 1))
    final_plot <- plot_grid(top_row, plot_whole, ncol = 1, align = 'v', rel_heights = c(1, bottom_height), axis = 'tb')
    
    # Prepare filename
    filename_base <- gsub("\n", "_", species_name)
    filename_base <- gsub(" ", "_", filename_base)
    
    # Add suffix for raster mode
    raster_suffix <- if(use_maxcell_inf) "_unconstrained" else "_constrained"
    
    # Save plots
    cat("Saving plots...\n")
    ggsave(filename = here("output", "output_figures_tables", paste0("fig_combined_", filename_base, raster_suffix, ".png")),
           plot = final_plot,
           width = 7,
           height = 5.5,
           dpi = 300,
           bg = "white",
           limitsize = FALSE)
    
    ggsave(filename = here("output", "output_figures_tables", paste0("fig_combined_", filename_base, raster_suffix, ".pdf")),
           plot = final_plot,
           width = 7,
           height = 5.5,
           device = cairo_pdf,
           limitsize = FALSE)
    
    cat(paste0("✓ Plots saved for ", species_name, "\n"))
    
    return(invisible(NULL))
  }
  
  ################################## Execute Based on Toggles ##################################
  
  cat("\n")
  cat("==========================================\n")
  cat("CORAL COVER MAP GENERATION\n")
  cat("==========================================\n")
  cat(paste0("Run mode: ", toupper(RUN_MODE), "\n"))
  cat(paste0("Raster constraint: ", if(UNCONSTRAINED_RASTER) "OFF (maxcell = Inf)" else "ON (default downsampling)", "\n"))
  
  if (RUN_MODE == "single") {
    cat(paste0("Selected species: ", SINGLE_SPECIES, "\n"))
    cat("==========================================\n")
    
    # Process single species
    process_species(SINGLE_SPECIES, use_maxcell_inf = UNCONSTRAINED_RASTER)
    
  } else if (RUN_MODE == "all") {
    cat(paste0("Processing all ", length(all_species), " species\n"))
    cat("==========================================\n")
    
    # Process all species
    for (species in all_species) {
      tryCatch({
        process_species(species, use_maxcell_inf = UNCONSTRAINED_RASTER)
      }, error = function(e) {
        cat(paste0("✗ ERROR processing ", species, ": ", e$message, "\n"))
      })
    }
    
    cat("\n==========================================\n")
    cat("ALL SPECIES PROCESSING COMPLETE\n")
    cat("==========================================\n")
    
  } else {
    stop("Invalid RUN_MODE. Must be 'single' or 'all'")
  }
  
  cat("\nScript completed successfully!\n")




  # # # can reference; but old
  # 
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
  # library(extrafont)
  # library(cowplot)
  # 
  # source(here("src/functions.R"))
  # 
  # ################################## Setup ##################################
  # 
  # load(here("output", "all_combined_data.rda"))
  # spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
  # bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
  # terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
  # 
  # load_spat_objects(directory = 'output/output_summarize_maps/')
  # source(here("src/functions.R"))
  # load(here("output/output_summarize_maps", "output_summarize_maps_workspace.Rdata"))
  # 
  # susc_rasters <- lapply(susc_rasters, function(x) {
  #   r <- terra::unwrap(x)
  #   terra::crs(r) <- spatial_metadata$bathy_final$crs
  #   return(r)
  # })
  # 
  # total_raster <- terra::unwrap(total_raster)
  # terra::crs(total_raster) <- spatial_metadata$bathy_final$crs
  # 
  # results_files <- list.files(here('output/output_maps'), pattern = '^results_light_.*\\.rds$', full.names = TRUE)
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
  # total_raster_latlon <- project(total_raster, "EPSG:4326")
  # domain_extent <- ext(-67.99, -63.99, 17.63, 18.83)
  # total_raster_cropped <- crop(total_raster_latlon, domain_extent)
  # 
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
  # ################################## Load saved OSM data ##################################
  # 
  # land_vect <- readRDS(here("output", "osm_land_vect.rds"))
  # island_df <- readRDS(here("output", "osm_island_labels.rds"))
  # 
  # ################################## Labels ##################################
  # 
  # labels_stt_stj <- data.frame(
  #   name = c("Saint Thomas", "Saint John", "Culebra", "Vieques", "BVI"),
  #   x = c(-64.93106, -64.75019, -65.28551, -65.34985, -64.65402),
  #   y = c(18.35037, 18.33783, 18.32049, 18.13872, 18.40760)
  # )
  # labels_stx <- data.frame(name = c("Saint Croix"), x = c(-64.77806), y = c(17.73595))
  # labels_whole <- data.frame(name = c("Puerto Rico"), x = c(-66.41416), y = c(18.23930))
  # 
  # get_region_labels <- function(region_name) {
  #   switch(region_name,
  #          "STT_STJ" = labels_stt_stj,
  #          "STX" = labels_stx,
  #          "Full_Domain" = labels_whole,
  #          data.frame(name = character(), x = numeric(), y = numeric()))
  # }
  # 
  # ################################## Regions ##################################
  # 
  # regions <- list(
  #   stt_stj = list(xlim = c(-65.39305, -64.65467), ylim = c(18.06478, 18.46411), name = "STT_STJ", title = "St. Thomas/St. John"),
  #   stx = list(xlim = c(-64.9474, -64.4094), ylim = c(17.6145, 17.8561), name = "STX", title = "St. Croix"),
  #   whole = list(xlim = c(-67.9895, -63.9886), ylim = c(17.6345, 18.8255), name = "Full_Domain", title = "Full Region")
  # )
  # 
  # ################################## Species ##################################
  # 
  # species_choice <- "total"  # Change this to your desired species
  # # species_choice <- "agaricia"
  # # species_choice <- "porites"
  # 
  # if (species_choice == "total") {
  #   selected_raster <- total_raster_latlon
  #   selected_raster_cropped <- total_raster_cropped
  #   species_name <- "Total Coral"
  # } else {
  #   selected_raster <- project(raster_list[[species_choice]], "EPSG:4326")
  #   selected_raster_cropped <- crop(selected_raster, domain_extent)
  #   species_name <- paste0(toupper(substring(species_choice, 1, 1)), substring(species_choice, 2))
  # }
  # 
  # max_coral_cover <- max(values(selected_raster_cropped), na.rm = TRUE) * 100
  # cat(paste0("Plotting species: ", species_name, " (max cover: ", round(max_coral_cover, 1), "%)\n"))
  # 
  # ################################## Region Plot Function ##################################
  # 
  # create_region_plot <- function(region, raster_data, land_data, max_value, show_legend = TRUE, 
  #                                textsize = 9, titlesize = 10) {
  #   region_labels <- get_region_labels(region$name)
  #   buffer <- 0.1
  #   region_extent <- ext(region$xlim[1]-buffer, region$xlim[2]+buffer, region$ylim[1]-buffer, region$ylim[2]+buffer)
  #   raster_region <- crop(raster_data, region_extent)
  #   
  #   ggplot() +
  #     geom_spatraster(data = raster_region * 100) +
  #     scale_fill_viridis_c(name = "Coral cover (%)", option = "turbo", limits = c(0, max_value), na.value = "transparent") +
  #     geom_spatvector(data = land_data, fill = "gray90", color = "black") +
  #     geom_text(data = region_labels, aes(x = x, y = y, label = name), size = 2.5, family = 'Georgia', color = 'black') +
  #     coord_sf(xlim = region$xlim, ylim = region$ylim, expand = FALSE) +
  #     theme_classic(base_family = "Georgia") +
  #     theme(
  #       legend.position = if(show_legend) "bottom" else "none",
  #       axis.title = element_blank(),
  #       axis.text = element_text(size = 6),
  #       plot.margin = margin(-10, -10, -10, -10, 'pt'),
  #       panel.spacing = unit(0, 'pt')
  #     )
  # }
  # 
  # ################################## Combine Plots (Tight Layout) ##################################
  # 
  # cat("Creating multi-panel plot...\n")
  # 
  # plot_stt_stj <- create_region_plot(regions$stt_stj, selected_raster_cropped, land_vect, max_coral_cover, FALSE)
  # plot_stx     <- create_region_plot(regions$stx, selected_raster_cropped, land_vect, max_coral_cover, FALSE)
  # plot_whole   <- create_region_plot(regions$whole, selected_raster_cropped, land_vect, max_coral_cover, TRUE)
  # 
  # top_row <- plot_grid(plot_stt_stj, plot_stx, nrow = 1, align = 'h', rel_widths = c(1, 1))
  # combined_plot <- plot_grid(top_row, plot_whole, ncol = 1, align = 'v', rel_heights = c(1, 1), axis = 'tb')
  # 
  # title <- ggdraw() + draw_label(species_name, fontface = 'bold', size = 14, hjust = 0.5, vjust = 0.9, fontfamily = 'Georgia')
  # final_plot <- plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.08, 1))
  # 
  # ################################## Save ##################################
  # 
  # # Display
  # print(final_plot)
  # 
  # # Save with proper species name in filename
  # ggsave(filename = here("output", paste0("fig_combined_", gsub(" ", "_", species_name), ".png")),
  #        plot = final_plot,
  #        width = 7,
  #        height = 5.5,
  #        dpi = 600,
  #        limitsize = FALSE)
  # 
  # ggsave(filename = here("output", paste0("fig_combined_", gsub(" ", "_", species_name), ".pdf")),
  #        plot = final_plot,
  #        width = 7,
  #        height = 5.5,
  #        device = cairo_pdf,
  #        limitsize = FALSE)
  # 
  # cat("Plots saved successfully!\n")
