  
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
  
  # Toggle: Raster constraint - TRUE for maxcell = Inf, FALSE for default downsampling
  UNCONSTRAINED_RASTER <- TRUE  # Options: TRUE or FALSE
  
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
  
  ################################## Create Richness Raster ##################################
  
  cat("Creating species richness raster...\n")
  
  # Create richness raster (count of species with cover > 0 at each pixel)
  richness_raster <- Reduce(`+`, lapply(raster_list, function(r) r > 0))
  
  # Restore CRS
  terra::crs(richness_raster) <- spatial_metadata$bathy_final$crs
  
  # Project to lat/lon
  richness_raster_latlon <- project(richness_raster, "EPSG:4326")
  
  # Crop to domain
  domain_extent <- ext(-67.99, -63.99, 17.63, 18.83)
  richness_raster_cropped <- crop(richness_raster_latlon, domain_extent)
  
  cat("Richness raster created\n")
  
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
      geom_spatraster(data = raster_region, maxcell = maxcell_param) +
      scale_fill_viridis_c(name = "Species richness", option = "viridis", limits = c(0, max_value), na.value = "transparent") +
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
  
  ################################## Main Processing ##################################
  
  cat("\n========================================\n")
  cat("Processing: Species Richness\n")
  cat(paste0("Raster mode: ", if(UNCONSTRAINED_RASTER) "Unconstrained (maxcell = Inf)" else "Constrained (default)", "\n"))
  cat("========================================\n")
  
  max_richness <- max(values(richness_raster_cropped), na.rm = TRUE)
  cat(paste0("Max species richness: ", round(max_richness, 0), " species\n"))
  
  # Create plots with individual margins
  plot_stt_stj <- create_region_plot(regions$stt_stj, richness_raster_cropped, land_vect, 
                                     max_richness, FALSE, textsize = 2.0, use_maxcell_inf = UNCONSTRAINED_RASTER) + 
    theme(plot.margin = margin(margin_stt_t, margin_stt_r, margin_stt_b, margin_stt_l, 'pt'),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
  
  plot_stx <- create_region_plot(regions$stx, richness_raster_cropped, land_vect, 
                                 max_richness, TRUE, textsize = 2.5, use_maxcell_inf = UNCONSTRAINED_RASTER) +
    scale_fill_viridis_c(name = "Richness", 
                         option = "viridis", 
                         limits = c(0, max_richness), 
                         breaks = c(0, max_richness),
                         labels = c("0", round(max_richness, 0)),
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
          panel.background = element_rect(fill = "white", color = NA)) +
    annotate("text", x = regions$stx$xlim[1] + 0.05, y = regions$stx$ylim[2] - 0.01,
             label = "Species richness", size = 4, fontface = "bold", family = "Georgia", hjust = 0)
  
  plot_whole <- create_region_plot(regions$whole, richness_raster_cropped, land_vect, 
                                   max_richness, FALSE, textsize = 3.0, use_maxcell_inf = UNCONSTRAINED_RASTER) +
    theme(plot.margin = margin(margin_whole_t, margin_whole_r, margin_whole_b, margin_whole_l, 'pt'),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
  
  # Combine plots
  top_row <- plot_grid(plot_stt_stj, plot_stx, nrow = 1, align = 'h', rel_widths = c(topleft_width, 1))
  final_plot <- plot_grid(top_row, plot_whole, ncol = 1, align = 'v', rel_heights = c(1, bottom_height), axis = 'tb')
  
  # Add suffix for raster mode
  raster_suffix <- if(UNCONSTRAINED_RASTER) "_unconstrained" else "_constrained"
  
  # Save plots
  cat("Saving plots...\n")
  ggsave(filename = here("output", "output_figures_tables", paste0("fig_combined_species_richness", raster_suffix, ".png")),
         plot = final_plot,
         width = 7.087,
         height = 5.5,
         dpi = 300,
         bg = "white",
         limitsize = FALSE)
  
  # ggsave(filename = here("output", "output_figures_tables", paste0("fig_combined_species_richness", raster_suffix, ".pdf")),
  #        plot = final_plot,
  #        width = 7.087,
  #        height = 5.5,
  #        device = cairo_pdf,
  #        limitsize = FALSE)
  
  cat("✓ Species richness map saved\n")
  cat("\nScript completed successfully!\n")