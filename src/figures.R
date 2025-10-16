 
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(terra)
  library(leaflet)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  
  source(here("src/functions.R"))
  
  ################################## load spatial metadata ##################################
  
  #load 650 m habitat grid
  load_spat_objects(directory = 'output/output_create_habitat_grid/')
  source(here("src/functions.R"))
  load(here('output', 'output_create_habitat_grid/create_habitat_grid_workspace.RData'))
  
  
  
  load(here("output", "all_combined_data.rda"))
  
  cat("Loading spatial metadata...\n")
  
  # Load the spatial metadata
  spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
  
  # Load the raster
  bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
  
  # Apply the stored CRS
  terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
  
  cat("CRS loaded successfully\n\n")
  
  
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
  
  
  ################################## create total cover map ##################################
  
  cat("Creating total cover map...\n")
  
  # Clamp total cover
  total_raster_clamped <- clamp(total_raster, lower = 0, upper = 0.5, values = TRUE)
  
  # # Create custom color palette matching the paper
  # # Blue → Cyan → Yellow → Orange/Red
  # paper_colors <- colorRampPalette(c("#0000FF", "#00BFFF", "#00FFFF", 
  #                                    "#FFFF00", "#FF8C00", "#FF4500"))(256)
  # 
  # # Use same paper-style palette for both
  # unified_pal <- colorNumeric(paper_colors, domain = c(0, 0.5), na.color = "transparent")
  # 
  # total_cover_map <- leaflet() %>%
  #   addTiles(group = "OpenStreetMap") %>%
  #   addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  #   addRasterImage(total_raster_clamped, colors = unified_pal, opacity = 0.7) %>%
  #   addCircleMarkers(data = psu_total_cover, ~lon, ~lat, radius = 8,
  #                    fillColor = ~unified_pal(total_cover_clamped), fillOpacity = 0.8,
  #                    color = "white", weight = 2, stroke = TRUE,
  #                    popup = ~paste0("<b>PSU:</b> ", PSU, "<br><b>Total Cover:</b> ", 
  #                                    round(total_cover, 2), "%")) %>%
  #   addLegend("bottomright", pal = unified_pal, values = values(total_raster_clamped),
  #             title = "Total Coral Cover (%)",
  #             labFormat = labelFormat(transform = function(x) x * 100, suffix = "%")) %>%
  #   addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
  #                    options = layersControlOptions(collapsed = FALSE))
  # 
  # total_cover_map
  # 
  # cat("\n=== TOTAL COVER MAP COMPLETE ===\n")
  # cat("Species included:", length(all_results), "\n")
  # cat("PSUs with data:", nrow(psu_total_cover), "\n")
  
  
  
  
  
  #experimenting
  cat("Creating total cover map...\n")
  
  # Clamp total cover
  total_raster_clamped <- clamp(total_raster, lower = 0, upper = 0.5, values = TRUE)
  # total_raster_clamped <- aggregate(total_raster_clamped, fact=2, fun=mean)
  
  # Create YlOrRd color palette from RColorBrewer
  library(RColorBrewer)
  library(viridis)
  # ylOrRd_colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(256)
  colormap <- viridis(256)
  # colormap <- viridis(256, option = 'turbo')
  
  unified_pal <- colorNumeric(colormap, domain = c(0, 0.5), na.color = "transparent")
  
  total_cover_map <- leaflet() %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addTiles(group = "OpenStreetMap") %>%
    addRasterImage(total_raster_clamped, colors = unified_pal, opacity = 1.0) %>% #0.7
    addCircleMarkers(data = psu_total_cover, ~lon, ~lat, radius = 8,
                     fillColor = ~unified_pal(total_cover_clamped), fillOpacity = 1.0, #0.8
                     color = "white", weight = 2, stroke = TRUE,
                     popup = ~paste0("<b>PSU:</b> ", PSU, "<br><b>Total Cover:</b> ", 
                                     round(total_cover, 2), "%")) %>%
    addLegend("bottomright", pal = unified_pal, values = values(total_raster_clamped),
              title = "Total Coral Cover (%)",
              labFormat = labelFormat(transform = function(x) x * 100, suffix = "%")) %>%
    addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
                     options = layersControlOptions(collapsed = FALSE))
  
  total_cover_map
  
  cat("\n=== TOTAL COVER MAP COMPLETE ===\n")
  cat("Species included:", length(all_results), "\n")
  cat("PSUs with data:", nrow(psu_total_cover), "\n")
  
  
  
  ################################## prepare data for ggplot ##################################
  
  
  # Display all palettes
  display.brewer.all(colorblindFriendly = TRUE)
  
  # Convert raster to data frame
  total_cover_df <- as.data.frame(total_raster_clamped, xy = TRUE)
  colnames(total_cover_df) <- c("Longitude", "Latitude", "Cover")
  
  # Convert to percentage
  total_cover_df$Cover_pct <- total_cover_df$Cover * 100
  
  # Remove NA values
  total_cover_df <- total_cover_df[!is.na(total_cover_df$Cover_pct), ]
  
  
  ################################## create publication quality plot ##################################
  
  # Define plot extent options
  # plot_extents = ext(280000, 310000, 2010000, 2060000) #for investigating drops
  # plot_extents = ext(280000, 310000, 2000000, 2040000) #for investigating south of STT
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
  
  cat("\nCreating publication quality plot...\n")
  
  # Transform observation points to match raster CRS (UTM)
  obs_points <- vect(psu_total_cover, geom = c("lon", "lat"), crs = "EPSG:4326")
  obs_points_utm <- project(obs_points, crs(total_raster_clamped))
  
  # Extract coordinates in UTM and clamp cover values to 0-50%
  psu_total_cover_utm <- as.data.frame(geom(obs_points_utm)[, c("x", "y")])
  psu_total_cover_utm$total_cover <- pmin(psu_total_cover$total_cover, 50)  # Clamp to 50%
  psu_total_cover_utm$PSU <- psu_total_cover$PSU
  
  # Create the plot
  p <- ggplot() +
    # Raster layer
    geom_raster(data = total_cover_df, 
                aes(x = Longitude, y = Latitude, fill = Cover_pct)) +
    
    # Add observation points (colored by cover, fixed size)
    geom_point(data = psu_total_cover_utm, 
               aes(x = x, y = y, fill = total_cover),
               shape = 21, color = "white", size = 1.4,
               stroke = 0.3, alpha = 0.9) +
    
    # viridis turbo color scale for both raster and points
    scale_fill_viridis_c(option = "turbo",
                         name = "Coral Cover (%)",
                         limits = c(0, 50),
                         breaks = seq(0, 50, 10)) +
    
    # Labels
    labs(title = "Total Coral Cover Across All Species",
         subtitle = "Model predictions (colored raster) with field observations (circles)",
         x = "Easting (m)",
         y = "Northing (m)") +
    
    # # Coordinate system
    # coord_equal() +
    
    # Coordinate system
    coord_equal(xlim = c(xmin(plot_extents), xmax(plot_extents)),
                ylim = c(ymin(plot_extents), ymax(plot_extents))) +
    
    # Clean theme for publication
    theme_bw(base_size = 12) +
    theme(
      # Panel
      panel.grid.major = element_line(color = "gray90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      
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
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      legend.key = element_rect(fill = "white"),
      legend.spacing.y = unit(0.3, "cm")
    )
  
  # Display the plot
  print(p)
  
  cat("\n✓ Publication-quality plot created\n")  
  
  
  # # Save high-resolution version
  # ggsave(here("output/output_maps/total_coral_cover_publication.png"), 
  #        plot = p, 
  #        width = 10, 
  #        height = 8, 
  #        dpi = 300, 
  #        bg = "white")
  # 
  # ggsave(here("output/output_maps/total_coral_cover_publication.pdf"), 
  #        plot = p, 
  #        width = 10, 
  #        height = 8, 
  #        device = "pdf")
  # 
  # cat("\n✓ Publication-quality plots saved:\n")
  # cat("  - PNG: output/output_maps/total_coral_cover_publication.png\n")
  # cat("  - PDF: output/output_maps/total_coral_cover_publication.pdf\n")
  
  
  ################################## predicted vs observed scatterplot ##################################
  
  cat("\nCreating predicted vs observed scatterplot...\n")
  
  # Extract predicted values at observation locations
  # First, create spatial points in lat/lon (WGS84)
  obs_points <- vect(psu_total_cover, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  # Transform to match raster CRS (UTM)
  obs_points_utm <- project(obs_points, crs(total_raster_clamped))
  
  # Extract values
  predicted_values <- extract(total_raster_clamped, obs_points_utm)
  
  # Combine observed and predicted
  comparison_df <- data.frame(
    PSU = psu_total_cover$PSU,
    observed = psu_total_cover$total_cover,
    predicted = predicted_values[,2] * 100  # Convert to percentage
  )
  
  # Check for NAs (observations that fall outside raster or on NA cells)
  n_total <- nrow(comparison_df)
  n_na <- sum(is.na(comparison_df$predicted))
  
  cat("  Total observation sites:", n_total, "\n")
  cat("  Sites with NA predictions (outside raster or on masked cells):", n_na, "\n")
  
  # Remove any NAs
  comparison_df <- comparison_df[complete.cases(comparison_df), ]
  
  cat("  Valid prediction/observation pairs:", nrow(comparison_df), "\n")
  
  if(nrow(comparison_df) < n_total * 0.5) {
    warning("More than 50% of observations fell on NA raster cells. Check raster coverage.")
  }
  
  # Calculate statistics
  correlation <- cor(comparison_df$observed, comparison_df$predicted)
  rmse <- sqrt(mean((comparison_df$observed - comparison_df$predicted)^2))
  mae <- mean(abs(comparison_df$observed - comparison_df$predicted))
  
  # Create scatterplot
  scatter_plot <- ggplot(comparison_df, aes(x = observed, y = predicted)) +
    geom_point(size = 3, alpha = 0.6, color = "#E74C3C") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "#3498DB", fill = "#3498DB", alpha = 0.2) +
    annotate("text", x = 5, y = 95, 
             label = paste0("r = ", round(correlation, 3), "\n",
                            "RMSE = ", round(rmse, 2), "%\n",
                            "MAE = ", round(mae, 2), "%"),
             hjust = 0, size = 4, fontface = "bold") +
    labs(title = "Predicted vs Observed Total Coral Cover",
         subtitle = paste0("n = ", nrow(comparison_df), " observation sites"),
         x = "Observed Cover (%)",
         y = "Predicted Cover (%)") +
    coord_equal(xlim = c(0, 50), ylim = c(0, 50)) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "gray90", size = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(face = "bold", size = 11),
      axis.ticks = element_line(color = "black")
    )
  
  print(scatter_plot)
  
  # # Save scatterplot
  # ggsave(here("output/output_maps/predicted_vs_observed_scatterplot.png"), 
  #        plot = scatter_plot, 
  #        width = 8, 
  #        height = 8, 
  #        dpi = 300, 
  #        bg = "white")
  # 
  # ggsave(here("output/output_maps/predicted_vs_observed_scatterplot.pdf"), 
  #        plot = scatter_plot, 
  #        width = 8, 
  #        height = 8, 
  #        device = "pdf")
  
  cat("\n✓ Scatterplot saved:\n")
  cat("  - PNG: output/output_maps/predicted_vs_observed_scatterplot.png\n")
  cat("  - PDF: output/output_maps/predicted_vs_observed_scatterplot.pdf\n")
  cat("\nModel Performance:\n")
  cat("  Correlation (r):", round(correlation, 3), "\n")
  cat("  RMSE:", round(rmse, 2), "%\n")
  cat("  MAE:", round(mae, 2), "%\n")
  
  
  # ################################## SD scatterplot ##################################
  # 
  # # Fit linear model and calculate SD
  # lm_fit <- lm(predicted ~ observed, data = comparison_df)
  # residual_sd <- sd(residuals(lm_fit))
  # 
  # # Create prediction data for smooth line with SD bands
  # pred_data <- data.frame(observed = seq(0, 50, length.out = 100))
  # pred_data$predicted <- predict(lm_fit, newdata = pred_data)
  # pred_data$lower_1sd <- pred_data$predicted - residual_sd
  # pred_data$upper_1sd <- pred_data$predicted + residual_sd
  # pred_data$lower_2sd <- pred_data$predicted - 2*residual_sd
  # pred_data$upper_2sd <- pred_data$predicted + 2*residual_sd
  # 
  # # Create scatterplot
  # scatter_plot <- ggplot(comparison_df, aes(x = observed, y = predicted)) +
  #   # # Add 2SD band (lighter)
  #   # geom_ribbon(data = pred_data, aes(y = predicted, ymin = lower_2sd, ymax = upper_2sd),
  #   #             fill = "#3498DB", alpha = 0.1) +
  #   # Add 1SD band (darker)
  #   geom_ribbon(data = pred_data, aes(y = predicted, ymin = lower_1sd, ymax = upper_1sd),
  #               fill = "#3498DB", alpha = 0.2) +
  #   # Regression line
  #   geom_line(data = pred_data, aes(x = observed, y = predicted), 
  #             color = "#3498DB", size = 1) +
  #   # Data points
  #   geom_point(size = 3, alpha = 0.6, color = "#E74C3C") +
  #   # 1:1 line
  #   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  #   # Stats annotation
  #   annotate("text", x = 5, y = 47, 
  #            label = paste0("r = ", round(correlation, 3), "\n",
  #                           "RMSE = ", round(rmse, 2), "%\n",
  #                           "MAE = ", round(mae, 2), "%\n",
  #                           "SD = ", round(residual_sd, 2), "%"),
  #            hjust = 0, size = 4, fontface = "bold") +
  #   labs(title = "Predicted vs Observed Total Coral Cover",
  #        subtitle = paste0("n = ", nrow(comparison_df), " observation sites | Shaded: ±1 SD (dark), ±2 SD (light)"),
  #        x = "Observed Cover (%)",
  #        y = "Predicted Cover (%)") +
  #   coord_equal(xlim = c(0, 50), ylim = c(0, 50)) +
  #   theme_bw(base_size = 12) +
  #   theme(
  #     panel.grid.major = element_line(color = "gray90", size = 0.3),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_rect(fill = "white"),
  #     panel.border = element_rect(color = "black", fill = NA, size = 1),
  #     plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
  #     plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
  #     axis.text = element_text(color = "black", size = 10),
  #     axis.title = element_text(face = "bold", size = 11),
  #     axis.ticks = element_line(color = "black")
  #   )
  # 
  # 
  # ################################## >0 cover obs scatter ##################################
  # 
  # 
  # # Filter to only observations with cover > 0
  # comparison_df_presence <- comparison_df %>%
  #   filter(observed > 0)
  # 
  # cat("  Observations with cover > 0:", nrow(comparison_df_presence), 
  #     "out of", nrow(comparison_df), "total\n")
  # 
  # # Recalculate statistics for presence-only data
  # correlation_presence <- cor(comparison_df_presence$observed, comparison_df_presence$predicted)
  # rmse_presence <- sqrt(mean((comparison_df_presence$observed - comparison_df_presence$predicted)^2))
  # mae_presence <- mean(abs(comparison_df_presence$observed - comparison_df_presence$predicted))
  # 
  # # Fit linear model and calculate SD
  # lm_fit_presence <- lm(predicted ~ observed, data = comparison_df_presence)
  # residual_sd_presence <- sd(residuals(lm_fit_presence))
  # 
  # # Create prediction data for smooth line with SD bands
  # pred_data_presence <- data.frame(observed = seq(0, 50, length.out = 100))
  # pred_data_presence$predicted <- predict(lm_fit_presence, newdata = pred_data_presence)
  # pred_data_presence$lower_1sd <- pred_data_presence$predicted - residual_sd_presence
  # pred_data_presence$upper_1sd <- pred_data_presence$predicted + residual_sd_presence
  # pred_data_presence$lower_2sd <- pred_data_presence$predicted - 2*residual_sd_presence
  # pred_data_presence$upper_2sd <- pred_data_presence$predicted + 2*residual_sd_presence
  # 
  # # Create scatterplot (presence only)
  # scatter_plot_presence <- ggplot(comparison_df_presence, aes(x = observed, y = predicted)) +
  #   # Add 2SD band (lighter)
  #   geom_ribbon(data = pred_data_presence, aes(y = predicted, ymin = lower_2sd, ymax = upper_2sd),
  #               fill = "#3498DB", alpha = 0.1) +
  #   # Add 1SD band (darker)
  #   geom_ribbon(data = pred_data_presence, aes(y = predicted, ymin = lower_1sd, ymax = upper_1sd),
  #               fill = "#3498DB", alpha = 0.2) +
  #   # Regression line
  #   geom_line(data = pred_data_presence, aes(x = observed, y = predicted), 
  #             color = "#3498DB", size = 1) +
  #   # Data points
  #   geom_point(size = 3, alpha = 0.6, color = "#E74C3C") +
  #   # 1:1 line
  #   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
  #   # Stats annotation
  #   annotate("text", x = 5, y = 47, 
  #            label = paste0("r = ", round(correlation_presence, 3), "\n",
  #                           "RMSE = ", round(rmse_presence, 2), "%\n",
  #                           "MAE = ", round(mae_presence, 2), "%\n",
  #                           "SD = ", round(residual_sd_presence, 2), "%"),
  #            hjust = 0, size = 4, fontface = "bold") +
  #   labs(title = "Predicted vs Observed Total Coral Cover (Presence Only)",
  #        subtitle = paste0("n = ", nrow(comparison_df_presence), " sites with observed cover > 0% | Shaded: ±1 SD (dark), ±2 SD (light)"),
  #        x = "Observed Cover (%)",
  #        y = "Predicted Cover (%)") +
  #   coord_equal(xlim = c(0, 50), ylim = c(0, 50)) +
  #   theme_bw(base_size = 12) +
  #   theme(
  #     panel.grid.major = element_line(color = "gray90", size = 0.3),
  #     panel.grid.minor = element_blank(),
  #     panel.background = element_rect(fill = "white"),
  #     panel.border = element_rect(color = "black", fill = NA, size = 1),
  #     plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
  #     plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
  #     axis.text = element_text(color = "black", size = 10),
  #     axis.title = element_text(face = "bold", size = 11),
  #     axis.ticks = element_line(color = "black")
  #   )
  # 
  # print(scatter_plot_presence)
  # 
  # # Save presence-only scatterplot
  # ggsave(here("output/output_maps/predicted_vs_observed_presence_only.png"), 
  #        plot = scatter_plot_presence, 
  #        width = 8, 
  #        height = 8, 
  #        dpi = 300, 
  #        bg = "white")
  # 
  # ggsave(here("output/output_maps/predicted_vs_observed_presence_only.pdf"), 
  #        plot = scatter_plot_presence, 
  #        width = 8, 
  #        height = 8, 
  #        device = "pdf")
  # 
  # cat("\n✓ Presence-only scatterplot saved\n")
  # cat("\nPresence-Only Model Performance:\n")
  # cat("  Correlation (r):", round(correlation_presence, 3), "\n")
  # cat("  RMSE:", round(rmse_presence, 2), "%\n")
  # cat("  MAE:", round(mae_presence, 2), "%\n")
  # cat("  SD:", round(residual_sd_presence, 2), "%\n")
  # 
  # 
