  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(terra)
  library(ggplot2)
  library(tidyterra)
  library(cmocean)
  library(patchwork)
  library(viridis)
  library(Cairo)
  library(extrafont)
  
  source(here("src/functions.R"))
  
  ################################## CONTROL TOGGLE ##################################
  
  # Create output directory if it doesn't exist
  output_dir <- here("output", "output_figures_tables")
  
  # Raster constraint - TRUE for maxcell = Inf, FALSE for default downsampling
  UNCONSTRAINED_RASTER <- TRUE  # Options: TRUE or FALSE
  
  # Set maxcell parameter based on toggle
  maxcell_param <- if(UNCONSTRAINED_RASTER) 5e6 else 5e5
  
  ################################## LEGEND SETTINGS ##################################
  
  # Colorbar dimensions
  colorbar_width <- 3.5
  colorbar_height <- 0.27
  
  # Colorbar spacing
  label_vjust <- 4.5      # Gap between colorbar and numbers
  title_vjust <- 1        # Gap between colorbar and title
  ticks_linewidth <- 0.5
  
  # Legend position and background
  legend_x <- 0.05        # Horizontal position (0.05 = far left)
  legend_y <- 0.001        # Vertical position (0.01 = very low)
  legend_bg_alpha <- 0    # Transparency of legend background (0 = fully transparent)
  
  # Legend text sizes
  legend_title_size <- 7
  legend_text_size <- 6
  
  # Legend margins
  legend_margin <- 2      # Margin around legend in points
  
  # Panel label settings
  panel_label_size <- 5
  panel_label_hjust <- -0.5
  panel_label_vjust <- 1.5
  panel_label_fontface <- "bold"
  
  # Land vector settings
  land_linewidth <- 0.05   # Line width for land vector borders
  
  ################################## Setup ##################################
  
  #load bathy_final
  # Load the spatial metadata
  spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
  #
  # Load the raster
  bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
  #
  # Apply the stored CRS
  terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs

  #load derived bathy rasters, and oceanographic rasters
  load_spat_objects(directory = 'output/output_calculate_bathy_rasters/')
  load_spat_objects(directory = 'output/output_calculate_ocean_rasters/')

  # Load land vector data
  land_vect <- readRDS(here("output", "osm_land_vect.rds"))

  # Project land vector to match raster CRS
  land_vect <- project(land_vect, crs(bathy_final))

  # Crop land vector to match raster extent
  raster_extent <- ext(bathy_final)
  land_vect <- crop(land_vect, raster_extent)
  
  ################################## Mask Oceanographic Rasters ##################################
  
  # Create a mask where bathy is deeper than -60m (set these areas to NA)
  bathy_mask <- bathy_final > -60
  
  # Apply mask to all oceanographic rasters
  mean_sst_raster <- mask(mean_sst_raster, bathy_mask, maskvalues = FALSE)
  dir_erddap_raster <- mask(dir_erddap_raster, bathy_mask, maskvalues = FALSE)
  mean_spm_raster <- mask(mean_spm_raster, bathy_mask, maskvalues = FALSE)
  TPI_terra <- mask(TPI_terra, bathy_mask, maskvalues = FALSE)
  mean_par_raster <- mask(mean_par_raster, bathy_mask, maskvalues = FALSE)
  
  ################################## Font Setup ##################################
  
  # Load fonts for use in plots
  extrafont::loadfonts(device = "win", quiet = TRUE)
  
  ################################## Calculate UTM Coordinate Limits ##################################
  
  # Convert geographic bounds to UTM
  geo_coords <- data.frame(
    lon = c(-67.97, -64.13192),
    lat = c(17.45, 18.84)  # Original min lat: 17.55 (reduced to create space for legend)
  )
  
  # Create a vector object in WGS84
  geo_vect <- vect(geo_coords, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  # Transform to UTM zone 20N (match bathy_final CRS)
  utm_vect <- project(geo_vect, crs(bathy_final))
  
  # Extract the coordinates
  utm_coords <- crds(utm_vect)
  
  # Get the bounding box
  xmin_utm <- min(utm_coords[, "x"])
  xmax_utm <- max(utm_coords[, "x"])
  ymin_utm <- min(utm_coords[, "y"])
  ymax_utm <- max(utm_coords[, "y"])
  
  ################################## Create Individual Plots ##################################
  
  # Common theme for all panels - legends inside, no titles, no axis elements
  common_theme <- theme_classic(base_family = "Georgia") +
    theme(
      legend.position = c(legend_x, legend_y),
      legend.justification = c(0, 0),  # Anchor at bottom-left of legend
      legend.direction = "horizontal",
      legend.key.width = unit(0.8, "cm"),
      legend.key.height = unit(0.2, "cm"),
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size),
      legend.title.position = "left",
      legend.background = element_rect(fill = alpha("white", legend_bg_alpha), color = NA),
      legend.box.background = element_rect(fill = alpha("white", legend_bg_alpha), color = NA),
      legend.margin = margin(legend_margin, legend_margin, legend_margin, legend_margin, 'pt'),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(2, 2, 2, 2, 'pt'),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # 1. Depth (Bathymetry) - with label
  plot_depth <- ggplot() +
    geom_spatraster(data = bathy_final, maxcell = maxcell_param) +
    scale_fill_gradientn(
      name = "meters",
      colors = rev(cmocean("deep")(256)),
      limits = c(-60, 0),
      breaks = c(-60, -45, -30, -15, 0),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = colorbar_width,
        barheight = colorbar_height,
        label.vjust = label_vjust,
        title.vjust = title_vjust,
        ticks.linewidth = ticks_linewidth,
        frame.colour = NA
      )
    ) +
    geom_spatvector(data = land_vect, fill = "gray90", color = "black", linewidth = land_linewidth) +
    annotate("text", x = -Inf, y = Inf, label = "a", 
             hjust = panel_label_hjust, vjust = panel_label_vjust, 
             size = panel_label_size, family = "Georgia", fontface = panel_label_fontface) +
    coord_sf(xlim = c(xmin_utm, xmax_utm), 
             ylim = c(ymin_utm, ymax_utm),
             expand = FALSE) +
    common_theme
  
  # 2. Mean SST - with label
  plot_sst <- ggplot() +
    geom_spatraster(data = mean_sst_raster, maxcell = maxcell_param) +
    scale_fill_gradientn(
      name = "°C",
      colors = cmocean("thermal")(256),
      limits = c(27.5, 28.2),
      breaks = c(27.5, 27.85, 28.2),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = colorbar_width,
        barheight = colorbar_height,
        label.vjust = label_vjust,
        title.vjust = title_vjust,
        ticks.linewidth = ticks_linewidth,
        frame.colour = NA
      )
    ) +
    geom_spatvector(data = land_vect, fill = "gray90", color = "black", linewidth = land_linewidth) +
    annotate("text", x = -Inf, y = Inf, label = "b", 
             hjust = panel_label_hjust, vjust = panel_label_vjust, 
             size = panel_label_size, family = "Georgia", fontface = panel_label_fontface) +
    coord_sf(xlim = c(xmin_utm, xmax_utm), 
             ylim = c(ymin_utm, ymax_utm),
             expand = FALSE) +
    common_theme
  
  # 3. Mean Direction - with label
  plot_dir <- ggplot() +
    geom_spatraster(data = dir_erddap_raster, maxcell = maxcell_param) +
    scale_fill_gradientn(
      name = "dir (°)",
      colors = cmocean("phase")(256),
      limits = c(0, 360),
      breaks = c(0, 120, 240, 360),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = colorbar_width,
        barheight = colorbar_height,
        label.vjust = label_vjust,
        title.vjust = title_vjust,
        ticks.linewidth = ticks_linewidth,
        frame.colour = NA
      )
    ) +
    geom_spatvector(data = land_vect, fill = "gray90", color = "black", linewidth = land_linewidth) +
    annotate("text", x = -Inf, y = Inf, label = "c", 
             hjust = panel_label_hjust, vjust = panel_label_vjust, 
             size = panel_label_size, family = "Georgia", fontface = panel_label_fontface) +
    coord_sf(xlim = c(xmin_utm, xmax_utm), 
             ylim = c(ymin_utm, ymax_utm),
             expand = FALSE) +
    common_theme
  
  # 4. Mean SPM - with label
  # Custom break scheme optimized for data distribution
  spm_breaks <- c(0.07, 0.4, 2, 10, 26)
  spm_log_breaks <- log10(spm_breaks)
  spm_values <- (spm_log_breaks - min(spm_log_breaks)) / (max(spm_log_breaks) - min(spm_log_breaks))
  
  plot_spm <- ggplot() +
    geom_spatraster(data = mean_spm_raster, maxcell = maxcell_param) +
    scale_fill_gradientn(
      name = "mg/m³",
      colors = cmocean("turbid")(256),
      values = spm_values,
      na.value = "transparent",
      trans = "log10",
      limits = c(0.07, 26),
      breaks = spm_breaks,
      labels = c("0.07", "0.4", "2", "10", "26"),
      guide = guide_colorbar(
        barwidth = colorbar_width,
        barheight = colorbar_height,
        label.vjust = label_vjust,
        title.vjust = title_vjust,
        ticks.linewidth = ticks_linewidth,
        frame.colour = NA
      )
    ) +
    geom_spatvector(data = land_vect, fill = "gray90", color = "black", linewidth = land_linewidth) +
    annotate("text", x = -Inf, y = Inf, label = "d", 
             hjust = panel_label_hjust, vjust = panel_label_vjust, 
             size = panel_label_size, family = "Georgia", fontface = panel_label_fontface) +
    coord_sf(xlim = c(xmin_utm, xmax_utm), 
             ylim = c(ymin_utm, ymax_utm),
             expand = FALSE) +
    common_theme
  
  # 5. TPI - with label
  # Symmetric divergent breaks centered on 0, with squish to handle out-of-bounds values
  tpi_breaks_full <- c(-50, -10, -1, 0, 1, 10, 50)
  tpi_values <- (tpi_breaks_full - min(tpi_breaks_full)) / (max(tpi_breaks_full) - min(tpi_breaks_full))
  tpi_breaks_display <- c(-50, -10, 10, 50)
  
  plot_tpi <- ggplot() +
    geom_spatraster(data = TPI_terra, maxcell = maxcell_param) +
    scale_fill_gradientn(
      name = "index",
      colors = cmocean("balance")(256),
      values = tpi_values,
      na.value = "transparent",
      limits = c(-50, 50),
      breaks = tpi_breaks_display,
      labels = c("-50", "-10", "10", "50"),
      oob = scales::squish,  # Squish out-of-bounds values to nearest color
      guide = guide_colorbar(
        barwidth = colorbar_width,
        barheight = colorbar_height,
        label.vjust = label_vjust,
        title.vjust = title_vjust,
        ticks.linewidth = ticks_linewidth,
        frame.colour = NA
      )
    ) +
    geom_spatvector(data = land_vect, fill = "gray90", color = "black", linewidth = land_linewidth) +
    annotate("text", x = -Inf, y = Inf, label = "e", 
             hjust = panel_label_hjust, vjust = panel_label_vjust, 
             size = panel_label_size, family = "Georgia", fontface = panel_label_fontface) +
    coord_sf(xlim = c(xmin_utm, xmax_utm), 
             ylim = c(ymin_utm, ymax_utm),
             expand = FALSE) +
    common_theme
  
  # 6. Mean PAR - with label
  plot_par <- ggplot() +
    geom_spatraster(data = mean_par_raster, maxcell = maxcell_param) +
    scale_fill_gradientn(
      name = "Einstein/m²/d",
      colors = cmocean("solar")(256),
      limits = c(37, 48),
      breaks = c(37, 42.5, 48),
      na.value = "transparent",
      guide = guide_colorbar(
        barwidth = colorbar_width,
        barheight = colorbar_height,
        label.vjust = label_vjust,
        title.vjust = title_vjust,
        ticks.linewidth = ticks_linewidth,
        frame.colour = NA
      )
    ) +
    geom_spatvector(data = land_vect, fill = "gray90", color = "black", linewidth = land_linewidth) +
    annotate("text", x = -Inf, y = Inf, label = "f", 
             hjust = panel_label_hjust, vjust = panel_label_vjust, 
             size = panel_label_size, family = "Georgia", fontface = panel_label_fontface) +
    coord_sf(xlim = c(xmin_utm, xmax_utm), 
             ylim = c(ymin_utm, ymax_utm),
             expand = FALSE) +
    common_theme
  
  ################################## Combine Panels ##################################
  
  # Create 3x2 grid using patchwork
  enviro_multipanel <- (plot_depth + plot_sst) /
    (plot_dir + plot_spm) /
    (plot_tpi + plot_par)
  
  # enviro_multipanel
  
  ################################## Save Plots ##################################
  
  # Save PNG
  ggsave(
    filename = here(output_dir, "fig2.png"),
    plot = enviro_multipanel,
    width = 7,
    height = 4,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  
  # # Save PDF
  # ggsave(
  #   filename = here(output_dir, "environmental_rasters_constrained.pdf"),
  #   plot = enviro_multipanel,
  #   width = 7,
  #   height = 4,
  #   device = cairo_pdf,
  #   limitsize = FALSE
  # )
  