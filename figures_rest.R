# .rs.restartR(clean = TRUE)
rm(list=ls())

library(here)
library(dplyr)
library(ggplot2)
library(terra)
library(tidyterra)
library(viridis)
library(cowplot)
library(cmocean)
library(extrafont)
library(ggnewscale)
source(here("src/functions.R"))

# Display colorblind-friendly ColorBrewer palettes (comment out when not needed)
# RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)

# Load fonts (run loadfonts() if Georgia not available)
# loadfonts(device = "all", quiet = TRUE)

################################## CONTROL TOGGLE ##################################

# Raster constraint - TRUE for maxcell = Inf, FALSE for default downsampling
UNCONSTRAINED_RASTER <- FALSE  # Options: TRUE or FALSE

# Set maxcell parameter based on toggle
maxcell_param <- if(UNCONSTRAINED_RASTER) Inf else 5e5

# Depth contour settings
show_depth_contour <- TRUE  # Toggle for -30m depth contour
contour_depth <- -30  # Depth for contour line
contour_color <- "gray30"  # Contour line color
contour_width <- 0.3  # Contour line width
contour_alpha <- 0.5  # Contour line transparency (0-1)

# Color palette for datasets (colorblind-friendly with high contrast against deep blue bathymetry)
# Based on well-tested colorblind-safe palette with maximum separation
dataset_colors <- c(
  "DCRMP" = "#8B4513",    # Dark Brown
  "DeepLion" = "#F0E442", # Bright Yellow
  "NCRMP" = "#009E73",    # Bluish Green/Teal
  "NODICE" = "#FFFFFF",   # White (with black border)
  "PRCRMP" = "#9400D3",   # Dark Violet/Purple (distinct from all others)
  "SESAP" = "#CC79A7",    # Rose Pink
  "TCRMP" = "#333333"     # Dark Gray (almost black)
)

# Point offset settings (for visibility over bathymetry)
use_point_offset <- TRUE  # Toggle for offset points with lines
offset_distance <- 0.0075  # Distance to offset points (in degrees)

# Margin settings (from reference script)
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

# Panel dimensions
bottom_height <- 1
topleft_width <- 1.2

# Legend positions
legend_depth_x <- 0.9  # Shifted left from 1.0
legend_depth_y <- 0.0  # Bottom of top-left panel
legend_dataset_x <- 0.0  # Left side of bottom panel
legend_dataset_y <- 1.0  # Top of bottom panel

# Text sizes
textsize_stt <- 2.0
textsize_stx <- 2.5
textsize_whole <- 3.0

# Point sizes
pointsize_stt <- 0.8
pointsize_stx <- 0.8
pointsize_whole <- 0.6

# Point alpha (transparency)
point_alpha <- 0.7

################################## Setup ##################################

load(here("output", "all_combined_data.rda"))
land_vect <- readRDS(here("output", "osm_land_vect.rds"))

# Get unique PSU locations
psu_locations <- combined_benthic_data_averaged_psu %>%
  select(PSU, lat, lon, dataset) %>%
  distinct(PSU, .keep_all = TRUE) %>%
  mutate(dataset = case_when(
    dataset == "NCRMP_benthic" ~ "NCRMP",
    dataset == "TCRMP_benthic" ~ "TCRMP",
    TRUE ~ dataset
  ))

# Load and prepare bathymetry
spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
bathy_latlon <- project(bathy_final, "EPSG:4326")

# Define island labels
labels_stt_stj <- data.frame(
  name = c("Saint Thomas", "Saint John", "Culebra", "Vieques"),
  x = c(-64.93106, -64.75019, -65.28551, -65.34985),
  y = c(18.3485, 18.33596, 18.318, 18.13872)
)
labels_stx <- data.frame(name = c("Saint Croix"), x = c(-64.77806), y = c(17.73595))
labels_whole <- data.frame(name = c("Puerto Rico"), x = c(-66.41416), y = c(18.23930))

# Define regions
regions <- list(
  stt_stj = list(xlim = c(-65.39305, -64.65467), ylim = c(18.06478, 18.46411)),
  stx = list(xlim = c(-64.9474, -64.4094), ylim = c(17.6145, 17.8561)),
  whole = list(xlim = c(-67.9895, -63.9886), ylim = c(17.6145, 18.8255))
)

################################## Create Multi-panel Map ##################################

# Prepare offset data if enabled
if(use_point_offset) {
  psu_locations <- psu_locations %>%
    mutate(
      lat_offset = lat + offset_distance,
      lon_offset = lon
    )
}

# Top left: St. Thomas/St. John
bathy_stt <- crop(bathy_latlon, ext(regions$stt_stj$xlim[1], regions$stt_stj$xlim[2], 
                                    regions$stt_stj$ylim[1], regions$stt_stj$ylim[2]))
plot_stt <- ggplot() +
  geom_spatraster(data = bathy_stt, maxcell = maxcell_param) +
  scale_fill_gradientn(name = "Depth (m)", colors = rev(cmocean("deep")(256)), limits = c(-60, 0), 
                       breaks = c(-60, 0), labels = c("-60", "0"), na.value = "transparent") +
  {if(show_depth_contour) geom_contour(data = as.data.frame(bathy_stt, xy = TRUE, na.rm = TRUE), 
                                       aes(x = x, y = y, z = .data[[names(bathy_stt)]]), 
                                       breaks = contour_depth, color = contour_color, 
                                       linewidth = contour_width, alpha = contour_alpha)} +
  new_scale_fill() +
  geom_spatvector(data = land_vect, fill = "tan", color = "black", linewidth = 0.3) +
  geom_text(data = labels_stt_stj, aes(x = x, y = y, label = name), size = textsize_stt, family = "Georgia", color = "black") +
  {if(use_point_offset) geom_segment(data = psu_locations, 
                                     aes(x = lon, y = lat, xend = lon_offset, yend = lat_offset),
                                     linewidth = 0.2, color = "gray30")} +
  geom_point(data = psu_locations, 
             aes(x = if(use_point_offset) lon_offset else lon, 
                 y = if(use_point_offset) lat_offset else lat, 
                 fill = dataset), 
             size = pointsize_stt, alpha = point_alpha, shape = 21, color = "black", stroke = 0.5) +
  scale_fill_manual(name = "Dataset:", values = dataset_colors) +
  coord_sf(xlim = regions$stt_stj$xlim, ylim = regions$stt_stj$ylim, expand = FALSE) +
  theme_classic(base_family = "Georgia") +
  theme(legend.position = c(legend_depth_x, legend_depth_y),
        legend.justification = c("right", "bottom"),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(0.35, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.spacing.x = unit(0.1, "cm"),
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(margin_stt_t, margin_stt_r, margin_stt_b, margin_stt_l, 'pt'),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  guides(fill = guide_colorbar(order = 1), 
         fill = guide_legend(override.aes = list(size = 3), order = 2))

# Top right: St. Croix
bathy_stx <- crop(bathy_latlon, ext(regions$stx$xlim[1], regions$stx$xlim[2], 
                                    regions$stx$ylim[1], regions$stx$ylim[2]))
plot_stx <- ggplot() +
  geom_spatraster(data = bathy_stx, maxcell = maxcell_param) +
  scale_fill_gradientn(name = "Depth", colors = rev(cmocean("deep")(256)), limits = c(-60, 0), na.value = "transparent", guide = "none") +
  {if(show_depth_contour) geom_contour(data = as.data.frame(bathy_stx, xy = TRUE, na.rm = TRUE), 
                                       aes(x = x, y = y, z = .data[[names(bathy_stx)]]), 
                                       breaks = contour_depth, color = contour_color, 
                                       linewidth = contour_width, alpha = contour_alpha)} +
  new_scale_fill() +
  geom_spatvector(data = land_vect, fill = "tan", color = "black", linewidth = 0.3) +
  geom_text(data = labels_stx, aes(x = x, y = y, label = name), size = textsize_stx, family = "Georgia", color = "black") +
  {if(use_point_offset) geom_segment(data = psu_locations, 
                                     aes(x = lon, y = lat, xend = lon_offset, yend = lat_offset),
                                     linewidth = 0.2, color = "gray30")} +
  geom_point(data = psu_locations, 
             aes(x = if(use_point_offset) lon_offset else lon, 
                 y = if(use_point_offset) lat_offset else lat, 
                 fill = dataset), 
             size = pointsize_stx, alpha = point_alpha, shape = 21, color = "black", stroke = 0.5) +
  scale_fill_manual(name = "Dataset:", values = dataset_colors) +
  coord_sf(xlim = regions$stx$xlim, ylim = regions$stx$ylim, expand = FALSE) +
  theme_classic(base_family = "Georgia") +
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(margin_stx_t, margin_stx_r, margin_stx_b, margin_stx_l, 'pt'),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Bottom: Full region
bathy_whole <- crop(bathy_latlon, ext(regions$whole$xlim[1], regions$whole$xlim[2], 
                                      regions$whole$ylim[1], regions$whole$ylim[2]))
plot_whole <- ggplot() +
  geom_spatraster(data = bathy_whole, maxcell = maxcell_param) +
  scale_fill_gradientn(name = "Depth", colors = rev(cmocean("deep")(256)), limits = c(-60, 0), na.value = "transparent", guide = "none") +
  {if(show_depth_contour) geom_contour(data = as.data.frame(bathy_whole, xy = TRUE, na.rm = TRUE), 
                                       aes(x = x, y = y, z = .data[[names(bathy_whole)]]), 
                                       breaks = contour_depth, color = contour_color, 
                                       linewidth = contour_width, alpha = contour_alpha)} +
  new_scale_fill() +
  geom_spatvector(data = land_vect, fill = "tan", color = "black", linewidth = 0.3) +
  geom_text(data = labels_whole, aes(x = x, y = y, label = name), size = textsize_whole, family = "Georgia", color = "black") +
  {if(use_point_offset) geom_segment(data = psu_locations, 
                                     aes(x = lon, y = lat, xend = lon_offset, yend = lat_offset),
                                     linewidth = 0.2, color = "gray30")} +
  geom_point(data = psu_locations, 
             aes(x = if(use_point_offset) lon_offset else lon, 
                 y = if(use_point_offset) lat_offset else lat, 
                 fill = dataset), 
             size = pointsize_whole, alpha = point_alpha, shape = 21, color = "black", stroke = 0.5) +
  scale_fill_manual(name = "Dataset:", values = dataset_colors) +
  coord_sf(xlim = regions$whole$xlim, ylim = regions$whole$ylim, expand = FALSE) +
  theme_classic(base_family = "Georgia") +
  theme(legend.position = c(legend_dataset_x, legend_dataset_y), 
        legend.justification = c("left", "top"),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_text(size = 8, margin = margin(r = 2)),
        legend.text = element_text(size = 7, margin = margin(l = 1, r = 2)),
        legend.spacing.x = unit(0.02, "cm"),
        legend.key.spacing.x = unit(0.02, "cm"),
        legend.key.width = unit(0.4, "cm"),
        axis.title = element_blank(),
        axis.text = element_text(size = 6),
        plot.margin = margin(margin_whole_t, margin_whole_r, margin_whole_b, margin_whole_l, 'pt'),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  guides(fill = guide_legend(override.aes = list(size = 3), nrow = 1, byrow = TRUE))

# Combine panels
top_row <- plot_grid(plot_stt, plot_stx, nrow = 1, align = 'h', rel_widths = c(topleft_width, 1))
map_multipanel <- plot_grid(top_row, plot_whole, ncol = 1, align = 'v', rel_heights = c(1, bottom_height))

# print(map_multipanel)

# Save plots
ggsave(filename = here("output", "output_figures_tables", "fig1.png"),
       plot = map_multipanel,
       width = 7,
       height = 5.5,
       dpi = 300,
       bg = "white",
       limitsize = FALSE)

ggsave(filename = here("output", "output_figures_tables", "fig1.pdf"),
       plot = map_multipanel,
       width = 7,
       height = 5.5,
       device = cairo_pdf,
       limitsize = FALSE)