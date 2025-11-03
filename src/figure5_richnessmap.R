# .rs.restartR(clean = TRUE)
# rm(list=ls())

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
UNCONSTRAINED_RASTER <- FALSE  # Options: TRUE or FALSE

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

# cat("Creating species richness raster...\n")
# 
# # Create richness raster (count of species with cover > 0 at each pixel)
# richness_raster <- Reduce(`+`, lapply(raster_list, function(r) r > 0))
# 
# # Restore CRS
# terra::crs(richness_raster) <- spatial_metadata$bathy_final$crs
# 
# # Project to lat/lon
# richness_raster_latlon <- project(richness_raster, "EPSG:4326")
# 
# # Crop to domain
# domain_extent <- ext(-67.99, -63.99, 17.63, 18.83)
# richness_raster_cropped <- crop(richness_raster_latlon, domain_extent)
# 
# cat("Richness raster created\n")

################################## Create Dominant Species Raster ##################################

cat("Creating dominant species raster...\n")

# Stack all species rasters
species_stack <- rast(raster_list)
names(species_stack) <- names(raster_list)

# Create a mask where ANY species has NA
na_mask <- any(is.na(species_stack))

# Find which species has maximum cover at each pixel
dominant_species <- which.max(species_stack)

# Set pixels where all species have 0 cover to value 0 (for gray)
total_cover <- sum(species_stack)
dominant_species[total_cover == 0] <- 0

# Apply NA mask - set to NA anywhere ANY species has NA
dominant_species[na_mask] <- NA

# Restore CRS
terra::crs(dominant_species) <- spatial_metadata$bathy_final$crs

# Project to lat/lon
dominant_species_latlon <- project(dominant_species, "EPSG:4326", method = "near")

# Crop to domain
domain_extent <- ext(-67.99, -63.99, 17.63, 18.83)
dominant_species_cropped <- crop(dominant_species_latlon, domain_extent)

# Rename species to desired codes
species_codes <- c("AGAR", "CNAT", "DCYL", "DSTO", "DLAB", "EFAS", 
                   "MADR", "MEAN", "MCAV", "MYCE", "ORBI", "PORI", 
                   "PSEU", "SIDE", "SINT")

# Define display names and desired order
display_names <- c("Agaricia", "Madracis", "Porites", "Siderastrea", "S. intersepta", 
                   "M. cavernosa", "Orbicella", "C. natans", "D. cylindrus", "D. stokesii", 
                   "D. labyrinthiformis", "E. fastigiata", "Meandrina", "Mycetophyllia", "Pseudodiploria")

# Susceptibility groups
sus_groups <- c("LS", "LS", "LS", "LS", "LS",
                "MS", "MS",
                "HS", "HS", "HS", "HS", "HS", "HS", "HS", "HS")

# Create mapping from original raster_list order to display order indices
# raster_list: agaricia, colpophyllia, dendrogyra, dichocoenia, diploria, eusmilia, madracis, meandrina, montastraea, mycetophyllia, orbicella, porites, pseudodiploria, siderastrea, stephanocoenia
display_order <- c(1, 8, 9, 10, 11, 12, 2, 13, 6, 14, 7, 3, 15, 4, 5)

# Find which species are actually dominant somewhere
dominant_values <- unique(values(dominant_species_cropped))
dominant_values <- dominant_values[!is.na(dominant_values)]
dominant_indices <- sort(unique(dominant_values[dominant_values > 0]))  # Exclude 0

# Subset to only species that are dominant somewhere
species_codes_present <- species_codes[dominant_indices]

n_species <- length(species_codes_present)

cat(paste0("Dominant species raster created with ", n_species, " species\n"))

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

create_dominant_plot <- function(region, raster_data, land_data, species_codes, dominant_indices, show_legend = TRUE,
                                 textsize = 2.5, use_maxcell_inf = TRUE) {
  region_labels <- get_region_labels(region$name)
  buffer <- 0.1
  region_extent <- ext(region$xlim[1]-buffer, region$xlim[2]+buffer, region$ylim[1]-buffer, region$ylim[2]+buffer)
  raster_region <- crop(raster_data, region_extent)
  
  maxcell_param <- if(use_maxcell_inf) Inf else 5e5
  
  # Create color palette for species present
  n_species <- length(species_codes)
  species_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_species)
  
  ggplot() +
    geom_spatraster(data = raster_region, maxcell = maxcell_param) +
    scale_fill_gradientn(name = "",
                         colors = c("gray85", species_colors),
                         values = scales::rescale(c(0, dominant_indices)),
                         breaks = dominant_indices,
                         labels = species_codes,
                         na.value = "transparent") +
    guides(fill = guide_legend(nrow = 3, ncol = 5, byrow = TRUE)) +
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

################################## Main Processing - Richness ##################################

# cat("\n========================================\n")
# cat("Processing: Species Richness\n")
# cat(paste0("Raster mode: ", if(UNCONSTRAINED_RASTER) "Unconstrained (maxcell = Inf)" else "Constrained (default)", "\n"))
# cat("========================================\n")
# 
# max_richness <- max(values(richness_raster_cropped), na.rm = TRUE)
# cat(paste0("Max species richness: ", round(max_richness, 0), " species\n"))
# 
# # Create plots with individual margins
# plot_stt_stj <- create_region_plot(regions$stt_stj, richness_raster_cropped, land_vect, 
#                                    max_richness, FALSE, textsize = 2.0, use_maxcell_inf = UNCONSTRAINED_RASTER) + 
#   theme(plot.margin = margin(margin_stt_t, margin_stt_r, margin_stt_b, margin_stt_l, 'pt'),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         plot.background = element_rect(fill = "white", color = NA),
#         panel.background = element_rect(fill = "white", color = NA))
# 
# plot_stx <- create_region_plot(regions$stx, richness_raster_cropped, land_vect, 
#                                max_richness, TRUE, textsize = 2.5, use_maxcell_inf = UNCONSTRAINED_RASTER) +
#   scale_fill_viridis_c(name = "Richness", 
#                        option = "viridis", 
#                        limits = c(0, max_richness), 
#                        breaks = c(0, max_richness),
#                        labels = c("0", round(max_richness, 0)),
#                        na.value = "transparent") +
#   theme(plot.margin = margin(margin_stx_t, margin_stx_r, margin_stx_b, margin_stx_l, 'pt'),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         legend.position = c(legend_x, legend_y),
#         legend.justification = c("right", "bottom"),
#         legend.direction = "horizontal",
#         legend.background = element_rect(fill = "transparent", color = NA),
#         legend.box.background = element_rect(fill = "transparent", color = NA),
#         legend.key = element_rect(fill = "transparent", color = NA),
#         legend.key.size = unit(0.35, "cm"),
#         legend.title = element_text(size = 8),
#         legend.text = element_text(size = 7),
#         plot.background = element_rect(fill = "white", color = NA),
#         panel.background = element_rect(fill = "white", color = NA)) +
#   annotate("text", x = regions$stx$xlim[1] + 0.05, y = regions$stx$ylim[2] - 0.01,
#            label = "Species richness", size = 4, fontface = "bold", family = "Georgia", hjust = 0)
# 
# plot_whole <- create_region_plot(regions$whole, richness_raster_cropped, land_vect, 
#                                  max_richness, FALSE, textsize = 3.0, use_maxcell_inf = UNCONSTRAINED_RASTER) +
#   theme(plot.margin = margin(margin_whole_t, margin_whole_r, margin_whole_b, margin_whole_l, 'pt'),
#         plot.background = element_rect(fill = "white", color = NA),
#         panel.background = element_rect(fill = "white", color = NA))
# 
# # Combine plots
# top_row <- plot_grid(plot_stt_stj, plot_stx, nrow = 1, align = 'h', rel_widths = c(topleft_width, 1))
# final_plot <- plot_grid(top_row, plot_whole, ncol = 1, align = 'v', rel_heights = c(1, bottom_height), axis = 'tb')
# 
# # Add suffix for raster mode
# raster_suffix <- if(UNCONSTRAINED_RASTER) "_unconstrained" else "_constrained"
# 
# # Save plots
# cat("Saving richness plots...\n")
# ggsave(filename = here("output", "output_figures_tables", paste0("fig_combined_species_richness", raster_suffix, ".png")),
#        plot = final_plot,
#        width = 7.087,
#        height = 5.5,
#        dpi = 300,
#        bg = "white",
#        limitsize = FALSE)
# 
# cat("✓ Species richness map saved\n")

################################## Main Processing - Dominant Species ##################################

cat("\n========================================\n")
cat("Processing: Dominant Species\n")
cat("========================================\n")

# Create dominant species plots
plot_dom_stt_stj <- create_dominant_plot(regions$stt_stj, dominant_species_cropped, land_vect,
                                         species_codes_present, dominant_indices, FALSE, textsize = 2.0, use_maxcell_inf = UNCONSTRAINED_RASTER) +
  theme(plot.margin = margin(margin_stt_t, margin_stt_r, margin_stt_b, margin_stt_l, 'pt'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

plot_dom_stx <- create_dominant_plot(regions$stx, dominant_species_cropped, land_vect,
                                     species_codes_present, dominant_indices, TRUE, textsize = 2.5, use_maxcell_inf = UNCONSTRAINED_RASTER) +
  theme(plot.margin = margin(margin_stx_t, margin_stx_r, margin_stx_b, margin_stx_l, 'pt'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

plot_dom_whole <- create_dominant_plot(regions$whole, dominant_species_cropped, land_vect,
                                       species_codes_present, dominant_indices, FALSE, textsize = 3.0, use_maxcell_inf = UNCONSTRAINED_RASTER) +
  theme(plot.margin = margin(margin_whole_t, margin_whole_r, margin_whole_b, margin_whole_l, 'pt'),
        legend.position = c(0.35, 0.85),
        legend.justification = c("center", "bottom"),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0.98), color = "black", linewidth = 0.3),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.spacing.x = unit(0.1, "cm"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Combine dominant species plots
top_row_dom <- plot_grid(plot_dom_stt_stj, plot_dom_stx, nrow = 1, align = 'h', rel_widths = c(topleft_width, 1))
final_plot_dom <- plot_grid(top_row_dom, plot_dom_whole, ncol = 1, align = 'v', rel_heights = c(1, bottom_height), axis = 'tb')

# Add suffix for raster mode
raster_suffix <- if(UNCONSTRAINED_RASTER) "_unconstrained" else "_constrained"

# Save dominant species plots
cat("Saving dominant species plots...\n")
ggsave(filename = here("output", "output_figures_tables", paste0("fig_combined_dominant_species", raster_suffix, ".png")),
       plot = final_plot_dom,
       width = 7.087,
       height = 5.5,
       dpi = 300,
       bg = "white",
       limitsize = FALSE)

cat("✓ Dominant species map saved\n")

################################## Calculate Species Extent ##################################

cat("\n========================================\n")
cat("Calculating species extent (km²)...\n")
cat("========================================\n")

# Grid cell area in km²
cell_area_km2 <- 0.05 * 0.05  # 50m x 50m = 0.0025 km²

# Create depth masks from bathymetry
# Shallow: < 30m depth (bathymetry values are negative, so > -30)
# Mesophotic: 30-60m depth (between -60 and -30)
shallow_mask <- bathy_final > -30
mesophotic_mask <- (bathy_final <= -30) & (bathy_final >= -60)

# Calculate extent for each species
species_extent <- data.frame(
  Taxon = display_names,
  Sus = sus_groups,
  Shallow = numeric(length(species_codes)),
  Mesophotic = numeric(length(species_codes)),
  Total = numeric(length(species_codes))
)

for (i in 1:length(raster_list)) {
  species_raster <- raster_list[[i]]
  
  # Map to display order
  display_idx <- display_order[i]
  
  # Shallow extent (< 30m, cover-weighted)
  shallow_extent <- sum(values(species_raster * shallow_mask) * cell_area_km2, na.rm = TRUE)
  species_extent$Shallow[display_idx] <- shallow_extent
  
  # Mesophotic extent (30-60m, cover-weighted)
  mesophotic_extent <- sum(values(species_raster * mesophotic_mask) * cell_area_km2, na.rm = TRUE)
  species_extent$Mesophotic[display_idx] <- mesophotic_extent
  
  # Total extent (cover-weighted)
  total_extent <- sum(values(species_raster) * cell_area_km2, na.rm = TRUE)
  species_extent$Total[display_idx] <- total_extent
}

# Calculate totals
# For cover-weighted: sum all cover across all species at each pixel
total_cover_stack <- Reduce(`+`, raster_list)
total_extent_km2 <- sum(values(total_cover_stack) * cell_area_km2, na.rm = TRUE)
total_shallow_km2 <- sum(values(total_cover_stack * shallow_mask) * cell_area_km2, na.rm = TRUE)
total_mesophotic_km2 <- sum(values(total_cover_stack * mesophotic_mask) * cell_area_km2, na.rm = TRUE)

# Keep species in original order (as listed in display_names)
# species_extent already in the correct order from the loop

################################## Print Table ##################################

# Print cover-weighted table
cat("\nTable: Spatial extent (km²) predicted by the combined SDM-SAM hurdle model.\n")
cat("Mesophotic depths are 30-60 m.\n\n")
cat("COVER-WEIGHTED EXTENT:\n")

# Calculate percentages based on TOTAL coral extent (not by depth zone)
species_extent$Shallow_pct <- (species_extent$Shallow / total_extent_km2) * 100
species_extent$Mesophotic_pct <- (species_extent$Mesophotic / total_extent_km2) * 100
species_extent$Total_pct <- (species_extent$Total / total_extent_km2) * 100

# Calculate group totals
ls_indices <- which(species_extent$Sus == "LS")
ms_indices <- which(species_extent$Sus == "MS")
hs_indices <- which(species_extent$Sus == "HS")

ls_total <- data.frame(
  Taxon = "Total", Sus = "",
  Shallow = sum(species_extent$Shallow[ls_indices]),
  Mesophotic = sum(species_extent$Mesophotic[ls_indices]),
  Total = sum(species_extent$Total[ls_indices]),
  Shallow_pct = (sum(species_extent$Shallow[ls_indices]) / total_extent_km2) * 100,
  Mesophotic_pct = (sum(species_extent$Mesophotic[ls_indices]) / total_extent_km2) * 100,
  Total_pct = (sum(species_extent$Total[ls_indices]) / total_extent_km2) * 100
)

ms_total <- data.frame(
  Taxon = "Total", Sus = "",
  Shallow = sum(species_extent$Shallow[ms_indices]),
  Mesophotic = sum(species_extent$Mesophotic[ms_indices]),
  Total = sum(species_extent$Total[ms_indices]),
  Shallow_pct = (sum(species_extent$Shallow[ms_indices]) / total_extent_km2) * 100,
  Mesophotic_pct = (sum(species_extent$Mesophotic[ms_indices]) / total_extent_km2) * 100,
  Total_pct = (sum(species_extent$Total[ms_indices]) / total_extent_km2) * 100
)

hs_total <- data.frame(
  Taxon = "Total", Sus = "",
  Shallow = sum(species_extent$Shallow[hs_indices]),
  Mesophotic = sum(species_extent$Mesophotic[hs_indices]),
  Total = sum(species_extent$Total[hs_indices]),
  Shallow_pct = (sum(species_extent$Shallow[hs_indices]) / total_extent_km2) * 100,
  Mesophotic_pct = (sum(species_extent$Mesophotic[hs_indices]) / total_extent_km2) * 100,
  Total_pct = (sum(species_extent$Total[hs_indices]) / total_extent_km2) * 100
)

grand_total <- data.frame(
  Taxon = "Grand Total", Sus = "",
  Shallow = total_shallow_km2,
  Mesophotic = total_mesophotic_km2,
  Total = total_extent_km2,
  Shallow_pct = (total_shallow_km2 / total_extent_km2) * 100,
  Mesophotic_pct = (total_mesophotic_km2 / total_extent_km2) * 100,
  Total_pct = 100.0
)

# Combine with group totals inserted
species_with_totals <- rbind(
  species_extent[1:5, ],
  ls_total,
  species_extent[6:7, ],
  ms_total,
  species_extent[8:15, ],
  hs_total,
  grand_total
)

# Format table with area and percentage columns
table1 <- data.frame(
  Taxon = species_with_totals$Taxon,
  Sus = species_with_totals$Sus,
  Shallow_Area = sprintf("%.2f", species_with_totals$Shallow),
  Shallow_Pct = sprintf("%.2f", species_with_totals$Shallow_pct),
  Mesophotic_Area = sprintf("%.2f", species_with_totals$Mesophotic),
  Mesophotic_Pct = sprintf("%.2f", species_with_totals$Mesophotic_pct),
  Total_Area = sprintf("%.2f", species_with_totals$Total),
  Total_Pct = sprintf("%.2f", species_with_totals$Total_pct)
)

# Rename columns for display
names(table1) <- c("Taxon", "Sus.", "Area (km²)", "%", "Area (km²)", "%", "Area (km²)", "%")

print(table1, row.names = FALSE)

# Print total available seafloor area
cat("\n========================================\n")
cat("Total available seafloor area:\n")
cat("========================================\n")

# Count non-NA cells in first species raster (all should have same NA pattern)
total_valid_cells <- sum(!is.na(values(raster_list[[1]])))
total_seafloor_km2 <- total_valid_cells * cell_area_km2

# Count by depth zone
shallow_valid_cells <- sum(values(!is.na(raster_list[[1]]) & shallow_mask), na.rm = TRUE)
mesophotic_valid_cells <- sum(values(!is.na(raster_list[[1]]) & mesophotic_mask), na.rm = TRUE)

shallow_seafloor_km2 <- shallow_valid_cells * cell_area_km2
mesophotic_seafloor_km2 <- mesophotic_valid_cells * cell_area_km2

cat(sprintf("Shallow (< 30m): %.1f km²\n", shallow_seafloor_km2))
cat(sprintf("Mesophotic (30-60m): %.1f km²\n", mesophotic_seafloor_km2))
cat(sprintf("Total: %.1f km²\n", total_seafloor_km2))
cat(sprintf("Total number of valid cells: %d\n", total_valid_cells))

# Count cells with non-zero vs zero coral cover (predicted)
nonzero_cells <- sum(values(total_cover_stack) > 0, na.rm = TRUE)
zero_cells <- sum(values(total_cover_stack) == 0, na.rm = TRUE)

nonzero_km2 <- nonzero_cells * cell_area_km2
zero_km2 <- zero_cells * cell_area_km2

nonzero_pct <- (nonzero_km2 / total_seafloor_km2) * 100
zero_pct <- (zero_km2 / total_seafloor_km2) * 100

# Print observation summary
cat("\n========================================\n")
cat("Field observation coverage summary:\n")
cat("========================================\n")

# Get all species names
all_species <- names(all_results)

# Calculate observed total cover by PSU
psu_total_cover <- combined_benthic_data_averaged %>%
  mutate(genus = tolower(sub(" .*", "", spp))) %>%
  filter(genus %in% all_species) %>%
  group_by(PSU) %>%
  summarise(lat = first(lat), 
            lon = first(lon), 
            total_cover = sum(cover, na.rm = TRUE), 
            .groups = 'drop')

# Create spatial points in lat/lon (WGS84)
obs_points <- vect(psu_total_cover, geom = c("lon", "lat"), crs = "EPSG:4326")

# Transform to match raster CRS (UTM)
obs_points_utm <- project(obs_points, crs(total_cover_stack))

# Extract cell indices for observation locations
obs_cells <- cellFromXY(total_cover_stack, crds(obs_points_utm))
obs_cells_df <- data.frame(
  cell = obs_cells,
  total_cover = psu_total_cover$total_cover
)

# Get unique cells and their properties
unique_obs_cells_df <- obs_cells_df %>%
  group_by(cell) %>%
  summarise(max_cover = max(total_cover), .groups = 'drop')

unique_obs_cells <- nrow(unique_obs_cells_df)

# Extract depth values for observed cells
unique_obs_cells_df$depth <- values(bathy_final)[unique_obs_cells_df$cell]

# Classify by depth and cover
unique_obs_cells_df <- unique_obs_cells_df %>%
  mutate(
    is_shallow = ifelse(!is.na(depth), depth > -30, NA),
    is_mesophotic = ifelse(!is.na(depth), depth <= -30 & depth >= -60, NA),
    has_coral = max_cover > 0
  )

# Count by categories (excluding NAs for depth stratification)
obs_nonzero <- sum(unique_obs_cells_df$has_coral)
obs_zero <- sum(!unique_obs_cells_df$has_coral)

obs_nonzero_shallow <- sum(unique_obs_cells_df$has_coral & unique_obs_cells_df$is_shallow, na.rm = TRUE)
obs_nonzero_meso <- sum(unique_obs_cells_df$has_coral & unique_obs_cells_df$is_mesophotic, na.rm = TRUE)

obs_zero_shallow <- sum(!unique_obs_cells_df$has_coral & unique_obs_cells_df$is_shallow, na.rm = TRUE)
obs_zero_meso <- sum(!unique_obs_cells_df$has_coral & unique_obs_cells_df$is_mesophotic, na.rm = TRUE)

# Calculate km² and percentages
obs_cells_km2 <- unique_obs_cells * cell_area_km2
obs_cells_pct <- (obs_cells_km2 / total_seafloor_km2) * 100

obs_nonzero_km2 <- obs_nonzero * cell_area_km2
obs_zero_km2 <- obs_zero * cell_area_km2
obs_nonzero_pct <- (obs_nonzero_km2 / total_seafloor_km2) * 100
obs_zero_pct <- (obs_zero_km2 / total_seafloor_km2) * 100

obs_nonzero_shallow_km2 <- obs_nonzero_shallow * cell_area_km2
obs_nonzero_meso_km2 <- obs_nonzero_meso * cell_area_km2
obs_zero_shallow_km2 <- obs_zero_shallow * cell_area_km2
obs_zero_meso_km2 <- obs_zero_meso * cell_area_km2

obs_nonzero_shallow_pct <- (obs_nonzero_shallow_km2 / total_seafloor_km2) * 100
obs_nonzero_meso_pct <- (obs_nonzero_meso_km2 / total_seafloor_km2) * 100
obs_zero_shallow_pct <- (obs_zero_shallow_km2 / total_seafloor_km2) * 100
obs_zero_meso_pct <- (obs_zero_meso_km2 / total_seafloor_km2) * 100

cat(sprintf("\nCells with field observations: %d (%.1f km², %.2f%% of seafloor)\n", 
            unique_obs_cells, obs_cells_km2, obs_cells_pct))
cat(sprintf("  - With non-zero coral cover: %d (%.1f km², %.2f%% of seafloor)\n", 
            obs_nonzero, obs_nonzero_km2, obs_nonzero_pct))
cat(sprintf("    • Shallow (0-29m): %d (%.1f km², %.2f%% of seafloor)\n", 
            obs_nonzero_shallow, obs_nonzero_shallow_km2, obs_nonzero_shallow_pct))
cat(sprintf("    • Mesophotic (30-60m): %d (%.1f km², %.2f%% of seafloor)\n", 
            obs_nonzero_meso, obs_nonzero_meso_km2, obs_nonzero_meso_pct))
cat(sprintf("  - With zero coral cover: %d (%.1f km², %.2f%% of seafloor)\n", 
            obs_zero, obs_zero_km2, obs_zero_pct))
cat(sprintf("    • Shallow (0-29m): %d (%.1f km², %.2f%% of seafloor)\n", 
            obs_zero_shallow, obs_zero_shallow_km2, obs_zero_shallow_pct))
cat(sprintf("    • Mesophotic (30-60m): %d (%.1f km², %.2f%% of seafloor)\n", 
            obs_zero_meso, obs_zero_meso_km2, obs_zero_meso_pct))
cat(sprintf("  - Total PSUs: %d\n", nrow(psu_total_cover)))

cat(sprintf("\nCells with model predictions: %d (%.1f km², %.2f%% of seafloor)\n", 
            total_valid_cells, total_seafloor_km2, 100.0))
cat(sprintf("  - With non-zero predicted coral cover: %d (%.1f km², %.2f%% of seafloor)\n", 
            nonzero_cells, nonzero_km2, nonzero_pct))
cat(sprintf("  - With zero predicted coral cover: %d (%.1f km², %.2f%% of seafloor)\n", 
            zero_cells, zero_km2, zero_pct))

cat("\n✓ Species extent calculations complete\n")
cat("\nScript completed successfully!\n")