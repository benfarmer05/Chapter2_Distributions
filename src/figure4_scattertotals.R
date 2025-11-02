# .rs.restartR(clean = TRUE)
# rm(list=ls())

library(here)
library(terra)
library(dplyr)
library(ggplot2)
library(patchwork)
library(extrafont)

################################## CONTROL SETTINGS ##################################

# Text sizes
plot_title_size <- 10
axis_title_size <- 9
axis_text_size <- 8
annotation_text_size <- 3
panel_label_size <- 5

# Plot margins
plot_margin_pts <- 2

# Panel label position (inside plot) - bottom right
panel_label_x <- 0.95
panel_label_y <- 0.05

# Point sizes
scatter_point_size <- 0.8

################################## Setup ##################################

# extrafont::loadfonts(device = "win", quiet = TRUE)
# 
# # Load combined benthic data
# load(here("output", "all_combined_data.rda"))
# 
# # Load spatial metadata
# cat("Loading spatial metadata...\n")
# spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
# bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
# terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs

################################## Read model outputs ##################################

# cat("Reading model outputs...\n")
# 
# # Get all results files
# results_files <- list.files(here('output/output_maps'), 
#                             pattern = '^results_light_.*\\.rds$', 
#                             full.names = TRUE)
# 
# # Read all results and unwrap rasters
# all_results <- lapply(results_files, function(f) {
#   result <- readRDS(f)
#   result$raster <- terra::unwrap(result$raster)
#   result$raster_raw <- terra::unwrap(result$raster_raw)
#   return(result)
# })
# names(all_results) <- gsub('results_light_|\\.rds', '', basename(results_files))
# 
# cat("Loaded", length(all_results), "species models\n")

################################## Restore CRS and sum rasters ##################################

# cat("Restoring CRS and calculating total cover...\n")
# 
# # Extract all rasters and restore CRS to each
# raster_list <- lapply(all_results, function(x) {
#   r <- x$raster
#   terra::crs(r) <- spatial_metadata$bathy_final$crs
#   return(r)
# })
# 
# # Stack all rasters and sum them
# total_raster <- Reduce(`+`, raster_list)
# 
# # Create richness raster (count of species with cover > 0 at each pixel)
# richness_raster <- Reduce(`+`, lapply(raster_list, function(r) r > 0))
# 
# cat("Total and richness rasters created\n")

################################## Calculate observed total cover by PSU ##################################

cat("Calculating observed total cover...\n")

all_species <- names(all_results)

# Sum observed cover across all species
psu_total_cover <- combined_benthic_data_averaged %>%
  mutate(genus = tolower(sub(" .*", "", spp))) %>%
  filter(genus %in% all_species) %>%
  group_by(PSU) %>%
  summarise(lat = first(lat), 
            lon = first(lon), 
            total_cover = sum(cover, na.rm = TRUE), 
            .groups = 'drop') %>%
  mutate(total_cover_prop = total_cover / 100)

cat("PSU total cover calculated\n")

################################## Calculate observed species richness by PSU ##################################

cat("Calculating observed species richness...\n")

# Count number of species present (cover > 0) at each PSU
# Include sites with richness = 0 (no species present)
psu_richness <- combined_benthic_data_averaged %>%
  mutate(genus = tolower(sub(" .*", "", spp))) %>%
  filter(genus %in% all_species) %>%
  group_by(PSU) %>%
  summarise(lat = first(lat),
            lon = first(lon),
            observed_richness = sum(cover > 0),  # Count species with cover > 0
            .groups = 'drop')

cat("PSU richness calculated\n")

################################## Extract predicted values ##################################

cat("Extracting predicted values at observation locations...\n")

# Create spatial points in lat/lon (WGS84)
obs_points <- vect(psu_total_cover, geom = c("lon", "lat"), crs = "EPSG:4326")

# Transform to match raster CRS (UTM)
obs_points_utm <- project(obs_points, crs(total_raster))

# Extract total cover values
predicted_cover <- extract(total_raster, obs_points_utm)

# Extract richness values
predicted_richness <- extract(richness_raster, obs_points_utm)

# Combine observed and predicted for cover
comparison_cover <- data.frame(
  PSU = psu_total_cover$PSU,
  observed = psu_total_cover$total_cover,
  predicted = predicted_cover[,2] * 100  # Convert to percentage
)

# Combine observed and predicted for richness
comparison_richness <- psu_richness %>%
  left_join(data.frame(PSU = psu_total_cover$PSU, 
                       predicted_richness = predicted_richness[,2]), 
            by = "PSU")

# Remove NAs
comparison_cover <- comparison_cover[complete.cases(comparison_cover), ]
comparison_richness <- comparison_richness[complete.cases(comparison_richness), ]

cat("Valid pairs - Cover:", nrow(comparison_cover), "| Richness:", nrow(comparison_richness), "\n")

################################## Calculate statistics ##################################

# OCCUPIED SITES ONLY (observed > 0 AND predicted > 0)
comparison_cover_occupied <- comparison_cover %>% filter(observed > 0 & predicted > 0)
comparison_richness_occupied <- comparison_richness %>% filter(observed_richness > 0 & predicted_richness > 0)

# Cover statistics (occupied)
cor_test_cover <- cor.test(comparison_cover_occupied$observed, comparison_cover_occupied$predicted)
r2_cover_occ <- cor_test_cover$estimate^2
p_cover_occ <- cor_test_cover$p.value
rmse_cover_occ <- sqrt(mean((comparison_cover_occupied$observed - comparison_cover_occupied$predicted)^2))
mae_cover_occ <- mean(abs(comparison_cover_occupied$observed - comparison_cover_occupied$predicted))

# Richness statistics (occupied)
cor_test_richness <- cor.test(comparison_richness_occupied$observed_richness, comparison_richness_occupied$predicted_richness)
r2_richness_occ <- cor_test_richness$estimate^2
p_richness_occ <- cor_test_richness$p.value
rmse_richness_occ <- sqrt(mean((comparison_richness_occupied$observed_richness - comparison_richness_occupied$predicted_richness)^2))
mae_richness_occ <- mean(abs(comparison_richness_occupied$observed_richness - comparison_richness_occupied$predicted_richness))

# Print R² values and p-values with full precision
cat("\n=== R² Values and P-values (full precision) ===\n")
cat("OCCUPIED SITES ONLY:\n")
cat("  Cover R²:    ", format(r2_cover_occ, digits = 10), "\n")
cat("  Cover p:     ", format(p_cover_occ, digits = 10), "\n")
cat("  Richness R²: ", format(r2_richness_occ, digits = 10), "\n")
cat("  Richness p:  ", format(p_richness_occ, digits = 10), "\n")
cat("  Difference:  ", format(abs(r2_cover_occ - r2_richness_occ), digits = 10), "\n")

################################## Theme ##################################

common_theme <- theme_classic(base_family = "Georgia") +
  theme(
    plot.title = element_text(size = plot_title_size, face = "bold", hjust = 0.5),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    plot.margin = margin(plot_margin_pts, plot_margin_pts, plot_margin_pts, plot_margin_pts, 'pt'),
    legend.position = "right",
    legend.title = element_text(size = axis_text_size),
    legend.text = element_text(size = axis_text_size - 1)
  )

# Colorblind-friendly colors
col_points <- "#E74C3C"  # Red for points

# Helper function to add panel labels
add_panel_label <- function(label, plot_max_x, plot_max_y) {
  annotate("text", 
           x = plot_max_x * 0.95, 
           y = plot_max_y * 0.95,
           label = label,
           hjust = 1, vjust = 1,
           size = panel_label_size,
           fontface = "bold",
           family = "Georgia")
}

################################## Create plots ##################################

cat("\nCreating scatterplots...\n")

# Calculate plot maxima - round up to nearest multiple of 5
plot_max_cover_occ <- ceiling(max(c(comparison_cover_occupied$observed, comparison_cover_occupied$predicted)) / 5) * 5
plot_max_richness_occ <- ceiling(max(c(comparison_richness_occupied$observed_richness, comparison_richness_occupied$predicted_richness)) / 5) * 5

# Create axis breaks (4 values total including 0)
cover_breaks <- seq(0, plot_max_cover_occ, length.out = 4)
richness_breaks <- seq(0, plot_max_richness_occ, length.out = 4)

cat("Cover plot max:", plot_max_cover_occ, "Breaks:", cover_breaks, "\n")
cat("Richness plot max:", plot_max_richness_occ, "Breaks:", richness_breaks, "\n")

################################## SCATTER VERSION ##################################

# Panel A: Total Cover (occupied sites only) - scatter
plot_cover_occ_scatter <- ggplot(comparison_cover_occupied, aes(x = observed, y = predicted)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_jitter(alpha = 0.4, size = scatter_point_size, color = col_points, width = 0.3, height = 0.3) +
  labs(x = "Observed cover (%)",
       y = "Predicted cover (%)") +
  annotate("label", x = plot_max_cover_occ * 0.05, y = plot_max_cover_occ * 0.95, 
           label = sprintf("R² = %.2f\nRMSE = %.2f%%\nMAE = %.2f%%\nN = %d", 
                           r2_cover_occ, rmse_cover_occ, mae_cover_occ, nrow(comparison_cover_occupied)),
           hjust = 0, vjust = 1, size = annotation_text_size, family = "Georgia",
           fill = "white", alpha = 0.6, label.size = 0) +
  add_panel_label("a", plot_max_cover_occ, plot_max_cover_occ) +
  scale_x_continuous(breaks = cover_breaks, limits = c(0, plot_max_cover_occ), expand = c(0, 0)) +
  scale_y_continuous(breaks = cover_breaks, limits = c(0, plot_max_cover_occ), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  common_theme +
  theme(legend.position = "none")

# Panel B: Species Richness (occupied sites only) - scatter
plot_richness_occ_scatter <- ggplot(comparison_richness_occupied, aes(x = observed_richness, y = predicted_richness)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_jitter(alpha = 0.5, size = scatter_point_size, color = col_points, width = 0.4, height = 0.4) +
  labs(x = "Observed richness",
       y = "Predicted richness") +
  annotate("label", x = plot_max_richness_occ * 0.05, y = plot_max_richness_occ * 0.95, 
           label = sprintf("R² = %.2f\nRMSE = %.2f\nMAE = %.2f\nN = %d", 
                           r2_richness_occ, rmse_richness_occ, mae_richness_occ, nrow(comparison_richness_occupied)),
           hjust = 0, vjust = 1, size = annotation_text_size, family = "Georgia",
           fill = "white", alpha = 0.6, label.size = 0) +
  add_panel_label("b", plot_max_richness_occ, plot_max_richness_occ) +
  scale_x_continuous(breaks = richness_breaks, limits = c(0, plot_max_richness_occ), expand = c(0, 0)) +
  scale_y_continuous(breaks = richness_breaks, limits = c(0, plot_max_richness_occ), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  common_theme +
  theme(legend.position = "none")

# Combine scatter plots (2x1 layout)
combined_figure_scatter <- plot_cover_occ_scatter + plot_richness_occ_scatter

################################## BIN2D VERSION ##################################

# Panel A: Total Cover (occupied sites only) - bin2d
# Build plot layer to extract actual binned data
temp_plot_a <- ggplot(comparison_cover_occupied, aes(x = observed, y = predicted)) +
  geom_bin2d(bins = 20)

# Extract the actual bin counts from the plot and filter out zeros
temp_build_a <- ggplot_build(temp_plot_a)
actual_counts_a <- temp_build_a$data[[1]]$count
actual_counts_a_nonzero <- actual_counts_a[actual_counts_a > 0]
count_max_a <- max(actual_counts_a_nonzero, na.rm = TRUE)
count_min_a <- min(actual_counts_a_nonzero, na.rm = TRUE)

# Panel B: Species Richness (occupied sites only) - bin2d with integer bins
# Create integer-aligned bins for discrete count data
richness_occ_breaks_x <- seq(floor(min(comparison_richness_occupied$observed_richness)) - 0.5, 
                             ceiling(max(comparison_richness_occupied$observed_richness)) + 0.5, 
                             by = 1)
richness_occ_breaks_y <- seq(floor(min(comparison_richness_occupied$predicted_richness)) - 0.5, 
                             ceiling(max(comparison_richness_occupied$predicted_richness)) + 0.5, 
                             by = 1)

# Build plot layer to extract actual binned data
temp_plot_b <- ggplot(comparison_richness_occupied, aes(x = observed_richness, y = predicted_richness)) +
  geom_bin2d(breaks = list(x = richness_occ_breaks_x, y = richness_occ_breaks_y))

# Extract the actual bin counts from the plot and filter out zeros
temp_build_b <- ggplot_build(temp_plot_b)
actual_counts_b <- temp_build_b$data[[1]]$count
actual_counts_b_nonzero <- actual_counts_b[actual_counts_b > 0]
count_max_b <- max(actual_counts_b_nonzero, na.rm = TRUE)
count_min_b <- min(actual_counts_b_nonzero, na.rm = TRUE)

# Use the maximum count across both panels for unified legend
count_max_unified <- max(count_max_a, count_max_b)
count_min_unified <- min(count_min_a, count_min_b)

# Create 4 breaks for legend (min, two intermediate, max)
legend_breaks <- round(seq(count_min_unified, count_max_unified, length.out = 4))

cat("Panel A (Cover) - Min count:", count_min_a, "Max count:", count_max_a, "\n")
cat("Panel B (Richness) - Min count:", count_min_b, "Max count:", count_max_b, "\n")
cat("Unified scale - Min:", count_min_unified, "Max:", count_max_unified, "\n")

plot_cover_occ_bin2d <- ggplot(comparison_cover_occupied, aes(x = observed, y = predicted)) +
  geom_bin2d(bins = 20) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "#fee5d9", high = "#a50f15", name = "count",
                      limits = c(count_min_unified, count_max_unified), 
                      breaks = legend_breaks,
                      na.value = "transparent") +
  labs(x = "Observed cover (%)",
       y = "Predicted cover (%)") +
  annotate("label", x = plot_max_cover_occ * 0.05, y = plot_max_cover_occ * 0.95, 
           label = sprintf("R² = %.2f\nRMSE = %.2f%%\nMAE = %.2f%%\nN = %d", 
                           r2_cover_occ, rmse_cover_occ, mae_cover_occ, nrow(comparison_cover_occupied)),
           hjust = 0, vjust = 1, size = annotation_text_size, family = "Georgia",
           fill = "white", alpha = 0.6, label.size = 0) +
  add_panel_label("a", plot_max_cover_occ, plot_max_cover_occ) +
  scale_x_continuous(breaks = cover_breaks, limits = c(0, plot_max_cover_occ), expand = c(0, 0)) +
  scale_y_continuous(breaks = cover_breaks, limits = c(0, plot_max_cover_occ), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  common_theme +
  theme(legend.position = "none")

plot_richness_occ_bin2d <- ggplot(comparison_richness_occupied, aes(x = observed_richness, y = predicted_richness)) +
  geom_bin2d(breaks = list(x = richness_occ_breaks_x, y = richness_occ_breaks_y)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "#fee5d9", high = "#a50f15", name = "count",
                      limits = c(count_min_unified, count_max_unified), 
                      breaks = legend_breaks,
                      na.value = "transparent") +
  labs(x = "Observed richness",
       y = "Predicted richness") +
  annotate("label", x = plot_max_richness_occ * 0.05, y = plot_max_richness_occ * 0.95, 
           label = sprintf("R² = %.2f\nRMSE = %.2f\nMAE = %.2f\nN = %d", 
                           r2_richness_occ, rmse_richness_occ, mae_richness_occ, nrow(comparison_richness_occupied)),
           hjust = 0, vjust = 1, size = annotation_text_size, family = "Georgia",
           fill = "white", alpha = 0.6, label.size = 0) +
  add_panel_label("b", plot_max_richness_occ, plot_max_richness_occ) +
  scale_x_continuous(breaks = richness_breaks, limits = c(0, plot_max_richness_occ), expand = c(0, 0)) +
  scale_y_continuous(breaks = richness_breaks, limits = c(0, plot_max_richness_occ), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  common_theme

# Combine bin2d plots (2x1 layout)
combined_figure_bin2d <- plot_cover_occ_bin2d + plot_richness_occ_bin2d

################################## Print and Save ##################################

print(combined_figure_scatter)
print(combined_figure_bin2d)

output_dir <- here("output", "output_figures_tables")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  filename = file.path(output_dir, "cover_richness_scatterplots_occupied.png"),
  plot = combined_figure_scatter, 
  width = 7.087,
  height = 3.5, 
  dpi = 300, 
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "cover_richness_bin2d_occupied.png"),
  plot = combined_figure_bin2d, 
  width = 7.087,
  height = 3.5, 
  dpi = 300, 
  bg = "white"
)

cat("\n✓ Figures saved:\n")
cat("  - cover_richness_scatterplots_occupied.png (occupied sites scatter)\n")
cat("  - cover_richness_bin2d_occupied.png (occupied sites bin2d with legends)\n")

cat("\n=== OCCUPIED SITES ONLY ===")
cat("\nTotal Cover Performance:\n")
cat("  R²:", round(r2_cover_occ, 3), "\n")
cat("  p-value:", format(p_cover_occ, scientific = TRUE, digits = 3), "\n")
cat("  RMSE:", round(rmse_cover_occ, 2), "%\n")
cat("  MAE:", round(mae_cover_occ, 2), "%\n")
cat("\nSpecies Richness Performance:\n")
cat("  R²:", round(r2_richness_occ, 3), "\n")
cat("  p-value:", format(p_richness_occ, scientific = TRUE, digits = 3), "\n")
cat("  RMSE:", round(rmse_richness_occ, 2), "\n")
cat("  MAE:", round(mae_richness_occ, 2), "\n")