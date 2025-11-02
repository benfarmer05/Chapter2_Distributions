  
# .rs.restartR(clean = TRUE)
rm(list=ls())

library(here)
library(mgcv)
library(ggplot2)
library(patchwork)
library(extrafont)
library(isotone)
library(tidyr)

source(here("src/functions.R"))

################################## CONTROL SETTINGS ##################################

# Text sizes
plot_title_size <- 10
axis_title_size <- 9
axis_text_size <- 8
legend_title_size <- 7
legend_text_size <- 7
legend_key_width <- 0.4  # cm
legend_key_height <- 0.4  # cm
annotation_text_size <- 3
panel_label_size <- 5

# Plot margins
plot_margin_pts <- 2

# Panel label position (inside plot) - bottom right
panel_label_x <- 0.95
panel_label_y <- 0.05

# Line and point sizes
calib_curve_line_size <- 0.8
calib_curve_point_size <- 0.8
density_line_size_raw <- 0.8
density_line_size_observed <- 1.0
density_line_size_calibrated <- 0.8
scatter_point_size <- 0.8

################################## Setup ##################################

extrafont::loadfonts(device = "win", quiet = TRUE)
load(here("output", "all_combined_data.rda"))

# EXPANDED TO FOUR SPECIES
species_list <- c("orbicella", "dendrogyra", "porites", "pseudodiploria")

all_species_models <- lapply(species_list, function(sp) {
  readRDS(here("output", "output_GAMs", paste0(sp, "_models.rds")))
})
names(all_species_models) <- species_list

cat("Loaded models for:", paste(species_list, collapse=", "), "\n")

################################## Calibration Functions ##################################

# Helper function for quantile calibration
calibrate_quantiles <- function(predictions, observations, n_quantiles = 100) {
  probs <- seq(0, 1, length.out = n_quantiles)
  pred_q <- quantile(predictions, probs)
  obs_q <- quantile(observations, probs)
  
  correction_fn <- approxfun(pred_q, obs_q, method = "linear", rule = 2)
  
  list(correction_function = correction_fn)
}

################################## Cross-Validation ##################################

calibration_results <- list()

for(species in species_list) {
  cat("\nProcessing", species, "...")
  
  presence_model <- all_species_models[[species]]$presence_model
  abundance_model <- all_species_models[[species]]$abundance_model
  model_data <- all_species_models[[species]]$model_data
  
  all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  abundance_train <- complete_data[complete_data$cover > 0, ]
  abundance_train$cover_prop <- abundance_train$cover / 100
  
  # In-sample predictions from full model
  train_abundance_pred <- predict(abundance_model, newdata = abundance_train, type = "response")
  
  # Calibration on full dataset
  cal_full <- calibrate_quantiles(train_abundance_pred, abundance_train$cover_prop, n_quantiles = 100)
  calibrated_full <- cal_full$correction_function(train_abundance_pred)
  
  # Cross-validation for out-of-sample predictions
  n <- nrow(abundance_train)
  fold_ids <- sample(rep(1:5, length.out = n))
  cv_results <- data.frame(observed = numeric(n), raw = numeric(n), calibrated = numeric(n))
  
  for(fold in 1:5) {
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    abundance_fit <- gam(formula(abundance_model), data = abundance_train[train_idx, ], family = gaussian())
    test_pred <- predict(abundance_fit, newdata = abundance_train[test_idx, ], type = "response")
    train_pred_fold <- predict(abundance_fit, newdata = abundance_train[train_idx, ], type = "response")
    train_obs_fold <- abundance_train$cover_prop[train_idx]
    
    cal_cv <- calibrate_quantiles(train_pred_fold, train_obs_fold, n_quantiles = 100)
    
    cv_results$observed[test_idx] <- abundance_train$cover_prop[test_idx]
    cv_results$raw[test_idx] <- test_pred
    cv_results$calibrated[test_idx] <- cal_cv$correction_function(test_pred)
  }
  
  calibration_results[[species]] <- list(
    # In-sample (full model)
    raw_pred_full = train_abundance_pred,
    calibrated_full = calibrated_full,
    observed = abundance_train$cover_prop,
    correction_function = cal_full$correction_function,
    # Out-of-sample (CV)
    cv_results = cv_results,
    fold_ids = fold_ids
  )
  cat(" done\n")
}

################################## Theme ##################################

common_theme <- theme_classic(base_family = "Georgia") +
  theme(
    plot.title = element_text(size = plot_title_size, face = "bold.italic", hjust = 0.5),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = legend_text_size),
    legend.key.size = unit(legend_key_width, "cm"),
    plot.margin = margin(plot_margin_pts, plot_margin_pts, plot_margin_pts, plot_margin_pts, 'pt')
  )

# Colorblind-friendly colors
col_raw <- "#E69F00"      # Orange
col_obs <- "#0072B2"      # Blue
col_calib <- "#009E73"    # Green

# Helper function to add panel labels inside plots (bottom right)
add_panel_label <- function(label, plot_max_x, plot_max_y) {
  annotate("text", 
           x = plot_max_x * panel_label_x, 
           y = plot_max_y * panel_label_y,
           label = label,
           hjust = 1, vjust = 0,
           size = panel_label_size,
           fontface = "bold",
           family = "Georgia")
}

################################################################################
##                              FIGURE 1                                      ##
################################################################################

cat("\nCreating figure...\n")

# Panel label counter
panel_labels <- letters[1:20]
label_idx <- 1

# Row 1: Calibration curves
calib_plots <- list()
for(i in seq_along(species_list)) {
  species <- species_list[i]
  data <- calibration_results[[species]]
  pred_grid <- seq(min(data$raw_pred), max(data$raw_pred), length.out = 200)
  curve_df <- data.frame(raw = pred_grid, calibrated = data$correction_function(pred_grid))
  plot_max <- max(c(curve_df$raw, curve_df$calibrated))
  
  p <- ggplot(curve_df, aes(x = raw, y = calibrated)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_line(color = col_calib, linewidth = calib_curve_line_size) +
    labs(title = stringr::str_to_title(species),
         x = "Raw",
         y = "Calibrated") +
    coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
    add_panel_label(panel_labels[label_idx], plot_max, plot_max) +
    common_theme + 
    theme(legend.position = "none")
  
  label_idx <- label_idx + 1
  calib_plots[[species]] <- p
}

# Row 2: Distributions (Raw, Observed, Calibrated)
dist_plots <- list()
for(i in seq_along(species_list)) {
  species <- species_list[i]
  data <- calibration_results[[species]]
  dist_df <- data.frame(
    value = c(data$raw_pred_full, data$observed, data$calibrated_full),
    type = rep(c("Raw", "Observed", "Calibrated"), each = length(data$observed))
  )
  dist_df$type <- factor(dist_df$type, levels = c("Raw", "Observed", "Calibrated"))
  
  # Get plot limits for label positioning
  dens_data <- density(dist_df$value)
  x_max <- max(dist_df$value)
  y_max <- max(dens_data$y)
  
  p <- ggplot(dist_df, aes(x = value, color = type, fill = type, linetype = type, linewidth = type)) +
    geom_density(alpha = 0.15) +
    scale_color_manual(values = c("Raw" = col_raw, "Observed" = col_obs, "Calibrated" = col_calib)) +
    scale_fill_manual(values = c("Raw" = col_raw, "Observed" = col_obs, "Calibrated" = col_calib)) +
    scale_linetype_manual(values = c("Raw" = "solid", "Observed" = "dashed", "Calibrated" = "solid")) +
    scale_linewidth_manual(values = c("Raw" = density_line_size_raw, "Observed" = density_line_size_observed, "Calibrated" = density_line_size_calibrated)) +
    labs(title = NULL,
         x = "Observed cover",
         y = "Density") +
    add_panel_label(panel_labels[label_idx], x_max, y_max) +
    common_theme +
    theme(aspect.ratio = 1)
  
  label_idx <- label_idx + 1
  dist_plots[[species]] <- p
}

# Row 3: Raw scatter (in-sample)
scatter_raw <- list()
for(i in seq_along(species_list)) {
  species <- species_list[i]
  data <- calibration_results[[species]]
  r2_raw <- cor(data$raw_pred_full, data$observed)^2
  plot_max <- max(c(data$raw_pred_full, data$observed))
  
  p <- ggplot(data.frame(predicted = data$raw_pred_full, observed = data$observed), aes(x = predicted, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_point(alpha = 0.4, size = scatter_point_size, color = col_raw) +
    labs(title = NULL,
         x = "Predicted",
         y = "Observed") +
    annotate("label", x = plot_max * 0.05, y = plot_max * 0.95, 
             label = sprintf("R² = %.2f\nN = %d", r2_raw, length(data$observed)),
             hjust = 0, vjust = 1, size = annotation_text_size, family = "Georgia",
             fill = "white", alpha = 0.6, label.size = 0) +
    add_panel_label(panel_labels[label_idx], plot_max, plot_max) +
    coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
    common_theme + 
    theme(legend.position = "none")
  
  label_idx <- label_idx + 1
  scatter_raw[[species]] <- p
}

# Row 4: Calibrated scatter (in-sample)
scatter_calib <- list()
for(i in seq_along(species_list)) {
  species <- species_list[i]
  data <- calibration_results[[species]]
  r2_calib <- cor(data$calibrated_full, data$observed)^2
  plot_max <- max(c(data$calibrated_full, data$observed))
  
  p <- ggplot(data.frame(predicted = data$calibrated_full, observed = data$observed), aes(x = predicted, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_point(alpha = 0.4, size = scatter_point_size, color = col_calib) +
    labs(title = NULL,
         x = "Predicted",
         y = "Observed") +
    annotate("label", x = plot_max * 0.05, y = plot_max * 0.95, 
             label = sprintf("R² = %.2f\nN = %d", r2_calib, length(data$observed)),
             hjust = 0, vjust = 1, size = annotation_text_size, family = "Georgia",
             fill = "white", alpha = 0.6, label.size = 0) +
    add_panel_label(panel_labels[label_idx], plot_max, plot_max) +
    coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
    common_theme + 
    theme(legend.position = "none")
  
  label_idx <- label_idx + 1
  scatter_calib[[species]] <- p
}

# Row 5: Calibrated scatter (out-of-sample)
scatter_cv <- list()
for(i in seq_along(species_list)) {
  species <- species_list[i]
  data <- calibration_results[[species]]
  chosen_fold <- sample(1:5, 1)
  test_idx <- which(data$fold_ids == chosen_fold)
  test_df <- data.frame(
    predicted = data$cv_results$calibrated[test_idx],
    observed = data$cv_results$observed[test_idx]
  )
  plot_max <- max(c(test_df$predicted, test_df$observed))
  r2_cv <- if(nrow(test_df) > 1) cor(test_df$predicted, test_df$observed)^2 else NA
  
  p <- ggplot(test_df, aes(x = predicted, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_point(alpha = 0.6, size = scatter_point_size, color = col_calib) +
    labs(title = NULL,
         x = "Predicted",
         y = "Observed") +
    annotate("label", x = plot_max * 0.05, y = plot_max * 0.95,
             label = sprintf("R² = %s\nN = %d", 
                             ifelse(is.na(r2_cv), "NA", sprintf("%.2f", r2_cv)), nrow(test_df)),
             hjust = 0, vjust = 1, size = annotation_text_size, family = "Georgia",
             fill = "white", alpha = 0.6, label.size = 0) +
    add_panel_label(panel_labels[label_idx], plot_max, plot_max) +
    coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
    common_theme + 
    theme(legend.position = "none")
  
  label_idx <- label_idx + 1
  scatter_cv[[species]] <- p
}

# Combine all rows with proper legend handling
figure1 <- (wrap_plots(calib_plots, ncol = 4) /
              wrap_plots(dist_plots, ncol = 4) /
              wrap_plots(scatter_raw, ncol = 4) /
              wrap_plots(scatter_calib, ncol = 4) /
              wrap_plots(scatter_cv, ncol = 4)) +
  plot_layout(heights = c(1, 1, 1, 1, 1), guides = "collect") &
  theme(legend.position = "bottom")

print(figure1)

################################## Save ##################################

output_dir <- here("output", "output_figures_tables")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  filename = file.path(output_dir, "calibration_figure.png"),
  plot = figure1, 
  width = 7.087,
  height = 8, 
  dpi = 300, 
  bg = "white"
)

cat("\n✓ Figure saved: calibration_figure.png\n")


# # VERSION WITH MORE TESTING; LOOKED AT OTHER TYPES OF CORRECTION
#   #
#   # .rs.restartR(clean = TRUE)
#   rm(list=ls())
#   
#   library(here)
#   library(mgcv)
#   library(ggplot2)
#   library(patchwork)
#   library(extrafont)
#   library(isotone)
#   library(tidyr)
#   
#   source(here("src/functions.R"))
#   
#   ################################## Setup ##################################
#   
#   extrafont::loadfonts(device = "win", quiet = TRUE)
#   load(here("output", "all_combined_data.rda"))
#   
#   species_list <- c("orbicella", "pseudodiploria", "porites")
#   
#   all_species_models <- lapply(species_list, function(sp) {
#     file_path <- here("output", "output_GAMs", paste0(sp, "_models.rds"))
#     if(file.exists(file_path)) {
#       readRDS(file_path)
#     } else {
#       warning(paste("Model file not found for:", sp))
#       NULL
#     }
#   })
#   names(all_species_models) <- species_list
#   
#   cat("Loaded models for:\n")
#   for(sp in species_list) {
#     if(!is.null(all_species_models[[sp]])) {
#       cat("  ✓", sp, "\n")
#     }
#   }
#   
#   ################################## Calibration Functions ##################################
#   
#   calibrate_quantiles <- function(predictions, observations, n_quantiles = 100) {
#     probs <- seq(0, 1, length.out = n_quantiles)
#     pred_q <- quantile(predictions, probs)
#     obs_q <- quantile(observations, probs)
#     correction_fn <- approxfun(pred_q, obs_q, method = "linear", rule = 2)
#     list(correction_function = correction_fn,
#          pred_quantiles = pred_q,
#          obs_quantiles = obs_q)
#   }
#   
#   calibrate_quantiles_smooth <- function(predictions, observations, n_quantiles = 20) {
#     probs <- seq(0, 1, length.out = n_quantiles)
#     pred_q <- quantile(predictions, probs)
#     obs_q <- quantile(observations, probs)
#     smooth_fit <- smooth.spline(pred_q, obs_q, spar = 0.5)
#     correction_fn <- function(x) predict(smooth_fit, x)$y
#     list(correction_function = correction_fn,
#          pred_quantiles = pred_q,
#          obs_quantiles = obs_q)
#   }
#   
#   calibrate_linear <- function(predictions, observations) {
#     fit <- lm(observations ~ predictions)
#     correction_fn <- function(x) predict(fit, newdata = data.frame(predictions = x))
#     list(correction_function = correction_fn,
#          intercept = coef(fit)[1],
#          slope = coef(fit)[2])
#   }
#   
#   calibrate_isotonic <- function(predictions, observations) {
#     ord <- order(predictions)
#     pred_sorted <- predictions[ord]
#     obs_sorted <- observations[ord]
#     iso_fit <- gpava(pred_sorted, obs_sorted)
#     correction_fn <- approxfun(pred_sorted, iso_fit$x, method = "linear", rule = 2)
#     list(correction_function = correction_fn,
#          fitted_values = iso_fit$x)
#   }
#   
#   calibrate_beta <- function(predictions, observations) {
#     epsilon <- 1e-6
#     pred_bounded <- pmax(pmin(predictions, 1 - epsilon), epsilon)
#     obs_bounded <- pmax(pmin(observations, 1 - epsilon), epsilon)
#     logit_pred <- log(pred_bounded / (1 - pred_bounded))
#     logit_obs <- log(obs_bounded / (1 - obs_bounded))
#     fit <- lm(logit_obs ~ logit_pred)
#     correction_fn <- function(x) {
#       x_bounded <- pmax(pmin(x, 1 - epsilon), epsilon)
#       logit_x <- log(x_bounded / (1 - x_bounded))
#       logit_calib <- predict(fit, newdata = data.frame(logit_pred = logit_x))
#       exp(logit_calib) / (1 + exp(logit_calib))
#     }
#     list(correction_function = correction_fn,
#          intercept = coef(fit)[1],
#          slope = coef(fit)[2])
#   }
#   
#   ################################## Cross-Validated Comparison ##################################
#   
#   calibration_results <- list()
#   
#   for(species in species_list) {
#     
#     cat("\n\n========================================\n")
#     cat("Processing", toupper(species), "\n")
#     cat("========================================\n")
#     
#     presence_model <- all_species_models[[species]]$presence_model
#     abundance_model <- all_species_models[[species]]$abundance_model
#     model_data <- all_species_models[[species]]$model_data
#     
#     all_vars <- unique(c(all.vars(formula(presence_model)), 
#                          all.vars(formula(abundance_model))))
#     complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
#     
#     abundance_train <- complete_data[complete_data$cover > 0, ]
#     abundance_train$cover_prop <- abundance_train$cover / 100
#     
#     train_abundance_pred <- predict(abundance_model, 
#                                     newdata = abundance_train, 
#                                     type = "response")
#     
#     cat("\nData summary:\n")
#     cat("  N observations:", length(train_abundance_pred), "\n")
#     cat("  Raw pred range:", round(range(train_abundance_pred), 4), "\n")
#     cat("  Observed range:", round(range(abundance_train$cover_prop), 4), "\n")
#     
#     # make CV randomized each run: DO NOT set a fixed seed here
#     n <- nrow(abundance_train)
#     # ensure approx 5-fold (each ~20%) by sampling rep(1:5)
#     fold_ids <- sample(rep(1:5, length.out = n))
#     
#     cv_results <- data.frame(
#       observed = numeric(n),
#       raw = numeric(n),
#       quantile = numeric(n),
#       quantile_smooth = numeric(n),
#       linear = numeric(n),
#       isotonic = numeric(n),
#       beta = numeric(n)
#     )
#     
#     cat("\nCross-validation:\n")
#     for(fold in 1:5) {
#       cat("  Fold", fold, "...")
#       
#       test_idx <- which(fold_ids == fold)
#       train_idx <- which(fold_ids != fold)
#       
#       abundance_fit <- gam(formula(abundance_model), 
#                            data = abundance_train[train_idx, ], 
#                            family = gaussian())
#       
#       test_pred <- predict(abundance_fit, 
#                            newdata = abundance_train[test_idx, ], 
#                            type = "response")
#       test_obs <- abundance_train$cover_prop[test_idx]
#       
#       train_pred_fold <- predict(abundance_fit, 
#                                  newdata = abundance_train[train_idx, ], 
#                                  type = "response")
#       train_obs_fold <- abundance_train$cover_prop[train_idx]
#       
#       cal_quantile <- calibrate_quantiles(train_pred_fold, train_obs_fold)
#       cal_quantile_smooth <- calibrate_quantiles_smooth(train_pred_fold, train_obs_fold)
#       cal_linear <- calibrate_linear(train_pred_fold, train_obs_fold)
#       cal_isotonic <- calibrate_isotonic(train_pred_fold, train_obs_fold)
#       cal_beta <- calibrate_beta(train_pred_fold, train_obs_fold)
#       
#       cv_results$observed[test_idx] <- test_obs
#       cv_results$raw[test_idx] <- test_pred
#       cv_results$quantile[test_idx] <- cal_quantile$correction_function(test_pred)
#       cv_results$quantile_smooth[test_idx] <- cal_quantile_smooth$correction_function(test_pred)
#       cv_results$linear[test_idx] <- cal_linear$correction_function(test_pred)
#       cv_results$isotonic[test_idx] <- cal_isotonic$correction_function(test_pred)
#       cv_results$beta[test_idx] <- cal_beta$correction_function(test_pred)
#       
#       cat(" done\n")
#     }
#     
#     calc_metrics <- function(pred, obs) {
#       rmse <- sqrt(mean((pred - obs)^2))
#       mae <- mean(abs(pred - obs))
#       r2 <- cor(pred, obs)^2
#       bias <- mean(pred - obs)
#       c(RMSE = rmse, MAE = mae, R2 = r2, Bias = bias)
#     }
#     
#     metrics <- data.frame(
#       Raw = calc_metrics(cv_results$raw, cv_results$observed),
#       Quantile = calc_metrics(cv_results$quantile, cv_results$observed),
#       Quantile_Smooth = calc_metrics(cv_results$quantile_smooth, cv_results$observed),
#       Linear = calc_metrics(cv_results$linear, cv_results$observed),
#       Isotonic = calc_metrics(cv_results$isotonic, cv_results$observed),
#       Beta = calc_metrics(cv_results$beta, cv_results$observed)
#     )
#     
#     cat("\nCross-validated metrics:\n")
#     print(round(metrics, 4))
#     
#     cal_quantile_full <- calibrate_quantiles(train_abundance_pred, abundance_train$cover_prop)
#     cal_quantile_smooth_full <- calibrate_quantiles_smooth(train_abundance_pred, abundance_train$cover_prop)
#     cal_linear_full <- calibrate_linear(train_abundance_pred, abundance_train$cover_prop)
#     cal_isotonic_full <- calibrate_isotonic(train_abundance_pred, abundance_train$cover_prop)
#     cal_beta_full <- calibrate_beta(train_abundance_pred, abundance_train$cover_prop)
#     
#     # save fold_ids and cv_results so we can re-use the exact CV test sets later
#     calibration_results[[species]] <- list(
#       cv_results = cv_results,
#       metrics = metrics,
#       raw_pred = train_abundance_pred,
#       observed = abundance_train$cover_prop,
#       quantile_full = cal_quantile_full,
#       quantile_smooth_full = cal_quantile_smooth_full,
#       linear_full = cal_linear_full,
#       isotonic_full = cal_isotonic_full,
#       beta_full = cal_beta_full,
#       fold_ids = fold_ids  # store the fold assignment for this species
#     )
#   }
#   
#   ################################## Common Theme ##################################
#   
#   common_theme <- theme_classic(base_family = "Georgia") +
#     theme(
#       plot.title = element_text(size = 10, face = "bold.italic", hjust = 0.5),
#       axis.title = element_text(size = 9),
#       axis.text = element_text(size = 8),
#       legend.position = "right",
#       legend.title = element_blank(),
#       legend.text = element_text(size = 7),
#       legend.key.size = unit(0.4, "cm"),
#       plot.margin = margin(5, 5, 5, 5, 'pt')
#     )
#   
#   ################################################################################
#   ##                       FIGURE 1: ALL METHODS COMPARISON                     ##
#   ################################################################################
#   
#   cat("\n\nCreating Figure 1: All methods comparison...\n")
#   
#   calib_curve_plots <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     pred_grid <- seq(min(data$raw_pred), max(data$raw_pred), length.out = 200)
#     curve_df <- data.frame(
#       raw = pred_grid,
#       Quantile = data$quantile_full$correction_function(pred_grid),
#       Quantile_Smooth = data$quantile_smooth_full$correction_function(pred_grid),
#       Linear = data$linear_full$correction_function(pred_grid),
#       Isotonic = data$isotonic_full$correction_function(pred_grid),
#       Beta = data$beta_full$correction_function(pred_grid)
#     )
#     curve_long <- pivot_longer(curve_df, 
#                                cols = c(Quantile, Quantile_Smooth, Linear, Isotonic, Beta),
#                                names_to = "Method", values_to = "calibrated")
#     plot_max <- max(c(curve_df$raw, curve_df$Quantile, curve_df$Quantile_Smooth, 
#                       curve_df$Linear, curve_df$Isotonic, curve_df$Beta))
#     p <- ggplot(curve_long, aes(x = raw, y = calibrated, color = Method)) +
#       geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
#                   color = "gray50", linewidth = 0.5) +
#       geom_line(linewidth = 1) +
#       scale_color_manual(values = c("Quantile" = "#d6604d", "Quantile_Smooth" = "#f4a582",
#                                     "Linear" = "#4393c3", "Isotonic" = "#5aae61",
#                                     "Beta" = "#9970ab")) +
#       labs(title = stringr::str_to_title(species), x = "Raw Predicted Cover",
#            y = "Calibrated Cover") +
#       coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
#       common_theme
#     calib_curve_plots[[species]] <- p
#   }
#   calib_curve_panel <- wrap_plots(calib_curve_plots, ncol = 3, guides = "collect")
#   
#   metrics_plots <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     metrics_long <- data.frame(
#       Metric = rep(rownames(data$metrics), ncol(data$metrics)),
#       Method = rep(colnames(data$metrics), each = nrow(data$metrics)),
#       Value = as.vector(as.matrix(data$metrics))
#     )
#     p_rmse <- ggplot(metrics_long[metrics_long$Metric == "RMSE", ], 
#                      aes(x = Method, y = Value, fill = Method)) +
#       geom_col() +
#       scale_fill_manual(values = c("Raw" = "#999999", "Quantile" = "#d6604d",
#                                    "Quantile_Smooth" = "#f4a582", "Linear" = "#4393c3",
#                                    "Isotonic" = "#5aae61", "Beta" = "#9970ab")) +
#       labs(title = stringr::str_to_title(species), y = "RMSE") +
#       theme_classic(base_family = "Georgia") +
#       theme(axis.title.x = element_blank(),
#             axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
#             legend.position = "none",
#             plot.title = element_text(size = 10, face = "bold.italic", hjust = 0.5))
#     metrics_plots[[species]] <- p_rmse
#   }
#   metrics_panel <- wrap_plots(metrics_plots, ncol = 3)
#   
#   scatter_plots <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     cv <- data$cv_results
#     scatter_long <- data.frame(
#       observed = rep(cv$observed, 6),
#       predicted = c(cv$raw, cv$quantile, cv$quantile_smooth, cv$linear, cv$isotonic, cv$beta),
#       method = rep(c("Raw", "Quantile", "Quantile_Smooth", "Linear", "Isotonic", "Beta"), 
#                    each = nrow(cv))
#     )
#     scatter_long$method <- factor(scatter_long$method, 
#                                   levels = c("Raw", "Quantile", "Quantile_Smooth", "Linear", "Isotonic", "Beta"))
#     plot_max <- max(c(scatter_long$predicted, scatter_long$observed))
#     p <- ggplot(scatter_long, aes(x = predicted, y = observed, color = method)) +
#       geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
#                   color = "gray50", linewidth = 0.5) +
#       geom_point(alpha = 0.3, size = 0.8) +
#       scale_color_manual(values = c("Raw" = "#999999", "Quantile" = "#d6604d",
#                                     "Quantile_Smooth" = "#f4a582", "Linear" = "#4393c3",
#                                     "Isotonic" = "#5aae61", "Beta" = "#9970ab")) +
#       facet_wrap(~method, ncol = 6) +
#       labs(title = stringr::str_to_title(species), x = "Predicted Cover",
#            y = "Observed Cover") +
#       coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
#       theme_classic(base_family = "Georgia") +
#       theme(legend.position = "none",
#             plot.title = element_text(size = 10, face = "bold.italic", hjust = 0.5),
#             strip.text = element_text(size = 8, face = "bold"),
#             strip.background = element_rect(fill = "gray90", color = NA),
#             axis.title = element_text(size = 8), axis.text = element_text(size = 7))
#     # add sample-size annotation for the aggregated scatter (use total CV n)
#     p <- p + annotate("text", x = plot_max * 0.05, y = plot_max * 0.85,
#                       label = paste0("n = ", nrow(cv)), hjust = 0, vjust = 1,
#                       size = 3, family = "Georgia")
#     scatter_plots[[species]] <- p
#   }
#   scatter_panel <- wrap_plots(scatter_plots, ncol = 1)
#   
#   figure1 <- calib_curve_panel / metrics_panel / scatter_panel +
#     plot_layout(heights = c(1, 0.8, 2.5)) +
#     plot_annotation(tag_levels = 'a',
#                     theme = theme(plot.tag = element_text(size = 10, face = "bold", family = "Georgia")))
#   
#   ################################################################################
#   ##              FIGURE 2: SMOOTHED QUANTILE DETAILED VIEW                     ##
#   ################################################################################
#   
#   cat("\nCreating Figure 2: Smoothed quantile detailed view...\n")
#   
#   smooth_calib_plots <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     pred_grid <- seq(min(data$raw_pred), max(data$raw_pred), length.out = 200)
#     curve_df <- data.frame(raw = pred_grid,
#                            calibrated = data$quantile_smooth_full$correction_function(pred_grid))
#     plot_max <- max(c(curve_df$raw, curve_df$calibrated))
#     p <- ggplot(curve_df, aes(x = raw, y = calibrated)) +
#       geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
#                   color = "gray50", linewidth = 0.5) +
#       geom_line(color = "#f4a582", linewidth = 1.2) +
#       geom_point(color = "#f4a582", size = 1.5, alpha = 0.6) +
#       labs(title = stringr::str_to_title(species),
#            x = "Raw Predicted Cover (proportion)", y = "Calibrated Cover (proportion)") +
#       coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
#       common_theme + theme(legend.position = "none")
#     smooth_calib_plots[[species]] <- p
#   }
#   
#   raw_vs_obs_plots <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     cv <- data$cv_results
#     dist_df <- data.frame(value = c(cv$raw, cv$observed),
#                           type = rep(c("Raw Predicted", "Observed"), each = nrow(cv)))
#     dist_df$type <- factor(dist_df$type, levels = c("Raw Predicted", "Observed"))
#     p <- ggplot(dist_df, aes(x = value, color = type, fill = type)) +
#       geom_density(alpha = 0.2, linewidth = 1.2) +
#       scale_fill_manual(values = c("Raw Predicted" = "#d6604d", "Observed" = "#2d2d2d")) +
#       scale_color_manual(values = c("Raw Predicted" = "#d6604d", "Observed" = "#2d2d2d")) +
#       labs(title = stringr::str_to_title(species), x = "Cover (proportion)", y = "Density") +
#       common_theme
#     raw_vs_obs_plots[[species]] <- p
#   }
#   
#   raw_vs_calib_smooth_plots <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     cv <- data$cv_results
#     dist_df <- data.frame(value = c(cv$raw, cv$quantile_smooth),
#                           type = rep(c("Raw Predicted", "Calibrated"), each = nrow(cv)))
#     dist_df$type <- factor(dist_df$type, levels = c("Raw Predicted", "Calibrated"))
#     p <- ggplot(dist_df, aes(x = value, color = type, fill = type)) +
#       geom_density(alpha = 0.2, linewidth = 1.2) +
#       scale_fill_manual(values = c("Raw Predicted" = "#d6604d", "Calibrated" = "#f4a582")) +
#       scale_color_manual(values = c("Raw Predicted" = "#d6604d", "Calibrated" = "#f4a582")) +
#       labs(title = stringr::str_to_title(species), x = "Cover (proportion)", y = "Density") +
#       common_theme
#     raw_vs_calib_smooth_plots[[species]] <- p
#   }
#   
#   scatter_raw_plots <- list()
#   scatter_calib_smooth_plots <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     cv <- data$cv_results
#     r2_raw <- cor(cv$raw, cv$observed)^2
#     r2_calib <- cor(cv$quantile_smooth, cv$observed)^2
#     plot_max <- max(c(cv$raw, cv$quantile_smooth, cv$observed))
#     
#     p_raw <- ggplot(data.frame(predicted = cv$raw, observed = cv$observed), 
#                     aes(x = predicted, y = observed)) +
#       geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
#                   color = "gray50", linewidth = 0.5) +
#       geom_point(alpha = 0.4, size = 1.2, color = "#d6604d") +
#       labs(title = stringr::str_to_title(species),
#            x = "Predicted Cover (proportion)", y = "Observed Cover (proportion)") +
#       annotate("text", x = plot_max * 0.05, y = plot_max * 0.95, 
#                label = sprintf("R² = %.3f", r2_raw),
#                hjust = 0, vjust = 1, size = 3, family = "Georgia") +
#       # sample size annotation for this scatter
#       annotate("text", x = plot_max * 0.05, y = plot_max * 0.85,
#                label = paste0("n = ", nrow(cv)), hjust = 0, vjust = 1,
#                size = 3, family = "Georgia") +
#       coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
#       common_theme + theme(legend.position = "none")
#     scatter_raw_plots[[species]] <- p_raw
#     
#     p_calib <- ggplot(data.frame(predicted = cv$quantile_smooth, observed = cv$observed), 
#                       aes(x = predicted, y = observed)) +
#       geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
#                   color = "gray50", linewidth = 0.5) +
#       geom_point(alpha = 0.4, size = 1.2, color = "#f4a582") +
#       labs(title = stringr::str_to_title(species),
#            x = "Predicted Cover (proportion)", y = "Observed Cover (proportion)") +
#       annotate("text", x = plot_max * 0.05, y = plot_max * 0.95, 
#                label = sprintf("R² = %.3f", r2_calib),
#                hjust = 0, vjust = 1, size = 3, family = "Georgia") +
#       # sample size annotation
#       annotate("text", x = plot_max * 0.05, y = plot_max * 0.85,
#                label = paste0("n = ", nrow(cv)), hjust = 0, vjust = 1,
#                size = 3, family = "Georgia") +
#       coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
#       common_theme + theme(legend.position = "none")
#     scatter_calib_smooth_plots[[species]] <- p_calib
#   }
#   
#   # create out-of-sample bottom row by selecting one of the 5 folds produced earlier
#   outsample_scatter_smooth <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     # reuse the stored fold assignment so we pick one of the actual CV folds
#     fold_ids <- data$fold_ids
#     chosen_fold <- sample(1:5, 1)
#     test_idx <- which(fold_ids == chosen_fold)
#     test_df <- data.frame(predicted = data$cv_results$quantile_smooth[test_idx],
#                           observed = data$cv_results$observed[test_idx])
#     plot_max <- max(c(test_df$predicted, test_df$observed))
#     r2_out <- if(nrow(test_df) > 1) cor(test_df$predicted, test_df$observed)^2 else NA
#     p_out <- ggplot(test_df, aes(x = predicted, y = observed)) +
#       geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
#       geom_point(alpha = 0.6, size = 1.4, color = "#f4a582") +
#       labs(title = stringr::str_to_title(species), x = "Predicted Cover (proportion)", y = "Observed Cover (proportion)") +
#       annotate("text", x = plot_max * 0.05, y = plot_max * 0.95,
#                label = sprintf("R² = %s", ifelse(is.na(r2_out), "NA", sprintf("%.3f", r2_out))),
#                hjust = 0, vjust = 1, size = 3, family = "Georgia") +
#       annotate("text", x = plot_max * 0.05, y = plot_max * 0.85,
#                label = paste0("n = ", nrow(test_df)), hjust = 0, vjust = 1,
#                size = 3, family = "Georgia") +
#       coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
#       common_theme + theme(legend.position = "none")
#     outsample_scatter_smooth[[species]] <- p_out
#   }
#   
#   figure2 <- wrap_plots(smooth_calib_plots, ncol = 3) /
#     wrap_plots(raw_vs_obs_plots, ncol = 3, guides = "collect") /
#     wrap_plots(raw_vs_calib_smooth_plots, ncol = 3, guides = "collect") /
#     wrap_plots(scatter_raw_plots, ncol = 3) /
#     wrap_plots(scatter_calib_smooth_plots, ncol = 3) /
#     wrap_plots(outsample_scatter_smooth, ncol = 3) +
#     plot_layout(heights = c(1, 1, 1, 1, 1, 1)) +
#     plot_annotation(tag_levels = 'a',
#                     theme = theme(plot.tag = element_text(size = 10, face = "bold", family = "Georgia")))
#   
#   ################################################################################
#   ##              FIGURE 3: REGULAR QUANTILE DETAILED VIEW                      ##
#   ################################################################################
#   
#   cat("\nCreating Figure 3: Regular quantile detailed view...\n")
#   
#   regular_calib_plots <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     pred_grid <- seq(min(data$raw_pred), max(data$raw_pred), length.out = 200)
#     curve_df <- data.frame(raw = pred_grid,
#                            calibrated = data$quantile_full$correction_function(pred_grid))
#     plot_max <- max(c(curve_df$raw, curve_df$calibrated))
#     p <- ggplot(curve_df, aes(x = raw, y = calibrated)) +
#       geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
#                   color = "gray50", linewidth = 0.5) +
#       geom_line(color = "#d6604d", linewidth = 1.2) +
#       geom_point(color = "#d6604d", size = 1.5, alpha = 0.6) +
#       labs(title = stringr::str_to_title(species),
#            x = "Raw Predicted Cover (proportion)", y = "Calibrated Cover (proportion)") +
#       coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
#       common_theme + theme(legend.position = "none")
#     regular_calib_plots[[species]] <- p
#   }
#   
#   raw_vs_calib_regular_plots <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     cv <- data$cv_results
#     dist_df <- data.frame(value = c(cv$raw, cv$quantile),
#                           type = rep(c("Raw Predicted", "Calibrated"), each = nrow(cv)))
#     dist_df$type <- factor(dist_df$type, levels = c("Raw Predicted", "Calibrated"))
#     p <- ggplot(dist_df, aes(x = value, color = type, fill = type)) +
#       geom_density(alpha = 0.2, linewidth = 1.2) +
#       scale_fill_manual(values = c("Raw Predicted" = "#d6604d", "Calibrated" = "#d6604d")) +
#       scale_color_manual(values = c("Raw Predicted" = "#d6604d", "Calibrated" = "#d6604d")) +
#       labs(title = stringr::str_to_title(species), x = "Cover (proportion)", y = "Density") +
#       common_theme
#     raw_vs_calib_regular_plots[[species]] <- p
#   }
#   
#   scatter_calib_regular_plots <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     cv <- data$cv_results
#     r2_calib <- cor(cv$quantile, cv$observed)^2
#     plot_max <- max(c(cv$raw, cv$quantile, cv$observed))
#     
#     p_calib <- ggplot(data.frame(predicted = cv$quantile, observed = cv$observed), 
#                       aes(x = predicted, y = observed)) +
#       geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
#                   color = "gray50", linewidth = 0.5) +
#       geom_point(alpha = 0.4, size = 1.2, color = "#d6604d") +
#       labs(title = stringr::str_to_title(species),
#            x = "Predicted Cover (proportion)", y = "Observed Cover (proportion)") +
#       annotate("text", x = plot_max * 0.05, y = plot_max * 0.95, 
#                label = sprintf("R² = %.3f", r2_calib),
#                hjust = 0, vjust = 1, size = 3, family = "Georgia") +
#       annotate("text", x = plot_max * 0.05, y = plot_max * 0.85,
#                label = paste0("n = ", nrow(cv)), hjust = 0, vjust = 1,
#                size = 3, family = "Georgia") +
#       coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
#       common_theme + theme(legend.position = "none")
#     scatter_calib_regular_plots[[species]] <- p_calib
#   }
#   
#   # out-of-sample bottom row for regular quantile: reuse chosen fold per species
#   outsample_scatter_regular <- list()
#   for(species in species_list) {
#     data <- calibration_results[[species]]
#     fold_ids <- data$fold_ids
#     chosen_fold <- sample(1:5, 1)
#     test_idx <- which(fold_ids == chosen_fold)
#     test_df <- data.frame(predicted = data$cv_results$quantile[test_idx],
#                           observed = data$cv_results$observed[test_idx])
#     plot_max <- max(c(test_df$predicted, test_df$observed))
#     r2_out <- if(nrow(test_df) > 1) cor(test_df$predicted, test_df$observed)^2 else NA
#     p_out <- ggplot(test_df, aes(x = predicted, y = observed)) +
#       geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
#       geom_point(alpha = 0.6, size = 1.4, color = "#d6604d") +
#       labs(title = stringr::str_to_title(species), x = "Predicted Cover (proportion)", y = "Observed Cover (proportion)") +
#       annotate("text", x = plot_max * 0.05, y = plot_max * 0.95,
#                label = sprintf("R² = %s", ifelse(is.na(r2_out), "NA", sprintf("%.3f", r2_out))),
#                hjust = 0, vjust = 1, size = 3, family = "Georgia") +
#       annotate("text", x = plot_max * 0.05, y = plot_max * 0.85,
#                label = paste0("n = ", nrow(test_df)), hjust = 0, vjust = 1,
#                size = 3, family = "Georgia") +
#       coord_fixed(ratio = 1, xlim = c(0, plot_max), ylim = c(0, plot_max)) +
#       common_theme + theme(legend.position = "none")
#     outsample_scatter_regular[[species]] <- p_out
#   }
#   
#   figure3 <- wrap_plots(regular_calib_plots, ncol = 3) /
#     wrap_plots(raw_vs_obs_plots, ncol = 3, guides = "collect") /
#     wrap_plots(raw_vs_calib_regular_plots, ncol = 3, guides = "collect") /
#     wrap_plots(scatter_raw_plots, ncol = 3) /
#     wrap_plots(scatter_calib_regular_plots, ncol = 3) /
#     wrap_plots(outsample_scatter_regular, ncol = 3) +
#     plot_layout(heights = c(1, 1, 1, 1, 1, 1)) +
#     plot_annotation(tag_levels = 'a',
#                     theme = theme(plot.tag = element_text(size = 10, face = "bold", family = "Georgia")))
#   
#   ################################## Display and Save ##################################
#   
#   print(figure1)
#   print(figure2)
#   print(figure3)
#   
#   output_dir <- here("output", "output_figures_tables")
#   if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
#   
#   ggsave(filename = file.path(output_dir, "calibration_comparison_all_methods.png"),
#          plot = figure1, width = 14, height = 14, dpi = 300, bg = "white")
#   cat("\n✓ Figure 1 saved: calibration_comparison_all_methods.png\n")
#   
#   ggsave(filename = file.path(output_dir, "calibration_smoothed_quantile_detailed.png"),
#          plot = figure2, width = 10, height = 12, dpi = 300, bg = "white")
#   cat("✓ Figure 2 saved: calibration_smoothed_quantile_detailed.png\n")
#   
#   ggsave(filename = file.path(output_dir, "calibration_regular_quantile_detailed.png"),
#          plot = figure3, width = 10, height = 12, dpi = 300, bg = "white")
#   cat("✓ Figure 3 saved: calibration_regular_quantile_detailed.png\n")
