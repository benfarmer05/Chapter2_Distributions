# Species toggle - change this to switch between species
species <- "orbicella"
# species <- "agaricia"
# species <- "porites"
# species <- "colpophyllia"
# species <- "madracis"

# Dynamically construct model variable names based on species
presence_model_name <- paste0(species, "_gam_presence_binom")
abundance_model_name <- paste0(species, "_gam_abundance_beta")

# Get the actual model objects
presence_model <- get(presence_model_name)
abundance_model <- get(abundance_model_name)

# Get model data for calibration
data_name <- paste0(species, "_model_data")
all_objects <- ls(envir = .GlobalEnv)
model_data <- if(data_name %in% all_objects) {
  get(data_name, envir = .GlobalEnv)
} else {
  presence_model$model
}

# ============================================================================
# SECTION 1: OPTIMIZE PRESENCE THRESHOLD (5-FOLD CV)
# ============================================================================

# # OPTION 1: OPTIMIZE BY PREVALENCE (comment out the below block and use this instead)
# #
# cat("Optimizing presence threshold via 5-fold cross-validation...\n")
# # Prepare complete data
# all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
# complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
# # 5-fold cross-validation
# # set.seed(42)
# n <- nrow(complete_data)
# fold_ids <- sample(rep(1:5, length.out = n))
# test_thresholds <- seq(0.01, 0.99, by = 0.001)
# prevalence_diffs <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
# for(fold in 1:5) {
#   cat("\n--- Fold", fold, "---\n")
# 
#   # Split data
#   test_idx <- which(fold_ids == fold)
#   train_data <- complete_data[-test_idx, ]
#   test_data <- complete_data[test_idx, ]
# 
#   cat("  Training size:", nrow(train_data), "| Test size:", nrow(test_idx), "\n")
# 
#   # Fit presence model on training data
#   presence_fit <- gam(formula(presence_model), data = train_data, family = binomial())
# 
#   # Predict on test data
#   presence_prob <- predict(presence_fit, newdata = test_data, type = "response")
#   actual_binary <- as.numeric(test_data$cover > 0)
#   actual_prevalence <- mean(actual_binary)
# 
#   cat("  Actual prevalence:", round(actual_prevalence, 4), "\n")
# 
#   # Calculate prevalence difference for each threshold
#   for(i in seq_along(test_thresholds)) {
#     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
#     prevalence_diffs[fold, i] <- abs(mean(predicted_binary) - mean(actual_binary))
#   }
# 
#   # Find best threshold for this fold
#   fold_best_threshold <- test_thresholds[which.min(prevalence_diffs[fold, ])]
#   fold_best_prev_diff <- min(prevalence_diffs[fold, ])
#   predicted_prevalence <- mean(as.numeric(presence_prob > fold_best_threshold))
# 
#   cat("  Best threshold for this fold:", fold_best_threshold, "\n")
#   cat("  Predicted prevalence at best threshold:", round(predicted_prevalence, 4), "\n")
#   cat("  Minimum prevalence difference:", round(fold_best_prev_diff, 6), "\n")
# }
# # Average across folds and find optimal threshold
# cat("\n--- Aggregating Results ---\n")
# mean_prevalence_diffs <- colMeans(prevalence_diffs)
# optimal_threshold <- test_thresholds[which.min(mean_prevalence_diffs)]
# min_mean_prev_diff <- min(mean_prevalence_diffs)
# cat("Mean prevalence difference across folds: ", round(min_mean_prev_diff, 6), "\n")
# cat("Optimal threshold (minimizes prevalence difference): ", optimal_threshold, "\n\n")





# OPTION 2: OPTIMIZE BY TSS (comment out the above block and use this instead)
cat("Optimizing presence threshold via 5-fold cross-validation...\n")
# Prepare complete data
all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
# 5-fold cross-validation
set.seed(42)
n <- nrow(complete_data)
fold_ids <- sample(rep(1:5, length.out = n))
test_thresholds <- seq(0.01, 0.99, by = 0.001)
tss_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
for(fold in 1:5) {
  cat("\n--- Fold", fold, "---\n")

  # Split data
  test_idx <- which(fold_ids == fold)
  train_data <- complete_data[-test_idx, ]
  test_data <- complete_data[test_idx, ]

  cat("  Training size:", nrow(train_data), "| Test size:", nrow(test_idx), "\n")

  # Fit presence model on training data
  presence_fit <- gam(formula(presence_model), data = train_data, family = binomial())

  # Predict on test data
  presence_prob <- predict(presence_fit, newdata = test_data, type = "response")
  actual_binary <- as.numeric(test_data$cover > 0)

  # Calculate TSS for each threshold
  for(i in seq_along(test_thresholds)) {
    predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
    TP <- sum(actual_binary == 1 & predicted_binary == 1)
    TN <- sum(actual_binary == 0 & predicted_binary == 0)
    FP <- sum(actual_binary == 0 & predicted_binary == 1)
    FN <- sum(actual_binary == 1 & predicted_binary == 0)
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    tss_values[fold, i] <- sensitivity + specificity - 1
  }

  # Find best threshold for this fold
  fold_best_threshold <- test_thresholds[which.max(tss_values[fold, ])]
  fold_best_tss <- max(tss_values[fold, ])

  cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  cat("  Maximum TSS:", round(fold_best_tss, 6), "\n")
}
# Average across folds and find optimal threshold
cat("\n--- Aggregating Results ---\n")
mean_tss_values <- colMeans(tss_values)
optimal_threshold <- test_thresholds[which.max(mean_tss_values)]
max_mean_tss <- max(mean_tss_values)
cat("Mean TSS across folds: ", round(max_mean_tss, 6), "\n")
cat("Optimal threshold (maximizes TSS): ", optimal_threshold, "\n\n")










# # Manual override option (comment out to use optimized threshold)
# optimal_threshold <- 0.5

# Evaluate optimal threshold on entire dataset
cat("\n=== PERFORMANCE ON ENTIRE DATASET ===\n")
full_presence_fit <- gam(formula(presence_model), data = complete_data, family = binomial())
full_presence_prob <- predict(full_presence_fit, type = "response")
full_actual_binary <- as.numeric(complete_data$cover > 0)
full_predicted_binary <- as.numeric(full_presence_prob > optimal_threshold)

# Calculate metrics
TP <- sum(full_actual_binary == 1 & full_predicted_binary == 1)
TN <- sum(full_actual_binary == 0 & full_predicted_binary == 0)
FP <- sum(full_actual_binary == 0 & full_predicted_binary == 1)
FN <- sum(full_actual_binary == 1 & full_predicted_binary == 0)

sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)
TSS <- sensitivity + specificity - 1
actual_prevalence <- mean(full_actual_binary)
predicted_prevalence <- mean(full_predicted_binary)
prevalence_diff <- abs(predicted_prevalence - actual_prevalence)

library(pROC)
full_auc <- as.numeric(auc(full_actual_binary, full_presence_prob))

cat("Threshold:", optimal_threshold, "\n")
cat("Sensitivity:", round(sensitivity, 3), "| Specificity:", round(specificity, 3), "\n")
cat("TSS:", round(TSS, 3), "| AUC:", round(full_auc, 3), "\n")
cat("Actual prevalence:", round(actual_prevalence, 4), "| Predicted prevalence:", round(predicted_prevalence, 4), "\n")
cat("Prevalence difference:", round(prevalence_diff, 6), "\n\n")


# ============================================================================
# SECTION 2: TRAIN ABUNDANCE CALIBRATION FUNCTION
# ============================================================================

cat("Training abundance calibration function...\n")

# Helper function for quantile calibration
calibrate_quantiles <- function(predictions, observations, n_quantiles = 100) {
  probs <- seq(0, 1, length.out = n_quantiles)
  pred_q <- quantile(predictions, probs)
  obs_q <- quantile(observations, probs)
  
  correction_fn <- approxfun(pred_q, obs_q, method = "linear", rule = 2)
  
  list(correction_function = correction_fn)
}

# Get abundance training data (only where cover > 0)
abundance_train <- complete_data[complete_data$cover > 0, ]
abundance_train$cover_prop <- abundance_train$cover / 100

# Predict on training data
train_abundance_pred <- predict(abundance_model, newdata = abundance_train, type = "response")

# Create calibration function
cal_result <- calibrate_quantiles(train_abundance_pred, abundance_train$cover_prop, n_quantiles = 100)
abundance_correction_fn <- cal_result$correction_function

cat("Abundance calibration function trained\n\n")

# ============================================================================
# SECTION 3: PREDICT ON FULL SPATIAL GRID
# ============================================================================

cat("Preparing spatial predictions...\n")

# Manually specify required variables
presence_vars <- all.vars(formula(presence_model))[-1]
abundance_vars <- all.vars(formula(abundance_model))[-1]
required_vars <- unique(c(presence_vars, abundance_vars))

# Subset raster stack to only required variables
required_vars <- gsub("^depth_bathy$", "depth", required_vars)
env_subset <- env_complex[[required_vars]]

# Convert subset raster to dataframe
env_df <- as.data.frame(env_subset, xy = TRUE, na.rm = TRUE)

# Rename depth column to match model expectation
names(env_df)[names(env_df) == "depth"] <- "depth_bathy"

# Sample for testing (1/100th of data)
sample_size <- ceiling(nrow(env_df) / 100)
# set.seed(123)
sample_indices <- sample(nrow(env_df), sample_size)
env_df_sample <- env_df[sample_indices, ]

cat("Species:", species, "\n")
cat("Full dataset size:", nrow(env_df), "cells\n")
cat("Sample size (1/100th):", nrow(env_df_sample), "cells\n\n")

# ============================================================================
# SECTION 4: MAKE PREDICTIONS AND APPLY CALIBRATION
# ============================================================================

# Presence predictions
cat("Predicting presence...\n")
start_time <- Sys.time()
presence_prob_sample <- predict(presence_model, newdata = env_df_sample, type = "response")
presence_time <- difftime(Sys.time(), start_time, units = "secs")
cat("Sample presence prediction:", round(presence_time, 2), "seconds\n")

# Apply optimized threshold
presence_binary_sample <- ifelse(presence_prob_sample > optimal_threshold, 1, 0)

# Abundance predictions
cat("Predicting abundance...\n")
start_time <- Sys.time()
abundance_pred_sample <- predict(abundance_model, newdata = env_df_sample, type = "response")
abundance_time <- difftime(Sys.time(), start_time, units = "secs")
cat("Sample abundance prediction:", round(abundance_time, 2), "seconds\n")

# Apply abundance calibration
abundance_pred_calibrated <- abundance_correction_fn(abundance_pred_sample)

# Create hurdle predictions with calibrated abundance
hurdle_pred_sample <- ifelse(presence_binary_sample == 1, abundance_pred_calibrated, 0)

# Add predictions to sample dataframe
env_df_sample$presence_prob <- presence_prob_sample
env_df_sample$presence_binary <- presence_binary_sample
env_df_sample$abundance_pred_raw <- abundance_pred_sample
env_df_sample$abundance_pred_calibrated <- abundance_pred_calibrated
env_df_sample$hurdle_pred <- hurdle_pred_sample

# Time estimates
total_sample_time <- as.numeric(presence_time + abundance_time)
estimated_full_time <- total_sample_time * 100

cat("\nSample results summary:\n")
cat("Presence threshold used:", optimal_threshold, "\n")
cat("Presence probability range:", round(range(presence_prob_sample), 4), "\n")
cat("Abundance prediction range (raw):", round(range(abundance_pred_sample), 6), "\n")
cat("Abundance prediction range (calibrated):", round(range(abundance_pred_calibrated), 6), "\n")
cat("Hurdle prediction range:", round(range(hurdle_pred_sample), 6), "\n")

cat("\nTime estimates:\n")
cat("Sample time:", round(total_sample_time, 2), "seconds\n")
cat("Estimated full dataset time:", round(estimated_full_time, 1), "seconds")
if(estimated_full_time > 60) {
  cat(" (", round(estimated_full_time/60, 1), " minutes)", sep="")
}
cat("\n\n")

# ============================================================================
# SECTION 5: SAMPLE VISUALIZATIONS
# ============================================================================

library(ggplot2)
library(viridis)

# Presence probability map
ggplot(env_df_sample, aes(x = x, y = y, color = presence_prob)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(name = "Presence\nProbability") +
  coord_equal() +
  theme_minimal() +
  ggtitle(paste(stringr::str_to_title(species), "Presence Probability (Sample)"))

# Hurdle prediction map (with calibrated abundance)
ggplot(env_df_sample, aes(x = x, y = y, color = hurdle_pred)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(name = "Predicted\nCover", limits = c(0, 0.1)) +
  coord_equal() +
  theme_minimal() +
  ggtitle(paste("Hurdle Model Prediction (Calibrated) -", 
                stringr::str_to_title(species), 
                "\nThreshold:", optimal_threshold, "(Sample)"))

# Histogram of non-zero predictions
hist(env_df_sample$hurdle_pred[env_df_sample$hurdle_pred > 0 & env_df_sample$hurdle_pred < 1],
     main = "Distribution of Non-Zero Hurdle Predictions (Sample)",
     xlab = "Predicted Cover")

# # ============================================================================
# # SECTION 6: FULL SPATIAL PREDICTIONS (SERIAL - DEFAULT)
# # ============================================================================
# # Comment out this section if using parallelized version below (Section 6.5)
# 
# cat("\n=== RUNNING FULL MODEL PREDICTIONS (SERIAL) ===\n")
# cat("Warning: This may take 15-20+ minutes depending on dataset size\n\n")
# 
# start_time <- Sys.time()
# 
# # Full presence prediction
# cat("Predicting presence on full dataset...\n")
# cat("  Started at:", format(Sys.time(), "%H:%M:%S"), "\n")
# presence_prob_full <- predict(presence_model, newdata = env_df, type = "response")
# cat("  Finished at:", format(Sys.time(), "%H:%M:%S"), "\n")
# presence_binary_full <- ifelse(presence_prob_full > optimal_threshold, 1, 0)
# 
# # Full abundance prediction
# cat("Predicting abundance on full dataset...\n")
# cat("  Started at:", format(Sys.time(), "%H:%M:%S"), "\n")
# abundance_pred_full <- predict(abundance_model, newdata = env_df, type = "response")
# cat("  Finished at:", format(Sys.time(), "%H:%M:%S"), "\n")
# 
# # Apply calibration to full abundance predictions
# cat("Applying calibration to abundance predictions...\n")
# abundance_pred_full_calibrated <- abundance_correction_fn(abundance_pred_full)
# 
# # Full hurdle prediction
# hurdle_pred_full <- ifelse(presence_binary_full == 1, abundance_pred_full_calibrated, 0)

# ============================================================================
# SECTION 6.5: FULL SPATIAL PREDICTIONS (PARALLELIZED - OPTIONAL)
# ============================================================================
# Uncomment this section (and comment out Section 6) for faster predictions
# Requires: parallel package
# Expected speedup: 2-4x on typical home machines

library(parallel)

cat("\n=== RUNNING FULL MODEL PREDICTIONS (PARALLELIZED) ===\n")
cat("Warning: This may take 5-10 minutes depending on cores and dataset size\n\n")

start_time <- Sys.time()

# Set up parallel cluster
n_cores <- detectCores() - 1  # Leave one core free
cat("Using", n_cores, "cores\n\n")
cl <- makeCluster(n_cores)

# Split data into chunks
chunks <- split(1:nrow(env_df), cut(1:nrow(env_df), n_cores))

# ---- Presence prediction ----
cat("Predicting presence on full dataset (parallel)...\n")
cat("  Started at:", format(Sys.time(), "%H:%M:%S"), "\n")

# Export to workers
clusterExport(cl, c("presence_model", "env_df"), envir = environment())
clusterEvalQ(cl, library(mgcv))

# Parallel prediction
presence_prob_list <- parLapply(cl, chunks, function(idx) {
  predict(presence_model, newdata = env_df[idx, ], type = "response")
})
presence_prob_full <- unlist(presence_prob_list)

cat("  Finished at:", format(Sys.time(), "%H:%M:%S"), "\n")
presence_binary_full <- ifelse(presence_prob_full > optimal_threshold, 1, 0)

# ---- Abundance prediction ----
cat("Predicting abundance on full dataset (parallel)...\n")
cat("  Started at:", format(Sys.time(), "%H:%M:%S"), "\n")

# Export abundance model
clusterExport(cl, c("abundance_model"), envir = environment())

# Parallel prediction
abundance_pred_list <- parLapply(cl, chunks, function(idx) {
  predict(abundance_model, newdata = env_df[idx, ], type = "response")
})
abundance_pred_full <- unlist(abundance_pred_list)

cat("  Finished at:", format(Sys.time(), "%H:%M:%S"), "\n")

# Clean up cluster
stopCluster(cl)

# Apply calibration to full abundance predictions
cat("Applying calibration to abundance predictions...\n")
abundance_pred_full_calibrated <- abundance_correction_fn(abundance_pred_full)

# Full hurdle prediction
hurdle_pred_full <- ifelse(presence_binary_full == 1, abundance_pred_full_calibrated, 0)

# ============================================================================
# SECTION 6.9: PREDICTION SUMMARY (COMMON TO BOTH METHODS)
# ============================================================================

full_time <- difftime(Sys.time(), start_time, units = "secs")
cat("\nFull prediction time:", round(full_time, 1), "seconds (", 
    round(full_time/60, 1), "minutes )\n")

# Add predictions to full dataframe
env_df$presence_prob <- presence_prob_full
env_df$presence_binary <- presence_binary_full
env_df$abundance_pred_raw <- abundance_pred_full
env_df$abundance_pred_calibrated <- abundance_pred_full_calibrated
env_df$hurdle_pred <- hurdle_pred_full

cat("\nFull prediction summary:\n")
cat("Presence probability range:", round(range(presence_prob_full), 4), "\n")
cat("Predicted prevalence:", round(mean(presence_binary_full), 4), "\n")
cat("Hurdle prediction range:", round(range(hurdle_pred_full), 6), "\n")
cat("Number of predicted presences:", sum(presence_binary_full), "\n\n")

# ============================================================================
# SECTION 7: CONVERT TO RASTER AND CREATE INTERACTIVE MAP
# ============================================================================

library(terra)
library(leaflet)

cat("Converting predictions to raster...\n")

# Convert to raster
hurdle_raster <- rast(env_df[, c("x", "y", "hurdle_pred")], type = "xyz")

# Apply shallow mask (if bathy_final exists)
if(exists("bathy_final")) {
  shallow_mask <- bathy_final > -60
  hurdle_raster <- resample(hurdle_raster, bathy_final)
  hurdle_raster <- mask(hurdle_raster, shallow_mask, maskvalues = FALSE)
}

# Clamp values for visualization
cover_clamp_val = 5 #in percent
cover_clamp_val = cover_clamp_val / 100
hurdle_raster_clamped <- clamp(hurdle_raster, lower = 0, upper = cover_clamp_val, values = TRUE)

cat("Creating interactive leaflet map...\n")

# Create leaflet map
pal <- colorNumeric("viridis", values(hurdle_raster_clamped), na.color = "transparent")

interactive_map <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(hurdle_raster_clamped, colors = pal, opacity = 0.7) %>%
  addLegend("bottomright", pal = pal, values = values(hurdle_raster_clamped),
            title = paste(stringr::str_to_title(species), "Cover (%)")) %>%
  addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                   options = layersControlOptions(collapsed = FALSE))

# Display the map
interactive_map

cat("\n=== PREDICTION COMPLETE ===\n")
cat("Threshold used:", optimal_threshold, "\n")
cat("Calibration applied: Yes\n")
cat("Full dataset predicted:", nrow(env_df), "cells\n")



#WITH CORAL COVER VERSION (VIRIDIS)
psu_cover_data <- combined_benthic_data_averaged %>%
  filter(grepl(species, spp, ignore.case = TRUE)) %>%
  group_by(PSU) %>%
  summarise(lat = first(lat), lon = first(lon), 
            avg_cover = mean(cover, na.rm = TRUE), .groups = 'drop')

# Clamp cover data to match cover_clamp_val (e.g., 10% = 0.1 in proportion)
psu_cover_data <- psu_cover_data %>%
  mutate(avg_cover_clamped = pmin(avg_cover / 100, cover_clamp_val))  # Convert to proportion AND clamp

# Use viridis palette for both, 0-cover_clamp_val proportion scale
unified_pal <- colorNumeric("viridis", 
                            domain = c(0, cover_clamp_val), na.color = "transparent")

interactive_map <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addRasterImage(hurdle_raster_clamped, colors = unified_pal, opacity = 0.7) %>%
  addCircleMarkers(data = psu_cover_data, ~lon, ~lat, radius = 8,
                   fillColor = ~unified_pal(avg_cover_clamped), fillOpacity = 0.8,  # Already in proportion
                   color = "white", weight = 2, stroke = TRUE,
                   popup = ~paste0("<b>PSU:</b> ", PSU, "<br><b>Cover:</b> ", round(avg_cover, 2), "%")) %>%
  addLegend("bottomright", pal = unified_pal, values = values(hurdle_raster_clamped),
            title = paste(stringr::str_to_title(species), "Cover (%)"),
            labFormat = labelFormat(transform = function(x) x * 100, suffix = "%")) %>%
  addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                   options = layersControlOptions(collapsed = FALSE))
interactive_map