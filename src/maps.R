  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(mgcv)
  library(pROC)
  library(terra)
  library(ggplot2)
  library(parallel)
  library(leaflet)
  library(tidyverse)
  
  source(here("src/functions.R"))
  
  ################################## setup ##################################
  
  load(here("output", "all_combined_data.rda"))
  
  
  #load bathy_final
  # Load the spatial metadata
  spatial_metadata <- readRDS(here('output/output_import_merge_rasters_higher-res/spatial_metadata.rds'))
  #
  # Load the raster
  bathy_final <- readRDS(here('output/output_import_merge_rasters_higher-res/bathy_final.rds'))
  #
  # Apply the stored CRS
  terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
  
  
  load_spat_objects(directory = 'output/output_GAMs/')
  source(here("src/functions.R"))
  load(here("output/output_GAMs", "output_GAMs_workspace.Rdata"))
  
  
  #load in the species-specific GAM model output
  #   NOTE - this will overwrite any species model data that were in 
  #             'output_GAMs_workspace'. that is fine, as I just want to make sure
  #             that the models are easily update-able.
  #Load all species models
  
  species_list <- c("agaricia", "colpophyllia", "dendrogyra", "dichocoenia",
                    "diploria", "eusmilia", "madracis", "meandrina",
                    "montastraea", "mycetophyllia", "orbicella", "porites",
                    "pseudodiploria", "siderastrea", "solenastrea")
  
  # Load all species models into a list
  all_species_models <- lapply(species_list, function(sp) {
    file_path <- here("output", "output_GAMs", paste0(sp, "_models.rds"))
    if(file.exists(file_path)) {
      readRDS(file_path)
    } else {
      warning(paste("Model file not found for:", sp))
      NULL
    }
  })
  names(all_species_models) <- species_list
  
  # Quick check of what loaded
  cat("Loaded models for:\n")
  for(sp in species_list) {
    if(!is.null(all_species_models[[sp]])) {
      cat("  ✓", sp, "\n")
    } else {
      cat("  ✗", sp, "(missing)\n")
    }
  }
  

  
  
  #save information for exporting new objects later
  existing_objects <- ls(envir = .GlobalEnv)
  
  
  ################################## select species ##################################
  
  # Species toggle - change this to switch between species
  # species <- "agaricia" # N = 651. best is target prev. ratio (OPTION 6C). looks good
  # species <- "colpophyllia" # N = 116. best is sens. w/ constraints (OPTION 8). a bit overpredicted but good
  # species <- "dendrogyra" # N = 35. best is sens. w/ EVEN MORE constraints (OPTION 9) - looks AWESOME!!
  # species <- "dichocoenia" # N = 76. best is sens. w/ MORE constraints (OPTION 9), could even use EVEN MORE constraints to ratchet down overprediction
  # species <- "diploria" # N = 120. best is sens. w/ constraints (OPTION 8). a bit overpredicted but good overall
  # species <- "eusmilia" # N = 39. best is sens. w/ EVEN MORE constraints (OPTION 9). looks excellent
  # species <- "madracis" # N = 131. best is sens. w/ constraints (OPTION 8). great!!
  # species <- "meandrina" # N = 193. best is sens. w/ constraints (OPTION 8). prob a bit overpredicted
  # species <- "montastraea" # N = 490. best is target prev. ratio (OPTION 6C). looks good! maybe a bit overpredicted
  # species <- "mycetophyllia" # N = 28. best is sens. w/ EVEN MORE constraints (OPTION 9). pretty good, slightly overpredicted
  # species <- "orbicella" # N = 795. best is target prev. ratio (OPTION 6C), though overstates extent a bit most likely. Kappa was useful for narrowing toward an ideal prev. ratio in the first place (useful for all common species)
  # species <- "porites" # N = 1016. best is target prev. ratio (OPTION 6C). looks good! maybe underestimates a bit...not sure (definitely in MCD though)
  # species <- "pseudodiploria" # N = 401. best is target prev. ratio (OPTION 6C), but really misses mesophotic pstrig (tried sens. w/ LESS constraint to handle slight underestimation but this just overstated extent in the wrong places)
  species <- "siderastrea" # N = 845. best is target prev. ratio (OPTION 6C). looks good! probably slightly overstated
  # species <- "solenastrea" # N = 21. tried sens. w/ EVEN MORE constraints (OPTION 9), but predictions are much too overestimated. may drop this species
  
  # Extract models for this species
  presence_model <- all_species_models[[species]]$presence_model
  abundance_model <- all_species_models[[species]]$abundance_model
  model_data <- all_species_models[[species]]$model_data
  
  
  # # Dynamically construct model variable names based on species
  # presence_model_name <- paste0(species, "_gam_presence_binom")
  # abundance_model_name <- paste0(species, "_gam_abundance_beta")
  # 
  # # Get the actual model objects
  # presence_model <- get(presence_model_name)
  # abundance_model <- get(abundance_model_name)
  # 
  # # Get model data for calibration
  # data_name <- paste0(species, "_model_data")
  # all_objects <- ls(envir = .GlobalEnv)
  # model_data <- if(data_name %in% all_objects) {
  #   get(data_name, envir = .GlobalEnv)
  # } else {
  #   presence_model$model
  # }
  
  
  
  
  
  ###                                                    ###
  ###                                                    ###
  ###     BELOW: options for optimizing P/A threshold    ###
  ###                                                    ###
  ###                                                    ###
  
  # ################################## OPTION 1: prevalence ##################################
  # 
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
  # 
  # 
  # 
  # ################################## OPTION 2: TSS ##################################
  # 
  # # OPTION 2: OPTIMIZE BY TSS (comment out the above block and use this instead)
  # cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # test_thresholds <- seq(0.0001, 0.99, by = 0.00001)
  # tss_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
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
  # 
  #   # Calculate TSS for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  #     TP <- sum(actual_binary == 1 & predicted_binary == 1)
  #     TN <- sum(actual_binary == 0 & predicted_binary == 0)
  #     FP <- sum(actual_binary == 0 & predicted_binary == 1)
  #     FN <- sum(actual_binary == 1 & predicted_binary == 0)
  #     sensitivity <- TP / (TP + FN)
  #     specificity <- TN / (TN + FP)
  #     tss_values[fold, i] <- sensitivity + specificity - 1
  #   }
  # 
  #   # Find best threshold for this fold
  #   fold_best_threshold <- test_thresholds[which.max(tss_values[fold, ])]
  #   fold_best_tss <- max(tss_values[fold, ])
  # 
  #   cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #   cat("  Maximum TSS:", round(fold_best_tss, 6), "\n")
  # }
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_tss_values <- colMeans(tss_values)
  # optimal_threshold <- test_thresholds[which.max(mean_tss_values)]
  # max_mean_tss <- max(mean_tss_values)
  # cat("Mean TSS across folds: ", round(max_mean_tss, 6), "\n")
  # cat("Optimal threshold (maximizes TSS): ", optimal_threshold, "\n\n")
  # 
  # ################################## OPTION 3: Kappa (maybe best for orbicella to set parameters) ##################################
  # 
  # cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # kappa_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
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
  # 
  #   # Calculate Kappa for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  #     TP <- sum(actual_binary == 1 & predicted_binary == 1)
  #     TN <- sum(actual_binary == 0 & predicted_binary == 0)
  #     FP <- sum(actual_binary == 0 & predicted_binary == 1)
  #     FN <- sum(actual_binary == 1 & predicted_binary == 0)
  # 
  #     n_total <- length(actual_binary)
  #     observed_accuracy <- (TP + TN) / n_total
  #     expected_accuracy <- ((TP + FN) * (TP + FP) + (TN + FP) * (TN + FN)) / (n_total^2)
  #     kappa_values[fold, i] <- (observed_accuracy - expected_accuracy) / (1 - expected_accuracy)
  #   }
  # 
  #   # Find best threshold for this fold
  #   fold_best_threshold <- test_thresholds[which.max(kappa_values[fold, ])]
  #   fold_best_kappa <- max(kappa_values[fold, ], na.rm = TRUE)
  # 
  #   cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #   cat("  Maximum Kappa:", round(fold_best_kappa, 6), "\n")
  # }
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_kappa_values <- colMeans(kappa_values, na.rm = TRUE)
  # optimal_threshold <- test_thresholds[which.max(mean_kappa_values)]
  # max_mean_kappa <- max(mean_kappa_values, na.rm = TRUE)
  # cat("Mean Kappa across folds: ", round(max_mean_kappa, 6), "\n")
  # cat("Optimal threshold (maximizes Kappa): ", optimal_threshold, "\n\n")
  # 
  # 
  # 
  # 
  # 
  # 
  # ################################## OPTION 4: F2-score ##################################
  # 
  # cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # f2_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
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
  # 
  #   # Calculate F2-score for each threshold (weights recall/sensitivity 2x more than precision)
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  #     TP <- sum(actual_binary == 1 & predicted_binary == 1)
  #     FP <- sum(actual_binary == 0 & predicted_binary == 1)
  #     FN <- sum(actual_binary == 1 & predicted_binary == 0)
  # 
  #     precision <- if((TP + FP) > 0) TP / (TP + FP) else 0
  #     recall <- if((TP + FN) > 0) TP / (TP + FN) else 0
  #     beta <- 2  # Weights recall (sensitivity) MORE than precision
  #     f2_values[fold, i] <- if((precision + recall) > 0) {
  #       (1 + beta^2) * (precision * recall) / ((beta^2 * precision) + recall)
  #     } else 0
  #   }
  # 
  #   # Find best threshold for this fold
  #   fold_best_threshold <- test_thresholds[which.max(f2_values[fold, ])]
  #   fold_best_f2 <- max(f2_values[fold, ], na.rm = TRUE)
  # 
  #   cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #   cat("  Maximum F2-score:", round(fold_best_f2, 6), "\n")
  # }
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_f2_values <- colMeans(f2_values, na.rm = TRUE)
  # optimal_threshold <- test_thresholds[which.max(mean_f2_values)]
  # max_mean_f2 <- max(mean_f2_values, na.rm = TRUE)
  # cat("Mean F2-score across folds: ", round(max_mean_f2, 6), "\n")
  # cat("Optimal threshold (maximizes F2-score): ", optimal_threshold, "\n\n")
  # 
  # 
  # 
  # 
  # ################################## OPTION 5: balance sens. & spec. ##################################
  # 
  # cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # # 5-fold cross-validation
  # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # balanced_acc_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
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
  # 
  #   # Calculate Balanced Accuracy for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  #     TP <- sum(actual_binary == 1 & predicted_binary == 1)
  #     TN <- sum(actual_binary == 0 & predicted_binary == 0)
  #     FP <- sum(actual_binary == 0 & predicted_binary == 1)
  #     FN <- sum(actual_binary == 1 & predicted_binary == 0)
  # 
  #     sensitivity <- TP / (TP + FN)
  #     specificity <- TN / (TN + FP)
  #     # Balanced Accuracy = average of sensitivity and specificity
  #     balanced_acc_values[fold, i] <- (sensitivity + specificity) / 2
  #   }
  # 
  #   # Find best threshold for this fold
  #   fold_best_threshold <- test_thresholds[which.max(balanced_acc_values[fold, ])]
  #   fold_best_balanced_acc <- max(balanced_acc_values[fold, ], na.rm = TRUE)
  # 
  #   cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #   cat("  Maximum Balanced Accuracy:", round(fold_best_balanced_acc, 6), "\n")
  # }
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_balanced_acc_values <- colMeans(balanced_acc_values, na.rm = TRUE)
  # optimal_threshold <- test_thresholds[which.max(mean_balanced_acc_values)]
  # max_mean_balanced_acc <- max(mean_balanced_acc_values, na.rm = TRUE)
  # cat("Mean Balanced Accuracy across folds: ", round(max_mean_balanced_acc, 6), "\n")
  # cat("Optimal threshold (maximizes Balanced Accuracy): ", optimal_threshold, "\n\n")
  # 
  # 
  # 
  # ################################## OPTION 6a: target sens. & spec. ##################################
  # 
  # cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # 
  # # Target values from well-performing model (e.g., Orbicella)
  # target_sensitivity <- 0.692
  # target_specificity <- 0.731
  # 
  # cat("Target sensitivity:", target_sensitivity, "| Target specificity:", target_specificity, "\n\n")
  # 
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # target_distance <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # 
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
  # 
  #   # Calculate distance from target for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  #     TP <- sum(actual_binary == 1 & predicted_binary == 1)
  #     TN <- sum(actual_binary == 0 & predicted_binary == 0)
  #     FP <- sum(actual_binary == 0 & predicted_binary == 1)
  #     FN <- sum(actual_binary == 1 & predicted_binary == 0)
  # 
  #     sensitivity <- TP / (TP + FN)
  #     specificity <- TN / (TN + FP)
  # 
  #     # Euclidean distance from target (equal weighting for both metrics)
  #     sens_diff <- (sensitivity - target_sensitivity)^2
  #     spec_diff <- (specificity - target_specificity)^2
  #     target_distance[fold, i] <- sqrt(sens_diff + spec_diff)
  #   }
  # 
  #   # Find best threshold for this fold (minimum distance from target)
  #   fold_best_threshold <- test_thresholds[which.min(target_distance[fold, ])]
  #   fold_min_distance <- min(target_distance[fold, ], na.rm = TRUE)
  # 
  #   cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #   cat("  Minimum distance from target:", round(fold_min_distance, 6), "\n")
  # }
  # 
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_target_distance <- colMeans(target_distance, na.rm = TRUE)
  # optimal_threshold <- test_thresholds[which.min(mean_target_distance)]
  # min_mean_distance <- min(mean_target_distance, na.rm = TRUE)
  # cat("Mean distance from target across folds: ", round(min_mean_distance, 6), "\n")
  # cat("Optimal threshold (minimizes distance from target): ", optimal_threshold, "\n\n")
  # 
  # 
  # ################################## OPTION 6b: target sens.:spec. ratio ##################################
  # 
  # cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # 
  # # Target ratio from well-performing model (e.g., Orbicella: 0.692/0.731 = 0.9467)
  # target_ratio <- 0.692 / 0.731
  # 
  # cat("Target sensitivity:specificity ratio:", round(target_ratio, 4), "\n\n")
  # 
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # ratio_diff <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # 
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
  # 
  #   # Calculate ratio difference for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  #     TP <- sum(actual_binary == 1 & predicted_binary == 1)
  #     TN <- sum(actual_binary == 0 & predicted_binary == 0)
  #     FP <- sum(actual_binary == 0 & predicted_binary == 1)
  #     FN <- sum(actual_binary == 1 & predicted_binary == 0)
  # 
  #     sensitivity <- TP / (TP + FN)
  #     specificity <- TN / (TN + FP)
  # 
  #     # Calculate ratio (avoid division by zero)
  #     current_ratio <- if(specificity > 0) sensitivity / specificity else NA
  # 
  #     # Absolute difference from target ratio
  #     ratio_diff[fold, i] <- abs(current_ratio - target_ratio)
  #   }
  # 
  #   # Find best threshold for this fold (minimum ratio difference)
  #   fold_best_threshold <- test_thresholds[which.min(ratio_diff[fold, ])]
  #   fold_min_diff <- min(ratio_diff[fold, ], na.rm = TRUE)
  # 
  #   cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #   cat("  Minimum ratio difference:", round(fold_min_diff, 6), "\n")
  # }
  # 
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_ratio_diff <- colMeans(ratio_diff, na.rm = TRUE)
  # optimal_threshold <- test_thresholds[which.min(mean_ratio_diff)]
  # min_mean_ratio_diff <- min(mean_ratio_diff, na.rm = TRUE)
  # cat("Mean ratio difference across folds: ", round(min_mean_ratio_diff, 6), "\n")
  # cat("Optimal threshold (minimizes ratio difference from target): ", optimal_threshold, "\n\n")
  # 
  # 
  ################################## OPTION 6c: target prev. ratio ##################################

  cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # Prepare complete data
  all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  complete_data <- model_data[complete.cases(model_data[, all_vars]), ]

  # Target prevalence ratio from well-calibrated model (e.g., Orbicella with Kappa)
  # 0.4412 / 0.398 = 1.109 (slight overprediction)
  target_prev_ratio <- 1.109

  cat("Target predicted:actual prevalence ratio:", round(target_prev_ratio, 3), "\n")
  cat("(Target allows ~", round((target_prev_ratio - 1) * 100, 1), "% overprediction)\n\n", sep="")

  # 5-fold cross-validation
  # set.seed(42)
  n <- nrow(complete_data)
  fold_ids <- sample(rep(1:5, length.out = n))
  test_thresholds <- seq(0.01, 0.99, by = 0.001)
  prev_ratio_diff <- matrix(NA, nrow = 5, ncol = length(test_thresholds))

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
    actual_prev <- mean(actual_binary)

    # Calculate prevalence ratio difference for each threshold
    for(i in seq_along(test_thresholds)) {
      predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
      predicted_prev <- mean(predicted_binary)

      # Calculate ratio (avoid division by zero)
      current_ratio <- if(actual_prev > 0) predicted_prev / actual_prev else NA

      # Absolute difference from target ratio
      prev_ratio_diff[fold, i] <- abs(current_ratio - target_prev_ratio)
    }

    # Find best threshold for this fold
    fold_best_threshold <- test_thresholds[which.min(prev_ratio_diff[fold, ])]
    fold_min_diff <- min(prev_ratio_diff[fold, ], na.rm = TRUE)

    # Get actual metrics at this threshold
    fold_pred_binary <- as.numeric(presence_prob > fold_best_threshold)
    fold_pred_prev <- mean(fold_pred_binary)
    fold_actual_ratio <- fold_pred_prev / actual_prev

    cat("  Best threshold for this fold:", fold_best_threshold, "\n")
    cat("  Achieved ratio:", round(fold_actual_ratio, 3),
        "| Target:", round(target_prev_ratio, 3), "\n")
  }

  # Average across folds and find optimal threshold
  cat("\n--- Aggregating Results ---\n")
  mean_prev_ratio_diff <- colMeans(prev_ratio_diff, na.rm = TRUE)
  optimal_threshold <- test_thresholds[which.min(mean_prev_ratio_diff)]
  min_mean_diff <- min(mean_prev_ratio_diff, na.rm = TRUE)
  cat("Mean prevalence ratio difference across folds: ", round(min_mean_diff, 6), "\n")
  cat("Optimal threshold (minimizes difference from target ratio): ", optimal_threshold, "\n\n")

  # ################################## OPTION 6d: target higher prev. ratio ##################################
  # 
  # cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # 
  # # Target prevalence ratio from well-calibrated model (e.g., Orbicella with Kappa)
  # # 0.4412 / 0.398 = 1.109 (slight overprediction)
  # target_prev_ratio <- 1.60
  # 
  # cat("Target predicted:actual prevalence ratio:", round(target_prev_ratio, 3), "\n")
  # cat("(Target allows ~", round((target_prev_ratio - 1) * 100, 1), "% overprediction)\n\n", sep="")
  # 
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # prev_ratio_diff <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # 
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
  #   actual_prev <- mean(actual_binary)
  # 
  #   # Calculate prevalence ratio difference for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  #     predicted_prev <- mean(predicted_binary)
  # 
  #     # Calculate ratio (avoid division by zero)
  #     current_ratio <- if(actual_prev > 0) predicted_prev / actual_prev else NA
  # 
  #     # Absolute difference from target ratio
  #     prev_ratio_diff[fold, i] <- abs(current_ratio - target_prev_ratio)
  #   }
  # 
  #   # Find best threshold for this fold
  #   fold_best_threshold <- test_thresholds[which.min(prev_ratio_diff[fold, ])]
  #   fold_min_diff <- min(prev_ratio_diff[fold, ], na.rm = TRUE)
  # 
  #   # Get actual metrics at this threshold
  #   fold_pred_binary <- as.numeric(presence_prob > fold_best_threshold)
  #   fold_pred_prev <- mean(fold_pred_binary)
  #   fold_actual_ratio <- fold_pred_prev / actual_prev
  # 
  #   cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #   cat("  Achieved ratio:", round(fold_actual_ratio, 3),
  #       "| Target:", round(target_prev_ratio, 3), "\n")
  # }
  # 
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_prev_ratio_diff <- colMeans(prev_ratio_diff, na.rm = TRUE)
  # optimal_threshold <- test_thresholds[which.min(mean_prev_ratio_diff)]
  # min_mean_diff <- min(mean_prev_ratio_diff, na.rm = TRUE)
  # cat("Mean prevalence ratio difference across folds: ", round(min_mean_diff, 6), "\n")
  # cat("Optimal threshold (minimizes difference from target ratio): ", optimal_threshold, "\n\n")
  # 
  # 
  # 
  # ################################## OPTION 6e: prevalence w/ constraint ##################################
  # 
  # cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # 
  # # Calculate actual prevalence
  # actual_prevalence_full <- mean(complete_data$cover > 0)
  # 
  # # Set prevalence-adaptive target ratio (always ≥ 1.0)
  # if(actual_prevalence_full >= 0.25) {
  #   target_prev_ratio <- 1.08  # 8% overprediction for common species
  #   cat("Species prevalence:", round(actual_prevalence_full, 3), "(common)\n")
  # } else if(actual_prevalence_full >= 0.15) {
  #   target_prev_ratio <- 1.12  # 12% overprediction for moderately common
  #   cat("Species prevalence:", round(actual_prevalence_full, 3), "(moderately common)\n")
  # } else if(actual_prevalence_full >= 0.08) {
  #   target_prev_ratio <- 1.25  # 25% overprediction for moderately rare
  #   cat("Species prevalence:", round(actual_prevalence_full, 3), "(moderately rare)\n")
  # } else {
  #   target_prev_ratio <- 2.00  # or could try 50% overprediction for rare species
  #   cat("Species prevalence:", round(actual_prevalence_full, 3), "(rare)\n")
  # }
  # 
  # cat("Target predicted:actual prevalence ratio:", round(target_prev_ratio, 3), "\n")
  # cat("(Allows ~", round((target_prev_ratio - 1) * 100, 0), "% overprediction)\n\n", sep="")
  # 
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # prev_ratio_diff <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # 
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
  #   actual_prev <- mean(actual_binary)
  # 
  #   # Calculate prevalence ratio difference for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  #     predicted_prev <- mean(predicted_binary)
  # 
  #     # Calculate ratio (avoid division by zero)
  #     current_ratio <- if(actual_prev > 0) predicted_prev / actual_prev else NA
  # 
  #     # CONSTRAINT: Penalize ratios < 1.0 (underprediction) heavily
  #     if(!is.na(current_ratio) && current_ratio < 1.0) {
  #       # Add large penalty for underprediction
  #       prev_ratio_diff[fold, i] <- abs(current_ratio - target_prev_ratio) + (1.0 - current_ratio) * 10
  #     } else {
  #       # Normal calculation for overprediction
  #       prev_ratio_diff[fold, i] <- abs(current_ratio - target_prev_ratio)
  #     }
  #   }
  # 
  #   # Find best threshold for this fold
  #   fold_best_threshold <- test_thresholds[which.min(prev_ratio_diff[fold, ])]
  #   fold_min_diff <- min(prev_ratio_diff[fold, ], na.rm = TRUE)
  # 
  #   # Get actual metrics at this threshold
  #   fold_pred_binary <- as.numeric(presence_prob > fold_best_threshold)
  #   fold_pred_prev <- mean(fold_pred_binary)
  #   fold_actual_ratio <- fold_pred_prev / actual_prev
  # 
  #   cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #   cat("  Achieved ratio:", round(fold_actual_ratio, 3),
  #       "| Target:", round(target_prev_ratio, 3), "\n")
  # }
  # 
  # # Average across folds and find optimal threshold
  # # NOTE - make sure this is in the right spot!!
  # cat("\n--- Aggregating Results ---\n")
  # mean_prev_ratio_diff <- colMeans(prev_ratio_diff, na.rm = TRUE)
  # optimal_threshold <- test_thresholds[which.min(mean_prev_ratio_diff)]
  # min_mean_diff <- min(mean_prev_ratio_diff, na.rm = TRUE)
  # cat("Mean prevalence ratio difference across folds: ", round(min_mean_diff, 6), "\n")
  # cat("Optimal threshold (minimizes difference from target ratio): ", optimal_threshold, "\n\n")
  # 
  # 
  # ################################## OPTION 7: sensitivity ##################################
  # 
  # cat("Optimizing presence threshold to maximize sensitivity via 5-fold cross-validation...\n")
  # 
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # 
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.001, 0.99, by = 0.0001)
  # sensitivity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # 
  # for(fold in 1:5) {
  #   cat("\n--- Fold", fold, "---\n")
  # 
  #   # Split data
  #   test_idx <- which(fold_ids == fold)
  #   train_data <- complete_data[-test_idx, ]
  #   test_data <- complete_data[test_idx, ]
  #   cat("  Training size:", nrow(train_data), "| Test size:", nrow(test_idx), "\n")
  # 
  #   # Fit presence model on training data
  #   presence_fit <- gam(formula(presence_model), data = train_data, family = binomial())
  # 
  #   # Predict on test data
  #   presence_prob <- predict(presence_fit, newdata = test_data, type = "response")
  #   actual_binary <- as.numeric(test_data$cover > 0)
  # 
  #   # Calculate Sensitivity for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  # 
  #     TP <- sum(actual_binary == 1 & predicted_binary == 1)
  #     FN <- sum(actual_binary == 1 & predicted_binary == 0)
  # 
  #     sensitivity <- TP / (TP + FN)
  #     sensitivity_values[fold, i] <- sensitivity
  #   }
  # 
  #   # Find best threshold for this fold
  #   fold_best_threshold <- test_thresholds[which.max(sensitivity_values[fold, ])]
  #   fold_best_sensitivity <- max(sensitivity_values[fold, ], na.rm = TRUE)
  #   cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #   cat("  Maximum Sensitivity:", round(fold_best_sensitivity, 6), "\n")
  # }
  # 
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_sensitivity_values <- colMeans(sensitivity_values, na.rm = TRUE)
  # optimal_threshold <- test_thresholds[which.max(mean_sensitivity_values)]
  # max_mean_sensitivity <- max(mean_sensitivity_values, na.rm = TRUE)
  # 
  # cat("Mean Sensitivity across folds: ", round(max_mean_sensitivity, 6), "\n")
  # cat("Optimal threshold (maximizes Sensitivity): ", optimal_threshold, "\n\n")
  # 
  # ################################## OPTION 8: sens. w/ constraints ##################################
  # 
  # cat("Optimizing presence threshold to maximize sensitivity (max 0.5) with specificity >= 0.8...\n")
  # 
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # 
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # sensitivity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # specificity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # 
  # for(fold in 1:5) {
  #   cat("\n--- Fold", fold, "---\n")
  # 
  #   # Split data
  #   test_idx <- which(fold_ids == fold)
  #   train_data <- complete_data[-test_idx, ]
  #   test_data <- complete_data[test_idx, ]
  #   cat("  Training size:", nrow(train_data), "| Test size:", nrow(test_idx), "\n")
  # 
  #   # Fit presence model on training data
  #   presence_fit <- gam(formula(presence_model), data = train_data, family = binomial())
  # 
  #   # Predict on test data
  #   presence_prob <- predict(presence_fit, newdata = test_data, type = "response")
  #   actual_binary <- as.numeric(test_data$cover > 0)
  # 
  #   # Calculate Sensitivity and Specificity for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  # 
  #     TP <- sum(actual_binary == 1 & predicted_binary == 1)
  #     TN <- sum(actual_binary == 0 & predicted_binary == 0)
  #     FP <- sum(actual_binary == 0 & predicted_binary == 1)
  #     FN <- sum(actual_binary == 1 & predicted_binary == 0)
  # 
  #     sensitivity <- TP / (TP + FN)
  #     specificity <- TN / (TN + FP)
  # 
  #     sensitivity_values[fold, i] <- sensitivity
  #     specificity_values[fold, i] <- specificity
  #   }
  # 
  #   # Find best threshold for this fold (with constraints)
  #   # Only consider thresholds where specificity >= 0.8 and sensitivity <= 0.5
  #   valid_thresholds <- which(specificity_values[fold, ] >= 0.8 & sensitivity_values[fold, ] <= 0.5)
  # 
  #   if(length(valid_thresholds) > 0) {
  #     fold_best_idx <- valid_thresholds[which.max(sensitivity_values[fold, valid_thresholds])]
  #     fold_best_threshold <- test_thresholds[fold_best_idx]
  #     fold_best_sensitivity <- sensitivity_values[fold, fold_best_idx]
  #     fold_best_specificity <- specificity_values[fold, fold_best_idx]
  # 
  #     cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #     cat("  Sensitivity:", round(fold_best_sensitivity, 4), "| Specificity:", round(fold_best_specificity, 4), "\n")
  #   } else {
  #     cat("  WARNING: No threshold met constraints (spec >= 0.8, sens <= 0.5) for this fold\n")
  #   }
  # }
  # 
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_sensitivity_values <- colMeans(sensitivity_values, na.rm = TRUE)
  # mean_specificity_values <- colMeans(specificity_values, na.rm = TRUE)
  # 
  # # Apply constraints: specificity >= 0.8 and sensitivity <= 0.5
  # valid_thresholds <- which(mean_specificity_values >= 0.8 & mean_sensitivity_values <= 0.5)
  # 
  # if(length(valid_thresholds) > 0) {
  #   optimal_idx <- valid_thresholds[which.max(mean_sensitivity_values[valid_thresholds])]
  #   optimal_threshold <- test_thresholds[optimal_idx]
  #   optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
  #   optimal_specificity <- mean_specificity_values[optimal_idx]
  # 
  #   cat("Mean Sensitivity across folds: ", round(optimal_sensitivity, 4), "\n")
  #   cat("Mean Specificity across folds: ", round(optimal_specificity, 4), "\n")
  #   cat("Optimal threshold: ", optimal_threshold, "\n")
  # 
  #   if(optimal_sensitivity >= 0.4) {
  #     cat("✓ Sensitivity >= 0.4 target achieved\n\n")
  #   } else {
  #     cat("⚠ Warning: Sensitivity < 0.4 (target not achieved)\n\n")
  #   }
  # } else {
  #   cat("ERROR: No threshold met constraints (spec >= 0.8, sens <= 0.5)\n")
  #   cat("Falling back to threshold that maximizes balanced accuracy...\n")
  #   optimal_idx <- which.max((mean_sensitivity_values + mean_specificity_values) / 2)
  #   optimal_threshold <- test_thresholds[optimal_idx]
  #   optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
  #   optimal_specificity <- mean_specificity_values[optimal_idx]
  # 
  #   cat("Fallback threshold: ", optimal_threshold, "\n")
  #   cat("Sensitivity: ", round(optimal_sensitivity, 4), "| Specificity: ", round(optimal_specificity, 4), "\n\n")
  # }
  # 
  # ################################## OPTION 9: sens. w/ MORE constraint ##################################
  # 
  # cat("Optimizing presence threshold to maximize sensitivity (max 0.5) with specificity >= 0.90...\n")
  # 
  # upper_constrant = 0.90
  # 
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # 
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # sensitivity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # specificity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # 
  # for(fold in 1:5) {
  #   cat("\n--- Fold", fold, "---\n")
  # 
  #   # Split data
  #   test_idx <- which(fold_ids == fold)
  #   train_data <- complete_data[-test_idx, ]
  #   test_data <- complete_data[test_idx, ]
  #   cat("  Training size:", nrow(train_data), "| Test size:", nrow(test_idx), "\n")
  # 
  #   # Fit presence model on training data
  #   presence_fit <- gam(formula(presence_model), data = train_data, family = binomial())
  # 
  #   # Predict on test data
  #   presence_prob <- predict(presence_fit, newdata = test_data, type = "response")
  #   actual_binary <- as.numeric(test_data$cover > 0)
  # 
  #   # Calculate Sensitivity and Specificity for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  # 
  #     TP <- sum(actual_binary == 1 & predicted_binary == 1)
  #     TN <- sum(actual_binary == 0 & predicted_binary == 0)
  #     FP <- sum(actual_binary == 0 & predicted_binary == 1)
  #     FN <- sum(actual_binary == 1 & predicted_binary == 0)
  # 
  #     sensitivity <- TP / (TP + FN)
  #     specificity <- TN / (TN + FP)
  # 
  #     sensitivity_values[fold, i] <- sensitivity
  #     specificity_values[fold, i] <- specificity
  #   }
  # 
  #   # Find best threshold for this fold (with constraints)
  #   # Only consider thresholds where specificity >= 0.8 and sensitivity <= 0.5
  #   valid_thresholds <- which(specificity_values[fold, ] >= upper_constrant & sensitivity_values[fold, ] <= 0.5)
  # 
  #   if(length(valid_thresholds) > 0) {
  #     fold_best_idx <- valid_thresholds[which.max(sensitivity_values[fold, valid_thresholds])]
  #     fold_best_threshold <- test_thresholds[fold_best_idx]
  #     fold_best_sensitivity <- sensitivity_values[fold, fold_best_idx]
  #     fold_best_specificity <- specificity_values[fold, fold_best_idx]
  # 
  #     cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #     cat("  Sensitivity:", round(fold_best_sensitivity, 4), "| Specificity:", round(fold_best_specificity, 4), "\n")
  #   } else {
  #     cat("  WARNING: No threshold met constraints (spec >= 0.8, sens <= 0.5) for this fold\n")
  #   }
  # }
  # 
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_sensitivity_values <- colMeans(sensitivity_values, na.rm = TRUE)
  # mean_specificity_values <- colMeans(specificity_values, na.rm = TRUE)
  # 
  # # Apply constraints: specificity >= 0.8 and sensitivity <= 0.5
  # valid_thresholds <- which(mean_specificity_values >= upper_constrant & mean_sensitivity_values <= 0.5)
  # 
  # if(length(valid_thresholds) > 0) {
  #   optimal_idx <- valid_thresholds[which.max(mean_sensitivity_values[valid_thresholds])]
  #   optimal_threshold <- test_thresholds[optimal_idx]
  #   optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
  #   optimal_specificity <- mean_specificity_values[optimal_idx]
  # 
  #   cat("Mean Sensitivity across folds: ", round(optimal_sensitivity, 4), "\n")
  #   cat("Mean Specificity across folds: ", round(optimal_specificity, 4), "\n")
  #   cat("Optimal threshold: ", optimal_threshold, "\n")
  # 
  #   if(optimal_sensitivity >= 0.4) {
  #     cat("✓ Sensitivity >= 0.4 target achieved\n\n")
  #   } else {
  #     cat("⚠ Warning: Sensitivity < 0.4 (target not achieved)\n\n")
  #   }
  # } else {
  #   cat("ERROR: No threshold met constraints (spec >= 0.8, sens <= 0.5)\n")
  #   cat("Falling back to threshold that maximizes balanced accuracy...\n")
  #   optimal_idx <- which.max((mean_sensitivity_values + mean_specificity_values) / 2)
  #   optimal_threshold <- test_thresholds[optimal_idx]
  #   optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
  #   optimal_specificity <- mean_specificity_values[optimal_idx]
  # 
  #   cat("Fallback threshold: ", optimal_threshold, "\n")
  #   cat("Sensitivity: ", round(optimal_sensitivity, 4), "| Specificity: ", round(optimal_specificity, 4), "\n\n")
  # }
  # 
  # 
  # ################################## OPTION 9: sens. w/ EVEN MORE constraint ##################################
  # 
  # cat("Optimizing presence threshold to maximize sensitivity (max 0.5) with specificity >= 0.95...\n")
  # 
  # upper_constrant = 0.95
  # 
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # 
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # sensitivity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # specificity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # 
  # for(fold in 1:5) {
  #   cat("\n--- Fold", fold, "---\n")
  # 
  #   # Split data
  #   test_idx <- which(fold_ids == fold)
  #   train_data <- complete_data[-test_idx, ]
  #   test_data <- complete_data[test_idx, ]
  #   cat("  Training size:", nrow(train_data), "| Test size:", nrow(test_idx), "\n")
  # 
  #   # Fit presence model on training data
  #   presence_fit <- gam(formula(presence_model), data = train_data, family = binomial())
  # 
  #   # Predict on test data
  #   presence_prob <- predict(presence_fit, newdata = test_data, type = "response")
  #   actual_binary <- as.numeric(test_data$cover > 0)
  # 
  #   # Calculate Sensitivity and Specificity for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  # 
  #     TP <- sum(actual_binary == 1 & predicted_binary == 1)
  #     TN <- sum(actual_binary == 0 & predicted_binary == 0)
  #     FP <- sum(actual_binary == 0 & predicted_binary == 1)
  #     FN <- sum(actual_binary == 1 & predicted_binary == 0)
  # 
  #     sensitivity <- TP / (TP + FN)
  #     specificity <- TN / (TN + FP)
  # 
  #     sensitivity_values[fold, i] <- sensitivity
  #     specificity_values[fold, i] <- specificity
  #   }
  # 
  #   # Find best threshold for this fold (with constraints)
  #   # Only consider thresholds where specificity >= 0.8 and sensitivity <= 0.5
  #   valid_thresholds <- which(specificity_values[fold, ] >= upper_constrant & sensitivity_values[fold, ] <= 0.5)
  # 
  #   if(length(valid_thresholds) > 0) {
  #     fold_best_idx <- valid_thresholds[which.max(sensitivity_values[fold, valid_thresholds])]
  #     fold_best_threshold <- test_thresholds[fold_best_idx]
  #     fold_best_sensitivity <- sensitivity_values[fold, fold_best_idx]
  #     fold_best_specificity <- specificity_values[fold, fold_best_idx]
  # 
  #     cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #     cat("  Sensitivity:", round(fold_best_sensitivity, 4), "| Specificity:", round(fold_best_specificity, 4), "\n")
  #   } else {
  #     cat("  WARNING: No threshold met constraints (spec >= 0.8, sens <= 0.5) for this fold\n")
  #   }
  # }
  # 
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_sensitivity_values <- colMeans(sensitivity_values, na.rm = TRUE)
  # mean_specificity_values <- colMeans(specificity_values, na.rm = TRUE)
  # 
  # # Apply constraints: specificity >= 0.8 and sensitivity <= 0.5
  # valid_thresholds <- which(mean_specificity_values >= upper_constrant & mean_sensitivity_values <= 0.5)
  # 
  # if(length(valid_thresholds) > 0) {
  #   optimal_idx <- valid_thresholds[which.max(mean_sensitivity_values[valid_thresholds])]
  #   optimal_threshold <- test_thresholds[optimal_idx]
  #   optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
  #   optimal_specificity <- mean_specificity_values[optimal_idx]
  # 
  #   cat("Mean Sensitivity across folds: ", round(optimal_sensitivity, 4), "\n")
  #   cat("Mean Specificity across folds: ", round(optimal_specificity, 4), "\n")
  #   cat("Optimal threshold: ", optimal_threshold, "\n")
  # 
  #   if(optimal_sensitivity >= 0.4) {
  #     cat("✓ Sensitivity >= 0.4 target achieved\n\n")
  #   } else {
  #     cat("⚠ Warning: Sensitivity < 0.4 (target not achieved)\n\n")
  #   }
  # } else {
  #   cat("ERROR: No threshold met constraints (spec >= 0.8, sens <= 0.5)\n")
  #   cat("Falling back to threshold that maximizes balanced accuracy...\n")
  #   optimal_idx <- which.max((mean_sensitivity_values + mean_specificity_values) / 2)
  #   optimal_threshold <- test_thresholds[optimal_idx]
  #   optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
  #   optimal_specificity <- mean_specificity_values[optimal_idx]
  # 
  #   cat("Fallback threshold: ", optimal_threshold, "\n")
  #   cat("Sensitivity: ", round(optimal_sensitivity, 4), "| Specificity: ", round(optimal_specificity, 4), "\n\n")
  # }
  # 
  # ################################## OPTION 10: sens. w/ LESS constraint ##################################
  # 
  # cat("Optimizing presence threshold to maximize sensitivity (max 0.6) with specificity >= 0.80...\n")
  # 
  # upper_constrant = 0.80
  # sens_upper_constraint = 0.60
  # 
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # 
  # # 5-fold cross-validation
  # # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
  # sensitivity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # specificity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
  # 
  # for(fold in 1:5) {
  #   cat("\n--- Fold", fold, "---\n")
  # 
  #   # Split data
  #   test_idx <- which(fold_ids == fold)
  #   train_data <- complete_data[-test_idx, ]
  #   test_data <- complete_data[test_idx, ]
  #   cat("  Training size:", nrow(train_data), "| Test size:", nrow(test_idx), "\n")
  # 
  #   # Fit presence model on training data
  #   presence_fit <- gam(formula(presence_model), data = train_data, family = binomial())
  # 
  #   # Predict on test data
  #   presence_prob <- predict(presence_fit, newdata = test_data, type = "response")
  #   actual_binary <- as.numeric(test_data$cover > 0)
  # 
  #   # Calculate Sensitivity and Specificity for each threshold
  #   for(i in seq_along(test_thresholds)) {
  #     predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
  # 
  #     TP <- sum(actual_binary == 1 & predicted_binary == 1)
  #     TN <- sum(actual_binary == 0 & predicted_binary == 0)
  #     FP <- sum(actual_binary == 0 & predicted_binary == 1)
  #     FN <- sum(actual_binary == 1 & predicted_binary == 0)
  # 
  #     sensitivity <- TP / (TP + FN)
  #     specificity <- TN / (TN + FP)
  # 
  #     sensitivity_values[fold, i] <- sensitivity
  #     specificity_values[fold, i] <- specificity
  #   }
  # 
  #   # Find best threshold for this fold (with constraints)
  #   # Only consider thresholds where specificity >= 0.8 and sensitivity <= 0.5
  #   valid_thresholds <- which(specificity_values[fold, ] >= upper_constrant & sensitivity_values[fold, ] <= sens_upper_constraint)
  # 
  #   if(length(valid_thresholds) > 0) {
  #     fold_best_idx <- valid_thresholds[which.max(sensitivity_values[fold, valid_thresholds])]
  #     fold_best_threshold <- test_thresholds[fold_best_idx]
  #     fold_best_sensitivity <- sensitivity_values[fold, fold_best_idx]
  #     fold_best_specificity <- specificity_values[fold, fold_best_idx]
  # 
  #     cat("  Best threshold for this fold:", fold_best_threshold, "\n")
  #     cat("  Sensitivity:", round(fold_best_sensitivity, 4), "| Specificity:", round(fold_best_specificity, 4), "\n")
  #   } else {
  #     cat("  WARNING: No threshold met constraints (spec >= 0.8, sens <= 0.5) for this fold\n")
  #   }
  # }
  # 
  # # Average across folds and find optimal threshold
  # cat("\n--- Aggregating Results ---\n")
  # mean_sensitivity_values <- colMeans(sensitivity_values, na.rm = TRUE)
  # mean_specificity_values <- colMeans(specificity_values, na.rm = TRUE)
  # 
  # # Apply constraints: specificity >= 0.8 and sensitivity <= 0.5
  # valid_thresholds <- which(mean_specificity_values >= upper_constrant & mean_sensitivity_values <= sens_upper_constraint)
  # 
  # if(length(valid_thresholds) > 0) {
  #   optimal_idx <- valid_thresholds[which.max(mean_sensitivity_values[valid_thresholds])]
  #   optimal_threshold <- test_thresholds[optimal_idx]
  #   optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
  #   optimal_specificity <- mean_specificity_values[optimal_idx]
  # 
  #   cat("Mean Sensitivity across folds: ", round(optimal_sensitivity, 4), "\n")
  #   cat("Mean Specificity across folds: ", round(optimal_specificity, 4), "\n")
  #   cat("Optimal threshold: ", optimal_threshold, "\n")
  # 
  #   if(optimal_sensitivity >= 0.4) {
  #     cat("✓ Sensitivity >= 0.4 target achieved\n\n")
  #   } else {
  #     cat("⚠ Warning: Sensitivity < 0.4 (target not achieved)\n\n")
  #   }
  # } else {
  #   cat("ERROR: No threshold met constraints (spec >= 0.8, sens <= 0.5)\n")
  #   cat("Falling back to threshold that maximizes balanced accuracy...\n")
  #   optimal_idx <- which.max((mean_sensitivity_values + mean_specificity_values) / 2)
  #   optimal_threshold <- test_thresholds[optimal_idx]
  #   optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
  #   optimal_specificity <- mean_specificity_values[optimal_idx]
  # 
  #   cat("Fallback threshold: ", optimal_threshold, "\n")
  #   cat("Sensitivity: ", round(optimal_sensitivity, 4), "| Specificity: ", round(optimal_specificity, 4), "\n\n")
  # }
  # 
  # 
  ################################## in-sample performance ##################################
  
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
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  TSS <- sensitivity + specificity - 1
  actual_prevalence <- mean(full_actual_binary)
  predicted_prevalence <- mean(full_predicted_binary)
  prevalence_diff <- abs(predicted_prevalence - actual_prevalence)
  
  
  full_auc <- as.numeric(auc(full_actual_binary, full_presence_prob))
  
  cat("Threshold:", optimal_threshold, "\n")
  cat("Sensitivity:", round(sensitivity, 3), "| Specificity:", round(specificity, 3), "\n")
  cat("Accuracy:", round(accuracy, 3), "| TSS:", round(TSS, 3), "| AUC:", round(full_auc, 3), "\n")  # MODIFIED THIS LINE
  cat("Actual prevalence:", round(actual_prevalence, 4), "| Predicted prevalence:", round(predicted_prevalence, 4), "\n")
  cat("Prevalence difference:", round(prevalence_diff, 6), "\n\n")
  
  
  
  
  
  
  ################################## train abundance calibration ##################################
  
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
  
  
  ################################## subset predictors across grid ##################################
  
  cat("Preparing spatial predictions (chunked to save memory)...\n")
  
  # Manually specify required variables
  presence_vars <- all.vars(formula(presence_model))[-1]
  abundance_vars <- all.vars(formula(abundance_model))[-1]
  required_vars <- unique(c(presence_vars, abundance_vars))
  
  # Subset raster stack to only required variables
  required_vars <- gsub("^depth_bathy$", "depth", required_vars)
  env_subset <- env_complex[[required_vars]]
  
  # Get dimensions
  n_cells <- ncell(env_subset)
  chunk_size <- 8000000  # Adjust based on your RAM (50k cells = ~5-10 MB per chunk)
  n_chunks <- ceiling(n_cells / chunk_size)
  
  cat("Total cells:", n_cells, "\n")
  cat("Processing in", n_chunks, "chunks of", chunk_size, "cells\n\n")
  
  # Process chunks and combine
  env_df_list <- vector("list", n_chunks)
  
  for(i in 1:n_chunks) {
    start_row <- (i - 1) * chunk_size + 1
    end_row <- min(i * chunk_size, n_cells)
    
    # Extract values for this chunk
    chunk_vals <- values(env_subset, mat = TRUE)[start_row:end_row, , drop = FALSE]
    chunk_coords <- xyFromCell(env_subset, start_row:end_row)
    
    # Combine coordinates and values
    chunk_df <- data.frame(chunk_coords, chunk_vals)
    
    # Remove NAs
    chunk_df <- chunk_df[complete.cases(chunk_df), ]
    
    env_df_list[[i]] <- chunk_df
    
    if(i %% 10 == 0) cat("  Processed chunk", i, "of", n_chunks, "\n")
  }
  
  # Combine all chunks
  env_df <- do.call(rbind, env_df_list)
  rm(env_df_list)  # Free memory
  gc()  # Garbage collection
  
  # Rename depth column to match model expectation
  names(env_df)[names(env_df) == "depth"] <- "depth_bathy"
  
  cat("Final dataset size:", nrow(env_df), "cells\n\n")
  
  # Sample for testing (1/100th of data)
  sample_size <- ceiling(nrow(env_df) / 100)
  # set.seed(123)
  sample_indices <- sample(nrow(env_df), sample_size)
  env_df_sample <- env_df[sample_indices, ]
  
  cat("Species:", species, "\n")
  cat("Full dataset size:", nrow(env_df), "cells\n")
  cat("Sample size (1/100th):", nrow(env_df_sample), "cells\n\n")
  
  # ################################## small-sample test (optional) ##################################
  # 
  # # Presence predictions
  # cat("Predicting presence...\n")
  # start_time <- Sys.time()
  # presence_prob_sample <- predict(presence_model, newdata = env_df_sample, type = "response")
  # presence_time <- difftime(Sys.time(), start_time, units = "secs")
  # cat("Sample presence prediction:", round(presence_time, 2), "seconds\n")
  # 
  # # Apply optimized threshold
  # presence_binary_sample <- ifelse(presence_prob_sample > optimal_threshold, 1, 0)
  # 
  # # Abundance predictions
  # cat("Predicting abundance...\n")
  # start_time <- Sys.time()
  # abundance_pred_sample <- predict(abundance_model, newdata = env_df_sample, type = "response")
  # abundance_time <- difftime(Sys.time(), start_time, units = "secs")
  # cat("Sample abundance prediction:", round(abundance_time, 2), "seconds\n")
  # 
  # # Apply abundance calibration
  # abundance_pred_calibrated <- abundance_correction_fn(abundance_pred_sample)
  # 
  # # Create hurdle predictions with calibrated abundance
  # hurdle_pred_sample <- ifelse(presence_binary_sample == 1, abundance_pred_calibrated, 0)
  # 
  # # Add predictions to sample dataframe
  # env_df_sample$presence_prob <- presence_prob_sample
  # env_df_sample$presence_binary <- presence_binary_sample
  # env_df_sample$abundance_pred_raw <- abundance_pred_sample
  # env_df_sample$abundance_pred_calibrated <- abundance_pred_calibrated
  # env_df_sample$hurdle_pred <- hurdle_pred_sample
  # 
  # # Time estimates
  # total_sample_time <- as.numeric(presence_time + abundance_time)
  # estimated_full_time <- total_sample_time * 100
  # 
  # cat("\nSample results summary:\n")
  # cat("Presence threshold used:", optimal_threshold, "\n")
  # cat("Presence probability range:", round(range(presence_prob_sample), 4), "\n")
  # cat("Abundance prediction range (raw):", round(range(abundance_pred_sample), 6), "\n")
  # cat("Abundance prediction range (calibrated):", round(range(abundance_pred_calibrated), 6), "\n")
  # cat("Hurdle prediction range:", round(range(hurdle_pred_sample), 6), "\n")
  # 
  # cat("\nTime estimates:\n")
  # cat("Sample time:", round(total_sample_time, 2), "seconds\n")
  # cat("Estimated full dataset time:", round(estimated_full_time, 1), "seconds")
  # if(estimated_full_time > 60) {
  #   cat(" (", round(estimated_full_time/60, 1), " minutes)", sep="")
  # }
  # cat("\n\n")
  # 
  # 
  # ################################## sample visualizations (optional) ##################################
  # 
  # # Presence probability map
  # ggplot(env_df_sample, aes(x = x, y = y, color = presence_prob)) +
  #   geom_point(size = 0.5) +
  #   scale_color_viridis_c(name = "Presence\nProbability") +
  #   coord_equal() +
  #   theme_minimal() +
  #   ggtitle(paste(stringr::str_to_title(species), "Presence Probability (Sample)"))
  # 
  # # Histogram of non-zero predictions
  # hist(env_df_sample$hurdle_pred[env_df_sample$hurdle_pred > 0 & env_df_sample$hurdle_pred < 1],
  #      main = "Distribution of Non-Zero Hurdle Predictions (Sample)",
  #      xlab = "Predicted Cover")
  # 
  # # Hurdle prediction map (with calibrated abundance)
  # ggplot(env_df_sample, aes(x = x, y = y, color = hurdle_pred)) +
  #   geom_point(size = 0.5) +
  #   scale_color_viridis_c(name = "Predicted\nCover", limits = c(0, 0.1)) +
  #   coord_equal() +
  #   theme_minimal() +
  #   ggtitle(paste("Hurdle Model Prediction (Calibrated) -", 
  #                 stringr::str_to_title(species), 
  #                 "\nThreshold:", optimal_threshold, "(Sample)"))
  # 
  # 
  ################################## full predictions (parallel) ##################################
  
  cat("\n=== RUNNING FULL MODEL PREDICTIONS (PARALLELIZED) ===\n")
  cat("Warning: This may take 5-10 minutes depending on cores and dataset size\n\n")
  
  start_time <- Sys.time()
  
  # Set up parallel cluster
  n_cores <- detectCores() - 1  # Leave one core free
  cat("Using", n_cores, "cores\n\n")
  cl <- makeCluster(n_cores)
  
  # Split data into chunks
  chunks <- split(1:nrow(env_df), cut(1:nrow(env_df), n_cores))
  
  
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
  

  ################################## prediction summary ##################################
  
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
  
  max_predicted_prev = round(max(hurdle_pred_full), 6)
  
  
  ################################## maps ##################################
  
  cat("Converting predictions to raster...\n")
  
  cover_clamp_val = max_predicted_prev * 100 #in percent
  
  # Convert to raster
  hurdle_raster <- rast(env_df[, c("x", "y", "hurdle_pred")], type = "xyz")
  
  # Apply shallow mask (if bathy_final exists)
  if(exists("bathy_final")) {
    shallow_mask <- bathy_final > -60
    hurdle_raster <- resample(hurdle_raster, bathy_final)
    hurdle_raster <- mask(hurdle_raster, shallow_mask, maskvalues = FALSE)
  }
  
  # Clamp values for visualization
  cover_clamp_val = cover_clamp_val / 100
  hurdle_raster_clamped <- clamp(hurdle_raster, lower = 0, upper = cover_clamp_val, values = TRUE)
  
  cat("Creating interactive leaflet map...\n")
  
  # Create custom color palette matching the paper
  # Blue → Cyan → Yellow → Orange/Red
  paper_colors <- colorRampPalette(c("#0000FF", "#00BFFF", "#00FFFF", 
                                     "#FFFF00", "#FF8C00", "#FF4500"))(256)
  
  # Create leaflet map
  pal <- colorNumeric(paper_colors, values(hurdle_raster_clamped), na.color = "transparent")
  
  # interactive_map <- leaflet() %>%
  #   addTiles(group = "OpenStreetMap") %>%
  #   addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  #   addRasterImage(hurdle_raster_clamped, colors = pal, opacity = 0.7) %>%
  #   addLegend("bottomright", pal = pal, values = values(hurdle_raster_clamped),
  #             title = paste(stringr::str_to_title(species), "Cover (%)"),
  #             labFormat = labelFormat(transform = function(x) x * 100, suffix = "%")) %>%
  #   addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
  #                    options = layersControlOptions(collapsed = FALSE))
  # 
  # # Display the map
  # interactive_map
  # cat("\n=== PREDICTION COMPLETE ===\n")
  # cat("Threshold used:", optimal_threshold, "\n")
  # cat("Calibration applied: Yes\n")
  # cat("Full dataset predicted:", nrow(env_df), "cells\n")
  
  # WITH CORAL COVER
  psu_cover_data <- combined_benthic_data_averaged %>%
    filter(grepl(species, spp, ignore.case = TRUE)) %>%
    group_by(PSU) %>%
    summarise(lat = first(lat), lon = first(lon), 
              avg_cover = mean(cover, na.rm = TRUE), .groups = 'drop')
  
  # Clamp cover data to match cover_clamp_val
  psu_cover_data <- psu_cover_data %>%
    mutate(avg_cover_clamped = pmin(avg_cover / 100, cover_clamp_val))
  
  # Use same paper-style palette for both
  unified_pal <- colorNumeric(paper_colors, 
                              domain = c(0, cover_clamp_val), na.color = "transparent")
  
  interactive_map <- leaflet() %>%
    addTiles(group = "OpenStreetMap") %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addRasterImage(hurdle_raster_clamped, colors = unified_pal, opacity = 0.7) %>%
    addCircleMarkers(data = psu_cover_data, ~lon, ~lat, radius = 8,
                     fillColor = ~unified_pal(avg_cover_clamped), fillOpacity = 0.8,
                     color = "white", weight = 2, stroke = TRUE,
                     popup = ~paste0("<b>PSU:</b> ", PSU, "<br><b>Cover:</b> ", round(avg_cover, 2), "%")) %>%
    addLegend("bottomright", pal = unified_pal, values = values(hurdle_raster_clamped),
              title = paste(stringr::str_to_title(species), "Cover (%)"),
              labFormat = labelFormat(transform = function(x) x * 100, suffix = "%")) %>%
    addLayersControl(baseGroups = c("OpenStreetMap", "Satellite"),
                     options = layersControlOptions(collapsed = FALSE))
  interactive_map
  
  
  ################################## save specific output ##################################
  
  
  predictions_df = env_df[, c("presence_prob", "presence_binary", "abundance_pred_raw", "abundance_pred_calibrated", "hurdle_pred")]
  
  # # heavyweight save; includes raster stack
  # output_list <- list(
  #   map = interactive_map,
  #   raster = hurdle_raster_clamped,
  #   predictions = env_df,
  #   psu_data = psu_cover_data,
  #   threshold = optimal_threshold,
  #   calibration_fn = abundance_correction_fn,
  #   max_cover = cover_clamp_val
  # )
  # 
  # saveRDS(output_list, file = here('output/output_maps', paste0('results_', species, '.rds')))
  
  # lightweight save; includes simply just what is needed for plotting
  output_list <- list(
    map = interactive_map,
    raster = terra::wrap(hurdle_raster_clamped),
    raster_raw = terra::wrap(hurdle_raster),
    predictions = predictions_df,
    threshold = optimal_threshold,
    max_cover = cover_clamp_val
  )

  saveRDS(output_list, file = here('output/output_maps', paste0('results_light_', species, '.rds')))
  
  
  # #to read:
  # results <- readRDS("results_sppX.rds")
  # results$map  # Display map
  # results$predictions  # Access predictions dataframe
  # results$raster  # Access raster
  
  ################################## Save objects/workspace ##################################
  
  # save_new_objects("output/output_maps", existing_objects)
  