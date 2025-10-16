  
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
  
  ################################## TOGGLE: RUN ALL SPECIES ##################################
  
  # Set to TRUE to run all species automatically, FALSE to run single species
  run_all_species <- TRUE
  
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
  
  # Species toggle - change this to switch between species (only used if run_all_species = FALSE)
  # species <- "agaricia" # N = 651. best is target prev. ratio (OPTION 6C). looks good
  # species <- "colpophyllia" # N = 116. best is sens. w/ constraints (OPTION 8). a bit overpredicted but good
  # species <- "dendrogyra" # N = 35. best is sens. w/ EVEN MORE constraints (OPTION 9 @ 0.95) - looks AWESOME!!
  # species <- "dichocoenia" # N = 76. best is sens. w/ MORE constraints (OPTION 9 @ 0.90), could even use EVEN MORE constraints to ratchet down overprediction
  # species <- "diploria" # N = 120. best is sens. w/ constraints (OPTION 8). a bit overpredicted but good overall
  # species <- "eusmilia" # N = 39. best is sens. w/ EVEN MORE constraints (OPTION 9 @ 0.95). looks excellent
  # species <- "madracis" # N = 131. best is sens. w/ constraints (OPTION 8). great!!
  # species <- "meandrina" # N = 193. best is sens. w/ constraints (OPTION 8). prob a bit overpredicted
  # species <- "montastraea" # N = 490. best is target prev. ratio (OPTION 6C). looks good! maybe a bit overpredicted
  # species <- "mycetophyllia" # N = 28. best is sens. w/ EVEN MORE constraints (OPTION 9 @ 0.95). pretty good, slightly overpredicted
  species <- "orbicella" # N = 795. best is target prev. ratio (OPTION 6C), though overstates extent a bit most likely. Kappa was useful for narrowing toward an ideal prev. ratio in the first place (useful for all common species)
  # species <- "porites" # N = 1016. best is target prev. ratio (OPTION 6C). looks good! maybe underestimates a bit...not sure (definitely in MCD though)
  # species <- "pseudodiploria" # N = 401. best is target prev. ratio (OPTION 6C), but really misses mesophotic pstrig (tried sens. w/ LESS constraint to handle slight underestimation but this just overstated extent in the wrong places)
  # species <- "siderastrea" # N = 845. best is target prev. ratio (OPTION 6C). looks good! probably slightly overstated
  # species <- "solenastrea" # N = 21. tried sens. w/ EVEN MORE constraints (OPTION 9 @ 0.95), but predictions are much too overestimated. may drop this species
  
  # Determine which species to run
  species_to_run <- if(run_all_species) {
    species_list[species_list != "solenastrea"]  # Exclude solenastrea as noted it may be dropped
  } else {
    species
  }
  
  cat("\n=== RUNNING PREDICTIONS FOR:", paste(species_to_run, collapse=", "), "===\n\n")
  
  
  ################################## MAIN LOOP ##################################
  
  for(current_species in species_to_run) {
    
    cat("\n\n")
    cat("================================================================================\n")
    cat("=== PROCESSING SPECIES:", toupper(current_species), "===\n")
    cat("================================================================================\n\n")
    
    species <- current_species
    
    # Extract models for this species
    presence_model <- all_species_models[[species]]$presence_model
    abundance_model <- all_species_models[[species]]$abundance_model
    model_data <- all_species_models[[species]]$model_data
    
    
    
    ###                                                    ###
    ###                                                    ###
    ###     BELOW: options for optimizing P/A threshold    ###
    ###                                                    ###
    ###                                                    ###
    
    ################################## OPTION 6c: target prev. ratio ##################################
    
    # Species using this option: agaricia, montastraea, orbicella, porites, pseudodiploria, siderastrea
    if(species %in% c("agaricia", "montastraea", "orbicella", "porites", "pseudodiploria", "siderastrea")) {
      
      cat("Using OPTION 6C: Target prevalence ratio\n")
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
      
    }
    
    
    ################################## OPTION 8: sens. w/ constraints ##################################
    
    # Species using this option: colpophyllia, diploria, madracis, meandrina
    if(species %in% c("colpophyllia", "diploria", "madracis", "meandrina")) {
      
      cat("Using OPTION 8: Sensitivity with constraints (specificity >= 0.8)\n")
      cat("Optimizing presence threshold to maximize sensitivity (max 0.5) with specificity >= 0.8...\n")
      
      # Prepare complete data
      all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
      complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
      
      # 5-fold cross-validation
      # set.seed(42)
      n <- nrow(complete_data)
      fold_ids <- sample(rep(1:5, length.out = n))
      test_thresholds <- seq(0.01, 0.99, by = 0.001)
      sensitivity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
      specificity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
      
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
        
        # Calculate Sensitivity and Specificity for each threshold
        for(i in seq_along(test_thresholds)) {
          predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
          
          TP <- sum(actual_binary == 1 & predicted_binary == 1)
          TN <- sum(actual_binary == 0 & predicted_binary == 0)
          FP <- sum(actual_binary == 0 & predicted_binary == 1)
          FN <- sum(actual_binary == 1 & predicted_binary == 0)
          
          sensitivity <- TP / (TP + FN)
          specificity <- TN / (TN + FP)
          
          sensitivity_values[fold, i] <- sensitivity
          specificity_values[fold, i] <- specificity
        }
        
        # Find best threshold for this fold (with constraints)
        # Only consider thresholds where specificity >= 0.8 and sensitivity <= 0.5
        valid_thresholds <- which(specificity_values[fold, ] >= 0.8 & sensitivity_values[fold, ] <= 0.5)
        
        if(length(valid_thresholds) > 0) {
          fold_best_idx <- valid_thresholds[which.max(sensitivity_values[fold, valid_thresholds])]
          fold_best_threshold <- test_thresholds[fold_best_idx]
          fold_best_sensitivity <- sensitivity_values[fold, fold_best_idx]
          fold_best_specificity <- specificity_values[fold, fold_best_idx]
          
          cat("  Best threshold for this fold:", fold_best_threshold, "\n")
          cat("  Sensitivity:", round(fold_best_sensitivity, 4), "| Specificity:", round(fold_best_specificity, 4), "\n")
        } else {
          cat("  WARNING: No threshold met constraints (spec >= 0.8, sens <= 0.5) for this fold\n")
        }
      }
      
      # Average across folds and find optimal threshold
      cat("\n--- Aggregating Results ---\n")
      mean_sensitivity_values <- colMeans(sensitivity_values, na.rm = TRUE)
      mean_specificity_values <- colMeans(specificity_values, na.rm = TRUE)
      
      # Apply constraints: specificity >= 0.8 and sensitivity <= 0.5
      valid_thresholds <- which(mean_specificity_values >= 0.8 & mean_sensitivity_values <= 0.5)
      
      if(length(valid_thresholds) > 0) {
        optimal_idx <- valid_thresholds[which.max(mean_sensitivity_values[valid_thresholds])]
        optimal_threshold <- test_thresholds[optimal_idx]
        optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
        optimal_specificity <- mean_specificity_values[optimal_idx]
        
        cat("Mean Sensitivity across folds: ", round(optimal_sensitivity, 4), "\n")
        cat("Mean Specificity across folds: ", round(optimal_specificity, 4), "\n")
        cat("Optimal threshold: ", optimal_threshold, "\n")
        
        if(optimal_sensitivity >= 0.4) {
          cat("✓ Sensitivity >= 0.4 target achieved\n\n")
        } else {
          cat("⚠ Warning: Sensitivity < 0.4 (target not achieved)\n\n")
        }
      } else {
        cat("ERROR: No threshold met constraints (spec >= 0.8, sens <= 0.5)\n")
        cat("Falling back to threshold that maximizes balanced accuracy...\n")
        optimal_idx <- which.max((mean_sensitivity_values + mean_specificity_values) / 2)
        optimal_threshold <- test_thresholds[optimal_idx]
        optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
        optimal_specificity <- mean_specificity_values[optimal_idx]
        
        cat("Fallback threshold: ", optimal_threshold, "\n")
        cat("Sensitivity: ", round(optimal_sensitivity, 4), "| Specificity: ", round(optimal_specificity, 4), "\n\n")
      }
      
    }
    
    
    ################################## OPTION 9 @ 0.90: sens. w/ MORE constraint ##################################
    
    # Species using this option: dichocoenia
    if(species %in% c("dichocoenia")) {
      
      cat("Using OPTION 9 @ 0.90: Sensitivity with MORE constraints (specificity >= 0.90)\n")
      cat("Optimizing presence threshold to maximize sensitivity (max 0.5) with specificity >= 0.90...\n")
      
      upper_constrant = 0.90
      
      # Prepare complete data
      all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
      complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
      
      # 5-fold cross-validation
      # set.seed(42)
      n <- nrow(complete_data)
      fold_ids <- sample(rep(1:5, length.out = n))
      test_thresholds <- seq(0.01, 0.99, by = 0.001)
      sensitivity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
      specificity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
      
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
        
        # Calculate Sensitivity and Specificity for each threshold
        for(i in seq_along(test_thresholds)) {
          predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
          
          TP <- sum(actual_binary == 1 & predicted_binary == 1)
          TN <- sum(actual_binary == 0 & predicted_binary == 0)
          FP <- sum(actual_binary == 0 & predicted_binary == 1)
          FN <- sum(actual_binary == 1 & predicted_binary == 0)
          
          sensitivity <- TP / (TP + FN)
          specificity <- TN / (TN + FP)
          
          sensitivity_values[fold, i] <- sensitivity
          specificity_values[fold, i] <- specificity
        }
        
        # Find best threshold for this fold (with constraints)
        # Only consider thresholds where specificity >= upper_constraint and sensitivity <= 0.5
        valid_thresholds <- which(specificity_values[fold, ] >= upper_constrant & sensitivity_values[fold, ] <= 0.5)
        
        if(length(valid_thresholds) > 0) {
          fold_best_idx <- valid_thresholds[which.max(sensitivity_values[fold, valid_thresholds])]
          fold_best_threshold <- test_thresholds[fold_best_idx]
          fold_best_sensitivity <- sensitivity_values[fold, fold_best_idx]
          fold_best_specificity <- specificity_values[fold, fold_best_idx]
          
          cat("  Best threshold for this fold:", fold_best_threshold, "\n")
          cat("  Sensitivity:", round(fold_best_sensitivity, 4), "| Specificity:", round(fold_best_specificity, 4), "\n")
        } else {
          cat("  WARNING: No threshold met constraints for this fold\n")
        }
      }
      
      # Average across folds and find optimal threshold
      cat("\n--- Aggregating Results ---\n")
      mean_sensitivity_values <- colMeans(sensitivity_values, na.rm = TRUE)
      mean_specificity_values <- colMeans(specificity_values, na.rm = TRUE)
      
      # Apply constraints
      valid_thresholds <- which(mean_specificity_values >= upper_constrant & mean_sensitivity_values <= 0.5)
      
      if(length(valid_thresholds) > 0) {
        optimal_idx <- valid_thresholds[which.max(mean_sensitivity_values[valid_thresholds])]
        optimal_threshold <- test_thresholds[optimal_idx]
        optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
        optimal_specificity <- mean_specificity_values[optimal_idx]
        
        cat("Mean Sensitivity across folds: ", round(optimal_sensitivity, 4), "\n")
        cat("Mean Specificity across folds: ", round(optimal_specificity, 4), "\n")
        cat("Optimal threshold: ", optimal_threshold, "\n")
        
        if(optimal_sensitivity >= 0.4) {
          cat("✓ Sensitivity >= 0.4 target achieved\n\n")
        } else {
          cat("⚠ Warning: Sensitivity < 0.4 (target not achieved)\n\n")
        }
      } else {
        cat("ERROR: No threshold met constraints\n")
        cat("Falling back to threshold that maximizes balanced accuracy...\n")
        optimal_idx <- which.max((mean_sensitivity_values + mean_specificity_values) / 2)
        optimal_threshold <- test_thresholds[optimal_idx]
        optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
        optimal_specificity <- mean_specificity_values[optimal_idx]
        
        cat("Fallback threshold: ", optimal_threshold, "\n")
        cat("Sensitivity: ", round(optimal_sensitivity, 4), "| Specificity: ", round(optimal_specificity, 4), "\n\n")
      }
      
    }
    
    
    ################################## OPTION 9 @ 0.95: sens. w/ EVEN MORE constraint ##################################
    
    # Species using this option: dendrogyra, eusmilia, mycetophyllia
    if(species %in% c("dendrogyra", "eusmilia", "mycetophyllia")) {
      
      cat("Using OPTION 9 @ 0.95: Sensitivity with EVEN MORE constraints (specificity >= 0.95)\n")
      cat("Optimizing presence threshold to maximize sensitivity (max 0.5) with specificity >= 0.95...\n")
      
      upper_constrant = 0.95
      
      # Prepare complete data
      all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
      complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
      
      # 5-fold cross-validation
      # set.seed(42)
      n <- nrow(complete_data)
      fold_ids <- sample(rep(1:5, length.out = n))
      test_thresholds <- seq(0.01, 0.99, by = 0.001)
      sensitivity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
      specificity_values <- matrix(NA, nrow = 5, ncol = length(test_thresholds))
      
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
        
        # Calculate Sensitivity and Specificity for each threshold
        for(i in seq_along(test_thresholds)) {
          predicted_binary <- as.numeric(presence_prob > test_thresholds[i])
          
          TP <- sum(actual_binary == 1 & predicted_binary == 1)
          TN <- sum(actual_binary == 0 & predicted_binary == 0)
          FP <- sum(actual_binary == 0 & predicted_binary == 1)
          FN <- sum(actual_binary == 1 & predicted_binary == 0)
          
          sensitivity <- TP / (TP + FN)
          specificity <- TN / (TN + FP)
          
          sensitivity_values[fold, i] <- sensitivity
          specificity_values[fold, i] <- specificity
        }
        
        # Find best threshold for this fold (with constraints)
        valid_thresholds <- which(specificity_values[fold, ] >= upper_constrant & sensitivity_values[fold, ] <= 0.5)
        
        if(length(valid_thresholds) > 0) {
          fold_best_idx <- valid_thresholds[which.max(sensitivity_values[fold, valid_thresholds])]
          fold_best_threshold <- test_thresholds[fold_best_idx]
          fold_best_sensitivity <- sensitivity_values[fold, fold_best_idx]
          fold_best_specificity <- specificity_values[fold, fold_best_idx]
          
          cat("  Best threshold for this fold:", fold_best_threshold, "\n")
          cat("  Sensitivity:", round(fold_best_sensitivity, 4), "| Specificity:", round(fold_best_specificity, 4), "\n")
        } else {
          cat("  WARNING: No threshold met constraints for this fold\n")
        }
      }
      
      # Average across folds and find optimal threshold
      cat("\n--- Aggregating Results ---\n")
      mean_sensitivity_values <- colMeans(sensitivity_values, na.rm = TRUE)
      mean_specificity_values <- colMeans(specificity_values, na.rm = TRUE)
      
      # Apply constraints
      valid_thresholds <- which(mean_specificity_values >= upper_constrant & mean_sensitivity_values <= 0.5)
      
      if(length(valid_thresholds) > 0) {
        optimal_idx <- valid_thresholds[which.max(mean_sensitivity_values[valid_thresholds])]
        optimal_threshold <- test_thresholds[optimal_idx]
        optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
        optimal_specificity <- mean_specificity_values[optimal_idx]
        
        cat("Mean Sensitivity across folds: ", round(optimal_sensitivity, 4), "\n")
        cat("Mean Specificity across folds: ", round(optimal_specificity, 4), "\n")
        cat("Optimal threshold: ", optimal_threshold, "\n")
        
        if(optimal_sensitivity >= 0.4) {
          cat("✓ Sensitivity >= 0.4 target achieved\n\n")
        } else {
          cat("⚠ Warning: Sensitivity < 0.4 (target not achieved)\n\n")
        }
      } else {
        cat("ERROR: No threshold met constraints\n")
        cat("Falling back to threshold that maximizes balanced accuracy...\n")
        optimal_idx <- which.max((mean_sensitivity_values + mean_specificity_values) / 2)
        optimal_threshold <- test_thresholds[optimal_idx]
        optimal_sensitivity <- mean_sensitivity_values[optimal_idx]
        optimal_specificity <- mean_specificity_values[optimal_idx]
        
        cat("Fallback threshold: ", optimal_threshold, "\n")
        cat("Sensitivity: ", round(optimal_sensitivity, 4), "| Specificity: ", round(optimal_specificity, 4), "\n\n")
      }
      
    }
    
    
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
    cat("Accuracy:", round(accuracy, 3), "| TSS:", round(TSS, 3), "| AUC:", round(full_auc, 3), "\n")
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
    
    if(!run_all_species) {
      print(interactive_map)  # Only display map if running single species
    }
    
    
    ################################## save specific output ##################################
    
    
    predictions_df = env_df[, c("presence_prob", "presence_binary", "abundance_pred_raw", "abundance_pred_calibrated", "hurdle_pred")]
    
    # lightweight save with terra::wrap() to preserve raster integrity
    output_list <- list(
      map = interactive_map,
      raster = terra::wrap(hurdle_raster_clamped),
      raster_raw = terra::wrap(hurdle_raster),
      predictions = predictions_df,
      threshold = optimal_threshold,
      max_cover = cover_clamp_val
    )
    
    saveRDS(output_list, file = here('output/output_maps', paste0('results_light_', species, '.rds')))
    
    cat("\n✓ Saved results for", species, "\n")
    
  }  # End of species loop
  
  cat("\n\n")
  cat("================================================================================\n")
  cat("=== ALL PROCESSING COMPLETE ===\n")
  cat("================================================================================\n")