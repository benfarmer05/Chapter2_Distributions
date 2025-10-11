  # Species toggle - change this to switch between species
  # species <- "orbicella"
  # species <- "agaricia"
  species <- "porites"
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
  
  
  
  
  
  # # OPTION 2: OPTIMIZE BY TSS (comment out the above block and use this instead)
  # cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # # 5-fold cross-validation
  # set.seed(42)
  # n <- nrow(complete_data)
  # fold_ids <- sample(rep(1:5, length.out = n))
  # test_thresholds <- seq(0.01, 0.99, by = 0.001)
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
  
  
  
  
  
  
  
  
  
  # # OPTION 4: OPTIMIZE BY KAPPA (accounts for chance agreement and class imbalance)
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

  
  
  
  
  
  
  
  # # OPTION 5: OPTIMIZE BY F2-SCORE (favors sensitivity over precision - good for rare species)
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
  
  
  
  
  
  
  # # OPTION 6: OPTIMIZE BY BALANCED ACCURACY (equal weight to sensitivity and specificity)
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
  
  
  
  
  
  # # OPTION 7A: OPTIMIZE TO TARGET SENSITIVITY AND SPECIFICITY (based on best-performing species)
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
  
  
  # # OPTION 7B: OPTIMIZE TO TARGET SENSITIVITY:SPECIFICITY RATIO
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
  
  
  # # OPTION 7C: OPTIMIZE TO TARGET PREDICTED:ACTUAL PREVALENCE RATIO
  # cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # # Prepare complete data
  # all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  # complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  # 
  # # Target prevalence ratio from well-calibrated model (e.g., Orbicella with Kappa)
  # # 0.4412 / 0.398 = 1.109 (slight overprediction)
  # target_prev_ratio <- 1.109
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
  
  
  
  # OPTION 7E: PREVALENCE-ADAPTIVE RATIO WITH MINIMUM CONSTRAINT
  cat("Optimizing presence threshold via 5-fold cross-validation...\n")
  # Prepare complete data
  all_vars <- unique(c(all.vars(formula(presence_model)), all.vars(formula(abundance_model))))
  complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
  
  # Calculate actual prevalence
  actual_prevalence_full <- mean(complete_data$cover > 0)
  
  # Set prevalence-adaptive target ratio (always ≥ 1.0)
  if(actual_prevalence_full >= 0.25) {
    target_prev_ratio <- 1.08  # 8% overprediction for common species
    cat("Species prevalence:", round(actual_prevalence_full, 3), "(common)\n")
  } else if(actual_prevalence_full >= 0.15) {
    target_prev_ratio <- 1.12  # 12% overprediction for moderately common
    cat("Species prevalence:", round(actual_prevalence_full, 3), "(moderately common)\n")
  } else if(actual_prevalence_full >= 0.08) {
    target_prev_ratio <- 1.25  # 25% overprediction for moderately rare
    cat("Species prevalence:", round(actual_prevalence_full, 3), "(moderately rare)\n")
  } else {
    target_prev_ratio <- 1.50  # 50% overprediction for rare species
    cat("Species prevalence:", round(actual_prevalence_full, 3), "(rare)\n")
  }
  
  cat("Target predicted:actual prevalence ratio:", round(target_prev_ratio, 3), "\n")
  cat("(Allows ~", round((target_prev_ratio - 1) * 100, 0), "% overprediction)\n\n", sep="")
  
  # 5-fold cross-validation
  set.seed(42)
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
      
      # CONSTRAINT: Penalize ratios < 1.0 (underprediction) heavily
      if(!is.na(current_ratio) && current_ratio < 1.0) {
        # Add large penalty for underprediction
        prev_ratio_diff[fold, i] <- abs(current_ratio - target_prev_ratio) + (1.0 - current_ratio) * 10
      } else {
        # Normal calculation for overprediction
        prev_ratio_diff[fold, i] <- abs(current_ratio - target_prev_ratio)
      }
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
  cover_clamp_val = 17 #in percent
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