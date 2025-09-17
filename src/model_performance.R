 
   # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(terra)
  library(mgcv)
  library(cmocean)
  library(ggplot2)
  library(gridExtra)
  library(tidyverse)
  library(gratia)
  library(tictoc)
  library(pROC)
  
  # library(gam.hp)
  
  source(here("src/functions.R"))
  
  
  ################################## TESTING, AND SHOWING PROMISE ##################################
  
  # COMPLETE BETA vs GAMMA HURDLE VALIDATION FUNCTION - Copy and paste this entire block
  
  validate_hurdle_simple <- function(species_name, train_percent = 80, seed = 300, presence_threshold = 0.5, 
                                     auto_plot = TRUE, use_beta = FALSE) {
    
    library(mgcv)
    library(pROC)
    
    set.seed(seed)
    
    # Get models based on distribution choice
    all_objects <- ls(envir = .GlobalEnv)
    species_gam_objects <- grep(paste0("^", species_name, "_gam_"), all_objects, value = TRUE)
    
    presence_model_name <- grep("_presence_binom$", species_gam_objects, value = TRUE)
    
    if(use_beta) {
      abundance_model_name <- grep("_abundance_beta$", species_gam_objects, value = TRUE)
      cat("Using BETA distribution for abundance model\n")
    } else {
      abundance_model_name <- grep("_abundance_gamma$", species_gam_objects, value = TRUE)
      cat("Using GAMMA distribution for abundance model\n")
    }
    
    if(length(presence_model_name) == 0 || length(abundance_model_name) == 0) {
      missing_models <- c()
      if(length(presence_model_name) == 0) missing_models <- c(missing_models, paste0(species_name, "_gam_presence_binom"))
      if(use_beta && length(abundance_model_name) == 0) missing_models <- c(missing_models, paste0(species_name, "_gam_abundance_beta"))
      if(!use_beta && length(abundance_model_name) == 0) missing_models <- c(missing_models, paste0(species_name, "_gam_abundance_gamma"))
      stop("Could not find required models: ", paste(missing_models, collapse = ", "))
    }
    
    presence_model <- get(presence_model_name, envir = .GlobalEnv)
    abundance_model <- get(abundance_model_name, envir = .GlobalEnv)
    
    # Get data
    data_name <- paste0(species_name, "_model_data")
    if(data_name %in% all_objects) {
      model_data <- get(data_name, envir = .GlobalEnv)
    } else {
      model_data <- presence_model$model
    }
    
    # Get complete data
    presence_vars <- all.vars(formula(presence_model))
    abundance_vars <- all.vars(formula(abundance_model))
    all_vars <- unique(c(presence_vars, abundance_vars))
    complete_data <- model_data[complete.cases(model_data[, all_vars]), ]
    
    # Train-test split
    n <- nrow(complete_data)
    train_size <- floor((train_percent / 100) * n)
    train_idx <- sample(1:n, size = train_size)
    
    train_data <- complete_data[train_idx, ]
    test_data <- complete_data[-train_idx, ]
    
    cat("Training on", nrow(train_data), "observations\n")
    cat("Testing on", nrow(test_data), "observations\n")
    cat("Train/test split:", train_percent, "/", 100 - train_percent, "\n")
    cat("Presence threshold:", presence_threshold, "\n")
    cat("Random seed:", seed, "\n")
    cat("Abundance model distribution:", ifelse(use_beta, "Beta", "Gamma"), "\n")
    
    # Fit models
    presence_fit <- gam(formula(presence_model), data = train_data, family = binomial())
    
    # Prepare abundance training data and fit abundance model
    abundance_train_data <- train_data[train_data$cover > 0, ]
    
    if(use_beta) {
      # For beta models, need to work with proportions
      abundance_train_data$cover_prop <- abundance_train_data$cover / 100
      abundance_formula_adj <- update(formula(abundance_model), cover_prop ~ .)
      abundance_fit <- gam(abundance_formula_adj, data = abundance_train_data, family = betar())
    } else {
      # For gamma models, work with original cover scale
      abundance_fit <- gam(formula(abundance_model), data = abundance_train_data, family = Gamma(link = "log"))
    }
    
    # Make predictions
    presence_prob <- predict(presence_fit, newdata = test_data, type = "response")
    
    if(use_beta) {
      # Beta model predictions are on 0-1 scale, convert to percentage
      abundance_pred_prop <- predict(abundance_fit, newdata = test_data, type = "response")
      abundance_pred <- abundance_pred_prop * 100  # Convert back to 0-100 scale
    } else {
      # Gamma model predictions are already on 0-100 scale
      abundance_pred <- predict(abundance_fit, newdata = test_data, type = "response")
    }
    
    # HURDLE MODEL: Binary presence threshold, then multiply
    presence_binary <- ifelse(presence_prob > presence_threshold, 1, 0)
    hurdle_predictions <- presence_binary * abundance_pred  # Only abundance where predicted present
    
    actual <- test_data$cover
    actual_binary <- ifelse(actual > 0, 1, 0)
    
    # Sample sizes
    n_total <- length(actual)
    n_actual_present <- sum(actual > 0)
    n_predicted_present <- sum(presence_binary == 1)
    n_both_present <- sum(actual > 0 & presence_binary == 1)
    n_actual_absent <- sum(actual == 0)
    
    # Overall metrics
    overall_r2 <- cor(actual, hurdle_predictions, use = "complete.obs")^2
    overall_rmse <- sqrt(mean((actual - hurdle_predictions)^2))
    overall_mae <- mean(abs(actual - hurdle_predictions))
    
    # Absence prediction metrics (renamed from "zero" to "absence")
    actual_absences <- sum(actual == 0)
    predicted_absences <- sum(hurdle_predictions == 0)
    correct_absences <- sum(actual == 0 & hurdle_predictions == 0)
    absence_accuracy <- if(actual_absences > 0) correct_absences / actual_absences else NA
    
    # Presence prediction accuracy metrics
    actual_presences <- sum(actual > 0)
    predicted_presences <- sum(presence_binary == 1)
    correct_presences <- sum(actual > 0 & presence_binary == 1)  # True positives
    false_positives <- sum(actual == 0 & presence_binary == 1)   # Predicted present but actually absent
    false_negatives <- sum(actual > 0 & presence_binary == 0)   # Predicted absent but actually present
    true_negatives <- sum(actual == 0 & presence_binary == 0)   # Correctly predicted absent
    
    # Calculate standard classification metrics
    presence_accuracy <- (correct_presences + true_negatives) / n_total  # Overall accuracy
    presence_sensitivity <- if(actual_presences > 0) correct_presences / actual_presences else NA  # True positive rate
    presence_specificity <- if(actual_absences > 0) true_negatives / actual_absences else NA  # True negative rate
    presence_precision <- if(predicted_presences > 0) correct_presences / predicted_presences else NA  # Positive predictive value
    
    # Method A: All actual presence sites
    method_A_r2 <- NA
    method_A_mae <- NA
    if(n_actual_present > 5) {
      actual_present <- actual[actual > 0]
      pred_present <- abundance_pred[actual > 0]
      method_A_r2 <- cor(actual_present, pred_present, use = "complete.obs")^2
      method_A_mae <- mean(abs(actual_present - pred_present))
    }
    
    # Method B: Both actual and predicted present
    method_B_r2 <- NA
    method_B_mae <- NA
    if(n_both_present > 5) {
      mask <- actual > 0 & presence_binary == 1
      actual_both <- actual[mask]
      pred_both <- abundance_pred[mask]
      method_B_r2 <- cor(actual_both, pred_both, use = "complete.obs")^2
      method_B_mae <- mean(abs(actual_both - pred_both))
    }
    
    # Method C: All predicted present sites (regardless of actual)
    method_C_r2 <- NA
    method_C_mae <- NA
    if(n_predicted_present > 5) {
      mask_C <- presence_binary == 1
      actual_pred_present <- actual[mask_C]
      pred_pred_present <- abundance_pred[mask_C]
      method_C_r2 <- cor(actual_pred_present, pred_pred_present, use = "complete.obs")^2
      method_C_mae <- mean(abs(actual_pred_present - pred_pred_present))
    }
    
    # Presence AUC
    presence_auc <- as.numeric(auc(actual > 0, presence_prob))
    
    # Print results
    cat("\n=== RESULTS ===\n")
    cat("Sample sizes:\n")
    cat("  Total test sites:", n_total, "\n")
    cat("  Actual presence sites:", n_actual_present, "\n")
    cat("  Predicted presence sites:", n_predicted_present, "\n")
    cat("  Both actual AND predicted present:", n_both_present, "\n")
    cat("  Actual absence sites:", n_actual_absent, "\n\n")
    
    cat("Overall performance (all sites including absences):\n")
    cat("  R²:", round(overall_r2, 3), "\n")
    cat("  RMSE:", round(overall_rmse, 3), "\n")
    cat("  MAE:", round(overall_mae, 3), "\n\n")
    
    cat("Absence prediction performance:\n")
    cat("  Actual absences:", actual_absences, "/ Predicted absences:", predicted_absences, "\n")
    cat("  Correctly predicted absences:", correct_absences, "\n")
    cat("  Absence prediction accuracy:", round(absence_accuracy, 3), "\n\n")
    
    cat("Presence prediction performance:\n")
    cat("  Overall accuracy:", round(presence_accuracy, 3), "\n")
    cat("  Sensitivity (true positive rate):", round(presence_sensitivity, 3), "\n")
    cat("  Specificity (true negative rate):", round(presence_specificity, 3), "\n")
    cat("  Precision (positive predictive value):", round(presence_precision, 3), "\n")
    cat("  AUC:", round(presence_auc, 3), "\n")
    cat("  Presence threshold used:", presence_threshold, "\n\n")
    
    cat("Abundance model - Method A (all actual presence sites):\n")
    if(!is.na(method_A_r2)) {
      cat("  R²:", round(method_A_r2, 3), "\n")
      cat("  MAE:", round(method_A_mae, 3), "\n")
      cat("  Sites used:", n_actual_present, "\n\n")
    } else {
      cat("  Cannot calculate (insufficient sites)\n\n")
    }
    
    cat("Abundance model - Method B (both actual and predicted present):\n")
    if(!is.na(method_B_r2)) {
      cat("  R²:", round(method_B_r2, 3), "\n")
      cat("  MAE:", round(method_B_mae, 3), "\n")
      cat("  Sites used:", n_both_present, "\n\n")
    } else {
      cat("  Cannot calculate (insufficient sites)\n\n")
    }
    
    cat("Abundance model - Method C (all predicted present sites):\n")
    if(!is.na(method_C_r2)) {
      cat("  R²:", round(method_C_r2, 3), "\n")
      cat("  MAE:", round(method_C_mae, 3), "\n")
      cat("  Sites used:", n_predicted_present, "\n\n")
    } else {
      cat("  Cannot calculate (insufficient sites)\n\n")
    }
    
    # Create test results for mapping
    test_results <- test_data
    test_results$predicted <- hurdle_predictions
    test_results$presence_prob <- presence_prob
    test_results$presence_binary <- presence_binary
    test_results$abundance_pred <- abundance_pred
    test_results$residuals <- actual - hurdle_predictions
    test_results$actual_binary <- actual_binary
    
    result <- list(
      test_results = test_results,
      species = species_name,
      response_variable = "cover",
      train_percent = train_percent,
      seed = seed,
      presence_threshold = presence_threshold,
      use_beta = use_beta,
      abundance_distribution = ifelse(use_beta, "Beta", "Gamma"),
      n_total = n_total,
      n_actual_present = n_actual_present,
      n_predicted_present = n_predicted_present,
      n_both_present = n_both_present,
      overall_r2 = overall_r2,
      overall_mae = overall_mae,
      actual_absences = actual_absences,
      predicted_absences = predicted_absences,
      correct_absences = correct_absences,
      absence_accuracy = absence_accuracy,
      presence_accuracy = presence_accuracy,
      presence_sensitivity = presence_sensitivity,
      presence_specificity = presence_specificity,
      presence_precision = presence_precision,
      method_A_r2 = method_A_r2,
      method_A_mae = method_A_mae,
      method_B_r2 = method_B_r2,
      method_B_mae = method_B_mae,
      method_C_r2 = method_C_r2,
      method_C_mae = method_C_mae,
      presence_auc = presence_auc
    )
    
    # Auto-plot if requested
    if(auto_plot) {
      plot_hurdle_complete(result)
    }
    
    return(result)
  }
  
  # Complete plotting function with all 8 plots and shared axis limits
  plot_hurdle_complete <- function(result) {
    library(ggplot2)
    library(viridis)
    library(gridExtra)
    
    plot_data <- result$test_results
    
    # Find coordinate columns
    coord_cols <- c("lon", "longitude", "lat", "latitude")
    available_coords <- coord_cols[coord_cols %in% colnames(plot_data)]
    
    if(length(available_coords) < 2) {
      stop("Could not find coordinate columns. Available: ", paste(colnames(plot_data), collapse = ", "))
    }
    
    x_col <- available_coords[grepl("lon", available_coords)][1]
    y_col <- available_coords[grepl("lat", available_coords)][1]
    
    cat("Using coordinates:", x_col, "and", y_col, "\n")
    
    # Common theme for maps
    map_theme <- theme_minimal() +
      theme(
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 9, hjust = 0.5),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      )
    
    # Calculate shared color scale for cover values
    actual_values <- plot_data$cover
    predicted_values <- plot_data$predicted
    shared_min <- min(c(actual_values, predicted_values), na.rm = TRUE)
    shared_max <- max(c(actual_values, predicted_values), na.rm = TRUE)
    
    
    # SPATIAL MAPS (Top two rows)
    # Plot 1: Actual cover values
    p1 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = "cover")) +
      geom_point(size = 1.5, alpha = 0.7) +
      scale_color_viridis_c(name = "Actual", option = "plasma", limits = c(shared_min, shared_max)) +
      labs(title = "Actual Cover", x = "Longitude", y = "Latitude") +
      map_theme
    
    # Plot 2: Predicted cover values
    p2 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = "predicted")) +
      geom_point(size = 1.5, alpha = 0.7) +
      scale_color_viridis_c(name = "Predicted", option = "plasma", limits = c(shared_min, shared_max)) +
      labs(title = "Predicted Cover", x = "Longitude", y = "Latitude") +
      map_theme
    
    # Plot 3: Observed presence/absence (binary)
    p3 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = "factor(actual_binary)")) +
      geom_point(size = 1.5, alpha = 0.7) +
      scale_color_manual(name = "Observed", values = c("0" = "red", "1" = "blue"), 
                         labels = c("0" = "Absent", "1" = "Present")) +
      labs(title = "Observed Presence/Absence", x = "Longitude", y = "Latitude") +
      map_theme
    
    # Plot 4: Predicted presence probability (continuous)
    p4 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = "presence_prob")) +
      geom_point(size = 1.5, alpha = 0.7) +
      scale_color_viridis_c(name = "P(Presence)", option = "viridis") +
      labs(title = "Predicted Presence Probability", x = "Longitude", y = "Latitude") +
      map_theme
    
    # Plot 5: Predicted presence/absence (binary, based on threshold)
    p5 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = "factor(presence_binary)")) +
      geom_point(size = 1.5, alpha = 0.7) +
      scale_color_manual(name = "Predicted", values = c("0" = "red", "1" = "blue"), 
                         labels = c("0" = "Absent", "1" = "Present")) +
      labs(title = paste("Predicted P/A (threshold =", result$presence_threshold, ")"), 
           x = "Longitude", y = "Latitude") +
      map_theme
    
    # ABUNDANCE SCATTER PLOTS (Bottom row) - WITH SHARED AXIS LIMITS
    # Prepare data for abundance scatter plots
    # Method A: All actual presence sites
    method_A_data <- plot_data[plot_data$cover > 0, ]
    
    # Method B: Both actual and predicted present
    method_B_data <- plot_data[plot_data$cover > 0 & plot_data$presence_binary == 1, ]
    
    # Method C: All predicted present sites
    method_C_data <- plot_data[plot_data$presence_binary == 1, ]
    
    # Calculate shared axis limits for all three scatter plots
    all_actual_values <- c(method_A_data$cover, method_B_data$cover, method_C_data$cover)
    all_pred_values <- c(method_A_data$abundance_pred, method_B_data$abundance_pred, method_C_data$abundance_pred)
    
    scatter_min <- min(c(all_actual_values, all_pred_values), na.rm = TRUE)
    scatter_max <- max(c(all_actual_values, all_pred_values), na.rm = TRUE)
    
    # Add a small buffer to axis limits
    axis_buffer <- (scatter_max - scatter_min) * 0.05
    axis_min <- scatter_min - axis_buffer
    axis_max <- scatter_max + axis_buffer
    
    # Plot 6: Method A scatter plot
    p6 <- ggplot(method_A_data, aes(x = abundance_pred, y = cover)) +
      geom_point(alpha = 0.6, color = "darkgreen", size = 1) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      xlim(axis_min, axis_max) +
      ylim(axis_min, axis_max) +
      labs(title = paste("Method A: Actual Present\n(R² =", round(result$method_A_r2, 3), ", n =", nrow(method_A_data), ")"),
           x = "Abundance Prediction", y = "Actual Cover") +
      theme_minimal() +
      theme(plot.title = element_text(size = 8))
    
    # Plot 7: Method B scatter plot  
    p7 <- ggplot(method_B_data, aes(x = abundance_pred, y = cover)) +
      geom_point(alpha = 0.6, color = "darkblue", size = 1) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      xlim(axis_min, axis_max) +
      ylim(axis_min, axis_max) +
      labs(title = paste("Method B: Both Actual & Predicted Present\n(R² =", round(result$method_B_r2, 3), ", n =", nrow(method_B_data), ")"),
           x = "Abundance Prediction", y = "Actual Cover") +
      theme_minimal() +
      theme(plot.title = element_text(size = 8))
    
    # Plot 8: Method C scatter plot
    p8 <- ggplot(method_C_data, aes(x = abundance_pred, y = cover)) +
      geom_point(alpha = 0.6, color = "darkorange", size = 1) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      xlim(axis_min, axis_max) +
      ylim(axis_min, axis_max) +
      labs(title = paste("Method C: Predicted Present\n(R² =", round(result$method_C_r2, 3), ", n =", nrow(method_C_data), ")"),
           x = "Abundance Prediction", y = "Actual Cover") +
      theme_minimal() +
      theme(plot.title = element_text(size = 8))
    
    # Plot 9: Overall scatter plot (all sites, using hurdle predictions)
    p9 <- ggplot(plot_data, aes(x = predicted, y = cover)) +
      geom_point(alpha = 0.6, color = "purple", size = 1) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      labs(title = paste("Overall Scatter\n(R² =", round(result$overall_r2, 3), ", n =", nrow(plot_data), ")"),
           x = "Hurdle Prediction", y = "Actual Cover") +
      theme_minimal() +
      theme(plot.title = element_text(size = 8))
    
    
    # # Arrange all 8 plots in 3x3 grid (with one empty space)
    # grid.arrange(p1, p2, p3, p4, p5, 
    #              p6, p7, p8, 
    #              ncol = 3, nrow = 3,
    #              top = paste("Hurdle Validation Results -", result$species, 
    #                          "(", result$abundance_distribution, "abundance model)"),
    #              layout_matrix = rbind(c(1,2,3), c(4,5,6), c(7,8,NA)))
    # Arrange all 9 plots in 3x3 grid (fills last empty spot)
    grid.arrange(p1, p2, p3, 
                 p4, p5, p6, 
                 p7, p8, p9,
                 ncol = 3, nrow = 3,
                 top = paste("Hurdle Validation Results -", result$species, 
                             "(", result$abundance_distribution, " abundance model)"))
    
    # return(invisible(list(actual_cover = p1, predicted_cover = p2, 
    #                       observed_binary = p3, predicted_prob = p4, predicted_binary = p5,
    #                       method_A_scatter = p6, method_B_scatter = p7, method_C_scatter = p8)))
    return(invisible(list(actual_cover = p1, predicted_cover = p2, 
                          observed_binary = p3, predicted_prob = p4, predicted_binary = p5,
                          method_A_scatter = p6, method_B_scatter = p7, method_C_scatter = p8)))
    
  }
  
  # Test function to compare thresholds (updated to include beta parameter)
  test_thresholds <- function(species_name, thresholds = c(0.3, 0.4, 0.5, 0.6, 0.7), seed = 300, use_beta = FALSE) {
    results <- list()
    
    cat("Testing presence thresholds for", species_name, "with", ifelse(use_beta, "Beta", "Gamma"), "abundance model:\n")
    
    for(thresh in thresholds) {
      cat("\n--- Threshold:", thresh, "---\n")
      result <- validate_hurdle_simple(species_name, seed = seed, presence_threshold = thresh, 
                                       auto_plot = FALSE, use_beta = use_beta)
      results[[as.character(thresh)]] <- result
    }
    
    # Summary table
    cat("\n=== THRESHOLD COMPARISON SUMMARY ===\n")
    summary_df <- data.frame(
      Threshold = thresholds,
      Predicted_Present = sapply(results, function(x) x$n_predicted_present),
      Both_Present = sapply(results, function(x) x$n_both_present),
      Absence_Accuracy = sapply(results, function(x) round(x$absence_accuracy, 3)),
      Presence_Accuracy = sapply(results, function(x) round(x$presence_accuracy, 3)),
      Overall_R2 = sapply(results, function(x) round(x$overall_r2, 3)),
      Method_A_R2 = sapply(results, function(x) round(x$method_A_r2, 3)),
      Method_B_R2 = sapply(results, function(x) round(x$method_B_r2, 3)),
      Method_C_R2 = sapply(results, function(x) round(x$method_C_r2, 3))
    )
    
    print(summary_df)
    
    return(results)
  }
  
  # Compare Beta vs Gamma function
  compare_beta_gamma <- function(species_name, train_percent = 80, seed = 300, presence_threshold = 0.5) {
    
    cat("=== COMPARING BETA vs GAMMA ABUNDANCE MODELS ===\n")
    cat("Species:", species_name, "\n")
    cat("Seed:", seed, "\n\n")
    
    # Run gamma model
    cat("Running GAMMA abundance model...\n")
    gamma_result <- validate_hurdle_simple(species_name, train_percent = train_percent, 
                                           seed = seed, presence_threshold = presence_threshold,
                                           auto_plot = FALSE, use_beta = FALSE)
    
    # Run beta model  
    cat("\nRunning BETA abundance model...\n")
    beta_result <- validate_hurdle_simple(species_name, train_percent = train_percent,
                                          seed = seed, presence_threshold = presence_threshold, 
                                          auto_plot = FALSE, use_beta = TRUE)
    
    # Comparison summary
    cat("\n=== COMPARISON SUMMARY ===\n")
    comparison_df <- data.frame(
      Model = c("Gamma", "Beta"),
      Overall_R2 = c(round(gamma_result$overall_r2, 3), round(beta_result$overall_r2, 3)),
      Method_A_R2 = c(round(gamma_result$method_A_r2, 3), round(beta_result$method_A_r2, 3)),
      Method_B_R2 = c(round(gamma_result$method_B_r2, 3), round(beta_result$method_B_r2, 3)),
      Method_C_R2 = c(round(gamma_result$method_C_r2, 3), round(beta_result$method_C_r2, 3)),
      Overall_MAE = c(round(gamma_result$overall_mae, 3), round(beta_result$overall_mae, 3)),
      Method_A_MAE = c(round(gamma_result$method_A_mae, 3), round(beta_result$method_A_mae, 3)),
      Presence_AUC = c(round(gamma_result$presence_auc, 3), round(beta_result$presence_auc, 3))
    )
    
    print(comparison_df)
    
    # Highlight better performing model for each metric
    cat("\n=== PERFORMANCE WINNERS ===\n")
    cat("Overall R²:", ifelse(gamma_result$overall_r2 > beta_result$overall_r2, "Gamma", "Beta"), 
        "(", round(max(gamma_result$overall_r2, beta_result$overall_r2), 3), ")\n")
    cat("Method A R²:", ifelse(gamma_result$method_A_r2 > beta_result$method_A_r2, "Gamma", "Beta"),
        "(", round(max(gamma_result$method_A_r2, beta_result$method_A_r2, na.rm = TRUE), 3), ")\n")
    cat("Overall MAE:", ifelse(gamma_result$overall_mae < beta_result$overall_mae, "Gamma", "Beta"),
        "(", round(min(gamma_result$overall_mae, beta_result$overall_mae), 3), ")\n")
    
    
    
    return(list(gamma = gamma_result, beta = beta_result, comparison = comparison_df))
  }
  
  # USAGE EXAMPLES:
  
  # Basic usage with default gamma
  # result_gamma <- validate_hurdle_simple("orbicella")
  
  # Test beta model
  # result_beta <- validate_hurdle_simple("orbicella", use_beta = TRUE)
  
  # Compare both models directly
  # comparison <- compare_beta_gamma("orbicella")
  
  # Custom parameters with beta
  result <- validate_hurdle_simple("orbicella", train_percent = 80, presence_threshold = 0.5, use_beta = TRUE) # 0.38
  
  # Test thresholds with beta model
  # threshold_results_beta <- test_thresholds("orbicella", use_beta = TRUE)
  
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
  
  # #load derived bathy rasters, and oceanographic rasters
  # load_spat_objects(directory = 'output/output_calculate_bathy_rasters/')
  # load_spat_objects(directory = 'output/output_calculate_ocean_rasters/')
  
  #pull habitat grid
  source(here("src/functions.R"))
  load(here('output', 'output_create_habitat_grid/create_habitat_grid_workspace.RData'))
  
  #load GAM output
  load(here('output', 'output_GAMs/output_GAMs_workspace.RData'))
  source(here("src/functions.R"))
  load_spat_objects(directory = 'output/output_GAMs/')
  
  # FIGURING THIS OUT
  
  # #load GAMs
  # load_spat_objects(directory = 'output/output_calculate_bathy_rasters/')
  # load_spat_objects(directory = 'output/output_calculate_ocean_rasters/')
  
  
  # # Load fitted GAM models
  # load_all_files <- TRUE #Toggle: set to TRUE to load all .rds files, FALSE to load only GAM models
  # #
  # if(load_all_files) {
  #   # Load all .rds files in the folder
  #   model_files <- list.files(here("output", "output_GAMs"), 
  #                             pattern = "\\.rds$", 
  #                             full.names = TRUE)
  # } else {
  #   # Load only GAM model files
  #   model_files <- list.files(here("output", "output_GAMs"), 
  #                             pattern = "_gam_(abundance_gamma|presence_binom)\\.rds$", 
  #                             full.names = TRUE)
  # }
  # 
  # for(file in model_files) {
  #   # Extract the base filename without extension
  #   model_name <- tools::file_path_sans_ext(basename(file))
  #   
  #   # Load the model and assign to global environment
  #   assign(model_name, readRDS(file), envir = .GlobalEnv)
  # }  
  
  
  
  
  # FIGURING THIS OUT
  
  
  
  
  
  #save information for exporting new objects later
  existing_objects <- ls(envir = .GlobalEnv)
  
  
  
  ################################## GAM Train-Test Split Function ##################################
  
  gam_train_test_split <- function(gam_model, gam_results_entry = NULL, train_percent = 80, 
                                   seed = 123, return_details = FALSE, use_hurdle = FALSE) {
    
    # Load required libraries
    library(mgcv)
    library(pROC)
    library(dplyr)
    library(ggplot2)
    library(viridis)
    library(gridExtra)
    
    set.seed(seed)
    
    if(use_hurdle && !is.null(gam_results_entry)) {
      # Use two-part hurdle model approach
      return(gam_hurdle_validation(gam_results_entry, train_percent, seed, return_details))
    }
    
    # Original single-model approach
    # Extract model information and use model's own data
    model_formula <- formula(gam_model)
    model_family <- family(gam_model)
    model_data <- gam_model$model
    
    # If full gam_results entry provided, use that data for coordinates
    if(!is.null(gam_results_entry) && !is.null(gam_results_entry$data)) {
      full_data <- gam_results_entry$data
      # Match the model data with full data to get coordinates
      # Assuming the model data is a subset of the full data
      model_data_with_coords <- merge(model_data, full_data, by = intersect(names(model_data), names(full_data)))
    } else {
      model_data_with_coords <- model_data
    }
    
    # Get response variable name
    response_var <- all.vars(model_formula)[1]
    
    # Remove rows with missing values for the variables in the formula
    formula_vars <- all.vars(model_formula)
    complete_data <- model_data_with_coords[complete.cases(model_data_with_coords[, formula_vars]), ]
    
    cat("Using", nrow(complete_data), "complete observations from model's data\n")
    cat("Model family:", model_family$family, "\n")
    cat("Train/test split:", train_percent, "/", 100 - train_percent, "\n")
    cat("Response variable:", response_var, "\n")
    
    if(nrow(complete_data) < 20) {
      warning("Insufficient data for validation (< 20 observations)")
      return(NULL)
    }
    
    # Determine if this is predicting binary outcomes (0/1 only) vs continuous
    unique_values <- unique(complete_data[[response_var]])
    unique_values <- unique_values[!is.na(unique_values)]
    is_binary <- length(unique_values) == 2 && all(unique_values %in% c(0, 1))
    
    # Calculate appropriate metrics
    calculate_metrics <- function(actual, predicted) {
      if(is_binary) {
        # For binary outcomes - calculate AUC and classification metrics
        if(length(unique(actual)) > 1) {
          auc_val <- as.numeric(auc(actual, predicted))
          pred_class <- ifelse(predicted > 0.5, 1, 0)
          accuracy <- mean(actual == pred_class)
          sensitivity <- if(sum(actual == 1) > 0) sum(actual == 1 & pred_class == 1) / sum(actual == 1) else NA
          specificity <- if(sum(actual == 0) > 0) sum(actual == 0 & pred_class == 0) / sum(actual == 0) else NA
          
          return(list(
            auc = auc_val,
            accuracy = accuracy,
            sensitivity = sensitivity,
            specificity = specificity
          ))
        } else {
          return(list(auc = NA, accuracy = NA, sensitivity = NA, specificity = NA))
        }
      } else {
        # For continuous outcomes - calculate RMSE, MAE, R², and Pearson correlation
        rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
        mae <- mean(abs(actual - predicted), na.rm = TRUE)
        
        # Calculate R²
        if(var(actual, na.rm = TRUE) > 0) {
          r_squared <- cor(actual, predicted, use = "complete.obs")^2
        } else {
          r_squared <- NA
        }
        
        # Calculate Pearson correlation coefficient and test significance
        cor_test <- cor.test(actual, predicted, method = "pearson")
        pearson_r <- cor_test$estimate
        pearson_p <- cor_test$p.value
        
        # Determine significance level
        if(is.na(pearson_p)) {
          significance <- "Cannot determine"
        } else if(pearson_p < 0.001) {
          significance <- "***"
        } else if(pearson_p < 0.01) {
          significance <- "**"
        } else if(pearson_p < 0.05) {
          significance <- "*"
        } else {
          significance <- "ns"
        }
        
        return(list(
          rmse = rmse,
          mae = mae,
          r_squared = r_squared,
          pearson_r = pearson_r,
          pearson_p = pearson_p,
          significance = significance
        ))
      }
    }
    
    # Train-test split
    n <- nrow(complete_data)
    train_size <- floor((train_percent / 100) * n)
    train_idx <- sample(1:n, size = train_size)
    
    train_data <- complete_data[train_idx, ]
    test_data <- complete_data[-train_idx, ]
    
    cat("Training observations:", nrow(train_data), "\n")
    cat("Testing observations:", nrow(test_data), "\n")
    
    # Fit model on training data
    tryCatch({
      fit <- gam(model_formula, data = train_data, family = model_family)
      
      # Make predictions on test data
      predictions <- predict(fit, newdata = test_data, type = "response")
      actual <- test_data[[response_var]]
      
      # Calculate metrics
      metrics <- calculate_metrics(actual, predictions)
      
      # Print key results
      cat("\n=== RESULTS ===\n")
      if(is_binary) {
        cat("AUC:", round(metrics$auc, 3), "\n")
        cat("Accuracy:", round(metrics$accuracy, 3), "\n")
      } else {
        cat("R²:", round(metrics$r_squared, 3), "\n")
        cat("Pearson r:", round(metrics$pearson_r, 3), "\n")
        cat("P-value:", format(metrics$pearson_p, scientific = TRUE, digits = 3), "\n")
        cat("Significance:", metrics$significance, "\n")
        cat("RMSE:", round(metrics$rmse, 3), "\n")
        
        # Provide interpretation
        if(metrics$significance == "***") {
          cat("Interpretation: Correlation is highly significant (P < 0.001)\n")
        } else if(metrics$significance == "**") {
          cat("Interpretation: Correlation is significant (P < 0.01)\n")
        } else if(metrics$significance == "*") {
          cat("Interpretation: Correlation is significant (P < 0.05)\n")
        } else if(metrics$significance == "ns") {
          cat("Interpretation: Correlation is not significant (P ≥ 0.05)\n")
        }
      }
      
      # Create results dataframe with spatial coordinates for mapping
      test_results <- test_data
      test_results$predicted <- predictions
      test_results$residuals <- actual - predictions
      
      result <- list(
        n_total = nrow(complete_data),
        n_train = nrow(train_data),
        n_test = nrow(test_data),
        train_percent = train_percent,
        family = model_family$family,
        response_variable = response_var,
        metrics = metrics,
        test_results = test_results,  # Always include this for mapping
        is_binary = is_binary
      )
      
      if(return_details) {
        result$predictions <- predictions
        result$actual <- actual
        result$train_data <- train_data
        result$test_data <- test_data
        result$fitted_model <- fit
        result$formula <- model_formula
      }
      
      return(result)
      
    }, error = function(e) {
      cat("Error in model fitting:", e$message, "\n")
      return(NULL)
    })
  }
  
  ################################## Two-Part Hurdle Model Validation ##################################
  
  gam_hurdle_validation <- function(gam_results_entry, train_percent = 80, seed = 123, return_details = FALSE) {
    
    set.seed(seed)
    
    # Check if both presence and abundance models exist
    if(is.null(gam_results_entry$presence) || is.null(gam_results_entry$abundance)) {
      stop("Both presence and abundance models required for hurdle validation")
    }
    
    # Get the full dataset
    full_data <- gam_results_entry$data
    
    # Get formulas from both models
    presence_formula <- formula(gam_results_entry$presence)
    abundance_formula <- formula(gam_results_entry$abundance)
    
    # Get all variables needed
    presence_vars <- all.vars(presence_formula)
    abundance_vars <- all.vars(abundance_formula)
    all_vars <- unique(c(presence_vars, abundance_vars))
    
    # Filter to complete cases for all variables
    complete_data <- full_data[complete.cases(full_data[, all_vars]), ]
    
    cat("=== HURDLE MODEL VALIDATION ===\n")
    cat("Using", nrow(complete_data), "complete observations\n")
    cat("Train/test split:", train_percent, "/", 100 - train_percent, "\n")
    
    if(nrow(complete_data) < 20) {
      warning("Insufficient data for validation (< 20 observations)")
      return(NULL)
    }
    
    # Train-test split
    n <- nrow(complete_data)
    train_size <- floor((train_percent / 100) * n)
    train_idx <- sample(1:n, size = train_size)
    
    train_data <- complete_data[train_idx, ]
    test_data <- complete_data[-train_idx, ]
    
    cat("Training observations:", nrow(train_data), "\n")
    cat("Testing observations:", nrow(test_data), "\n")
    
    tryCatch({
      # Fit presence model on training data
      presence_fit <- gam(presence_formula, data = train_data, family = binomial())
      
      # Fit abundance model on training data (only where cover > 0)
      abundance_train_data <- train_data[train_data$cover > 0, ]
      abundance_fit <- gam(abundance_formula, data = abundance_train_data, family = Gamma(link = "log"))
      
      # Make predictions on test data
      # Step 1: Predict presence probability
      presence_prob <- predict(presence_fit, newdata = test_data, type = "response")
      
      # Step 2: Predict abundance (conditional on presence)
      abundance_pred <- predict(abundance_fit, newdata = test_data, type = "response")
      
      # Step 3: Combine using hurdle approach with threshold
      # Method 1: Expected value = P(presence) * E(abundance | presence)
      hurdle_predictions_raw <- presence_prob * abundance_pred
      
      # Apply threshold for more realistic zero predictions
      presence_threshold <- 0.5  # Can be adjusted
      probability_threshold <- 0.1  # Very low probability predictions become zero
      
      # Method A: Use presence classification (0.5 threshold)
      presence_binary <- ifelse(presence_prob > presence_threshold, 1, 0)
      classification_predictions <- presence_binary * abundance_pred
      
      # Method B: Use probability threshold (very low probabilities become zero)
      hurdle_predictions_threshold <- ifelse(presence_prob < probability_threshold, 0, hurdle_predictions_raw)
      
      # Use the classification approach as default (more zeros)
      hurdle_predictions <- classification_predictions
      
      actual <- test_data$cover
      
      # Calculate metrics for both approaches
      calculate_hurdle_metrics <- function(actual, predicted) {
        rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
        mae <- mean(abs(actual - predicted), na.rm = TRUE)
        
        if(var(actual, na.rm = TRUE) > 0) {
          r_squared <- cor(actual, predicted, use = "complete.obs")^2
        } else {
          r_squared <- NA
        }
        
        cor_test <- cor.test(actual, predicted, method = "pearson")
        pearson_r <- cor_test$estimate
        pearson_p <- cor_test$p.value
        
        if(is.na(pearson_p)) {
          significance <- "Cannot determine"
        } else if(pearson_p < 0.001) {
          significance <- "***"
        } else if(pearson_p < 0.01) {
          significance <- "**"
        } else if(pearson_p < 0.05) {
          significance <- "*"
        } else {
          significance <- "ns"
        }
        
        # Calculate how many zeros are correctly predicted
        actual_zeros <- sum(actual == 0)
        predicted_zeros <- sum(predicted == 0)
        correct_zeros <- sum(actual == 0 & predicted == 0)
        
        return(list(
          rmse = rmse,
          mae = mae,
          r_squared = r_squared,
          pearson_r = pearson_r,
          pearson_p = pearson_p,
          significance = significance,
          actual_zeros = actual_zeros,
          predicted_zeros = predicted_zeros,
          correct_zeros = correct_zeros
        ))
      }
      
      # Use the expected value approach as primary
      metrics <- calculate_hurdle_metrics(actual, hurdle_predictions)
      
      # Print results
      cat("\n=== HURDLE MODEL RESULTS ===\n")
      cat("R²:", round(metrics$r_squared, 3), "\n")
      cat("Pearson r:", round(metrics$pearson_r, 3), "\n")
      cat("P-value:", format(metrics$pearson_p, scientific = TRUE, digits = 3), "\n")
      cat("Significance:", metrics$significance, "\n")
      cat("RMSE:", round(metrics$rmse, 3), "\n")
      cat("Actual zeros:", metrics$actual_zeros, "/ Predicted zeros:", metrics$predicted_zeros, "\n")
      cat("Correctly predicted zeros:", metrics$correct_zeros, "\n")
      cat("Presence threshold used:", presence_threshold, "\n")
      
      # Provide interpretation
      if(metrics$significance == "***") {
        cat("Interpretation: Correlation is highly significant (P < 0.001)\n")
      } else if(metrics$significance == "**") {
        cat("Interpretation: Correlation is significant (P < 0.01)\n")
      } else if(metrics$significance == "*") {
        cat("Interpretation: Correlation is significant (P < 0.05)\n")
      } else if(metrics$significance == "ns") {
        cat("Interpretation: Correlation is not significant (P ≥ 0.05)\n")
      }
      
      # Create results dataframe
      test_results <- test_data
      test_results$predicted <- hurdle_predictions
      test_results$presence_prob <- presence_prob
      test_results$abundance_pred <- abundance_pred
      test_results$residuals <- actual - hurdle_predictions
      
      result <- list(
        n_total = nrow(complete_data),
        n_train = nrow(train_data),
        n_test = nrow(test_data),
        train_percent = train_percent,
        family = "Hurdle (Binomial + Gamma)",
        response_variable = "cover",
        metrics = metrics,
        test_results = test_results,
        is_binary = FALSE,
        model_type = "hurdle"
      )
      
      if(return_details) {
        result$presence_fit <- presence_fit
        result$abundance_fit <- abundance_fit
        result$presence_predictions <- presence_prob
        result$abundance_predictions <- abundance_pred
        result$hurdle_predictions <- hurdle_predictions
        result$hurdle_predictions_raw <- hurdle_predictions_raw
        result$hurdle_predictions_threshold <- hurdle_predictions_threshold
        result$presence_threshold <- presence_threshold
        result$probability_threshold <- probability_threshold
      }
      
      return(result)
      
    }, error = function(e) {
      cat("Error in hurdle model fitting:", e$message, "\n")
      return(NULL)
    })
  }
  
  ################################## Mapping Function ##################################
  
  plot_validation_maps <- function(validation_result, coord_system = "latlon") {
    
    if(is.null(validation_result) || is.null(validation_result$test_results)) {
      stop("No test results found for mapping")
    }
    
    plot_data <- validation_result$test_results
    
    # Print available columns for debugging
    cat("Available columns in test results:\n")
    print(colnames(plot_data))
    
    # Determine coordinate columns with fallbacks
    if(coord_system == "latlon") {
      x_candidates <- c("lon", "longitude", "Longitude", "LON")
      y_candidates <- c("lat", "latitude", "Latitude", "LAT")
      x_lab <- "Longitude"
      y_lab <- "Latitude"
    } else {
      x_candidates <- c("x_utm", "X_utm", "utm_x", "UTM_X")
      y_candidates <- c("y_utm", "Y_utm", "utm_y", "UTM_Y")
      x_lab <- "UTM X"
      y_lab <- "UTM Y"
    }
    
    # Find the actual column names
    x_col <- x_candidates[x_candidates %in% colnames(plot_data)]
    y_col <- y_candidates[y_candidates %in% colnames(plot_data)]
    
    if(length(x_col) == 0 || length(y_col) == 0) {
      cat("Could not find coordinate columns. Trying the other coordinate system...\n")
      
      # Try the other coordinate system
      if(coord_system == "latlon") {
        x_candidates <- c("x_utm", "X_utm", "utm_x", "UTM_X")
        y_candidates <- c("y_utm", "Y_utm", "utm_y", "UTM_Y")
        x_lab <- "UTM X"
        y_lab <- "UTM Y"
      } else {
        x_candidates <- c("lon", "longitude", "Longitude", "LON")
        y_candidates <- c("lat", "latitude", "Latitude", "LAT")
        x_lab <- "Longitude"
        y_lab <- "Latitude"
      }
      
      x_col <- x_candidates[x_candidates %in% colnames(plot_data)]
      y_col <- y_candidates[y_candidates %in% colnames(plot_data)]
    }
    
    if(length(x_col) == 0 || length(y_col) == 0) {
      stop("Could not find any coordinate columns. Available columns: ", paste(colnames(plot_data), collapse = ", "))
    }
    
    # Use first match
    x_col <- x_col[1]
    y_col <- y_col[1]
    
    cat("Using coordinates:", x_col, "and", y_col, "\n")
    
    # Common theme for maps
    map_theme <- theme_minimal() +
      theme(
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)
      )
    
    # Get response variable name
    response_var <- validation_result$response_variable
    
    # Calculate shared color scale limits
    actual_values <- plot_data[[response_var]]
    predicted_values <- plot_data$predicted
    shared_min <- min(c(actual_values, predicted_values), na.rm = TRUE)
    shared_max <- max(c(actual_values, predicted_values), na.rm = TRUE)
    
    cat("Shared color scale range:", round(shared_min, 2), "to", round(shared_max, 2), "\n")
    
    # Plot 1: Actual cover values
    p1 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = response_var)) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_viridis_c(name = "Actual", option = "plasma", limits = c(shared_min, shared_max)) +
      labs(title = paste("Actual", response_var, "(Test Set)"), 
           x = x_lab, y = y_lab) +
      map_theme
    
    # Plot 2: Predicted cover values
    p2 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = "predicted")) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_viridis_c(name = "Predicted", option = "plasma", limits = c(shared_min, shared_max)) +
      labs(title = paste("Predicted", response_var), 
           x = x_lab, y = y_lab) +
      map_theme
    
    # Plot 3: Residuals (actual - predicted)
    p3 <- ggplot(plot_data, aes_string(x = x_col, y = y_col, color = "residuals")) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_gradient2(name = "Residual", 
                            low = "red", mid = "white", high = "blue",
                            midpoint = 0) +
      labs(title = "Residuals (Actual - Predicted)", 
           x = x_lab, y = y_lab) +
      map_theme
    
    # Plot 4: Predicted vs Actual scatter plot
    p4 <- ggplot(plot_data, aes_string(x = "predicted", y = response_var)) +
      geom_point(alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      labs(title = "Predicted vs Actual", 
           x = paste("Predicted", response_var), 
           y = paste("Actual", response_var)) +
      theme_minimal()
    
    # Arrange plots
    grid.arrange(p1, p2, p3, p4, ncol = 2,
                 top = paste("Validation Results -", 
                             validation_result$family, "model"))
    
    # Return the individual plots as well
    invisible(list(actual = p1, predicted = p2, residuals = p3, scatter = p4))
  }
  
  ################################## Usage Examples ##################################
  
  sppofinterest = "Orbicella"
  
  # Example 1: Standard single-model validation
  result <- gam_train_test_split(gam_results[[sppofinterest]]$complex,
                                 gam_results_entry = gam_results[[sppofinterest]],
                                 train_percent = 80)
  plot_validation_maps(result)
  
  
  # Example 2: Hurdle model validation (uses presence + abundance models)
  result_hurdle <- gam_train_test_split(gam_model = NULL,  # Not used for hurdle
                                        gam_results_entry = gam_results[[sppofinterest]],
                                        train_percent = 80,
                                        use_hurdle = TRUE, seed = sample(1:10000, 1))
  plot_validation_maps(result_hurdle)
  
  # Example 3: Compare single vs hurdle approaches
  # single_result <- gam_train_test_split(gam_results[["Agaricia"]]$complex, 
  #                                      gam_results_entry = gam_results[["Agaricia"]])
  # hurdle_result <- gam_train_test_split(gam_model = NULL, 
  #                                      gam_results_entry = gam_results[["Agaricia"]], 
  #                                      use_hurdle = TRUE)
  # 
  # cat("Single model R²:", round(single_result$metrics$r_squared, 3), "\n")
  # cat("Hurdle model R²:", round(hurdle_result$metrics$r_squared, 3), "\n")
  # cat("Single model predicted zeros:", sum(single_result$test_results$predicted == 0), "\n")
  # cat("Hurdle model predicted zeros:", hurdle_result$metrics$predicted_zeros, "\n")
  
  # Example 4: Generate summary table for multiple species
  # results_summary <- data.frame(
  #   Species = character(),
  #   Model_Type = character(),
  #   Pearson_r = numeric(),
  #   P_value = numeric(),
  #   Significance = character(),
  #   R_squared = numeric(),
  #   RMSE = numeric(),
  #   stringsAsFactors = FALSE
  # )
  # 
  # for(sp in names(gam_results)) {
  #   for(model_type in c("complex", "reduced")) {
  #     result <- gam_train_test_split(gam_results[[sp]][[model_type]], 
  #                                   gam_results_entry = gam_results[[sp]], 
  #                                   train_percent = 80)
  #     if(!is.null(result) && !result$is_binary) {
  #       results_summary <- rbind(results_summary, data.frame(
  #         Species = sp,
  #         Model_Type = model_type,
  #         Pearson_r = round(result$metrics$pearson_r, 3),
  #         P_value = result$metrics$pearson_p,
  #         Significance = result$metrics$significance,
  #         R_squared = round(result$metrics$r_squared, 3),
  #         RMSE = round(result$metrics$rmse, 3)
  #       ))
  #     }
  #   }
  # }
  # print(results_summary)
  
  # ################################## make predictions ##################################
  # 
  # # NOTE - this is sparse right now compared to all available rasters, for computational
  # #         reasons. can return to this
  # # First, create your stack properly
  # env_stack <- env_complex
  # 
  # # Check geometry - compareGeom works on individual rasters, not the stack
  # # Let's check that all layers match the first one (bathy_final)
  # terra::compareGeom(bathy_final, aspect_terra, stopOnError = FALSE)
  # terra::compareGeom(bathy_final, slope_terra, stopOnError = FALSE)
  # # etc... or loop through them:
  # 
  # for(i in 2:nlyr(env_stack)) {
  #   cat("Checking layer", i, ":", names(env_stack)[i], "\n")
  #   print(terra::compareGeom(env_stack[[1]], env_stack[[i]], stopOnError = FALSE))
  # }
  # 
  # # Alternative: just check basic properties
  # terra::res(env_stack)  # Should all be the same
  # terra::ext(env_stack)  # Should all be the same
  # 
  # # Create prediction grid
  # prediction_grid <- as.data.frame(env_stack, xy = TRUE)
  # 
  # # # Remove rows with any NA values
  # # prediction_grid <- prediction_grid[complete.cases(prediction_grid), ]
  # 
  # # Check your grid
  # dim(prediction_grid)
  # head(prediction_grid)
  # 
  # # Better names for your variables
  # names(prediction_grid) <- c("x", "y", "depth", "aspect", "slope", "complexity",
  #                             "TPI", "VRM", "planform_curv", "max_Hsig", "mean_dir",
  #                             "mean_Hsig", "mean_SST", "range_SST")
  # 
  # # Convert prediction grid UTM coordinates back to lat/lon for the spatial term
  # prediction_coords_utm <- vect(cbind(prediction_grid$x, prediction_grid$y),
  #                               crs = crs(bathy_final))
  # prediction_coords_latlon <- project(prediction_coords_utm, "EPSG:4326")
  # coords_df <- as.data.frame(geom(prediction_coords_latlon)[, c("x", "y")])
  # 
  # # Add lat/lon to prediction grid
  # prediction_grid$lon <- coords_df$x
  # prediction_grid$lat <- coords_df$y
  # 
  # # Remove any rows with NAs in the required variables
  # # pred_vars <- c("bathymetry", "TPI", "VRM", "lon", "lat")
  # pred_vars <- c("depth", "TPI", "slope", "complexity", "planform_curv",
  #                "range_SST", "mean_SST", "mean_dir", "max_Hsig", "mean_Hsig", "lon", "lat")
  # prediction_grid_clean <- prediction_grid[complete.cases(prediction_grid[, pred_vars]), ]
  # 
  # cat("Prediction grid size:", nrow(prediction_grid_clean), "cells\n")
  # 
  # 
  # # # Make predictions with the spatial model
  # # predictions <- predict(gam_spatial, prediction_grid_clean, type = "response")
  # #
  # # # Add predictions to the grid
  # # prediction_grid_clean$predicted_cover <- predictions
  # #
  # # # Check prediction range
  # # summary(prediction_grid_clean$predicted_cover)
  # 
  # 
  # # # Simple approach with time estimates
  # # start_time <- Sys.time()
  # #
  # # # Test prediction on small subset to estimate time
  # # test_subset <- prediction_grid_clean[1:1000, ]
  # # test_start <- Sys.time()
  # # test_pred <- predict(gam_complex, test_subset, type = "response")
  # # test_time <- as.numeric(difftime(Sys.time(), test_start, units = "secs"))
  # #
  # # # Estimate total time
  # # estimated_total <- (test_time / 1000) * nrow(prediction_grid_clean)
  # # cat("Estimated total time:", round(estimated_total/60, 1), "minutes\n")
  # #
  # # # Now do full prediction
  # # predictions <- predict(gam_spatial, prediction_grid_clean, type = "response")
  # # cat("Actual time:", round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes\n")
  # 
  # 
  # # Resample your environmental stack to 200m resolution
  # env_simple_200m <- aggregate(env_stack, fact = 4, fun = "mean")  # 4x aggregation = 200m
  # 
  # # Check the new resolution
  # res(env_simple_200m)  # Should show 200, 200
  # 
  # # Create prediction grid from 200m rasters
  # prediction_grid_200m <- as.data.frame(env_simple_200m, xy = TRUE)
  # names(prediction_grid_200m) <- c("x", "y", "bathymetry", "slope")
  # 
  # # Add TPI and VRM at 200m
  # env_complex_200m <- aggregate(env_complex, fact = 4, fun = "mean")
  # complex_grid_200m <- as.data.frame(env_complex_200m, xy = TRUE)
  # names(complex_grid_200m) <- c("x", "y", "bathymetry", "slope", "TPI", "roughness", "VRM",
  #                               "max_curv", "mean_curv", "planform_curv", "profile_curv")
  # 
  # # Keep just the variables you need
  # prediction_grid_200m <- complex_grid_200m[, c("x", "y", "bathymetry", "TPI", "VRM")]
  # 
  # # Remove NAs
  # prediction_grid_200m <- prediction_grid_200m[complete.cases(prediction_grid_200m), ]
  # 
  # cat("200m grid size:", nrow(prediction_grid_200m), "cells\n")
  # cat("50m grid size was:", nrow(prediction_grid_clean), "cells\n")
  # cat("Reduction factor:", round(nrow(prediction_grid_clean)/nrow(prediction_grid_200m), 1), "\n")
  # 
  # # Convert UTM to lat/lon for the spatial term
  # prediction_coords_utm_200m <- vect(cbind(prediction_grid_200m$x, prediction_grid_200m$y),
  #                                    crs = crs(bathy_final))
  # prediction_coords_latlon_200m <- project(prediction_coords_utm_200m, "EPSG:4326")
  # coords_df_200m <- as.data.frame(geom(prediction_coords_latlon_200m)[, c("x", "y")])
  # 
  # prediction_grid_200m$lon <- coords_df_200m$x
  # prediction_grid_200m$lat <- coords_df_200m$y
  # 
  # cat("Final 200m grid size:", nrow(prediction_grid_200m), "cells\n")
  # 
  # # Step 4: Estimate time before full prediction
  # start_time <- Sys.time()
  # 
  # # Test prediction on small subset to estimate time
  # test_subset_200m <- prediction_grid_200m[1:1000, ]
  # test_start <- Sys.time()
  # test_pred_200m <- predict(gam_spatial, test_subset_200m, type = "response")
  # test_time <- as.numeric(difftime(Sys.time(), test_start, units = "secs"))
  # 
  # # Estimate total time
  # estimated_total_200m <- (test_time / 1000) * nrow(prediction_grid_200m)
  # cat("200m grid size:", nrow(prediction_grid_200m), "cells\n")
  # cat("Test time for 1000 cells:", round(test_time, 2), "seconds\n")
  # cat("Estimated total time:", round(estimated_total_200m/60, 1), "minutes\n")
  # 
  # # Compare to original estimate
  # original_cells <- nrow(prediction_grid_clean)
  # reduction_factor <- original_cells / nrow(prediction_grid_200m)
  # cat("Reduction factor:", round(reduction_factor, 1), "x smaller\n")
  # cat("Original estimate was 61 minutes\n")
  # cat("New estimate should be ~", round(61/reduction_factor, 1), "minutes\n")
  # 
  # # Ask user if they want to proceed
  # cat("\nProceed with full prediction? (y/n)\n")
  # user_input <- readline()
  # 
  # if(tolower(user_input) == "y") {
  #   cat("Starting full prediction...\n")
  # 
  #   # Full prediction with timing
  #   full_start <- Sys.time()
  #   predictions_200m <- predict(gam_spatial, prediction_grid_200m, type = "response")
  #   prediction_grid_200m$predicted_cover <- predictions_200m
  #   full_end <- Sys.time()
  # 
  #   cat("Actual prediction time:", round(difftime(full_end, full_start, units = "mins"), 1), "minutes\n")
  #   cat("Prediction complete!\n")
  # 
  # } else {
  #   cat("Prediction cancelled. Consider further reducing grid size or using sampling.\n")
  # }
  # 
  # # Convert to raster
  # pred_raster_200m <- rast(prediction_grid_200m[, c("x", "y", "predicted_cover")],
  #                          crs = crs(bathy_final))
  # 
  # # Plot
  # plot(pred_raster_200m,
  #      main = "Predicted Coral Cover (%) - 200m resolution",
  #      col = viridis(100))
  # 
  # # Add survey points
  # points(model_data_complex$x_utm, model_data_complex$y_utm,
  #        pch = 21,
  #        bg = heat.colors(10)[cut(model_data_complex$cover, breaks = 10)],
  #        cex = 0.8)
  # 
  # 
  # ################################## maps ##################################
  # 
  # 
  # # library(viridis)
  # library(sf)
  # 
  # # Convert predictions back to raster
  # pred_raster <- rast(prediction_grid_clean[, c("x", "y", "predicted_cover")], 
  #                     crs = crs(bathy_final))
  # 
  # # Create a nice map
  # # Option 1: Using terra plot
  # plot(pred_raster, 
  #      main = "Predicted Coral Cover (%)",
  #      col = viridis(100),
  #      range = c(0, max(prediction_grid_clean$predicted_cover, na.rm = TRUE)))
  # 
  # # Add your actual survey points
  # points(model_data_complex$x_utm, model_data_complex$y_utm, 
  #        pch = 21, 
  #        bg = heat.colors(10)[cut(model_data_complex$cover, breaks = 10)],
  #        cex = 0.8)
  # 
  # # Add legend for points
  # legend("topright", 
  #        title = "Observed Cover",
  #        legend = c("0-10%", "10-20%", "20-30%", "30%+"),
  #        pch = 21,
  #        pt.bg = heat.colors(4),
  #        cex = 0.8)
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  
  
  
  ################################## Save objects/workspace ##################################
  
  # #updated way to handle saving of new objects
  # save_new_objects("output/GAMs", existing_objects)
