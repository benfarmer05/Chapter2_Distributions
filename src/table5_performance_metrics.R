# .rs.restartR(clean = TRUE)
rm(list=ls())

library(here)
library(mgcv)
library(dplyr)
library(ggplot2)
library(patchwork)

################################## Setup ##################################

cat("========================================\n")
cat("Abundance Model Cross-Validation Summary\n")
cat("========================================\n\n")

# Load combined data
load(here("output", "all_combined_data.rda"))

cat("Reading model outputs...\n")

# Get all results files
results_files <- list.files(here('output/output_maps'), 
                            pattern = '^results_light_.*\\.rds$', 
                            full.names = TRUE)

# Read all results
all_results <- lapply(results_files, function(f) {
  result <- readRDS(f)
  return(result)
})
names(all_results) <- gsub('results_light_|\\.rds', '', basename(results_files))

cat("Loaded", length(all_results), "species models\n\n")

existing_objects <- ls(envir = .GlobalEnv)

################################## Helper Functions ##################################

# Calculate regression metrics
calculate_regression_metrics <- function(observed, predicted) {
  valid_idx <- complete.cases(observed, predicted)
  observed <- observed[valid_idx]
  predicted <- predicted[valid_idx]
  
  r <- cor(observed, predicted, method = "pearson")
  r_squared <- r^2
  
  ss_res <- sum((observed - predicted)^2)
  ss_tot <- sum((observed - mean(observed))^2)
  R2 <- 1 - (ss_res / ss_tot)
  
  rmse <- sqrt(mean((observed - predicted)^2))
  mae <- mean(abs(observed - predicted))
  
  return(c(r = r, r2 = r_squared, R2 = R2, RMSE = rmse, MAE = mae))
}

################################## Cross-Validation ##################################

cat("Running 5-fold cross-validation on abundance models...\n\n")

# Calibration function
calibrate_quantiles <- function(predictions, observations, n_quantiles = 100) {
  probs <- seq(0, 1, length.out = n_quantiles)
  pred_q <- quantile(predictions, probs)
  obs_q <- quantile(observations, probs)
  correction_fn <- approxfun(pred_q, obs_q, method = "linear", rule = 2)
  list(correction_function = correction_fn)
}

# Initialize results storage
cv_results_raw_list <- list()
cv_results_calib_list <- list()
overall_results_raw_list <- list()
overall_results_calib_list <- list()
scatter_data_list <- list()  # Store data for scatter plots

for(species in names(all_results)) {
  cat("Processing", species, "...\n")
  
  model_file <- here("output", "output_GAMs", paste0(species, "_models.rds"))
  if(!file.exists(model_file)) {
    cat("  Warning: Model file not found, skipping\n")
    next
  }
  
  species_models <- readRDS(model_file)
  abundance_model <- species_models$abundance_model
  model_data <- species_models$model_data
  
  all_vars <- all.vars(formula(abundance_model))
  
  if("cover_prop" %in% names(model_data)) {
    presence_data <- model_data[model_data$cover_prop > 0, ]
  } else {
    model_data$cover_prop <- model_data$cover / 100
    presence_data <- model_data[model_data$cover_prop > 0, ]
  }
  
  complete_data <- presence_data[complete.cases(presence_data[, all_vars]), ]
  
  cat("  Data: N =", nrow(complete_data), 
      "| Mean cover:", round(mean(complete_data$cover_prop) * 100, 2), "%\n")
  
  if(nrow(complete_data) < 20) {
    cat("  Warning: Too few observations (< 20), skipping\n")
    next
  }
  
  # ========== OVERALL METRICS (full dataset) ==========
  overall_pred_raw <- predict(abundance_model, newdata = complete_data, type = "response")
  
  # Raw metrics
  overall_metrics_raw <- calculate_regression_metrics(complete_data$cover_prop, overall_pred_raw)
  overall_results_raw_list[[species]] <- overall_metrics_raw
  
  # Calibrated metrics
  cal_full <- calibrate_quantiles(overall_pred_raw, complete_data$cover_prop, n_quantiles = 100)
  overall_pred_calibrated <- cal_full$correction_function(overall_pred_raw)
  overall_metrics_calib <- calculate_regression_metrics(complete_data$cover_prop, overall_pred_calibrated)
  overall_results_calib_list[[species]] <- overall_metrics_calib
  
  # ========== CROSS-VALIDATION WITH NESTED CALIBRATION ==========
  n <- nrow(complete_data)
  fold_ids <- sample(rep(1:5, length.out = n))
  
  cv_metrics_raw <- matrix(NA, nrow = 5, ncol = 5)
  cv_metrics_calib <- matrix(NA, nrow = 5, ncol = 5)
  colnames(cv_metrics_raw) <- colnames(cv_metrics_calib) <- c("r", "r2", "R2", "RMSE", "MAE")
  
  # Store one random fold for plotting
  chosen_fold <- sample(1:5, 1)
  
  for(fold in 1:5) {
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)
    
    if(length(test_idx) < 5) {
      cat("  Warning: Fold", fold, "test set too small, skipping fold\n")
      next
    }
    
    abundance_fit <- tryCatch({
      gam(formula(abundance_model), 
          data = complete_data[train_idx, ], 
          family = betar())
    }, error = function(e) {
      cat("  Warning: Fold", fold, "model fitting failed\n")
      return(NULL)
    })
    
    if(is.null(abundance_fit)) next
    
    test_pred_raw <- predict(abundance_fit, 
                             newdata = complete_data[test_idx, ], 
                             type = "response")
    
    train_pred_fold <- predict(abundance_fit,
                               newdata = complete_data[train_idx, ],
                               type = "response")
    train_obs_fold <- complete_data$cover_prop[train_idx]
    
    # Calibration
    cal_cv <- calibrate_quantiles(train_pred_fold, train_obs_fold, n_quantiles = 100)
    test_pred_calibrated <- cal_cv$correction_function(test_pred_raw)
    test_obs <- complete_data$cover_prop[test_idx]
    
    # Calculate metrics for raw
    cv_metrics_raw[fold, ] <- calculate_regression_metrics(test_obs, test_pred_raw)
    
    # Calculate metrics for calibrated
    cv_metrics_calib[fold, ] <- calculate_regression_metrics(test_obs, test_pred_calibrated)
    
    # Store scatter data for chosen fold
    if(fold == chosen_fold) {
      scatter_data_list[[species]] <- list(
        overall_raw = data.frame(observed = complete_data$cover_prop, 
                                 predicted = overall_pred_raw, 
                                 type = "Overall_Raw"),
        overall_calib = data.frame(observed = complete_data$cover_prop, 
                                   predicted = overall_pred_calibrated, 
                                   type = "Overall_Calib"),
        test_raw = data.frame(observed = test_obs, 
                              predicted = test_pred_raw, 
                              type = "Test_Raw"),
        test_calib = data.frame(observed = test_obs, 
                                predicted = test_pred_calibrated, 
                                type = "Test_Calib")
      )
    }
  }
  
  # Average across folds
  cv_results_raw_list[[species]] <- colMeans(cv_metrics_raw, na.rm = TRUE)
  cv_results_calib_list[[species]] <- colMeans(cv_metrics_calib, na.rm = TRUE)
  
  cat("  CV complete\n")
}

cat("\n✓ Cross-validation complete for all species\n\n")

################################## Format Results ##################################

# Create summary tables
cv_summary_raw <- do.call(rbind, cv_results_raw_list)
cv_summary_raw_df <- as.data.frame(cv_summary_raw)
cv_summary_raw_df$Species <- rownames(cv_summary_raw_df)

cv_summary_calib <- do.call(rbind, cv_results_calib_list)
cv_summary_calib_df <- as.data.frame(cv_summary_calib)
cv_summary_calib_df$Species <- rownames(cv_summary_calib_df)

overall_summary_raw <- do.call(rbind, overall_results_raw_list)
overall_summary_raw_df <- as.data.frame(overall_summary_raw)
overall_summary_raw_df$Species <- rownames(overall_summary_raw_df)

overall_summary_calib <- do.call(rbind, overall_results_calib_list)
overall_summary_calib_df <- as.data.frame(overall_summary_calib)
overall_summary_calib_df$Species <- rownames(overall_summary_calib_df)

# Species mapping
species_codes <- c("agaricia", "madracis", "porites", "siderastrea", "stephanocoenia",
                   "montastraea", "orbicella",
                   "colpophyllia", "dendrogyra", "dichocoenia", "diploria",
                   "eusmilia", "meandrina", "mycetophyllia", "pseudodiploria")

display_names <- c("Agaricia", "Madracis", "Porites", "Siderastrea", "S. intersepta",
                   "M. cavernosa", "Orbicella",
                   "C. natans", "D. cylindrus", "D. stokesii", "D. labyrinthiformis",
                   "E. fastigiata", "Meandrina", "Mycetophyllia", "Pseudodiploria")

sus_groups <- c("LS", "LS", "LS", "LS", "LS", "MS", "MS", 
                "HS", "HS", "HS", "HS", "HS", "HS", "HS", "HS")

species_mapping <- data.frame(code = species_codes, display = display_names, 
                              sus = sus_groups, stringsAsFactors = FALSE)

# Format tables
format_table <- function(cv_df, overall_df) {
  cv_df <- cv_df %>%
    mutate(Species_lower = tolower(Species)) %>%
    left_join(species_mapping, by = c("Species_lower" = "code")) %>%
    select(Taxon = display, Sus = sus, r, r2, R2, RMSE, MAE) %>%
    arrange(match(Taxon, display_names))
  
  overall_df <- overall_df %>%
    mutate(Species_lower = tolower(Species)) %>%
    left_join(species_mapping, by = c("Species_lower" = "code")) %>%
    select(Taxon = display, Sus = sus, r, r2, R2, RMSE, MAE) %>%
    arrange(match(Taxon, display_names))
  
  combined <- cv_df %>%
    left_join(overall_df, by = "Taxon", suffix = c("_cv", "_all")) %>%
    select(Taxon, r_All = r_all, r_Test = r_cv, r2_All = r2_all, r2_Test = r2_cv,
           R2_All = R2_all, R2_Test = R2_cv, RMSE_All = RMSE_all, RMSE_Test = RMSE_cv,
           MAE_All = MAE_all, MAE_Test = MAE_cv)
  
  return(combined)
}

combined_raw <- format_table(cv_summary_raw_df, overall_summary_raw_df)
combined_calib <- format_table(cv_summary_calib_df, overall_summary_calib_df)

################################## Print Results ##################################

print_table <- function(combined, title) {
  display_table <- combined
  display_table[, 2:11] <- lapply(display_table[, 2:11], function(x) {
    if(all(abs(x) < 1, na.rm = TRUE)) sprintf("%.4f", x) else sprintf("%.2f", x)
  })
  
  cat("\n========================================\n")
  cat(title, "\n")
  cat("========================================\n\n")
  cat(sprintf("%-20s %17s %17s %17s %17s %17s\n",
              "Taxon", "r", "r²", "R²", "RMSE", "MAE"))
  cat(sprintf("%-20s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n",
              "", "All", "Test", "All", "Test", "All", "Test", "All", "Test", "All", "Test"))
  cat(paste(rep("-", 135), collapse = ""), "\n")
  
  for(i in 1:nrow(display_table)) {
    cat(sprintf("%-20s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n",
                display_table$Taxon[i], display_table$r_All[i], display_table$r_Test[i],
                display_table$r2_All[i], display_table$r2_Test[i],
                display_table$R2_All[i], display_table$R2_Test[i],
                display_table$RMSE_All[i], display_table$RMSE_Test[i],
                display_table$MAE_All[i], display_table$MAE_Test[i]))
  }
}

print_table(combined_raw, "RAW PREDICTIONS - MODEL PERFORMANCE")
print_table(combined_calib, "CALIBRATED PREDICTIONS - MODEL PERFORMANCE")

cat("\n✓ Analysis complete\n")
cat("\nNote: 'All' = full dataset; 'Test' = 5-fold CV average\n")
cat("      r = Pearson correlation; r² = r-squared; R² = coefficient of determination\n")

# Save tables
write.csv(combined_raw, here("output", "output_figures_tables", "cv_abundance_raw_table.csv"), row.names = FALSE)
write.csv(combined_calib, here("output", "output_figures_tables", "cv_abundance_calibrated_table.csv"), row.names = FALSE)
cat("\n✓ Tables saved\n")

################################## Scatter Plots ##################################

cat("\n========================================\n")
cat("CREATING SCATTER PLOTS\n")
cat("========================================\n\n")

available_species <- names(scatter_data_list)
set.seed(123)
selected_species <- sample(available_species, min(3, length(available_species)))

cat("Creating scatter plots for:", paste(selected_species, collapse=", "), "\n\n")

for(sp in selected_species) {
  data_list <- scatter_data_list[[sp]]
  display_name <- species_mapping$display[species_mapping$code == sp]
  
  # Get metrics
  metrics_overall_raw <- overall_results_raw_list[[sp]]
  metrics_overall_calib <- overall_results_calib_list[[sp]]
  metrics_test_raw <- cv_results_raw_list[[sp]]
  metrics_test_calib <- cv_results_calib_list[[sp]]
  
  # Create 4 subplots
  plot_max_overall <- max(c(data_list$overall_raw$observed, data_list$overall_raw$predicted,
                            data_list$overall_calib$predicted))
  plot_max_test <- max(c(data_list$test_raw$observed, data_list$test_raw$predicted,
                         data_list$test_calib$predicted))
  
  p1 <- ggplot(data_list$overall_raw, aes(x = predicted, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_point(alpha = 0.4, size = 1, color = "#E69F00") +
    labs(title = "Overall - Raw", x = "Predicted", y = "Observed") +
    annotate("text", x = plot_max_overall * 0.05, y = plot_max_overall * 0.95,
             label = sprintf("r²=%.2f\nR²=%.2f\nN=%d", 
                             metrics_overall_raw["r2"], metrics_overall_raw["R2"], 
                             nrow(data_list$overall_raw)),
             hjust = 0, vjust = 1, size = 2.5) +
    coord_fixed(ratio = 1, xlim = c(0, plot_max_overall), ylim = c(0, plot_max_overall)) +
    theme_classic(base_size = 9)
  
  p2 <- ggplot(data_list$overall_calib, aes(x = predicted, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_point(alpha = 0.4, size = 1, color = "#009E73") +
    labs(title = "Overall - Calibrated", x = "Predicted", y = "Observed") +
    annotate("text", x = plot_max_overall * 0.05, y = plot_max_overall * 0.95,
             label = sprintf("r²=%.2f\nR²=%.2f\nN=%d", 
                             metrics_overall_calib["r2"], metrics_overall_calib["R2"],
                             nrow(data_list$overall_calib)),
             hjust = 0, vjust = 1, size = 2.5) +
    coord_fixed(ratio = 1, xlim = c(0, plot_max_overall), ylim = c(0, plot_max_overall)) +
    theme_classic(base_size = 9)
  
  p3 <- ggplot(data_list$test_raw, aes(x = predicted, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_point(alpha = 0.5, size = 1.2, color = "#E69F00") +
    labs(title = "Test Fold - Raw", x = "Predicted", y = "Observed") +
    annotate("text", x = plot_max_test * 0.05, y = plot_max_test * 0.95,
             label = sprintf("r²=%.2f\nR²=%.2f\nN=%d", 
                             metrics_test_raw["r2"], metrics_test_raw["R2"],
                             nrow(data_list$test_raw)),
             hjust = 0, vjust = 1, size = 2.5) +
    coord_fixed(ratio = 1, xlim = c(0, plot_max_test), ylim = c(0, plot_max_test)) +
    theme_classic(base_size = 9)
  
  p4 <- ggplot(data_list$test_calib, aes(x = predicted, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_point(alpha = 0.5, size = 1.2, color = "#009E73") +
    labs(title = "Test Fold - Calibrated", x = "Predicted", y = "Observed") +
    annotate("text", x = plot_max_test * 0.05, y = plot_max_test * 0.95,
             label = sprintf("r²=%.2f\nR²=%.2f\nN=%d", 
                             metrics_test_calib["r2"], metrics_test_calib["R2"],
                             nrow(data_list$test_calib)),
             hjust = 0, vjust = 1, size = 2.5) +
    coord_fixed(ratio = 1, xlim = c(0, plot_max_test), ylim = c(0, plot_max_test)) +
    theme_classic(base_size = 9)
  
  combined_plot <- (p1 | p2) / (p3 | p4) +
    plot_annotation(title = display_name, 
                    theme = theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 11)))
  
  print(combined_plot)
}

cat("\n✓ Scatter plots created\n")